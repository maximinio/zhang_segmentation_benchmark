"""Native learned segmentation only; common RNA measurement is a separate job."""
from pathlib import Path
import os,sys,json,time,random,hashlib,importlib.metadata as md
cfg=json.loads(Path(sys.argv[1]).read_text());v=next(x for x in cfg['learned_variants'] if x['name']==sys.argv[2]);d=Path(cfg['stage']);out=d/'native'/v['name'];out.mkdir(parents=True,exist_ok=True)
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='4'
for key,sub in [('CELLPOSE_LOCAL_MODELS_PATH','cellpose'),('MICROSAM_CACHEDIR','microsam'),('INSTANSEG_BIOIMAGEIO_PATH','instanseg'),('HF_HOME','huggingface'),('TORCH_HOME','torch')]:os.environ[key]=str(d/'models'/sub)
os.environ['QT_QPA_PLATFORM']='offscreen'
import numpy as np
import torch
import imageio.v3 as iio
from PIL import Image
torch.set_num_threads(4);random.seed(cfg['seed']);np.random.seed(cfg['seed']);torch.manual_seed(cfg['seed']);assert torch.cuda.is_available();start_load=time.monotonic();model=None
if v['name'].startswith('Cellpose'):
 from cellpose import models
 if v['name']=='Cellpose3':model=models.CellposeModel(gpu=True,model_type='cyto3')
 else:
  assert v['model'] in models.MODEL_NAMES
  model=models.CellposeModel(gpu=True,pretrained_model=v['model'],use_bfloat16=torch.cuda.is_bf16_supported())
  assert Path(model.pretrained_model).name==v['model'],'Refuse silent fallback to another Cellpose model'
elif v['name'].startswith('InstanSeg'):
 from instanseg import InstanSeg
 model=InstanSeg(v['model'],device='cuda',image_reader='skimage.io')
else:
 from micro_sam.automatic_segmentation import get_predictor_and_segmenter,automatic_instance_segmentation
 from micro_sam.instance_segmentation import TiledAutomaticMaskGenerator,TiledInstanceSegmentationWithDecoder
 kwargs={'points_per_batch':8} if v['mode']=='amg' else {}
 model=get_predictor_and_segmenter(model_type=v['model'],device='cuda',segmentation_mode=v['mode'],is_tiled=True,**kwargs)
 assert isinstance(model[1],TiledAutomaticMaskGenerator if v['mode']=='amg' else TiledInstanceSegmentationWithDecoder)
 if v['mode']=='amg':model[0].model.mask_threshold=0.
load_seconds=time.monotonic()-start_load
for s in cfg['samples']:
 p=d/'inputs'/f'{s}.png';rgb=np.asarray(Image.open(p).convert('RGB'));g=json.loads((d/'inputs'/f'{s}.geometry.json').read_text());dest=out/f'{s}.native.npy';assert not dest.exists();print('LEARNED_NATIVE_SECTION_START',v['name'],s,flush=True);start=time.monotonic();torch.cuda.reset_peak_memory_stats()
 if v['name']=='Cellpose3':
  image=np.asarray(Image.open(d/'inputs'/f'{s}.hematoxylin_for_cyto3.png'));labels=model.eval(image,channels=[0,0],diameter=v['diameter_px'],batch_size=1)[0]
 elif v['name'].startswith('Cellpose'):labels=model.eval(rgb,channel_axis=-1,diameter=v['diameter_px'],batch_size=1)[0]
 elif v['name'].startswith('InstanSeg'):
  target='cells' if v['model']=='fluorescence_nuclei_and_cells' else 'nuclei'
  labels=model.eval_medium_image(torch.from_numpy(rgb.transpose(2,0,1).copy()).float(),pixel_size=g['um_per_pixel'],target=target,return_image_tensor=False,rescale_output=True,tile_size=512,batch_size=1);labels=labels.cpu().numpy().squeeze()
 else:
  image=np.asarray(Image.open(p).convert('L' if v['mode']=='amg' else 'RGB'));input_path=out/f'{s}.model_input.tif';iio.imwrite(input_path,image)
  kwargs={'pred_iou_thresh':.88,'stability_score_thresh':.95} if v['mode']=='amg' else {'output_mode':'instance_segmentation'}
  tile=(1200,1200) if v['mode']=='amg' else (384,384);halo=(128,128) if v['mode']=='amg' else (64,64)
  labels=automatic_instance_segmentation(*model,str(input_path),embedding_path=str(out/f'{s}.embeddings.zarr'),ndim=2,tile_shape=tile,halo=halo,verbose=True,**kwargs)
 seconds=time.monotonic()-start;labels=np.asarray(labels);assert labels.ndim==2 and labels.shape==rgb.shape[:2] and np.all(np.isfinite(labels)) and np.all(labels>=0) and np.all(labels==np.floor(labels))
 np.save(dest,labels);assert np.array_equal(np.load(dest),labels)
 q={'sample':s,'method':v['name'],'configuration':v,'seconds_inference':seconds,'seconds_model_load':load_seconds,'peak_gpu_allocated_mb':torch.cuda.max_memory_allocated()/2**20,'native_mask_readback_exact':True,'native_mask_md5':hashlib.md5(dest.read_bytes()).hexdigest(),'shape':list(labels.shape),'model_output_not_thresholded_or_expanded':True,'versions':{n:md.version(n) for n in (['cellpose'] if v['name'].startswith('Cellpose') else ['instanseg-torch'] if v['name'].startswith('InstanSeg') else ['micro-sam'])}}
 q['inference_geometry']={'um_per_pixel':g['um_per_pixel'],'native_image_shape':list(rgb.shape),'rescale_output_to_original_frame':True if v['name'].startswith('InstanSeg') else None}
 q['cellpose_precision']='bfloat16' if torch.cuda.is_bf16_supported() else 'float32' if v['name'].startswith('Cellpose') and v['name']!='Cellpose3' else None
 if v['name'].startswith('Cellpose'):q['actual_checkpoint']=str(model.pretrained_model)
 if v['name'].startswith('InstanSeg'):q['native_parameters']={'tile_size':512,'batch_size':1,'target':target,'pixel_size':g['um_per_pixel'],'rescale_output':True}
 elif not v['name'].startswith('Cellpose'):q['native_parameters']={'tile_shape':tile,'halo':halo,**kwargs}
 (out/f'{s}.native.json').write_text(json.dumps(q,indent=2));print('LEARNED_NATIVE_SECTION_SAVED',v['name'],s,'positive_pixels',int(np.sum(labels>0)),flush=True);del labels,rgb;torch.cuda.empty_cache()
print('LEARNED_NATIVE_BOTH_SECTIONS_FINISHED',v['name'],flush=True)
