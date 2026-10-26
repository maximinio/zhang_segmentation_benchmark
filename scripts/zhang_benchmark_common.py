"""Common full-section RNA accounting and honest image-only benchmark diagnostics."""
from pathlib import Path
import json,time,hashlib,resource
import numpy as np
import pandas as pd
import tifffile
from skimage import io,color,filters,morphology,segmentation
from scipy import ndimage as ndi

def load(cfg,sample):
 d=Path(cfg['stage'])/'inputs';g=json.loads((d/f'{sample}.geometry.json').read_text());rgb=io.imread(d/f'{sample}.png')[...,:3];bins=pd.read_csv(d/f'{sample}.bins.csv')
 assert rgb.dtype==np.uint8 and rgb.shape[:2]==tuple(g['shape'][:2]);assert len(bins)==g['total_bins'] and int(bins.nUMI.sum())==g['total_umi']
 return rgb,bins,g

def hematoxylin_and_tissue(rgb):
 a=rgb.astype(np.float32)/255;h=color.rgb2hed(a)[...,0];tissue=(a.mean(2)<.94)&((a.max(2)-a.min(2))>.025)
 tissue=morphology.remove_small_objects(morphology.binary_closing(tissue,morphology.disk(3)),min_size=100)
 return h,tissue

def prepare_inputs(cfg):
 d=Path(cfg['stage'])/'inputs';rows=[]
 for s in cfg['samples']:
  rgb,b,g=load(cfg,s);h,t=hematoxylin_and_tissue(rgb)
  for name,a in [('hematoxylin',h),('common_tissue',t)]:
   p=d/f'{s}.{name}.npy';assert not p.exists();np.save(p,a);assert np.array_equal(np.load(p),a)
  q=float(np.percentile(h[t],99));gray=np.uint8(255*np.clip(h/max(q,np.finfo(float).eps),0,1));io.imsave(d/f'{s}.hematoxylin_for_cyto3.png',gray,check_contrast=False)
  rows.append({'sample':s,'hematoxylin_p99':q,'common_tissue_pixels':int(t.sum()),'same_mask_for_all_methods':True})
 (d/'prepared.qa.json').write_text(json.dumps({'valid':True,'rows':rows},indent=2))

def quantiles(a):
 return dict(zip(['p10','p50','p90','p99','max'],map(float,np.quantile(a,[.1,.5,.9,.99,1])))) if len(a) else {k:None for k in ['p10','p50','p90','p99','max']}

def save_result(cfg,sample,variant,labels,seconds,extra=None):
 started=time.monotonic();out=Path(cfg['output'])/'segmentations'/variant['name'];out.mkdir(parents=True,exist_ok=True)
 rgb,b,g=load(cfg,sample);lab=np.asarray(labels);assert lab.ndim==2 and lab.shape==rgb.shape[:2] and np.all(np.isfinite(lab)) and np.all(lab>=0) and np.all(lab==np.floor(lab))
 lab=lab.astype(np.uint32);p=out/f'{sample}.labels.tiff';assert not p.exists();tifffile.imwrite(p,lab,compression='zlib');assert np.array_equal(tifffile.imread(p),lab)
 px=b.px.to_numpy();py=b.py.to_numpy();x=np.floor(px).astype(int);y=np.floor(py).astype(int);inside=(x>=0)&(x<lab.shape[1])&(y>=0)&(y<lab.shape[0]);a=np.zeros(len(b),np.uint32);a[inside]=lab[y[inside],x[inside]]
 umi=b.nUMI.to_numpy(dtype=np.int64);bylabel=np.bincount(a,weights=umi,minlength=int(lab.max())+1);nbin=np.bincount(a,minlength=len(bylabel));areas=np.bincount(lab.ravel(),minlength=len(bylabel));assert np.all(bylabel==np.floor(bylabel)) and int(bylabel.sum())==int(umi.sum())
 np.save(out/f'{sample}.bin_labels.npy',a);tissue=np.load(Path(cfg['stage'])/'inputs'/f'{sample}.common_tissue.npy')
 # Same original lattice corners for all methods. This is a boundary-ambiguity
 # diagnostic at the available image resolution, not molecule localisation.
 cross=np.zeros(len(b),bool);basis=np.asarray(g['basis'])
 for r,c in [(-.5,-.5),(-.5,.5),(.5,-.5),(.5,.5)]:
  delta=basis@np.array([r,c]);cx=np.floor(px+delta[0]).astype(int);cy=np.floor(py+delta[1]).astype(int);valid=(cx>=0)&(cx<lab.shape[1])&(cy>=0)&(cy<lab.shape[0]);corner=np.zeros(len(b),np.uint32);corner[valid]=lab[cy[valid],cx[valid]];cross|=corner!=a
 object_ids=np.flatnonzero(areas);object_ids=object_ids[object_ids>0];n=len(object_ids);object_umi=bylabel[object_ids];object_areas=areas[object_ids];object_bins=nbin[object_ids];assigned=int(object_umi.sum());inside_umi=int(umi[inside].sum())
 row={'sample':sample,'method':variant['name'],'target':variant.get('target','expanded nuclear candidates' if variant.get('expansion_um',0)>0 else 'nuclear candidates'),
 'model':variant.get('model','classical watershed'),'input':variant.get('input','hematoxylin from original H&E'),'expansion_um':variant.get('expansion_um',0),
 'source_umi':int(umi.sum()),'source_umi_inside_image':inside_umi,'assigned_umi':assigned,'unassigned_umi':int(bylabel[0]),'assigned_umi_pct':100*assigned/int(umi.sum()),'assigned_umi_of_inside_image_pct':100*assigned/inside_umi,
 'total_bins':len(b),'assigned_bins':int((a>0).sum()),'assigned_bins_pct':100*float(np.mean(a>0)),'outside_image_bins':int((~inside).sum()),'outside_image_umi':int(umi[~inside].sum()),
 'objects':n,'objects_with_rna':int(np.sum(object_umi>0)),'objects_zero_rna':int(np.sum(object_umi==0)),'zero_rna_objects_pct':100*float(np.mean(object_umi==0)) if n else None,
 'area_um2':quantiles(object_areas*g['um_per_pixel']**2),'umi_per_object':quantiles(object_umi),'bins_per_object':quantiles(object_bins),
 'image_area_covered_pct':100*float(np.mean(lab>0)),'common_tissue_area_covered_pct':100*float(np.sum((lab>0)&tissue)/np.sum(tissue)),'corner_discordant_bins_pct':100*float(np.mean(cross)),'corner_discordant_umi_pct':100*int(umi[cross].sum())/int(umi.sum()),
 'seconds_inference':seconds,'seconds_common_metrics':time.monotonic()-started,'peak_process_rss_mb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'umi_conserved_exact':True,'mask_readback_exact':True,
 'label_md5':hashlib.md5(p.read_bytes()).hexdigest(),'image_md5':g['image_md5'],'source_md5':cfg['source_md5'],'status':'zero_recovery' if n==0 else 'computed_not_biologically_validated','boundary_accuracy':'not measurable without independent ground truth','molecule_positions_available':False,'configuration':variant,'extra':extra or {}}
 (out/f'{sample}.metrics.json').write_text(json.dumps(row,indent=2));print('BENCHMARK_SECTION_MEASURED',json.dumps({k:row[k] for k in ['sample','method','objects','assigned_umi','assigned_umi_pct','zero_rna_objects_pct','status']}),flush=True)
 return row

if __name__=='__main__':
 import sys
 cfg=json.loads(Path(sys.argv[1]).read_text())
 if sys.argv[2]=='prepare':prepare_inputs(cfg)
 elif sys.argv[2]=='measure':
  name,sample=sys.argv[3:5];v=next(x for x in cfg['learned_variants'] if x['name']==name);d=Path(cfg['stage'])/'native'/name;q=json.loads((d/f'{sample}.native.json').read_text());assert json.loads((d/'inference.exit.json').read_text())['exit_code']==0
  save_result(cfg,sample,v,np.load(d/f'{sample}.native.npy'),q['seconds_inference'],q)
 else:raise ValueError('Unknown benchmark action')
