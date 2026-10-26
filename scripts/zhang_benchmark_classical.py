from pathlib import Path
import sys,json,time,os
for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[n]='4'
import numpy as np
from scipy import ndimage as ndi
from skimage import filters,morphology,segmentation,feature
import cv2
from zhang_benchmark_common import load,hematoxylin_and_tissue,save_result
cfg=json.loads(Path(sys.argv[1]).read_text());name=sys.argv[2];v=next(x for x in cfg['classical_variants'] if x['name']==name);cv2.setNumThreads(4)
for s in cfg['samples']:
 rgb,b,g=load(cfg,s);um=g['um_per_pixel'];start=time.monotonic();h,tissue=hematoxylin_and_tissue(rgb)
 if v['backend']=='zhang':
  sm=filters.gaussian(h,sigma=.7/um,preserve_range=True);threshold=filters.threshold_otsu(sm[tissue]);fg=(sm>threshold)&tissue;minimum=max(3,round(12/um**2));fg=morphology.remove_small_objects(fg,min_size=minimum)
  fg=ndi.binary_fill_holes(fg) if v['holes']=='all' else morphology.remove_small_holes(fg,area_threshold=minimum);distance=ndi.distance_transform_edt(fg);spacing=max(1,round(v['spacing_um']/um));peaks=feature.peak_local_max(distance,min_distance=spacing,labels=fg,exclude_border=False)
 else:
  sm=filters.gaussian(h,sigma=1,preserve_range=True);threshold=filters.threshold_otsu(sm[tissue]);fg=(sm>threshold)&tissue;fg=morphology.remove_small_objects(fg,min_size=15);distance=ndi.distance_transform_edt(fg);spacing=v['spacing_px'];peaks=feature.peak_local_max(distance,min_distance=spacing,labels=fg)
 markers=np.zeros(fg.shape,np.int32);markers[tuple(peaks.T)]=np.arange(1,len(peaks)+1)
 if v['backend']=='opencv':
  markers=markers+1;markers[fg&(markers==1)]=0;labels=cv2.watershed(rgb.copy(),markers);labels[labels<=1]=0;labels[labels>1]-=1;labels[~fg]=0
 else:labels=segmentation.watershed(-distance,markers,mask=fg)
 if v['expansion_um']>0:labels=segmentation.expand_labels(labels,distance=v['expansion_um']/um);labels[~tissue]=0
 save_result(cfg,s,v,labels,time.monotonic()-start,{'effective_seed_spacing_px':spacing,'hematoxylin_threshold':float(threshold),'holes':v['holes'],'no_cell_identity_inferred':True})
print('CLASSICAL_METHOD_EXIT_SUCCESS',name,flush=True)
