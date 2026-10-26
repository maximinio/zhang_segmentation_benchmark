"""Recolor the exact original S1-S8 raster layout with native RCTD argmax calls.

No new segmentation, fitting, threshold, probability normalization or uncertainty
analysis. The original Python renderer is reused to preserve the user's layout.
"""
from pathlib import Path
import csv,json,hashlib,sys,subprocess,shutil,time,os
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from skimage import io
from PIL import Image

COLORS={'Myeloid':'#D55E00','CD4T':'#0072B2','CD8T':'#56B4E9','ILC':'#CC79A7',
 'B_Plasma':'#882255','Epithelial':'#E69F00','Mesenchymal':'#009E73',
 'Endothelial':'#44AA99','Neural':'#AA4499','Not fitted':'#B3B3B3','Tied maximum':'#222222'}

def md5(p):
 h=hashlib.md5()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()

def render(archive,work,out):
 cfg=json.loads((archive/'provenance/zhang_benchmark_config.json').read_text())
 rois=json.loads((archive/'provenance/visual_fields.json').read_text())
 variants=cfg['learned_variants']+cfg['classical_variants'];assert len(variants)==20
 summary={x['job']:x for x in csv.DictReader((archive/'tables/rctd_full_major_summary.csv').open())}
 font=subprocess.check_output(['fc-match','Arial','--format','%{family}'],text=True).split(',')[0]
 assert font=='Arial';plt.rcParams['font.family']=font
 names=list(COLORS);cmap=ListedColormap(list(COLORS.values()));audit=[]
 for sample in cfg['samples']:
  rgb=io.imread(work/f'{sample}.png');geom=json.loads((archive/f'provenance/{sample}.geometry.json').read_text())
  assert md5(work/f'{sample}.png')==geom['image_md5']
  for page in range(4):
   group=variants[page*5:page*5+5]
   fig,axes=plt.subplots(5,3,figsize=(16,22),layout='constrained',squeeze=False)
   for row,v in enumerate(group):
    method=v['name'];job=f'{method}__{sample}';entry=summary[job]
    maskpath=archive/'segmentations'/method/f'{sample}.labels.tiff'
    metrics=json.loads((archive/'segmentations'/method/f'{sample}.metrics.json').read_text())
    assert md5(maskpath)==metrics['label_md5']
    mask=tifffile.imread(maskpath);assert mask.shape==rgb.shape[:2]
    records=list(csv.DictReader((work/'calls'/f'{job}.csv').open()))
    labels=np.unique(mask[mask>0]);assert set(labels)=={int(x['source_label']) for x in records}
    lookup=np.zeros(int(mask.max())+1,dtype=np.uint8)
    for x in records:lookup[int(x['source_label'])]=names.index(x['display_major'])+1
    assert len(records)==int(entry['input_regions'])
    categories=lookup[mask]
    count_fitted=sum(x['display_major']!='Not fitted' for x in records)
    assert count_fitted==int(entry['native_retained_regions'])
    for col,ax in enumerate(axes[row]):
     if col==0:image=rgb;field=categories;fieldname='Full section'
     else:
      q=rois[sample][col-1];sl=np.s_[q['y']:q['y']+q['size'],q['x']:q['x']+q['size']]
      image=rgb[sl];field=categories[sl];fieldname=q['label']
     ax.imshow(image,interpolation='nearest');ax.set_axis_off()
     # Direct pixel-to-label lookup, not coordinate fitting or rerasterization.
     ax.imshow(np.ma.masked_where(field==0,field-1),cmap=cmap,vmin=0,vmax=len(names)-1,interpolation='nearest',alpha=.70)
     scale_um=1000 if col==0 else 100;length=scale_um/geom['um_per_pixel'];margin=.04*image.shape[1];yb=.94*image.shape[0]
     ax.plot([margin,margin+length],[yb,yb],color='black',linewidth=4)
     ax.plot([margin,margin+length],[yb,yb],color='white',linewidth=2)
     ax.text(margin,yb-.02*image.shape[0],f'{scale_um} µm',color='black',fontsize=9,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
     title=method.replace('_',' ')+' · '+fieldname
     if col==0:title+='\n'+(f'{count_fitted:,} fitted / {len(records):,} regions' if records else 'No regions — no RCTD fit')
     ax.set_title(title,fontsize=11)
    audit.append({'job':job,'regions':len(records),'native_fitted':count_fitted,'ties':sum(x['display_major']=='Tied maximum' for x in records),'mask_md5':metrics['label_md5'],'bundle_md5':entry['md5']})
    del mask,categories,lookup,records
    print('PANELS_READY',job,flush=True)
   fig.suptitle(f'Zhang {sample} · RCTD major display call · matched S1–S8 fields · page {page+1}\nColor = largest native full-mode coefficient; not calibrated probability. Uncertainty analysis deferred.',fontsize=14)
   # Native constrained-layout canvas allocation reserves a legend margin;
   # Matplotlib 3.6 does not yet implement loc='outside lower center'.
   fig.set_layout_engine('constrained',rect=(0,.035,1,.965))
   fig.legend(handles=[Patch(facecolor=COLORS[n],label=n.replace('_',' / ')) for n in names],loc='lower center',ncol=6,fontsize=10,frameon=False)
   target=out/f'{sample}_major_calls_{page+1:02d}.png'
   fig.savefig(target,dpi=300);plt.close(fig)
   with Image.open(target) as im:assert im.size==(4800,6600) and abs(im.info['dpi'][0]-300)<.1
   print('FIGURE_WRITTEN',target,target.stat().st_size,flush=True)
 qa={'status':'RENDERED_PENDING_VISUAL_QA','dpi':300,'font':font,'layout':'Original 5 methods x 3 views; 16 x 22 inches; 8 pages',
  'fields':rois,'major_palette':COLORS,'display_rule':'Unique maximum native full-mode coefficient; exact ties neutral; no clipping or renormalization',
  'no_new_rctd_or_segmentation':True,'uncertainty_analysis_performed':False,'all_40_scopes':audit,
  'figures':[{ 'path':str(p),'md5':md5(p),'bytes':p.stat().st_size} for p in sorted(out.glob('*.png'))]}
 assert len(audit)==40 and len(qa['figures'])==8
 (out/'render.qa.json').write_text(json.dumps(qa,indent=2))

def worker(archive,work,out):
 code=work/'frozen_code';provenance=out/'provenance';provenance.mkdir()
 for source in code.iterdir():shutil.copy2(source,provenance/source.name)
 commands=[['/usr/bin/Rscript',str(code/'zhang_major_call_supplement.R'),str(archive),str(work)]]
 for command in commands:
  start=time.time();p=subprocess.Popen(command,stdout=open(work/'extract.log','w'),stderr=subprocess.STDOUT)
  (work/'child.json').write_text(json.dumps({'pid':p.pid,'argv':command}))
  rc=p.wait();(provenance/'extract.exit.json').write_text(json.dumps({'argv':command,'pid':p.pid,'started':start,'ended':time.time(),'exit_code':rc}))
  if rc:raise RuntimeError('Native extraction failed; no figures accepted')
 render(archive,work,out)
 shutil.copy2(work/'extract.log',provenance/'extract.log')
 # All per-region calls are deterministic display intermediates, not new fits.
 shutil.rmtree(work/'calls')
 for s in ['CD1','CD2']:(work/f'{s}.png').unlink()

if __name__=='__main__':
 archive=Path(sys.argv[2]);work=Path(sys.argv[3]);out=archive/'figures'/'supplementary'/'rctd'
 if sys.argv[1]=='launch':
  assert not out.exists() and not work.exists();work.mkdir(parents=True);out.mkdir()
  code=work/'frozen_code';code.mkdir()
  for name in ['zhang_major_call_supplement.py','zhang_major_call_supplement.R']:
   shutil.copy2(Path(__file__).with_name(name),code/name)
  p=subprocess.Popen(['/usr/bin/python3','-u',str(code/Path(__file__).name),'worker',str(archive),str(work)],stdout=open(work/'render.log','w'),stderr=subprocess.STDOUT,start_new_session=True)
  (work/'supervisor.json').write_text(json.dumps({'pid':p.pid,'cmdline':str(code/Path(__file__).name)}));print('SUPERVISOR',p.pid,flush=True)
 else:
  start=time.time();rc=0
  try:worker(archive,work,out)
  except Exception:rc=1;raise
  finally:(work/'exit.json').write_text(json.dumps({'pid':os.getpid(),'started':start,'ended':time.time(),'exit_code':rc}))
