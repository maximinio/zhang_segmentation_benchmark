"""Matched H&E panels and full-section metrics; never call coverage accuracy."""
from pathlib import Path
import json,sys,math,hashlib,subprocess
import numpy as np
import pandas as pd
from skimage import io,filters,segmentation
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def report(cfg,partial=False):
 d=Path(cfg['stage']);out=Path(cfg['output']);figdir=out/'figures/supplementary/segmentation';figdir.mkdir(parents=True,exist_ok=True);(out/'figures/main').mkdir(parents=True,exist_ok=True);(out/'tables').mkdir(exist_ok=True);variants=cfg['learned_variants']+cfg['classical_variants'];rows=[];completed=[]
 font=subprocess.check_output(['fc-match','Arial','--format','%{family}'],text=True).split(',')[0];plt.rcParams['font.family']=font
 for v in variants:
  p=d/'jobs'/v['name']/'method.exit.json';state='pending' if not p.exists() else 'failed' if json.loads(p.read_text())['exit_code'] else 'completed'
  if state=='completed':
   completed.append(v)
   for s in cfg['samples']:
    r=json.loads((out/'segmentations'/v['name']/f'{s}.metrics.json').read_text());assert r['umi_conserved_exact'] and r['mask_readback_exact'];rows.append(r)
  else:
   for s in cfg['samples']:rows.append({'sample':s,'method':v['name'],'status':state,'target':v.get('target','expanded nuclear candidates' if v.get('expansion_um',0)>0 else 'nuclear candidates'),'configuration':v})
 flat=[]
 for r in rows:
  a={k:v for k,v in r.items() if not isinstance(v,(dict,list))}
  for name in ['area_um2','umi_per_object','bins_per_object']:
   for k,v in r.get(name,{}).items():a[name+'_'+k]=v
  a['configuration_json']=json.dumps(r['configuration']);flat.append(a)
  if 'extra' in r:
   a['peak_inference_gpu_allocated_mb']=r['extra'].get('peak_gpu_allocated_mb')
   a['seconds_model_load']=r['extra'].get('seconds_model_load')
   a['peak_process_rss_scope']='common RNA measurement only' if r['extra'].get('native_mask_readback_exact') else 'segmentation plus common RNA measurement'
 pd.DataFrame(flat).to_csv(out/('partial_metrics.csv' if partial else 'tables/segmentation_metrics.csv'),index=False,encoding='utf-8-sig')
 roipath=d/'inputs/visual_fields.json'
 if roipath.exists():rois=json.loads(roipath.read_text())
 else:
  rois={}
  for s in cfg['samples']:
   h=np.load(d/'inputs'/f'{s}.hematoxylin.npy');t=np.load(d/'inputs'/f'{s}.common_tissue.npy');threshold=filters.threshold_otsu(h[t]);candidates=[];size=384
   for y in range(0,h.shape[0]-size+1,size):
    for x in range(0,h.shape[1]-size+1,size):
     sl=np.s_[y:y+size,x:x+size]
     if t[sl].mean()>=.5:candidates.append((float((h[sl]>threshold).mean()),y,x))
   assert candidates;candidates.sort();selected=[candidates[-1],candidates[len(candidates)//2]];rois[s]=[{'label':label,'x':x,'y':y,'size':size,'source_hematoxylin_density':density} for label,(density,y,x) in zip(['Dense field','Medium-density field'],selected)]
  roipath.write_text(json.dumps(rois,indent=2))
 active=completed if partial else variants
 for sample in cfg['samples']:
  rgb=io.imread(d/'inputs'/f'{sample}.png');g=json.loads((d/'inputs'/f'{sample}.geometry.json').read_text())
  for page in range(math.ceil(len(active)/5)):
   group=active[5*page:5*page+5];fig,axes=plt.subplots(len(group),3,figsize=(16,4.4*len(group)),layout='constrained',squeeze=False)
   for r,v in enumerate(group):
    entry=next(z for z in rows if z['sample']==sample and z['method']==v['name']);p=out/'segmentations'/v['name']/f'{sample}.labels.tiff';labels=None
    if entry['status'] not in ['pending','failed']:
     assert hashlib.md5(p.read_bytes()).hexdigest()==entry['label_md5'];labels=tifffile.imread(p)
    for c,ax in enumerate(axes[r]):
     if c==0:image=rgb;mask=labels;field='Full section'
     else:
      q=rois[sample][c-1];sl=np.s_[q['y']:q['y']+q['size'],q['x']:q['x']+q['size']];image=rgb[sl];mask=labels[sl] if labels is not None else None;field=q['label']
     ax.imshow(image,interpolation='nearest');ax.set_axis_off()
     if mask is not None:
      edges=segmentation.find_boundaries(mask,mode='outer');ax.imshow(np.ma.masked_where(~edges,edges),cmap=ListedColormap(['#00ef70']),vmin=0,vmax=1,interpolation='nearest',alpha=.85)
     scale_um=1000 if c==0 else 100;length=scale_um/g['um_per_pixel'];margin=.04*image.shape[1];yb=.94*image.shape[0]
     ax.plot([margin,margin+length],[yb,yb],color='black',linewidth=4);ax.plot([margin,margin+length],[yb,yb],color='white',linewidth=2)
     ax.text(margin,yb-.02*image.shape[0],f'{scale_um} µm',color='black',fontsize=9,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
     title=v['name'].replace('_',' ')+' · '+field
     if c==0:title+='\n'+(f"{entry['objects']:,} regions · {entry['assigned_umi_pct']:.2f}% original UMI" if labels is not None else entry['status'].upper())
     ax.set_title(title,fontsize=11)
   fig.suptitle(f"Zhang {sample} · same original H&E and matched fields · page {page+1}\nAll detected boundaries shown; nuclei and expanded/cellular targets differ. RNA coverage is not segmentation accuracy.",fontsize=14)
   fig.savefig(figdir/f"{'partial_' if partial else ''}{sample}_comparison_{page+1:02d}.png",dpi=300);plt.close(fig)
 summary=[]
 for v in variants:
  pair=[r for r in rows if r['method']==v['name']];good=all('assigned_umi' in r for r in pair);denom=sum(r.get('source_umi',0) for r in pair);a=pair[0];b=pair[1]
  summary.append([v['name'].replace('_',' '),'Nuclei +5um' if v.get('expansion_um',0)>0 else v.get('target','Nuclei'),f"{a['assigned_umi_pct']:.2f}" if good else a['status'],f"{b['assigned_umi_pct']:.2f}" if good else b['status'],f"{100*sum(r['assigned_umi'] for r in pair)/denom:.2f}" if good else '—',f"{a['objects']:,} / {b['objects']:,}" if good else '—',f"{a['area_um2']['p50']:.1f} / {b['area_um2']['p50']:.1f}" if good and a['objects'] and b['objects'] else '—',f"{a['zero_rna_objects_pct']:.1f} / {b['zero_rna_objects_pct']:.1f}" if good and a['objects'] and b['objects'] else '—',f"{a['corner_discordant_umi_pct']:.1f} / {b['corner_discordant_umi_pct']:.1f}" if good else '—'])
 fig,ax=plt.subplots(figsize=(20,12),layout='constrained');ax.axis('off');table=ax.table(cellText=summary,colLabels=['Method','Target','UMI% CD1','UMI% CD2','UMI% pooled','Regions\nCD1 / CD2','Median area um²\nCD1 / CD2','No assigned UMI\nregions% CD1 / CD2','Bin-corner discordance\nUMI% CD1 / CD2'],cellLoc='left',loc='center',colWidths=[.21,.18,.06,.06,.06,.13,.10,.10,.10]);table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,2)
 ax.set_title('Zhang full-section benchmark · original RNA denominators\nNuclear, cellular and expanded targets differ; no boundary accuracy or automatic winner.\nCorner discordance diagnoses 8um-bin assignment ambiguity, not the fraction of wrongly assigned molecules.',fontsize=14,pad=20);fig.savefig((out/'figures/main')/('partial_metric_summary.png' if partial else 'metric_summary.png'),dpi=300);plt.close(fig)
 q={'all_required_methods_complete':len(completed)==len(variants),'completed_methods':len(completed),'required_methods':len(variants),'full_section_metrics':True,'font':font,'dpi':300,'visual_fields_fixed_independent_of_method':rois,'no_accuracy_claim':True,'status':'PARTIAL' if partial else 'COMPUTED_AWAITING_VISUAL_REVIEW'};(out/('partial_report.qa.json' if partial else 'report.qa.json')).write_text(json.dumps(q,indent=2));print('BENCHMARK_REPORT_WRITTEN',q['status'],len(completed),len(variants),flush=True)
 if not partial:assert q['all_required_methods_complete'],'Final report cannot omit failed/pending methods'

if __name__=='__main__':report(json.loads(Path(sys.argv[1]).read_text()),partial='--partial' in sys.argv[2:])
