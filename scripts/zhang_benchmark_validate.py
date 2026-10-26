"""Independent full-section accounting audit, not biological validation."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
import pandas as pd
import tifffile

def md5(p):
 h=hashlib.md5()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

cfg=json.loads(Path(sys.argv[1]).read_text());d=Path(cfg['stage']);out=Path(cfg['output'])
variants=cfg['learned_variants']+cfg['classical_variants']
assert len({v['name'] for v in variants})==len(variants)
assert md5(Path(cfg['source']))==cfg['source_md5']
rows=[]
for s in cfg['samples']:
 b=pd.read_csv(d/'inputs'/f'{s}.bins.csv');g=json.loads((d/'inputs'/f'{s}.geometry.json').read_text())
 umi=b.nUMI.to_numpy(np.int64);x=np.floor(b.px).astype(int);y=np.floor(b.py).astype(int)
 for v in variants:
  name=v['name'];assert json.loads((d/'jobs'/name/'method.exit.json').read_text())['exit_code']==0,name
  folder=out/'segmentations'/name;q=json.loads((folder/f'{s}.metrics.json').read_text());p=folder/f'{s}.labels.tiff';lab=tifffile.imread(p)
  assert md5(p)==q['label_md5'];assert lab.shape==tuple(g['shape'][:2])
  if name in {z['name'] for z in cfg['learned_variants']}:
   native=d/'native'/name;assert json.loads((native/'inference.exit.json').read_text())['exit_code']==0
   nq=json.loads((native/f'{s}.native.json').read_text());npth=native/f'{s}.native.npy'
   assert md5(npth)==nq['native_mask_md5'];assert np.array_equal(lab,np.load(npth)),name
  inside=(x>=0)&(y>=0)&(x<lab.shape[1])&(y<lab.shape[0]);assigned=np.zeros(len(b),np.uint32);assigned[inside]=lab[y[inside],x[inside]]
  assert np.array_equal(assigned,np.load(folder/f'{s}.bin_labels.npy'))
  captured=int(umi[assigned>0].sum());lost=int(umi[assigned==0].sum());n=len(np.unique(lab[lab>0]))
  assert captured==q['assigned_umi'] and lost==q['unassigned_umi'] and n==q['objects']
  assert captured+lost==q['source_umi']==g['total_umi'];assert q['source_md5']==cfg['source_md5']
  assert abs(q['assigned_umi_pct']-100*captured/(captured+lost))<1e-10
  rows.append({'sample':s,'method':name,'assigned_umi':captured,'unassigned_umi':lost,'objects':n,'label_md5':q['label_md5']})
  print('INDEPENDENT_ACCOUNTING_PASS',s,name,captured,lost,n,flush=True)
summary=pd.read_csv(out/'tables/segmentation_metrics.csv');assert len(summary)==len(rows)
assert not summary.status.isin(['pending','failed']).any()
assert all((out/'figures/supplementary/segmentation'/f'{s}_comparison_{i:02d}.png').is_file() for s in cfg['samples'] for i in range(1,(len(variants)+4)//5+1))
assert (out/'figures/main/metric_summary.png').is_file()
qa={'computational_validation_passed':True,'biological_accuracy_validated':False,'visual_review_required':True,'source_md5':cfg['source_md5'],'required_methods':len(variants),'verified_full_sections':len(rows),'whole_bin_centre_assignment_not_molecule_positions':True,'rows':rows}
(out/'independent_validation.json').write_text(json.dumps(qa,indent=2));print('BENCHMARK_ALL_METHODS_ACCOUNTING_VALIDATED',len(rows),flush=True)
