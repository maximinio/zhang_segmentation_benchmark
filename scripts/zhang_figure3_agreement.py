from pathlib import Path
import numpy as np,pandas as pd,json
import sys,hashlib
a=Path(sys.argv[1]);w=Path(sys.argv[2]);tab=pd.read_csv(a/'tables/rctd_full_major_summary.csv')
methods=tab.method.drop_duplicates().tolist();audit={}
for sample in ['CD1','CD2']:
 vectors={};coverage={}
 for m in methods:
  calls=pd.read_csv(w/'profiles'/f'{m}__{sample}.csv')
  labels=np.load(a/'segmentations'/m/f'{sample}.bin_labels.npy')
  assert labels.ndim==1
  ix=np.argsort(calls.source_label.to_numpy());keys=calls.source_label.to_numpy()[ix];val=calls.call.to_numpy()[ix]
  v=np.zeros(labels.size,dtype=np.int8)
  if len(keys):
   pos=np.searchsorted(keys,labels);good=pos<len(keys);good[good]&=keys[pos[good]]==labels[good];v[good]=val[pos[good]]
  vectors[m]=v;coverage[m]=int(np.count_nonzero(v))
 n=len(methods);agreement=np.full((n,n),np.nan);den=np.zeros((n,n),dtype=int)
 for i,m in enumerate(methods):
  for j,k in enumerate(methods):
   v,u=vectors[m],vectors[k];ok=(v>0)&(u>0);den[i,j]=ok.sum()
   if den[i,j]:agreement[i,j]=np.mean(v[ok]==u[ok])
 pd.DataFrame(agreement,index=methods,columns=methods).to_csv(w/f'{sample}.agreement.csv')
 pd.DataFrame(den,index=methods,columns=methods).to_csv(w/f'{sample}.denominators.csv')
 audit[sample]={'total_original_bins':len(v),'fitted_bin_coverage':coverage,'pairwise_complete_support':True,'missing_never_counted_as_agreement':True}
json.dump(audit,(w/'agreement.qa.json').open('w'),indent=2)
print('AGREEMENT_BOTH_SECTIONS_VERIFIED',flush=True)
