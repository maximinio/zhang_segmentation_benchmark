"""Lossless polygonization of the user's accepted labels; no segmentation rerun."""
from pathlib import Path
import sys,json,hashlib,csv,time
import numpy as np
import tifffile
from rasterio.features import shapes
from affine import Affine
from shapely.geometry import shape,MultiPolygon,mapping

def md5(p):
 h=hashlib.md5()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def polygonize(mask):
 """Keep exact pixel corners, all parts and holes, connectivity 4, no simplification."""
 parts={}
 for g,value in shapes(mask.astype(np.int32),mask=mask>0,connectivity=4,transform=Affine.identity()):
  label=int(value);parts.setdefault(label,[]).append(shape(g))
 result={}
 for label,ps in parts.items():
  g=ps[0] if len(ps)==1 else MultiPolygon(ps)
  if not g.is_valid:raise ValueError(f'Invalid exact source polygon {label}')
  result[label]=g
 return result

def self_test():
 a=np.array([[1,1,1,0,2],[1,0,1,2,0],[1,1,1,0,0]],dtype=np.int32)
 g=polygonize(a)
 assert g[1].area==8 and len(g[1].interiors)==1
 assert g[2].area==2 and len(g[2].geoms)==2
 print('POLYGONIZER_HOLES_MULTIPART_EXACT_TEST_PASS',flush=True)

def main(config):
 cfg=json.loads(Path(config).read_text());root=Path(cfg['archive']);out=Path(cfg['stage'])/'geometry'
 out.mkdir(parents=True,exist_ok=False);qa={}
 for sample in cfg['samples']:
  start=time.time();source=root/'segmentations'/cfg['method']/f'{sample}.labels.tiff'
  assert md5(source)==cfg['expected'][sample]['mask_md5']
  mask=tifffile.imread(source);assert mask.ndim==2 and np.issubdtype(mask.dtype,np.integer)
  labels,areas=np.unique(mask[mask>0],return_counts=True)
  assert len(labels)==cfg['expected'][sample]['polygons']
  print('POLYGONIZING_ACCEPTED_MASK',sample,len(labels),flush=True)
  gs=polygonize(mask);assert set(gs)==set(map(int,labels))
  target=out/f'{sample}.polygons.geojson'
  with target.open('w') as f:
   f.write('{"type":"FeatureCollection","features":[')
   for i,(label,area) in enumerate(zip(labels,areas)):
    g=gs[int(label)];assert g.area==int(area)
    if i:f.write(',')
    json.dump({'type':'Feature','properties':{'label':int(label),'barcodes':f'{sample}_region_{int(label):07d}','area_hires_px2':int(area)},'geometry':mapping(g)},f,separators=(',',':'))
   f.write(']}')
  archived=np.load(root/'segmentations'/cfg['method']/f'{sample}.bin_labels.npy',allow_pickle=False)
  assignment=out/f'{sample}.bin_assignment.csv';n=0;assigned_umi=0;source_umi=0
  with (root/'provenance'/f'{sample}.bins.csv').open(newline='') as src,assignment.open('w',newline='') as dest:
   reader=csv.DictReader(src);writer=csv.writer(dest);writer.writerow(['barcode','label','x_fullres','y_fullres','inside_image'])
   for n,row in enumerate(reader,1):
    x=float(row['px']);y=float(row['py']);inside=0<=x<mask.shape[1] and 0<=y<mask.shape[0]
    label=int(mask[int(np.floor(y)),int(np.floor(x))]) if inside else 0
    assert label==int(archived[n-1]);umi=int(float(row['nUMI']));source_umi+=umi
    if label>0:assigned_umi+=umi
    writer.writerow([row['cell'],label,row['x'],row['y'],int(inside)])
  assert n==len(archived)
  assert source_umi==cfg['expected'][sample]['source_umi'] and assigned_umi==cfg['expected'][sample]['assigned_umi']
  qa[sample]={'polygons':len(labels),'pixel_areas_exact':True,'all_parts_and_holes_preserved':True,'source_mask_md5':md5(source),'geometry_md5':md5(target),'bin_assignment_md5':md5(assignment),'bin_assignments_identical_to_benchmark':True,'source_bins':n,'source_umi':source_umi,'assigned_umi':assigned_umi,'seconds':time.time()-start}
  print('EXACT_POLYGONIZATION_PASS',sample,qa[sample],flush=True)
 (Path(cfg['stage'])/'polygonization.qa.json').write_text(json.dumps(qa,indent=2))

if __name__=='__main__':
 self_test()
 if len(sys.argv)>1 and sys.argv[1]!='--self-test':main(sys.argv[1])
