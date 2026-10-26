"""Read accepted masks and assignments; never alter or re-segment a source mask."""
from pathlib import Path
import sys,json,hashlib,csv
import numpy as np
import pandas as pd
import tifffile

def digest(p,kind='md5'):
    h=hashlib.new(kind)
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def run(config):
    cfg=json.loads(Path(config).read_text());root=Path(cfg['archive']);stage=Path(cfg['stage'])
    out=stage/'assignments';out.mkdir(exist_ok=False)
    methods=sorted(p for p in (root/'segmentations').iterdir() if p.is_dir())
    assert len(methods)==20
    jobs=[]
    for method in methods:
        for sample in cfg['samples']:
            maskfile=method/f'{sample}.labels.tiff';assignment=method/f'{sample}.bin_labels.npy'
            metrics=json.loads((method/f'{sample}.metrics.json').read_text())
            assert digest(maskfile)==metrics['label_md5'] and metrics['source_md5']==cfg['source_md5']
            mask=tifffile.imread(maskfile);a=np.load(assignment,allow_pickle=False)
            bins=pd.read_csv(root/'provenance'/f'{sample}.bins.csv')
            assert a.ndim==1 and len(a)==len(bins) and len(bins)==metrics['total_bins']
            assert np.issubdtype(mask.dtype,np.integer) and mask.min()>=0
            px=np.floor(bins.px.to_numpy()).astype(np.int64);py=np.floor(bins.py.to_numpy()).astype(np.int64)
            inside=(px>=0)&(py>=0)&(px<mask.shape[1])&(py<mask.shape[0]);exact=np.zeros(len(a),dtype=a.dtype)
            exact[inside]=mask[py[inside],px[inside]]
            assert np.array_equal(exact,a),f'Assignment differs from accepted mask: {method.name}/{sample}'
            assert int(bins.nUMI.to_numpy()[a>0].sum())==metrics['assigned_umi']
            areas=np.bincount(mask.ravel());labels=np.flatnonzero(areas);labels=labels[labels>0]
            assert len(labels)==metrics['objects'] and mask.max()<np.iinfo(np.int32).max
            sx=np.bincount(mask.ravel(),weights=np.tile(np.arange(mask.shape[1],dtype=np.float64),mask.shape[0]),minlength=len(areas))
            sy=np.bincount(mask.ravel(),weights=np.repeat(np.arange(mask.shape[0],dtype=np.float64),mask.shape[1]),minlength=len(areas))
            job=f'{method.name}__{sample}';folder=out/job;folder.mkdir()
            labelsfile=folder/'labels.i32';a.astype('<i4').tofile(labelsfile)
            regionfile=folder/'regions.csv'
            with regionfile.open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['region_id','source_label','area_hires_px2','centroid_x_hires','centroid_y_hires'])
                for label in labels:w.writerow([f'{job}__{label}',int(label),int(areas[label]),float(sx[label]/areas[label]),float(sy[label]/areas[label])])
            jobs.append(dict(job=job,method=method.name,sample=sample,mask=str(maskfile),mask_md5=metrics['label_md5'],assignment=str(assignment),assignment_md5=digest(assignment),labels=str(labelsfile),labels_md5=digest(labelsfile),regions=str(regionfile),regions_md5=digest(regionfile),regions_n=int(len(labels)),source_bins=int(len(a)),source_umi=metrics['source_umi'],assigned_umi=metrics['assigned_umi'],target=metrics['target'],mask_shape=list(mask.shape)))
            print('ACCEPTED_MASK_ASSIGNMENTS_EXACT',job,len(labels),metrics['assigned_umi'],flush=True)
    assert len(jobs)==40
    (stage/'jobs.json').write_text(json.dumps({'jobs':jobs,'all_40_source_assignments_exact':True,'no_segmentation_performed':True},indent=2))
    print('ALL_40_MASKS_PREPARED',flush=True)

if __name__=='__main__':run(sys.argv[1])
