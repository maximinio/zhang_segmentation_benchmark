"""Independent frozen benchmark runner; sequential full fits, receipts, no Drive writes."""
from pathlib import Path
import os,sys,json,subprocess,time,fcntl,hashlib,csv,shutil

def md5(p):
    h=hashlib.md5()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def run(config,start_phase='prepare_masks'):
    cfg=json.loads(Path(config).read_text());d=Path(cfg['stage']);code=Path(__file__).parent;dest=Path(cfg['destination'])
    assert cfg['rctd_mode']=='full' and cfg['reference_labels']=='major_cluster' and cfg['no_drive_upload']
    with (d/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);(d/'processing.pid').write_text(str(os.getpid()))
        dest.mkdir(exist_ok=True);(dest/'results').mkdir(exist_ok=True)
        phases=[('prepare_masks',['/usr/bin/python3',str(code/'zhang_all_masks_major_full_prepare.py'),config]),('prepare_counts',['/usr/bin/Rscript',str(code/'zhang_all_masks_major_full_counts.R'),config])]
        if start_phase not in ['prepare_masks','prepare_counts']:
            manifest=json.loads((d/'jobs.json').read_text());jobs=manifest['jobs']
        else:jobs=None
        active=start_phase=='prepare_masks'
        def execute(name,cmd):
            (d/'status').write_text(name);started=time.time();print('START',name,time.strftime('%Y-%m-%dT%H:%M:%S%z'),flush=True)
            with (d/f'{name}.log').open('x') as log:
                p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
                (d/'child.json').write_text(json.dumps({'phase':name,'pid':p.pid,'parent':os.getpid(),'argv':cmd,'started':started},indent=2));rc=p.wait()
            (d/f'{name}.exit.json').write_text(json.dumps({'exit_code':rc,'started':started,'ended':time.time(),'argv':cmd},indent=2))
            if rc:
                (d/'status').write_text('FAILED_'+name);raise RuntimeError(f'{name} failed with actual exit {rc}; successful independent scopes preserved')
            print('CLEAN_EXIT',name,flush=True)
        for name,cmd in phases:
            if name==start_phase:active=True
            if active:execute(name,cmd)
            else:assert json.loads((d/f'{name}.exit.json').read_text())['exit_code']==0
        jobs=json.loads((d/'jobs.json').read_text())['jobs']
        active=start_phase in ['prepare_masks','prepare_counts']
        for j in jobs:
            name='fit_'+j['job']
            if name==start_phase:active=True
            if active:execute(name,['/usr/bin/Rscript',str(code/'zhang_all_masks_major_full_fit.R'),config,j['job']])
            else:assert json.loads((d/f'{name}.exit.json').read_text())['exit_code']==0
            qa=json.loads((dest/'results'/j['job']/'qa.json').read_text());assert qa['readback_exact'] and md5(qa['path'])==qa['md5']
            cp=d/'counts'/f"{j['job']}.counts.rds"
            if cp.exists():cp.unlink() # input counts are now included verbatim in the verified final bundle
        rows=[]
        for j in jobs:
            qa=json.loads((dest/'results'/j['job']/'qa.json').read_text());assert qa['readback_exact'] and md5(qa['path'])==qa['md5']
            nq=qa.pop('numeric_qa');qa['method']=j['method'];qa['sample']=j['sample'];qa['cell_types']=';'.join(qa.get('cell_types') or []);qa.update(nq);rows.append(qa)
        fields=list(dict.fromkeys(k for row in rows for k in row))
        with (dest/'summary.csv').open('x',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
        provenance=Path(cfg['archive'])/'provenance/rctd';provenance.mkdir(exist_ok=False)
        shutil.copytree(code,provenance/'frozen_code')
        for name in ['source_revision.json','jobs.json','reference.qa.json','processing.log']:
            p=d/name
            if p.exists():shutil.copy2(p,provenance/name)
        for p in d.glob('*.exit.json'):shutil.copy2(p,provenance/p.name)
        for p in d.glob('*.log'):
            if p.name!='processing.log':shutil.copy2(p,provenance/p.name)
        shutil.copytree(d/'counts',provenance/'counts_QA')
        # These are task-created conversion intermediates, already represented
        # exactly in final bundles and backed by durable accepted source masks.
        for j in jobs:
            folder=d/'assignments'/j['job']
            assert folder.resolve().is_relative_to((d/'assignments').resolve())
            for n in ['labels.i32','regions.csv']:(folder/n).unlink()
            folder.rmdir()
        (d/'assignments').rmdir()
        (d/'reference.rds').unlink()
        complete={'status':'ALL_40_SCOPES_VERIFIED_PENDING_INTERPRETATION','jobs':len(rows),'full_fits':sum(r['status']=='NATIVE_FULL_FIT_READBACK_VERIFIED' for r in rows),'zero_region_scopes':sum(r['status']=='ZERO_REGIONS_NO_RCTD_FIT' for r in rows),'no_drive_upload':True,'all_native_containers_readback_and_md5_verified':True,'boundary_accuracy_not_measured':True}
        (dest/'completion.json').write_text(json.dumps(complete,indent=2));(d/'status').write_text(complete['status']);print(json.dumps(complete),flush=True)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('config');p.add_argument('--start-phase',default='prepare_masks');a=p.parse_args();run(str(Path(a.config).resolve()),a.start_phase)
