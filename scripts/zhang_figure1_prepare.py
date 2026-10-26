"""Read the accepted display mask; no segmentation or RCTD is recomputed."""
from pathlib import Path
import sys,json,importlib.util
import numpy as np
import tifffile
from shapely.geometry import mapping

archive=Path(sys.argv[1]);work=Path(sys.argv[2])
method='Zhang_revised_s3p5_expand5';sample='CD2'
spec=importlib.util.spec_from_file_location('lossless',Path(__file__).with_name('polygonize_zhang_selected.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
source=archive/'segmentations'/method/f'{sample}.labels.tiff'
metrics=json.loads((source.parent/f'{sample}.metrics.json').read_text())
assert module.md5(source)==metrics['label_md5']
mask=tifffile.imread(source)
roi=json.loads((archive/'provenance/visual_fields.json').read_text())[sample][0]
x,y,size=roi['x'],roi['y'],roi['size']
labels=np.unique(mask[y:y+size,x:x+size]);labels=labels[labels>0]
selected=np.where(np.isin(mask,labels),mask,0)
polygons=module.polygonize(selected)
assert set(polygons)==set(map(int,labels))
ids,area=np.unique(selected[selected>0],return_counts=True)
features=[]
for label,n in zip(ids,area):
    geometry=polygons[int(label)]
    assert geometry.area==int(n)
    features.append({'type':'Feature','properties':{'label':int(label),'area_hires_px2':int(n)},'geometry':mapping(geometry)})
(work/'field.polygons.geojson').write_text(json.dumps({'type':'FeatureCollection','features':features},separators=(',',':')))
(work/'field_geometry.qa.json').write_text(json.dumps({'method':method,'sample':sample,'roi':roi,'mask_md5':metrics['label_md5'],
    'regions':len(features),'exact_pixel_area':True,'all_parts_and_holes_preserved':True,'no_simplification':True,'no_resegmentation':True},indent=2))
print('EXACT_ACCEPTED_D_FIELD_READY',len(features),flush=True)
