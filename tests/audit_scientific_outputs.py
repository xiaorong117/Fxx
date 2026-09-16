#!/usr/bin/env python3
import csv,json,math,sys,xml.etree.ElementTree as ET
from pathlib import Path

root=Path(sys.argv[1]);segments=int(sys.argv[2])
def rows(name):
    with (root/name).open(newline='') as stream:return list(csv.DictReader(stream,delimiter='\t'))
manifest=json.loads((root/'run_manifest.json').read_text());steps=int(manifest['accepted_steps'])
g=rows('global_history.tsv');s=rows('segment_history.tsv');b=rows('balance_history.tsv');geom=rows('segment_geometry.tsv')
if len(g)!=steps or len(b)!=steps or len(s)!=steps*segments or len(geom)!=segments:
    raise SystemExit('scientific row counts mismatch')
if not (root/'BALANCE_README.md').is_file():raise SystemExit('missing BALANCE_README')
for index,row in enumerate(g,1):
    part=[x for x in s if int(x['step'])==index]
    for local,global_key in [('segment_gas_volume_m3','gas_volume_m3'),('segment_free_gas_mol','free_gas_mol'),
                             ('segment_dissolved_gas_mol','dissolved_gas_mol'),
                             ('segment_cumulative_interphase_transfer_mol','cumulative_interphase_transfer_mol')]:
        if not math.isclose(sum(float(x[local]) for x in part),float(row[global_key]),rel_tol=1e-11,abs_tol=1e-28):
            raise SystemExit(f'segment/global invariant failed: {local}')
    if int(b[index-1]['component_balance_pass'])!=1 or int(b[index-1]['liquid_balance_pass'])!=1:
        raise SystemExit('balance PASS failed')
if manifest['linear_backend']['selected']!='amgx' or not manifest['linear_backend']['amgx_compiled']:
    raise SystemExit('manifest does not prove AMGX execution')
if manifest['termination']['actual_reason']!='right_censored':raise SystemExit('maximum-time case was not right_censored')
cp=json.loads((root/'checkpoint.json').read_text())
if cp['format']!='bubble_checkpoint_v6' or len(cp['cumulative_interphase_transfer_by_node_mol'])!=len(cp['state']):
    raise SystemExit('checkpoint v6 cumulative state missing')
entries=list(ET.parse(root/'states.pvd').getroot().find('Collection'))
if not entries or any(not x.attrib['file'].endswith('.vtp') for x in entries):
    raise SystemExit('production-format PVD does not reference VTP')
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'relative-vis-deps'))
import vtk
reader=vtk.vtkXMLPolyDataReader();reader.SetFileName(str(root/entries[-1].attrib['file']));reader.Update()
data=reader.GetOutput()
if reader.GetErrorCode() or data.GetNumberOfPoints()!=5 or data.GetNumberOfCells()!=4:
    raise SystemExit('compressed VTP readback topology failed')
for field in ('normalized_x','segment_id','local_Sg','initial_Sg','free_gas_moles_lost',
              'cumulative_interphase_transfer','adjacent_constriction_ratio_p50'):
    if data.GetPointData().GetArray(field) is None:raise SystemExit(f'compressed VTP lacks {field}')
print('PASS AMGX scientific histories, invariants, right censoring, checkpoint continuity state')
