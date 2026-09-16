#!/usr/bin/env python3
"""Read back the small relative-pressure visualization delivery."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

output=Path(sys.argv[1]).resolve()
script=Path(sys.argv[2]).resolve()
spec=importlib.util.spec_from_file_location('relative_visualize',script)
rv=importlib.util.module_from_spec(spec);spec.loader.exec_module(rv)

audit=json.loads((output/'visualization_audit.json').read_text())
if not audit.get('xml_readback_passed') or audit.get('local_edge_count',0)<1:
    raise SystemExit('visualization XML/local-region audit did not pass')
if audit.get('paraview_render_validated') not in (True,False):
    raise SystemExit('ParaView validation state is not explicit')

required_points={'node_id','source_type','boundary_flag','is_virtual_boundary','component_id',
    'Pl_relative_Pa','Pg_relative_Pa','Pl_absolute_Pa','Pg_absolute_Pa','C','Sw','Sg',
    'active_gas','gas_moles','gas_volume','bubble_equivalent_radius'}
required_cells={'liquid_flux_m3_s','advective_component_flux_mol_s',
                'diffusive_component_flux_mol_s','is_boundary_edge'}
entries=[]
for dataset in ET.parse(output/'states.pvd').getroot().find('Collection'):
    time=float(dataset.attrib['timestep']);path=output/dataset.attrib['file'];data=rv.read(path)
    points={data.GetPointData().GetArrayName(i) for i in range(data.GetPointData().GetNumberOfArrays())}
    cells={data.GetCellData().GetArrayName(i) for i in range(data.GetCellData().GetNumberOfArrays())}
    if not required_points<=points or not required_cells<=cells:
        raise SystemExit(f'missing required VTP fields at {time}')
    entries.append((time,path))
if not entries or entries[0][0]!=0 or entries[-1][0]<=0:
    raise SystemExit('PVD physical times are invalid')

for name,flag in [('inlet_boundary.vtp',1),('outlet_boundary.vtp',2)]:
    data=rv.read(output/name);pd=data.GetPointData()
    for field in ('node_id','boundary_flag','radius','connected_real_node_id'):
        if pd.GetArray(field) is None: raise SystemExit(f'{name} lacks {field}')
    if any(int(pd.GetArray('boundary_flag').GetTuple1(i))!=flag for i in range(data.GetNumberOfPoints())):
        raise SystemExit(f'{name} has wrong boundary flag')
connections=rv.read(output/'boundary_connections.vtp')
if connections.GetNumberOfCells()!=2 or connections.GetNumberOfPoints()!=4:
    raise SystemExit('boundary connection geometry readback mismatch')

local=rv.read(output/('local_'+entries[0][1].name))
ids={int(local.GetPointData().GetArray('node_id').GetTuple1(i)) for i in range(local.GetNumberOfPoints())}
if ids!=set(audit['local_node_ids']) or local.GetNumberOfCells()!=audit['local_edge_count']:
    raise SystemExit('local remapping left unused points or lost connectivity')

# Explicitly exercise the empty-bubble fallback and empty-local failure path.
data=rv.read(entries[0][1]);active=data.GetPointData().GetArray('active_gas')
for i in range(active.GetNumberOfTuples()): active.SetTuple1(i,0)
import numpy as np
virtual=np.array([bool(data.GetPointData().GetArray('is_virtual_boundary').GetTuple1(i))
                  for i in range(data.GetNumberOfPoints())])
component=np.array([int(data.GetPointData().GetArray('component_id').GetTuple1(i))
                    for i in range(data.GetNumberOfPoints())])
if rv.choose_local_region(data,virtual,component,audit['main_flow_component'])[1] != \
        'fallback_main_flow_real_node_no_active_bubble':
    raise SystemExit('empty-bubble fallback was not selected')
try: rv.line_subset(data,[])
except RuntimeError as error:
    if 'no internal edge' not in str(error): raise
else: raise SystemExit('empty local edge set was not handled explicitly')

rows=list(csv.DictReader((output/'diagnostics.tsv').open(),delimiter='\t'))
if not rows or not {'free_mol_total','free_mol_main_flow_component',
                    'gas_volume_total','gas_volume_main_flow_component'}<=set(rows[0]):
    raise SystemExit('full-network/main-component histories are not distinct')

manifest=json.loads((output/'run_manifest.json').read_text())
for name in ('states.pvd','local_states.pvd','inlet_boundary.vtp','outlet_boundary.vtp',
             'boundary_connections.vtp','visualization_audit.json','history.png','load_paraview.py'):
    path=(output/name).resolve();expected=hashlib.sha256(path.read_bytes()).hexdigest()
    if manifest['outputs'].get(str(path))!=expected:
        raise SystemExit(f'postprocessing manifest hash mismatch: {name}')
print(f'PASS relative visualization readback: {len(entries)} frames; ParaView rendered={audit["paraview_render_validated"]}')
