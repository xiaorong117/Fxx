#!/usr/bin/env python3
"""Regression for one-component derived-PNM exclusion and provenance."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


project = Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory(prefix="bubble_boundaryless_filter_") as temporary:
    root = Path(temporary); pore = root / "pore.txt"; throat = root / "throat.txt"
    config = root / "config.json"; output = root / "out"
    pores = np.array([
        [1,0,1,0,0,0,2,3,.03,10,0,.5],
        [2,0,2,1,0,0,2,3,.03,10,0,1],
        [3,0,3,10,0,0,2,3,.03,10,0,.5],
        [4,2,3,11,0,0,2,3,.03,10,1,1],
    ],float)
    edges = np.array([
        [1,1,2,1,1,1,.5,1,1,.03,1],
        [2,3,4,0,1,1,1,1,1,.03,0],
    ],float)
    np.savetxt(pore,pores);np.savetxt(throat,edges)
    config.write_text(json.dumps({
        "geometry":{"input_length_to_m":1e-6,"input_volume_to_m3":1e-18},
        "physics":{"sigma_N_m":.07197,"contact_angle_rad":0,"qin_b":6.83,
                   "Z":1,"T_K":298.15},
        "pressure_reference":{"background_abs_Pa":101325},
        "boundary":{"Pl_in_rel_Pa":100,"Pl_out_rel_Pa":0},
        "initial":{"C_initial_mol_m3":.001},
        "active_set":{"Sw_min":1e-6,"Sg_off":1e-8,"ng_off_mol":1e-30},
    }))
    before=(pore.read_bytes(),throat.read_bytes())
    subprocess.run([sys.executable,str(project/'scripts/exclude_boundaryless_component.py'),
                    '--pore',str(pore),'--throat',str(throat),'--config',str(config),
                    '--component-node-id','1','--output-directory',str(output)],check=True,
                   stdout=subprocess.DEVNULL)
    audit=json.loads((output/'boundaryless_component_freeze_audit.json').read_text())
    assert audit['excluded_old_node_ids']==[1,2]
    assert audit['excluded_boundary_flag_set']==[0]
    assert audit['dynamic_node_count']==2 and audit['dynamic_edge_count']==1
    derived=np.loadtxt(output/'PNM_new_pore_virtual_dynamic.txt')
    assert derived[:,0].tolist()==[1,2]
    assert int(derived[1,2])==1  # virtual source_id remapped with its real pore
    assert (pore.read_bytes(),throat.read_bytes())==before
    bad=root/'bad'
    failed=subprocess.run([sys.executable,str(project/'scripts/exclude_boundaryless_component.py'),
                           '--pore',str(pore),'--throat',str(throat),'--config',str(config),
                           '--component-node-id','3','--output-directory',str(bad)],
                          stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    assert failed.returncode!=0
print('PASS case-specific boundaryless component exclusion')
