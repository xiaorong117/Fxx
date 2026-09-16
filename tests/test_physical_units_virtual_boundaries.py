#!/usr/bin/env python3
"""Audit physical units and one-to-one virtual-boundary generation."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


source = Path(sys.argv[1]).resolve()
root = Path(sys.argv[2]).resolve()
unit_script = source / "scripts/audit_pnm_physical_units.py"
virtual_script = source / "scripts/generate_virtual_boundary_network.py"
pnextract_source = Path("/workspace/pnflow/src/pnm/pnextract/blockNet_write_cnm.cpp")
old_pore = Path("/workspace/SinglePhase/Project/kong/fina_test/reversed_full_pore.txt")
old_throat = Path("/workspace/SinglePhase/Project/kong/fina_test/reversed_full_throat.txt")

raw_names = [
    "PNM_pb.txt", "PNM_pt.txt", "PNM_new_pore.txt",
    "PNM_new_throat.txt", "PNM_new_pore_connect.txt",
    "convert_throats_to_pores.py", "trans-react.cpp",
]
before = {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in raw_names}
expected_baseline = {
    "PNM_pb.txt": "72e32ad727f31912529f2b21a9c90bcdf4d9079a22ed07288fca891974f18282",
    "PNM_pt.txt": "24ac31276f29073fa8f2c0cfab4d6ea7ffd35eb4a79b7d0d9bbf2e950ca0055a",
    "PNM_new_pore.txt": "08d1bfc68ffdf678aeac3ebfc42c5dcb4a725213139c6f5316a9923c2cd9c5d8",
    "PNM_new_throat.txt": "c6a5fa83b2993715159a5afc35c37d1e1a530b368273bedb9b13fdb1e3f1c53a",
    "PNM_new_pore_connect.txt": "d819271381d01217661679c52db3673482cd2b7ed7b880eb7b55dc0384de0561",
    "convert_throats_to_pores.py": "ed0df60db78db0f67b94183fcf1eacbae6b27f89a9bc603b82a3e9c8367e5547",
    "trans-react.cpp": "d70e9c44eb54e51418791a24d587ef3e12b75314bafcb6e1e3fcd0e92269444a",
}
if before != expected_baseline:
    raise SystemExit("raw/reference PNM or source baseline hash mismatch")

with tempfile.TemporaryDirectory(prefix="bubble_virtual_boundary_") as temporary:
    out = Path(temporary)
    unit_audit = out / "unit.json"
    subprocess.run([
        sys.executable, str(unit_script), "--root", str(root),
        "--pnextract-source", str(pnextract_source), "--output", str(unit_audit),
    ], check=True)
    unit = json.loads(unit_audit.read_text())
    if (unit["image_dimensions_voxels"] != [1200, 1200, 1200] or
            unit["physical_domain_um"] != [9000.0, 9000.0, 9000.0] or
            not math.isclose(1200 * unit["voxel_size_um"], 9000.0) or
            unit["coordinate_min_um"] != [3.75, 3.75, 3.75] or
            unit["coordinate_max_um"] != [8996.3, 8996.3, 8996.3] or
            unit["volume_quantum_um3"] != 421.875):
        raise SystemExit("physical-image/PNM unit evidence mismatch")
    scales = unit["si_conversions"]
    if (scales["input_length_to_m"] != 1e-6 or
            scales["input_area_to_m2"] != scales["input_length_to_m"] ** 2 or
            not math.isclose(scales["input_volume_to_m3"],
                             scales["input_length_to_m"] ** 3, rel_tol=1e-15)):
        raise SystemExit("SI length/area/volume conversion mismatch")

    pore = out / "pore.txt"
    throat = out / "throat.txt"
    connect = out / "connect.txt"
    audit_path = out / "virtual.json"
    subprocess.run([
        sys.executable, str(virtual_script),
        "--base-pore", str(root / "PNM_new_pore.txt"),
        "--base-throat", str(root / "PNM_new_throat.txt"),
        "--base-connect", str(root / "PNM_new_pore_connect.txt"),
        "--output-pore", str(pore), "--output-throat", str(throat),
        "--output-connect", str(connect), "--audit-output", str(audit_path),
        "--axis", "x", "--extension-voxels", "1", "--voxel-size-um", "7.5",
        "--old-model-pore", str(old_pore), "--old-model-throat", str(old_throat),
    ], check=True)
    base_pores = np.loadtxt(root / "PNM_new_pore.txt")
    derived_pores = np.loadtxt(pore)
    base_edges = np.loadtxt(root / "PNM_new_throat.txt")
    derived_edges = np.loadtxt(throat)
    audit = json.loads(audit_path.read_text())
    virtual = derived_pores[:, 1] == 2
    real = ~virtual
    flags = derived_pores[:, 10].astype(int)
    if (np.count_nonzero(virtual & (flags == 1)) != 724 or
            np.count_nonzero(virtual & (flags == 2)) != 785 or
            np.any(flags[real] != 0) or np.any(~np.isin(flags[virtual], (1, 2)))):
        raise SystemExit("boundary flags were not transferred one-to-one")
    degree = np.zeros(len(derived_pores), dtype=int)
    for a, b in derived_edges[:, 1:3].astype(int):
        degree[a - 1] += 1
        degree[b - 1] += 1
    if np.any(degree[virtual] != 1):
        raise SystemExit("virtual node coordination is not exactly one")
    boundary_edges = derived_edges[len(base_edges):]
    if (len(boundary_edges) != np.count_nonzero(virtual) or
            np.any(boundary_edges[:, 5] != 7.5) or
            np.any(boundary_edges[:, 3] != 0) or
            audit["new_zero_length_boundary_edge_count"] != 0 or
            audit["virtual_nodes_in_nontrivial_zero_distance_cluster"] != 0):
        raise SystemExit("generated boundary-edge topology/length is invalid")
    if (not math.isclose(audit["connected_real_virtual_clearance_min_um"], 7.5,
                         abs_tol=2e-12) or
            not math.isclose(audit["connected_real_virtual_clearance_max_um"], 7.5,
                             abs_tol=2e-12)):
        raise SystemExit("connected virtual-real clearance is not 7.5 um")
    source = derived_pores[virtual, 2].astype(int) - 1
    direction = np.sign(derived_pores[virtual, 3] - base_pores[source, 3])
    if np.any(direction[flags[virtual] == 1] != -1) or np.any(direction[flags[virtual] == 2] != 1):
        raise SystemExit("virtual inlet/outlet extension direction is incorrect")
    if audit["isolated_former_boundary_count"] != 35:
        raise SystemExit("isolated former boundary inventory mismatch")

after = {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in raw_names}
if after != before:
    raise SystemExit("raw/reference PNM or source file changed during generation")
print("PASS physical um units and 1509 one-to-one virtual boundaries")
