#!/usr/bin/env python3
"""Audit the deterministic virtual-boundary CLI run."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path


output = Path(sys.argv[1])
with (output / "diagnostics.tsv").open(newline="") as stream:
    rows = list(csv.DictReader(stream, delimiter="\t"))
if len(rows) != 5:
    raise SystemExit("virtual CLI did not accept five steps")
for row in rows:
    if (float(row["residual_scaled_inf"]) > 1e-8 or
            float(row["relative_update_norm"]) > 1e-9 or
            float(row["epsilon_N_adjusted"]) > 1e-10 or
            float(row["epsilon_V_adjusted"]) > 1e-10):
        raise SystemExit("virtual CLI convergence/conservation audit failed")

manifest = json.loads((output / "run_manifest.json").read_text())
network = manifest["network"]
if (network["control_volumes"] != 5 or network["half_edges"] != 4 or
        network["virtual_boundary_nodes"] != 2 or network["boundary_edges"] != 2 or
        network["base_node_count"] != 3 or network["base_edge_count"] != 2):
    raise SystemExit("virtual CLI network manifest mismatch")
if manifest["virtual_boundary"]["inlet_virtual_nodes"] != 1 or \
        manifest["virtual_boundary"]["outlet_virtual_nodes"] != 1:
    raise SystemExit("virtual CLI inlet/outlet count mismatch")
if set(manifest.get("input_provenance", {})) != {
        "raw_pore_file", "raw_throat_file", "base_pore_file",
        "base_throat_file", "base_connect_file", "input_geometry_audit_file",
        "virtual_boundary_geometry_audit_file"}:
    raise SystemExit("virtual CLI provenance hashes are incomplete")

checkpoint = json.loads((output / "checkpoint.json").read_text())
if len(checkpoint["state"]) != 5 or len(checkpoint["state_extended_precision"]) != 5:
    raise SystemExit("virtual source_type=2 state was not checkpointed")
for node in (3, 4):
    state = checkpoint["state"][node]
    if state[2] != 1.0 or state[4]:
        raise SystemExit("virtual reservoir acquired storage/gas state")

vtk = (output / "state_000005.vtk").read_text()
for field in ("node_id", "source_type", "boundary_flag", "is_virtual_boundary",
              "active_gas", "Pl", "Pg", "Sw", "C", "is_boundary_edge",
              "edge_length", "edge_length_voxels", "source_old_throat_id"):
    if f"SCALARS {field} " not in vtk:
        raise SystemExit(f"virtual CLI VTK lacks {field}")
for name in ("run.log", "input_geometry_unit_audit.json",
             "virtual_boundary_geometry_audit.json", "degenerate_geometry_audit.json"):
    if not (output / name).is_file():
        raise SystemExit(f"virtual CLI output missing {name}")
print("PASS virtual-boundary CLI flux/conservation/checkpoint/VTK audit")
