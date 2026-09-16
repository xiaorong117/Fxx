#!/usr/bin/env python3
"""Regression test for common 9-mm batch geometry and flow calculations."""

from __future__ import annotations

import importlib.util
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


source = Path(sys.argv[1]).resolve()
spec = importlib.util.spec_from_file_location(
    "prepare_batch35_inputs", source / "scripts/prepare_batch35_inputs.py")
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

domain_m = 0.009
porosity = 0.30
flow = module.flow_metrics(porosity * domain_m ** 3)
if not math.isclose(flow["domain_m"], domain_m, rel_tol=0.0, abs_tol=1e-16):
    raise SystemExit("batch domain is not 9 mm")
if not math.isclose(flow["q_peclet"], 1.0e-8, rel_tol=2e-15):
    raise SystemExit("target-Peclet flow is not 1e-8 m3/s")
if not math.isclose(flow["achieved_peclet"], module.TARGET_SAMPLE_PECLET,
                    rel_tol=2e-15):
    raise SystemExit("sample-scale Peclet calculation mismatch")

with tempfile.TemporaryDirectory(prefix="bubble_batch_scale_") as temporary:
    root = Path(temporary)
    raw_pore = root / "pore.txt"
    raw_throat = root / "throat.txt"
    out = root / "derived"
    pores = np.array([
        [1.0, 2.0, 3.0, 4.0, 5.0, 0.03, 6.0, 0.0, 0.4],
        [7.0, 8.0, 9.0, 10.0, 11.0, 0.04, 12.0, 0.0, 0.6],
    ])
    throats = np.array([
        [1.0, 2.0, 2.0, 3.0, 0.05, 4.0, 0.0, 5.0, 6.0, 0.5],
        [2.0, 1.0, 2.0, 3.0, 0.05, 0.0, 0.0, 5.0, 6.0, 0.5],
    ])
    np.savetxt(raw_pore, pores)
    np.savetxt(raw_throat, throats)
    subprocess.run([
        sys.executable, str(source / "scripts/convert_batch_pnm_full_geometry.py"),
        "--raw-pore", str(raw_pore), "--raw-throat", str(raw_throat),
        "--output-directory", str(out), "--geometry-scale", "2.0",
        "--normalize-zero-volume-saturation",
    ], check=True, stdout=subprocess.DEVNULL)
    converted_pores = np.loadtxt(out / "PNM_new_pore_base.txt")
    converted_edges = np.loadtxt(out / "PNM_new_throat_base.txt")
    if not np.allclose(converted_pores[0, 3:7], [2.0, 4.0, 6.0, 8.0]):
        raise SystemExit("pore coordinates/radius were not length-scaled")
    if not math.isclose(converted_pores[0, 7], 20.0):
        raise SystemExit("pore area was not squared-scaled")
    if not math.isclose(converted_pores[0, 9], 48.0):
        raise SystemExit("pore volume was not cubed-scaled")
    if not np.allclose(converted_edges[:, 5], [10.0, 12.0, 10.0, 12.0]):
        raise SystemExit("half-edge lengths were not scaled")
    expected_throat_geometry = np.array([
        [4.0, 12.0, 0.05, 32.0], [4.0, 12.0, 0.05, 32.0],
        [4.0, 12.0, 0.05, 0.0], [4.0, 12.0, 0.05, 0.0],
    ])
    if not np.allclose(converted_edges[:, 7:11], expected_throat_geometry):
        raise SystemExit("throat radius/area/shape/volume scaling mismatch")
    audit = json.loads((out / "conversion_audit.json").read_text())
    normalization = audit["geometry_normalization"]
    if (normalization["isotropic_length_scale"] != 2.0 or
            normalization["area_scale"] != 4.0 or
            normalization["volume_scale"] != 8.0):
        raise SystemExit("conversion scaling provenance is incomplete")
    zero_measure = audit["zero_measure_saturation_normalization"]
    if zero_measure["source_throat_ids"] != [2]:
        raise SystemExit("zero-measure saturation correction was not audited")
    if not np.allclose(converted_edges[2:, 6], 1.0):
        raise SystemExit("zero-volume throat remained incorrectly marked with gas")

print("PASS common 9-mm batch flow and dimensional geometry normalization")
