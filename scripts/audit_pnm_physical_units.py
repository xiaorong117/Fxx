#!/usr/bin/env python3
"""Machine-readable audit that the supplied pnextract geometry is in micrometres."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pnextract-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--image-dimensions", type=int, nargs=3,
                        default=[1200, 1200, 1200])
    parser.add_argument("--voxel-size-um", type=float, default=7.5)
    args = parser.parse_args()

    root = args.root.resolve()
    source = args.pnextract_source.resolve()
    pore_path = root / "PNM_pb.txt"
    throat_path = root / "PNM_pt.txt"
    pores = np.loadtxt(pore_path)
    throats = np.loadtxt(throat_path)

    image_dimensions = list(args.image_dimensions)
    voxel_size_um = args.voxel_size_um
    require(all(value > 0 for value in image_dimensions) and voxel_size_um > 0.0,
            "image dimensions and voxel size must be positive")
    physical_domain_um = [value * voxel_size_um for value in image_dimensions]
    coordinate_min_um = [float(pores[:, axis].min()) for axis in range(3)]
    coordinate_max_um = [float(pores[:, axis].max()) for axis in range(3)]
    require(all(value >= 0.0 for value in coordinate_min_um), "negative PNM coordinate")
    require(all(value <= domain + 1e-9 for value, domain in
                zip(coordinate_max_um, physical_domain_um)),
            "PNM coordinate exceeds the configured image domain")
    require(all(value >= 0.99 * domain for value, domain in
                zip(coordinate_max_um, physical_domain_um)),
            "PNM coordinate range does not span the known image domain")

    quantum_um3 = voxel_size_um ** 3
    volume = np.concatenate((pores[:, 6], throats[:, 5]))
    nearest = np.rint(volume / quantum_um3) * quantum_um3
    relative_quantization_error = np.abs(volume - nearest) / np.maximum(
        np.abs(nearest), quantum_um3
    )
    maximum_quantization_error = float(relative_quantization_error.max())
    require(maximum_quantization_error <= 5.0e-5,
            "volume values are inconsistent with 7.5^3 um3 quantization")

    text = source.read_text(encoding="utf-8", errors="replace")
    evidence = {
        "voxel_size_binding": "const double dx = cg.vxlSize;",
        "throat_radius": "tr.radius()*dx",
        "throat_length": "t_lengthP1toP2s[ti]*dx",
        "throat_volume": "tr.volumn*dx*dx*dx",
        "pore_coordinates": "por.mb->fi*dx",
        "pore_volume": "por.volumn*dx*dx*dx",
        "pore_radius": "por.radius()*dx",
    }
    for label, snippet in evidence.items():
        require(snippet in text, f"pnextract source evidence missing: {label}")

    result = {
        "format": "bubble_pnm_physical_unit_audit_v1",
        "conclusion": "consistent_with_physical_um_not_voxel_indices_or_nm",
        "input_basis": "physical",
        "input_length_unit": "um",
        "image_dimensions_voxels": image_dimensions,
        "voxel_size_um": voxel_size_um,
        "physical_domain_um": physical_domain_um,
        "coordinate_min_um": coordinate_min_um,
        "coordinate_max_um": coordinate_max_um,
        "volume_quantum_um3": quantum_um3,
        "volume_quantization_max_relative_error": maximum_quantization_error,
        "volume_quantization_rounding_limit": 5.0e-5,
        "si_conversions": {
            "input_length_to_m": 1.0e-6,
            "input_area_to_m2": 1.0e-12,
            "input_volume_to_m3": 1.0e-18,
            "area_equals_length_squared": True,
            "volume_equals_length_cubed": True,
        },
        "inputs": {
            str(pore_path): sha256(pore_path),
            str(throat_path): sha256(throat_path),
        },
        "pnextract_source": {
            "path": str(source),
            "sha256": sha256(source),
            "matched_write_expressions": evidence,
        },
    }
    require(math.isclose(result["si_conversions"]["input_area_to_m2"],
                         result["si_conversions"]["input_length_to_m"] ** 2,
                         rel_tol=1e-15), "area scale mismatch")
    require(math.isclose(result["si_conversions"]["input_volume_to_m3"],
                         result["si_conversions"]["input_length_to_m"] ** 3,
                         rel_tol=1e-15), "volume scale mismatch")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "conclusion": result["conclusion"],
        "coordinate_min_um": coordinate_min_um,
        "coordinate_max_um": coordinate_max_um,
        "volume_quantization_max_relative_error": maximum_quantization_error,
    }))


if __name__ == "__main__":
    main()
