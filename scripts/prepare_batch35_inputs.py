#!/usr/bin/env python3
"""Prepare immutable, audited solver inputs for all PNM-0.6/0.8 gradations."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np


IMAGE_DIMENSIONS = (1200, 1200, 1200)
VOXEL_SIZE_UM = 7.5
TARGET_DOMAIN_UM = IMAGE_DIMENSIONS[0] * VOXEL_SIZE_UM
TARGET_SAMPLE_PECLET = 555.5555555555555
DIFFUSIVITY_M2_S = 2.0e-9
VISCOSITY_PA_S = 8.9e-4
DENSITY_KG_M3 = 997.0
SURFACE_TENSION_N_M = 0.07197
MAX_CAPILLARY_NUMBER = 1.0e-5


def source_domain_um(directory: Path) -> float:
    """Read the physical cube size recorded by pnextract itself."""
    header = (directory / "default_pore.out").open(
        encoding="utf-8", errors="replace").readline().split()
    if len(header) < 4:
        raise RuntimeError(f"{directory}: malformed default_pore.out header")
    dimensions_um = [float(value) * 1.0e6 for value in header[1:4]]
    if max(dimensions_um) - min(dimensions_um) > 1.0e-6 * max(dimensions_um):
        raise RuntimeError(f"{directory}: non-cubic pnextract domain {dimensions_um}")
    return float(sum(dimensions_um) / 3.0)


def flow_metrics(control_volume_m3: float) -> dict[str, float]:
    domain_m = TARGET_DOMAIN_UM * 1.0e-6
    area_m2 = domain_m * domain_m
    porosity = control_volume_m3 / domain_m ** 3
    if not 0.0 < porosity < 1.0:
        raise RuntimeError(f"invalid control-volume porosity {porosity}")
    q_peclet = TARGET_SAMPLE_PECLET * DIFFUSIVITY_M2_S * area_m2 / domain_m
    q_capillary = (MAX_CAPILLARY_NUMBER * porosity * area_m2 *
                   SURFACE_TENSION_N_M / VISCOSITY_PA_S)
    recommended_q = min(q_peclet, q_capillary)
    achieved_peclet = (recommended_q * domain_m /
                        (area_m2 * DIFFUSIVITY_M2_S))
    interstitial_velocity = recommended_q / (porosity * area_m2)
    return {
        "domain_m": domain_m,
        "area_m2": area_m2,
        "porosity": porosity,
        "q_peclet": q_peclet,
        "q_capillary": q_capillary,
        "recommended_q": recommended_q,
        "achieved_peclet": achieved_peclet,
        "interstitial_velocity": interstitial_velocity,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    temporary.replace(path)


def slug(text: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_.-")
    return value or "unnamed"


def command(log: Path, arguments: list[str]) -> None:
    with log.open("w") as stream:
        completed = subprocess.run(arguments, stdout=stream,
                                   stderr=subprocess.STDOUT, check=False)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}); see {log}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    input_root = args.input_root.resolve()
    output_root = args.output_root.resolve()
    project = args.project.resolve()
    if output_root.exists():
        raise SystemExit(f"refusing to overwrite {output_root}")
    output_root.mkdir(parents=True)

    complete: list[tuple[str, str, Path]] = []
    parents: set[str] = set()
    for pore in input_root.rglob("PNM_pb.txt"):
        directory = pore.parent
        if directory.name not in ("PNM-0.6", "PNM-0.8"):
            continue
        parent = directory.parent.name
        parents.add(parent)
        if (directory / "PNM_pt.txt").is_file():
            complete.append((parent, directory.name, directory))
    complete.sort(key=lambda item: (item[0], item[1]))
    missing = [{"parent_case": parent, "variant": variant}
               for parent in sorted(parents)
               for variant in ("PNM-0.6", "PNM-0.8")
               if not any(item[0] == parent and item[1] == variant
                          for item in complete)]
    variant_counts = {
        variant: sum(item[1] == variant for item in complete)
        for variant in ("PNM-0.6", "PNM-0.8")
    }
    if (len(parents) != 35 or len(complete) != 70 or missing or
            variant_counts != {"PNM-0.6": 35, "PNM-0.8": 35}):
        raise SystemExit(
            "expected 35 parents and complete 35+35 variants, got "
            f"{len(parents)} parents, {len(complete)} complete, "
            f"counts={variant_counts}, missing={missing}")

    results = []
    try:
        for index, (parent, variant, source) in enumerate(complete, 1):
            variant_slug = "p06" if variant == "PNM-0.6" else "p08"
            case_slug = slug(parent)
            destination = output_root / case_slug / variant_slug
            atomic_json(output_root / "progress.json", {
                "stage": "prepare_batch35_inputs", "status": "running",
                "completed": index - 1, "total": len(complete),
                "case": parent, "variant": variant,
                "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            })
            raw_pore = source / "PNM_pb.txt"
            raw_throat = source / "PNM_pt.txt"
            pores = np.loadtxt(raw_pore)
            throats = np.loadtxt(raw_throat)
            coordinate_max = pores[:, :3].max(axis=0)
            source_domain = source_domain_um(source)
            if not np.all(coordinate_max <= source_domain + 1e-6) or not np.all(
                    coordinate_max >= 0.99 * source_domain):
                raise RuntimeError(
                    f"{source}: coordinate range inconsistent with source domain "
                    f"{source_domain} um")
            geometry_scale = TARGET_DOMAIN_UM / source_domain
            control_volume_m3 = float(
                pores[:, 6].sum() + throats[:, 5].sum()) * geometry_scale ** 3 * 1e-18
            flow = flow_metrics(control_volume_m3)
            domain_m = flow["domain_m"]
            area_m2 = flow["area_m2"]
            porosity = flow["porosity"]
            q_peclet = flow["q_peclet"]
            q_capillary = flow["q_capillary"]
            recommended_q = flow["recommended_q"]
            achieved_peclet = flow["achieved_peclet"]
            interstitial_velocity = flow["interstitial_velocity"]
            throat_radius_p50_m = (float(np.median(throats[:, 2])) *
                                   geometry_scale * 1e-6)
            pore_reynolds = (DENSITY_KG_M3 * interstitial_velocity *
                             2.0 * throat_radius_p50_m / VISCOSITY_PA_S)
            capillary_number = (VISCOSITY_PA_S * interstitial_velocity /
                                SURFACE_TENSION_N_M)

            command(output_root / f"{case_slug}_{variant_slug}_convert.log", [
                sys.executable, str(project / "scripts/convert_batch_pnm_full_geometry.py"),
                "--raw-pore", str(raw_pore), "--raw-throat", str(raw_throat),
                "--output-directory", str(destination),
                "--geometry-scale", str(geometry_scale),
                "--normalize-zero-volume-saturation"])
            command(destination / "virtual_boundary.log", [
                sys.executable, str(project / "scripts/generate_virtual_boundary_network.py"),
                "--base-pore", str(destination / "PNM_new_pore_base.txt"),
                "--base-throat", str(destination / "PNM_new_throat_base.txt"),
                "--base-connect", str(destination / "PNM_new_connect_base.txt"),
                "--output-pore", str(destination / "PNM_new_pore_virtual.txt"),
                "--output-throat", str(destination / "PNM_new_throat_virtual.txt"),
                "--output-connect", str(destination / "PNM_new_connect_virtual.txt"),
                "--audit-output", str(destination / "virtual_boundary_geometry_audit.json"),
                "--axis", "x", "--extension-voxels", "1",
                "--voxel-size-um", str(VOXEL_SIZE_UM),
                "--old-model-pore", "/workspace/SinglePhase/Project/kong/fina_test/reversed_full_pore.txt",
                "--old-model-throat", "/workspace/SinglePhase/Project/kong/fina_test/reversed_full_throat.txt"])
            command(destination / "strict_geometry.log", [
                sys.executable, str(project / "scripts/audit_strict_geometry.py"),
                str(destination), "--pore-file", "PNM_new_pore_virtual.txt",
                "--throat-file", "PNM_new_throat_virtual.txt", "--report",
                str(destination / "strict_geometry_audit.tsv")])
            conversion = json.loads((destination / "conversion_audit.json").read_text())
            virtual = json.loads((destination / "virtual_boundary_geometry_audit.json").read_text())
            record = {
                "case_index": index,
                "parent_case": parent,
                "case_slug": case_slug,
                "variant": variant,
                "variant_slug": variant_slug,
                "source_directory": str(source),
                "derived_directory": str(destination),
                "source_pnextract_domain_um": [source_domain] * 3,
                "geometry_length_scale_to_effective_domain": geometry_scale,
                "nominal_domain_um": [TARGET_DOMAIN_UM] * 3,
                "image_dimensions_voxels": list(IMAGE_DIMENSIONS),
                "voxel_size_um": VOXEL_SIZE_UM,
                "geometry_basis": "user_confirmed_1200_cubed_at_7.5_um",
                "source_metadata_conflict_preserved": not np.isclose(
                    source_domain, TARGET_DOMAIN_UM),
                "zero_measure_saturation_normalization": conversion[
                    "zero_measure_saturation_normalization"],
                "control_volume_m3": control_volume_m3,
                "control_volume_porosity": porosity,
                "target_sample_peclet": TARGET_SAMPLE_PECLET,
                "q_for_target_peclet_m3_s": q_peclet,
                "q_capillary_limit_m3_s": q_capillary,
                "recommended_inlet_flow_m3_s": recommended_q,
                "achieved_sample_peclet": achieved_peclet,
                "estimated_pore_reynolds_p50": pore_reynolds,
                "estimated_capillary_number": capillary_number,
                "raw_pore_count": conversion["raw_pore_count"],
                "raw_throat_count": conversion["raw_throat_count"],
                "derived_node_count": virtual["derived_node_count"],
                "derived_edge_count": virtual["derived_edge_count"],
                "inlet_virtual_count": virtual["inlet_virtual_count"],
                "outlet_virtual_count": virtual["outlet_virtual_count"],
                "raw_pore_sha256": sha256(raw_pore),
                "raw_throat_sha256": sha256(raw_throat),
                "derived_pore_sha256": sha256(destination / "PNM_new_pore_virtual.txt"),
                "derived_throat_sha256": sha256(destination / "PNM_new_throat_virtual.txt"),
                "derived_connect_sha256": sha256(destination / "PNM_new_connect_virtual.txt"),
                "status": "passed",
            }
            normalization = conversion["geometry_normalization"]
            expected_volume_scale = geometry_scale ** 3
            pore_ratio = (normalization["derived_pore_volume_sum_um3"] /
                          normalization["source_pore_volume_sum_um3"])
            throat_ratio = (normalization["derived_throat_volume_sum_um3"] /
                            normalization["source_throat_volume_sum_um3"])
            if (not np.isclose(pore_ratio, expected_volume_scale, rtol=2e-14) or
                    not np.isclose(throat_ratio, expected_volume_scale,
                                   rtol=2e-14)):
                raise RuntimeError(f"{source}: dimensional geometry scaling audit failed")
            geometry_audit = {
                "format": "bubble_batch_geometry_normalization_audit_v1",
                "status": "passed",
                "policy": "normalize each immutable raw PNM to user-confirmed "
                          "1200^3 voxels at 7.5 um",
                "source_pnextract_domain_um": [source_domain] * 3,
                "target_image_dimensions_voxels": list(IMAGE_DIMENSIONS),
                "target_voxel_size_um": VOXEL_SIZE_UM,
                "target_physical_domain_um": [TARGET_DOMAIN_UM] * 3,
                "isotropic_length_scale": geometry_scale,
                "area_scale": geometry_scale ** 2,
                "volume_scale": expected_volume_scale,
                "source_metadata_conflict": not np.isclose(
                    source_domain, TARGET_DOMAIN_UM),
                "conversion_audit": str(destination / "conversion_audit.json"),
                "raw_files_unchanged": True,
            }
            atomic_json(destination / "input_geometry_audit.json", geometry_audit)
            atomic_json(destination / "batch_input_manifest.json", record)
            results.append(record)
    except Exception as error:
        atomic_json(output_root / "progress.json", {
            "stage": "prepare_batch35_inputs", "status": "failed",
            "completed": len(results), "total": len(complete),
            "error": str(error),
            "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        })
        raise

    summary = {
        "format": "bubble_batch35_input_inventory_v1",
        "status": "passed", "parent_case_count": len(parents),
        "complete_variant_count": len(results), "missing_inputs": missing,
        "variant_counts": variant_counts,
        "geometry_policy": "all effective samples are 1200^3 voxels at 7.5 um",
        "effective_physical_domain_um": [TARGET_DOMAIN_UM] * 3,
        "flow_rule": "min(target sample Pe_L, Ca cap)",
        "target_sample_peclet": TARGET_SAMPLE_PECLET,
        "maximum_capillary_number": MAX_CAPILLARY_NUMBER,
        "cases": results,
    }
    atomic_json(output_root / "batch_input_inventory.json", summary)
    columns = ["case_index", "parent_case", "variant", "nominal_domain_um",
               "voxel_size_um", "control_volume_porosity",
               "recommended_inlet_flow_m3_s", "achieved_sample_peclet",
               "estimated_pore_reynolds_p50", "estimated_capillary_number",
               "derived_node_count", "derived_edge_count", "status"]
    with (output_root / "batch_input_inventory.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in results:
            item = dict(row)
            item["nominal_domain_um"] = item["nominal_domain_um"][0]
            writer.writerow({key: item[key] for key in columns})
    atomic_json(output_root / "progress.json", {
        "stage": "prepare_batch35_inputs", "status": "passed",
        "completed": len(results), "total": len(complete),
        "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
