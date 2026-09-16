#!/usr/bin/env python3
"""Generate immutable per-network production configs from an audited inventory."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--config-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--pardiso-threads", type=int, default=8)
    parser.add_argument("--maximum-output-bytes", type=int,
                        default=4 * 1024 ** 3)
    parser.add_argument("--t-end-s", type=float, default=1.0e6)
    args = parser.parse_args()

    inventory_path = args.inventory.resolve()
    template_path = args.template.resolve()
    config_root = args.config_root.resolve()
    output_root = args.output_root.resolve()
    if config_root.exists():
        raise SystemExit(f"refusing to overwrite config root {config_root}")
    if args.pardiso_threads <= 0 or args.maximum_output_bytes <= 0:
        raise SystemExit("thread count and output byte cap must be positive")
    inventory = json.loads(inventory_path.read_text())
    if (inventory.get("status") != "passed" or
            inventory.get("variant_counts") != {"PNM-0.6": 35, "PNM-0.8": 35} or
            inventory.get("missing_inputs")):
        raise SystemExit("input inventory is not a complete audited 35+35 set")
    template = json.loads(template_path.read_text())
    config_root.mkdir(parents=True)
    entries = []
    for record in inventory["cases"]:
        variant = record["variant_slug"]
        run_id = f'{record["case_slug"]}_{variant}'
        derived = Path(record["derived_directory"]).resolve()
        output = output_root / run_id
        config = copy.deepcopy(template)
        config["case"] = {
            "name": run_id,
            "kind": "reference_physics_with_unvalidated_closures",
            "scientific_parameters_confirmed": False,
            "gas_component": "effective_dry_air",
            "liquid_composition": "pure_water",
            "notes": ("35-gradation comparison; immutable raw PNM normalized to the "
                      "user-confirmed common 1200^3 at 7.5 um geometry; global inlet "
                      "flow equation, degassed inlet, PARDISO, and conservative batch "
                      "microbubble retirement."),
        }
        config["io"].update({
            "pore_file": str(derived / "PNM_new_pore_virtual.txt"),
            "throat_file": str(derived / "PNM_new_throat_virtual.txt"),
            "connect_audit_file": str(derived / "PNM_new_connect_virtual.txt"),
            "input_geometry_audit_file": str(derived / "input_geometry_audit.json"),
            "virtual_boundary_geometry_audit_file": str(
                derived / "virtual_boundary_geometry_audit.json"),
            "output_directory": str(output),
            "checkpoint_interval_steps": 100,
            "vtk_interval_steps": 100000000,
            "vtk_output_times_s": [0.0],
            "checkpoint_output_times_s": [],
        })
        config["geometry"].update({
            "voxel_size_um": 7.5,
            "image_dimensions_voxels": [1200, 1200, 1200],
            "physical_domain_um": [9000, 9000, 9000],
        })
        config["input_provenance"] = {
            "raw_pore_file": record["source_directory"] + "/PNM_pb.txt",
            "raw_throat_file": record["source_directory"] + "/PNM_pt.txt",
            "base_pore_file": str(derived / "PNM_new_pore_base.txt"),
            "base_throat_file": str(derived / "PNM_new_throat_base.txt"),
            "base_connect_file": str(derived / "PNM_new_connect_base.txt"),
            "raw_pore_sha256": record["raw_pore_sha256"],
            "raw_throat_sha256": record["raw_throat_sha256"],
            "batch_input_manifest": str(derived / "batch_input_manifest.json"),
            "geometry_normalization_scale": record[
                "geometry_length_scale_to_effective_domain"],
        }
        config["boundary"]["inlet_flow_m3_s"] = record[
            "recommended_inlet_flow_m3_s"]
        config["linear"]["backend"] = "pardiso"
        config["linear"]["pardiso_threads"] = args.pardiso_threads
        config["time"]["t_end_s"] = args.t_end_s
        config["termination"]["maximum_physical_time_s"] = args.t_end_s
        config["termination"]["maximum_output_bytes"] = args.maximum_output_bytes
        config["active_set"]["batch_microbubble_retirement"] = {
            "enabled": True,
            "maximum_individual_free_gas_fraction": 1.0e-10,
            "maximum_batch_free_gas_fraction_per_step": 1.0e-8,
        }
        config["source_notes"] = {
            "geometry": ("Derived-only isotropic normalization to 9000 um; source "
                         f'pnextract domain was {record["source_pnextract_domain_um"]} um.'),
            "physics": ("25 C air-water configuration; mass_transfer_multiplier=1. "
                        "Henry coefficient, kL, trapped-bubble immobility and other "
                        "closures remain unconfirmed experimentally."),
            "boundary": ("One global inlet-flow equation; case-specific Q selected by "
                         "target Pe_L=555.555 with Ca<=1e-5; outlet Pl_rel=0; C_in=0 "
                         "and C_backflow=0; outlet diffusive zero gradient."),
            "initial": "Fresh saturation distribution from this PNM variant; starts at t=0.",
            "adaptive_dt": ("Frozen validated controller and conservative microbubble "
                            "batch-retirement thresholds; no tolerance relaxation."),
        }
        path = config_root / f"{run_id}.json"
        path.write_text(json.dumps(config, indent=2) + "\n")
        entries.append({
            "run_id": run_id,
            "parent_case": record["parent_case"],
            "variant": record["variant"],
            "config": str(path),
            "config_sha256": sha256(path),
            "output_directory": str(output),
            "inlet_flow_m3_s": record["recommended_inlet_flow_m3_s"],
            "sample_peclet": record["achieved_sample_peclet"],
            "capillary_number": record["estimated_capillary_number"],
            "node_count": record["derived_node_count"],
            "edge_count": record["derived_edge_count"],
        })
    index = {
        "format": "bubble_batch35_production_configs_v1",
        "status": "prepared",
        "inventory": str(inventory_path),
        "inventory_sha256": sha256(inventory_path),
        "template": str(template_path),
        "template_sha256": sha256(template_path),
        "output_root": str(output_root),
        "pardiso_threads_per_job": args.pardiso_threads,
        "maximum_output_bytes_per_case": args.maximum_output_bytes,
        "case_count": len(entries),
        "variant_counts": {"PNM-0.6": 35, "PNM-0.8": 35},
        "cases": entries,
    }
    index_path = config_root / "batch_cases.json"
    index_path.write_text(json.dumps(index, indent=2) + "\n")
    print(index_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
