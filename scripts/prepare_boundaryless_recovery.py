#!/usr/bin/env python3
"""Prepare reproducible configs for the one approved boundaryless recovery."""

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
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--derived-input", type=Path, required=True)
    parser.add_argument("--config-directory", type=Path, required=True)
    parser.add_argument("--verification-root", type=Path, required=True)
    parser.add_argument("--production-output", type=Path, required=True)
    args = parser.parse_args()
    base_path = args.base_config.resolve(); derived = args.derived_input.resolve()
    config_dir = args.config_directory.resolve()
    if config_dir.exists():
        raise SystemExit(f"refusing to overwrite {config_dir}")
    config_dir.mkdir(parents=True)
    base = json.loads(base_path.read_text())

    def make(name: str, output: Path, end: float, filtered: bool,
             fixed_dt: bool = False) -> Path:
        value = copy.deepcopy(base)
        value["case"]["name"] = name
        value["case"]["kind"] = (
            "reference_physics_with_unvalidated_closures" if "production" in name
            else "numerical_smoke_only")
        value["case"]["notes"] = (
            "Case-specific boundaryless-component recovery. One audited three-node "
            "component is frozen at t=0 and excluded from dynamics; all boundary-connected "
            "components retain the original equations, physics and tolerances."
            if filtered else
            "Short unchanged-network comparator for boundaryless-component recovery; not production data.")
        value["io"]["output_directory"] = str(output.resolve())
        value["io"]["checkpoint_interval_steps"] = 1
        value["io"]["vtk_interval_steps"] = 100000000
        value["io"]["vtk_output_times_s"] = [0.0]
        value["time"]["t_end_s"] = end
        value["termination"]["maximum_physical_time_s"] = end
        if fixed_dt:
            value["time"].update({
                "dt_initial_s": 0.001, "dt_min_s": 0.001,
                "dt_max_s": 0.001, "growth_factor": 1.0,
                "fast_newton_iterations": 0})
            value["time_step_control"].update({
                "enabled": True, "dt_min_s": 0.001,
                "dt_max_s": 0.001, "growth_factor": 1.0})
        if filtered:
            value["io"]["pore_file"] = str(derived / "PNM_new_pore_virtual_dynamic.txt")
            value["io"]["throat_file"] = str(derived / "PNM_new_throat_virtual_dynamic.txt")
            value["io"]["connect_audit_file"] = str(derived / "PNM_new_connect_virtual_dynamic.txt")
            value["io"]["input_geometry_audit_file"] = str(
                derived / "boundaryless_component_freeze_audit.json")
            value["io"]["virtual_boundary_geometry_audit_file"] = str(
                derived / "filtered_virtual_boundary_audit.json")
            value["input_provenance"]["boundaryless_component_freeze_audit"] = str(
                derived / "boundaryless_component_freeze_audit.json")
            value["input_provenance"]["parent_failed_config"] = str(base_path)
        path = config_dir / f"{name}.json"
        path.write_text(json.dumps(value, indent=2) + "\n")
        return path

    verification = args.verification_root.resolve()
    configs = [
        make("85A15D_p08_original_equivalence_t0p02",
             verification / "original_equivalence_t0p02", 0.02, False, True),
        make("85A15D_p08_filtered_equivalence_t0p02",
             verification / "filtered_equivalence_t0p02", 0.02, True, True),
        make("85A15D_p08_filtered_cross_failure_pilot_t0p6",
             verification / "filtered_cross_failure_pilot_t0p6", 0.6, True),
        make("85A15D_p08_boundaryless_frozen_production",
             args.production_output.resolve(), 1.0e6, True),
    ]
    manifest = {
        "format": "bubble_boundaryless_recovery_configs_v1",
        "status": "prepared",
        "base_config": str(base_path), "base_config_sha256": sha256(base_path),
        "derived_input": str(derived),
        "freeze_audit_sha256": sha256(derived / "boundaryless_component_freeze_audit.json"),
        "configs": [{"path": str(path), "sha256": sha256(path)} for path in configs],
    }
    (config_dir / "recovery_configs_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
