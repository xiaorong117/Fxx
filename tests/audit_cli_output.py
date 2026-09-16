#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

root = Path(sys.argv[1])
required = {
    "step", "time_s", "dt_s", "newton_iterations", "linear_iterations",
    "residual_scaled_inf", "state_scaled_inf", "update_scaled_inf",
    "relative_update_norm", "update_tolerance_limit", "update_criterion_ratio",
    "residual_criterion_passed", "update_criterion_passed",
    "Rl_scaled_inf", "Rd_scaled_inf", "Rg_scaled_inf", "Rc_scaled_inf",
    "Rl_max_node_id", "Rd_max_node_id", "Rg_max_node_id", "Rc_max_node_id",
    "boundary_liquid_outflow_m3_s", "dissolved_mol", "free_mol",
    "boundary_component_outflow_mol_s", "epsilon_N_equation", "epsilon_V_equation",
    "epsilon_N_adjusted", "epsilon_V_adjusted", "active_gas_count",
    "threshold_moles_transferred", "active_set_mole_ledger",
    "active_set_liquid_volume_ledger", "retired_bubble_count",
    "batched_retired_bubble_count", "batched_retired_moles",
    "batched_retired_gas_volume_m3", "batch_aggregate_free_gas_fraction",
    "event_timestep_localized",
}
with (root / "diagnostics.tsv").open(newline="") as stream:
    table = list(csv.DictReader(stream, delimiter="\t"))
if len(table) != 5 or not required.issubset(table[0]):
    raise SystemExit("TSV field or accepted-step count audit failed")
for row in table:
    for key, value in row.items():
        if key not in {"step", "dt_change_reason", "termination_reason"} and not math.isfinite(float(value)):
            raise SystemExit(f"TSV contains non-finite value: {key}")
    if float(row["epsilon_N_equation"]) > 1e-10 or float(row["epsilon_V_equation"]) > 1e-10:
        raise SystemExit("small CLI conservation audit failed")
    if (float(row["residual_scaled_inf"]) > 1e-8 or
            float(row["relative_update_norm"]) > 1e-9 or
            float(row["update_scaled_inf"]) > float(row["update_tolerance_limit"]) or
            float(row["update_criterion_ratio"]) > 1.0 or
            int(row["residual_criterion_passed"]) != 1 or
            int(row["update_criterion_passed"]) != 1):
        raise SystemExit("small CLI accepted a step without both nonlinear convergence gates")
manifest = json.loads((root / "run_manifest.json").read_text())
if manifest["exit_status"] != 0 or manifest["accepted_steps"] != 5:
    raise SystemExit("manifest status/step audit failed")
batch_manifest = manifest.get("active_set", {}).get("batch_microbubble_retirement", {})
if batch_manifest.get("enabled") is not False or batch_manifest.get("event_log") != "gas_disappearance_events.tsv":
    raise SystemExit("manifest batch-retirement configuration missing")
expected_network = {
    "control_volumes": 4, "half_edges": 3,
    "effective_zero_distance_clusters": 4, "production_positive_edges": 3,
    "zero_distance_constraints": 0, "massless_junctions": 0,
    "nontrivial_zero_distance_clusters": 0,
    "throats_with_two_zero_halves": 0,
    "ignored_boundary_to_boundary_edges": 0, "connected_components": 1,
    "boundaryless_components": 0, "isolated_pores": 0,
    "final_pressure_gauges": 0,
    "virtual_boundary_nodes": 0, "boundary_edges": 0,
    "base_node_count": 4, "base_edge_count": 3,
}
if manifest["network"] != expected_network:
    raise SystemExit("manifest network audit failed")
if "union-find zero-distance elimination" not in manifest.get(
        "degenerate_geometry_formulation", ""):
    raise SystemExit("manifest lacks degenerate-geometry formulation")
if len(manifest["inputs"]) != 4 or not manifest["run_command"] or not manifest.get("run_argv"):
    raise SystemExit("manifest reproducibility fields missing")
build = manifest.get("build", {})
if (build.get("cxx_standard") != 17 or build.get("build_type") not in {"Release", "Debug"}
        or len(build.get("project_source_tree_sha256", "")) != 64):
    raise SystemExit("manifest build/source provenance missing")
for name, expected in manifest["outputs"].items():
    path = Path(name)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit(f"output checksum mismatch: {path}")
vtk = (root / "state_000005.vtk").read_text()
if " nan" in vtk.lower() or " inf" in vtk.lower():
    raise SystemExit("VTK contains NaN/Inf")
for field in ("Pl", "Pg", "Sw", "Sg", "C", "n_g", "ntr", "active_gas",
              "node_id", "source_type", "source_id", "boundary_flag",
              "is_virtual_boundary", "pressure_gauge",
              "zero_distance_cluster", "massless_junction",
              "Q", "F_adv", "F_diff",
              "source_old_throat_id", "is_boundary_edge", "edge_length",
              "edge_length_voxels", "side", "zero_distance_constraint"):
    if f"SCALARS {field} " not in vtk:
        raise SystemExit(f"VTK field missing: {field}")
for field in ("normalized_x", "segment_id", "local_Sg", "initial_Sg",
              "free_gas_moles_lost", "cumulative_interphase_transfer",
              "adjacent_constriction_ratio_p50",
              "batch_microbubble_retirement_enabled",
              "batch_individual_free_gas_fraction",
              "batch_total_free_gas_fraction_per_step"):
    if f"SCALARS {field} " not in vtk:
        raise SystemExit(f"scientific VTK field missing: {field}")
checkpoint = json.loads((root / "checkpoint.json").read_text())
if checkpoint.get("format") != "bubble_checkpoint_v6" or not all(
    key in checkpoint for key in (
        "next_dt_s", "consecutive_fast_steps", "provenance",
        "state_extended_precision", "batch_microbubble_retirement_enabled",
        "batch_individual_free_gas_fraction",
        "batch_total_free_gas_fraction_per_step",
    )
):
    raise SystemExit("checkpoint lacks adaptive or extended-precision continuation state")
if len(checkpoint["state_extended_precision"]) != len(checkpoint["state"]):
    raise SystemExit("checkpoint extended-precision state count mismatch")
for legacy, precise in zip(checkpoint["state"], checkpoint["state_extended_precision"]):
    if len(precise) != 4 or not all(isinstance(value, str) for value in precise):
        raise SystemExit("malformed checkpoint extended-precision state")
    if not all(math.isfinite(float(precise[index])) for index in (0, 2, 3)):
        raise SystemExit("non-finite checkpoint extended-precision active state")
    if legacy[4] and not math.isfinite(float(precise[1])):
        raise SystemExit("non-finite checkpoint extended-precision gas pressure")
if len(checkpoint.get("cumulative_interphase_transfer_by_node_mol", [])) != len(checkpoint["state"]):
    raise SystemExit("checkpoint lacks per-node cumulative transfer")

with (root / "global_history.tsv").open(newline="") as stream:
    global_rows = list(csv.DictReader(stream, delimiter="\t"))
with (root / "segment_history.tsv").open(newline="") as stream:
    segment_rows = list(csv.DictReader(stream, delimiter="\t"))
with (root / "segment_geometry.tsv").open(newline="") as stream:
    geometry_rows = list(csv.DictReader(stream, delimiter="\t"))
with (root / "balance_history.tsv").open(newline="") as stream:
    balance_rows = list(csv.DictReader(stream, delimiter="\t"))
if len(global_rows) != 5 or len(balance_rows) != 5 or len(geometry_rows) != 20 or len(segment_rows) != 100:
    raise SystemExit("scientific history row counts are incomplete")
required_global = {"main_real_gas_saturation", "cumulative_interphase_transfer_mol",
    "free_gas_moles_lost_mol", "assembly_seconds", "linear_setup_seconds",
    "linear_solve_seconds", "cpu_linear_relative_residual", "retired_bubble_count",
    "batched_retired_bubble_count", "batch_aggregate_free_gas_fraction"}
if not required_global.issubset(global_rows[0]):
    raise SystemExit("global scientific fields are incomplete")
with (root / "gas_disappearance_events.tsv").open(newline="") as stream:
    event_reader = csv.DictReader(stream, delimiter="\t")
    if not {"event_node_id", "segment_id", "event_mode",
            "gas_moles_transferred_to_dissolved",
            "batch_aggregate_free_gas_fraction"}.issubset(event_reader.fieldnames or []):
        raise SystemExit("gas-disappearance long-table schema is incomplete")
for row in balance_rows:
    if int(row["component_balance_pass"]) != 1 or int(row["liquid_balance_pass"]) != 1:
        raise SystemExit("readable balance history failed")
for step in range(1, 6):
    rows = [r for r in segment_rows if int(r["step"]) == step]
    gas = sum(float(r["segment_gas_volume_m3"]) for r in rows)
    free = sum(float(r["segment_free_gas_mol"]) for r in rows)
    dissolved = sum(float(r["segment_dissolved_gas_mol"]) for r in rows)
    global_row = global_rows[step - 1]
    for actual, expected in ((gas, float(global_row["gas_volume_m3"])),
                             (free, float(global_row["free_gas_mol"])),
                             (dissolved, float(global_row["dissolved_gas_mol"]))):
        if not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-30):
            raise SystemExit("segment/global inventory invariant failed")
degenerate = json.loads((root / "degenerate_geometry_audit.json").read_text())
if (degenerate.get("format") != "bubble_degenerate_geometry_audit_v1" or
        degenerate.get("effective_cluster_count") != 4 or
        any(value != 0.0 for key, value in
            degenerate.get("aggregation_invariants", {}).items()
            if key.endswith("relative_error"))):
    raise SystemExit("degenerate geometry audit is incomplete")
print("PASS CLI TSV/VTK/checkpoint/manifest/checksum audit")
