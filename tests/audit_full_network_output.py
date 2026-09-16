#!/usr/bin/env python3
"""Audit the five-step full-network degenerate-geometry smoke run."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path


root = Path(sys.argv[1])
output = Path(sys.argv[2])

expected_raw_hashes = {
    "PNM_pb.txt": "72e32ad727f31912529f2b21a9c90bcdf4d9079a22ed07288fca891974f18282",
    "PNM_pt.txt": "24ac31276f29073fa8f2c0cfab4d6ea7ffd35eb4a79b7d0d9bbf2e950ca0055a",
}
for name, expected in expected_raw_hashes.items():
    actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit(f"raw input changed: {name}")

manifest = json.loads((output / "run_manifest.json").read_text())
if (manifest.get("exit_status") != 0 or manifest.get("accepted_steps") != 5 or
        manifest.get("final_time_s") != 5e-4):
    raise SystemExit("full-network manifest status/time/step mismatch")
network = manifest.get("network", {})
expected_network = {
    "control_volumes": 148144,
    "half_edges": 225771,
    "base_node_count": 146635,
    "base_edge_count": 224262,
    "virtual_boundary_nodes": 1509,
    "boundary_edges": 1509,
    "effective_zero_distance_clusters": 147630,
    "production_positive_edges": 225257,
    "zero_distance_constraints": 514,
    "nontrivial_zero_distance_clusters": 219,
    "massless_junctions": 2,
    "throats_with_two_zero_halves": 0,
    "connected_components": 711,
    "boundaryless_components": 672,
    "isolated_pores": 631,
    "final_pressure_gauges": 672,
    "ignored_boundary_to_boundary_edges": 0,
}
if network != expected_network:
    raise SystemExit(f"full-network topology mismatch: {network}")
if manifest.get("degenerate_geometry_formulation") != (
        "raw geometry unchanged; union-find zero-distance elimination plus massless junctions"):
    raise SystemExit("manifest lacks exact degenerate-geometry formulation")
configuration = manifest.get("configuration", {})
geometry = configuration.get("geometry", {})
if (geometry.get("input_basis") != "physical" or
        geometry.get("input_length_unit") != "um" or
        geometry.get("voxel_size_um") != 7.5 or
        geometry.get("image_dimensions_voxels") != [1200, 1200, 1200] or
        geometry.get("physical_domain_um") != [9000, 9000, 9000] or
        geometry.get("input_length_to_m") != 1e-6 or
        geometry.get("input_area_to_m2") != 1e-12 or
        geometry.get("input_volume_to_m3") != 1e-18):
    raise SystemExit("full-network physical geometry metadata mismatch")
virtual = manifest.get("virtual_boundary", {})
if (virtual.get("enabled") is not True or virtual.get("boundary_axis") != "x" or
        virtual.get("extension_voxels") != 1.0 or virtual.get("voxel_size_um") != 7.5 or
        virtual.get("boundary_throat_length_um") != 7.5 or
        virtual.get("inlet_virtual_nodes") != 724 or
        virtual.get("outlet_virtual_nodes") != 785):
    raise SystemExit("full-network virtual-boundary manifest mismatch")
if set(manifest.get("input_provenance", {})) != {
        "raw_pore_file", "raw_throat_file", "base_pore_file", "base_throat_file",
        "base_connect_file", "input_geometry_audit_file",
        "virtual_boundary_geometry_audit_file"}:
    raise SystemExit("full-network original/derived provenance is incomplete")
nonlinear = configuration.get("nonlinear", {})
residual_tolerance = nonlinear.get("residual_tolerance")
update_tolerance = nonlinear.get("update_tolerance")
if residual_tolerance != 1e-8 or update_tolerance != 1e-9:
    raise SystemExit("full-network run did not use the required nonlinear tolerances")
for name, expected in manifest.get("outputs", {}).items():
    path = Path(name)
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise SystemExit(f"output checksum mismatch: {path}")

with (output / "diagnostics.tsv").open(newline="") as stream:
    rows = list(csv.DictReader(stream, delimiter="\t"))
if len(rows) != 5 or [int(row["step"]) for row in rows] != [1, 2, 3, 4, 5]:
    raise SystemExit("full-network diagnostics do not contain exactly five accepted steps")
previous_free = math.inf
block_names = ("Rl", "Rd", "Rg", "Rc")
for row in rows:
    numeric = {key: float(value) for key, value in row.items()
               if key not in {"step", "dt_change_reason", "termination_reason"}}
    if not all(math.isfinite(value) for value in numeric.values()):
        raise SystemExit("full-network diagnostics contain NaN/Inf")
    if numeric["epsilon_N_equation"] > 1e-6 or numeric["epsilon_V_equation"] > 1e-6:
        raise SystemExit("full-network conservation hard limit exceeded")
    if numeric["epsilon_N_adjusted"] > 1e-6 or numeric["epsilon_V_adjusted"] > 1e-6:
        raise SystemExit("full-network adjusted conservation hard limit exceeded")
    if numeric["residual_scaled_inf"] > residual_tolerance:
        raise SystemExit("full-network total nonlinear residual tolerance exceeded")
    if numeric["relative_update_norm"] > update_tolerance:
        raise SystemExit("full-network relative update tolerance exceeded")
    if numeric["update_scaled_inf"] > numeric["update_tolerance_limit"]:
        raise SystemExit("full-network scaled update limit exceeded")
    if (numeric["update_criterion_ratio"] > 1.0 or
            numeric["residual_criterion_passed"] != 1.0 or
            numeric["update_criterion_passed"] != 1.0):
        raise SystemExit("full-network step lacks both convergence gates")
    if max(numeric[f"{name}_scaled_inf"] for name in block_names) > residual_tolerance:
        raise SystemExit("full-network residual block tolerance exceeded")
    if any(int(row[f"{name}_max_node_id"]) <= 0 for name in block_names):
        raise SystemExit("full-network residual maximum lacks an original node ID")
    if numeric["free_mol"] > previous_free + 1e-12 * max(previous_free, 1e-30):
        raise SystemExit("free gas increased")
    previous_free = numeric["free_mol"]

convergence = manifest.get("nonlinear_convergence_audit", {})
if (convergence.get("residual_tolerance") != residual_tolerance or
        convergence.get("update_tolerance") != update_tolerance or
        convergence.get("all_steps_passed_residual_and_update_criteria") is not True or
        len(convergence.get("steps", [])) != 5):
    raise SystemExit("manifest nonlinear convergence audit is incomplete")
max_residual_row = max(rows, key=lambda row: float(row["residual_scaled_inf"]))
if (convergence.get("maximum_residual_scaled_inf") !=
        float(max_residual_row["residual_scaled_inf"]) or
        convergence.get("maximum_residual_step") != int(max_residual_row["step"])):
    raise SystemExit("manifest total residual maximum disagrees with diagnostics")
max_update_row = max(rows, key=lambda row: float(row["update_criterion_ratio"]))
if (convergence.get("maximum_update_criterion_ratio") !=
        float(max_update_row["update_criterion_ratio"]) or
        convergence.get("maximum_update_criterion_ratio_step") != int(max_update_row["step"])):
    raise SystemExit("manifest update maximum disagrees with diagnostics")
for name in block_names:
    maximum = max(rows, key=lambda row: float(row[f"{name}_scaled_inf"]))
    expected = {
        "scaled_inf": float(maximum[f"{name}_scaled_inf"]),
        "step": int(maximum["step"]),
        "node_id": int(maximum[f"{name}_max_node_id"]),
    }
    if convergence.get("block_maxima", {}).get(name) != expected:
        raise SystemExit(f"manifest {name} maximum disagrees with diagnostics")

degenerate = json.loads((output / "degenerate_geometry_audit.json").read_text())
if (degenerate.get("format") != "bubble_degenerate_geometry_audit_v1" or
        degenerate.get("zero_distance_edge_count") != 514 or
        degenerate.get("massless_junction_node_ids") != [32123, 140815] or
        len(degenerate.get("zero_distance_edges", [])) != 514 or
        len(degenerate.get("nontrivial_clusters", [])) != 219):
    raise SystemExit("degenerate geometry audit mapping is incomplete")
if any(value > 1e-15 for key, value in degenerate["aggregation_invariants"].items()
       if key.endswith("relative_error")):
    raise SystemExit("cluster aggregation changed an initial invariant")

checkpoint = json.loads((output / "checkpoint.json").read_text())
if (checkpoint.get("format") != "bubble_checkpoint_v6" or
        checkpoint.get("step_number") != 5 or checkpoint.get("time_s") != 5e-4):
    raise SystemExit("full-network checkpoint status mismatch")
state = checkpoint["state"]
if len(state) != 148144:
    raise SystemExit("full-network checkpoint node count mismatch")
precise_state = checkpoint.get("state_extended_precision")
if not isinstance(precise_state, list) or len(precise_state) != len(state):
    raise SystemExit("full-network checkpoint lacks the extended-precision state")
for legacy, precise in zip(state, precise_state):
    if len(precise) != 4 or not all(isinstance(value, str) for value in precise):
        raise SystemExit("malformed full-network extended-precision checkpoint state")
    if not all(math.isfinite(float(precise[index])) for index in (0, 2, 3)):
        raise SystemExit("non-finite full-network extended-precision checkpoint state")
    if legacy[4] and not math.isfinite(float(precise[1])):
        raise SystemExit("non-finite full-network extended-precision gas pressure")
for cluster in degenerate["nontrivial_clusters"]:
    members = [state[node_id - 1] for node_id in cluster["member_node_ids"]]
    precise_members = [precise_state[node_id - 1] for node_id in cluster["member_node_ids"]]
    pl, concentration = members[0][0], members[0][3]
    if any(member[0] != pl or member[3] != concentration for member in members[1:]):
        raise SystemExit(f"cluster {cluster['cluster_id']} does not share Pl/C")
    precise_pl, precise_concentration = precise_members[0][0], precise_members[0][3]
    if any(member[0] != precise_pl or member[3] != precise_concentration
           for member in precise_members[1:]):
        raise SystemExit(f"cluster {cluster['cluster_id']} does not share precise Pl/C")
for node_id in (32123, 140815):
    row = state[node_id - 1]
    if row[2] != 1.0 or row[4]:
        raise SystemExit(f"massless junction {node_id} has gas/storage saturation state")
for node_id in range(146636, 148145):
    row = state[node_id - 1]
    if row[2] != 1.0 or row[4]:
        raise SystemExit(f"virtual reservoir {node_id} has gas/storage saturation state")

virtual_audit = json.loads((output / "virtual_boundary_geometry_audit.json").read_text())
if (virtual_audit.get("inlet_virtual_count") != 724 or
        virtual_audit.get("outlet_virtual_count") != 785 or
        virtual_audit.get("boundary_connection_count") != 1509 or
        virtual_audit.get("derived_node_count") != 148144 or
        virtual_audit.get("derived_edge_count") != 225771 or
        virtual_audit.get("new_zero_length_boundary_edge_count") != 0 or
        virtual_audit.get("virtual_nodes_in_nontrivial_zero_distance_cluster") != 0 or
        virtual_audit.get("virtual_node_coordination_distribution") != {"1": 1509} or
        virtual_audit.get("isolated_former_boundary_count") != 35 or
        virtual_audit.get("virtual_volume_counted_as_internal_storage") is not False or
        not math.isclose(virtual_audit.get("connected_real_virtual_clearance_min_um"),
                         7.5, abs_tol=2e-12)):
    raise SystemExit("full-network virtual-boundary geometry audit is incomplete")
isolated = virtual_audit.get("isolated_former_boundary_nodes", [])
if (virtual_audit.get("isolated_former_inlet_count") != 14 or
        virtual_audit.get("isolated_former_outlet_count") != 21 or
        [item["real_node_id"] for item in isolated if item["input_active_gas"]] != [18230] or
        any(item["base_component_size"] != 1 or item["derived_component_size"] != 2 or
            item["in_largest_sample_component"] or item["singular_equation_expected"]
            for item in isolated)):
    raise SystemExit("isolated former boundary classification is incomplete")
overlap = virtual_audit.get("overlap_counts", {})
if overlap != {
        "virtual_nonconnected_real": 96,
        "virtual_nonconnected_real_new_vs_source_layout": 75,
        "virtual_virtual": 1,
        "virtual_virtual_new_vs_real_boundary_layout": 1,
}:
    raise SystemExit("full-network virtual overlap classification changed")
unit_audit = json.loads((output / "input_geometry_unit_audit.json").read_text())
if (unit_audit.get("conclusion") !=
        "consistent_with_physical_um_not_voxel_indices_or_nm" or
        unit_audit.get("physical_domain_um") != [9000.0, 9000.0, 9000.0]):
    raise SystemExit("full-network copied unit audit is incomplete")
if not (output / "run.log").is_file():
    raise SystemExit("full-network run.log missing")

required_point = {
    "Pl", "Pg", "Sw", "Sg", "C", "n_g", "ntr", "active_gas",
    "node_id", "source_type", "source_id", "boundary_flag",
    "is_virtual_boundary", "pressure_gauge",
    "zero_distance_cluster", "massless_junction",
}
required_cell = {
    "Q", "F_adv", "F_diff", "source_old_throat_id", "side",
    "is_boundary_edge", "edge_length", "edge_length_voxels",
    "zero_distance_constraint",
}
seen: set[str] = set()
virtual_point_count = 0
source_type_two_count = 0
boundary_edge_count = 0
active_name: str | None = None
remaining = 0
index = 0
zero_edge_indices = {item["edge_id"] - 1 for item in degenerate["zero_distance_edges"]}
data_count = 0
with (output / "state_000005.vtk").open() as stream:
    for line in stream:
        if line.startswith(("POINT_DATA ", "CELL_DATA ")):
            data_count = int(line.split()[1])
        if line.startswith("SCALARS "):
            active_name = line.split()[1]
            seen.add(active_name)
            remaining = data_count
            index = 0
            continue
        if active_name is None or line.startswith("LOOKUP_TABLE"):
            continue
        value = float(line)
        if not math.isfinite(value):
            raise SystemExit(f"VTK contains NaN/Inf in {active_name}")
        if active_name == "Pl" and value <= 0:
            raise SystemExit("VTK contains non-positive liquid pressure")
        if active_name in {"C", "ntr"} and value < -1e-12:
            raise SystemExit(f"VTK contains negative {active_name}")
        if active_name == "Sw" and not (0 < value <= 1):
            raise SystemExit("VTK contains invalid Sw")
        if active_name in {"Q", "F_adv", "F_diff"} and index in zero_edge_indices and value != 0:
            raise SystemExit(f"zero-distance edge has assembled {active_name}")
        if active_name == "is_virtual_boundary" and value == 1.0:
            virtual_point_count += 1
        if active_name == "source_type" and value == 2.0:
            source_type_two_count += 1
        if active_name == "is_boundary_edge" and value == 1.0:
            boundary_edge_count += 1
        if index >= 224262 and active_name == "edge_length" and not math.isclose(
                value, 7.5e-6, rel_tol=1e-14):
            raise SystemExit("VTK boundary edge length is not 7.5 um")
        if index >= 224262 and active_name == "edge_length_voxels" and not math.isclose(
                value, 1.0, rel_tol=1e-14):
            raise SystemExit("VTK boundary edge length is not one voxel")
        if index >= 224262 and active_name == "source_old_throat_id" and value != 0.0:
            raise SystemExit("VTK generated boundary edge has a false source throat ID")
        index += 1
        remaining -= 1
        if remaining == 0:
            active_name = None
if not (required_point | required_cell).issubset(seen):
    raise SystemExit("full-network VTK lacks required fields")
if (virtual_point_count != 1509 or source_type_two_count != 1509 or
        boundary_edge_count != 1509):
    raise SystemExit("full-network VTK virtual node/edge classification count mismatch")

print("PASS full-network 5-step degenerate-geometry/output/conservation audit")
