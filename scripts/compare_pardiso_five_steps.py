#!/usr/bin/env python3
"""Compare existing Eigen and PARDISO five-step checkpoint sequences."""

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

eigen = Path(sys.argv[1])
pardiso = Path(sys.argv[2])
output = Path(sys.argv[3])


def table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


eigen_rows = table(eigen / "steps.tsv")
pardiso_rows = table(pardiso / "steps.tsv")
if len(eigen_rows) != 5 or len(pardiso_rows) != 5:
    raise SystemExit("both backends must have five accepted steps")

timing_fields = [
    "step_wall_seconds", "assembly_seconds", "csr_seconds", "setup_seconds",
    "symbolic_analysis_seconds", "numeric_factorization_seconds", "solve_seconds",
    "newton_iterations", "linear_iterations", "cpu_linear_relative_residual",
    "residual_scaled_inf", "relative_update_norm", "epsilon_N", "epsilon_V",
    "symbolic_analysis_calls", "numeric_factorization_calls",
]


def aggregate(rows):
    sums = {key: sum(float(row.get(key, 0.0)) for row in rows)
            for key in timing_fields}
    return {
        "steps": len(rows),
        "sum": sums,
        "mean": {key: value / len(rows) for key, value in sums.items()},
        "max": {key: max(float(row.get(key, 0.0)) for row in rows)
                for key in timing_fields},
    }


with (output / "pardiso_performance_comparison.tsv").open("w", newline="") as stream:
    fields = ["backend", "accepted_index", "global_step", "time_s",
              "accepted_dt_s", *timing_fields, "factor_nonzeros", "peak_memory_kib"]
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for backend, rows in (("eigen_bicgstab_ilut", eigen_rows),
                          ("pardiso", pardiso_rows)):
        for row in rows:
            writer.writerow({key: backend if key == "backend" else row.get(key, "0")
                             for key in fields})

eigen_summary = aggregate(eigen_rows)
pardiso_summary = aggregate(pardiso_rows)
eigen_manifest = json.loads((eigen / "run_manifest.json").read_text())
pardiso_manifest = json.loads((pardiso / "run_manifest.json").read_text())
eigen_checkpoint = json.loads((eigen / "checkpoint.json").read_text())
pardiso_checkpoint = json.loads((pardiso / "checkpoint.json").read_text())
if (eigen_checkpoint["step_number"], eigen_checkpoint["time_s"]) != (
        pardiso_checkpoint["step_number"], pardiso_checkpoint["time_s"]):
    raise SystemExit("final checkpoint time/step mismatch")

eigen_state = np.array([[float(value) if value is not None else np.nan
                         for value in row[:4]] for row in eigen_checkpoint["state"]])
pardiso_state = np.array([[float(value) if value is not None else np.nan
                           for value in row[:4]] for row in pardiso_checkpoint["state"]])
active = np.array([bool(row[4]) for row in eigen_checkpoint["state"]])
scales = np.array([100.0, 10000.0, 1.0, 1.0])
state_differences = {}
for column, name in enumerate(("Pl", "Pg", "Sw", "C")):
    mask = np.isfinite(eigen_state[:, column]) & np.isfinite(pardiso_state[:, column])
    difference = np.abs(eigen_state[mask, column] - pardiso_state[mask, column])
    denominator = np.maximum.reduce((np.abs(eigen_state[mask, column]),
                                     np.abs(pardiso_state[mask, column]),
                                     np.full(difference.shape, 1e-300)))
    state_differences[name] = {
        "maximum_absolute_difference": float(difference.max()),
        "maximum_relative_difference": float((difference / denominator).max()),
        "maximum_scaled_difference": float((difference / scales[column]).max()),
    }

config = json.loads((eigen / "resolved_config.json").read_text())
root = Path("/workspace/zz")
pore = np.loadtxt(root / config["io"]["pore_file"])
throat = np.loadtxt(root / config["io"]["throat_file"])
node_count = len(pore)
rows = np.r_[throat[:, 1].astype(int) - 1, throat[:, 2].astype(int) - 1]
columns = np.r_[throat[:, 2].astype(int) - 1, throat[:, 1].astype(int) - 1]
_, labels = connected_components(
    coo_matrix((np.ones(len(rows)), (rows, columns)), shape=(node_count, node_count)),
    directed=False)
virtual = pore[:, 1].astype(int) == 2
flags = pore[:, 10].astype(int)
candidates = set(labels[virtual & (flags == 1)]) & set(labels[virtual & (flags == 2)])
main_component = max(candidates,
                     key=lambda label: int(np.count_nonzero((labels == label) & ~virtual)))
volume = pore[:, 9] * float(config["geometry"]["input_volume_to_m3"])
real = (~virtual) & (volume > 0.0)
main_real = real & (labels == main_component)
gas_constant = 8.31446261815324
temperature = float(config["physics"]["T_K"])
compressibility = float(config["physics"]["Z"])
background = float(config["pressure_reference"]["background_abs_Pa"])


def inventory(state):
    pl, pg, sw, concentration = state.T
    del pl
    gas = active & real
    free = float(np.sum((background + pg[gas]) * volume[gas] * (1.0 - sw[gas]) /
                        (compressibility * gas_constant * temperature)))
    dissolved = float(np.sum(volume[real] * sw[real] * concentration[real]))
    saturation = float(np.sum(volume[main_real] * (1.0 - sw[main_real])) /
                       np.sum(volume[main_real]))
    return {"main_real_gas_saturation": saturation,
            "free_gas_mol": free, "dissolved_gas_mol": dissolved}


eigen_inventory = inventory(eigen_state)
pardiso_inventory = inventory(pardiso_state)
inventory_differences = {}
for key in eigen_inventory:
    absolute = abs(eigen_inventory[key] - pardiso_inventory[key])
    relative = absolute / max(abs(eigen_inventory[key]), abs(pardiso_inventory[key]), 1e-300)
    inventory_differences[key] = {"eigen": eigen_inventory[key],
                                  "pardiso": pardiso_inventory[key],
                                  "absolute_difference": absolute,
                                  "relative_difference": relative}

eigen_linear = (eigen_summary["sum"]["setup_seconds"] +
                eigen_summary["sum"]["solve_seconds"] +
                eigen_summary["sum"]["csr_seconds"])
pardiso_linear = (pardiso_summary["sum"]["setup_seconds"] +
                  pardiso_summary["sum"]["solve_seconds"] +
                  pardiso_summary["sum"]["csr_seconds"])
result = {
    "format": "bubble_eigen_pardiso_five_step_comparison_v1",
    "initial_checkpoint": pardiso_manifest["initial_checkpoint"],
    "initial_checkpoint_sha256": pardiso_manifest["initial_checkpoint_sha256"],
    "accepted_steps": 5,
    "pardiso_threads": pardiso_manifest["pardiso_threads"],
    "eigen": eigen_summary,
    "pardiso": pardiso_summary,
    "speedup_eigen_over_pardiso": {
        "process_wall": eigen_manifest["total_wall_seconds"] /
                        pardiso_manifest["total_wall_seconds"],
        "recorded_step_wall": eigen_summary["sum"]["step_wall_seconds"] /
                              pardiso_summary["sum"]["step_wall_seconds"],
        "complete_linear_stage": eigen_linear / pardiso_linear,
        "pure_solve": eigen_summary["sum"]["solve_seconds"] /
                      pardiso_summary["sum"]["solve_seconds"],
    },
    "state_differences": state_differences,
    "inventory_differences": inventory_differences,
    "compatibility": {
        "same_initial_checkpoint": eigen_manifest["initial_checkpoint_sha256"] ==
                                   pardiso_manifest["initial_checkpoint_sha256"],
        "same_final_time_and_step": True,
        "same_checkpoint_schema": sorted(eigen_checkpoint) == sorted(pardiso_checkpoint),
        "active_gas_mismatch_count": sum(
            bool(a[4]) != bool(b[4]) for a, b in
            zip(eigen_checkpoint["state"], pardiso_checkpoint["state"])),
    },
}
result["acceptance"] = {
    "state_scaled_below_1e-7": max(value["maximum_scaled_difference"]
                                    for value in state_differences.values()) <= 1e-7,
    "inventories_relative_below_1e-8": max(value["relative_difference"]
                                            for value in inventory_differences.values()) <= 1e-8,
    "both_nonlinear_gates_pass": max(eigen_summary["max"]["residual_scaled_inf"],
                                      pardiso_summary["max"]["residual_scaled_inf"]) <= 1e-8 and
                                 max(eigen_summary["max"]["relative_update_norm"],
                                     pardiso_summary["max"]["relative_update_norm"]) <= 1e-9,
    "both_conservation_pass": max(eigen_summary["max"]["epsilon_N"],
                                  pardiso_summary["max"]["epsilon_N"],
                                  eigen_summary["max"]["epsilon_V"],
                                  pardiso_summary["max"]["epsilon_V"]) <= 1e-6,
    "structure_reused_after_first_analysis":
        pardiso_summary["sum"]["symbolic_analysis_calls"] == 1.0,
}
(output / "pardiso_performance_comparison.json").write_text(
    json.dumps(result, indent=2) + "\n")
markdown = f"""# Eigen / PARDISO five-step comparison

| Metric | Eigen | PARDISO (16 threads) | Eigen/PARDISO speedup |
|---|---:|---:|---:|
| Five-step process wall [s] | {eigen_manifest['total_wall_seconds']:.6g} | {pardiso_manifest['total_wall_seconds']:.6g} | {result['speedup_eigen_over_pardiso']['process_wall']:.6g} |
| Recorded step wall [s] | {eigen_summary['sum']['step_wall_seconds']:.6g} | {pardiso_summary['sum']['step_wall_seconds']:.6g} | {result['speedup_eigen_over_pardiso']['recorded_step_wall']:.6g} |
| Complete linear stage [s] | {eigen_linear:.6g} | {pardiso_linear:.6g} | {result['speedup_eigen_over_pardiso']['complete_linear_stage']:.6g} |
| Pure solve [s] | {eigen_summary['sum']['solve_seconds']:.6g} | {pardiso_summary['sum']['solve_seconds']:.6g} | {result['speedup_eigen_over_pardiso']['pure_solve']:.6g} |
| Symbolic analyses | 0 | {pardiso_summary['sum']['symbolic_analysis_calls']:.0f} | |

All acceptance checks passed: `{all(result['acceptance'].values())}`.
"""
(output / "pardiso_performance_comparison.md").write_text(markdown)
if not all(result["acceptance"].values()):
    raise SystemExit("PARDISO five-step comparison failed acceptance")
print(json.dumps(result["speedup_eigen_over_pardiso"], indent=2))
