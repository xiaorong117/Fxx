#!/usr/bin/env python3
"""Audit short PARDISO pilots before any 70-case production launch."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    index = json.loads(args.index.resolve().read_text())
    results = []
    for item in index["cases"]:
        output = Path(item["output_directory"])
        config = json.loads(Path(item["config"]).read_text())
        with (output / "diagnostics.tsv").open(newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        manifest = json.loads((output / "run_manifest.json").read_text())
        checkpoint = json.loads((output / "checkpoint.json").read_text())
        failures = []
        if len(rows) < 5:
            failures.append(f"only {len(rows)} accepted steps")
        required = ["dt_s", "newton_iterations", "timestep_retries",
                    "linear_residual", "residual_scaled_inf",
                    "relative_update_norm", "epsilon_N_adjusted",
                    "epsilon_V_adjusted", "free_mol_total",
                    "dissolved_mol_total", "gas_volume_total",
                    "inlet_liquid_inflow", "inlet_flow_target_m3_s"]
        for field in required:
            if field not in rows[0]:
                failures.append(f"missing diagnostics field {field}")
        numeric = {field: [float(row[field]) for row in rows]
                   for field in required if field in rows[0]}
        if any(not math.isfinite(value) for values in numeric.values()
               for value in values):
            failures.append("non-finite diagnostics value")
        if any(row.get("residual_criterion_passed") != "1" or
               row.get("update_criterion_passed") != "1" for row in rows):
            failures.append("Newton dual acceptance gate failed")
        mass_limit = config["acceptance"]["mass_balance_hard_limit"]
        volume_limit = config["acceptance"]["volume_balance_hard_limit"]
        max_mass = max(abs(v) for v in numeric["epsilon_N_adjusted"])
        max_volume = max(abs(v) for v in numeric["epsilon_V_adjusted"])
        if max_mass > mass_limit:
            failures.append(f"mass conservation {max_mass} > {mass_limit}")
        if max_volume > volume_limit:
            failures.append(f"volume conservation {max_volume} > {volume_limit}")
        target = config["boundary"]["inlet_flow_m3_s"]
        max_flow_relative_error = max(
            abs(v - target) / target for v in numeric["inlet_liquid_inflow"])
        if max_flow_relative_error > 1.0e-8:
            failures.append(f"global inlet-flow equation relative error "
                            f"{max_flow_relative_error}")
        if (manifest.get("exit_status") != 0 or
                manifest.get("linear_backend", {}).get("selected") != "pardiso"):
            failures.append("successful PARDISO manifest missing")
        active = manifest.get("active_set", {})
        if not active.get("batch_microbubble_retirement", {}).get("enabled", False):
            failures.append("batch microbubble retirement not effective")
        last_time = float(rows[-1]["time_s"])
        if not math.isclose(checkpoint["time_s"], last_time,
                            rel_tol=0.0, abs_tol=1e-15):
            failures.append("checkpoint does not match last diagnostics state")
        if not list(output.glob("*.vtp")):
            failures.append("initial VTP is missing")
        results.append({
            "run_id": item["run_id"], "status": "passed" if not failures else "failed",
            "accepted_steps": len(rows), "final_time_s": last_time,
            "minimum_dt_s": min(numeric["dt_s"]),
            "maximum_dt_s": max(numeric["dt_s"]),
            "total_retries": sum(int(float(v)) for v in numeric["timestep_retries"]),
            "maximum_mass_balance_error": max_mass,
            "maximum_volume_balance_error": max_volume,
            "maximum_inlet_flow_relative_error": max_flow_relative_error,
            "failures": failures,
        })
    report = {
        "format": "bubble_batch35_pilot_audit_v1",
        "status": "passed" if all(row["status"] == "passed" for row in results)
                  else "failed",
        "pilots": results,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
