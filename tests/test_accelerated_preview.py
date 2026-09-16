#!/usr/bin/env python3
import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

binary = Path(sys.argv[1]).resolve()
source = Path(sys.argv[2]).resolve()
template = json.loads((source / "tests/data/accelerated_preview_small.json").read_text())

with tempfile.TemporaryDirectory(prefix="bubble_accelerated_preview_") as temporary:
    root = Path(temporary)
    output = root / "output"
    config = root / "config.json"
    template["io"]["output_directory"] = str(output)
    config.write_text(json.dumps(template) + "\n")
    run = subprocess.run([str(binary), "--config", str(config)], cwd=source.parent,
                         text=True, capture_output=True)
    if run.returncode:
        raise SystemExit(f"accelerated small run failed:\n{run.stdout}\n{run.stderr}")

    manifest = json.loads((output / "run_manifest.json").read_text())
    checkpoint = json.loads((output / "checkpoint.json").read_text())
    rows = list(csv.DictReader((output / "diagnostics.tsv").open(), delimiter="\t"))
    if not rows or manifest.get("termination_reason") != "residual_free_gas_fraction_reached":
        raise SystemExit("1% residual free-gas termination did not fire")
    initial = float(manifest["termination"]["initial_total_free_gas_moles"])
    free = [float(row["total_free_gas_moles"]) for row in rows]
    dissolved = [float(row["total_dissolved_moles"]) for row in rows]
    if free[-1] > 0.01 * initial or any(b > a * (1 + 1e-12) for a, b in zip(free, free[1:])):
        raise SystemExit("accelerated free gas is not monotonically decreasing to 1%")
    if any(b < a - 1e-25 for a, b in zip(dissolved, dissolved[1:])):
        raise SystemExit("accelerated dissolved inventory is not monotonic")
    if min(float(row["min_free_gas_moles"]) for row in rows) < 0:
        raise SystemExit("negative free gas appeared")
    if max(float(row["epsilon_N_adjusted"]) for row in rows) > 1e-6:
        raise SystemExit("component conservation failed")
    if max(float(row["epsilon_V_adjusted"]) for row in rows) > 1e-6:
        raise SystemExit("liquid-volume conservation failed")
    if any(float(row["mass_transfer_multiplier"]) != 100 for row in rows):
        raise SystemExit("diagnostics lost mass-transfer multiplier")
    if checkpoint.get("mass_transfer_multiplier") != 100:
        raise SystemExit("checkpoint lost mass-transfer multiplier")
    if checkpoint.get("gas_disappearance_events", 0) < 1:
        raise SystemExit("accelerated active-set disappearance event did not occur")
    if not any(state[5] for state in checkpoint["state"]):
        raise SystemExit("accepted disappearance event did not set permanently_inactive")
    if manifest.get("mass_transfer", {}).get("mass_transfer_multiplier") != 100:
        raise SystemExit("manifest lost mass-transfer multiplier")
    pvd = (output / "states.pvd").read_text()
    vtp = next(output.glob("state_*.vtp")).read_bytes()
    for field in (b"mass_transfer_multiplier", b"active_gas",
                  b"permanently_inactive", b"gas_moles", b"dissolved_moles"):
        if field not in vtp:
            raise SystemExit(f"VTP field missing: {field.decode()}")
    if "mass_transfer_multiplier=100" not in pvd:
        raise SystemExit("PVD lost mass-transfer metadata")
    print(f"PASS accelerated preview small termination steps={len(rows)} "
          f"final_fraction={free[-1]/initial:.9g}")
