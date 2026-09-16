#!/usr/bin/env python3
"""Exercise the production CLI restart path with an adaptive timestep."""

import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


binary = Path(sys.argv[1]).resolve()
source = Path(sys.argv[2]).resolve()
base = json.loads((source / "tests/data/small_cli.json").read_text())


def write_config(path: Path, output: Path, end: float) -> None:
    config = json.loads(json.dumps(base))
    for key in ("pore_file", "throat_file", "connect_audit_file"):
        config["io"][key] = str((source.parent / config["io"][key]).resolve())
    config["io"]["output_directory"] = str(output)
    config["time"].update({
        "t_end_s": end,
        "dt_initial_s": 1.0e-4,
        "dt_max_s": 4.0e-4,
        "growth_factor": 2.0,
        "fast_newton_iterations": 30,
    })
    path.write_text(json.dumps(config, indent=2) + "\n")


def run(config: Path, restart: Path | None = None) -> None:
    command = [str(binary), "--config", str(config)]
    if restart is not None:
        command.extend(["--restart", str(restart)])
    completed = subprocess.run(command, text=True, capture_output=True)
    if completed.returncode:
        raise SystemExit(
            f"CLI restart command failed ({completed.returncode}):\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )


def run_expect_failure(config: Path, restart: Path, expected: str) -> None:
    completed = subprocess.run(
        [str(binary), "--config", str(config), "--restart", str(restart)],
        text=True,
        capture_output=True,
    )
    if completed.returncode == 0 or expected not in completed.stdout + completed.stderr:
        raise SystemExit(
            "expected restart rejection was not observed: "
            f"expected={expected!r}, returncode={completed.returncode}, "
            f"stdout={completed.stdout!r}, stderr={completed.stderr!r}"
        )


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=2.0e-14, abs_tol=2.0e-16)


with tempfile.TemporaryDirectory(prefix="bubble_cli_restart_") as temporary:
    root = Path(temporary)
    continuous_output = root / "continuous"
    restarted_output = root / "restarted"
    continuous_config = root / "continuous.json"
    first_config = root / "first.json"
    restart_config = root / "restart.json"
    incompatible_config = root / "incompatible.json"
    relocated_config = root / "relocated.json"
    relocated_output = root / "relocated"
    write_config(continuous_config, continuous_output, 7.5e-4)
    write_config(first_config, restarted_output, 2.0e-4)
    write_config(restart_config, restarted_output, 7.5e-4)
    write_config(incompatible_config, restarted_output, 7.5e-4)
    write_config(relocated_config, relocated_output, 4.0e-4)
    incompatible = json.loads(incompatible_config.read_text())
    incompatible["physics"]["Hcp_mol_m3_Pa"] *= 2.0
    incompatible_config.write_text(json.dumps(incompatible, indent=2) + "\n")

    run(continuous_config)
    run(first_config)
    checkpoint_path = restarted_output / "checkpoint.json"
    saved_checkpoint_path = root / "checkpoint_at_step_2.json"
    saved_checkpoint_path.write_bytes(checkpoint_path.read_bytes())
    restart_hash = hashlib.sha256(checkpoint_path.read_bytes()).hexdigest()
    first_checkpoint = json.loads(checkpoint_path.read_text())
    if first_checkpoint.get("format") != "bubble_checkpoint_v6":
        raise SystemExit("adaptive checkpoint did not use v6 format")
    if first_checkpoint.get("step_number") != 2 or not close(
        first_checkpoint.get("next_dt_s", -1.0), 2.0e-4
    ):
        raise SystemExit("checkpoint did not retain the grown next timestep")
    for key in ("cumulative_liquid_inflow_m3", "cumulative_liquid_outflow_m3",
                "cumulative_interphase_transfer_mol",
                "cumulative_interphase_transfer_by_node_mol"):
        if key not in first_checkpoint:
            raise SystemExit(f"checkpoint v6 lacks {key}")

    run_expect_failure(
        incompatible_config, checkpoint_path, "checkpoint provenance does not match"
    )
    corrupt_checkpoint = root / "corrupt_reactivation.json"
    corrupt = json.loads(saved_checkpoint_path.read_text())
    corrupt["state"][2][1] = 102000.0
    corrupt["state"][2][2] = 0.5
    corrupt["state"][2][4] = True
    # v3 checkpoints restore the extended-precision state in preference to the
    # binary64 compatibility array, so corrupt both representations to reach
    # the intended irreversible-reactivation validation.
    corrupt["state_extended_precision"][2][1] = "102000"
    corrupt["state_extended_precision"][2][2] = "0.5"
    corrupt_checkpoint.write_text(json.dumps(corrupt) + "\n")
    run_expect_failure(
        restart_config, corrupt_checkpoint, "forbidden activation of an initially liquid node"
    )
    run(relocated_config, saved_checkpoint_path)
    relocated_diagnostics = (relocated_output / "diagnostics.tsv").read_text().splitlines()
    if len(relocated_diagnostics) != 2 or not relocated_diagnostics[0].startswith("step\t"):
        raise SystemExit("restart into a new output directory did not create a valid TSV header")
    for name in ("global_history.tsv", "segment_history.tsv", "balance_history.tsv"):
        lines = (relocated_output / name).read_text().splitlines()
        if len(lines) < 2 or not lines[0].startswith("step\t"):
            raise SystemExit(f"restart into new output lacks {name} header")
    run(restart_config, checkpoint_path)

    reference = json.loads((continuous_output / "checkpoint.json").read_text())
    restarted = json.loads(checkpoint_path.read_text())
    for key in (
        "time_s",
        "step_number",
        "cumulative_mole_ledger",
        "cumulative_volume_ledger",
        "next_dt_s",
        "consecutive_fast_steps",
    ):
        a, b = reference[key], restarted[key]
        if isinstance(a, float):
            if not close(a, b):
                raise SystemExit(f"restart metadata differs for {key}: {a} != {b}")
        elif a != b:
            raise SystemExit(f"restart metadata differs for {key}: {a} != {b}")
    for row, (a, b) in enumerate(zip(reference["state"], restarted["state"]), 1):
        if len(a) != len(b):
            raise SystemExit(f"restart state width differs at row {row}")
        for column, (x, y) in enumerate(zip(a, b), 1):
            if isinstance(x, bool):
                same = x == y
            else:
                same = close(x, y)
            if not same:
                raise SystemExit(
                    f"restart state differs at row {row}, column {column}: {x} != {y}"
                )

    reference_rows = (continuous_output / "diagnostics.tsv").read_text().splitlines()
    restarted_rows = (restarted_output / "diagnostics.tsv").read_text().splitlines()
    if reference_rows != restarted_rows:
        raise SystemExit("restarted adaptive diagnostics differ from uninterrupted run")
    for name in ("segment_history.tsv", "balance_history.tsv"):
        if (continuous_output / name).read_bytes() != (restarted_output / name).read_bytes():
            raise SystemExit(f"restarted {name} differs from uninterrupted run")
    import csv
    def history(path):
        with path.open(newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))
    continuous_global, restarted_global = history(continuous_output / "global_history.tsv"), history(restarted_output / "global_history.tsv")
    timing = {"assembly_seconds", "csr_seconds", "linear_upload_seconds", "linear_setup_seconds",
              "linear_solve_seconds", "linear_download_seconds"}
    if len(continuous_global) != len(restarted_global):
        raise SystemExit("restart global history row count changed")
    for a, b in zip(continuous_global, restarted_global):
        if {k:v for k,v in a.items() if k not in timing} != {k:v for k,v in b.items() if k not in timing}:
            raise SystemExit("restart global cumulative/scientific history is discontinuous")

    manifest = json.loads((restarted_output / "run_manifest.json").read_text())
    if manifest.get("accepted_steps") != 5 or manifest.get("restart", {}).get("sha256") != restart_hash:
        raise SystemExit("restart manifest lacks accepted-step or input-checksum provenance")

print("PASS adaptive CLI checkpoint/restart equivalence")
