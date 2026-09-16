#!/usr/bin/env python3
"""Verify unsafe numerical configurations fail before mesh loading."""

import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path


binary = Path(sys.argv[1]).resolve()
source = Path(sys.argv[2]).resolve()
base = json.loads((source / "tests/data/small_cli.json").read_text())

cases = [
    (("case", "gas_component"), "   ", "gas_component must be nonblank"),
    (("source_notes", "physics"), "", "source_notes.physics must be nonblank"),
    (("io", "pore_file"), "", "io paths must be nonempty"),
    (("physics", "contact_angle_rad"), 4.0, "contact_angle_rad"),
    (("time", "fast_newton_iterations"), -1, "iteration limits"),
    (("nonlinear", "armijo_c1"), 0.0, "armijo_c1"),
    (("nonlinear", "line_search_reduction"), 1.0, "line_search_reduction"),
    (("nonlinear", "line_search_min_lambda"), 0.0, "line_search_min_lambda"),
    (("linear", "backend"), "invalid_backend", "unsupported linear.backend"),
    (("linear", "relative_tolerance"), 0.0, "relative_tolerance"),
    (("linear", "ilut_drop_tolerance"), -1.0, "ilut_drop_tolerance"),
    (("linear", "ilut_fill_factor"), 0, "ilut_fill_factor"),
    (("active_set", "transfer_switch_tolerance_mol_m3"), -1.0,
     "transfer_switch_tolerance_mol_m3"),
    (("boundary", "inlet_flag"), 0, "distinct and nonzero"),
    (("geometry", "input_basis"), "", "input_basis must be explicit physical"),
    (("geometry", "input_area_to_m2"), 2.0,
     "area conversion must equal length conversion squared"),
    (("geometry", "input_volume_to_m3"), 2.0,
     "volume conversion must equal length conversion cubed"),
]

with tempfile.TemporaryDirectory(prefix="bubble_invalid_config_") as temporary:
    root = Path(temporary)
    for index, (keys, value, expected) in enumerate(cases):
        config = copy.deepcopy(base)
        config["io"]["output_directory"] = str(root / f"output_{index}")
        # If validation accidentally reaches M1, this missing path makes the
        # ordering failure unmistakable rather than exercising a valid mesh.
        config["io"]["pore_file"] = str(root / "must_not_be_read.txt")
        config[keys[0]][keys[1]] = value
        path = root / f"case_{index}.json"
        path.write_text(json.dumps(config) + "\n")
        completed = subprocess.run(
            [str(binary), "--config", str(path)], text=True, capture_output=True
        )
        combined = completed.stdout + completed.stderr
        if completed.returncode == 0 or expected not in combined:
            raise SystemExit(
                f"invalid config case {keys} was not rejected before mesh loading:\n{combined}"
            )

    batch_cases = [
        ({"enabled": False, "maximum_individual_free_gas_fraction": 0.0,
          "maximum_batch_free_gas_fraction_per_step": 1e-8},
         "batch retirement fractions"),
        ({"enabled": False, "maximum_individual_free_gas_fraction": 1e-6,
          "maximum_batch_free_gas_fraction_per_step": 1e-8},
         "aggregate limit must not be smaller"),
        ({"enabled": True, "maximum_individual_free_gas_fraction": 1e-10,
          "maximum_batch_free_gas_fraction_per_step": 1e-8},
         "requires time_step_control.enabled=true"),
    ]
    for offset, (batch, expected) in enumerate(batch_cases, start=len(cases)):
        config = copy.deepcopy(base)
        config["io"]["output_directory"] = str(root / f"output_{offset}")
        config["io"]["pore_file"] = str(root / "must_not_be_read.txt")
        config["active_set"]["batch_microbubble_retirement"] = batch
        path = root / f"case_{offset}.json"
        path.write_text(json.dumps(config) + "\n")
        completed = subprocess.run(
            [str(binary), "--config", str(path)], text=True, capture_output=True
        )
        combined = completed.stdout + completed.stderr
        if completed.returncode == 0 or expected not in combined:
            raise SystemExit(
                f"invalid batch-retirement config was not rejected before mesh loading:\n{combined}"
            )

print(f"PASS {len(cases) + len(batch_cases)} unsafe configurations rejected before mesh loading")
