#!/usr/bin/env python3
"""A post-step output failure must preserve accepted progress in its manifest."""

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


binary = Path(sys.argv[1]).resolve()
source = Path(sys.argv[2]).resolve()
base = json.loads((source / "tests/data/small_cli.json").read_text())

with tempfile.TemporaryDirectory(prefix="bubble_failure_progress_") as temporary:
    root = Path(temporary)
    output = root / "output"
    output.mkdir()
    # The initial VTK remains writable, while the first accepted-step VTK
    # deliberately collides with a directory.
    (output / "state_000001.vtk").mkdir()
    for key in ("pore_file", "throat_file", "connect_audit_file"):
        base["io"][key] = str((source.parent / base["io"][key]).resolve())
    base["io"]["output_directory"] = str(output)
    base["time"]["t_end_s"] = 2.0e-4
    config = root / "config.json"
    config.write_text(json.dumps(base, indent=2) + "\n")
    completed = subprocess.run(
        [str(binary), "--config", str(config)], text=True, capture_output=True
    )
    if completed.returncode == 0:
        raise SystemExit("deliberate post-step output failure unexpectedly succeeded")
    manifest = json.loads((output / "run_manifest.failed.json").read_text())
    if manifest.get("accepted_steps") != 1 or not math.isclose(
        manifest.get("failure_time_s", -1.0), 1.0e-4, rel_tol=0.0, abs_tol=1.0e-18
    ):
        raise SystemExit("failure manifest discarded already accepted progress")
    if not manifest.get("run_argv") or manifest.get("exit_status") != 1:
        raise SystemExit("failure manifest lacks command/status provenance")

print("PASS failure manifest preserves accepted progress")
