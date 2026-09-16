#!/usr/bin/env python3
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
path = root / "bubble_dissolution_solver/output/full_network_smoke/run_manifest.failed.json"
doc = json.loads(path.read_text())
if doc.get("exit_status") != 1 or doc.get("accepted_steps") != 0:
    raise SystemExit("failure manifest has wrong status/step count")
if doc.get("error") != "pore 32123 volume must be finite and positive":
    raise SystemExit(f"unexpected M1 failure: {doc.get('error')}")
if (not doc.get("run_command") or not doc.get("run_argv") or
        not doc.get("started_utc") or not doc.get("ended_utc")):
    raise SystemExit("failure manifest lacks reproducibility metadata")
if len(doc.get("build", {}).get("project_source_tree_sha256", "")) != 64:
    raise SystemExit("failure manifest lacks project source-tree provenance")
expected = {
    str((root / name).resolve()): hashlib.sha256((root / name).read_bytes()).hexdigest()
    for name in ("PNM_new_pore.txt", "PNM_new_throat.txt", "PNM_new_pore_connect.txt",
                 "bubble_dissolution_solver/configs/bubble_solver.smoke.json")
}
if doc.get("inputs") != expected:
    raise SystemExit("failure manifest input hashes mismatch")
print("PASS full-network M1 failure manifest audit")
