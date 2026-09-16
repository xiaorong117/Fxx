#!/usr/bin/env python3
"""Clone selected production configs into short, non-overlapping t=0 pilots."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--config-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--run-id", action="append", required=True)
    parser.add_argument("--t-end-s", type=float, default=0.002)
    args = parser.parse_args()
    index = json.loads(args.index.resolve().read_text())
    selected = {row["run_id"]: row for row in index["cases"]
                if row["run_id"] in args.run_id}
    if len(selected) != len(set(args.run_id)):
        raise SystemExit("one or more pilot run IDs are absent")
    config_root = args.config_root.resolve()
    if config_root.exists():
        raise SystemExit(f"refusing to overwrite {config_root}")
    config_root.mkdir(parents=True)
    rows = []
    for run_id in args.run_id:
        source = Path(selected[run_id]["config"])
        config = json.loads(source.read_text())
        pilot_id = run_id + "_pilot"
        output = args.output_root.resolve() / pilot_id
        config["case"]["name"] = pilot_id
        config["case"]["kind"] = "numerical_smoke_only"
        config["case"]["notes"] += " Short t=0 acceptance pilot; not production data."
        config["io"]["output_directory"] = str(output)
        config["io"]["checkpoint_interval_steps"] = 1
        config["time"]["t_end_s"] = args.t_end_s
        config["termination"]["maximum_physical_time_s"] = args.t_end_s
        config["termination"]["maximum_output_bytes"] = 1024 ** 3
        path = config_root / (pilot_id + ".json")
        path.write_text(json.dumps(config, indent=2) + "\n")
        rows.append({"run_id": pilot_id, "source_run_id": run_id,
                     "config": str(path), "output_directory": str(output)})
    pilot_index = {
        "format": "bubble_batch35_pilot_configs_v1", "status": "prepared",
        "case_count": len(rows), "output_root": str(args.output_root.resolve()),
        "cases": rows,
    }
    path = config_root / "pilot_cases.json"
    path.write_text(json.dumps(pilot_index, indent=2) + "\n")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
