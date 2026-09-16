#!/usr/bin/env python3
"""Temporarily hold queued case locks during a scheduler concurrency handoff."""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import time
from pathlib import Path


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--scheduler-progress", type=Path, required=True)
    parser.add_argument("--ready", type=Path, required=True)
    parser.add_argument("--release", type=Path, required=True)
    args = parser.parse_args()
    index = json.loads(args.index.resolve().read_text())
    progress = json.loads(args.scheduler_progress.resolve().read_text())
    running = set(progress.get("running", {}))
    handles = []
    locked = []
    for row in index["cases"]:
        if row["run_id"] in running:
            continue
        output = Path(row["output_directory"])
        output.mkdir(parents=True, exist_ok=True)
        handle = (output / "run.lock").open("w")
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f'already locked: {row["run_id"]}') from error
        handles.append(handle)
        locked.append(row["run_id"])
    atomic_json(args.ready, {
        "status": "ready", "pid": os.getpid(),
        "excluded_running_count": len(running),
        "locked_queued_count": len(locked),
        "excluded_running": sorted(running),
    })
    while not args.release.exists():
        time.sleep(1)
    for handle in handles:
        fcntl.flock(handle, fcntl.LOCK_UN)
        handle.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
