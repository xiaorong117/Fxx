#!/usr/bin/env python3
"""Bounded-concurrency, restart-safe runner for audited batch35 configs."""

from __future__ import annotations

import argparse
import asyncio
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


async def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=6)
    parser.add_argument("--threads-per-job", type=int, default=8)
    parser.add_argument("--only-run-id", action="append", default=[])
    args = parser.parse_args()
    index_path = args.index.resolve()
    binary = args.binary.resolve()
    index = json.loads(index_path.read_text())
    if (index.get("case_count") != len(index.get("cases", [])) or
            not index.get("cases") or args.jobs <= 0 or
            args.threads_per_job <= 0):
        raise SystemExit("invalid batch index or concurrency")
    selected = [row for row in index["cases"] if
                not args.only_run_id or row["run_id"] in args.only_run_id]
    if args.only_run_id and len(selected) != len(set(args.only_run_id)):
        raise SystemExit("one or more requested run IDs are absent from index")
    root = Path(index["output_root"]).resolve()
    root.mkdir(parents=True, exist_ok=True)
    driver_lock = (root / "batch.run.lock").open("w")
    try:
        fcntl.flock(driver_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        raise SystemExit("another batch scheduler already holds batch.run.lock") from error
    progress_path = root / "batch_progress.json"
    state = {
        "format": "bubble_batch35_scheduler_progress_v1",
        "status": "running",
        "binary": str(binary),
        "binary_sha256": sha256(binary),
        "index": str(index_path),
        "selected_count": len(selected),
        "jobs": args.jobs,
        "threads_per_job": args.threads_per_job,
        "running": {}, "completed": [], "failed": [],
    }

    def save() -> None:
        state["updated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        atomic_json(progress_path, state)

    semaphore = asyncio.Semaphore(args.jobs)
    stop_launching = asyncio.Event()

    async def run(row: dict) -> None:
        async with semaphore:
            if stop_launching.is_set():
                return
            run_id = row["run_id"]
            output = Path(row["output_directory"])
            output.mkdir(parents=True, exist_ok=True)
            case_lock = (output / "run.lock").open("w")
            try:
                fcntl.flock(case_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                state["failed"].append({
                    "run_id": run_id, "return_code": None,
                    "error": "run.lock is held by another process",
                    "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                })
                stop_launching.set()
                save()
                case_lock.close()
                return
            checkpoint = output / "checkpoint.json"
            command = [str(binary), "--config", row["config"]]
            restarted = checkpoint.is_file()
            if restarted:
                command += ["--restart", str(checkpoint)]
            log_path = output / "driver.log"
            log = log_path.open("a" if restarted else "w")
            env = os.environ.copy()
            env.update({
                "OMP_NUM_THREADS": str(args.threads_per_job),
                "MKL_NUM_THREADS": str(args.threads_per_job),
                "OPENBLAS_NUM_THREADS": "1",
            })
            process = await asyncio.create_subprocess_exec(
                *command, cwd=str(index_path.parent.parent.parent), env=env,
                stdout=log, stderr=asyncio.subprocess.STDOUT)
            state["running"][run_id] = {
                "pid": process.pid, "log": str(log_path), "restart": restarted,
                "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            }
            save()
            return_code = await process.wait()
            log.close()
            fcntl.flock(case_lock, fcntl.LOCK_UN)
            case_lock.close()
            state["running"].pop(run_id, None)
            result = {"run_id": run_id, "return_code": return_code,
                      "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
            if return_code == 0:
                state["completed"].append(result)
            else:
                state["failed"].append(result)
                stop_launching.set()
            save()

    save()
    await asyncio.gather(*(run(row) for row in selected))
    state["status"] = "failed" if state["failed"] else "passed"
    state["not_launched_after_failure"] = len(selected) - (
        len(state["completed"]) + len(state["failed"]))
    save()
    return 1 if state["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main()))
