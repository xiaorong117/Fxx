#!/usr/bin/env python3
"""Restricted control operations for production bubble-solver cases."""
from __future__ import annotations

import hashlib
import fcntl
import json
import os
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

PROJECT = Path(__file__).resolve().parents[1]
PRODUCTION_ROOT = (PROJECT / "output" / "production").resolve()
RUNNER = (PROJECT / "scripts" / "control_runner.py").resolve()
CASE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")


class ControlError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise ControlError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ControlError(f"JSON root is not an object: {path}")
    return value


def resolve_case(name: str, root: Path = PRODUCTION_ROOT) -> Path:
    root = root.resolve()
    if not CASE_NAME.fullmatch(name or ""):
        raise ControlError("invalid case identifier")
    case = (root / name).resolve()
    try:
        relative = case.relative_to(root)
    except ValueError as error:
        raise ControlError("case path escapes production root") from error
    if len(relative.parts) != 1 or not case.is_dir():
        raise ControlError("case is not a direct production child")
    return case


def solver_processes(case: Path) -> list[int]:
    completed = subprocess.run(["ps", "-eo", "pid=,args="], check=False,
                               capture_output=True, text=True, timeout=3)
    needles = (str(case),)
    result = []
    for line in completed.stdout.splitlines():
        match = re.match(r"\s*(\d+)\s+(.*)$", line)
        if not match or "bubble_solver" not in match.group(2):
            continue
        if any(needle in match.group(2) for needle in needles):
            result.append(int(match.group(1)))
    return result


def run_lock_held(case: Path) -> bool:
    fd = os.open(case / "run.lock", os.O_RDWR | os.O_CREAT, 0o600)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(fd, fcntl.LOCK_UN)
        return False
    finally:
        os.close(fd)


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def timestamps() -> tuple[str, str]:
    now = datetime.now(timezone.utc)
    return (now.isoformat().replace("+00:00", "Z"),
            now.astimezone(ZoneInfo("Asia/Shanghai")).isoformat())


def audit(case: Path, action: str, result: str,
          details: dict[str, Any] | None = None) -> None:
    utc, shanghai = timestamps()
    record = {"time_utc": utc, "time_shanghai": shanghai,
              "action": action, "case_name": case.name, "result": result,
              "details": details or {}}
    control = case / "control"
    control.mkdir(parents=True, exist_ok=True)
    with (control / "control_history.jsonl").open("a") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def configured_output(config: dict[str, Any]) -> Path:
    value = config.get("io", {}).get("output_directory")
    if not value:
        raise ControlError("configuration has no output_directory")
    path = Path(str(value))
    return path.resolve() if path.is_absolute() else (PROJECT.parent / path).resolve()


def find_config(case: Path) -> Path:
    state_path = case / "control" / "active_config.json"
    if state_path.is_file():
        candidate = Path(read_json(state_path).get("path", "")).resolve()
        if candidate.is_file() and configured_output(read_json(candidate)) == case:
            return candidate
    manifest_path = case / "run_manifest.json"
    if manifest_path.is_file():
        manifest = read_json(manifest_path)
        argv = manifest.get("run_argv", [])
        if isinstance(argv, list):
            for index, item in enumerate(argv[:-1]):
                if item == "--config":
                    candidate = Path(str(argv[index + 1])).resolve()
                    if candidate.is_file() and configured_output(read_json(candidate)) == case:
                        return candidate
    for candidate in sorted((PROJECT / "configs").glob("*.json")):
        try:
            config = read_json(candidate)
            if configured_output(config) == case:
                return candidate.resolve()
        except ControlError:
            continue
    raise ControlError("no validated configuration found for case")


def stage_settings(name: str, changes: dict[str, Any],
                   root: Path = PRODUCTION_ROOT) -> dict[str, Any]:
    case = resolve_case(name, root)
    allowed = {"dt_max_s", "checkpoint_interval_steps", "vtk_interval_steps"}
    if not changes or set(changes) - allowed:
        raise ControlError("only dt_max_s, checkpoint_interval_steps and vtk_interval_steps may be staged")
    source = find_config(case)
    config = read_json(source)
    before_hash = sha256(source)
    applied: dict[str, Any] = {}
    if "dt_max_s" in changes:
        value = float(changes["dt_max_s"])
        minimum = float(config.get("time_step_control", {}).get(
            "dt_min_s", config.get("time", {}).get("dt_min_s", 0)))
        if not (value >= minimum and value <= 1.0e6):
            raise ControlError("dt_max_s is outside the validated range")
        config["time"]["dt_max_s"] = value
        config.setdefault("time_step_control", {})["dt_max_s"] = value
        applied["dt_max_s"] = value
    for key in ("checkpoint_interval_steps", "vtk_interval_steps"):
        if key not in changes:
            continue
        value = int(changes[key])
        if value < 1 or value > 1_000_000_000:
            raise ControlError(f"{key} is outside the validated range")
        config["io"][key] = value
        applied[key] = value
    control = case / "control"
    control.mkdir(parents=True, exist_ok=True)
    existing = list(control.glob("config_revision_*.json"))
    revision = 1 + max([int(p.stem.rsplit("_", 1)[1]) for p in existing] or [0])
    output = control / f"config_revision_{revision:04d}.json"
    atomic_json(output, config)
    metadata = {"revision": revision, "path": str(output),
                "sha256": sha256(output), "source_path": str(source),
                "source_sha256": before_hash, "changes": applied,
                "effective": "next_restart"}
    atomic_json(control / "staged_config.json", metadata)
    audit(case, "stage_settings", "success", metadata)
    return metadata


def request_stop(name: str, root: Path = PRODUCTION_ROOT) -> dict[str, Any]:
    case = resolve_case(name, root)
    pids = solver_processes(case)
    if not run_lock_held(case):
        raise ControlError("no running solver process found")
    stop = case / "STOP"
    stop.write_text("requested by case dashboard\n")
    details = {"pids": pids, "stop_path": str(stop),
               "checkpoint_sha256_before": sha256(case / "checkpoint.json")
                   if (case / "checkpoint.json").is_file() else None}
    audit(case, "safe_stop", "requested", details)
    return details


def choose_restart_config(case: Path) -> tuple[Path, dict[str, Any] | None]:
    staged_path = case / "control" / "staged_config.json"
    if staged_path.is_file():
        staged = read_json(staged_path)
        candidate = Path(str(staged.get("path", ""))).resolve()
        if not candidate.is_file() or sha256(candidate) != staged.get("sha256"):
            raise ControlError("staged configuration hash mismatch")
        if configured_output(read_json(candidate)) != case:
            raise ControlError("staged configuration output mismatch")
        return candidate, staged
    return find_config(case), None


def restart_case(name: str, solver_binary: Path,
                 root: Path = PRODUCTION_ROOT) -> dict[str, Any]:
    case = resolve_case(name, root)
    if run_lock_held(case):
        raise ControlError("solver is already running")
    checkpoint = case / "checkpoint.json"
    if not checkpoint.is_file():
        raise ControlError("checkpoint.json is missing")
    checkpoint_data = read_json(checkpoint)
    if checkpoint_data.get("format") not in ("bubble_checkpoint_v5", "bubble_checkpoint_v6"):
        raise ControlError("unsupported checkpoint format")
    config, staged = choose_restart_config(case)
    solver_binary = solver_binary.resolve()
    if not solver_binary.is_file() or not os.access(solver_binary, os.X_OK):
        raise ControlError("fixed solver binary is missing or not executable")
    stop = case / "STOP"
    if stop.exists():
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        stop.rename(case / f"STOP.control_archive_{stamp}")
    session = "bubble_ctl_" + re.sub(r"[^A-Za-z0-9_-]", "_", name)[:80]
    if subprocess.run(["tmux", "has-session", "-t", session], check=False,
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0:
        raise ControlError("control tmux session already exists")
    command = shlex.join([sys.executable, str(RUNNER), "--binary", str(solver_binary),
                          "--config", str(config), "--checkpoint", str(checkpoint),
                          "--output", str(case), "--production-root", str(root.resolve())])
    completed = subprocess.run(["tmux", "new-session", "-d", "-s", session,
                                "-c", str(PROJECT.parent), command],
                               check=False, capture_output=True, text=True, timeout=5)
    if completed.returncode:
        raise ControlError("tmux start failed: " + completed.stderr.strip())
    metadata = {"tmux_session": session, "binary_path": str(solver_binary),
                "binary_sha256": sha256(solver_binary), "config_path": str(config),
                "config_sha256": sha256(config), "checkpoint_path": str(checkpoint),
                "checkpoint_sha256": sha256(checkpoint),
                "checkpoint_step": checkpoint_data.get("step_number"),
                "checkpoint_time_s": checkpoint_data.get("time_s"),
                "staged_revision": staged.get("revision") if staged else None}
    atomic_json(case / "control" / "active_config.json",
                {"path": str(config), "sha256": sha256(config),
                 "revision": metadata["staged_revision"]})
    staged_path = case / "control" / "staged_config.json"
    if staged_path.exists():
        staged_path.rename(case / "control" / "last_applied_config.json")
    audit(case, "restart", "started", metadata)
    return metadata


def control_state(name: str, root: Path = PRODUCTION_ROOT) -> dict[str, Any]:
    case = resolve_case(name, root)
    staged = read_json(case / "control" / "staged_config.json") \
        if (case / "control" / "staged_config.json").is_file() else None
    return {"running": run_lock_held(case),
            "checkpoint_exists": (case / "checkpoint.json").is_file(),
            "staged": staged}
