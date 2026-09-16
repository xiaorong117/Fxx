#!/usr/bin/env python3
"""Read-only case dashboard for the bubble dissolution production runs.

The service deliberately restricts all case reads to one resolved production
root.  It never reads environment variables, command history, or credentials.
Use ``--snapshot`` for a machine-readable audit without starting the server.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import hmac
import html
import json
import math
import os
import re
import subprocess
import threading
import time
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

import case_control

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = (PROJECT / "output" / "production").resolve()


def sha256(path: Path) -> str | None:
    try:
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1 << 20), b""):
                h.update(block)
        return h.hexdigest()
    except (OSError, ValueError):
        return None


@lru_cache(maxsize=512)
def _sha256_immutable(path_text: str, size: int, modified_ns: int) -> str | None:
    del size, modified_ns
    return sha256(Path(path_text))


def sha256_immutable(path: Path) -> str | None:
    try:
        stat = path.stat()
        return _sha256_immutable(str(path.resolve()), stat.st_size, stat.st_mtime_ns)
    except OSError:
        return None


def embedded_source_hash(executable: Path) -> str | None:
    try:
        hashes = {x.decode() for x in re.findall(
            rb"(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])", executable.read_bytes())}
    except OSError:
        return None
    dependencies = {
        "d88119e194e6ba96d5fdfd22212e3e101a5a95038b27c2c967bb8051927bbe2e",
        "dc7f7d9a31ec7dc1824a55cefe96feaa25db1be46f0a8be84c6c870fc7628afe",
        "57821f5579438289c048ab699d6bec173badf05319c212e464f6612fa14634e9",
    }
    candidates = sorted(hashes - dependencies)
    return candidates[0] if len(candidates) == 1 else None


@lru_cache(maxsize=16)
def virtual_counts(path_text: str, inlet_flag: int, outlet_flag: int) -> tuple[int, int]:
    inlet = outlet = 0
    try:
        for line in Path(path_text).read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 11 and int(fields[1]) == 2:
                flag = int(fields[10])
                inlet += flag == inlet_flag; outlet += flag == outlet_flag
    except (OSError, ValueError):
        pass
    return inlet, outlet


@lru_cache(maxsize=16)
def pore_segment_map(path_text: str, segments: int) -> dict[int, int]:
    real = []
    try:
        for line in Path(path_text).read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 10 and int(fields[1]) != 2 and float(fields[9]) > 0:
                real.append((int(fields[0]), float(fields[3])))
    except (OSError, ValueError):
        return {}
    if not real or segments <= 0:
        return {}
    xmin = min(x for _, x in real); xmax = max(x for _, x in real)
    span = xmax - xmin
    return {node: 1 + min(segments - 1, max(0, int(((x - xmin) / span if span else 0.0) * segments)))
            for node, x in real}


@lru_cache(maxsize=16)
def initial_segment_saturation(path_text: str, segments: int) -> dict[int, float]:
    mapping = pore_segment_map(path_text, segments)
    gas = {sid: 0.0 for sid in range(1, segments + 1)}
    volume = {sid: 0.0 for sid in range(1, segments + 1)}
    try:
        for line in Path(path_text).read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            node = int(fields[0]); sid = mapping.get(node)
            if sid is None or len(fields) < 12:
                continue
            v = float(fields[9]); sw = float(fields[11])
            volume[sid] += v
            gas[sid] += v * (1.0 - sw)
    except (OSError, ValueError, IndexError):
        return {}
    return {sid: gas[sid] / volume[sid] if volume[sid] > 0 else 0.0
            for sid in range(1, segments + 1)}


@lru_cache(maxsize=16)
def initial_segment_inventory(path_text: str, segments: int,
                              sg_off: float, volume_scale: float) -> dict[int, dict[str, float]]:
    mapping = pore_segment_map(path_text, segments)
    result = {sid: {"initial_segment_bubble_count": 0,
                    "initial_segment_gas_volume_m3": 0.0}
              for sid in range(1, segments + 1)}
    try:
        for line in Path(path_text).read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split(); node = int(fields[0]); sid = mapping.get(node)
            if sid is None or len(fields) < 12:
                continue
            volume = float(fields[9]) * volume_scale
            sg = 1.0 - float(fields[11])
            result[sid]["initial_segment_gas_volume_m3"] += volume * sg
            if sg > sg_off:
                result[sid]["initial_segment_bubble_count"] += 1
    except (OSError, ValueError, IndexError):
        return {}
    return result


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text())
        return value if isinstance(value, dict) else None
    except (OSError, ValueError, UnicodeError):
        return None


@lru_cache(maxsize=512)
def _read_json_immutable(path_text: str, size: int,
                         modified_ns: int) -> dict[str, Any] | None:
    del size, modified_ns
    return read_json(Path(path_text))


def read_config(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    try:
        stat = path.stat()
        return _read_json_immutable(str(path.resolve()), stat.st_size, stat.st_mtime_ns)
    except OSError:
        return None


def last_row(path: Path) -> dict[str, str] | None:
    rows = tail_rows(path, 1)
    return rows[-1] if rows else None


def tail_rows(path: Path, count: int = 10,
              byte_limit: int = 1 << 20) -> list[dict[str, str]]:
    """Read complete trailing TSV rows without loading an unbounded history."""
    try:
        with path.open("rb") as stream:
            header = stream.readline().decode("utf-8").rstrip("\r\n")
            stream.seek(0, os.SEEK_END)
            end = stream.tell()
            start = max(0, end - byte_limit)
            stream.seek(start)
            text = stream.read().decode("utf-8", errors="strict")
        lines = text.splitlines()
        if start and lines:
            lines = lines[1:]
        lines = [line for line in lines if line.strip()][-count:]
        if not header or not lines:
            return []
        return list(csv.DictReader([header, *lines], delimiter="\t"))
    except (OSError, csv.Error, UnicodeError):
        return []


def history_series(path: Path, limit: int = 360) -> list[dict[str, float]]:
    """Read a compact, downsampled set of global-history observables."""
    wanted = {
        "time_s", "free_gas_mol", "dissolved_gas_mol", "total_component_mol",
        "main_real_gas_saturation", "mean_Sg", "dt_s", "newton_iterations",
        "linear_iterations", "residual_scaled_inf", "relative_update_norm",
        "epsilon_N_adjusted", "epsilon_V_adjusted", "active_bubble_count",
    }
    try:
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
    except (OSError, csv.Error, UnicodeError):
        return []
    if len(rows) > limit:
        stride = max(1, (len(rows) - 1) // (limit - 1))
        selected = rows[::stride]
        if selected[-1] is not rows[-1]:
            selected.append(rows[-1])
    else:
        selected = rows
    result = []
    for row in selected:
        item: dict[str, float] = {}
        for key in wanted:
            try:
                item[key] = float(row[key])
            except (KeyError, TypeError, ValueError):
                continue
        if "time_s" in item:
            result.append(item)
    return result


def conservation_series(path: Path, limit: int = 360) -> list[dict[str, float]]:
    """Read solver-accepted conservation errors from diagnostics.tsv."""
    wanted = ("time_s", "epsilon_N_adjusted", "epsilon_V_adjusted",
              "residual_scaled_inf", "relative_update_norm")
    try:
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
    except (OSError, csv.Error, UnicodeError):
        return []
    if len(rows) > limit:
        stride = max(1, (len(rows) - 1) // (limit - 1))
        selected = rows[::stride]
        if selected[-1] is not rows[-1]:
            selected.append(rows[-1])
    else:
        selected = rows
    result = []
    for row in selected:
        item = {}
        for key in wanted:
            try:
                value = float(row[key])
                if math.isfinite(value):
                    item[key] = value
            except (KeyError, TypeError, ValueError):
                continue
        if "time_s" in item:
            result.append(item)
    return result


def latest_segment_profile(path: Path) -> list[dict[str, float]]:
    try:
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
    except (OSError, csv.Error, UnicodeError):
        return []
    if not rows:
        return []
    latest_step = rows[-1].get("step")
    selected = [row for row in rows if row.get("step") == latest_step]
    result = []
    for row in selected:
        try:
            result.append({
                "segment_id": float(row["segment_id"]),
                "segment_mean_Sg": float(row["segment_mean_Sg"]),
                "segment_gas_volume_m3": float(row["segment_gas_volume_m3"]),
                "segment_active_bubble_count": float(row["segment_active_bubble_count"]),
            })
        except (KeyError, TypeError, ValueError):
            continue
    return sorted(result, key=lambda x: x["segment_id"])


def disappearance_events(case_dir: Path) -> list[dict[str, Any]]:
    event_path = case_dir / "gas_disappearance_events.tsv"
    if event_path.is_file():
        try:
            with event_path.open(newline="") as stream:
                return [{"time_s": float(row["estimated_crossing_time_s"]),
                         "accepted_state_time_s": float(row["accepted_state_time_s"]),
                         "event_node_id": float(row["event_node_id"]),
                         "event_dt_s": float(row["accepted_dt_s"]),
                         "event_mode": row.get("event_mode", "unknown")}
                        for row in csv.DictReader(stream, delimiter="\t")]
        except (OSError, csv.Error, ValueError, TypeError, KeyError):
            return []
    path = case_dir / "diagnostics.tsv"
    try:
        with path.open(newline="") as stream:
            rows = csv.DictReader(stream, delimiter="\t")
            result = []
            for row in rows:
                if str(row.get("gas_disappearance_event", "0")) != "1":
                    continue
                result.append({"time_s": float(row["time_s"]),
                               "event_node_id": float(row.get("event_node_id", 0)),
                               "event_dt_s": float(row.get("event_dt_s", 0))})
            return result
    except (OSError, csv.Error, ValueError, TypeError):
        return []


def bubble_process_table() -> list[dict[str, Any]]:
    """Inspect all bubble_solver processes once per dashboard snapshot."""
    found: list[dict[str, Any]] = []
    try:
        completed = subprocess.run(
            ["ps", "-eo", "pid=,etimes=,lstart=,args="],
            check=False, capture_output=True, text=True, timeout=2,
        )
    except (OSError, subprocess.SubprocessError):
        return found
    for line in completed.stdout.splitlines():
        match = re.match(r"\s*(\d+)\s+(\d+)\s+(\w+\s+\w+\s+\d+\s+[\d:]+\s+\d{4})\s+(.*)$", line)
        if not match or "case_dashboard" in match.group(4):
            continue
        pid = int(match.group(1))
        cmd = match.group(4)
        if "bubble_solver" not in cmd:
            continue
        exe = None
        try:
            exe = os.readlink(f"/proc/{pid}/exe")
        except OSError:
            pass
        found.append({"pid": pid, "elapsed_seconds": int(match.group(2)),
                      "started_local": match.group(3), "argv": cmd, "executable": exe})
    return found


def process_snapshot(output: Path, config_path: Path | None = None,
                     table: list[dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    """Match an output/config against one already-bounded process inventory."""
    needles = [str(output)]
    if config_path:
        needles.append(str(config_path))
    return [item for item in (table if table is not None else bubble_process_table())
            if any(needle in item.get("argv", "") for needle in needles)]


def find_config(case_dir: Path, manifest: dict[str, Any] | None) -> Path | None:
    argv = manifest.get("run_argv", []) if manifest else []
    if isinstance(argv, list):
        for index, value in enumerate(argv[:-1]):
            if value == "--config":
                candidate = Path(str(argv[index + 1])).resolve()
                if candidate.is_file():
                    return candidate
    configs = PROJECT / "configs"
    for candidate in sorted(configs.rglob("*.json"), reverse=True):
        value = read_json(candidate)
        if not value:
            continue
        case = value.get("case", {})
        output = value.get("io", {}).get("output_directory")
        if output and (PROJECT.parent / str(output)).resolve() == case_dir.resolve():
            return candidate.resolve()
        if case.get("name") == (manifest or {}).get("case_name") and case_dir.name == case.get("name"):
            return candidate.resolve()
    return None


def safe_child(root: Path, item: Path) -> bool:
    try:
        item.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def configured_output(config: dict[str, Any]) -> Path | None:
    value = config.get("io", {}).get("output_directory")
    if not value:
        return None
    path = Path(str(value))
    return path.resolve() if path.is_absolute() else (PROJECT.parent / path).resolve()


def config_catalog(root: Path) -> dict[Path, Path]:
    """Return one newest config for every output confined below production root."""
    result: dict[Path, Path] = {}
    for candidate in (PROJECT / "configs").rglob("*.json"):
        config = read_config(candidate)
        output = configured_output(config or {})
        if output is None or not safe_child(root, output):
            continue
        prior = result.get(output)
        if prior is None or (candidate.stat().st_mtime, str(candidate)) > (
                prior.stat().st_mtime, str(prior)):
            result[output] = candidate.resolve()
    return result


def scheduler_statuses(root: Path) -> tuple[dict[str, str], list[dict[str, Any]]]:
    statuses: dict[str, str] = {}
    schedulers = []
    for path in root.rglob("batch_progress.json"):
        value = read_json(path)
        if not value or not safe_child(root, path):
            continue
        index_path = Path(str(value.get("index", "")))
        index = read_json(index_path) if index_path.is_file() else None
        if value.get("status") == "running" and index:
            for item in index.get("cases", []):
                run_id = item.get("run_id", "")
                if run_id:
                    statuses.setdefault(run_id, "queued")
        for run_id in value.get("running", {}):
            statuses[run_id] = "running"
        for item in value.get("completed", []):
            statuses[item.get("run_id", "")] = "completed"
        for item in value.get("failed", []):
            statuses[item.get("run_id", "")] = "failed"
        schedulers.append({
            "path": str(path), "status": value.get("status"),
            "selected_count": value.get("selected_count"),
            "running": len(value.get("running", {})),
            "completed": len(value.get("completed", [])),
            "failed": len(value.get("failed", [])),
            "jobs": value.get("jobs"),
            "updated_utc": value.get("updated_utc"),
        })
    return statuses, schedulers


def console_tail(root: Path, case_name: str, line_limit: int = 80,
                 byte_limit: int = 131072) -> dict[str, Any]:
    """Return a bounded log tail for one production descendant only."""
    parts = Path(case_name or "").parts
    if (not parts or Path(case_name).is_absolute() or
            any(not re.fullmatch(r"[A-Za-z0-9_.-]{1,160}", part)
                for part in parts)):
        raise ValueError("invalid case name")
    resolved_root = root.resolve()
    case_dir = (resolved_root / case_name).resolve()
    if not safe_child(resolved_root, case_dir) or not case_dir.is_dir():
        raise ValueError("case is outside the production root or missing")
    line_limit = max(10, min(int(line_limit), 200))
    log_path = case_dir / "run.log"
    if not log_path.is_file():
        log_path = case_dir / "driver.log"
    try:
        size = log_path.stat().st_size
        with log_path.open("rb") as stream:
            offset = max(0, size - byte_limit)
            stream.seek(offset)
            payload = stream.read(byte_limit)
        text = payload.decode("utf-8", errors="replace")
        lines = text.splitlines()
        if offset and lines:
            lines = lines[1:]
        lines = lines[-line_limit:]
        modified = log_path.stat().st_mtime
        result = {
            "case_name": case_name,
            "source": log_path.name,
            "size_bytes": size,
            "line_limit": line_limit,
            "lines": lines,
            "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                         time.gmtime(modified)),
            "age_seconds": max(0.0, time.time() - modified),
            "available": True,
        }
        session_file = case_dir / "dashboard_tmux_session.txt"
        try:
            session = session_file.read_text().strip()
            if not re.fullmatch(r"[A-Za-z0-9_.-]{1,128}", session):
                raise ValueError("invalid tmux session metadata")
            target = f"{session}:live.0"
            pane_command = subprocess.run(
                ["tmux", "display-message", "-p", "-t", target,
                 "#{pane_start_command}"], check=True, capture_output=True,
                text=True, timeout=2).stdout.strip()
            if "tail " not in pane_command or str(log_path) not in pane_command:
                raise ValueError("tmux live pane is not the approved run.log tail")
            captured = subprocess.run(
                ["tmux", "capture-pane", "-p", "-J", "-S", "-200", "-t", target],
                check=True, capture_output=True, text=True, timeout=2).stdout.splitlines()
            while captured and not captured[-1].strip():
                captured.pop()
            result["lines"] = captured[-line_limit:]
            result["source"] = f"tmux:{session}:live"
            result["tmux_session"] = session
            result["tmux_window"] = "live"
        except (OSError, ValueError, subprocess.SubprocessError):
            pass
        return result
    except OSError:
        return {"case_name": case_name, "source": log_path.name, "lines": [],
                "available": False, "line_limit": line_limit}


def warnings_for(case: dict[str, Any]) -> list[str]:
    warnings: list[str] = []
    if case["status"] == "running" and not case["processes"]:
        warnings.append("进程不存在但状态仍为 running")
    if case["manifest"] is None and case["status"] != "running":
        warnings.append("manifest 缺失")
    if case["checkpoint"] is None:
        warnings.append("checkpoint 缺失或不可读")
    if case["termination_events"] is None and case["status"] != "running":
        warnings.append("termination_events 缺失")
    row = case.get("last_diagnostics") or {}
    try:
        dt = float(row.get("dt_s", "nan")); dt_min = float(case["effective"].get("dt_min_s", "nan"))
        if dt_min > 0 and dt <= 10 * dt_min:
            warnings.append("dt 接近 dt_min")
        retries = int(float(row.get("timestep_retries", "0")))
        if retries > 0:
            warnings.append("最近时间步发生 retries")
        hard = int(case["effective"].get("hard_newton_iterations") or 12)
        if int(float(row.get("newton_iterations", "0"))) > hard:
            warnings.append("Newton 迭代次数过高")
        if float(row.get("epsilon_N_adjusted", "0")) > float(case["effective"].get("mass_balance_hard_limit", "inf")):
            warnings.append("摩尔守恒误差超过容差")
        if float(row.get("epsilon_V_adjusted", "0")) > float(case["effective"].get("volume_balance_hard_limit", "inf")):
            warnings.append("液体体积守恒误差超过容差")
    except (ValueError, TypeError):
        warnings.append("diagnostics 数值字段无法解析")
    history = case.get("recent_diagnostics") or []
    try:
        if len(history) >= 5 and all(float(x.get("dt_s", 0)) <= 1.01 * float(history[0].get("dt_s", 0)) for x in history[-5:]):
            warnings.append("最近5个接受步 dt 未增长")
        if len(history) >= 3 and all(int(float(x.get("timestep_retries", 0))) > 0 for x in history[-3:]):
            warnings.append("retries 连续发生")
        free = [float(x["total_free_gas_moles"]) for x in history if x.get("total_free_gas_moles")]
        if any(b > a * (1 + 1e-10) for a, b in zip(free, free[1:])):
            warnings.append("自由气体摩尔数增加")
        total = [float(x["total_gas_moles"]) for x in history if x.get("total_gas_moles")]
        if total and (max(total) - min(total)) / max(abs(total[0]), 1e-300) > 1e-3:
            warnings.append("total gas 库存变化超过0.1%；请结合边界通量核对")
    except (ValueError, TypeError, KeyError):
        pass
    if case["status"] == "running" and case.get("diagnostics_age_seconds", 0) > 600:
        warnings.append("diagnostics 超过10分钟未更新")
    try:
        interval = int(case["effective"].get("checkpoint_interval_steps") or 1)
        if case.get("checkpoint") and row and int(case["checkpoint"].get("step_number", 0)) + 2 * interval < int(row.get("step", 0)):
            warnings.append("checkpoint 落后于 diagnostics")
    except (ValueError, TypeError):
        pass
    if case["kind"] == "accelerated_preview" and case["scientific_interpretation"] == "scientific baseline":
        warnings.append("kLa×100 预览被标记为 scientific baseline")
    if case.get("config_sha256") and case.get("manifest_config_sha256") and case["config_sha256"] != case["manifest_config_sha256"]:
        warnings.append("配置哈希与启动时记录不一致")
    if (case.get("executable_sha256") and case.get("executable_path_current_sha256") and
            case["executable_sha256"] != case["executable_path_current_sha256"]):
        warnings.append("运行中二进制映像与当前可执行文件哈希不一致")
    return warnings


def snapshot_case(case_dir: Path, root: Path,
                  config_path_override: Path | None = None,
                  process_table: list[dict[str, Any]] | None = None,
                  scheduler_status: str | None = None) -> dict[str, Any]:
    manifest = read_json(case_dir / "run_manifest.json")
    checkpoint = read_json(case_dir / "checkpoint.json")
    config_path = config_path_override or find_config(case_dir, manifest)
    config = read_config(config_path)
    row = last_row(case_dir / "diagnostics.tsv")
    recent = tail_rows(case_dir / "diagnostics.tsv", 10)
    run_values: dict[str, Any] = {}
    try:
        for line in (case_dir / "run.log").read_text(errors="replace").splitlines()[:8]:
            if line.startswith("started_utc="):
                run_values["started_utc"] = line.split("=", 1)[1]
            elif line.startswith("config="):
                run_values["config_path"] = line.split("=", 1)[1]
            elif line.startswith("nodes="):
                for key, value in re.findall(r"(nodes|edges|virtual_nodes|boundary_edges)=(\d+)", line):
                    run_values[key] = int(value)
    except OSError:
        pass
    config_hash = sha256(config_path) if config_path else None
    manifest_inputs = (manifest or {}).get("inputs", {})
    manifest_config_hash = None
    if config_path:
        for key, value in manifest_inputs.items():
            if Path(str(key)).name == config_path.name and isinstance(value, str):
                manifest_config_hash = value
    term = (manifest or {}).get("termination", {})
    physics = (config or {}).get("physics", {})
    time_cfg = (config or {}).get("time_step_control", {})
    time_legacy = (config or {}).get("time", {})
    batch_cfg = (config or {}).get("active_set", {}).get(
        "batch_microbubble_retirement", {})
    effective = {
        "backend": (manifest or {}).get("linear_backend", {}).get("selected", (config or {}).get("linear", {}).get("backend")),
        "pardiso_threads": (manifest or {}).get("linear_backend", {}).get(
            "pardiso_threads", (config or {}).get("linear", {}).get("pardiso_threads")),
        "dt_initial_s": time_legacy.get("dt_initial_s"),
        "dt_min_s": time_cfg.get("dt_min_s", time_legacy.get("dt_min_s")),
        "dt_max_s": time_cfg.get("dt_max_s", time_legacy.get("dt_max_s")),
        "growth_factor": time_cfg.get("growth_factor", time_legacy.get("growth_factor")),
        "shrink_factor": time_cfg.get("shrink_factor", time_legacy.get("cut_factor")),
        "adaptive_dt_enabled": time_cfg.get("enabled", False),
        "residual_tolerance": (config or {}).get("nonlinear", {}).get("residual_tolerance"),
        "update_tolerance": (config or {}).get("nonlinear", {}).get("update_tolerance"),
        "max_newton_iterations": (config or {}).get("nonlinear", {}).get("max_iterations"),
        "max_linear_iterations": (config or {}).get("linear", {}).get("max_iterations"),
        "mass_balance_hard_limit": (config or {}).get("acceptance", {}).get("mass_balance_hard_limit"),
        "volume_balance_hard_limit": (config or {}).get("acceptance", {}).get("volume_balance_hard_limit"),
        "checkpoint_interval_steps": (config or {}).get("io", {}).get("checkpoint_interval_steps"),
        "vtk_interval_steps": (config or {}).get("io", {}).get("vtk_interval_steps"),
        "residual_free_gas_fraction": term.get("residual_free_gas_fraction", (config or {}).get("termination", {}).get("residual_free_gas_fraction")),
        "batch_microbubble_retirement_enabled": batch_cfg.get("enabled", False),
        "batch_individual_free_gas_fraction": batch_cfg.get(
            "maximum_individual_free_gas_fraction", 1e-10),
        "batch_total_free_gas_fraction_per_step": batch_cfg.get(
            "maximum_batch_free_gas_fraction_per_step", 1e-8),
    }
    if checkpoint:
        effective["current_dt_s"] = checkpoint.get("next_dt_s")
    if row:
        effective["current_dt_s"] = row.get("dt_s")
    effective.update({
        "hard_newton_iterations": time_cfg.get("hard_newton_iterations"),
        "easy_newton_iterations": time_cfg.get("easy_newton_iterations"),
        "normal_newton_iterations": time_cfg.get("normal_newton_iterations"),
        "concentration_cap_enabled": time_cfg.get("concentration_change_limit_enabled"),
        "gas_timescale_fraction": time_cfg.get("gas_timescale_fraction"),
        "maximum_physical_time_s": (config or {}).get("termination", {}).get("maximum_physical_time_s", time_legacy.get("t_end_s")),
        "maximum_accepted_steps": (config or {}).get("termination", {}).get("maximum_accepted_steps"),
    })
    output_path = case_dir.resolve()
    configured_output = (config or {}).get("io", {}).get("output_directory")
    config_matches_directory = bool(configured_output) and (
        (PROJECT.parent / str(configured_output)).resolve() == output_path)
    processes = process_snapshot(output_path,
                                 config_path if config_matches_directory else None,
                                 process_table)
    status = "running" if processes else (
        "completed" if manifest and manifest.get("exit_status") == 0 else
        "failed" if ((manifest and manifest.get("exit_status", 0) != 0) or
                     (case_dir / "run_manifest.failed.json").is_file()) else
        scheduler_status or ("interrupted" if case_dir.exists() else "queued"))
    case_name = (manifest or {}).get("case_name") or (config or {}).get("case", {}).get("name") or case_dir.name
    multiplier = physics.get("mass_transfer_multiplier", (manifest or {}).get("mass_transfer", {}).get("mass_transfer_multiplier", 1.0))
    kind = "accelerated_preview" if float(multiplier or 1) != 1.0 or "accelerated" in case_name.lower() else "scientific baseline"
    manifest_view = None
    if manifest:
        manifest_view = {
            key: manifest.get(key) for key in (
                "case_name", "case_kind", "started_utc", "ended_utc", "termination_reason",
                "accepted_steps", "final_time_s", "elapsed_seconds", "exit_status", "build",
                "linear_backend", "termination", "network", "mass_transfer", "run_argv",
            ) if key in manifest
        }
        if "inputs" in manifest:
            manifest_view["inputs"] = manifest["inputs"]
        if "configuration" in manifest:
            manifest_view["configuration"] = manifest["configuration"]
    checkpoint_view = None
    if checkpoint:
        checkpoint_view = {key: checkpoint.get(key) for key in (
            "format", "time_s", "step_number", "next_dt_s", "consecutive_fast_steps",
            "mass_transfer_multiplier", "gas_disappearance_events",
            "batch_microbubble_retirement_enabled",
            "batch_individual_free_gas_fraction",
            "batch_total_free_gas_fraction_per_step",
            "last_gas_disappearance_node_id", "last_gas_disappearance_time_s",
            "background_abs_Pa", "provenance",
        ) if key in checkpoint}
    network = (manifest or {}).get("network", {})
    geometry_audit = read_json(case_dir / "degenerate_geometry_audit.json")
    global_row = last_row(case_dir / "global_history.tsv") or {}
    series = history_series(case_dir / "global_history.tsv")
    conservation = conservation_series(case_dir / "diagnostics.tsv")
    segment_profile = latest_segment_profile(case_dir / "segment_history.tsv")
    events = disappearance_events(case_dir)
    segment_volume = 0.0
    try:
        with (case_dir / "segment_geometry.tsv").open(newline="") as stream:
            segment_volume = sum(float(x["pore_volume_m3"]) for x in csv.DictReader(stream, delimiter="\t"))
    except (OSError, ValueError, KeyError, csv.Error):
        pass
    pnm = []
    for key in ("pore_file", "throat_file", "connect_audit_file"):
        value = (config or {}).get("io", {}).get(key)
        if value:
            path = (PROJECT.parent / value).resolve() if not Path(value).is_absolute() else Path(value).resolve()
            pnm.append({"role": key, "path": str(path),
                        "sha256": sha256_immutable(path)})
    pore_path = next((Path(item["path"]) for item in pnm if item["role"] == "pore_file"), None)
    inlet_virtual = outlet_virtual = None
    if pore_path:
        inlet_virtual, outlet_virtual = virtual_counts(
            str(pore_path), int((config or {}).get("boundary", {}).get("inlet_flag", 1)),
            int((config or {}).get("boundary", {}).get("outlet_flag", 2)))
    segment_count = int((config or {}).get("analysis", {}).get("x_segments", 20))
    node_segments = pore_segment_map(str(pore_path), segment_count) if pore_path else {}
    initial_sg = initial_segment_saturation(str(pore_path), segment_count) if pore_path else {}
    initial_inventory = initial_segment_inventory(
        str(pore_path), segment_count,
        float((config or {}).get("active_set", {}).get("Sg_off", 1.0e-8)),
        float((config or {}).get("geometry", {}).get("input_volume_to_m3", 1.0))) if pore_path else {}
    event_summary = {sid: {"event_count": 0, "last_event_time_s": None}
                     for sid in range(1, segment_count + 1)}
    for event in events:
        sid = node_segments.get(int(event["event_node_id"]))
        if sid is None:
            continue
        event["segment_id"] = sid
        event_summary[sid]["event_count"] += 1
        event_summary[sid]["last_event_time_s"] = event["time_s"]
    segment_spacetime = []
    for item in segment_profile:
        sid = int(item["segment_id"])
        current_sg = item.get("segment_mean_Sg")
        initial_value = initial_sg.get(sid)
        segment_spacetime.append({**item, **event_summary.get(sid, {}),
                                  **initial_inventory.get(sid, {}),
                                  "initial_segment_mean_Sg": initial_value,
                                  "segment_Sg_change": current_sg - initial_value
                                      if current_sg is not None and initial_value is not None else None})
    parameters = []
    def parameter(name: str, configured: Any, effective_value: Any, unit: str, source: str, classification: str) -> None:
        conflict = configured is not None and effective_value is not None and str(configured) != str(effective_value)
        parameters.append({"name": name, "configured": configured, "effective": effective_value,
                           "unit": unit, "source": source, "classification": classification,
                           "conflict": conflict})
    parameter("backend", (config or {}).get("linear", {}).get("backend"), effective["backend"], "", "config / manifest / runtime", "assumed")
    parameter("PARDISO threads", (config or {}).get("linear", {}).get("pardiso_threads"),
              effective["pardiso_threads"], "threads", "config / manifest / runtime", "assumed")
    parameter("batch microbubble retirement", batch_cfg.get("enabled", False),
              effective["batch_microbubble_retirement_enabled"], "", "config / manifest", "assumed")
    parameter("batch individual gas fraction", batch_cfg.get("maximum_individual_free_gas_fraction"),
              effective["batch_individual_free_gas_fraction"], "1", "config / manifest", "assumed")
    parameter("batch aggregate gas fraction", batch_cfg.get("maximum_batch_free_gas_fraction_per_step"),
              effective["batch_total_free_gas_fraction_per_step"], "1 per accepted step", "config / manifest", "assumed")
    parameter("dt_max_s", time_cfg.get("dt_max_s", time_legacy.get("dt_max_s")), effective["dt_max_s"], "s", "config / checkpoint", "assumed")
    parameter("mass_transfer_multiplier", physics.get("mass_transfer_multiplier", 1.0),
              (checkpoint or {}).get("mass_transfer_multiplier", multiplier), "1", "config / checkpoint / diagnostics", "preview-only" if kind == "accelerated_preview" else "assumed")
    for name, key, unit, classification in (
        ("temperature", "T_K", "K", "assumed"), ("water density", "rho_l_kg_m3", "kg m^-3", "assumed"),
        ("water viscosity", "mu_l_Pa_s", "Pa s", "assumed"), ("gas molar mass", "molar_mass_kg_mol", "kg mol^-1", "assumed"),
        ("Henry coefficient", "Hcp_mol_m3_Pa", "mol m^-3 Pa^-1", "assumed"),
        ("molecular diffusion", "Dm_m2_s", "m^2 s^-1", "assumed"),
        ("surface tension", "sigma_N_m", "N m^-1", "assumed"), ("contact angle", "contact_angle_rad", "rad", "unconfirmed"),
        ("mass-transfer coefficient", "kL_m_s", "m s^-1", "unconfirmed"),
    ):
        parameter(name, physics.get(key), physics.get(key), unit, "config / physical_parameter_audit", classification)
    parameter("residual tolerance", (config or {}).get("nonlinear", {}).get("residual_tolerance"), effective["residual_tolerance"], "1", "config / manifest", "assumed")
    parameter("update tolerance", (config or {}).get("nonlinear", {}).get("update_tolerance"), effective["update_tolerance"], "1", "config / manifest", "assumed")
    parameter("gas constant", None, 8.31446261815324, "J mol^-1 K^-1", "compiled solver constant", "assumed")
    parameter("gas-liquid interfacial area", None, "(36*pi)^(1/3) * Vg^(2/3)", "m^2", "solver / physical_parameter_audit", "unconfirmed")
    parameter("capillary pressure", None, "Young-Laplace spherical bubble", "Pa", "solver", "unconfirmed")
    parameter("relative permeability/conductance", physics.get("kr_exponent_m"), "kr=Sw^m; throat Poiseuille conductance", "", "config / solver", "unconfirmed")
    case = {
        "case_name": case_name, "run_id": (manifest or {}).get("run_id", case_dir.name),
        "description": (config or {}).get("case", {}).get("notes", ""), "kind": kind,
        "scientific_interpretation": "非真实传质时间预测" if kind == "accelerated_preview" else "scientific baseline",
        "status": status, "directory": str(output_path), "config_path": str(config_path) if config_path else None,
        "relative_path": str(output_path.relative_to(root.resolve())),
        "control_allowed": output_path.parent == root.resolve(),
        "scheduler_status": scheduler_status,
        "config_sha256": config_hash, "manifest_config_sha256": manifest_config_hash,
        "source_sha256": (manifest or {}).get("build", {}).get("project_source_tree_sha256"),
        "executable": (manifest or {}).get("run_argv", [None])[0] if manifest else None,
        "executable_sha256": None, "started_utc": (manifest or {}).get("started_utc", run_values.get("started_utc")),
        "last_update": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime((case_dir / "diagnostics.tsv").stat().st_mtime)) if (case_dir / "diagnostics.tsv").exists() else None,
        "manifest": manifest_view, "checkpoint": checkpoint_view, "checkpoint_sha256": sha256(case_dir / "checkpoint.json"),
        "last_diagnostics": row, "termination_events": read_json(case_dir / "termination_events.json"),
        "config": config, "physics": physics, "effective": effective,
        "processes": processes, "geometry_audit": geometry_audit,
        "physical_audit": read_json(case_dir / "physical_parameter_audit.json"), "warnings": [],
        "recent_diagnostics": recent, "global_latest": global_row, "parameters": parameters,
        "history_series": series,
        "conservation_series": conservation,
        "segment_profile": segment_profile,
        "disappearance_events": events,
        "segment_spacetime": segment_spacetime,
        "pnm_inputs": pnm, "segment_pore_volume_m3": segment_volume,
        "geometry_network": {
            "image_dimensions_voxels": (config or {}).get("geometry", {}).get("image_dimensions_voxels"),
            "voxel_size_um": (config or {}).get("geometry", {}).get("voxel_size_um"),
            "flow_axis": (config or {}).get("virtual_boundary", {}).get("boundary_axis"),
            "nodes": network.get("control_volumes", (geometry_audit or {}).get("raw_node_count", run_values.get("nodes"))),
            "edges": network.get("half_edges", (geometry_audit or {}).get("raw_edge_count", run_values.get("edges"))),
            "virtual_boundary_nodes": network.get("virtual_boundary_nodes", run_values.get("virtual_nodes")),
            "inlet_virtual_nodes": inlet_virtual,
            "outlet_virtual_nodes": outlet_virtual,
            "boundary_edges": network.get("boundary_edges", run_values.get("boundary_edges")),
            "main_real_pore_volume_m3": global_row.get("main_real_pore_volume_m3"),
            "control_volume_total_m3": segment_volume,
            "initial_bubbles": (read_json(case_dir / "physical_parameter_audit.json") or {}).get("initial_active_bubble_count"),
            "active_bubbles": row.get("active_gas_count") if row else None,
            "segments": (config or {}).get("analysis", {}).get("x_segments"),
        },
        "boundary_runtime": {
            "inlet_boundary_type": (config or {}).get("boundary", {}).get("inlet_type", "pressure"),
            "configured_inlet_flow_m3_s": (config or {}).get("boundary", {}).get("inlet_flow_m3_s"),
            "outlet_boundary_type": (config or {}).get("boundary", {}).get("outlet_type", "pressure"),
            "pressure_type": "Dirichlet",
            "inlet_relative_pressure_Pa": (config or {}).get("boundary", {}).get("Pl_in_rel_Pa"),
            "outlet_relative_pressure_Pa": (config or {}).get("boundary", {}).get("Pl_out_rel_Pa"),
            "background_absolute_pressure_Pa": (config or {}).get("pressure_reference", {}).get("background_abs_Pa"),
            "concentration_type": (config or {}).get("boundary", {}).get("species_mode"),
            "copy_adjacent_concentration": (config or {}).get("boundary", {}).get("species_mode") == "copy_adjacent_concentration",
            "diffusive_boundary_flux_zero": (config or {}).get("boundary", {}).get("outlet_diffusive_condition") == "zero_gradient",
            "advective_boundary": "upwind, enabled",
            "flow_axis": (config or {}).get("virtual_boundary", {}).get("boundary_axis"),
            "inlet_flow_m3_s": row.get("inlet_liquid_inflow") if row else None,
            "outlet_flow_m3_s": row.get("outlet_liquid_outflow") if row else None,
            "throughflow_m3_s": row.get("throughflow_rate") if row else None,
            "cumulative_makeup_liquid_m3": row.get("cumulative_net_liquid_inflow") if row else None,
            "inlet_reverse_flow_edges": row.get("inlet_reverse_flow_edges") if row else None,
            "inlet_reverse_flow_m3_s": row.get("inlet_reverse_flow_m3_s") if row else None,
            "outlet_backflow_edges": row.get("outlet_backflow_edges") if row else None,
            "outlet_backflow_m3_s": row.get("outlet_backflow_m3_s") if row else None,
        },
    }
    diagnostics_path = case_dir / "diagnostics.tsv"
    case["diagnostics_age_seconds"] = max(0.0, time.time() - diagnostics_path.stat().st_mtime) if diagnostics_path.exists() else None
    case["wall_clock_runtime_seconds"] = max(
        (p.get("elapsed_seconds", 0) for p in processes),
        default=(manifest or {}).get("elapsed_seconds"))
    bubble_process = next((p for p in processes if p.get("executable") and "bubble_solver" in p["executable"]), None)
    if bubble_process:
        case["executable"] = bubble_process["executable"].replace(" (deleted)", "")
        case["executable_sha256"] = sha256(Path(f"/proc/{bubble_process['pid']}/exe"))
        case["executable_path_current_sha256"] = sha256(Path(case["executable"]))
        if not case["source_sha256"]:
            case["source_sha256"] = embedded_source_hash(Path(f"/proc/{bubble_process['pid']}/exe"))
    if case["executable"] and not bubble_process:
        case["executable_sha256"] = sha256(Path(case["executable"]))
        if not case["source_sha256"]:
            case["source_sha256"] = embedded_source_hash(Path(case["executable"]))
    case["warnings"] = warnings_for(case)
    return case


def snapshot_case_brief(case_dir: Path, root: Path, config_path: Path | None,
                        process_table: list[dict[str, Any]],
                        scheduler_status: str | None) -> dict[str, Any]:
    """Cheap list-row snapshot; expensive histories are loaded by /api/case."""
    config = read_config(config_path)
    manifest_path = case_dir / "run_manifest.json"
    manifest_exists = manifest_path.is_file()
    checkpoint_exists = (case_dir / "checkpoint.json").is_file()
    row = last_row(case_dir / "diagnostics.tsv")
    output = case_dir.resolve()
    configured = configured_output(config or {})
    processes = process_snapshot(
        output, config_path if configured == output else None, process_table)
    termination = read_json(case_dir / "termination_events.json")
    term_reason = (termination or {}).get("termination_reason")
    if processes:
        status = "running"
    elif ((case_dir / "run_manifest.failed.json").is_file() or
          scheduler_status == "failed"):
        status = "failed"
    elif manifest_exists:
        status = "interrupted" if term_reason in (
            "safe_stop_requested", "wall_clock_limit") else "completed"
    elif scheduler_status:
        status = scheduler_status
    elif case_dir.exists() and any(case_dir.iterdir()):
        status = "interrupted"
    else:
        status = "queued"
    case_name = ((config or {}).get("case", {}).get("name") or case_dir.name)
    physics = (config or {}).get("physics", {})
    multiplier = physics.get("mass_transfer_multiplier", 1.0)
    kind = ("accelerated_preview" if float(multiplier or 1.0) != 1.0 or
            "accelerated" in case_name.lower() else "scientific baseline")
    time_cfg = (config or {}).get("time_step_control", {})
    legacy_time = (config or {}).get("time", {})
    checkpoint_view = ({"format": "present", "time_s": (row or {}).get("time_s"),
                        "step_number": (row or {}).get("step"),
                        "next_dt_s": (row or {}).get("next_dt_s")}
                       if checkpoint_exists else None)
    last_update = None
    diagnostics = case_dir / "diagnostics.tsv"
    if diagnostics.is_file():
        last_update = time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                    time.gmtime(diagnostics.stat().st_mtime))
    initial_free = None
    if termination:
        initial_free = termination.get("initial_total_free_gas_moles")
    retained = None
    try:
        current_free = float((row or {}).get("total_free_gas_moles",
                                             (row or {}).get("free_mol_total")))
        if initial_free and float(initial_free) > 0:
            retained = current_free / float(initial_free)
    except (TypeError, ValueError):
        pass
    return {
        "case_name": case_name,
        "run_id": case_name,
        "relative_path": str(output.relative_to(root.resolve())),
        "directory": str(output),
        "status": status,
        "scheduler_status": scheduler_status,
        "variant": ("PNM-0.6" if case_name.endswith("_p06") else
                    "PNM-0.8" if case_name.endswith("_p08") else None),
        "kind": kind,
        "scientific_interpretation": (
            "非真实传质时间预测" if kind == "accelerated_preview"
            else "scientific baseline"),
        "config_path": str(config_path) if config_path else None,
        "checkpoint": checkpoint_view,
        "last_diagnostics": row,
        "last_update": last_update,
        "free_gas_retained_fraction": retained,
        "processes": processes,
        "wall_clock_runtime_seconds": max(
            (item.get("elapsed_seconds", 0) for item in processes), default=None),
        "effective": {
            "backend": (config or {}).get("linear", {}).get("backend"),
            "current_dt_s": ((row or {}).get("dt_s") or
                             (checkpoint_view or {}).get("next_dt_s")),
            "dt_min_s": time_cfg.get("dt_min_s", legacy_time.get("dt_min_s")),
            "dt_max_s": time_cfg.get("dt_max_s", legacy_time.get("dt_max_s")),
        },
        "mass_transfer_multiplier": multiplier,
        "control_allowed": output.parent == root.resolve(),
    }


def snapshot(root: Path) -> dict[str, Any]:
    root = root.resolve()
    if not root.is_dir():
        raise SystemExit(f"production root is not a directory: {root}")
    catalog = config_catalog(root)
    status_by_id, schedulers = scheduler_statuses(root)
    process_table = bubble_process_table()
    actual_dirs: set[Path] = set()
    for marker in ("diagnostics.tsv", "checkpoint.json", "run.log",
                   "driver.log", "run_manifest.json", "run_manifest.failed.json"):
        for path in root.rglob(marker):
            if safe_child(root, path.parent):
                actual_dirs.add(path.parent.resolve())
    case_dirs = sorted(set(catalog) | actual_dirs)
    cases = []
    for case_dir in case_dirs:
        config_path = catalog.get(case_dir)
        config = read_config(config_path)
        run_id = (config or {}).get("case", {}).get("name", case_dir.name)
        cases.append(snapshot_case_brief(
            case_dir, root, config_path, process_table,
            status_by_id.get(run_id)))
    priority = {"failed": 0, "running": 1, "interrupted": 2,
                "queued": 3, "completed": 4}
    cases.sort(key=lambda item: (priority.get(item["status"], 9),
                                 item["case_name"]))
    counts = {name: sum(case["status"] == name for case in cases)
              for name in priority}
    variants = {
        variant: {name: sum(case.get("variant") == variant and
                            case["status"] == name for case in cases)
                  for name in priority}
        for variant in ("PNM-0.6", "PNM-0.8")
    }
    return {"production_root": str(root),
            "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "summary": {"total": len(cases), "status_counts": counts,
                        "variant_status_counts": variants,
                        "schedulers": schedulers},
            "cases": cases}


HTML = r'''<!doctype html><meta charset="utf-8"><title>Bubble dissolution cases</title>
<style>:root{--ink:#172b3a;--muted:#6b7c8f;--line:#e6ebf0;--panel:#fff;--bg:#f3f6f9;--blue:#1769aa;--red:#a33a3a;--green:#19754a;--purple:#7950a8;--orange:#c7791f}*{box-sizing:border-box}body{font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;margin:0;background:var(--bg);color:var(--ink)}header{padding:22px 30px;background:linear-gradient(120deg,#102e49,#1a5276);color:#fff;box-shadow:0 2px 8px #102e4933}header h2{margin:0 0 4px;font-size:22px}header #root{opacity:.78;font:12px ui-monospace,monospace;overflow-wrap:anywhere}main{display:grid;grid-template-columns:340px minmax(0,1fr);gap:18px;padding:18px;max-width:1900px;margin:auto}.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:18px;margin-bottom:14px;box-shadow:0 2px 10px #2030400b}.card h2,.card h3{margin:0 0 14px}.card h3{font-size:16px;color:#21445e}.case-list{position:sticky;top:14px;max-height:calc(100vh - 40px);overflow:auto}.case{display:block;width:100%;text-align:left;padding:12px;border:1px solid transparent;border-radius:8px;background:#fff;color:var(--ink);cursor:pointer;margin:6px 0}.case:hover{border-color:#a9c8df;background:#f5fbff}.case.sel{border-color:#4d9bd0;background:#eaf5fd}.case-name{font-weight:650;display:block;overflow-wrap:anywhere}.case-meta{display:block;margin-top:3px;color:var(--muted);font-size:12px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}.badge{display:inline-block;padding:3px 8px;border-radius:999px;font-size:11px;font-weight:650;margin:0 5px 8px 0;background:#e8eef3;color:#415568}.badge.running{background:#e2f5e9;color:var(--green)}.badge.queued{background:#edf1f5;color:#596b79}.badge.failed{background:#fde8e8;color:var(--red)}.badge.completed{background:#e7f3ff;color:var(--blue)}.badge.accelerated_preview{background:#fff0d6;color:#8a5311}.badge.scientific{background:#e8eefb;color:#315b9a}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:0 18px}.kv{padding:9px 0;border-bottom:1px solid var(--line);min-width:0}.k{display:block;color:var(--muted);font-size:12px;margin-bottom:2px}.v{display:block;font:12px ui-monospace,SFMono-Regular,Menlo,monospace;overflow-wrap:anywhere}.src{display:block;color:#94a1ad;font-size:10px;margin-top:2px}.warn{color:var(--red);background:#fff3f3;border-left:4px solid #d86a6a;padding:8px 10px;margin:6px 0;border-radius:4px}.table-wrap{overflow:auto}.charts{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}.chart-title{font-weight:650;color:#21445e;margin-bottom:6px}.chart-wrap{min-width:0}.chart-wrap svg{display:block;width:100%;height:220px;overflow:visible}.chart-grid{stroke:var(--line);stroke-width:1}.chart-axis{stroke:#9aa9b5;stroke-width:1}.chart-label{fill:var(--muted);font-size:11px}.chart-line{fill:none;stroke-width:2;vector-effect:non-scaling-stroke}.chart-legend{display:flex;flex-wrap:wrap;gap:10px;margin-top:5px;color:var(--muted);font-size:11px}.chart-legend span{white-space:nowrap}.chart-legend i{display:inline-block;width:18px;height:3px;vertical-align:middle;margin-right:4px}.chart-empty{height:220px;display:grid;place-items:center;color:var(--muted);border:1px dashed var(--line)}table{border-collapse:collapse;width:100%;font-size:12px;min-width:700px}td,th{padding:8px 7px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}th{color:var(--muted);font-weight:600;position:sticky;top:0;background:#fff}.conflict{background:#fff1df}.batch-summary{display:grid;grid-template-columns:repeat(3,1fr);gap:6px;margin:12px 0}.summary-cell{background:#f4f8fb;border:1px solid var(--line);border-radius:7px;padding:7px;text-align:center}.summary-cell b{display:block;font-size:17px;color:#21445e}.summary-cell span{font-size:10px;color:var(--muted)}.case-filter{display:grid;grid-template-columns:1fr 1fr;gap:7px;margin:10px 0}.case-filter input{grid-column:1/-1}.case-filter input,.case-filter select{width:100%;padding:7px;border:1px solid #cdd8e0;border-radius:6px;background:#fff}.list-count{font-size:11px;color:var(--muted);margin-bottom:6px}@media(max-width:900px){main{grid-template-columns:1fr;padding:10px}.case-list{position:static;max-height:360px}.grid{grid-template-columns:repeat(auto-fit,minmax(160px,1fr))}.charts{grid-template-columns:1fr}}
@keyframes pulse{0%,100%{opacity:1}50%{opacity:.35}}.live-dot{display:inline-block;width:9px;height:9px;border-radius:50%;background:var(--green);margin-right:5px}.live-dot.on{animation:pulse 1.4s infinite}.runtime-strip{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px;margin:10px 0}.runtime-item{padding:10px;border:1px solid var(--line);border-radius:8px;background:#f8fbfd;min-width:0}.runtime-item .label{display:block;color:var(--muted);font-size:11px}.runtime-item .value{display:block;font:14px ui-monospace,monospace;font-weight:600;white-space:nowrap}.runtime-item .context{display:block;color:var(--muted);font-size:11px;margin-top:2px;white-space:nowrap}.progress-track{height:6px;background:#e7edf1;border-radius:4px;overflow:hidden;margin-top:6px}.progress-bar{height:100%;background:var(--blue);transition:width .5s ease}
.data-group th{padding-top:18px;color:#21445e;background:#f4f8fb;font-weight:600;border-bottom:2px solid #dce6ed}
.control-note{color:var(--muted)}.control-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px}.control-grid label{display:grid;gap:5px;color:var(--muted);font-size:12px}.control-grid input{width:100%;padding:8px;border:1px solid var(--line);border-radius:6px;background:#fff;color:var(--ink);font:13px ui-monospace,monospace}.control-actions{display:flex;flex-wrap:wrap;gap:8px;margin-top:14px}.control-actions button{width:auto;border:1px solid #9cb4c5;border-radius:6px;padding:8px 12px;background:#f6fafc;color:var(--ink);cursor:pointer}.control-actions button:first-child{border-color:#d27b7b;color:#9c2929}.control-actions button:disabled{opacity:.45;cursor:not-allowed}.control-result{white-space:pre-wrap;max-height:180px;overflow:auto;padding:9px;background:#f4f7f9;border-radius:6px;color:var(--muted);font-size:11px}.control-result.ok{color:var(--green)}
.console-card{padding:0;overflow:hidden;border-color:#213544}.console-head{display:flex;align-items:center;justify-content:space-between;gap:12px;padding:13px 16px;background:#152633;color:#dceaf2}.console-title{display:flex;align-items:center;gap:8px;font-weight:650}.console-actions{display:flex;align-items:center;gap:8px}.console-actions button{border:1px solid #486171;border-radius:6px;background:#203847;color:#dceaf2;padding:5px 10px;cursor:pointer}.console-actions button:hover{background:#29495c}.console-screen{margin:0;height:360px;overflow:auto;padding:14px 16px;background:#07131b;color:#c7d8e2;font:13px/1.55 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;white-space:pre-wrap;overflow-wrap:anywhere;tab-size:2}.console-line{display:block;min-height:1.55em}.console-line.step{color:#c7e9d7}.console-line.retry{color:#ffd28a}.console-line.event{color:#d7b9ff;background:#5d3a7b33}.console-line.error{color:#ff9e9e;background:#8f232333}.console-line.meta{color:#8fb9d3}.console-cursor{display:inline-block;width:8px;height:1.1em;margin-left:5px;vertical-align:-2px;background:#5ed294;animation:pulse 1s infinite}.console-foot{display:flex;flex-wrap:wrap;justify-content:space-between;gap:8px;padding:9px 16px;background:#10212c;color:#91a8b6;font-size:12px}.console-empty{color:#91a8b6}.event-summary{display:flex;flex-wrap:wrap;gap:8px;margin:0 0 14px}.event-pill{padding:7px 10px;border-radius:7px;background:#f5f1fb;color:#604080;border:1px solid #dbcdeb;font-size:12px}.event-pill.exact{background:#fff2f2;color:#933;border-color:#efcccc}.event-line{vector-effect:non-scaling-stroke;stroke-width:1.5}.event-dot{stroke:#fff;stroke-width:1;vector-effect:non-scaling-stroke}
</style>
<header><h2>气泡溶解算例监控 / Bubble dissolution case monitor</h2><div id="root"></div></header><main><section><div class="card case-list"><b>全部生产算例 / All production cases</b><div id="batch-summary" class="batch-summary"></div><div class="case-filter"><input id="case-search" placeholder="搜索算例 / Search"><select id="variant-filter"><option value="">全部级配 / Both variants</option><option>PNM-0.6</option><option>PNM-0.8</option></select><select id="status-filter"><option value="">全部状态 / All status</option><option value="running">运行中 / Running</option><option value="queued">排队中 / Queued</option><option value="completed">已完成 / Completed</option><option value="failed">失败 / Failed</option><option value="interrupted">已中断 / Interrupted</option></select></div><div id="list-count" class="list-count"></div><div id="cases">正在读取 / Loading…</div></div></section><section id="detail"><div class="card">请选择一个算例 / Select a case.</div></section></main>
<script>
let data; let selectedDirectory=null; let selectedRelative=null; let selectedCase=null; let controlToken=''; let consolePaused=false; let consoleAutoScroll=true; const esc=s=>String(s??'—').replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));
function fmt(v,key=''){if(v===null||v===undefined||v==='')return '—';if(typeof v==='boolean')return v?'yes':'no';if(Array.isArray(v))return v.map(x=>fmt(x)).join(' × ');const s=String(v);if(!/^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/.test(s))return s;const n=Number(s);if(!Number.isFinite(n))return s;if(/(?:^|_)(?:step|id|count|iterations|events)$/.test(key)&&Number.isInteger(n))return String(n);if(n===0)return '0';const a=Math.abs(n);if(a>=1e4||a<1e-3)return n.toExponential(3);if(a>=100)return n.toFixed(2);if(a>=1)return n.toFixed(3);return n.toPrecision(4);}
function shanghaiTime(v){if(!v)return '—';const d=v instanceof Date?v:new Date(v);if(Number.isNaN(d.valueOf()))return String(v);return new Intl.DateTimeFormat('zh-CN',{timeZone:'Asia/Shanghai',year:'numeric',month:'2-digit',day:'2-digit',hour:'2-digit',minute:'2-digit',second:'2-digit',hour12:false}).format(d);}
const statusZh=s=>({running:'运行中',queued:'排队中',completed:'已完成',failed:'失败',interrupted:'已中断'}[s]||s||'未知');
const kindZh=s=>({accelerated_preview:'加速预览',scientific_baseline:'科研基准','scientific baseline':'科研基准'}[s]||s||'未分类');
const statusBi=s=>({running:'运行中 / Running',queued:'排队中 / Queued',completed:'已完成 / Completed',failed:'失败 / Failed',interrupted:'已中断 / Interrupted'}[s]||`${statusZh(s)} / ${s||'Unknown'}`);
const kindBi=s=>({accelerated_preview:'加速预览 / Accelerated preview',scientific_baseline:'科研基准 / Scientific baseline','scientific baseline':'科研基准 / Scientific baseline'}[s]||`${kindZh(s)} / ${s||'Unclassified'}`);
const classZh=s=>({confirmed:'已确认 / confirmed',assumed:'假设值 / assumed','preview-only':'仅用于预览 / preview-only',unconfirmed:'尚未确认 / unconfirmed'}[s]||s||'未分类 / unclassified');
const sourceZh=s=>String(s||'').replaceAll('physical_parameter_audit','物性参数审计 physical_parameter_audit').replaceAll('geometry audit','几何审计 geometry audit').replaceAll('config','配置文件 config').replaceAll('checkpoint','检查点 checkpoint').replaceAll('diagnostics','诊断记录 diagnostics').replaceAll('manifest','运行清单 manifest').replaceAll('runtime','运行时 runtime').replaceAll('solver','求解器 solver');
const uiZh={'Runtime curves':'运行历程曲线 / Runtime curves','Gas inventory':'气体摩尔量 / Gas inventory','Gas saturation':'气体饱和度 / Gas saturation','Accepted time step':'已接受时间步长 / Accepted time step','Conservation error':'守恒误差 / Conservation error','Flow-direction segment profile':'沿流向分段统计 / Flow-direction segment profile','Latest segment gas saturation':'最新分段气体饱和度 / Latest segment gas saturation','Geometry & network':'几何与网络 / Geometry and network','Physics':'物理参数 / Physics','Boundary & numerical parameters':'边界与数值参数 / Boundary and numerical parameters','Geometry, boundary and provenance':'几何、边界与溯源 / Geometry, boundary and provenance','PNM inputs':'PNM输入文件 / PNM inputs','Solver heartbeat':'求解器运行心跳 / Solver heartbeat','Segment':'区段 / Segment','Current Sg':'当前气体饱和度 / Current Sg','Active bubbles':'活跃气泡数 / Active bubbles','Gas volume [m³]':'气体体积 / Gas volume [m³]','Disappearance events':'气泡消失事件数 / Disappearance events','Last event time [s]':'最近事件时间 / Last event time [s]'};
Object.assign(uiZh,{'配置值与实际运行值':'配置值与实际运行值 / Configured and effective values','最新运行数据':'最新运行数据 / Latest runtime data','参数（中文解释 / 原字段）':'参数 / Parameter（中文解释 · raw field）','配置值':'配置值 / Configured','实际值':'实际值 / Effective','单位 / 来源':'单位与来源 / Unit and source','可信度分类':'可信度分类 / Confidence class','数值':'数值 / Value','数据来源':'数据来源 / Data source','求解器状态':'求解器状态 / Solver status','已接受时间步':'已接受时间步 / Accepted step','物理模拟时间':'物理模拟时间 / Simulation time','真实墙钟耗时':'真实墙钟耗时 / Wall-clock time','接受步长 / Newton迭代':'接受步长与Newton迭代 / Accepted dt and Newton iterations','最新接受数据（上海时间）':'最新接受数据（上海时间）/ Latest accepted data (Shanghai)','面板刷新时间（上海时间）':'面板刷新时间（上海时间）/ Dashboard refresh (Shanghai)'});
Object.assign(uiZh,{role:'用途 / Role',path:'路径 / Path','SHA-256':'哈希 / SHA-256','分段气体饱和度 Sg':'分段气体饱和度 / Segment gas saturation Sg','沿流向归一化区段':'沿流向归一化区段 / Normalized x segment'});
function translateUI(){document.querySelectorAll('h3,.chart-title,th,.runtime-item .label,svg text').forEach(e=>{if(uiZh[e.textContent.trim()])e.textContent=uiZh[e.textContent.trim()];});const sh=document.querySelector('#segment-table-body')?.closest('table')?.querySelector('thead tr');if(sh)sh.innerHTML='<th>区段 / Segment</th><th>初始 Sg / Initial Sg</th><th>当前 Sg / Current Sg</th><th>Sg 变化量 / Change</th><th>初始气泡数 / Initial bubbles</th><th>当前气泡数 / Active bubbles</th><th>初始气体体积 / Initial gas volume [m³]</th><th>当前气体体积 / Current gas volume [m³]</th><th>消失事件数 / Events</th><th>最近事件时间 / Last event [s]</th>';const pl=document.querySelector('#chart-profile .chart-legend');if(pl)pl.innerHTML=pl.innerHTML.replace('初始 Sg','初始 Sg / Initial Sg').replace('当前 Sg','当前 Sg / Current Sg').replace('个区段','个区段 / segments');}
const fieldInfo={
status:['运行状态','求解器当前是否运行、完成、中断或失败'],backend:['线性求解后端','实际使用的Eigen或AMGX线性求解器'],dt_max_s:['最大时间步长','自适应时间步允许增长到的上限'],temperature:['温度','液体和气体物性采用的绝对温度'],"water density":['水密度','液相密度'],"water viscosity":['水动力黏度','液相动力黏度'],"gas molar mass":['气体摩尔质量','有效空气组分的摩尔质量'],"Henry coefficient":['Henry溶解系数','气液平衡浓度与绝对气压的比例系数'],"molecular diffusion":['分子扩散系数','溶解气体在水中的分子扩散系数'],"surface tension":['表面张力','气液界面表面张力'],"contact angle":['接触角','毛细压力模型采用的接触角'],"mass-transfer coefficient":['传质系数 kL','未乘倍率前的气液界面传质系数'],"residual tolerance":['残量容差','Newton缩放残量接受阈值'],"update tolerance":['更新量容差','Newton相对状态增量接受阈值'],"gas constant":['通用气体常数','理想气体状态方程采用的R'],"gas-liquid interfacial area":['气液界面面积模型','根据局部气体体积计算界面面积'],"capillary pressure":['毛细压力模型','连接气液压力差与局部饱和度'],"relative permeability/conductance":['相对渗透率与导流模型','饱和度对液相导流能力的影响'],step:['接受步编号','已通过全部收敛与守恒检查的时间步编号'],time_s:['物理时间','模型内部累计的物理模拟时间'],physical_time_s:['物理时间','模型内部累计的物理模拟时间'],dt_s:['接受时间步长','本次成功接受的物理时间增量'],requested_dt_s:['请求时间步长','自适应控制器最初请求尝试的步长'],next_dt_s:['下一步建议步长','自适应控制器给下一次求解的建议值'],timestep_retries:['重试次数','当前接受步在成功前缩小时间步重算的次数'],newton_iterations:['Newton迭代次数','当前接受步的非线性迭代次数'],linear_iterations:['线性迭代次数','当前接受步累计线性求解迭代数'],residual_scaled_inf:['缩放残量无穷范数','Newton方程缩放后最大残量'],relative_update_norm:['相对更新范数','Newton状态增量相对于当前状态的大小'],epsilon_N_adjusted:['组分摩尔守恒误差','计入active-set账本后的相对摩尔守恒误差'],epsilon_V_adjusted:['液体体积守恒误差','计入active-set账本后的相对体积守恒误差'],overall_Sg:['总体气体饱和度','真实样品控制体的体积加权平均气体饱和度'],main_connected_component_Sg:['主贯通域气体饱和度','入口—出口主连通分量的体积加权气体饱和度'],total_free_gas_moles:['自由气体摩尔数','当前仍处于气泡相中的气体总摩尔数'],total_dissolved_moles:['溶解气体摩尔数','当前液相中溶解气体总摩尔数'],total_gas_moles:['内部气体组分库存','自由气体与液相溶解气体库存之和'],active_gas_count:['活跃气泡数','当前仍包含自由气相的控制体数量'],gas_disappearance_events:['累计气泡消失事件','已接受的active-set气泡消失事件总数'],min_Sg:['最小局部气体饱和度','所有活跃气泡中的最小Sg'],transfer_rate_mol_s:['界面传质速率','当前所有气液界面传质源项之和'],inlet_liquid_inflow:['入口净流入量','所有入口边界喉液体通量的有符号总和'],outlet_liquid_outflow:['出口净流出量','所有出口边界喉液体通量的有符号总和'],throughflow_rate:['贯穿流量','入口与出口主流量的代表值'],net_liquid_inflow:['净补液流量','入口流量减出口流量，用于填充气泡收缩空间'],inlet_reverse_flow_edges:['入口反流喉数','从多孔介质内部流向入口储层的边界喉数量'],inlet_reverse_flow_m3_s:['入口反流量','所有入口反流喉流量绝对值之和'],outlet_backflow_edges:['出口回流喉数','从出口储层流回样品的边界喉数量'],outlet_backflow_m3_s:['出口回流量','所有出口回流喉流量绝对值之和'],inlet_pressure_relative_Pa:['入口相对压力','入口压力相对于背景绝对压力的值'],inlet_pressure_absolute_Pa:['入口绝对压力','背景绝对压力加入口相对压力'],inlet_flow_target_m3_s:['目标入口流量','Neumann流量边界指定的入口总流量'],mass_transfer_multiplier:['传质倍率','对原始气液传质系数统一施加的倍率'],segment_mean_Sg:['分段平均气体饱和度','该归一化x区段内的体积加权平均Sg'],segment_active_bubble_count:['分段活跃气泡数','该区段内仍含自由气相的控制体数量']};
Object.assign(fieldInfo,{backend:['线性求解后端 / Linear backend','实际使用的Eigen、AMGX或PARDISO线性求解器'],'PARDISO threads':['PARDISO线程数 / PARDISO threads','PARDISO数值分解与回代使用的CPU线程数']});
Object.assign(fieldInfo,{'batch microbubble retirement':['微小气泡批量退休 / Batch microbubble retirement','是否允许质量可忽略的阈值穿越气泡在同一接受步批量退出active-set'],'batch individual gas fraction':['单泡批量上限 / Individual batch limit','单个候选气泡相对于全局自由气体的最大比例'],'batch aggregate gas fraction':['单步批量总上限 / Aggregate batch limit','同一接受步全部批量候选的最大总自由气体比例']});
const diagnosticsZh={linear_residual:'线性方程真实相对残量',state_scaled_inf:'状态量缩放无穷范数',update_scaled_inf:'状态更新缩放无穷范数',update_tolerance_limit:'更新量接受上限',update_criterion_ratio:'更新判据比值',residual_criterion_passed:'残量判据是否通过',update_criterion_passed:'更新量判据是否通过',Rl_scaled_inf:'液体体积残量 Rl',Rd_scaled_inf:'溶解组分残量 Rd',Rg_scaled_inf:'自由气体残量 Rg',Rc_scaled_inf:'毛细闭合残量 Rc',Rl_max_node_id:'最大 Rl 所在节点',Rd_max_node_id:'最大 Rd 所在节点',Rg_max_node_id:'最大 Rg 所在节点',Rc_max_node_id:'最大 Rc 所在节点',boundary_liquid_outflow_m3_s:'边界液体净流出率',dissolved_mol:'溶解气体摩尔库存',free_mol:'自由气体摩尔库存',boundary_component_outflow_mol_s:'边界组分净流出率',epsilon_N_equation:'方程层摩尔守恒误差',epsilon_V_equation:'方程层体积守恒误差',threshold_moles_transferred:'阈值切换转移摩尔数',active_set_mole_ledger:'活跃集摩尔账本',active_set_liquid_volume_ledger:'活跃集液体体积账本',free_mol_total:'自由气体总摩尔数',free_mol_main_flow_component:'主贯通域自由气体摩尔数',free_mol_boundary_connected:'边界连通区自由气体摩尔数',free_mol_boundaryless:'非边界连通区自由气体摩尔数',gas_volume_total:'气体总体积',gas_volume_main_flow_component:'主贯通域气体体积',gas_volume_boundary_connected:'边界连通区气体体积',gas_volume_boundaryless:'非边界连通区气体体积',dissolved_mol_total:'溶解气体总摩尔数',boundary_advective_component_outflow:'边界对流组分净流出率',boundary_diffusive_component_outflow:'边界扩散组分净流出率',cumulative_net_liquid_inflow:'累计净补液体积',liquid_volume_change:'本步液体体积增量',cumulative_liquid_volume_change:'累计液体体积增量',wrong_direction_boundary_edges:'反向边界喉总数',internal_pressure_min_Pa:'内部液相相对压力最小值',internal_pressure_max_Pa:'内部液相相对压力最大值',maximum_dissolution_driving_force:'最大溶解驱动力',active_main_flow:'主贯通域活跃气泡数',active_boundary_connected:'边界连通区活跃气泡数',active_boundaryless:'非边界连通区活跃气泡数',cumulative_boundary_component_outflow_mol:'累计边界组分净流出量',inlet_advective_component_inflow:'入口对流组分流入率',outlet_adv_component_outflow:'出口组分流出率',outlet_advective_component_outflow:'出口对流组分流出率',cumulative_outlet_advective_component_outflow_mol:'累计出口对流组分流出量',cumulative_component_inflow_mol:'累计组分流入量',cumulative_component_outflow_mol:'累计组分流出量',actual_liquid_inflow:'所有边界实际液体流入率',actual_liquid_outflow:'所有边界实际液体流出率',Pl_relative_min_node_id:'液相相对压力最小节点',Pl_relative_max_node_id:'液相相对压力最大节点',Pg_relative_min_Pa:'气相相对压力最小值',Pg_relative_max_Pa:'气相相对压力最大值',Pg_relative_min_node_id:'气相相对压力最小节点',Pg_relative_max_node_id:'气相相对压力最大节点',dt_change_reason:'时间步调整原因',min_free_gas_moles:'单孔最小自由气体摩尔数',gas_timescale_cap_s:'气体时间尺度步长上限',concentration_timescale_cap_s:'浓度变化步长上限',gas_disappearance_event:'本步是否发生气泡消失',event_node_id:'气泡消失事件节点',event_time_s:'气泡消失事件物理时间',event_dt_s:'气泡消失事件步长',event_gas_moles_before:'事件前自由气体摩尔数',event_gas_moles_after:'事件后自由气体摩尔数',concentration_limit_node_id:'浓度限制控制节点',concentration_limit_C_mol_m3:'限制节点当前浓度',concentration_limit_dC_dt_mol_m3_s:'限制节点浓度变化率',concentration_reference_mol_m3:'浓度变化参考尺度',concentration_tau_s:'局部浓度特征时间',termination_reason:'终止原因'};
Object.assign(diagnosticsZh,{retired_bubble_count:'本步退休气泡数 / Retired bubbles this step',batched_retired_bubble_count:'本步批量退休气泡数 / Batched retirements this step',batched_retired_moles:'批量转入溶解相摩尔数 / Batched moles transferred',batched_retired_gas_volume_m3:'批量退休气体体积 / Batched retired gas volume',batch_aggregate_free_gas_fraction:'本步批量自由气体占比 / Aggregate batch gas fraction',event_timestep_localized:'本步是否执行精确事件定位 / Exact event localization used'});
Object.assign(fieldInfo,{'run ID':['运行编号 / Run ID','输出目录对应的运行标识'],output:['输出目录 / Output directory','本次运行的结果目录'],config:['配置文件 / Configuration file','本次运行实际使用的JSON配置'],'config SHA-256':['配置文件哈希 / Config SHA-256','配置文件内容哈希'],'source SHA-256':['源码聚合哈希 / Source tree SHA-256','构建时源码树哈希'],'executable SHA-256':['可执行文件哈希 / Executable SHA-256','求解器二进制哈希'],checkpoint:['检查点格式 / Checkpoint format','重启文件格式'],'checkpoint step/time':['检查点步数与时间 / Checkpoint step and time','最近检查点状态'],'checkpoint SHA-256':['检查点哈希 / Checkpoint SHA-256','检查点文件哈希'],'开始时间（上海）':['开始时间 / Started time','上海时区启动时间'],'真实墙钟耗时':['真实墙钟耗时 / Wall-clock elapsed time','实际计算耗时'],'最新数据时间（上海）':['最新数据时间 / Latest data time','上海时区最后数据时间'],termination:['终止状态 / Termination status','停止或完成原因'],'image dimensions':['原始图像尺寸 / Image dimensions','PNM来源图像体素尺寸'],'voxel size':['体素边长 / Voxel size','单个体素物理边长'],'flow axis':['流动方向 / Flow axis','样品主流动轴'],nodes:['节点数 / Node count','控制体总数'],edges:['孔喉边数 / Edge count','计算网络边数'],'virtual inlet nodes':['虚拟入口节点数 / Virtual inlet nodes','入口虚拟节点数'],'virtual outlet nodes':['虚拟出口节点数 / Virtual outlet nodes','出口虚拟节点数'],'boundary edges':['边界连接数 / Boundary edges','样品与虚拟储层连接数'],'control volume total':['真实控制体总体积 / Total control-volume volume','排除虚拟节点后的体积'],'main pore volume':['主贯通域孔隙体积 / Main connected pore volume','入口出口主连通域体积'],'initial bubbles':['初始气泡数 / Initial bubble count','初始含自由气体控制体数'],'active bubbles':['当前活跃气泡数 / Active bubble count','当前含自由气体控制体数'],segments:['分段数量 / Segment count','归一化x分段数'],'liquid / gas':['液体与气体组分 / Liquid and gas','算例流体组分'],'background pressure':['背景绝对压力 / Background absolute pressure','相对压力基准'],diffusivity:['分子扩散系数 / Molecular diffusivity','溶解气体扩散系数'],kL:['液膜传质系数 / Mass-transfer coefficient kL','倍率施加前传质系数'],'initial concentration':['初始溶解浓度 / Initial dissolved concentration','样品初始液相浓度'],interpretation:['物理解释范围 / Interpretation scope','结果适用范围'],'pressure BC':['压力边界条件 / Pressure boundary condition','入口出口压力条件'],'concentration BC':['浓度边界条件 / Concentration boundary condition','入口与回流浓度规则'],'diffusive BC':['扩散边界条件 / Diffusive boundary condition','边界扩散通量规则'],'adaptive dt':['自适应时间步 / Adaptive time stepping','是否启用自适应步长'],'dt initial':['初始时间步长 / Initial dt','运行起始步长'],'dt current':['当前时间步长 / Current dt','最近接受步长'],'dt min/max':['时间步上下限 / Minimum and maximum dt','允许步长范围'],'growth/shrink':['步长增减因子 / Growth and shrink factors','步长调整倍率'],'Newton tolerance':['Newton收敛容差 / Newton tolerances','残量与更新双门槛'],'checkpoint interval':['检查点间隔 / Checkpoint interval','检查点保存间隔'],'VTK interval':['VTK输出间隔 / VTK output interval','状态文件输出间隔'],'residual free gas fraction':['残余自由气体比例 / Residual free-gas fraction','自由气体终止阈值']});
const info=k=>fieldInfo[k]||[diagnosticsZh[k]||String(k).replaceAll('_',' '),diagnosticsZh[k]?`原始字段：${k}`:'尚未提供中文说明'];
const groupStarts={step:'时间步与非线性收敛 / Time stepping and nonlinear convergence',boundary_liquid_outflow_m3_s:'守恒与气体库存 / Conservation and gas inventory',boundary_advective_component_outflow:'边界流动与组分输运 / Boundary flow and component transport',internal_pressure_min_Pa:'压力范围与气泡状态 / Pressure range and bubble state',next_dt_s:'自适应时间步与消失事件 / Adaptive time stepping and disappearance events',physical_time_s:'科学统计与运行配置 / Scientific statistics and runtime configuration'};
function runtimeRows(d){let out='';for(const[k,v]of Object.entries(d)){if(groupStarts[k])out+=`<tr class="data-group"><th colspan="3">${groupStarts[k]}</th></tr>`;out+=`<tr data-tooltip="${esc(info(k)[1])}"><td>${esc(label(k))}</td><td>${esc(fmt(v,k))}</td><td>诊断记录 diagnostics.tsv</td></tr>`;}return out;}
const label=k=>`${info(k)[0]} · ${k}`;
const kv=(k,v,source='runtime')=>`<div class="kv" data-tooltip="${esc(info(k)[1])}"><span class="k">${esc(label(k))}</span><br><span class="v">${esc(fmt(v,k))}</span> <small>(${esc(sourceZh(source))})</small></div>`;
const colors=['var(--blue)','var(--orange)','var(--purple)','var(--green)'];
function drawChart(id,rows,defs,yLabel){const box=document.getElementById(id);if(!box)return;const svg=box.querySelector('svg');if(!rows.length){box.innerHTML='<div class="chart-empty">暂无 global_history.tsv 数据</div>';return;}const w=Math.max(360,box.clientWidth||640),h=220,L=58,R=16,T=16,B=32;svg.setAttribute('viewBox',`0 0 ${w} ${h}`);svg.innerHTML='';const xs=rows.map(r=>r.time_s).filter(Number.isFinite),vals=defs.flatMap(d=>rows.map(r=>r[d.key]).filter(Number.isFinite));let xmin=Math.min(...xs),xmax=Math.max(...xs);if(xmin===xmax)xmax=xmin||1;let ymin=Math.min(...vals),ymax=Math.max(...vals);if(ymin===ymax){const e=Math.max(Math.abs(ymin)*.05,1);ymin-=e;ymax+=e;}const pad=(ymax-ymin)*.08;ymin-=pad;ymax+=pad;const X=x=>L+(x-xmin)/(xmax-xmin)*(w-L-R),Y=y=>h-B-(y-ymin)/(ymax-ymin)*(h-T-B),ns='http://www.w3.org/2000/svg';const line=(tag,a)=>{const e=document.createElementNS(ns,tag);for(const[k,v]of Object.entries(a))e.setAttribute(k,v);svg.appendChild(e);return e;};for(let i=0;i<5;i++){const y=T+i*(h-T-B)/4;line('line',{x1:L,x2:w-R,y1:y,y2:y,class:'chart-grid'});const t=line('text',{x:L-8,y:y+4,'text-anchor':'end',class:'chart-label'});t.textContent=(ymax-(y-T)/(h-T-B)*(ymax-ymin)).toExponential(2);}line('line',{x1:L,x2:w-R,y1:h-B,y2:h-B,class:'chart-axis'});line('line',{x1:L,x2:L,y1:T,y2:h-B,class:'chart-axis'});for(const[i,x]of[xmin,(xmin+xmax)/2,xmax].entries()){const t=line('text',{x:X(x),y:h-10,'text-anchor':i===0?'start':i===2?'end':'middle',class:'chart-label'});t.textContent=x.toPrecision(4)+' s';}const yl=line('text',{x:14,y:(h+T-B)/2,transform:`rotate(-90 14 ${(h+T-B)/2})`,'text-anchor':'middle',class:'chart-label'});yl.textContent=yLabel;for(const[i,d]of defs.entries()){const pts=rows.filter(r=>Number.isFinite(r[d.key]));if(!pts.length)continue;line('polyline',{points:pts.map(r=>`${X(r.time_s)},${Y(r[d.key])}`).join(' '),class:'chart-line',stroke:colors[i%colors.length]});}box.querySelector('.chart-legend').innerHTML=defs.map((d,i)=>`<span><i style="background:${colors[i%colors.length]}"></i>${esc(d.label)}</span>`).join('');}
function markEvents(id,events,rows){const box=document.getElementById(id),svg=box&&box.querySelector('svg');if(!box||!svg||!events.length||!rows.length)return;const frame=svg.viewBox.baseVal,labels=box.querySelector('.chart-legend'),batch=events.filter(e=>e.event_mode==='batched_microbubble'),exact=events.filter(e=>e.event_mode!=='batched_microbubble');labels.innerHTML+=`<span><i style="background:var(--purple)"></i>批量退休 ${batch.length} / Batched</span><span><i style="background:var(--red)"></i>精确定位 ${exact.length} / Exact</span>`;const ns='http://www.w3.org/2000/svg',times=rows.map(e=>e.time_s),maxTime=Math.max(...times),minTime=Math.min(...times),visible=events.filter(e=>e.time_s>=minTime&&e.time_s<=maxTime),stride=Math.max(1,Math.ceil(visible.length/40));for(const e of visible.filter((_,i)=>i%stride===0)){const batched=e.event_mode==='batched_microbubble',color=batched?'var(--purple)':'var(--red)',x=58+(e.time_s-minTime)/Math.max(maxTime-minTime,1e-30)*(frame.width-74),line=document.createElementNS(ns,'line');line.setAttribute('x1',x);line.setAttribute('x2',x);line.setAttribute('y1',16);line.setAttribute('y2',frame.height-32);line.setAttribute('class','event-line');line.setAttribute('stroke',color);line.setAttribute('stroke-opacity','.58');line.setAttribute('stroke-dasharray',batched?'2 4':'6 3');const title=document.createElementNS(ns,'title');title.textContent=`${batched?'微小气泡批量退休 / Batched retirement':'精确消失事件 / Exact event'} · node ${e.event_node_id} · t=${Number(e.time_s).toPrecision(6)} s`;line.appendChild(title);svg.appendChild(line);const dot=document.createElementNS(ns,'circle');dot.setAttribute('cx',x);dot.setAttribute('cy',18);dot.setAttribute('r',batched?3:5);dot.setAttribute('fill',color);dot.setAttribute('class','event-dot');svg.appendChild(dot);}}
function addProfileTicks(c){const svg=document.querySelector('#chart-profile svg'),rows=c.segment_spacetime||[];if(!svg||!rows.length)return;const w=svg.viewBox.baseVal.width,h=svg.viewBox.baseVal.height,L=50,R=16,T=16,B=32,ns='http://www.w3.org/2000/svg';const values=rows.flatMap(r=>[r.segment_mean_Sg,r.initial_segment_mean_Sg]).filter(Number.isFinite),ymax=Math.max(...values,1e-12)*1.12;for(let i=0;i<5;i++){const y=T+i*(h-T-B)/4,t=document.createElementNS(ns,'text');t.setAttribute('x',L-7);t.setAttribute('y',y+4);t.setAttribute('text-anchor','end');t.setAttribute('class','chart-label');t.textContent=(ymax*(1-i/4)).toFixed(3);svg.appendChild(t);}for(const [i,sid] of [rows[0].segment_id,rows[Math.floor((rows.length-1)/2)].segment_id,rows[rows.length-1].segment_id].entries()){const x=L+(sid-1)/Math.max(rows.length-1,1)*(w-L-R),t=document.createElementNS(ns,'text');t.setAttribute('x',x);t.setAttribute('y',h-18);t.setAttribute('text-anchor',i===0?'start':i===2?'end':'middle');t.setAttribute('class','chart-label');t.textContent=String(sid);svg.appendChild(t);}}
function fillSegmentTable(c){const table=document.getElementById('segment-table-body');if(!table)return;table.innerHTML=(c.segment_spacetime||[]).map(r=>`<tr><td>${fmt(r.segment_id,'id')}</td><td>${fmt(r.initial_segment_mean_Sg,'Sg')}</td><td>${fmt(r.segment_mean_Sg,'Sg')}</td><td>${fmt(r.segment_Sg_change,'Sg')}</td><td>${fmt(r.initial_segment_bubble_count,'count')}</td><td>${fmt(r.segment_active_bubble_count,'count')}</td><td>${fmt(r.initial_segment_gas_volume_m3,'volume')}</td><td>${fmt(r.segment_gas_volume_m3,'volume')}</td><td>${fmt(r.event_count,'count')}</td><td>${fmt(r.last_event_time_s,'time_s')}</td></tr>`).join('');}
function drawCharts(c){const r=c.history_series||[],q=c.conservation_series||[],events=c.disappearance_events||[],batch=events.filter(e=>e.event_mode==='batched_microbubble'),exact=events.filter(e=>e.event_mode!=='batched_microbubble'),summary=document.getElementById('event-summary');if(summary)summary.innerHTML=`<span class="event-pill">微小气泡批量退休 / Batched: <b>${batch.length}</b></span><span class="event-pill exact">精确事件定位 / Exact: <b>${exact.length}</b></span><span class="event-pill">最近事件 / Last: <b>${events.length?fmt(events[events.length-1].time_s,'time_s')+' s':'—'}</b></span>`;drawChart('chart-inventory',r,[{key:'free_gas_mol',label:'自由气体 / Free gas'},{key:'dissolved_gas_mol',label:'溶解气体 / Dissolved gas'}],'摩尔数 / Moles [mol]');drawChart('chart-saturation',r,[{key:'mean_Sg',label:'总体 Sg / Overall Sg'},{key:'main_real_gas_saturation',label:'主贯通域 Sg / Main-domain Sg'}],'气体饱和度 / Gas saturation');drawChart('chart-timestep',r,[{key:'dt_s',label:'已接受 dt / Accepted dt'}],'时间步长 / Time step [s]');drawChart('chart-quality',q,[{key:'epsilon_N_adjusted',label:'摩尔守恒误差 / Molar balance εN'},{key:'epsilon_V_adjusted',label:'体积守恒误差 / Volume balance εV'}],'相对守恒误差 / Relative balance error');drawProfile(c);fillSegmentTable(c);addProfileTicks(c);for(const id of ['chart-inventory','chart-saturation','chart-timestep'])markEvents(id,events,r);}
function drawProfile(c){const box=document.getElementById('chart-profile'),rows=c.segment_spacetime||c.segment_profile||[];if(!box||!rows.length){if(box)box.innerHTML='<div class="chart-empty">暂无分段历史数据</div>';return;}const svg=box.querySelector('svg'),w=Math.max(360,box.clientWidth||640),h=220,L=50,R=16,T=16,B=32;svg.setAttribute('viewBox',`0 0 ${w} ${h}`);svg.innerHTML='';const allSg=rows.flatMap(r=>[r.segment_mean_Sg,r.initial_segment_mean_Sg]).filter(Number.isFinite),ymax=Math.max(...allSg,1e-12)*1.12,X=i=>L+(i-1)/Math.max(rows.length-1,1)*(w-L-R),Y=v=>h-B-v/ymax*(h-T-B),ns='http://www.w3.org/2000/svg';for(let i=0;i<5;i++){const y=T+i*(h-T-B)/4,e=document.createElementNS(ns,'line');e.setAttribute('x1',L);e.setAttribute('x2',w-R);e.setAttribute('y1',y);e.setAttribute('y2',y);e.setAttribute('class','chart-grid');svg.appendChild(e);}for(const [key,color] of [['initial_segment_mean_Sg',colors[1]],['segment_mean_Sg',colors[0]]]){const path=document.createElementNS(ns,'polyline');path.setAttribute('points',rows.filter(r=>Number.isFinite(r[key])).map(r=>`${X(r.segment_id)},${Y(r[key])}`).join(' '));path.setAttribute('class','chart-line');path.setAttribute('stroke',color);svg.appendChild(path);}const yl=document.createElementNS(ns,'text');yl.setAttribute('x',12);yl.setAttribute('y',(h+T-B)/2);yl.setAttribute('transform',`rotate(-90 12 ${(h+T-B)/2})`);yl.setAttribute('text-anchor','middle');yl.setAttribute('class','chart-label');yl.textContent='分段气体饱和度 Sg';svg.appendChild(yl);const xl=document.createElementNS(ns,'text');xl.setAttribute('x',w/2);xl.setAttribute('y',h-6);xl.setAttribute('text-anchor','middle');xl.setAttribute('class','chart-label');xl.textContent='沿流向归一化区段';svg.appendChild(xl);box.querySelector('.chart-legend').innerHTML=`<span><i style="background:${colors[1]}"></i>初始 Sg</span><span><i style="background:${colors[0]}"></i>当前 Sg</span><span>${rows.length} 个区段</span>`;const table=document.getElementById('segment-table-body');if(table)table.innerHTML=rows.map(r=>`<tr><td>${fmt(r.segment_id,'id')}</td><td>${fmt(r.initial_segment_mean_Sg,'Sg')}</td><td>${fmt(r.segment_mean_Sg,'Sg')}</td><td>${fmt(r.segment_Sg_change,'Sg')}</td><td>${fmt(r.segment_active_bubble_count,'count')}</td><td>${fmt(r.segment_gas_volume_m3,'volume')}</td><td>${fmt(r.event_count,'count')}</td><td>${fmt(r.last_event_time_s,'time_s')}</td></tr>`).join('');}
function updateRuntime(c){const el=document.getElementById('runtime-strip');if(!el)return;const d=c.last_diagnostics||{}, running=c.status==='running';const t=d.physical_time_s||d.time_s||c.checkpoint?.time_s||'—';const pct=Math.min(100,Math.max(0,(Number(t)/(Number(c.effective?.maximum_physical_time_s)||1))*100));const age=Number(c.diagnostics_age_seconds);const ageText=Number.isFinite(age)?`${Math.max(0,Math.round(age))} 秒前`:'未知';el.innerHTML=`<div class="runtime-item"><span class="label">求解器状态</span><span class="value"><span class="live-dot ${running?'on':''}"></span>${esc(c.status)}</span></div><div class="runtime-item"><span class="label">已接受时间步</span><span class="value">${esc(fmt(d.step||c.checkpoint?.step_number||'—','step'))}</span></div><div class="runtime-item"><span class="label">物理模拟时间</span><span class="value">${esc(fmt(t,'time_s'))} s</span><div class="progress-track"><div class="progress-bar" style="width:${pct}%"></div></div></div><div class="runtime-item"><span class="label">真实墙钟耗时</span><span class="value">${esc(fmt(c.wall_clock_runtime_seconds,'time_s'))} s</span></div><div class="runtime-item"><span class="label">接受步长 / Newton迭代</span><span class="value">${esc(fmt(d.dt_s||'—','dt_s'))} s / ${esc(fmt(d.newton_iterations||'—','iterations'))}</span></div><div class="runtime-item"><span class="label">最新接受数据（上海时间）</span><span class="value">${esc(shanghaiTime(c.last_update))}</span><span class="context">${esc(ageText)}</span></div><div class="runtime-item"><span class="label">面板刷新时间（上海时间）</span><span class="value">${esc(shanghaiTime(new Date()))}</span><span class="context">Asia/Shanghai · UTC+8</span></div>`;}
function ensureCharts(){if(document.getElementById('chart-inventory'))return;const s='<div class="card"><h3>Runtime curves</h3><div id="event-summary" class="event-summary"></div><div class="charts"><div class="chart-wrap" id="chart-inventory"><div class="chart-title">Gas inventory</div><svg role="img" aria-label="Free and dissolved gas over time"><title>Gas inventory</title></svg><div class="chart-legend"></div></div><div class="chart-wrap" id="chart-saturation"><div class="chart-title">Gas saturation</div><svg role="img" aria-label="Gas saturation over time"><title>Gas saturation</title></svg><div class="chart-legend"></div></div><div class="chart-wrap" id="chart-timestep"><div class="chart-title">Accepted time step</div><svg role="img" aria-label="Accepted time step over time"><title>Accepted time step</title></svg><div class="chart-legend"></div></div><div class="chart-wrap" id="chart-quality"><div class="chart-title">Conservation error</div><svg role="img" aria-label="Mass and volume conservation error over time"><title>Conservation error</title></svg><div class="chart-legend"></div></div></div><h3>Flow-direction segment profile</h3><div class="chart-wrap" id="chart-profile"><div class="chart-title">Latest segment gas saturation</div><svg role="img" aria-label="Latest gas saturation by normalized x segment"><title>Latest segment gas saturation</title></svg><div class="chart-legend"></div></div><div class="table-wrap"><table><thead><tr><th>Segment</th><th>Current Sg</th><th>Active bubbles</th><th>Gas volume [m³]</th><th>Disappearance events</th><th>Last event time [s]</th></tr></thead><tbody id="segment-table-body"></tbody></table></div></div>';document.getElementById('detail').insertAdjacentHTML('afterbegin',s);}
async function controlPost(path,body){const result=document.getElementById('control-result');result.textContent='正在提交 / Submitting…';try{const response=await fetch(path,{method:'POST',headers:{'Content-Type':'application/json','X-Control-Token':controlToken},body:JSON.stringify(body)});const payload=await response.json();result.textContent=payload.ok?'操作成功 / Success\n'+JSON.stringify(payload.result,null,2):'操作失败 / Failed: '+payload.error;result.className=payload.ok?'control-result ok':'control-result warn';return payload;}catch(error){result.textContent='请求失败 / Request failed: '+error;result.className='control-result warn';}}
function ensureControl(c){if(!c.control_allowed||!data?.control?.enabled||document.getElementById('case-control'))return;const running=c.status==='running',q=c.effective||{},io=c.config?.io||{};const html=`<div class="card" id="case-control"><h3>运行控制 / Run control</h3><p class="control-note">仅允许安全停止、checkpoint恢复和低风险运行参数修订。批量调度算例在此面板中保持只读，避免绕过队列锁。</p><div class="control-grid"><label>控制令牌 / Control token<input id="control-token" type="password" autocomplete="off" value="${esc(controlToken)}"></label><label>最大时间步 / dt_max_s [s]<input id="ctl-dtmax" type="number" min="0.000001" max="1000000" step="any" value="${esc(q.dt_max_s??'')}"></label><label>检查点间隔 / checkpoint_interval_steps<input id="ctl-checkpoint" type="number" min="1" step="1" value="${esc(io.checkpoint_interval_steps??'')}"></label><label>VTK间隔 / vtk_interval_steps<input id="ctl-vtk" type="number" min="1" step="1" value="${esc(io.vtk_interval_steps??'')}"></label></div><div class="control-actions"><button id="ctl-stop" ${running?'':'disabled'}>安全停止 / Safe stop</button><button id="ctl-stage">暂存参数 / Stage settings</button><button id="ctl-restart" ${running?'disabled':''}>从检查点恢复 / Restart</button><button disabled>复制新物理算例 / Fork new physical case（尚未开放）</button></div><pre id="control-result" class="control-result">控制模式已启用；操作将写入 control_history.jsonl。</pre></div>`;document.getElementById('detail').insertAdjacentHTML('afterbegin',html);document.getElementById('control-token').oninput=e=>controlToken=e.target.value;document.getElementById('ctl-stop').onclick=async()=>{controlToken=document.getElementById('control-token').value;if(confirm(`确认安全停止 ${c.case_name}？\nConfirm safe stop?`))await controlPost('/api/control/stop',{case_name:c.run_id});};document.getElementById('ctl-stage').onclick=async()=>{controlToken=document.getElementById('control-token').value;const changes={dt_max_s:Number(document.getElementById('ctl-dtmax').value),checkpoint_interval_steps:Number(document.getElementById('ctl-checkpoint').value),vtk_interval_steps:Number(document.getElementById('ctl-vtk').value)};if(confirm(`参数仅在下次restart生效。\nStage settings for ${c.case_name}?`))await controlPost('/api/control/stage',{case_name:c.run_id,changes});};document.getElementById('ctl-restart').onclick=async()=>{controlToken=document.getElementById('control-token').value;if(confirm(`从最新checkpoint恢复 ${c.case_name}？\nRestart from latest checkpoint?`))await controlPost('/api/control/restart',{case_name:c.run_id});};}
function consoleClass(line){if(/(?:ERROR|FAILED|NOT_CONVERGED|NaN|Inf|dt_min_failure)/i.test(line))return'error';if(/(?:rejected_dt_s=|retries=[1-9]|retry)/i.test(line))return'retry';if(/(?:gas_disappearance_event=1|event_node_id=[1-9]|disappearance|batched_retired_bubbles=[1-9])/i.test(line))return'event';if(/^step=\d+/.test(line))return'step';if(/^(?:started_utc|config=|nodes=|elapsed_seconds=)/.test(line))return'meta';return'';}
function ensureConsole(c){if(document.getElementById('solver-console'))return;const html=`<div class="card console-card" id="solver-console"><div class="console-head"><div class="console-title"><span class="live-dot ${c.status==='running'?'on':''}"></span>实时求解终端 / Live tmux solver console</div><div class="console-actions"><button id="console-pause">${consolePaused?'继续 / Resume':'暂停 / Pause'}</button><button id="console-scroll">${consoleAutoScroll?'自动滚动：开 / Auto-scroll: on':'自动滚动：关 / Auto-scroll: off'}</button></div></div><pre class="console-screen" id="console-screen" aria-live="polite"><span class="console-empty">正在连接专用 tmux live 窗口… / Connecting…</span></pre><div class="console-foot"><span id="console-source">只读取受限的 run.log 投射窗格</span><span id="console-updated">—</span></div></div>`;document.getElementById('detail').insertAdjacentHTML('afterbegin',html);document.getElementById('console-pause').onclick=()=>{consolePaused=!consolePaused;document.getElementById('console-pause').textContent=consolePaused?'继续 / Resume':'暂停 / Pause';if(!consolePaused)refreshConsole();};document.getElementById('console-scroll').onclick=()=>{consoleAutoScroll=!consoleAutoScroll;document.getElementById('console-scroll').textContent=consoleAutoScroll?'自动滚动：开 / Auto-scroll: on':'自动滚动：关 / Auto-scroll: off';};refreshConsole();}
async function refreshConsole(){if(consolePaused||!selectedRelative||!data)return;const requestDirectory=selectedDirectory;if(!selectedCase||!document.getElementById('console-screen'))return;try{const response=await fetch(`/api/console?case=${encodeURIComponent(selectedRelative)}&lines=100`,{cache:'no-store'});const payload=await response.json();if(!response.ok)throw new Error(payload.error||response.statusText);if(selectedDirectory!==requestDirectory)return;const screen=document.getElementById('console-screen'),source=document.getElementById('console-source'),updated=document.getElementById('console-updated');if(!screen)return;const lines=payload.lines||[];screen.innerHTML=lines.length?lines.map(line=>`<span class="console-line ${consoleClass(line)}">${esc(line)}</span>`).join('')+`<span class="console-cursor" title="等待下一条求解输出"></span>`:'<span class="console-empty">tmux live窗口尚无输出 / Waiting for solver output…</span>';if(source)source.textContent=`来源 / Source: ${payload.source||'run.log'} · ${fmt(payload.size_bytes||0)} bytes`;if(updated)updated.textContent=`更新 / Updated: ${shanghaiTime(payload.updated_utc)} · ${fmt(payload.age_seconds||0)} s ago`;if(consoleAutoScroll)screen.scrollTop=screen.scrollHeight;}catch(error){const screen=document.getElementById('console-screen');if(screen&&selectedDirectory===requestDirectory)screen.innerHTML=`<span class="console-line error">终端投射暂不可用 / Console unavailable: ${esc(error)}</span>`;}}
function ensureRuntime(){if(document.getElementById('runtime-strip'))return;document.getElementById('detail').insertAdjacentHTML('afterbegin','<div class="card"><h3>Solver heartbeat</h3><div id="runtime-strip" class="runtime-strip" aria-live="polite"></div></div>');}
function render(c){selectedDirectory=c.directory;selectedRelative=c.relative_path;selectedCase=c;document.querySelectorAll('button.case').forEach(b=>b.classList.toggle('sel',b.dataset.path===c.relative_path)); let m=c.manifest||{}, p=c.physics||{}, q=c.effective||{}, d=c.last_diagnostics||{}, cp=c.checkpoint||{}, g=c.geometry_network||{}; const warn=(c.warnings||[]).map(x=>`<div class="warn">⚠ ${esc(x)}</div>`).join(''); const charts='<div class="card"><h3>Runtime curves</h3><div class="charts"><div class="chart-wrap" id="chart-inventory"><div class="chart-title">Gas inventory</div><svg role="img" aria-label="Free and dissolved gas over time"><title>Gas inventory</title></svg><div class="chart-legend"></div></div><div class="chart-wrap" id="chart-saturation"><div class="chart-title">Gas saturation</div><svg role="img" aria-label="Gas saturation over time"><title>Gas saturation</title></svg><div class="chart-legend"></div></div><div class="chart-wrap" id="chart-timestep"><div class="chart-title">Accepted time step</div><svg role="img" aria-label="Accepted time step over time"><title>Accepted time step</title></svg><div class="chart-legend"></div></div><div class="chart-wrap" id="chart-quality"><div class="chart-title">Conservation error</div><svg role="img" aria-label="Mass and volume conservation error over time"><title>Conservation error</title></svg><div class="chart-legend"></div></div></div></div>';
document.getElementById('detail').innerHTML=`<div class="card"><h2>${esc(c.case_name)}</h2><span class="badge ${esc(c.status)}">${esc(statusBi(c.status))}</span><span class="badge ${esc(c.kind)}">${esc(kindBi(c.kind))}</span><span class="badge">后端 / Backend：${esc(q.backend)}</span><p>${esc(c.description)}</p><p><b>科学解释 / Scientific interpretation：</b> ${esc(c.scientific_interpretation)}</p>${warn}<div class="grid">${kv('status',statusBi(c.status))}${kv('run ID',c.run_id)}${kv('output',c.directory)}${kv('config',c.config_path)}${kv('config SHA-256',c.config_sha256)}${kv('source SHA-256',c.source_sha256)}${kv('executable SHA-256',c.executable_sha256)}${kv('checkpoint',cp.format)}${kv('checkpoint step/time',`${cp.step_number??'—'} / ${cp.time_s??'—'} s`)}${kv('checkpoint SHA-256',c.checkpoint_sha256)}${kv('开始时间（上海）',shanghaiTime(c.started_utc),'manifest / run.log')}${kv('真实墙钟耗时',`${fmt(c.wall_clock_runtime_seconds,'time_s')} s`)}${kv('最新数据时间（上海）',shanghaiTime(c.last_update),'diagnostics.tsv')}${kv('termination',m.termination_reason||m.termination?.actual_reason||c.termination_events?.termination_reason)}</div></div>
<div class="card"><h3>Geometry & network</h3><div class="grid">${kv('image dimensions',g.image_dimensions_voxels,'config')} ${kv('voxel size',`${g.voxel_size_um??'—'} µm`,'config')} ${kv('flow axis',g.flow_axis,'config')} ${kv('nodes',g.nodes,'manifest/run.log')} ${kv('edges',g.edges,'manifest/run.log')} ${kv('virtual inlet nodes',g.inlet_virtual_nodes,'geometry/config')} ${kv('virtual outlet nodes',g.outlet_virtual_nodes,'geometry/config')} ${kv('boundary edges',g.boundary_edges,'manifest/run.log')} ${kv('control volume total',`${g.control_volume_total_m3??'—'} m³`,'segment_geometry.tsv')} ${kv('main pore volume',`${g.main_real_pore_volume_m3??'—'} m³`,'global_history.tsv')} ${kv('initial bubbles',g.initial_bubbles,'physical_parameter_audit')} ${kv('active bubbles',g.active_bubbles,'diagnostics.tsv')} ${kv('segments',g.segments,'config')}</div></div>
<div class="card"><h3>Physics</h3><div class="grid">${kv('liquid / gas',`${c.config?.case?.liquid_composition||'pure_water'} / ${c.config?.case?.gas_component||'air'}`)}${kv('temperature',`${p.T_K??'—'} K`)}${kv('background pressure',`${c.config?.pressure_reference?.background_abs_Pa??'—'} Pa`)}${kv('Henry coefficient',`${p.Hcp_mol_m3_Pa??'—'} mol m⁻³ Pa⁻¹`)}${kv('diffusivity',`${p.Dm_m2_s??'—'} m² s⁻¹`)}${kv('kL',`${p.kL_m_s??'—'} m s⁻¹`)}${kv('mass_transfer_multiplier',p.mass_transfer_multiplier,'config/preview-only')}${kv('initial concentration',`${c.config?.initial?.C_initial_mol_m3??'—'} mol m⁻³`)}${kv('interpretation','unconfirmed closures / preview-only','physical_parameter_audit')}</div></div>
<div class="card"><h3>Boundary & numerical parameters</h3><div class="grid">${kv('pressure BC',`${c.config?.boundary?.Pl_in_rel_Pa??'—'} / ${c.config?.boundary?.Pl_out_rel_Pa??'—'} Pa relative`)}${kv('concentration BC',c.config?.boundary?.species_mode)}${kv('diffusive BC',c.config?.boundary?.outlet_diffusive_condition)}${kv('adaptive dt',q.adaptive_dt_enabled,'config/runtime')}${kv('dt initial',q.dt_initial_s,'config')}${kv('dt current',q.current_dt_s,'checkpoint/diagnostics')}${kv('dt min/max',`${q.dt_min_s} / ${q.dt_max_s}`,'config')}${kv('growth/shrink',`${q.growth_factor} / ${q.shrink_factor}`,'config')}${kv('Newton tolerance',`${q.residual_tolerance} / ${q.update_tolerance}`,'config')}${kv('checkpoint interval',q.checkpoint_interval_steps,'config')}${kv('VTK interval',q.vtk_interval_steps,'config')}${kv('residual free gas fraction',q.residual_free_gas_fraction,'config/termination')}</div></div>
<div class="card"><h3>配置值与实际运行值</h3><table><tr><th>参数（中文解释 / 原字段）</th><th>配置值</th><th>实际值</th><th>单位 / 来源</th><th>可信度分类</th></tr>${(c.parameters||[]).map(x=>`<tr class="${x.conflict?'warn':''}" data-tooltip="${esc(info(x.name)[1])}"><td>${esc(label(x.name))}</td><td>${esc(fmt(x.configured,x.name))}</td><td>${esc(fmt(x.effective,x.name))}</td><td>${esc(x.unit)} / ${esc(sourceZh(x.source))}</td><td>${esc(classZh(x.classification))}${x.conflict?' ⚠ 数值冲突':''}</td></tr>`).join('')}</table></div>
<div class="card"><h3>Geometry, boundary and provenance</h3><div class="grid">${Object.entries(c.geometry_network||{}).map(([k,v])=>kv(k,v,'geometry audit / manifest / diagnostics')).join('')}${Object.entries(c.boundary_runtime||{}).map(([k,v])=>kv(k,v,'config / diagnostics')).join('')}</div><h4>PNM inputs</h4><table><tr><th>role</th><th>path</th><th>SHA-256</th></tr>${(c.pnm_inputs||[]).map(x=>`<tr><td>${esc(x.role)}</td><td>${esc(x.path)}</td><td>${esc(x.sha256)}</td></tr>`).join('')}</table></div>
<div class="card"><h3>最新运行数据</h3><table><tr><th>参数（中文解释 / 原字段）</th><th>数值</th><th>数据来源</th></tr>${runtimeRows(d)}</table></div>`;}
const originalRender=render; render=(c)=>{originalRender(c);ensureControl(c);ensureRuntime();ensureConsole(c);updateRuntime(c);ensureCharts();drawCharts(c);translateUI();};
async function loadDetail(brief){selectedRelative=brief.relative_path;selectedDirectory=brief.directory;document.getElementById('detail').innerHTML='<div class="card">正在读取算例详情 / Loading case details…</div>';try{const response=await fetch(`/api/case?case=${encodeURIComponent(brief.relative_path)}`,{cache:'no-store'});const full=await response.json();if(!response.ok)throw new Error(full.error||response.statusText);if(selectedRelative===brief.relative_path)render(full);}catch(error){document.getElementById('detail').innerHTML=`<div class="card warn">详情读取失败 / Detail unavailable: ${esc(error)}</div>`;}}
function renderCaseList(){if(!data)return;const search=document.getElementById('case-search').value.toLowerCase(),variant=document.getElementById('variant-filter').value,status=document.getElementById('status-filter').value;const visible=data.cases.filter(c=>(!search||c.case_name.toLowerCase().includes(search))&&(!variant||c.variant===variant)&&(!status||c.status===status));document.getElementById('list-count').textContent=`显示 ${visible.length} / ${data.cases.length} · Showing ${visible.length} of ${data.cases.length}`;document.getElementById('cases').innerHTML=visible.map(c=>`<button class="case ${c.relative_path===selectedRelative?'sel':''}" data-path="${esc(c.relative_path)}"><span class="case-name">${esc(c.case_name)}</span><span class="case-meta"><span class="badge ${esc(c.status)}">${esc(statusBi(c.status))}</span>${esc(c.variant||'—')} · step=${esc(fmt(c.last_diagnostics?.step||c.checkpoint?.step_number||'—','step'))} · t=${esc(fmt(c.last_diagnostics?.physical_time_s||c.last_diagnostics?.time_s||c.checkpoint?.time_s||'—','time_s'))} s · dt=${esc(fmt(c.last_diagnostics?.dt_s||c.effective?.current_dt_s||'—','dt_s'))} s</span></button>`).join('')||'<div class="console-empty">没有匹配算例 / No matching cases</div>';document.querySelectorAll('button.case').forEach(button=>button.onclick=()=>{const brief=data.cases.find(c=>c.relative_path===button.dataset.path);if(brief)loadDetail(brief);});}
function renderSummary(){const s=data?.summary||{},c=s.status_counts||{},v=s.variant_status_counts||{},count=status=>(v['PNM-0.6']?.[status]||0)+(v['PNM-0.8']?.[status]||0),p06=v['PNM-0.6']?Object.values(v['PNM-0.6']).reduce((a,b)=>a+b,0):0,p08=v['PNM-0.8']?Object.values(v['PNM-0.8']).reduce((a,b)=>a+b,0):0,batchFailed=count('failed'),otherFailed=Math.max(0,(c.failed||0)-batchFailed);document.getElementById('batch-summary').innerHTML=`<div class="summary-cell"><b>${p06+p08}</b><span>批量总数 / Batch</span></div><div class="summary-cell"><b>${count('running')}</b><span>批量运行 / Running</span></div><div class="summary-cell"><b>${count('queued')}</b><span>批量排队 / Queued</span></div><div class="summary-cell"><b>${count('completed')}</b><span>批量完成 / Completed</span></div><div class="summary-cell"><b>${batchFailed}</b><span>批量失败 / Batch failed</span></div><div class="summary-cell"><b>${otherFailed}</b><span>其他历史失败 / Historical</span></div><div class="summary-cell"><b>${p06} + ${p08}</b><span>PNM-0.6 + PNM-0.8</span></div>`;}
async function refreshAll(initial=false){try{const response=await fetch('/api/cases',{cache:'no-store'});const x=await response.json();data=x;document.getElementById('root').textContent=`${x.production_root} · 更新 / Updated ${shanghaiTime(x.generated_utc)}`;renderSummary();renderCaseList();const brief=(selectedRelative&&x.cases.find(c=>c.relative_path===selectedRelative))||(initial&&x.cases.find(c=>c.status==='running'&&c.variant))||(initial&&x.cases[0]);if(brief&&!document.getElementById('case-control')?.contains(document.activeElement))await loadDetail(brief);}catch(error){document.getElementById('list-count').textContent=`刷新失败 / Refresh failed: ${error}`;}}
document.getElementById('case-search').oninput=renderCaseList;document.getElementById('variant-filter').onchange=renderCaseList;document.getElementById('status-filter').onchange=renderCaseList;refreshAll(true);setInterval(()=>refreshAll(false),10000);
setInterval(refreshConsole,2000);
</script>'''


class ReusableThreadingHTTPServer(ThreadingHTTPServer):
    allow_reuse_address = True


class Handler(BaseHTTPRequestHandler):
    root: Path = DEFAULT_ROOT
    control_enabled = False
    control_token = ""
    solver_binary = Path("/workspace/zz/build-bubble-flow-release/bubble_solver")
    control_lock = threading.Lock()

    def send_json(self, status: int, value: dict[str, Any]) -> None:
        payload = json.dumps(value, ensure_ascii=False).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        try:
            self.wfile.write(payload)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def authorized(self) -> bool:
        if not self.control_enabled:
            return False
        supplied = self.headers.get("X-Control-Token", "")
        if not supplied or not hmac.compare_digest(supplied, self.control_token):
            return False
        origin = self.headers.get("Origin")
        if origin:
            origin_host = urlparse(origin).hostname
            request_host = urlparse("//" + self.headers.get("Host", "")).hostname
            loopback = {"127.0.0.1", "localhost", "::1"}
            if origin_host not in loopback or request_host not in loopback:
                return False
        return True

    def do_GET(self):  # noqa: N802
        request = urlparse(self.path)
        if request.path == "/api/cases":
            value = snapshot(self.root)
            value["control"] = {"enabled": self.control_enabled,
                "solver_binary": str(self.solver_binary) if self.control_enabled else None,
                "solver_binary_sha256": sha256(self.solver_binary) if self.control_enabled else None,
                "allowed_settings": ["dt_max_s", "checkpoint_interval_steps", "vtk_interval_steps"]}
            if self.control_enabled:
                for case in value["cases"]:
                    if not case.get("control_allowed"):
                        case["control_state"] = {"read_only": True}
                        continue
                    try:
                        case["control_state"] = case_control.control_state(case["run_id"], self.root)
                    except case_control.ControlError as error:
                        case["control_state"] = {"error": str(error)}
            self.send_json(200, value); return
        if request.path == "/api/case":
            try:
                query = parse_qs(request.query, keep_blank_values=True)
                relative = query.get("case", [""])[0]
                parts = Path(relative).parts
                if (not parts or Path(relative).is_absolute() or
                        any(not re.fullmatch(r"[A-Za-z0-9_.-]{1,160}", part)
                            for part in parts)):
                    raise ValueError("invalid case path")
                case_dir = (self.root / relative).resolve()
                if not safe_child(self.root, case_dir):
                    raise ValueError("case is outside production root")
                catalog = config_catalog(self.root)
                actual = case_dir.exists() and any(
                    (case_dir / marker).is_file() for marker in
                    ("diagnostics.tsv", "checkpoint.json", "run.log",
                     "driver.log", "run_manifest.json", "run_manifest.failed.json"))
                if case_dir not in catalog and not actual:
                    raise ValueError("unknown production case")
                status_by_id, _ = scheduler_statuses(self.root)
                config_path = catalog.get(case_dir)
                config = read_config(config_path)
                run_id = (config or {}).get("case", {}).get("name", case_dir.name)
                value = snapshot_case(case_dir, self.root, config_path,
                                      bubble_process_table(),
                                      status_by_id.get(run_id))
                self.send_json(200, value)
            except (ValueError, TypeError) as error:
                self.send_json(400, {"ok": False, "error": str(error)})
            return
        if request.path == "/api/console":
            try:
                query = parse_qs(request.query, keep_blank_values=True)
                case_name = query.get("case", [""])[0]
                line_limit = int(query.get("lines", ["80"])[0])
                self.send_json(200, console_tail(self.root, case_name, line_limit))
            except (ValueError, TypeError) as error:
                self.send_json(400, {"ok": False, "error": str(error)})
            return
        if request.path in ("/", "/index.html"):
            payload = HTML.encode(); self.send_response(200); self.send_header("Content-Type", "text/html; charset=utf-8"); self.send_header("Cache-Control", "no-store, no-cache, must-revalidate"); self.end_headers()
            try:
                self.wfile.write(payload)
            except (BrokenPipeError, ConnectionResetError):
                pass
            return
        self.send_error(404)

    def do_POST(self):  # noqa: N802
        if not self.authorized():
            self.send_json(403, {"ok": False, "error": "control disabled or authentication failed"}); return
        if self.headers.get("Content-Type", "").split(";", 1)[0] != "application/json":
            self.send_json(415, {"ok": False, "error": "application/json required"}); return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > 16384:
                raise case_control.ControlError("invalid request size")
            body = json.loads(self.rfile.read(length))
            case_name = str(body.get("case_name", ""))
            with self.control_lock:
                if self.path == "/api/control/stop":
                    result = case_control.request_stop(case_name, self.root)
                elif self.path == "/api/control/stage":
                    result = case_control.stage_settings(case_name, body.get("changes", {}), self.root)
                elif self.path == "/api/control/restart":
                    result = case_control.restart_case(case_name, self.solver_binary, self.root)
                else:
                    self.send_json(404, {"ok": False, "error": "unknown control endpoint"}); return
            self.send_json(200, {"ok": True, "result": result})
        except (ValueError, json.JSONDecodeError, case_control.ControlError) as error:
            try:
                if 'case_name' in locals():
                    case_control.audit(case_control.resolve_case(case_name, self.root),
                        self.path.rsplit('/', 1)[-1], "failed", {"error": str(error)})
            except case_control.ControlError:
                pass
            self.send_json(400, {"ok": False, "error": str(error)})
    def log_message(self, *_):
        return


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--production-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--snapshot", action="store_true")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--enable-control", action="store_true")
    parser.add_argument("--control-token-file", type=Path)
    parser.add_argument("--solver-binary", type=Path,
                        default=Path("/workspace/zz/build-bubble-flow-release/bubble_solver"))
    args = parser.parse_args()
    root = args.production_root.resolve()
    if root != DEFAULT_ROOT:
        parser.error(f"production root is fixed to {DEFAULT_ROOT}; refusing to scan another path")
    if not root.is_dir():
        parser.error(f"production root is not a directory: {root}")
    if args.snapshot:
        print(json.dumps(snapshot(root), ensure_ascii=False, indent=2)); return
    if args.enable_control:
        if args.host not in ("127.0.0.1", "localhost", "::1"):
            parser.error("control mode may only bind to a loopback address")
        if not args.control_token_file or not args.control_token_file.is_file():
            parser.error("--control-token-file is required in control mode")
        if args.control_token_file.stat().st_mode & 0o077:
            parser.error("control token file permissions must be 0600 or stricter")
        token = args.control_token_file.read_text().strip()
        if len(token) < 32:
            parser.error("control token must contain at least 32 characters")
        if not args.solver_binary.resolve().is_file() or not os.access(args.solver_binary.resolve(), os.X_OK):
            parser.error("--solver-binary is missing or not executable")
        Handler.control_enabled = True
        Handler.control_token = token
        Handler.solver_binary = args.solver_binary.resolve()
    Handler.root = root
    ReusableThreadingHTTPServer((args.host, args.port), Handler).serve_forever()


if __name__ == "__main__":
    main()
