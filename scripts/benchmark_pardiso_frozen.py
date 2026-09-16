#!/usr/bin/env python3
"""Benchmark oneMKL PARDISO phases on a frozen nonsymmetric linear system."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import pathlib
import resource
import subprocess
import sys
import time


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: pathlib.Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def worker(args: argparse.Namespace) -> int:
    # Keep this diagnostic below the priority of the ongoing production runs.
    os.nice(10)
    import numpy as np
    from scipy.io import mmread
    from scipy.sparse import csr_matrix
    from pypardiso import PyPardisoSolver

    started = dt.datetime.now(dt.timezone.utc).isoformat()
    load_started = time.perf_counter()
    matrix = csr_matrix(mmread(args.matrix), dtype=np.float64)
    rhs_loaded = mmread(args.rhs)
    rhs = np.asarray(rhs_loaded, dtype=np.float64).reshape(-1)
    matrix.sum_duplicates()
    matrix.sort_indices()
    if matrix.indices.dtype != np.int32:
        matrix.indices = matrix.indices.astype(np.int32)
    if matrix.indptr.dtype != np.int32:
        matrix.indptr = matrix.indptr.astype(np.int32)
    load_seconds = time.perf_counter() - load_started

    solver = PyPardisoSolver(mtype=11)
    zeros = np.zeros_like(rhs)
    result: dict[str, object] = {
        "format": "pardiso_frozen_worker_v1",
        "threads": args.threads,
        "started_utc": started,
        "matrix_rows": int(matrix.shape[0]),
        "matrix_nonzeros": int(matrix.nnz),
        "load_seconds": load_seconds,
        "mtype": 11,
    }
    try:
        phase_started = time.perf_counter()
        solver.set_phase(11)
        solver._call_pardiso(matrix, zeros)
        result["symbolic_analysis_seconds"] = time.perf_counter() - phase_started

        phase_started = time.perf_counter()
        solver.set_phase(22)
        solver._call_pardiso(matrix, zeros)
        result["numeric_factorization_seconds"] = time.perf_counter() - phase_started

        phase_started = time.perf_counter()
        solver.set_phase(33)
        solution = solver._call_pardiso(matrix, rhs)
        result["solve_seconds"] = time.perf_counter() - phase_started
        residual = matrix @ solution - rhs
        result["cpu_recomputed_relative_residual"] = float(
            np.linalg.norm(residual) / max(np.linalg.norm(rhs), 1e-300)
        )
        iparm = solver.get_iparms()
        result.update({
            "pardiso_iterative_refinement_steps": int(iparm[7]),
            "pardiso_symbolic_peak_memory_kib": int(iparm[15]),
            "pardiso_symbolic_permanent_memory_kib": int(iparm[16]),
            "pardiso_numeric_memory_kib": int(iparm[17]),
            "pardiso_factor_nonzeros": int(iparm[18]),
            "pardiso_factor_mflops": int(iparm[19]),
            "status": "passed",
        })
    except Exception as error:  # preserve every failed worker as evidence
        result.update({"status": "failed", "error": repr(error)})
    finally:
        try:
            solver.set_phase(-1)
            solver._call_pardiso(matrix, zeros)
        except Exception as cleanup_error:
            result["cleanup_error"] = repr(cleanup_error)
    result["maximum_resident_set_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result["completed_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    atomic_json(args.result, result)
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["status"] == "passed" else 2


def parent(args: argparse.Namespace) -> int:
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    matrix = args.matrix.resolve()
    rhs = args.rhs.resolve()
    threads = [int(value) for value in args.thread_list.split(",")]
    rows: list[dict[str, object]] = []
    for index, count in enumerate(threads, start=1):
        atomic_json(output / "progress.json", {
            "stage": "pardiso_frozen_thread_scaling",
            "status": "running",
            "threads": count,
            "completed": index - 1,
            "total": len(threads),
            "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        })
        result_path = output / f"threads_{count}.json"
        log_path = output / f"threads_{count}.log"
        environment = os.environ.copy()
        environment.update({
            "MKL_NUM_THREADS": str(count),
            "OMP_NUM_THREADS": str(count),
            "MKL_DYNAMIC": "FALSE",
            "OMP_DYNAMIC": "FALSE",
            "OMP_PROC_BIND": "close",
            "OMP_PLACES": "cores",
            "PYPARDISO_MKL_RT": "/opt/conda/envs/Mesh/lib/libmkl_rt.so.2",
            "LD_LIBRARY_PATH": "/opt/conda/envs/Mesh/lib:" + environment.get("LD_LIBRARY_PATH", ""),
        })
        command = [sys.executable, str(pathlib.Path(__file__).resolve()),
                   "--worker", "--matrix", str(matrix), "--rhs", str(rhs),
                   "--threads", str(count), "--result", str(result_path)]
        with log_path.open("w") as stream:
            try:
                completed = subprocess.run(command, env=environment, stdout=stream,
                                           stderr=subprocess.STDOUT, timeout=args.timeout_seconds,
                                           check=False)
                return_code = completed.returncode
            except subprocess.TimeoutExpired:
                return_code = 124
                atomic_json(result_path, {"status": "timeout", "threads": count,
                                          "timeout_seconds": args.timeout_seconds})
        row = json.loads(result_path.read_text())
        row["return_code"] = return_code
        row["log"] = str(log_path)
        rows.append(row)
        if row.get("status") != "passed":
            break

    summary = {
        "format": "pardiso_frozen_thread_scaling_v1",
        "matrix": str(matrix), "matrix_sha256": sha256(matrix),
        "rhs": str(rhs), "rhs_sha256": sha256(rhs),
        "python": sys.executable,
        "mkl_runtime": "/opt/conda/envs/Mesh/lib/libmkl_rt.so.2",
        "mkl_runtime_sha256": sha256(pathlib.Path("/opt/conda/envs/Mesh/lib/libmkl_rt.so.2")),
        "results": rows,
    }
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "progress.json", {
        "stage": "pardiso_frozen_thread_scaling",
        "status": "passed" if len(rows) == len(threads) and all(r.get("status") == "passed" for r in rows) else "failed",
        "completed": len(rows), "total": len(threads),
        "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    })
    return 0 if len(rows) == len(threads) and all(r.get("status") == "passed" for r in rows) else 2


def parse() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--matrix", type=pathlib.Path, required=True)
    parser.add_argument("--rhs", type=pathlib.Path, required=True)
    parser.add_argument("--result", type=pathlib.Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--thread-list", default="8,16,32,64")
    parser.add_argument("--timeout-seconds", type=int, default=900)
    args = parser.parse_args()
    if args.worker and args.result is None:
        parser.error("--worker requires --result")
    if not args.worker and args.output is None:
        parser.error("parent mode requires --output")
    return args


if __name__ == "__main__":
    arguments = parse()
    raise SystemExit(worker(arguments) if arguments.worker else parent(arguments))
