#!/usr/bin/env python3
"""Run a small, reproducible AMGX configuration set on the frozen Newton-0 system."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import pathlib
import subprocess


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=pathlib.Path, required=True)
    parser.add_argument("--build", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    project = args.project.resolve()
    build = args.build.resolve()
    output = args.output.resolve()
    result_dir = output / "frozen_results"
    if result_dir.exists():
        raise SystemExit(f"refusing to overwrite {result_dir}")
    result_dir.mkdir(parents=True)

    prefix = project / "output/verification/frozen_status3/newton0"
    matrix = pathlib.Path(str(prefix) + "_matrix_scaled.mtx")
    rhs = pathlib.Path(str(prefix) + "_rhs_scaled.mtx")
    executable = build / "solve_frozen_amgx"
    configs = [
        project / "configs/amgx_fgmres_block2_amg_nosolver.json",
        output / "configs/fgmres_amg_restart50.json",
        output / "configs/fgmres_amg_restart100.json",
        output / "configs/fgmres_dilu_restart50.json",
        output / "configs/fgmres_ilu1_restart50.json",
        output / "configs/pbicgstab_amg_dilu.json",
    ]
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    rows: list[dict[str, object]] = []
    for index, config in enumerate(configs, start=1):
        name = "reference_restart500" if index == 1 else config.stem
        atomic_json(output / "progress.json", {
            "stage": "frozen_configurations", "status": "running",
            "configuration": name, "completed": index - 1, "total": len(configs),
            "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        })
        wrapper = result_dir / f"{name}_wrapper.json"
        dump_prefix = result_dir / name
        log = result_dir / f"{name}.log"
        environment = os.environ.copy()
        environment["BUBBLE_AMGX_DUMP_PREFIX"] = str(dump_prefix)
        command = [str(executable), str(matrix), str(rhs), str(config), str(wrapper)]
        with log.open("w") as stream:
            completed = subprocess.run(command, env=environment, stdout=stream,
                                       stderr=subprocess.STDOUT, check=False)
        wrapped = json.loads(wrapper.read_text())
        internal_path = pathlib.Path(str(dump_prefix) + "_amgx_result.json")
        internal = json.loads(internal_path.read_text()) if internal_path.exists() else {}
        rows.append({
            "name": name,
            "config_path": str(config.resolve()),
            "config_sha256": sha256(config),
            "return_code": completed.returncode,
            "status": wrapped.get("status", "UNKNOWN"),
            "amgx_status_code": internal.get("status_code"),
            "amgx_status_name": internal.get("status_name"),
            "iterations": wrapped.get("iterations", internal.get("iterations")),
            "cpu_recomputed_relative_residual": wrapped.get("cpu_recomputed_relative_residual",
                                                               internal.get("cpu_recomputed_relative_residual")),
            "csr_seconds": wrapped.get("csr_seconds"),
            "upload_seconds": wrapped.get("upload_seconds"),
            "setup_seconds": wrapped.get("setup_seconds"),
            "solve_seconds": wrapped.get("solve_seconds"),
            "download_seconds": wrapped.get("download_seconds"),
            "error": wrapped.get("error"),
            "log": str(log),
            "residual_history": str(dump_prefix) + "_amgx_residual_history.tsv",
        })

    summary = {
        "format": "amgx_frozen_iteration_diagnosis_v1",
        "started_utc": started,
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "matrix": str(matrix), "matrix_sha256": sha256(matrix),
        "rhs": str(rhs), "rhs_sha256": sha256(rhs),
        "executable": str(executable), "executable_sha256": sha256(executable),
        "strict_acceptance": "AMGX_SOLVE_SUCCESS and CPU relative residual <= 1e-8",
        "results": rows,
    }
    atomic_json(output / "frozen_summary.json", summary)
    columns = ["name", "amgx_status_name", "iterations", "cpu_recomputed_relative_residual",
               "setup_seconds", "solve_seconds", "upload_seconds", "status", "return_code"]
    with (output / "frozen_summary.tsv").open("w") as stream:
        stream.write("\t".join(columns) + "\n")
        for row in rows:
            stream.write("\t".join("" if row.get(key) is None else str(row[key])
                                   for key in columns) + "\n")
    atomic_json(output / "progress.json", {
        "stage": "frozen_configurations", "status": "complete",
        "completed": len(configs), "total": len(configs),
        "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
