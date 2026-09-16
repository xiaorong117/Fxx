#!/usr/bin/env python3
"""Audit raw saturation-history completeness for all batch PNM variants."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def numeric_rows(path: Path, minimum_columns: int) -> int:
    count = 0
    for line in path.read_text(errors="strict").splitlines():
        fields = line.split()
        if len(fields) < minimum_columns:
            continue
        try:
            [float(value) for value in fields]
        except ValueError:
            continue
        count += 1
    return count


def declared_count(path: Path) -> int:
    return int(path.read_text(errors="strict").splitlines()[0].split()[0])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.input_root.resolve(); cases = []
    for pnm in sorted(root.rglob("PNM_pt.txt")):
        if pnm.parent.name not in ("PNM-0.6", "PNM-0.8"):
            continue
        directory = pnm.parent
        required = [directory / name for name in
                    ("default_throat.out", "default_throat_sw.out",
                     "default_pore.out", "default_pore_sw.out", "PNM_pb.txt")]
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            cases.append({"directory": str(directory), "status": "failed",
                          "missing": missing})
            continue
        throat_declared = declared_count(required[0])
        throat_history_rows = numeric_rows(required[1], 4)
        pore_declared = declared_count(required[2])
        pore_history_rows = numeric_rows(required[3], 6)
        throats = np.loadtxt(pnm)
        pores = np.loadtxt(directory / "PNM_pb.txt")
        zero_fraction = float(np.count_nonzero(throats[:, 9] == 0.0) / len(throats))
        failures = []
        if throat_history_rows != throat_declared:
            failures.append("default_throat_sw.out data-row count differs from default_throat.out")
        if pore_history_rows != pore_declared:
            failures.append("default_pore_sw.out data-row count differs from default_pore.out")
        if len(throats) > throat_declared:
            failures.append("PNM_pt internal row count exceeds declared raw throat count")
        if len(pores) != pore_declared:
            failures.append("PNM_pb row count differs from declared raw pore count")
        cases.append({
            "parent_case": directory.parent.name, "variant": directory.name,
            "directory": str(directory), "status": "failed" if failures else "passed",
            "failures": failures,
            "declared_throat_count": throat_declared,
            "throat_saturation_history_row_count": throat_history_rows,
            "missing_throat_saturation_history_rows": throat_declared - throat_history_rows,
            "PNM_pt_internal_throat_count": len(throats),
            "PNM_pt_zero_Sw_count": int(np.count_nonzero(throats[:, 9] == 0.0)),
            "PNM_pt_zero_Sw_fraction": zero_fraction,
            "declared_pore_count": pore_declared,
            "pore_saturation_history_row_count": pore_history_rows,
            "PNM_pb_count": len(pores),
            "sha256": {path.name: sha256(path) for path in [pnm, *required]},
        })
    failures = [case for case in cases if case["status"] == "failed"]
    report = {
        "format": "bubble_batch_saturation_source_audit_v1",
        "status": "failed" if failures else "passed",
        "case_count": len(cases), "passed_count": len(cases) - len(failures),
        "failed_count": len(failures),
        "conclusion": "incomplete upstream saturation histories cannot be repaired by zero fill",
        "failed_cases": failures, "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in
                      ("status", "case_count", "passed_count", "failed_count",
                       "failed_cases")}, ensure_ascii=False, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
