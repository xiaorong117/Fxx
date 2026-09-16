#!/usr/bin/env python3
"""Compare unchanged and filtered networks by immutable source identity."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def checkpoint(path: Path) -> dict:
    value = json.loads(path.read_text())
    if value.get("format") != "bubble_checkpoint_v6":
        raise RuntimeError("unexpected checkpoint format")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--filtered", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    original = checkpoint(args.original / "checkpoint.json")
    filtered = checkpoint(args.filtered / "checkpoint.json")
    if not math.isclose(original["time_s"], filtered["time_s"], rel_tol=0.0,
                        abs_tol=1e-15):
        raise RuntimeError("comparison times differ")
    mapping = []
    with args.mapping.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            mapping.append((int(row["old_node_id"]) - 1,
                            int(row["new_node_id"]) - 1))
    maxima = {"Pl_Pa": 0.0, "Pg_Pa": 0.0, "Sw": 0.0, "C_mol_m3": 0.0}
    maximum_nodes = {key: 0 for key in maxima}
    for old, new in mapping:
        a, b = original["state"][old], filtered["state"][new]
        for key, column in (("Pl_Pa", 0), ("Sw", 2), ("C_mol_m3", 3)):
            difference = abs(float(a[column]) - float(b[column]))
            if difference > maxima[key]: maxima[key], maximum_nodes[key] = difference, old + 1
        if bool(a[4]) != bool(b[4]) or bool(a[5]) != bool(b[5]):
            raise RuntimeError(f"active-set mismatch at old node {old+1}")
        if a[4]:
            difference = abs(float(a[1]) - float(b[1]))
            if difference > maxima["Pg_Pa"]:
                maxima["Pg_Pa"], maximum_nodes["Pg_Pa"] = difference, old + 1
    def last(path: Path) -> dict[str, str]:
        with (path / "diagnostics.tsv").open(newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))[-1]
    a, b = last(args.original), last(args.filtered)
    globals_ = {}
    for field in ("main_connected_component_Sg", "free_mol_main_flow_component",
                  "gas_volume_main_flow_component"):
        globals_[field] = abs(float(a[field]) - float(b[field]))
    tolerances = {"Pl_Pa": 1e-7, "Pg_Pa": 1e-7, "Sw": 1e-10,
                  "C_mol_m3": 1e-10}
    failures = [f"{key}={value}>{tolerances[key]}" for key, value in maxima.items()
                if value > tolerances[key]]
    failures += [f"{key}={value}>1e-10" for key, value in globals_.items()
                 if value > 1e-10]
    report = {
        "format": "bubble_boundaryless_recovery_equivalence_v1",
        "status": "passed" if not failures else "failed",
        "time_s": original["time_s"], "mapped_node_count": len(mapping),
        "maximum_absolute_state_differences": maxima,
        "maximum_difference_old_node_ids": maximum_nodes,
        "main_domain_absolute_differences": globals_,
        "tolerances": tolerances, "failures": failures,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
