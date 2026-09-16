#!/usr/bin/env python3
"""Audit degenerate PNM geometry without inventing replacement values."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def rows(path: Path):
    with path.open(encoding="utf-8-sig") as stream:
        for line_number, line in enumerate(stream, 1):
            text = line.strip()
            if text and not text.startswith("#"):
                yield line_number, [float(item) for item in text.split()]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", nargs="?", default=".", type=Path)
    parser.add_argument("--report", type=Path, help="write every defect as TSV")
    parser.add_argument("--pore-file", default="PNM_new_pore.txt")
    parser.add_argument("--throat-file", default="PNM_new_throat.txt")
    args = parser.parse_args()
    root = args.root
    failures: list[str] = []
    defects: list[dict[str, object]] = []
    categories: dict[str, list[int]] = {}
    pore_flags: list[int] = []
    pore_sw: list[float] = []
    pore_volume: list[float] = []
    pore_source_type: list[int] = []
    edge_endpoints: list[tuple[int, int]] = []
    edge_ids: list[int] = []
    edge_lengths: list[float] = []
    zero_distance: list[int] = []
    massless: list[int] = []

    def defect(kind: str, entity: str, entity_id: int, source_id: int,
               field: str, value: float, endpoint_a: object = "",
               endpoint_b: object = "", side: object = "") -> None:
        categories.setdefault(kind, []).append(entity_id)
        defects.append({
            "kind": kind, "entity": entity, "entity_id": entity_id,
            "source_id": source_id, "field": field, "value": value,
            "endpoint_a": endpoint_a, "endpoint_b": endpoint_b, "side": side,
        })

    pore_fields = {6: "radius", 7: "area", 8: "shape_factor"}
    edge_fields = {7: "radius", 8: "area", 9: "shape_factor"}
    for _, row in rows(root / args.pore_file):
        if len(row) != 12:
            failures.append("pore row does not have 12 columns")
            break
        if not all(math.isfinite(value) for value in row):
            failures.append(f"pore {int(row[0])} contains non-finite data")
        pore_flags.append(int(row[10]))
        pore_source_type.append(int(row[1]))
        pore_sw.append(row[11])
        pore_volume.append(row[9])
        for column, field in pore_fields.items():
            if row[column] <= 0:
                defect(f"non_positive_pore_{field}", "pore", int(row[0]),
                       int(row[2]), field, row[column])
        if row[9] < 0:
            defect("negative_pore_volume", "pore", int(row[0]), int(row[2]),
                   "volume", row[9])
        elif row[9] == 0 and row[11] < 1:
            defect("zero_volume_with_gas", "pore", int(row[0]), int(row[2]),
                   "volume", row[9])
        elif row[9] == 0 and int(row[1]) != 2:
            massless.append(int(row[0]))
    for _, row in rows(root / args.throat_file):
        if len(row) != 11:
            failures.append("edge row does not have 11 columns")
            break
        if not all(math.isfinite(value) for value in row):
            failures.append(f"edge {int(row[0])} contains non-finite data")
        edge_endpoints.append((int(row[1]) - 1, int(row[2]) - 1))
        edge_ids.append(int(row[0]))
        edge_lengths.append(row[5])
        for column, field in edge_fields.items():
            if row[column] <= 0:
                defect(f"non_positive_edge_{field}", "edge", int(row[0]),
                       int(row[3]), field, row[column], int(row[1]),
                       int(row[2]), int(row[4]))
        if row[5] < 0:
            defect("negative_edge_half_length", "edge", int(row[0]),
                   int(row[3]), "half_length", row[5], int(row[1]),
                   int(row[2]), int(row[4]))
        elif row[5] == 0:
            zero_distance.append(int(row[0]))
        if row[10] < 0:
            defect("negative_edge_source_volume_metadata", "edge", int(row[0]),
                   int(row[3]), "source_volume_metadata", row[10], int(row[1]),
                   int(row[2]), int(row[4]))
        elif row[10] == 0 and row[6] < 1:
            defect("zero_source_volume_with_gas", "edge", int(row[0]),
                   int(row[3]), "source_volume_metadata", row[10], int(row[1]),
                   int(row[2]), int(row[4]))
    adjacency = [[] for _ in pore_flags]
    for a, b in edge_endpoints:
        if 0 <= a < len(adjacency) and 0 <= b < len(adjacency):
            adjacency[a].append(b)
            adjacency[b].append(a)
    visited = [False] * len(adjacency)
    component_count = 0
    boundaryless_count = 0
    for root_node in range(len(adjacency)):
        if visited[root_node]:
            continue
        component_count += 1
        stack = [root_node]
        visited[root_node] = True
        has_boundary = False
        while stack:
            node = stack.pop()
            has_boundary = has_boundary or pore_flags[node] != 0
            for neighbor in adjacency[node]:
                if not visited[neighbor]:
                    visited[neighbor] = True
                    stack.append(neighbor)
        if not has_boundary:
            boundaryless_count += 1
    print(
        "INFO  topology: "
        f"components={component_count}, boundaryless_components={boundaryless_count}, "
        f"isolated_pores={sum(not neighbors for neighbors in adjacency)}"
    )
    parent = list(range(len(pore_flags)))

    def find(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a: int, b: int) -> None:
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a

    for (a, b), length in zip(edge_endpoints, edge_lengths):
        if length == 0:
            union(a, b)
    groups: dict[int, list[int]] = {}
    for node in range(len(parent)):
        groups.setdefault(find(node), []).append(node)
    cluster_id = {root: index + 1 for index, root in enumerate(groups)}
    external: dict[int, list[int]] = {root: [] for root in groups}
    zero_edges_by_cluster: dict[int, list[int]] = {root: [] for root in groups}
    for edge_id, (a, b), length in zip(edge_ids, edge_endpoints, edge_lengths):
        ca, cb = find(a), find(b)
        if length == 0:
            zero_edges_by_cluster[ca].append(edge_id)
        elif ca != cb:
            external[ca].append(edge_id)
            external[cb].append(edge_id)
    print(
        "PASS_WITH_CLASSIFICATION  degenerate geometry: "
        f"raw_nodes={len(pore_flags)}, effective_clusters={len(groups)}, "
        f"raw_edges={len(edge_ids)}, production_positive_edges="
        f"{sum(length > 0 for length in edge_lengths)}, "
        f"zero_distance_constraints={len(zero_distance)}, "
        f"massless_junctions={len(massless)}"
    )
    print(f"INFO  massless junction node IDs: {massless}")
    print(f"INFO  zero-distance edge IDs: {zero_distance[:20]}"
          f"{' ...' if len(zero_distance) > 20 else ''}")
    for root, members in groups.items():
        if len(members) > 1:
            print(
                f"INFO  zero-distance cluster {cluster_id[root]}: "
                f"members={[node + 1 for node in members]}, "
                f"zero_edges={zero_edges_by_cluster[root]}, "
                f"external_positive_edges={external[root]}"
            )
    internal = [flag == 0 for flag in pore_flags]
    raw_storage = math.fsum(v for v, keep in zip(pore_volume, internal) if keep)
    clustered_storage = math.fsum(
        math.fsum(pore_volume[node] for node in members if internal[node])
        for members in groups.values()
    )
    raw_liquid = math.fsum(v * sw for v, sw, keep in zip(pore_volume, pore_sw, internal)
                           if keep)
    clustered_liquid = math.fsum(
        math.fsum(pore_volume[node] * pore_sw[node] for node in members if internal[node])
        for members in groups.values()
    )
    raw_gas = math.fsum(v * (1 - sw) for v, sw, keep in zip(pore_volume, pore_sw, internal)
                        if keep)
    clustered_gas = math.fsum(
        math.fsum(pore_volume[node] * (1 - pore_sw[node]) for node in members
                  if internal[node])
        for members in groups.values()
    )
    tolerance = 32 * math.ulp(max(abs(raw_storage), 1.0))
    if (abs(raw_storage - clustered_storage) > tolerance or
            abs(raw_liquid - clustered_liquid) > tolerance or
            abs(raw_gas - clustered_gas) > tolerance):
        failures.append("cluster aggregation changed storage/liquid/gas volume")
    print(
        "INFO  aggregation invariants: "
        f"storage={raw_storage:.17e}/{clustered_storage:.17e}, "
        f"liquid={raw_liquid:.17e}/{clustered_liquid:.17e}, "
        f"gas={raw_gas:.17e}/{clustered_gas:.17e}"
    )
    for kind, ids in categories.items():
        failures.append(f"{kind}: count={len(ids)}, first_ids={ids[:20]}")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open("w", encoding="utf-8", newline="") as stream:
            columns = ["kind", "entity", "entity_id", "source_id", "field",
                       "value", "endpoint_a", "endpoint_b", "side"]
            writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(defects)
        print(f"INFO  wrote {len(defects)} defect records to {args.report}")
    if failures:
        for failure in failures:
            print("FAIL ", failure)
        print(f"STRICT GEOMETRY AUDIT FAILED: {len(failures)} category/categories")
        return 1
    print("DEGENERATE GEOMETRY AUDIT PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
