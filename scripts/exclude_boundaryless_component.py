#!/usr/bin/env python3
"""Create a derived PNM with one audited boundaryless component frozen out.

The raw/production network is never modified.  The excluded component's initial
inventory is written to an audit and can be added back in post-processing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import deque
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_table(path: Path, header: str, values: np.ndarray,
                integer_columns: int) -> None:
    formats = ["%d"] * integer_columns + ["%.17g"] * (
        values.shape[1] - integer_columns)
    np.savetxt(path, values, fmt=formats, delimiter="\t", header=header,
               comments="# ")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pore", type=Path, required=True)
    parser.add_argument("--throat", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--component-node-id", type=int, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    pore_path, throat_path = args.pore.resolve(), args.throat.resolve()
    config_path = args.config.resolve()
    output = args.output_directory.resolve()
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    original_hashes = {"pore": sha256(pore_path), "throat": sha256(throat_path),
                       "config": sha256(config_path)}
    pores = np.loadtxt(pore_path)
    edges = np.loadtxt(throat_path)
    if pores.ndim != 2 or pores.shape[1] != 12 or edges.ndim != 2 or edges.shape[1] != 11:
        raise RuntimeError("expected 12-column pore and 11-column throat tables")
    if args.component_node_id < 1 or args.component_node_id > len(pores):
        raise RuntimeError("component node ID lies outside network")
    adjacency: list[list[int]] = [[] for _ in range(len(pores))]
    for row in edges:
        a, b = int(row[1]) - 1, int(row[2]) - 1
        adjacency[a].append(b); adjacency[b].append(a)
    start = args.component_node_id - 1
    excluded = {start}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        for neighbor in adjacency[node]:
            if neighbor not in excluded:
                excluded.add(neighbor); queue.append(neighbor)
    flags = {int(pores[node, 10]) for node in excluded}
    if flags != {0}:
        raise RuntimeError(f"refusing to exclude a boundary-connected component: flags={flags}")
    excluded_edges = [k for k, row in enumerate(edges)
                      if int(row[1]) - 1 in excluded or int(row[2]) - 1 in excluded]
    if any((int(edges[k, 1]) - 1 in excluded) !=
           (int(edges[k, 2]) - 1 in excluded) for k in excluded_edges):
        raise RuntimeError("component traversal left a cut edge")

    keep_nodes = [node for node in range(len(pores)) if node not in excluded]
    old_to_new = {old: new for new, old in enumerate(keep_nodes)}
    kept_pores = pores[keep_nodes].copy()
    kept_pores[:, 0] = np.arange(1, len(kept_pores) + 1)
    # Virtual source_id points to a node ID, unlike source IDs for real nodes.
    for new, old in enumerate(keep_nodes):
        if int(kept_pores[new, 1]) == 2:
            old_source = int(pores[old, 2]) - 1
            if old_source not in old_to_new:
                raise RuntimeError("kept virtual node references excluded real node")
            kept_pores[new, 2] = old_to_new[old_source] + 1
    keep_edge_rows = [row for row in edges
                      if int(row[1]) - 1 not in excluded and
                      int(row[2]) - 1 not in excluded]
    kept_edges = np.asarray(keep_edge_rows, dtype=float).copy()
    kept_edges[:, 0] = np.arange(1, len(kept_edges) + 1)
    for row in kept_edges:
        row[1] = old_to_new[int(row[1]) - 1] + 1
        row[2] = old_to_new[int(row[2]) - 1] + 1

    pore_out = output / "PNM_new_pore_virtual_dynamic.txt"
    throat_out = output / "PNM_new_throat_virtual_dynamic.txt"
    connect_out = output / "PNM_new_connect_virtual_dynamic.txt"
    write_table(pore_out,
        "new_pore_id\tsource_type\tsource_id\tx\ty\tz\tradius\tarea\tshape_factor\tvolume\tboundary_flag\tSw",
        kept_pores, 3)
    write_table(throat_out,
        "new_throat_id\tpore_id_1\tpore_id_2\tsource_old_throat_id\tside\tlength\tSw\tradius\tarea\tshape_factor\tvolume",
        kept_edges, 5)
    new_neighbors: list[list[int]] = [[] for _ in range(len(kept_pores))]
    for row in kept_edges:
        a, b = int(row[1]), int(row[2])
        new_neighbors[a - 1].append(b); new_neighbors[b - 1].append(a)
    maximum = max(map(len, new_neighbors), default=0)
    with connect_out.open("w") as stream:
        stream.write("# pore_id\tcoordination_number" + "".join(
            f"\tconnection_pair_{i:03d}" for i in range(1, maximum + 1)) + "\n")
        for node_id, neighbors in enumerate(new_neighbors, 1):
            neighbors.sort()
            values = [str(node_id), str(len(neighbors))]
            values.extend(f"{node_id}-{neighbor}" for neighbor in neighbors)
            values.extend([""] * (maximum - len(neighbors)))
            stream.write("\t".join(values) + "\n")
    with (output / "old_to_new_node_id.tsv").open("w") as stream:
        stream.write("old_node_id\tnew_node_id\tsource_type\tsource_id\n")
        for old in keep_nodes:
            stream.write(f"{old+1}\t{old_to_new[old]+1}\t{int(pores[old,1])}\t{int(pores[old,2])}\n")
    with (output / "frozen_boundaryless_component.tsv").open("w") as stream:
        stream.write("old_node_id\tsource_type\tsource_id\tx_um\ty_um\tz_um\tvolume_um3\tinput_Sw\n")
        for node in sorted(excluded):
            row = pores[node]
            stream.write(f"{node+1}\t{int(row[1])}\t{int(row[2])}\t{row[3]:.17g}\t{row[4]:.17g}\t{row[5]:.17g}\t{row[9]:.17g}\t{row[11]:.17g}\n")

    config = json.loads(config_path.read_text())
    length_to_m = float(config["geometry"]["input_length_to_m"])
    volume_to_m3 = float(config["geometry"]["input_volume_to_m3"])
    physics = config["physics"]; active = config["active_set"]
    background = float(config["pressure_reference"]["background_abs_Pa"])
    pl_in = float(config["boundary"]["Pl_in_rel_Pa"])
    pl_out = float(config["boundary"]["Pl_out_rel_Pa"])
    initial_c = float(config["initial"]["C_initial_mol_m3"])
    xmin, xmax = float(pores[:, 3].min()), float(pores[:, 3].max())
    free_moles = dissolved_moles = liquid_volume = gas_volume = 0.0
    active_bubbles = 0
    initialized = []
    for node in sorted(excluded):
        row = pores[node]
        volume = float(row[9]) * volume_to_m3
        sw = float(row[11])
        if sw <= float(active["Sw_min"]):
            sw = min(100.0 * float(active["Sw_min"]),
                     0.5 * (float(active["Sw_min"]) + 1.0 - float(active["Sg_off"])))
        if volume == 0.0:
            sw = 1.0
        liquid_volume += volume * sw
        gas_volume += volume * (1.0 - sw)
        dissolved = volume * sw * initial_c
        dissolved_moles += dissolved
        ng = 0.0
        pg_rel = None
        if volume > 0.0 and 1.0 - sw > float(active["Sg_off"]):
            xi = (float(row[3]) - xmin) / (xmax - xmin)
            pl = pl_in + xi * (pl_out - pl_in)
            radius = float(row[6]) * length_to_m
            pc = (2.0 * float(physics["sigma_N_m"]) *
                  math.cos(float(physics["contact_angle_rad"])) /
                  (radius * -math.expm1(-float(physics["qin_b"]) * sw)))
            pg_rel = pl + pc
            pg_abs = background + pg_rel
            ng = (pg_abs * volume * (1.0 - sw) /
                  (float(physics["Z"]) * 8.31446261815324 * float(physics["T_K"])))
            if ng > float(active["ng_off_mol"]):
                free_moles += ng; active_bubbles += 1
            else:
                dissolved_moles += ng; ng = 0.0
        initialized.append({"old_node_id": node + 1, "initialized_Sw": sw,
                            "initialized_Pg_relative_Pa": pg_rel,
                            "initialized_free_gas_mol": ng,
                            "initialized_dissolved_gas_mol": dissolved})

    excluded_volume = float(pores[list(excluded), 9].sum()) * volume_to_m3
    original_volume = float(pores[pores[:, 1] != 2, 9].sum()) * volume_to_m3
    audit = {
        "format": "bubble_frozen_boundaryless_component_audit_v1",
        "status": "passed",
        "strategy": "exclude exactly one fully boundaryless component from dynamics and freeze its initial inventory for post-processing",
        "scope": "case-specific recovery only; no solver or other-case policy change",
        "requested_component_node_id": args.component_node_id,
        "excluded_old_node_ids": [node + 1 for node in sorted(excluded)],
        "excluded_component_node_count": len(excluded),
        "excluded_component_edge_count": len(excluded_edges),
        "excluded_boundary_flag_set": sorted(flags),
        "excluded_component_volume_m3": excluded_volume,
        "original_real_control_volume_m3": original_volume,
        "excluded_volume_fraction": excluded_volume / original_volume,
        "frozen_initial_liquid_volume_m3": liquid_volume,
        "frozen_initial_gas_volume_m3": gas_volume,
        "frozen_initial_dissolved_gas_mol": dissolved_moles,
        "frozen_initial_free_gas_mol": free_moles,
        "frozen_initial_total_component_mol": dissolved_moles + free_moles,
        "frozen_initial_active_bubble_count": active_bubbles,
        "initialized_frozen_nodes": initialized,
        "original_node_count": len(pores), "dynamic_node_count": len(kept_pores),
        "original_edge_count": len(edges), "dynamic_edge_count": len(kept_edges),
        "original_inputs": {"pore": str(pore_path), "throat": str(throat_path),
                            "config": str(config_path), "sha256": original_hashes},
        "derived_inputs": {
            "pore": str(pore_out), "pore_sha256": sha256(pore_out),
            "throat": str(throat_out), "throat_sha256": sha256(throat_out),
            "connect": str(connect_out), "connect_sha256": sha256(connect_out)},
        "raw_inputs_unchanged": original_hashes == {
            "pore": sha256(pore_path), "throat": sha256(throat_path),
            "config": sha256(config_path)},
        "scientific_reporting": "primary main-connected-domain result is dynamic; add frozen inventory only to explicitly labeled whole-network totals",
    }
    (output / "boundaryless_component_freeze_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
    virtual = kept_pores[:, 1].astype(int) == 2
    kept_flags = kept_pores[:, 10].astype(int)
    boundary_edges = sum(
        bool(virtual[int(row[1]) - 1]) != bool(virtual[int(row[2]) - 1])
        for row in kept_edges)
    virtual_audit = {
        "format": "bubble_filtered_virtual_boundary_audit_v1",
        "status": "passed",
        "parent_network_pore_file": str(pore_path),
        "parent_network_throat_file": str(throat_path),
        "filter_audit": str(output / "boundaryless_component_freeze_audit.json"),
        "derived_node_count": len(kept_pores),
        "derived_edge_count": len(kept_edges),
        "virtual_node_count": int(np.count_nonzero(virtual)),
        "inlet_virtual_count": int(np.count_nonzero(virtual & (kept_flags == 1))),
        "outlet_virtual_count": int(np.count_nonzero(virtual & (kept_flags == 2))),
        "boundary_edge_count": boundary_edges,
        "virtual_node_coordination_distribution": {
            str(value): int(sum(len(new_neighbors[i]) == value
                                for i in np.flatnonzero(virtual)))
            for value in sorted({len(new_neighbors[i]) for i in np.flatnonzero(virtual)})},
        "all_virtual_nodes_coordination_one": all(
            len(new_neighbors[i]) == 1 for i in np.flatnonzero(virtual)),
        "excluded_component_was_boundaryless": flags == {0},
        "derived_pore_sha256": sha256(pore_out),
        "derived_throat_sha256": sha256(throat_out),
        "derived_connect_sha256": sha256(connect_out),
    }
    if not virtual_audit["all_virtual_nodes_coordination_one"]:
        raise RuntimeError("filtered virtual node lost one-to-one topology")
    (output / "filtered_virtual_boundary_audit.json").write_text(
        json.dumps(virtual_audit, indent=2) + "\n")
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
