#!/usr/bin/env python3
"""Convert raw pnextract pore/throat tables to the solver's full 12/11-column graph."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_table(path: Path, header: str, rows: np.ndarray,
                integer_columns: int) -> None:
    formats = ["%d"] * integer_columns + ["%.17g"] * (
        rows.shape[1] - integer_columns)
    np.savetxt(path, rows, fmt=formats, delimiter="\t", header=header,
               comments="# ")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-pore", type=Path, required=True)
    parser.add_argument("--raw-throat", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--geometry-scale", type=float, default=1.0,
                        help="isotropic length scale applied to derived geometry")
    parser.add_argument("--normalize-zero-volume-saturation", action="store_true",
                        help="mark zero-volume controls Sw=1 in derived data")
    args = parser.parse_args()
    raw_pore = args.raw_pore.resolve()
    raw_throat = args.raw_throat.resolve()
    output = args.output_directory.resolve()
    scale = args.geometry_scale
    require(np.isfinite(scale) and scale > 0.0,
            "geometry scale must be finite and positive")
    require(not output.exists(), f"refusing to overwrite {output}")
    output.mkdir(parents=True)

    pores = np.loadtxt(raw_pore)
    throats = np.loadtxt(raw_throat)
    require(pores.ndim == 2 and pores.shape[1] == 9,
            "PNM_pb.txt must contain nine columns")
    require(throats.ndim == 2 and throats.shape[1] == 10,
            "PNM_pt.txt must contain ten columns")
    require(np.isfinite(pores).all() and np.isfinite(throats).all(),
            "raw PNM contains NaN or Inf")
    source_coordinate_min = pores[:, :3].min(axis=0).copy()
    source_coordinate_max = pores[:, :3].max(axis=0).copy()
    source_pore_volume_sum = float(pores[:, 6].sum())
    source_throat_volume_sum = float(throats[:, 5].sum())
    zero_pore_gas = np.flatnonzero((pores[:, 6] == 0.0) & (pores[:, 8] < 1.0))
    zero_throat_gas = np.flatnonzero((throats[:, 5] == 0.0) &
                                     (throats[:, 9] < 1.0))
    if (len(zero_pore_gas) or len(zero_throat_gas)) and not \
            args.normalize_zero_volume_saturation:
        raise RuntimeError(
            "zero-volume controls marked with gas require explicit "
            "--normalize-zero-volume-saturation")
    pores = pores.copy()
    throats = throats.copy()
    pores[zero_pore_gas, 8] = 1.0
    throats[zero_throat_gas, 9] = 1.0
    pores[:, 0:4] *= scale
    pores[:, 4] *= scale ** 2
    pores[:, 6] *= scale ** 3
    throats[:, 2] *= scale
    throats[:, 3] *= scale ** 2
    throats[:, 5] *= scale ** 3
    throats[:, 7:9] *= scale

    pore_count, throat_count = len(pores), len(throats)
    endpoints = throats[:, :2].astype(np.int64)
    require(np.array_equal(endpoints, throats[:, :2]),
            "throat endpoint IDs are non-integral")
    require(int(endpoints.min()) >= 1 and int(endpoints.max()) <= pore_count,
            "throat endpoint outside pore ID range")
    require(np.all(pores[:, 3:7] >= 0.0) and np.all(throats[:, 2:6] >= 0.0),
            "raw PNM contains negative geometry")
    require(np.all((pores[:, 8] >= 0.0) & (pores[:, 8] <= 1.0)) and
            np.all((throats[:, 9] >= 0.0) & (throats[:, 9] <= 1.0)),
            "raw saturation lies outside [0,1]")

    new_count = pore_count + throat_count
    throat_nodes = np.arange(pore_count + 1, new_count + 1, dtype=np.int64)
    new_pores = np.zeros((new_count, 12), dtype=np.float64)
    new_pores[:pore_count, 0] = np.arange(1, pore_count + 1)
    new_pores[:pore_count, 1] = 0
    new_pores[:pore_count, 2] = np.arange(1, pore_count + 1)
    new_pores[:pore_count, 3:12] = pores
    new_pores[pore_count:, 0] = throat_nodes
    new_pores[pore_count:, 1] = 1
    new_pores[pore_count:, 2] = np.arange(1, throat_count + 1)
    p1 = endpoints[:, 0] - 1
    p2 = endpoints[:, 1] - 1
    new_pores[pore_count:, 3:6] = 0.5 * (pores[p1, :3] + pores[p2, :3])
    new_pores[pore_count:, 6:10] = throats[:, 2:6]
    new_pores[pore_count:, 10] = 0
    new_pores[pore_count:, 11] = throats[:, 9]

    edges = np.zeros((2 * throat_count, 11), dtype=np.float64)
    edges[:, 0] = np.arange(1, 2 * throat_count + 1)
    edges[0::2, 1] = endpoints[:, 0]
    edges[0::2, 2] = throat_nodes
    edges[1::2, 1] = throat_nodes
    edges[1::2, 2] = endpoints[:, 1]
    edges[0::2, 3] = edges[1::2, 3] = np.arange(1, throat_count + 1)
    edges[0::2, 4] = 1
    edges[1::2, 4] = 2
    edges[0::2, 5] = throats[:, 7]
    edges[1::2, 5] = throats[:, 8]
    for column, source in ((6, 9), (7, 2), (8, 3), (9, 4), (10, 5)):
        edges[0::2, column] = throats[:, source]
        edges[1::2, column] = throats[:, source]

    pore_path = output / "PNM_new_pore_base.txt"
    throat_path = output / "PNM_new_throat_base.txt"
    connect_path = output / "PNM_new_connect_base.txt"
    write_table(pore_path,
        "new_pore_id\tsource_type(0=old_pore,1=old_throat)\tsource_id\t"
        "x\ty\tz\tradius\tarea\tshape_factor\tvolume\tboundary_flag\tSw",
        new_pores, 3)
    write_table(throat_path,
        "new_throat_id\tpore_id_1\tpore_id_2\tsource_old_throat_id\tside\t"
        "length\tSw\tradius\tarea\tshape_factor\tvolume", edges, 5)

    neighbors: list[list[int]] = [[] for _ in range(new_count)]
    for first, second in edges[:, 1:3].astype(np.int64):
        neighbors[first - 1].append(int(second))
        neighbors[second - 1].append(int(first))
    maximum_coordination = max(map(len, neighbors))
    with connect_path.open("w", encoding="utf-8") as stream:
        columns = ["pore_id", "coordination_number"] + [
            f"connection_pair_{index:03d}"
            for index in range(1, maximum_coordination + 1)]
        stream.write("# " + "\t".join(columns) + "\n")
        for node_id, connected in enumerate(neighbors, 1):
            connected.sort()
            row = [str(node_id), str(len(connected))]
            row.extend(f"{node_id}-{neighbor}" for neighbor in connected)
            row.extend([""] * (maximum_coordination - len(connected)))
            stream.write("\t".join(row) + "\n")

    audit = {
        "format": "bubble_batch_full_geometry_conversion_v1",
        "raw_pore": {"path": str(raw_pore), "sha256": sha256(raw_pore)},
        "raw_throat": {"path": str(raw_throat), "sha256": sha256(raw_throat)},
        "raw_pore_count": pore_count,
        "raw_throat_count": throat_count,
        "converted_control_count": new_count,
        "converted_half_edge_count": len(edges),
        "zero_half_length_count": int(np.count_nonzero(edges[:, 5] == 0.0)),
        "massless_control_count": int(np.count_nonzero(new_pores[:, 9] == 0.0)),
        "boundary_flag_counts": {
            str(flag): int(np.count_nonzero(pores[:, 7].astype(int) == flag))
            for flag in sorted(set(pores[:, 7].astype(int)))},
        "geometry_normalization": {
            "isotropic_length_scale": scale,
            "area_scale": scale ** 2,
            "volume_scale": scale ** 3,
            "source_coordinate_min_um": source_coordinate_min.tolist(),
            "source_coordinate_max_um": source_coordinate_max.tolist(),
            "derived_coordinate_min_um": pores[:, :3].min(axis=0).tolist(),
            "derived_coordinate_max_um": pores[:, :3].max(axis=0).tolist(),
            "source_pore_volume_sum_um3": source_pore_volume_sum,
            "source_throat_volume_sum_um3": source_throat_volume_sum,
            "derived_pore_volume_sum_um3": float(pores[:, 6].sum()),
            "derived_throat_volume_sum_um3": float(throats[:, 5].sum()),
        },
        "zero_measure_saturation_normalization": {
            "enabled": args.normalize_zero_volume_saturation,
            "policy": "set derived Sw=1 only where source control volume is exactly zero",
            "source_pore_ids": (zero_pore_gas + 1).tolist(),
            "source_throat_ids": (zero_throat_gas + 1).tolist(),
            "global_liquid_volume_change_um3": 0.0,
            "global_gas_volume_change_um3": 0.0,
        },
        "coordinate_min_um": pores[:, :3].min(axis=0).tolist(),
        "coordinate_max_um": pores[:, :3].max(axis=0).tolist(),
        "pore_saturation_range": [float(pores[:, 8].min()),
                                   float(pores[:, 8].max())],
        "throat_saturation_range": [float(throats[:, 9].min()),
                                     float(throats[:, 9].max())],
        "source_throat_volume_copied_to_both_half_edge_metadata": True,
        "source_throat_volume_counted_once_as_throat_control_volume": True,
        "outputs": {
            str(pore_path): sha256(pore_path),
            str(throat_path): sha256(throat_path),
            str(connect_path): sha256(connect_path),
        },
    }
    (output / "conversion_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pore_count": new_count, "edge_count": len(edges),
                      "maximum_coordination": maximum_coordination}))


if __name__ == "__main__":
    main()
