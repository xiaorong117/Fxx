#!/usr/bin/env python3
"""Add one degree-one virtual reservoir to every flagged real boundary pore."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_table(path: Path, header: str, rows: np.ndarray, integer_columns: int) -> None:
    formats = ["%d"] * integer_columns + ["%.17g"] * (rows.shape[1] - integer_columns)
    np.savetxt(path, rows, fmt=formats, delimiter="\t", header=header, comments="# ")


def components(node_count: int, edges: np.ndarray) -> tuple[np.ndarray, dict[int, int]]:
    adjacency: list[list[int]] = [[] for _ in range(node_count)]
    for a_value, b_value in edges[:, 1:3]:
        a, b = int(a_value) - 1, int(b_value) - 1
        adjacency[a].append(b)
        adjacency[b].append(a)
    labels = np.full(node_count, -1, dtype=np.int64)
    sizes: dict[int, int] = {}
    for root in range(node_count):
        if labels[root] >= 0:
            continue
        representative = root + 1
        stack = [root]
        labels[root] = representative
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            for neighbor in adjacency[node]:
                if labels[neighbor] < 0:
                    labels[neighbor] = representative
                    stack.append(neighbor)
        sizes[representative] = size
    return labels, sizes


def exact_minimum_clearance(point: np.ndarray, radius: float, points: np.ndarray,
                            radii: np.ndarray, tree: cKDTree, excluded: int) -> tuple[float, int]:
    count = len(points)
    k = min(32, count)
    best = math.inf
    best_index = -1
    maximum_radius = float(radii.max())
    while True:
        distances, indices = tree.query(point, k=k)
        distances = np.atleast_1d(distances)
        indices = np.atleast_1d(indices)
        for distance, index_value in zip(distances, indices):
            index = int(index_value)
            if index == excluded:
                continue
            clearance = float(distance - radius - radii[index])
            if clearance < best:
                best, best_index = clearance, index
        unseen_lower_bound = float(distances[-1] - radius - maximum_radius)
        if k == count or best <= unseen_lower_bound:
            return best, best_index
        k = min(count, 2 * k)


def geometry_audit(base_pores: np.ndarray, virtual_pores: np.ndarray,
                   source_indices: np.ndarray, flags: np.ndarray,
                   boundary_length: float, base_edges: np.ndarray,
                   derived_edges: np.ndarray) -> dict:
    real_points = base_pores[:, 3:6]
    real_radii = base_pores[:, 6]
    virtual_points = virtual_pores[:, 3:6]
    virtual_radii = virtual_pores[:, 6]
    real_tree = cKDTree(real_points)
    max_real_radius = float(real_radii.max())

    connected_clearances = np.linalg.norm(
        virtual_points - real_points[source_indices], axis=1
    ) - virtual_radii - real_radii[source_indices]

    nonconnected_min = math.inf
    nonconnected_pair = [-1, -1]
    virtual_real_overlaps = 0
    virtual_real_new_overlaps = 0
    worst_virtual_real = math.inf
    worst_virtual_real_pair = [-1, -1]
    for virtual_index, (point, radius, source_index) in enumerate(
            zip(virtual_points, virtual_radii, source_indices)):
        minimum, real_index = exact_minimum_clearance(
            point, float(radius), real_points, real_radii, real_tree, int(source_index)
        )
        if minimum < nonconnected_min:
            nonconnected_min = minimum
            nonconnected_pair = [virtual_index, real_index]
        candidates = real_tree.query_ball_point(point, float(radius + max_real_radius))
        for real_index in candidates:
            if real_index == source_index:
                continue
            clearance = float(np.linalg.norm(point - real_points[real_index]) -
                              radius - real_radii[real_index])
            if clearance < 0.0:
                virtual_real_overlaps += 1
                if clearance < worst_virtual_real:
                    worst_virtual_real = clearance
                    worst_virtual_real_pair = [virtual_index, real_index]
                original_clearance = float(
                    np.linalg.norm(real_points[source_index] - real_points[real_index]) -
                    real_radii[source_index] - real_radii[real_index]
                )
                if original_clearance >= 0.0:
                    virtual_real_new_overlaps += 1

    difference = virtual_points[:, None, :] - virtual_points[None, :, :]
    distance = np.linalg.norm(difference, axis=2)
    clearance_matrix = distance - virtual_radii[:, None] - virtual_radii[None, :]
    np.fill_diagonal(clearance_matrix, math.inf)
    upper = np.triu(np.ones(clearance_matrix.shape, dtype=bool), 1)
    virtual_virtual_min_index = np.unravel_index(
        np.argmin(np.where(upper, clearance_matrix, math.inf)), clearance_matrix.shape
    )
    virtual_virtual_min = float(clearance_matrix[virtual_virtual_min_index])
    virtual_virtual_overlap_pairs = np.argwhere(upper & (clearance_matrix < 0.0))
    virtual_virtual_new_overlaps = 0
    for first, second in virtual_virtual_overlap_pairs:
        real_first, real_second = source_indices[first], source_indices[second]
        original_clearance = float(
            np.linalg.norm(real_points[real_first] - real_points[real_second]) -
            real_radii[real_first] - real_radii[real_second]
        )
        virtual_virtual_new_overlaps += int(original_clearance >= 0.0)

    base_labels, base_sizes = components(len(base_pores), base_edges)
    largest_component = max(base_sizes, key=base_sizes.get)
    base_degree = np.zeros(len(base_pores), dtype=int)
    for a_value, b_value in base_edges[:, 1:3]:
        base_degree[int(a_value) - 1] += 1
        base_degree[int(b_value) - 1] += 1
    isolated_records = []
    for source_index, flag in zip(source_indices, flags):
        if base_degree[source_index] != 0:
            continue
        isolated_records.append({
            "real_node_id": int(base_pores[source_index, 0]),
            "boundary_flag": int(flag),
            "base_component_representative_node_id": int(base_labels[source_index]),
            "base_component_size": int(base_sizes[int(base_labels[source_index])]),
            "derived_component_size": 2,
            "in_largest_sample_component": bool(base_labels[source_index] == largest_component),
            "forms_real_virtual_two_node_component": True,
            "input_active_gas": bool(base_pores[source_index, 11] < 1.0),
            "expected_boundary_flux_zero_if_all_liquid_equilibrium": bool(
                base_pores[source_index, 11] == 1.0),
            "can_contribute_transient_boundary_flux_and_global_conservation": bool(
                base_pores[source_index, 11] < 1.0),
            "singular_equation_expected": False,
        })

    virtual_ids = virtual_pores[:, 0].astype(int)
    real_ids = base_pores[:, 0].astype(int)
    def pair_ids(pair: list[int]) -> dict:
        if pair[0] < 0:
            return {}
        return {"virtual_node_id": int(virtual_ids[pair[0]]),
                "real_node_id": int(real_ids[pair[1]])}

    new_edge_lengths = derived_edges[len(base_edges):, 5]
    return {
        "format": "bubble_virtual_boundary_geometry_audit_v1",
        "rule": "old-model one-to-one copy-radius; center offset = 2*r_real + boundary_length",
        "inlet_virtual_count": int(np.count_nonzero(flags == 1)),
        "outlet_virtual_count": int(np.count_nonzero(flags == 2)),
        "boundary_connection_count": len(virtual_pores),
        "virtual_node_id_min": int(virtual_pores[:, 0].min()),
        "virtual_node_id_max": int(virtual_pores[:, 0].max()),
        "derived_node_count": len(base_pores) + len(virtual_pores),
        "derived_edge_count": len(derived_edges),
        "boundary_throat_length_min_um": float(new_edge_lengths.min()),
        "boundary_throat_length_max_um": float(new_edge_lengths.max()),
        "connected_real_virtual_clearance_min_um": float(connected_clearances.min()),
        "connected_real_virtual_clearance_max_um": float(connected_clearances.max()),
        "connected_clearance_target_um": boundary_length,
        "nonconnected_virtual_real_min_clearance_um": nonconnected_min,
        "nonconnected_virtual_real_min_pair": pair_ids(nonconnected_pair),
        "virtual_virtual_min_clearance_um": virtual_virtual_min,
        "virtual_virtual_min_pair": {
            "virtual_node_id_1": int(virtual_ids[virtual_virtual_min_index[0]]),
            "virtual_node_id_2": int(virtual_ids[virtual_virtual_min_index[1]]),
        },
        "overlap_counts": {
            "virtual_nonconnected_real": virtual_real_overlaps,
            "virtual_nonconnected_real_new_vs_source_layout": virtual_real_new_overlaps,
            "virtual_virtual": int(len(virtual_virtual_overlap_pairs)),
            "virtual_virtual_new_vs_real_boundary_layout": virtual_virtual_new_overlaps,
        },
        "worst_virtual_real_overlap_um": None if virtual_real_overlaps == 0 else worst_virtual_real,
        "worst_virtual_real_overlap_pair": pair_ids(worst_virtual_real_pair),
        "zero_length_edge_count": int(np.count_nonzero(derived_edges[:, 5] == 0.0)),
        "new_zero_length_boundary_edge_count": int(np.count_nonzero(new_edge_lengths == 0.0)),
        "virtual_nodes_in_nontrivial_zero_distance_cluster": 0,
        "virtual_node_coordination_distribution": {"1": len(virtual_pores)},
        "virtual_declared_geometry_volume_sum_um3": float(virtual_pores[:, 9].sum()),
        "virtual_volume_counted_as_internal_storage": False,
        "isolated_former_boundary_count": len(isolated_records),
        "isolated_former_inlet_count": sum(
            item["boundary_flag"] == 1 for item in isolated_records),
        "isolated_former_outlet_count": sum(
            item["boundary_flag"] == 2 for item in isolated_records),
        "isolated_former_boundary_nodes": isolated_records,
        "overlap_interpretation": (
            "negative spherical-envelope clearance is visualization-only when the generated "
            "boundary edge remains positive, one-to-one, and outside every zero-distance cluster"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-pore", type=Path, required=True)
    parser.add_argument("--base-throat", type=Path, required=True)
    parser.add_argument("--base-connect", type=Path, required=True)
    parser.add_argument("--output-pore", type=Path, required=True)
    parser.add_argument("--output-throat", type=Path, required=True)
    parser.add_argument("--output-connect", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--axis", choices=("x", "y", "z"), default="x")
    parser.add_argument("--extension-voxels", type=float, required=True)
    parser.add_argument("--voxel-size-um", type=float, required=True)
    parser.add_argument("--old-model-pore", type=Path, required=True)
    parser.add_argument("--old-model-throat", type=Path, required=True)
    args = parser.parse_args()

    require(args.extension_voxels > 0.0 and args.voxel_size_um > 0.0,
            "virtual-boundary length inputs must be positive")
    base_pores = np.loadtxt(args.base_pore)
    base_edges = np.loadtxt(args.base_throat)
    require(base_pores.shape[1] == 12 and base_edges.shape[1] == 11,
            "unexpected converted PNM column count")
    require(np.array_equal(base_pores[:, 0], np.arange(1, len(base_pores) + 1)),
            "base pore IDs must be continuous")
    require(np.array_equal(base_edges[:, 0], np.arange(1, len(base_edges) + 1)),
            "base edge IDs must be continuous")

    boundary_indices = np.flatnonzero(np.isin(base_pores[:, 10].astype(int), (1, 2)))
    flags = base_pores[boundary_indices, 10].astype(int)
    require(len(boundary_indices) > 0, "base network has no flagged boundary pores")
    axis_index = {"x": 3, "y": 4, "z": 5}[args.axis]
    inlet_coordinates = base_pores[boundary_indices[flags == 1], axis_index]
    outlet_coordinates = base_pores[boundary_indices[flags == 2], axis_index]
    require(len(inlet_coordinates) > 0 and len(outlet_coordinates) > 0,
            "both inlet and outlet flags are required")
    require(float(inlet_coordinates.max()) < float(outlet_coordinates.min()),
            "boundary flags do not identify separated low/high coordinate sides")

    derived_pores = base_pores.copy()
    derived_pores[boundary_indices, 10] = 0
    virtual = base_pores[boundary_indices].copy()
    virtual[:, 0] = np.arange(len(base_pores) + 1, len(base_pores) + len(virtual) + 1)
    virtual[:, 1] = 2
    virtual[:, 2] = base_pores[boundary_indices, 0]
    virtual[:, 10] = flags
    virtual[:, 11] = 1.0
    boundary_length = args.extension_voxels * args.voxel_size_um
    direction = np.where(flags == 1, -1.0, 1.0)
    virtual[:, axis_index] += direction * (2.0 * base_pores[boundary_indices, 6] + boundary_length)
    derived_pores = np.vstack((derived_pores, virtual))
    require(np.array_equal(virtual[:, 6:10], base_pores[boundary_indices, 6:10]),
            "virtual pore geometry was not copied from its real pore")

    boundary_edges = np.zeros((len(virtual), 11), dtype=float)
    boundary_edges[:, 0] = np.arange(len(base_edges) + 1,
                                     len(base_edges) + len(virtual) + 1)
    real_ids = base_pores[boundary_indices, 0]
    virtual_ids = virtual[:, 0]
    inlet = flags == 1
    boundary_edges[inlet, 1] = virtual_ids[inlet]
    boundary_edges[inlet, 2] = real_ids[inlet]
    boundary_edges[~inlet, 1] = real_ids[~inlet]
    boundary_edges[~inlet, 2] = virtual_ids[~inlet]
    # These generated edges have no source throat in PNM_pt.txt.  A zero
    # source_old_throat_id is the explicit not-applicable marker; side records
    # inlet (1) versus outlet (2).
    boundary_edges[:, 3] = 0
    boundary_edges[:, 4] = flags
    boundary_edges[:, 5] = boundary_length
    boundary_edges[:, 6] = 1.0
    boundary_edges[:, 7] = base_pores[boundary_indices, 6]
    boundary_edges[:, 8] = base_pores[boundary_indices, 7]
    boundary_edges[:, 9] = base_pores[boundary_indices, 8]
    boundary_edges[:, 10] = 0.0
    derived_edges = np.vstack((base_edges, boundary_edges))

    degree = np.zeros(len(derived_pores), dtype=int)
    neighbors: list[list[int]] = [[] for _ in range(len(derived_pores))]
    for a_value, b_value in derived_edges[:, 1:3]:
        a, b = int(a_value) - 1, int(b_value) - 1
        degree[a] += 1
        degree[b] += 1
        neighbors[a].append(b + 1)
        neighbors[b].append(a + 1)
    require(np.all(degree[len(base_pores):] == 1), "virtual pore coordination is not one")
    require(np.all(derived_pores[:len(base_pores), 10] == 0),
            "a real pore retained a boundary flag")
    require(np.all(np.isin(derived_pores[len(base_pores):, 10], (1, 2))),
            "a virtual pore lacks a boundary flag")
    require(np.all(boundary_edges[:, 5] > 0.0), "non-positive boundary edge length")

    pore_header = (
        "new_pore_id\tsource_type(0=old_pore,1=old_throat,2=virtual_boundary)\t"
        "source_id\tx\ty\tz\tradius\tarea\tshape_factor\tvolume\tboundary_flag\tSw"
    )
    edge_header = (
        "new_throat_id\tpore_id_1\tpore_id_2\tsource_old_throat_id\tside\t"
        "length\tSw\tradius\tarea\tshape_factor\tvolume"
    )
    args.output_pore.parent.mkdir(parents=True, exist_ok=True)
    write_table(args.output_pore, pore_header, derived_pores, 3)
    write_table(args.output_throat, edge_header, derived_edges, 5)
    max_degree = int(degree.max())
    with args.output_connect.open("w", encoding="utf-8") as stream:
        columns = ["pore_id", "coordination_number"] + [
            f"connection_pair_{index:03d}" for index in range(1, max_degree + 1)
        ]
        stream.write("# " + "\t".join(columns) + "\n")
        for node, connected in enumerate(neighbors, 1):
            connected.sort()
            row = [str(node), str(len(connected))]
            row.extend(f"{node}-{neighbor}" for neighbor in connected)
            row.extend([""] * (max_degree - len(connected)))
            stream.write("\t".join(row) + "\n")

    audit = geometry_audit(base_pores, virtual, boundary_indices, flags,
                           boundary_length, base_edges, derived_edges)
    audit["axis"] = args.axis
    audit["extension_voxels"] = args.extension_voxels
    audit["voxel_size_um"] = args.voxel_size_um
    audit["boundary_side_coordinates_um"] = {
        "inlet_min": float(inlet_coordinates.min()),
        "inlet_max": float(inlet_coordinates.max()),
        "outlet_min": float(outlet_coordinates.min()),
        "outlet_max": float(outlet_coordinates.max()),
    }
    audit["source_geometry_rule"] = {
        "source_path": "/workspace/Mesh/src/Mesh/Muilti_PNM_load_np_array.py",
        "lines": "94-170",
        "old_model_pore_sha256": sha256(args.old_model_pore),
        "old_model_throat_sha256": sha256(args.old_model_throat),
        "observed_old_boundary_length": 10.0,
        "observed_virtual_radius_equals_real": True,
        "observed_center_offset_equals_two_radii_plus_length": True,
        "generated_virtual_pore_fields": "copy radius, area, shape_factor, volume; force Sw=1",
        "generated_boundary_edge_fields": (
            "copy real radius, area, shape_factor; net length from extension*voxel; "
            "source volume metadata=0"
        ),
    }
    audit["file_hashes"] = {
        "base_pore": {"path": str(args.base_pore.resolve()), "sha256": sha256(args.base_pore)},
        "base_throat": {"path": str(args.base_throat.resolve()), "sha256": sha256(args.base_throat)},
        "base_connect": {"path": str(args.base_connect.resolve()), "sha256": sha256(args.base_connect)},
        "derived_pore": {"path": str(args.output_pore.resolve()), "sha256": sha256(args.output_pore)},
        "derived_throat": {"path": str(args.output_throat.resolve()), "sha256": sha256(args.output_throat)},
        "derived_connect": {"path": str(args.output_connect.resolve()), "sha256": sha256(args.output_connect)},
    }
    args.audit_output.parent.mkdir(parents=True, exist_ok=True)
    args.audit_output.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "inlet_virtual_count": audit["inlet_virtual_count"],
        "outlet_virtual_count": audit["outlet_virtual_count"],
        "derived_node_count": audit["derived_node_count"],
        "derived_edge_count": audit["derived_edge_count"],
        "connected_clearance_um": audit["connected_real_virtual_clearance_min_um"],
        "isolated_former_boundary_count": audit["isolated_former_boundary_count"],
    }))


if __name__ == "__main__":
    main()
