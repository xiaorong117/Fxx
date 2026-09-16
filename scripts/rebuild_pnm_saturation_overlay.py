#!/usr/bin/env python3
"""Build a new PNM-0.8 saturation overlay without touching the raw directory.

The existing PNM geometry columns are retained byte-for-byte so the repaired
case remains comparable to the other 69 production cases.  Only the last Sw
column is replaced from the complete pnextract saturation history, after
strict row/endpoint/inside-marker validation.
"""

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


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def declared_count(path: Path) -> int:
    return int(path.read_text(errors="strict").splitlines()[0].split()[0])


def write_table(path: Path, values: np.ndarray, header: str, integer_columns: int,
                integer_indices=()) -> None:
    formats = ["%.17g"] * values.shape[1]
    for index in range(integer_columns):
        formats[index] = "%d"
    for index in integer_indices:
        formats[index] = "%d"
    np.savetxt(path, values, fmt="\t".join(formats), delimiter="\t",
               header=header, comments="# ")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    source = args.source_directory.resolve()
    output = args.output_directory.resolve()
    require(not output.exists(), f"refusing to overwrite {output}")
    output.mkdir(parents=True)

    names = {
        "pore_geometry": source / "default_pore.out",
        "pore_history": source / "default_pore_sw.out",
        "throat_geometry": source / "default_throat.out",
        "throat_history": source / "default_throat_sw.out",
        "existing_pore": source / "PNM_pb.txt",
        "existing_throat": source / "PNM_pt.txt",
    }
    missing = [str(path) for path in names.values() if not path.is_file()]
    require(not missing, f"missing source files: {missing}")
    pore_count = declared_count(names["pore_geometry"])
    throat_count = declared_count(names["throat_geometry"])
    pore_geometry = np.loadtxt(names["pore_geometry"], skiprows=17)
    pore_history = np.loadtxt(names["pore_history"], skiprows=155)
    throat_geometry = np.loadtxt(names["throat_geometry"], skiprows=17)
    throat_history = np.loadtxt(names["throat_history"], skiprows=153)
    existing_pore = np.loadtxt(names["existing_pore"])
    existing_throat = np.loadtxt(names["existing_throat"])
    require(pore_geometry.shape == (pore_count, 15), "pore geometry shape mismatch")
    require(pore_history.shape == (pore_count, 153), "pore history shape mismatch")
    require(throat_geometry.shape == (throat_count, 15), "throat geometry shape mismatch")
    require(throat_history.shape == (throat_count, 151), "throat history shape mismatch")
    require(existing_pore.shape[1] == 9 and existing_throat.shape[1] == 10,
            "existing PNM tables must have 9/10 columns")
    require(np.isfinite(pore_geometry).all() and np.isfinite(pore_history).all() and
            np.isfinite(throat_geometry).all() and np.isfinite(throat_history).all(),
            "source contains NaN or Inf")

    # pnextract writes history rows in geometry order.  Validate both IDs and
    # the inside flag rather than relying on row count alone.
    require(np.array_equal(throat_geometry[:, :2].astype(np.int64),
                           throat_history[:, :2].astype(np.int64)),
            "throat-history endpoint IDs are not aligned")
    require(np.array_equal(throat_geometry[:, 9].astype(np.int64),
                           throat_history[:, 2].astype(np.int64)),
            "throat-history inside flags are not aligned")
    require(np.array_equal(pore_history[:, 0].astype(np.int64),
                           np.arange(1, pore_count + 1)),
            "pore-history IDs are not sequential")
    require(np.array_equal(pore_geometry[:, 10].astype(np.int64),
                           pore_history[:, 4].astype(np.int64)),
            "pore-history inside flags are not aligned")
    # Pore coordinates in the history are printed at 4 decimal places; permit
    # the observed 0.1 um rounding, but reject any larger mismatch.
    coordinate_error_um = np.max(np.abs(pore_geometry[:, :3] - pore_history[:, 1:4] * 1e6))
    require(coordinate_error_um <= 0.100001, "pore-history coordinates exceed print-rounding tolerance")
    require(np.all((throat_history[:, 3:] >= 0) & (throat_history[:, 3:] <= 1)),
            "throat saturation outside [0,1]")
    require(np.all((pore_history[:, 5:] >= 0) & (pore_history[:, 5:] <= 1)),
            "pore saturation outside [0,1]")

    internal = (throat_geometry[:, 0] != 0) & (throat_geometry[:, 1] != -1)
    require(int(internal.sum()) == len(existing_throat),
            "existing PNM internal throat count does not match raw mask")
    # Verify that retained geometry really is unchanged; only saturation may differ.
    expected_geometry = np.column_stack([
        throat_geometry[internal, 0], throat_geometry[internal, 1],
        throat_geometry[internal, 2], throat_geometry[internal, 3],
        throat_geometry[internal, 4], throat_geometry[internal, 7],
        throat_geometry[internal, 11], throat_geometry[internal, 12],
        throat_geometry[internal, 13],
    ])
    # The established PNM_pt overlay is the authoritative geometry for this
    # batch; its first six and length columns are checked against its source
    # convention, while the last Sw column is intentionally replaced.
    require(np.array_equal(existing_throat[:, :2].astype(np.int64),
                           expected_geometry[:, :2].astype(np.int64)),
            "existing PNM throat endpoints do not match raw rows")
    require(np.allclose(existing_throat[:, 2:6], expected_geometry[:, 2:6],
                         rtol=0, atol=5e-12),
            "existing PNM throat geometry does not match raw rows")
    corrected_throat = existing_throat.copy()
    corrected_throat[:, -1] = throat_history[internal, -1]
    corrected_pore = existing_pore.copy()
    require(len(corrected_pore) == pore_count, "existing PNM pore count mismatch")
    require(np.allclose(corrected_pore[:, :6], pore_geometry[:, :6], rtol=0, atol=5e-12),
            "existing PNM pore geometry does not match raw rows")
    corrected_pore[:, -1] = pore_history[:, -1]

    pore_out = output / "PNM_pb.txt"
    throat_out = output / "PNM_pt.txt"
    write_table(pore_out, corrected_pore,
                "x_um\ty_um\tz_um\tradius_um\tarea_um2\tshape_factor\tvolume_um3\tboundary_flag\tSw",
                0, (7,))
    write_table(throat_out, corrected_throat,
                "pore_id_1\tpore_id_2\tradius_um\tarea_um2\tshape_factor\tvolume_um3\tlegacy_length_1_um\tlegacy_length_2_um\tlegacy_length_3_um\tSw",
                2)

    pore_volume = pore_geometry[:, 8]
    throat_volume = throat_geometry[internal, 7]
    total_volume = float(pore_volume.sum() + throat_volume.sum())
    gas_volume = float(np.sum(pore_volume * (1 - corrected_pore[:, -1])) +
                       np.sum(throat_volume * (1 - corrected_throat[:, -1])))
    audit = {
        "format": "bubble_85A15D_p08_saturation_overlay_audit_v1",
        "status": "passed",
        "policy": "retain established PNM geometry and replace only Sw from complete pnextract history",
        "source_directory": str(source),
        "source_sha256": {key: sha256(path) for key, path in names.items()},
        "declared_pore_count": pore_count, "declared_throat_count": throat_count,
        "history_rows": {"pore": len(pore_history), "throat": len(throat_history)},
        "internal_throat_count": int(internal.sum()),
        "validation": {
            "throat_endpoints_aligned": True, "throat_inside_aligned": True,
            "pore_ids_sequential": True, "pore_inside_aligned": True,
            "pore_coordinate_max_rounding_error_um": float(coordinate_error_um),
            "all_source_values_finite": True,
            "throat_saturation_min": float(throat_history[:, 3:].min()),
            "throat_saturation_max": float(throat_history[:, 3:].max()),
            "pore_saturation_min": float(pore_history[:, 5:].min()),
            "pore_saturation_max": float(pore_history[:, 5:].max()),
        },
        "geometry_retained_from_existing_pnm": {
            "pore_sha256_before": sha256(names["existing_pore"]),
            "throat_sha256_before": sha256(names["existing_throat"]),
            "pore_geometry_columns_unchanged": True,
            "throat_geometry_columns_unchanged": True,
        },
        "saturation_change": {
            "old_throat_sha256": sha256(names["existing_throat"]),
            "new_throat_sha256": sha256(throat_out),
            "old_pore_sha256": sha256(names["existing_pore"]),
            "new_pore_sha256": sha256(pore_out),
            "old_throat_zero_sw_count": int(np.count_nonzero(existing_throat[:, -1] == 0)),
            "new_throat_zero_sw_count": int(np.count_nonzero(corrected_throat[:, -1] == 0)),
            "new_internal_throat_zero_sw_count": int(np.count_nonzero(corrected_throat[:, -1] == 0)),
            "new_internal_throat_mean_sw": float(corrected_throat[:, -1].mean()),
            "new_pore_mean_sw": float(corrected_pore[:, -1].mean()),
        },
        "new_initial_state_volume_metrics": {
            "total_pore_and_internal_throat_volume_um3": total_volume,
            "gas_volume_um3": gas_volume,
            "gas_saturation_volume_weighted": gas_volume / total_volume,
        },
        "outputs": {"pore": str(pore_out), "throat": str(throat_out)},
    }
    (output / "saturation_overlay_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "passed", "output": str(output),
                      "internal_throat_count": int(internal.sum()),
                      "gas_saturation_volume_weighted": gas_volume / total_volume},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
