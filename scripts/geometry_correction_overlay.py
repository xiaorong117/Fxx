#!/usr/bin/env python3
"""Create or apply an explicit, provenance-bearing raw-geometry correction overlay.

The raw PNM files are never modified.  Applying an overlay requires exactly
one correction for every non-positive storage volume and half-throat length,
no corrections to already-positive fields, matching source-file SHA-256, a
finite positive replacement, and a non-empty source note.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


FIELDS = {
    "PNM_pb.txt": {"volume": 6},
    "PNM_pt.txt": {"volume": 5, "id1_half_length": 7, "id2_half_length": 8},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_lines(path: Path) -> tuple[list[str], list[list[str]]]:
    raw = path.read_text(encoding="utf-8-sig").splitlines()
    fields: list[list[str]] = []
    for line_number, line in enumerate(raw, 1):
        row = line.split()
        if not row:
            raise ValueError(f"{path}:{line_number}: blank rows are unsupported")
        try:
            [float(value) for value in row]
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: non-numeric raw row") from error
        fields.append(row)
    return raw, fields


def defect_keys(root: Path) -> tuple[dict[str, list[list[str]]], set[tuple[str, int, str]]]:
    tables: dict[str, list[list[str]]] = {}
    defects: set[tuple[str, int, str]] = set()
    expected_columns = {"PNM_pb.txt": 9, "PNM_pt.txt": 10}
    for filename, mapped_fields in FIELDS.items():
        _, rows = load_lines(root / filename)
        tables[filename] = rows
        for row_id, row in enumerate(rows, 1):
            if len(row) != expected_columns[filename]:
                raise ValueError(
                    f"{filename}:{row_id}: expected {expected_columns[filename]} columns"
                )
            for field, column in mapped_fields.items():
                value = float(row[column])
                if not math.isfinite(value):
                    raise ValueError(f"{filename}:{row_id}:{field} is non-finite")
                if value <= 0.0:
                    defects.add((filename, row_id, field))
    return tables, defects


def make_template(root: Path, destination: Path) -> None:
    _, defects = defect_keys(root)
    document = {
        "format": "bubble_geometry_correction_overlay_v1",
        "input_sha256": {
            name: sha256(root / name) for name in sorted(FIELDS)
        },
        "units": {
            "PNM_pb.txt.volume": "um3",
            "PNM_pt.txt.volume": "um3",
            "PNM_pt.txt.id1_half_length": "um",
            "PNM_pt.txt.id2_half_length": "um",
        },
        "corrections": [
            {
                "file": filename,
                "row_id": row_id,
                "field": field,
                "replacement_value": None,
                "source_note": None,
            }
            for filename, row_id, field in sorted(defects)
        ],
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE TEMPLATE {destination} with {len(defects)} required corrections")


def apply_overlay(root: Path, overlay_path: Path, output_dir: Path) -> None:
    tables, required = defect_keys(root)
    document = json.loads(overlay_path.read_text(encoding="utf-8"))
    if document.get("format") != "bubble_geometry_correction_overlay_v1":
        raise ValueError("unsupported overlay format")
    for filename in FIELDS:
        expected = sha256(root / filename)
        if document.get("input_sha256", {}).get(filename) != expected:
            raise ValueError(f"source checksum mismatch for {filename}")

    supplied: dict[tuple[str, int, str], tuple[float, str]] = {}
    for item in document.get("corrections", []):
        key = (item.get("file"), item.get("row_id"), item.get("field"))
        if key in supplied:
            raise ValueError(f"duplicate correction: {key}")
        if key not in required:
            raise ValueError(f"correction is not a current non-positive field: {key}")
        value = item.get("replacement_value")
        note = item.get("source_note")
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"replacement must be finite and positive: {key}")
        if not isinstance(note, str) or not note.strip():
            raise ValueError(f"source_note is required: {key}")
        supplied[key] = (float(value), note.strip())
    missing = sorted(required - supplied.keys())
    if missing:
        raise ValueError(f"overlay is missing {len(missing)} corrections; first={missing[:5]}")

    output_dir.mkdir(parents=True, exist_ok=False)
    ledger_corrections = []
    for filename, rows in tables.items():
        for row_id, row in enumerate(rows, 1):
            for field, column in FIELDS[filename].items():
                key = (filename, row_id, field)
                if key not in supplied:
                    continue
                old_text = row[column]
                value, note = supplied[key]
                row[column] = format(value, ".17g")
                ledger_corrections.append({
                    "file": filename, "row_id": row_id, "field": field,
                    "old_text": old_text, "replacement_value": value,
                    "source_note": note,
                })
        (output_dir / filename).write_text(
            "\n".join("\t".join(row) for row in rows) + "\n", encoding="utf-8"
        )
    ledger = {
        "format": "bubble_geometry_correction_ledger_v1",
        "overlay_sha256": sha256(overlay_path),
        "input_sha256": {name: sha256(root / name) for name in sorted(FIELDS)},
        "output_sha256": {name: sha256(output_dir / name) for name in sorted(FIELDS)},
        "correction_count": len(ledger_corrections),
        "corrections": ledger_corrections,
    }
    (output_dir / "geometry_correction_ledger.json").write_text(
        json.dumps(ledger, indent=2) + "\n", encoding="utf-8"
    )
    print(f"APPLIED {len(ledger_corrections)} corrections to {output_dir}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write-template", type=Path)
    action.add_argument("--overlay", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    try:
        if args.write_template:
            if args.output_dir:
                raise ValueError("--output-dir is only valid with --overlay")
            make_template(args.root, args.write_template)
        else:
            if not args.output_dir:
                raise ValueError("--output-dir is required with --overlay")
            apply_overlay(args.root, args.overlay, args.output_dir)
        return 0
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
