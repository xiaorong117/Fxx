#!/usr/bin/env python3
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
template = Path(sys.argv[2])
doc = json.loads(template.read_text())
if doc.get("format") != "bubble_geometry_correction_overlay_v1":
    raise SystemExit("wrong correction overlay format")
expected_hashes = {
    name: hashlib.sha256((root / name).read_bytes()).hexdigest()
    for name in ("PNM_pb.txt", "PNM_pt.txt")
}
if doc.get("input_sha256") != expected_hashes:
    raise SystemExit("correction template input hashes are stale")
items = doc.get("corrections", [])
if len(items) != 516 or len({(x["file"], x["row_id"], x["field"]) for x in items}) != 516:
    raise SystemExit("correction template count/uniqueness mismatch")
if any(x["replacement_value"] is not None or x["source_note"] is not None for x in items):
    raise SystemExit("unapproved values found in required-correction template")
counts = {}
for item in items:
    key = (item["file"], item["field"])
    counts[key] = counts.get(key, 0) + 1
expected_counts = {
    ("PNM_pb.txt", "volume"): 1,
    ("PNM_pt.txt", "volume"): 1,
    ("PNM_pt.txt", "id1_half_length"): 512,
    ("PNM_pt.txt", "id2_half_length"): 2,
}
if counts != expected_counts:
    raise SystemExit(f"unexpected raw correction categories: {counts}")
print("PASS optional correction overlay remains untouched: 516 null raw fields")
