#!/usr/bin/env python3
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

script = Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory(prefix="bubble_overlay_test_") as temp:
    root = Path(temp)
    pore = root / "PNM_pb.txt"
    throat = root / "PNM_pt.txt"
    pore.write_text("0 0 0 2 12 0.05 0 0 1\n")
    throat.write_text("1 1 1 3 0.05 0 0 0 2 1\n")
    originals = {
        pore.name: hashlib.sha256(pore.read_bytes()).hexdigest(),
        throat.name: hashlib.sha256(throat.read_bytes()).hexdigest(),
    }
    overlay = root / "overlay.json"
    subprocess.run(
        [sys.executable, str(script), "--root", str(root),
         "--write-template", str(overlay)], check=True
    )
    doc = json.loads(overlay.read_text())
    if len(doc["corrections"]) != 3:
        raise SystemExit("synthetic overlay defect count mismatch")
    replacements = {
        ("PNM_pb.txt", "volume"): 4.0,
        ("PNM_pt.txt", "volume"): 5.0,
        ("PNM_pt.txt", "id1_half_length"): 1.0,
    }
    for item in doc["corrections"]:
        item["replacement_value"] = replacements[(item["file"], item["field"])]
        item["source_note"] = "deterministic overlay unit-test fixture"
    overlay.write_text(json.dumps(doc, indent=2) + "\n")
    output = root / "derived"
    subprocess.run(
        [sys.executable, str(script), "--root", str(root), "--overlay", str(overlay),
         "--output-dir", str(output)], check=True
    )
    if any(hashlib.sha256((root / name).read_bytes()).hexdigest() != digest
           for name, digest in originals.items()):
        raise SystemExit("overlay modified raw input")
    corrected_pore = [float(x) for x in (output / "PNM_pb.txt").read_text().split()]
    corrected_throat = [float(x) for x in (output / "PNM_pt.txt").read_text().split()]
    if corrected_pore[6] != 4.0 or corrected_throat[5] != 5.0 or corrected_throat[7] != 1.0:
        raise SystemExit("overlay replacements were not written exactly")
    ledger = json.loads((output / "geometry_correction_ledger.json").read_text())
    if ledger["correction_count"] != 3 or ledger["input_sha256"] != originals:
        raise SystemExit("overlay ledger mismatch")
    for name, digest in ledger["output_sha256"].items():
        if hashlib.sha256((output / name).read_bytes()).hexdigest() != digest:
            raise SystemExit("overlay output checksum mismatch")
print("PASS geometry overlay end-to-end without raw mutation")
