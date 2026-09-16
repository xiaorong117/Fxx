#!/usr/bin/env python3
"""Validate that client/Codex delivery documentation is complete and secret-free."""

from __future__ import annotations

import re
import sys
from pathlib import Path


project = Path(sys.argv[1]).resolve()
manual = project / "docs" / "CLIENT_DELIVERY_GUIDE_ZH.md"
agents = project / "AGENTS.md"
for path in (manual, agents):
    if not path.is_file() or path.stat().st_size < 4000:
        raise SystemExit(f"missing or incomplete delivery document: {path}")

equations_doc = project / "docs" / "静止气泡溶解_液相对流扩散孔隙网络模型_全耦合控制方程_实现更新版_20260912.docx"
if not equations_doc.is_file() or equations_doc.stat().st_size < 40000:
    raise SystemExit(f"missing updated control-equations DOCX: {equations_doc}")

combined = manual.read_text() + "\n" + agents.read_text()
required = (
    "6ea4e9c1a3e00f3919e2415edeaa49e527bbd6f7c9042bbb98d7ee44695ebd04",
    "batch35_flow_pe555_pardiso_v3/batch_cases.json",
    "batch_progress.json", "diagnostics.tsv", "global_history.tsv",
    "segment_history.tsv", "gas_disappearance_events.tsv", "checkpoint.json",
    "run_manifest.failed.json", "mass_transfer_multiplier", "QIN-container",
    "SSH", "flock", "tmux", "PARDISO", "AMGX", "t=0", "未实验确认",
)
for marker in required:
    if marker not in combined:
        raise SystemExit(f"delivery documentation missing marker: {marker}")

for relative in ("PHYSICS_IMPLEMENTATION_ZH.md", "../AGENTS.md"):
    target = (manual.parent / relative).resolve()
    if not target.is_file():
        raise SystemExit(f"broken required documentation link: {relative}")

for forbidden in (
        r"-----BEGIN (?:OPENSSH|RSA|EC) PRIVATE KEY-----",
        r"X-Control-Token\s*[:=]\s*[A-Za-z0-9_-]{16,}",
        r"\.dashboard-control-token\s*[:=]\s*[A-Fa-f0-9]{32,}"):
    if re.search(forbidden, combined):
        raise SystemExit("delivery documentation contains credential material")

from zipfile import ZipFile
from xml.etree import ElementTree as ET
with ZipFile(equations_doc) as archive:
    for part in ("[Content_Types].xml", "_rels/.rels", "word/document.xml",
                 "word/styles.xml", "word/settings.xml", "docProps/core.xml"):
        ET.fromstring(archive.read(part))
    document_text = ''.join(ET.fromstring(archive.read("word/document.xml")).itertext())
    if "27. 守恒型微小气泡批量退休" not in document_text or \
            "30. checkpoint、输出与可复现性" not in document_text:
        raise SystemExit("updated DOCX appendix is incomplete")

print("PASS client and Codex delivery documentation")
