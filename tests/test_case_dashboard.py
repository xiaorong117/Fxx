#!/usr/bin/env python3
import json
import subprocess
import sys
from pathlib import Path

project = Path(sys.argv[1]).resolve()
root = project / "output" / "production"
result = subprocess.run([sys.executable, str(project / "scripts/case_dashboard.py"),
                         "--production-root", str(root), "--snapshot"],
                        check=True, capture_output=True, text=True)
payload = json.loads(result.stdout)
assert Path(payload["production_root"]).resolve() == root.resolve()
assert payload["cases"], "production case list is empty"
preview = next((x for x in payload["cases"] if x["case_name"] ==
                "bubble_dissolution_accelerated_preview_kLa100"), None)
assert preview is not None
assert preview["kind"] == "accelerated_preview"
assert preview["config_path"].endswith("bubble_solver.accelerated_preview.json")
assert preview["effective"]["backend"] == "eigen_bicgstab_ilut"
assert all(Path(case["directory"]).resolve().is_relative_to(root.resolve())
           for case in payload["cases"])
sys.path.insert(0, str(project / "scripts"))
import case_dashboard  # noqa: E402
preview_dir = Path(preview["directory"])
detail = case_dashboard.snapshot_case(
    preview_dir, root, Path(preview["config_path"]),
    case_dashboard.bubble_process_table())
assert detail["physics"]["mass_transfer_multiplier"] == 100
assert any(x["name"] == "mass_transfer_multiplier" and
           x["classification"] == "preview-only" for x in detail["parameters"])
assert detail["geometry_network"]["segments"] == 20
if (detail["last_diagnostics"] or {}).get("step"):
    assert detail["history_series"], "history series missing despite diagnostics"
batch = [case for case in payload["cases"] if case.get("variant")]
assert len(batch) == 70
assert sum(case["variant"] == "PNM-0.6" for case in batch) == 35
assert sum(case["variant"] == "PNM-0.8" for case in batch) == 35
batch_statuses = payload["summary"]["variant_status_counts"]
assert sum(batch_statuses[variant][status] for variant in ("PNM-0.6", "PNM-0.8")
           for status in ("running", "queued", "completed", "failed", "interrupted")) == 70
source_text = (project / "scripts/case_dashboard.py").read_text()
for marker in ("chart-inventory", "chart-saturation", "chart-timestep",
               "chart-quality", "setInterval", "/api/case", "batch-summary",
               "variant-filter", "status-filter", "loadDetail", "event-summary",
               "batched_microbubble", "批量失败 / Batch failed", "statusBi"):
    assert marker in source_text, f"dashboard chart marker missing: {marker}"
assert "document.querySelectorAll('.badge.running" not in source_text
print("PASS scalable 70-case dashboard and on-demand detail confinement")
