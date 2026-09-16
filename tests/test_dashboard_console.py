#!/usr/bin/env python3
import sys
import tempfile
from pathlib import Path

project = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(project / "scripts"))
import case_dashboard  # noqa: E402

with tempfile.TemporaryDirectory(prefix="bubble_console_test_") as temporary:
    root = Path(temporary) / "production"
    case = root / "case_01"
    case.mkdir(parents=True)
    lines = [f"step={index} time_s={index / 10}" for index in range(250)]
    (case / "run.log").write_text("\n".join(lines) + "\n")
    result = case_dashboard.console_tail(root, "case_01", 25, 4096)
    assert result["available"] is True
    assert result["source"] == "run.log"
    assert len(result["lines"]) == 25
    assert result["lines"][-1].startswith("step=249 ")
    assert result["lines"][0].startswith("step=225 ")
    for invalid in ("../case_01", "case_01/child", "", "/tmp"):
        try:
            case_dashboard.console_tail(root, invalid)
        except ValueError:
            pass
        else:
            raise SystemExit("console path traversal was accepted")
    nested = root / "batch" / "case_02"
    nested.mkdir(parents=True)
    (nested / "driver.log").write_text("queued\nstep=1 time_s=0.1\n")
    result = case_dashboard.console_tail(root, "batch/case_02", 10, 4096)
    assert result["available"] is True
    assert result["source"] == "driver.log"
    assert result["lines"][-1] == "step=1 time_s=0.1"

source = (project / "scripts/case_dashboard.py").read_text()
for marker in ("/api/console", "capture-pane", "solver-console",
               "Live tmux solver console", "setInterval(refreshConsole,2000)",
               "requestDirectory=selectedDirectory",
               "selectedDirectory!==requestDirectory"):
    assert marker in source, f"live console marker missing: {marker}"

print("PASS bounded live-console tail and production-root confinement")
