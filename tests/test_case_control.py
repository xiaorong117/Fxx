#!/usr/bin/env python3
import json
import sys
import tempfile
from pathlib import Path

project = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(project / "scripts"))
import case_control  # noqa: E402

with tempfile.TemporaryDirectory(prefix="bubble_control_test_") as temporary:
    root = Path(temporary) / "production"; root.mkdir()
    case = root / "case_01"; (case / "control").mkdir(parents=True)
    config = Path(temporary) / "base.json"
    original = {
        "io": {"output_directory": str(case), "checkpoint_interval_steps": 100,
               "vtk_interval_steps": 1000},
        "time": {"dt_min_s": 1e-6, "dt_max_s": 0.1},
        "time_step_control": {"dt_min_s": 1e-6, "dt_max_s": 0.1},
    }
    config.write_text(json.dumps(original))
    (case / "control" / "active_config.json").write_text(
        json.dumps({"path": str(config)}))

    staged = case_control.stage_settings("case_01", {
        "dt_max_s": 10.0,
        "checkpoint_interval_steps": 25,
        "vtk_interval_steps": 250,
    }, root)
    revised = json.loads(Path(staged["path"]).read_text())
    assert json.loads(config.read_text()) == original
    assert revised["time"]["dt_max_s"] == 10.0
    assert revised["time_step_control"]["dt_max_s"] == 10.0
    assert revised["io"]["checkpoint_interval_steps"] == 25
    assert revised["io"]["vtk_interval_steps"] == 250
    assert (case / "control" / "control_history.jsonl").is_file()

    for invalid in ("../case_01", "case_01/child", "", "x" * 129):
        try:
            case_control.resolve_case(invalid, root)
        except case_control.ControlError:
            pass
        else:
            raise SystemExit("path traversal or invalid identifier accepted")

    try:
        case_control.stage_settings("case_01", {"mass_transfer_multiplier": 2}, root)
    except case_control.ControlError:
        pass
    else:
        raise SystemExit("unsafe physical parameter was accepted for hot restart")

    try:
        case_control.restart_case("case_01", Path("/bin/true"), root)
    except case_control.ControlError as error:
        assert "checkpoint" in str(error)
    else:
        raise SystemExit("restart without checkpoint was accepted")

    (case / "checkpoint.json").write_text(json.dumps({
        "format": "bubble_checkpoint_v6", "step_number": 12, "time_s": 1.5
    }))
    restarted = case_control.restart_case("case_01", Path("/bin/true"), root)
    assert restarted["checkpoint_step"] == 12
    assert restarted["staged_revision"] == 1
    assert restarted["tmux_session"] == "bubble_ctl_case_01"
    assert (case / "control" / "active_config.json").is_file()
    assert (case / "control" / "last_applied_config.json").is_file()

print("PASS restricted dashboard control, revision audit and path confinement")
