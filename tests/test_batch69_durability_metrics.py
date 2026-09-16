#!/usr/bin/env python3
import importlib.util
import math
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/analyze_batch69_durability.py"
spec = importlib.util.spec_from_file_location("durability", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_identity_parser_preserves_unverified_e_label():
    parsed = module.parse_case_identity("sands-of-100B-mu-0.5-20kpa-e-0.52")
    assert parsed["gradation"] == "100B"
    assert parsed["mu"] == 0.5
    assert parsed["stress_kpa"] == 20
    assert parsed["void_ratio_label"] == "0.52"
    assert parsed["void_ratio_label_verified"] is False


def test_auc_and_interpolation_use_common_pv_window():
    x = np.array([1.0, 2.0, 5.0])
    retention = np.array([0.8, 0.6, 0.0])
    assert math.isclose(module.normalized_auc(x, retention, 5.0), 0.5, rel_tol=1e-12)
    values = module.interpolate([0, 1, 2], [1, .5, 0], [.5, 1.5, 3])
    assert np.allclose(values[:2], [.75, .25])
    assert math.isnan(values[2])


def test_entropy_and_gini_have_expected_limits():
    uniform = module.distribution_metrics([1, 1, 1, 1])
    assert math.isclose(uniform["gas_entropy_normalized"], 1.0)
    assert math.isclose(uniform["gas_effective_count_simpson"], 4.0)
    assert math.isclose(uniform["gas_gini"], 0.0)
    concentrated = module.distribution_metrics([1000, 1, 1, 1])
    assert concentrated["gas_entropy_normalized"] < uniform["gas_entropy_normalized"]
    assert concentrated["gas_gini"] > 0.7


def test_weibull_recovers_synthetic_parameters():
    pv = np.linspace(0.01, 20, 300)
    retention = module.weibull_retention(pv, 8.0, 1.4)
    fit = module.fit_weibull(pv, retention)
    assert math.isclose(fit["weibull_scale_pv"], 8.0, rel_tol=2e-3)
    assert math.isclose(fit["weibull_shape"], 1.4, rel_tol=2e-3)
    assert fit["weibull_r2"] > 0.9999


def test_plot_bundle_accepts_full_pv50_field_name():
    case = pd.DataFrame({
        "variant": ["PNM-0.6", "PNM-0.8"],
        "control_volume_porosity": [.2, .2],
        "pv_to_50pct_free_gas_loss": [7., 8.],
    })
    curves = pd.DataFrame({
        "variant": ["PNM-0.6", "PNM-0.6", "PNM-0.8", "PNM-0.8"],
        "run_id": ["a", "a", "b", "b"], "pv": [0., 1., 0., 1.],
        "free_gas_retained_fraction": [1., .8, 1., .85],
    })
    paired = pd.DataFrame({
        "parent_case": ["p"], "pv_to_50pct_free_gas_loss_p06": [7.],
        "pv_to_50pct_free_gas_loss_p08": [8.],
    })
    with tempfile.TemporaryDirectory() as directory:
        tmp_path = Path(directory)
        module.make_plots(case, curves, paired, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), tmp_path)
        assert (tmp_path / "01_free_gas_retention_common_pv.png").is_file()
        assert (tmp_path / "02_paired_pv50.png").is_file()
        assert (tmp_path / "03_porosity_vs_pv50.png").is_file()


def test_segment_spatial_moments_have_known_uniform_limits():
    result = module.segment_spatial_moments(np.ones(20), np.ones(20), np.ones(20))
    assert math.isclose(result["segment_gas_spatial_entropy_normalized"], 1.0)
    assert math.isclose(result["gas_center_x_normalized_20segment"], 0.5)
    assert result["gas_longitudinal_variance_normalized_20segment"] > 0


if __name__ == "__main__":
    for name, value in sorted(globals().copy().items()):
        if name.startswith("test_"):
            value()
    print("batch69 durability metric tests passed")
