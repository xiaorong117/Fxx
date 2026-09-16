#!/usr/bin/env python3
"""Reproducible, read-only durability analysis for the valid batch-35 runs.

The script never changes solver outputs.  It excludes inputs rejected by the
saturation-source audit, derives case/pair/geometry/segment metrics, and writes
all analysis products under a separate output directory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import struct
import sys
import zlib
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import friedmanchisquare, rankdata, spearmanr, wilcoxon


SCRIPT_VERSION = "1.0.0"
LOSS_LEVELS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90)
PV_GRID_POINTS = 201
SEGMENT_PV_TARGETS = (0.0, 1.0, 5.0, 10.0, 20.0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def finite(value, default=math.nan):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def parse_case_identity(parent_case: str) -> dict:
    pattern = re.compile(
        r"^sands-of-(?P<gradation>.+?)-mu-(?P<mu>[0-9.]+)-"
        r"(?P<stress>[0-9.]+)kpa-e-(?P<void_ratio_label>.*)$"
    )
    match = pattern.match(parent_case)
    if not match:
        raise ValueError(f"cannot parse parent case identity: {parent_case}")
    result = match.groupdict()
    result["mu"] = float(result["mu"])
    result["stress_kpa"] = float(result.pop("stress"))
    # Several folder names contain an empty or zero e suffix.  Preserve it as
    # an unverified label; measured PNM porosity is used in all calculations.
    result["void_ratio_label_verified"] = False
    return result


def strictly_increasing_last(x: np.ndarray, *ys: np.ndarray):
    """Keep the last row at duplicate abscissae and return sorted arrays."""
    order = np.argsort(x, kind="stable")
    x = np.asarray(x, dtype=float)[order]
    sorted_y = [np.asarray(y, dtype=float)[order] for y in ys]
    keep = np.r_[x[1:] != x[:-1], True]
    return (x[keep], *[y[keep] for y in sorted_y])


def interpolate(x, y, targets):
    x, y = strictly_increasing_last(np.asarray(x), np.asarray(y))
    targets = np.asarray(targets, dtype=float)
    result = np.full(targets.shape, np.nan)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 2:
        return result
    inside = (targets >= x[0]) & (targets <= x[-1])
    result[inside] = np.interp(targets[inside], x, y)
    return result


def normalized_auc(x, retention, horizon):
    if horizon <= 0:
        raise ValueError("AUC horizon must be positive")
    grid = np.linspace(0.0, horizon, PV_GRID_POINTS)
    values = interpolate(
        np.r_[0.0, np.asarray(x, dtype=float)],
        np.r_[1.0, np.asarray(retention, dtype=float)],
        grid,
    )
    if np.any(~np.isfinite(values)):
        return math.nan
    return float(np.trapezoid(values, grid) / horizon)


def weibull_retention(pv, scale, shape):
    return np.exp(-np.power(np.maximum(pv, 0.0) / scale, shape))


def fit_weibull(pv, retention):
    pv = np.asarray(pv, dtype=float)
    retention = np.asarray(retention, dtype=float)
    ok = np.isfinite(pv) & np.isfinite(retention) & (pv >= 0) & (retention > 0)
    pv, retention = pv[ok], retention[ok]
    if len(pv) < 8 or pv.max() <= 0:
        return {"weibull_scale_pv": math.nan, "weibull_shape": math.nan,
                "weibull_rmse": math.nan, "weibull_r2": math.nan}
    grid = np.linspace(0.0, pv.max(), PV_GRID_POINTS)
    observed = interpolate(np.r_[0.0, pv], np.r_[1.0, retention], grid)
    ok = np.isfinite(observed) & (observed > 0)
    grid, observed = grid[ok], observed[ok]
    try:
        p50 = grid[np.argmin(np.abs(observed - 0.5))]
        p0 = (max(float(p50), 1e-3), 1.0)
        popt, _ = curve_fit(
            weibull_retention,
            grid,
            observed,
            p0=p0,
            bounds=([1e-8, 0.05], [1e6, 10.0]),
            maxfev=20000,
        )
        predicted = weibull_retention(grid, *popt)
        residual = observed - predicted
        sse = float(np.sum(residual ** 2))
        sst = float(np.sum((observed - observed.mean()) ** 2))
        return {
            "weibull_scale_pv": float(popt[0]),
            "weibull_shape": float(popt[1]),
            "weibull_rmse": float(np.sqrt(np.mean(residual ** 2))),
            "weibull_r2": 1.0 - sse / sst if sst > 0 else math.nan,
        }
    except (RuntimeError, ValueError, FloatingPointError):
        return {"weibull_scale_pv": math.nan, "weibull_shape": math.nan,
                "weibull_rmse": math.nan, "weibull_r2": math.nan}


def weighted_quantile(values, weights, quantile):
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[ok], weights[ok]
    if not len(values):
        return math.nan
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    return float(np.interp(quantile * weights.sum(), cumulative, values))


def distribution_metrics(weights):
    weights = np.asarray(weights, dtype=float)
    weights = weights[np.isfinite(weights) & (weights > 0)]
    if not len(weights):
        return {key: math.nan for key in (
            "gas_entropy", "gas_entropy_normalized", "gas_effective_count_entropy",
            "gas_simpson_concentration", "gas_effective_count_simpson", "gas_gini",
            "largest_gas_unit_fraction", "top_1pct_gas_fraction",
            "top_5pct_gas_fraction", "top_10pct_gas_fraction")}
    shares = weights / weights.sum()
    entropy = float(-np.sum(shares * np.log(shares)))
    sorted_w = np.sort(weights)
    n = len(sorted_w)
    gini = float((2 * np.sum(np.arange(1, n + 1) * sorted_w) /
                  (n * sorted_w.sum())) - (n + 1) / n)

    def top_fraction(frac):
        count = max(1, int(math.ceil(n * frac)))
        return float(sorted_w[-count:].sum() / sorted_w.sum())

    simpson = float(np.sum(shares ** 2))
    return {
        "gas_entropy": entropy,
        "gas_entropy_normalized": entropy / math.log(n) if n > 1 else 0.0,
        "gas_effective_count_entropy": math.exp(entropy),
        "gas_simpson_concentration": simpson,
        "gas_effective_count_simpson": 1.0 / simpson,
        "gas_gini": gini,
        "largest_gas_unit_fraction": float(shares.max()),
        "top_1pct_gas_fraction": top_fraction(0.01),
        "top_5pct_gas_fraction": top_fraction(0.05),
        "top_10pct_gas_fraction": top_fraction(0.10),
    }


VTK_DTYPES = {
    "Float64": np.dtype("<f8"), "Float32": np.dtype("<f4"),
    "Int64": np.dtype("<i8"), "UInt64": np.dtype("<u8"),
    "Int32": np.dtype("<i4"), "UInt32": np.dtype("<u4"),
}


def read_vtp_arrays(path: Path, requested):
    """Read selected zlib-compressed appended arrays written by this solver."""
    payload = path.read_bytes()
    marker = b'<AppendedData encoding="raw">_'
    start = payload.index(marker) + len(marker)
    header = payload[:start].decode("utf-8", "ignore")
    arrays = {}
    for name in requested:
        pattern = (
            r'<DataArray type="([^"]+)" Name="' + re.escape(name) +
            r'" NumberOfComponents="(\d+)" format="appended" offset="(\d+)"'
        )
        match = re.search(pattern, header)
        if not match:
            raise ValueError(f"array {name!r} absent from {path}")
        type_name, components, offset = match.groups()
        dtype = VTK_DTYPES[type_name]
        pos = start + int(offset)
        blocks, block_size, last_size = struct.unpack_from("<QQQ", payload, pos)
        sizes = struct.unpack_from("<" + "Q" * blocks, payload, pos + 24)
        compressed = pos + 24 + 8 * blocks
        raw_parts = []
        cursor = compressed
        for size in sizes:
            raw_parts.append(zlib.decompress(payload[cursor:cursor + size]))
            cursor += size
        raw = b"".join(raw_parts)
        expected = (blocks - 1) * block_size + last_size if blocks else 0
        if len(raw) != expected:
            raise ValueError(f"decompressed byte count mismatch for {name} in {path}")
        values = np.frombuffer(raw, dtype=dtype)
        components = int(components)
        arrays[name] = values.reshape((-1, components)) if components > 1 else values
    return arrays


def initial_spatial_metrics(vtp_path: Path, domain_length_m: float):
    arrays = read_vtp_arrays(
        vtp_path,
        ("node_id", "gas_moles", "gas_volume", "active_gas", "segment_id",
         "normalized_x", "Points", "bubble_equivalent_radius"),
    )
    active = ((arrays["active_gas"] > 0.5) & (arrays["segment_id"] >= 1) &
              (arrays["gas_moles"] > 0))
    gas = arrays["gas_moles"][active]
    volume = arrays["gas_volume"][active]
    xyz = arrays["Points"][active]
    xnorm = arrays["normalized_x"][active]
    segments = arrays["segment_id"][active].astype(int)
    node_ids = arrays["node_id"][active].astype(np.int64)
    radii = arrays["bubble_equivalent_radius"][active]
    total = gas.sum()
    center = np.average(xyz, axis=0, weights=gas)
    variance = np.average((xyz - center) ** 2, axis=0, weights=gas)
    result = distribution_metrics(gas)
    result.update({
        "initial_real_active_gas_units": int(len(gas)),
        "initial_gas_moles_vtp": float(total),
        "initial_gas_volume_vtp_m3": float(volume.sum()),
        "initial_bubble_radius_p10_m": float(np.quantile(radii, 0.10)),
        "initial_bubble_radius_p50_m": float(np.quantile(radii, 0.50)),
        "initial_bubble_radius_p90_m": float(np.quantile(radii, 0.90)),
        "gas_center_x_normalized": float(np.average(xnorm, weights=gas)),
        "gas_longitudinal_variance_normalized": float(variance[0] / domain_length_m ** 2),
        "gas_transverse_variance_normalized": float(
            0.5 * (variance[1] + variance[2]) / domain_length_m ** 2),
    })
    return result, {
        "node_id": node_ids, "gas_moles": gas, "gas_volume": volume,
        "segment_id": segments, "radius_m": radii,
    }


def event_and_survival_metrics(events_path: Path, initial_nodes, time, pv, final_time, final_pv):
    events = pd.read_csv(events_path, sep="\t")
    if len(events):
        events["event_pv"] = interpolate(
            np.r_[0.0, time], np.r_[0.0, pv], events["estimated_crossing_time_s"].to_numpy())
        events = events.sort_values("estimated_crossing_time_s").drop_duplicates(
            subset="event_node_id", keep="first")
    event_map = events.set_index("event_node_id").to_dict("index") if len(events) else {}
    rows = []
    for node, moles, volume, segment, radius in zip(
        initial_nodes["node_id"], initial_nodes["gas_moles"],
        initial_nodes["gas_volume"], initial_nodes["segment_id"], initial_nodes["radius_m"]
    ):
        event = event_map.get(int(node))
        rows.append({
            "node_id": int(node), "segment_id": int(segment),
            "initial_gas_moles": float(moles), "initial_gas_volume_m3": float(volume),
            "initial_radius_m": float(radius), "event_observed": int(event is not None),
            "duration_s": float(event["estimated_crossing_time_s"] if event else final_time),
            "duration_pv": float(event["event_pv"] if event else final_pv),
            "event_mode": event["event_mode"] if event else "right_censored",
        })
    survival = pd.DataFrame(rows)
    observed = survival[survival.event_observed == 1]
    initial_weights = survival.initial_gas_moles.to_numpy()
    result = {
        "event_records": int(len(events)),
        "exact_localized_events": int((events.event_mode == "exact_localized").sum()) if len(events) else 0,
        "batched_microbubble_events": int((events.event_mode == "batched_microbubble").sum()) if len(events) else 0,
        "censored_gas_units_at_termination": int((survival.event_observed == 0).sum()),
        "gas_unit_event_fraction": float(survival.event_observed.mean()) if len(survival) else math.nan,
        "event_pv_p50_count_weighted": float(observed.duration_pv.median()) if len(observed) else math.nan,
        "event_pv_p50_initial_mole_weighted": weighted_quantile(
            observed.duration_pv, observed.initial_gas_moles, 0.50) if len(observed) else math.nan,
        "surviving_initial_gas_weight_fraction": float(
            initial_weights[survival.event_observed.to_numpy() == 0].sum() / initial_weights.sum()),
    }
    return result, survival, events


def load_global_history(path: Path):
    columns = [
        "step", "time_s", "injected_pore_volumes", "main_real_gas_saturation",
        "gas_volume_m3", "free_gas_mol", "dissolved_gas_mol", "active_bubble_count",
        "throughflow_m3_s", "gas_liquid_interfacial_area_m2",
    ]
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    use = [column for column in columns if column in header]
    return pd.read_csv(path, sep="\t", usecols=use)


def diagnostics_metrics(path: Path):
    columns = ["epsilon_N_adjusted", "epsilon_V_adjusted", "residual_scaled_inf",
               "relative_update_norm", "timestep_retries"]
    data = pd.read_csv(path, sep="\t", usecols=columns)
    return {
        "max_abs_epsilon_N": float(data.epsilon_N_adjusted.abs().max()),
        "max_abs_epsilon_V": float(data.epsilon_V_adjusted.abs().max()),
        "max_residual_scaled_inf": float(data.residual_scaled_inf.abs().max()),
        "max_relative_update_norm": float(data.relative_update_norm.abs().max()),
        "total_timestep_retries": int(data.timestep_retries.sum()),
    }


def geometry_metrics(path: Path):
    data = pd.read_csv(path, sep="\t")
    pore_w = data.real_pore_count.to_numpy(float)
    throat_w = np.maximum(data.internal_throat_count.to_numpy(float), 1.0)

    def weighted(column, weights):
        return float(np.average(data[column].to_numpy(float), weights=weights))

    def segment_cv(column):
        values = data[column].to_numpy(float)
        mean = values.mean()
        return float(values.std() / mean) if mean != 0 else math.nan

    return {
        "pore_radius_mean_m": weighted("pore_radius_mean_m", pore_w),
        "pore_radius_p50_segment_weighted_m": weighted("pore_radius_p50_m", pore_w),
        "throat_radius_mean_m": weighted("throat_radius_mean_m", throat_w),
        "throat_radius_p50_segment_weighted_m": weighted("throat_radius_p50_m", throat_w),
        "mean_coordination_number": weighted("mean_coordination_number", pore_w),
        "constriction_ratio_p50_segment_weighted": weighted("constriction_ratio_p50", throat_w),
        "hydraulic_conductance_sum_m3_Pa_s": float(data.hydraulic_conductance_sum_m3_Pa_s.sum()),
        "cross_segment_throat_count": int(data.cross_segment_throat_count.sum()),
        "segment_pore_radius_cv": segment_cv("pore_radius_mean_m"),
        "segment_throat_radius_cv": segment_cv("throat_radius_mean_m"),
        "segment_hydraulic_conductance_cv": segment_cv("hydraulic_conductance_sum_m3_Pa_s"),
    }


def load_selected_segment_states(path: Path, target_steps):
    columns = ["step", "segment_id", "segment_mean_Sg", "segment_gas_volume_m3",
               "segment_free_gas_mol", "segment_active_bubble_count"]
    selected = []
    target_steps = set(int(x) for x in target_steps)
    for chunk in pd.read_csv(path, sep="\t", usecols=columns, chunksize=250000):
        part = chunk[chunk.step.isin(target_steps)]
        if len(part):
            selected.append(part)
    if not selected:
        return pd.DataFrame(columns=columns)
    return pd.concat(selected, ignore_index=True)


def segment_spatial_moments(gas_moles, gas_volume, active_units):
    gas_moles = np.asarray(gas_moles, dtype=float)
    gas_volume = np.asarray(gas_volume, dtype=float)
    active_units = np.asarray(active_units, dtype=float)
    x = (np.arange(20, dtype=float) + 0.5) / 20.0
    total = float(np.sum(gas_moles))
    if total > 0:
        shares = gas_moles / total
        positive = shares > 0
        entropy = float(-np.sum(shares[positive] * np.log(shares[positive])) / math.log(20))
        center = float(np.sum(shares * x))
        variance = float(np.sum(shares * (x - center) ** 2))
    else:
        entropy = center = variance = math.nan
    mean = float(np.mean(gas_moles))
    return {
        "free_gas_moles_zeroth_moment": total,
        "gas_volume_zeroth_moment_m3": float(np.sum(gas_volume)),
        "active_gas_units": int(np.sum(active_units)),
        "segment_gas_spatial_entropy_normalized": entropy,
        "gas_center_x_normalized_20segment": center,
        "gas_longitudinal_variance_normalized_20segment": variance,
        "segment_gas_moles_cv": float(np.std(gas_moles) / mean) if mean > 0 else math.nan,
    }


def bootstrap_ci(values, statistic=np.median, seed=20260915, repetitions=10000):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return math.nan, math.nan
    rng = np.random.default_rng(seed)
    estimates = np.empty(repetitions)
    for i in range(repetitions):
        estimates[i] = statistic(rng.choice(values, size=len(values), replace=True))
    return tuple(float(x) for x in np.quantile(estimates, [0.025, 0.975]))


def paired_rank_biserial(differences):
    d = np.asarray(differences, dtype=float)
    d = d[np.isfinite(d) & (d != 0)]
    if not len(d):
        return math.nan
    ranks = rankdata(np.abs(d))
    return float((ranks[d > 0].sum() - ranks[d < 0].sum()) / ranks.sum())


def correlation_rows(frame, predictors, outcomes, analysis_unit):
    rows = []
    for predictor in predictors:
        for outcome in outcomes:
            subset = frame[[predictor, outcome]].replace([np.inf, -np.inf], np.nan).dropna()
            if len(subset) < 5:
                continue
            estimate, pvalue = spearmanr(subset[predictor], subset[outcome])
            rng = np.random.default_rng(20260915)
            boot = []
            values = subset.to_numpy()
            for _ in range(5000):
                sampled = values[rng.integers(0, len(values), len(values))]
                if np.unique(sampled[:, 0]).size < 2 or np.unique(sampled[:, 1]).size < 2:
                    continue
                boot.append(spearmanr(sampled[:, 0], sampled[:, 1]).statistic)
            low, high = np.quantile(boot, [0.025, 0.975]) if boot else (math.nan, math.nan)
            rows.append({
                "analysis_unit": analysis_unit, "predictor": predictor, "outcome": outcome,
                "n": len(subset), "spearman_rho": estimate, "p_value_unadjusted": pvalue,
                "bootstrap_ci95_low": low, "bootstrap_ci95_high": high,
            })
    result = pd.DataFrame(rows)
    if len(result):
        # Benjamini-Hochberg adjustment within this planned family.
        order = np.argsort(result.p_value_unadjusted.to_numpy())
        p = result.p_value_unadjusted.to_numpy()[order]
        adjusted = np.minimum.accumulate((p * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
        restored = np.empty_like(adjusted)
        restored[order] = np.minimum(adjusted, 1.0)
        result["p_value_fdr_bh"] = restored
    return result


def write_tsv(frame: pd.DataFrame, path: Path):
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.12g")


def configure_plots():
    cjk_font = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    if cjk_font.is_file():
        font_manager.fontManager.addfont(cjk_font)
        family = font_manager.FontProperties(fname=cjk_font).get_name()
    else:
        family = "DejaVu Sans"
    plt.rcParams.update({
        "font.family": family,
        "axes.unicode_minus": False, "figure.dpi": 140, "savefig.dpi": 220,
        "axes.spines.top": False, "axes.spines.right": False,
    })


def make_plots(case, curves, paired, correlations, segment, spatial_moments, figures: Path):
    configure_plots()
    figures.mkdir(parents=True, exist_ok=True)
    colors = {"PNM-0.6": "#2166ac", "PNM-0.8": "#b2182b"}

    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    for variant, group in curves.groupby("variant"):
        pivot = group.pivot(index="pv", columns="run_id", values="free_gas_retained_fraction")
        ax.plot(pivot.index, pivot, color=colors[variant], alpha=0.09, lw=0.7)
        median = pivot.median(axis=1)
        q1, q3 = pivot.quantile(0.25, axis=1), pivot.quantile(0.75, axis=1)
        ax.plot(pivot.index, median, color=colors[variant], lw=2.5, label=f"{variant} median / 中位数")
        ax.fill_between(pivot.index, q1, q3, color=colors[variant], alpha=0.18)
    ax.set(xlabel="Injected pore volumes / 注入孔隙体积 PV",
           ylabel="Free-gas retention $n_g/n_{g,0}$ / 自由气体保留率",
           xlim=(0, curves.pv.max()), ylim=(0, 1.03))
    ax.legend(frameon=False); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(figures / "01_free_gas_retention_common_pv.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    pv50_name = "pv_to_50pct_free_gas_loss"
    order = paired.sort_values(f"{pv50_name}_p06").parent_case.to_list()
    y = np.arange(len(order)); lookup = paired.set_index("parent_case").loc[order]
    for j, (_, row) in enumerate(lookup.iterrows()):
        ax.plot([row[f"{pv50_name}_p06"], row[f"{pv50_name}_p08"]],
                [j, j], color="#bdbdbd", lw=1)
    ax.scatter(lookup[f"{pv50_name}_p06"], y, color=colors["PNM-0.6"], s=24, label="PNM-0.6")
    ax.scatter(lookup[f"{pv50_name}_p08"], y, color=colors["PNM-0.8"], s=24, label="PNM-0.8")
    ax.set(xlabel="$PV_{50}$ / 50%自由气体损失所需PV", ylabel="Paired geometry / 配对几何")
    ax.set_yticks([]); ax.grid(axis="x", alpha=.2); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(figures / "02_paired_pv50.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    for variant, marker in [("PNM-0.6", "o"), ("PNM-0.8", "s")]:
        part = case[case.variant == variant]
        ax.scatter(part.control_volume_porosity, part[pv50_name], color=colors[variant], marker=marker,
                   alpha=.78, label=variant)
    ax.set(xlabel="PNM porosity / 孔隙率", ylabel="$PV_{50}$")
    ax.grid(alpha=.2); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(figures / "03_porosity_vs_pv50.png"); plt.close(fig)

    if len(correlations):
        primary = correlations[correlations.analysis_unit == "complete-pair geometry"]
        table = primary.pivot(index="predictor", columns="outcome", values="spearman_rho")
        fig, ax = plt.subplots(figsize=(7.8, 5.3))
        image = ax.imshow(table, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(len(table.columns)), table.columns, rotation=30, ha="right")
        ax.set_yticks(range(len(table.index)), table.index)
        for i in range(table.shape[0]):
            for j in range(table.shape[1]):
                ax.text(j, i, f"{table.iloc[i,j]:.2f}", ha="center", va="center", fontsize=8)
        ax.set_title("Spearman association / Spearman关联（34个完整配对几何）")
        fig.colorbar(image, ax=ax, label="Spearman $\\rho$")
        fig.tight_layout(); fig.savefig(figures / "04_geometry_durability_correlations.png"); plt.close(fig)

    if len(segment):
        part = segment[segment.target_pv == 5.0]
        table = part.groupby(["variant", "segment_id"]).event_fraction_of_initial_units.mean().unstack(0)
        fig, ax = plt.subplots(figsize=(7.3, 4.4))
        for variant in table.columns:
            ax.plot(table.index, table[variant], marker="o", ms=3, color=colors[variant], label=variant)
        ax.set(xlabel="Flow-direction segment / 沿流向区段",
               ylabel="Retired fraction by 5 PV / 5 PV前消失比例")
        ax.set_xticks(range(1, 21)); ax.grid(alpha=.2); ax.legend(frameon=False)
        fig.tight_layout(); fig.savefig(figures / "05_segment_event_fraction_at_5pv.png"); plt.close(fig)

    if len(spatial_moments):
        common = spatial_moments[spatial_moments.target_pv <= 5.0]
        metrics = [
            ("segment_gas_spatial_entropy_normalized", "Spatial entropy / 空间熵"),
            ("gas_center_x_normalized_20segment", "Gas centroid $x_c/L$ / 气体重心"),
            ("gas_longitudinal_variance_normalized_20segment", "Longitudinal variance / 纵向方差"),
        ]
        fig, axes = plt.subplots(1, 3, figsize=(12.8, 3.9))
        for ax, (metric, label) in zip(axes, metrics):
            for variant in ("PNM-0.6", "PNM-0.8"):
                part = common[common.variant == variant]
                grouped = part.groupby("target_pv")[metric]
                median, q1, q3 = grouped.median(), grouped.quantile(.25), grouped.quantile(.75)
                ax.plot(median.index, median, marker="o", color=colors[variant], label=variant)
                ax.fill_between(median.index, q1, q3, color=colors[variant], alpha=.18)
            ax.set(xlabel="Injected PV / 注入PV", ylabel=label)
            ax.grid(alpha=.2)
        axes[0].legend(frameon=False)
        fig.suptitle("Anna-style 20-segment gas spatial moments / 20段气体空间矩")
        fig.tight_layout(); fig.savefig(figures / "06_gas_spatial_moments_0_5pv.png"); plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--index", type=Path)
    parser.add_argument("--inventory", type=Path)
    parser.add_argument("--source-audit", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--common-pv-horizon", type=float, default=5.0)
    args = parser.parse_args()
    project = args.project_root.resolve()
    index_path = (args.index or project / "configs/batch35_flow_pe555_pardiso_v3/batch_cases.json").resolve()
    inventory_path = (args.inventory or project / "output/batch35_inputs/batch_input_inventory.json").resolve()
    audit_path = (args.source_audit or project / "output/verification/batch35_saturation_source_audit.json").resolve()
    output = (args.output or project / "output/analysis/batch69_bubble_durability").resolve()
    if output.exists() and any(output.iterdir()):
        raise SystemExit(f"refusing to overwrite non-empty analysis directory: {output}")
    output.mkdir(parents=True, exist_ok=True)

    index = json.loads(index_path.read_text())
    inventory = json.loads(inventory_path.read_text())
    source_audit = json.loads(audit_path.read_text())
    failed = {(x["parent_case"], x["variant"]) for x in source_audit["cases"] if x["status"] != "passed"}
    inventory_map = {(x["parent_case"], x["variant"]): x for x in inventory["cases"]}
    domain_length_m = float(inventory["effective_physical_domain_um"][0]) * 1e-6

    case_rows, curve_rows, all_survival, segment_rows, spatial_moment_rows = [], [], [], [], []
    input_records, exclusions = [], []
    for number, entry in enumerate(index["cases"], start=1):
        key = (entry["parent_case"], entry["variant"])
        if key in failed:
            exclusions.append({"run_id": entry["run_id"], "parent_case": key[0],
                               "variant": key[1], "reason": "failed saturation-source audit"})
            continue
        run = Path(entry["output_directory"]).resolve()
        required = [run / name for name in (
            "global_history.tsv", "diagnostics.tsv", "termination_events.json",
            "gas_disappearance_events.tsv", "segment_geometry.tsv", "state_000000.vtp",
            "run_manifest.json", "segment_history.tsv")]
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            raise FileNotFoundError(f"incomplete valid run {entry['run_id']}: {missing}")
        history = load_global_history(required[0])
        term = json.loads(required[2].read_text())
        manifest = json.loads(required[6].read_text())
        if term.get("termination_reason") != "main_real_gas_saturation_reached":
            raise ValueError(f"unexpected termination for {entry['run_id']}: {term.get('termination_reason')}")
        identity = parse_case_identity(entry["parent_case"])
        inv = inventory_map[key]
        time = history.time_s.to_numpy(float)
        pv = history.injected_pore_volumes.to_numpy(float)
        free = history.free_gas_mol.to_numpy(float)
        volume = history.gas_volume_m3.to_numpy(float)
        initial_free = float(term["initial_total_free_gas_moles"])
        retained = free / initial_free
        spatial, initial_nodes = initial_spatial_metrics(required[5], domain_length_m)
        initial_volume = spatial["initial_gas_volume_vtp_m3"]
        volume_retained = volume / initial_volume
        event_metrics, survival, events = event_and_survival_metrics(
            required[3], initial_nodes, time, pv, float(term["actual_final_time_s"]), float(pv[-1]))
        survival.insert(0, "run_id", entry["run_id"])
        survival.insert(1, "parent_case", entry["parent_case"])
        survival.insert(2, "variant", entry["variant"])
        all_survival.append(survival)

        event_lookup = {(event.get("kind"), float(event.get("target", -1))): event
                        for event in term.get("events", [])}
        row = {
            "run_id": entry["run_id"], "parent_case": entry["parent_case"],
            "variant": entry["variant"], **identity,
            "control_volume_m3": inv["control_volume_m3"],
            "control_volume_porosity": inv["control_volume_porosity"],
            "raw_pore_count": inv["raw_pore_count"], "raw_throat_count": inv["raw_throat_count"],
            "derived_node_count": inv["derived_node_count"], "derived_edge_count": inv["derived_edge_count"],
            "estimated_pore_reynolds_p50": inv["estimated_pore_reynolds_p50"],
            "estimated_capillary_number": inv["estimated_capillary_number"],
            "achieved_sample_peclet": inv["achieved_sample_peclet"],
            "initial_main_Sg": float(history.iloc[0].main_real_gas_saturation),
            "initial_free_gas_moles": initial_free,
            "initial_gas_volume_m3": spatial["initial_gas_volume_vtp_m3"],
            "initial_active_gas_units_global": int(history.iloc[0].active_bubble_count),
            "final_time_s": float(term["actual_final_time_s"]),
            "final_injected_pore_volumes": float(pv[-1]),
            "final_main_Sg": float(term["actual_final_main_real_gas_saturation"]),
            "final_free_gas_moles": float(term["current_total_free_gas_moles"]),
            "free_gas_retained_fraction_final": float(free[-1] / initial_free),
            "auc_free_gas_retention_5pv": normalized_auc(pv, retained, args.common_pv_horizon),
            "auc_gas_volume_retention_5pv": normalized_auc(pv, volume_retained, args.common_pv_horizon),
            "free_gas_retention_at_1pv": interpolate(np.r_[0, pv], np.r_[1, retained], [1])[0],
            "free_gas_retention_at_5pv": interpolate(np.r_[0, pv], np.r_[1, retained], [5])[0],
            "termination_reason": term["termination_reason"],
            "accepted_steps": int(history.iloc[-1].step),
            "actual_final_throughflow_m3_s": float(history.iloc[-1].throughflow_m3_s),
            **geometry_metrics(required[4]), **spatial, **event_metrics,
            **diagnostics_metrics(required[1]),
            **fit_weibull(pv, retained),
        }
        for loss in LOSS_LEVELS:
            event = event_lookup.get(("free_gas_loss_fraction", loss))
            label = int(round(loss * 100))
            row[f"time_to_{label}pct_free_gas_loss_s"] = finite(event and event.get("estimated_crossing_time_s"))
            row[f"pv_to_{label}pct_free_gas_loss"] = finite(event and event.get("estimated_crossing_injected_pore_volumes"))
            row[f"{label}pct_loss_right_censored"] = int(event is None)
        case_rows.append(row)

        common_grid = np.linspace(0.0, args.common_pv_horizon, PV_GRID_POINTS)
        rn = interpolate(np.r_[0, pv], np.r_[1, retained], common_grid)
        rv = interpolate(np.r_[0, pv], np.r_[1, volume_retained], common_grid)
        for x, a, b in zip(common_grid, rn, rv):
            curve_rows.append({"run_id": entry["run_id"], "parent_case": entry["parent_case"],
                               "variant": entry["variant"], "pv": x,
                               "free_gas_retained_fraction": a, "gas_volume_retained_fraction": b})

        target_step = {}
        for target in SEGMENT_PV_TARGETS:
            if target > 0 and target <= pv[-1]:
                target_step[target] = int(history.iloc[int(np.argmin(np.abs(pv - target)))].step)
        selected_segment = load_selected_segment_states(required[7], target_step.values())
        initial_segment = pd.DataFrame({
            "segment_id": np.arange(1, 21),
            "segment_free_gas_mol": [initial_nodes["gas_moles"][initial_nodes["segment_id"] == s].sum()
                                     for s in range(1, 21)],
            "segment_gas_volume_m3": [initial_nodes["gas_volume"][initial_nodes["segment_id"] == s].sum()
                                      for s in range(1, 21)],
            "segment_active_bubble_count": [(initial_nodes["segment_id"] == s).sum()
                                             for s in range(1, 21)],
        })
        initial_moment = segment_spatial_moments(
            initial_segment.segment_free_gas_mol, initial_segment.segment_gas_volume_m3,
            initial_segment.segment_active_bubble_count)
        spatial_moment_rows.append({
            "run_id": entry["run_id"], "parent_case": entry["parent_case"],
            "variant": entry["variant"], "gradation": identity["gradation"],
            "target_pv": 0.0, "actual_pv": 0.0, "source_step": 0, **initial_moment,
        })
        for target, step in target_step.items():
            state = selected_segment[selected_segment.step == step].sort_values("segment_id")
            if len(state) != 20:
                raise ValueError(f"expected 20 segments for {entry['run_id']} step {step}, got {len(state)}")
            moment = segment_spatial_moments(
                state.segment_free_gas_mol, state.segment_gas_volume_m3,
                state.segment_active_bubble_count)
            history_row = history[history.step == step].iloc[-1]
            spatial_moment_rows.append({
                "run_id": entry["run_id"], "parent_case": entry["parent_case"],
                "variant": entry["variant"], "gradation": identity["gradation"],
                "target_pv": target, "actual_pv": float(history_row.injected_pore_volumes),
                "source_step": step, **moment,
            })

        # Segment-level initial population and disappearance fraction at selected PVs.
        initial_by_segment = survival.groupby("segment_id").agg(
            initial_gas_units=("node_id", "size"), initial_gas_moles=("initial_gas_moles", "sum"),
            initial_gas_volume_m3=("initial_gas_volume_m3", "sum"),
        )
        observed = survival[survival.event_observed == 1]
        for target in SEGMENT_PV_TARGETS:
            if target > pv[-1]:
                continue
            happened = observed[observed.duration_pv <= target].groupby("segment_id").agg(
                disappeared_units=("node_id", "size"),
                disappeared_initial_moles=("initial_gas_moles", "sum"),
            )
            for segment_id in range(1, 21):
                initial = initial_by_segment.loc[segment_id] if segment_id in initial_by_segment.index else None
                gone = happened.loc[segment_id] if segment_id in happened.index else None
                n0 = int(initial.initial_gas_units) if initial is not None else 0
                m0 = float(initial.initial_gas_moles) if initial is not None else 0.0
                ng = int(gone.disappeared_units) if gone is not None else 0
                mg = float(gone.disappeared_initial_moles) if gone is not None else 0.0
                segment_rows.append({
                    "run_id": entry["run_id"], "parent_case": entry["parent_case"],
                    "variant": entry["variant"], "gradation": identity["gradation"],
                    "target_pv": target, "segment_id": segment_id,
                    "initial_gas_units": n0, "initial_gas_moles": m0,
                    "disappeared_units": ng, "disappeared_initial_moles": mg,
                    "event_fraction_of_initial_units": ng / n0 if n0 else math.nan,
                    "event_fraction_of_initial_moles": mg / m0 if m0 else math.nan,
                })

        input_records.append({
            "run_id": entry["run_id"],
            "files": {path.name: {"path": str(path), "sha256": sha256(path)} for path in required},
        })
        print(f"[{number:02d}/{len(index['cases'])}] {entry['run_id']}", flush=True)

    case = pd.DataFrame(case_rows).sort_values(["parent_case", "variant"])
    curves = pd.DataFrame(curve_rows)
    segment = pd.DataFrame(segment_rows)
    spatial_moments = pd.DataFrame(spatial_moment_rows)
    survival = pd.concat(all_survival, ignore_index=True)
    if len(case) != 69:
        raise ValueError(f"expected 69 valid cases, obtained {len(case)}")
    if case.parent_case.nunique() != 35:
        raise ValueError("expected 35 unique geometries")

    # Complete-pair comparison; the invalid 85A15D/mu=.5/100 kPa PNM-0.8 is not imputed.
    pair_wide = case.pivot(index="parent_case", columns="variant")
    complete_parents = [p for p in pair_wide.index
                        if pd.notna(pair_wide.loc[p, ("pv_to_50pct_free_gas_loss", "PNM-0.6")])
                        and pd.notna(pair_wide.loc[p, ("pv_to_50pct_free_gas_loss", "PNM-0.8")])]
    pair_metrics = [
        "initial_main_Sg", "initial_free_gas_moles", "initial_real_active_gas_units",
        "gas_entropy_normalized", "gas_gini", "gas_center_x_normalized",
        "pv_to_50pct_free_gas_loss", "auc_free_gas_retention_5pv",
        "final_injected_pore_volumes", "weibull_scale_pv", "weibull_shape",
    ]
    paired_rows = []
    for parent in complete_parents:
        base = {"parent_case": parent}
        for metric in pair_metrics:
            a = float(pair_wide.loc[parent, (metric, "PNM-0.6")])
            b = float(pair_wide.loc[parent, (metric, "PNM-0.8")])
            base[f"{metric}_p06"] = a; base[f"{metric}_p08"] = b
            base[f"{metric}_difference_p08_minus_p06"] = b - a
            base[f"{metric}_relative_difference"] = b / a - 1 if a != 0 else math.nan
        paired_rows.append(base)
    paired = pd.DataFrame(paired_rows)

    paired_stats = []
    for metric in pair_metrics:
        diff = paired[f"{metric}_difference_p08_minus_p06"].to_numpy(float)
        relative = paired[f"{metric}_relative_difference"].to_numpy(float)
        ci_low, ci_high = bootstrap_ci(diff)
        try:
            test = wilcoxon(diff, alternative="two-sided", zero_method="wilcox")
            statistic, pvalue = float(test.statistic), float(test.pvalue)
        except ValueError:
            statistic, pvalue = math.nan, math.nan
        paired_stats.append({
            "metric": metric, "n_pairs": len(diff), "median_p06": np.nanmedian(paired[f"{metric}_p06"]),
            "median_p08": np.nanmedian(paired[f"{metric}_p08"]), "median_difference_p08_minus_p06": np.nanmedian(diff),
            "median_relative_difference": np.nanmedian(relative), "bootstrap_median_difference_ci95_low": ci_low,
            "bootstrap_median_difference_ci95_high": ci_high, "wilcoxon_statistic": statistic,
            "wilcoxon_p_value_unadjusted": pvalue, "paired_rank_biserial": paired_rank_biserial(diff),
        })
    paired_stats = pd.DataFrame(paired_stats)

    # One row per complete-pair geometry prevents pseudo-replication.
    structure = [
        "control_volume_porosity", "raw_pore_count", "raw_throat_count",
        "estimated_pore_reynolds_p50", "estimated_capillary_number",
        "pore_radius_mean_m", "throat_radius_mean_m", "mean_coordination_number",
        "constriction_ratio_p50_segment_weighted", "hydraulic_conductance_sum_m3_Pa_s",
        "segment_pore_radius_cv", "segment_throat_radius_cv",
        "segment_hydraulic_conductance_cv",
    ]
    outcome = [
        "pv_to_50pct_free_gas_loss", "auc_free_gas_retention_5pv",
        "final_injected_pore_volumes", "weibull_scale_pv", "weibull_shape",
    ]
    complete_case = case[case.parent_case.isin(complete_parents)]
    geometry = complete_case.groupby("parent_case", as_index=False)[structure + outcome].mean()
    correlations = correlation_rows(geometry, structure, outcome, "complete-pair geometry")
    morphology_predictors = [
        "initial_main_Sg", "initial_free_gas_moles", "initial_real_active_gas_units",
        "gas_entropy_normalized", "gas_gini", "largest_gas_unit_fraction",
        "top_10pct_gas_fraction", "gas_center_x_normalized",
        "gas_longitudinal_variance_normalized", "gas_transverse_variance_normalized",
        "initial_bubble_radius_p50_m", "initial_bubble_radius_p90_m",
    ]
    morphology_outcomes = ["pv_to_50pct_free_gas_loss", "auc_free_gas_retention_5pv",
                           "weibull_scale_pv", "weibull_shape"]
    morphology_correlations = pd.concat([
        correlation_rows(case[case.variant == variant], morphology_predictors,
                         morphology_outcomes, f"case within {variant}")
        for variant in ("PNM-0.6", "PNM-0.8")
    ], ignore_index=True)

    # Planned blocked tests: friction at 100 kPa; stress at mu=0.5.
    factor_rows, factor_tests = [], []
    for block_name, factor, levels, mask in (
        ("friction_at_100kpa", "mu", [0.1, 0.3, 0.5], case.stress_kpa == 100),
        ("stress_at_mu_0p5", "stress_kpa", [20.0, 100.0, 200.0], case.mu == 0.5),
    ):
        for variant in ("PNM-0.6", "PNM-0.8"):
            part = case[mask & (case.variant == variant)]
            for metric in ("pv_to_50pct_free_gas_loss", "auc_free_gas_retention_5pv",
                           "final_injected_pore_volumes"):
                pivot = part.pivot(index="gradation", columns=factor, values=metric).dropna()
                arrays = [pivot[level].to_numpy() for level in levels]
                test = friedmanchisquare(*arrays) if len(pivot) >= 3 else None
                factor_tests.append({
                    "block": block_name, "variant": variant, "metric": metric,
                    "complete_gradation_blocks": len(pivot),
                    "friedman_chi_square": finite(test and test.statistic),
                    "friedman_p_value": finite(test and test.pvalue),
                })
                for level in levels:
                    values = part.loc[part[factor] == level, metric].dropna()
                    factor_rows.append({
                        "block": block_name, "variant": variant, "metric": metric,
                        "factor": factor, "level": level, "n": len(values),
                        "median": values.median(), "q25": values.quantile(.25),
                        "q75": values.quantile(.75), "mean": values.mean(),
                    })
    gradations = sorted(case.gradation.unique())
    for variant in ("PNM-0.6", "PNM-0.8"):
        part = case[case.variant == variant].copy()
        part["condition"] = part.mu.astype(str) + "_" + part.stress_kpa.astype(str)
        for metric in ("pv_to_50pct_free_gas_loss", "auc_free_gas_retention_5pv",
                       "final_injected_pore_volumes"):
            pivot = part.pivot(index="condition", columns="gradation", values=metric).dropna()
            arrays = [pivot[level].to_numpy() for level in gradations]
            test = friedmanchisquare(*arrays) if len(pivot) >= 3 else None
            factor_tests.append({
                "block": "gradation_across_matched_conditions", "variant": variant,
                "metric": metric, "complete_gradation_blocks": len(pivot),
                "friedman_chi_square": finite(test and test.statistic),
                "friedman_p_value": finite(test and test.pvalue),
            })
            for level in gradations:
                values = part.loc[part.gradation == level, metric].dropna()
                factor_rows.append({
                    "block": "gradation_across_matched_conditions", "variant": variant,
                    "metric": metric, "factor": "gradation", "level": level,
                    "n": len(values), "median": values.median(), "q25": values.quantile(.25),
                    "q75": values.quantile(.75), "mean": values.mean(),
                })
    factor_summary, factor_tests = pd.DataFrame(factor_rows), pd.DataFrame(factor_tests)
    if len(factor_tests):
        order = np.argsort(factor_tests.friedman_p_value.to_numpy())
        pvalues = factor_tests.friedman_p_value.to_numpy()[order]
        adjusted = np.minimum.accumulate(
            (pvalues * len(pvalues) / np.arange(1, len(pvalues) + 1))[::-1])[::-1]
        restored = np.empty_like(adjusted); restored[order] = np.minimum(adjusted, 1.0)
        factor_tests["friedman_p_value_fdr_bh"] = restored

    write_tsv(case, output / "case_metrics.tsv")
    write_tsv(curves, output / "retention_curves_common_0_5pv.tsv")
    write_tsv(paired, output / "paired_p06_p08.tsv")
    write_tsv(paired_stats, output / "paired_effects.tsv")
    write_tsv(geometry, output / "complete_pair_geometry_metrics.tsv")
    write_tsv(correlations, output / "geometry_durability_correlations.tsv")
    write_tsv(morphology_correlations, output / "initial_gas_morphology_correlations.tsv")
    write_tsv(factor_summary, output / "factor_level_summaries.tsv")
    write_tsv(factor_tests, output / "blocked_factor_tests.tsv")
    write_tsv(segment, output / "segment_event_metrics.tsv")
    write_tsv(spatial_moments, output / "gas_spatial_moments_20segments.tsv")
    with gzip.open(output / "gas_unit_survival.tsv.gz", "wt", newline="") as stream:
        survival.to_csv(stream, sep="\t", index=False, float_format="%.12g")
    make_plots(case, curves, paired, correlations, segment, spatial_moments, output / "figures")

    pv50 = case.pv_to_50pct_free_gas_loss
    pair_pv50 = paired_stats.set_index("metric").loc["pv_to_50pct_free_gas_loss"]
    geometry_key = correlations[
        (correlations.predictor == "control_volume_porosity") &
        (correlations.outcome == "pv_to_50pct_free_gas_loss")]
    morphology_keys = {}
    for variant in ("PNM-0.6", "PNM-0.8"):
        unit = f"case within {variant}"
        for predictor, outcome, label in (
            ("gas_entropy_normalized", "weibull_shape", "entropy_shape"),
            ("gas_longitudinal_variance_normalized", "pv_to_50pct_free_gas_loss", "spread_pv50"),
        ):
            morphology_keys[(variant, label)] = morphology_correlations[
                (morphology_correlations.analysis_unit == unit) &
                (morphology_correlations.predictor == predictor) &
                (morphology_correlations.outcome == outcome)].iloc[0]
    spatial_summary = spatial_moments[spatial_moments.target_pv.isin([0.0, 5.0])].groupby(
        ["variant", "target_pv"])[
            ["segment_gas_spatial_entropy_normalized", "gas_center_x_normalized_20segment",
             "gas_longitudinal_variance_normalized_20segment"]].median()
    event_counts = Counter(survival.event_mode)
    max_eps_n, max_eps_v = case.max_abs_epsilon_N.max(), case.max_abs_epsilon_V.max()
    report = f"""# 69个有效算例的气泡耐久性定量分析（第一版）

## 数据范围与质量控制

- 纳入：69个有效运行，35套几何；PNM-0.6为35个，PNM-0.8为34个。
- 排除：`sands-of-85A15D-mu-0.5-100kpa-e_p08`，原因是饱和度源审计失败；未插补。
- 34套几何具有完整的PNM-0.6/0.8配对。
- 69个运行均以`main_real_gas_saturation_reached`正常结束。
- 最大`|epsilon_N_adjusted|`={max_eps_n:.6g}；最大`|epsilon_V_adjusted|`={max_eps_v:.6g}。
- 公共曲线与AUC窗口为0–{args.common_pv_horizon:g} PV，所有有效算例均覆盖该范围。

## 核心定量结果

- `PV50`范围：{pv50.min():.4g}–{pv50.max():.4g} PV；中位数{pv50.median():.4g} PV。
- PNM-0.6的`PV50`中位数：{pair_pv50.median_p06:.4g} PV。
- PNM-0.8的`PV50`中位数：{pair_pv50.median_p08:.4g} PV。
- 34个配对中，PNM-0.8相对PNM-0.6的`PV50`中位相对变化：{100*pair_pv50.median_relative_difference:.3g}%。
- 配对差值的95% bootstrap区间：[{pair_pv50.bootstrap_median_difference_ci95_low:.4g}, {pair_pv50.bootstrap_median_difference_ci95_high:.4g}] PV；Wilcoxon p={pair_pv50.wilcoxon_p_value_unadjusted:.4g}。
- 完整配对几何层面，孔隙率与平均`PV50`的Spearman rho={geometry_key.iloc[0].spearman_rho:.4g}，
  95% geometry-bootstrap区间=[{geometry_key.iloc[0].bootstrap_ci95_low:.4g}, {geometry_key.iloc[0].bootstrap_ci95_high:.4g}]。
- 初始气体摩尔分布熵与Weibull形状参数的关联：PNM-0.6 rho={morphology_keys[("PNM-0.6", "entropy_shape")].spearman_rho:.4g}，
  PNM-0.8 rho={morphology_keys[("PNM-0.8", "entropy_shape")].spearman_rho:.4g}。
- 初始气体纵向分散度与`PV50`的关联：PNM-0.6 rho={morphology_keys[("PNM-0.6", "spread_pv50")].spearman_rho:.4g}，
  PNM-0.8 rho={morphology_keys[("PNM-0.8", "spread_pv50")].spearman_rho:.4g}。这些是未控制孔隙率/级配的探索性关联。
- 20段气体空间熵中位数从0到5 PV：PNM-0.6为
  {spatial_summary.loc[("PNM-0.6",0.0),"segment_gas_spatial_entropy_normalized"]:.4g}→{spatial_summary.loc[("PNM-0.6",5.0),"segment_gas_spatial_entropy_normalized"]:.4g}；
  PNM-0.8为{spatial_summary.loc[("PNM-0.8",0.0),"segment_gas_spatial_entropy_normalized"]:.4g}→{spatial_summary.loc[("PNM-0.8",5.0),"segment_gas_spatial_entropy_normalized"]:.4g}。
- 初始活动含气控制体生存表共有{len(survival):,}条；观察到{event_counts['batched_microbubble']:,}个批量退休事件、
  {event_counts['exact_localized']:,}个精确定位事件，其余{event_counts['right_censored']:,}个在终止时仍存活。

## 解释限制

1. `active_bubble_count`是活动含气控制体数，不等同于经三维连通聚类识别的真实气泡数。
2. 批量退休事件的穿越时间是接受步内估计值，不是重新求解的精确事件时刻。
3. 因此自由气体摩尔保留曲线、`PV50`和公共窗口AUC是主指标；数量生存曲线是敏感性指标。
4. 文件名中的`e`后缀不完整且未经统一审计，主分析使用PNM实际孔隙率，不使用文件名`e`作为孔隙比。
5. μ效应只在100 kPa子设计中检验；压力效应只在μ=0.5子设计中检验，不外推不存在的组合。
6. 当前是关联分析，不足以单独证明因果关系。Henry系数、kL、界面面积、毛细压力和相对导流闭合仍未实验标定。

## 文件说明

- `case_metrics.tsv`：69行算例级指标。
- `retention_curves_common_0_5pv.tsv`：公共PV网格上的气体保留曲线。
- `paired_p06_p08.tsv`、`paired_effects.tsv`：34组配对及其统计检验。
- `complete_pair_geometry_metrics.tsv`、`geometry_durability_correlations.tsv`：避免伪重复的几何层分析。
- `initial_gas_morphology_correlations.tsv`：分别在PNM-0.6/0.8内部计算的初始困气形态关联。
- `factor_level_summaries.tsv`、`blocked_factor_tests.tsv`：分块工况分析。
- `segment_event_metrics.tsv`：沿流向20段的消失比例。
- `gas_spatial_moments_20segments.tsv`：Anna式20段零阶矩、重心、纵向方差和空间熵。
- `gas_unit_survival.tsv.gz`：逐含气控制体生存/删失记录。
- `figures/`：论文图初稿。
- `analysis_manifest.json`：输入、排除项、方法和SHA-256来源。
"""
    (output / "ANALYSIS_REPORT_ZH.md").write_text(report)

    manifest = {
        "format": "bubble_batch69_durability_analysis_v1",
        "status": "passed", "script_version": SCRIPT_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_is_read_only": True,
        "case_count": len(case), "unique_geometry_count": int(case.parent_case.nunique()),
        "complete_pair_count": len(paired), "variant_counts": case.variant.value_counts().to_dict(),
        "common_pv_horizon": args.common_pv_horizon,
        "source_files": {
            "batch_index": {"path": str(index_path), "sha256": sha256(index_path)},
            "input_inventory": {"path": str(inventory_path), "sha256": sha256(inventory_path)},
            "saturation_source_audit": {"path": str(audit_path), "sha256": sha256(audit_path)},
            "analysis_script": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        },
        "excluded_cases": exclusions,
        "input_run_files": input_records,
        "definitions": {
            "primary_inventory": "free gas moles; volume and saturation are secondary because gas is compressible",
            "pv50": "injected pore volumes at 50% loss of initial total free gas moles",
            "auc": f"trapezoidal mean of n_g/n_g0 over common 0-{args.common_pv_horizon:g} PV window",
            "weibull": "n_g/n_g0 = exp(-(PV/scale)^shape), fitted on a uniform PV grid",
            "spatial_weight": "initial free-gas moles",
            "bubble_count_semantics": "active gas-bearing control volumes, not reconstructed connected gas clusters",
            "batched_event_time": "estimated within accepted backward-Euler step; not an exactly re-solved event time",
        },
        "software": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__},
        "outputs": {},
    }
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "analysis_manifest.json":
            manifest["outputs"][str(path.relative_to(output))] = sha256(path)
    (output / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": "passed", "output": str(output), "cases": len(case),
                      "pairs": len(paired)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
