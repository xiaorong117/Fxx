#!/usr/bin/env python3
"""Macro input -> pore structure -> trapped-gas pattern -> bubble durability.

Read-only, offline linkage analysis for the 35-geometry / 70-run bubble
dissolution batch.  It combines three already-existing artefact families:

1. the user-supplied gradation-curve workbook
   (/workspace/zz/级配曲线数据/PSD-zhangzhe2025.1.22.xlsx);
2. the DEM/triaxial design factors already encoded in the case names
   (mu = interparticle friction coefficient, overburden stress in kPa);
3. the completed 70-case durability statistics
   (output/analysis/batch70_bubble_durability_v1).

Nothing is re-simulated.  All new files are written to a fresh output
directory; raw PNM files, production runs and verification runs are not
touched.

Statistical notes
-----------------
* Gradation parameters are constants per gradation (7 unique values), so the
  effective sample size for gradation effects is 7, not 35/70.  Gradation-level
  correlations are therefore reported with exact blocked-permutation p-values
  (7! = 5040 relabellings) and labelled exploratory.
* mu and overburden stress are within-gradation design factors with 3 ordered
  levels each; ordered effects use the Page trend test for a randomised block
  design (7 gradations as blocks).
* All correlations are Spearman rank correlations.  Between/within
  decompositions are computed by centring inside gradation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata, spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DEFAULT_PSD = Path("/workspace/zz/级配曲线数据/PSD-zhangzhe2025.1.22.xlsx")
DEFAULT_BATCH70 = REPO / "output/analysis/batch70_bubble_durability_v1"
DEFAULT_OUT = REPO / "output/analysis/macro_linkage_v1"

# Gradations that actually have PNM production runs in this batch.
BATCH_GRADATIONS = ["100B", "100C", "100D", "25ABCD", "33ABC", "60A40D", "85A15D"]
# Condition blocks, in the fixed order used throughout.
BLOCKS = [(0.1, 100.0), (0.3, 100.0), (0.5, 20.0), (0.5, 100.0), (0.5, 200.0)]
RNG_SEED = 20260917
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


# --------------------------------------------------------------------------- #
# gradation-curve workbook
# --------------------------------------------------------------------------- #
def _colnum(ref: str) -> int:
    letters = re.match(r"([A-Z]+)", ref).group(1)
    value = 0
    for ch in letters:
        value = value * 26 + ord(ch) - 64
    return value


def parse_psd_workbook(path: Path) -> dict[str, list[tuple[float, float]]]:
    """Return {gradation: [(d_mm, percent_finer), ...]} from the Test-data block.

    The workbook stores, for every gradation, a measured "Test data" curve
    (nominal sieve in col 1, measured diameter in mm in col 2, percentage finer
    in col 3) next to a rounded "DEM simulation" copy (cols 6-8).  The measured
    curve is used because it carries the non-rounded diameters.
    """
    with zipfile.ZipFile(path) as zf:
        shared = [
            "".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t"))
            for si in ET.fromstring(zf.read("xl/sharedStrings.xml")).findall("m:si", NS)
        ]
        sheet = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
    grid: dict[int, dict[int, object]] = {}
    for row in sheet.find("m:sheetData", NS).findall("m:row", NS):
        cells: dict[int, object] = {}
        for cell in row.findall("m:c", NS):
            raw = cell.find("m:v", NS)
            if raw is None or raw.text is None:
                continue
            value: object = raw.text
            if cell.get("t") == "s":
                value = shared[int(raw.text)]
            cells[_colnum(cell.get("r"))] = value
        grid[int(row.get("r"))] = cells

    headers = [
        (rn, cells[2])
        for rn, cells in sorted(grid.items())
        if isinstance(cells.get(2), str) and cells[2].startswith("Test data-")
    ]
    curves: dict[str, list[tuple[float, float]]] = {}
    for idx, (start, title) in enumerate(headers):
        stop = headers[idx + 1][0] if idx + 1 < len(headers) else max(grid) + 1
        gradation = title.replace("Test data-", "")
        points: list[tuple[float, float]] = []
        for rn in range(start + 2, stop):
            cells = grid.get(rn, {})
            try:
                d_mm = float(cells[2])
                finer = float(cells[3])
            except (KeyError, TypeError, ValueError):
                continue
            points.append((d_mm, finer))
        points.sort()
        curves[gradation] = points
    return curves


def d_at_percent(curve: list[tuple[float, float]], percent: float) -> float:
    """Diameter (mm) at a given percentage finer, log-linear interpolation.

    Inside a perfectly flat segment the geometric mid-point of the flat run is
    returned; the percentile is otherwise unidentifiable there.
    """
    d = [p[0] for p in curve]
    q = [p[1] for p in curve]
    for i in range(len(curve) - 1):
        if q[i] <= percent <= q[i + 1]:
            if q[i + 1] == q[i]:
                j = i
                while j + 1 < len(curve) and q[j + 1] == q[i]:
                    j += 1
                return float(np.sqrt(d[i] * d[j]))
            frac = (percent - q[i]) / (q[i + 1] - q[i])
            return float(np.exp(np.log(d[i]) + frac * (np.log(d[i + 1]) - np.log(d[i]))))
    raise ValueError(f"percentile {percent} outside the curve")


def percent_at_d(curve: list[tuple[float, float]], d_target: float) -> float:
    d = [p[0] for p in curve]
    q = [p[1] for p in curve]
    for i in range(len(curve) - 1):
        if d[i] <= d_target <= d[i + 1]:
            frac = (np.log(d_target) - np.log(d[i])) / (np.log(d[i + 1]) - np.log(d[i]))
            return float(q[i] + frac * (q[i + 1] - q[i]))
    raise ValueError(f"diameter {d_target} outside the curve")


def plateau_gap(curve: list[tuple[float, float]], slope_threshold: float = 10.0) -> tuple[float, float]:
    """Return (gap_fraction_percent, plateau_size_ratio) for a gap-graded curve.

    A "plateau" is a run of consecutive segments that (i) lies inside the
    d10-d90 range, (ii) starts above 50 % passing, and (iii) keeps a cumulative
    slope below `slope_threshold` percentage points per decade.  The coarse-mode
    fraction is then 100 - (mean percentage finer over the longest plateau) and
    the size ratio is the diameter ratio across that plateau.  Narrow-graded
    curves have no such plateau and return (0, 1).
    """
    d = [p[0] for p in curve]
    q = [p[1] for p in curve]
    d10 = d_at_percent(curve, 10)
    d90 = d_at_percent(curve, 90)
    best = (0.0, 0.0)  # (log10 span, gap fraction)
    i = 0
    while i < len(curve) - 1:
        if d[i + 1] <= d[i]:
            i += 1
            continue
        if d[i] < d10 or d[i + 1] > d90 or q[i] < 50.0:
            i += 1
            continue
        slope = (q[i + 1] - q[i]) / np.log10(d[i + 1] / d[i])
        if slope <= slope_threshold:
            j = i
            while j < len(curve) - 2:
                if d[j + 2] > d90:
                    break
                nxt = (q[j + 2] - q[j + 1]) / np.log10(d[j + 2] / d[j + 1])
                if nxt > slope_threshold:
                    break
                j += 1
            span = np.log10(d[j + 1] / d[i])
            plateau_mean = float(np.mean(q[i : j + 2]))
            if span > best[0]:
                best = (float(span), float(100.0 - plateau_mean))
            i = j + 1
        else:
            i += 1
    return best[1], float(10 ** best[0]) if best[0] else 1.0


def gradation_parameters(curves: dict[str, list[tuple[float, float]]]) -> pd.DataFrame:
    rows = []
    for gradation, curve in curves.items():
        if len(curve) < 3:
            continue
        d10, d30, d50, d60, d90 = (d_at_percent(curve, p) for p in (10, 30, 50, 60, 90))
        gap_fraction, gap_ratio = plateau_gap(curve)
        rows.append(
            dict(
                gradation=gradation,
                n_points=len(curve),
                d10_mm=d10,
                d30_mm=d30,
                d50_mm=d50,
                d60_mm=d60,
                d90_mm=d90,
                dmax_mm=max(p[0] for p in curve),
                Cu=d60 / d10,
                Cc=d30**2 / (d10 * d60),
                span_d90_d10=d90 / d10,
                span_ratio=(d90 - d10) / d50,
                pct_passing_1mm=percent_at_d(curve, 1.0) if curve[0][0] < 1.0 < curve[-1][0] else np.nan,
                gap_fraction_pct=gap_fraction,
                gap_size_ratio=gap_ratio,
            )
        )
    return pd.DataFrame(rows).set_index("gradation").sort_index()


# --------------------------------------------------------------------------- #
# statistics helpers
# --------------------------------------------------------------------------- #
def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 5000, seed: int = RNG_SEED):
    rho, p = spearmanr(x, y)
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            boots[b] = np.nan
            continue
        boots[b] = spearmanr(x[idx], y[idx])[0]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return float(rho), float(p), float(lo), float(hi)


def blocked_permutation_pvalue(param_by_gradation: pd.Series, outcome: np.ndarray,
                               blocks: np.ndarray, gradations: np.ndarray,
                               n_perm: int = 5040, seed: int = RNG_SEED) -> tuple[float, float]:
    """Spearman rho and permutation p-value for a gradation-level parameter."""
    frame = pd.DataFrame({"outcome": outcome, "block": blocks, "gradation": gradations})
    frame["centered"] = frame["outcome"] - frame.groupby("block")["outcome"].transform("mean")
    block_mean = frame.groupby("gradation")["centered"].mean()
    values = param_by_gradation.reindex(block_mean.index).to_numpy(float)
    rho = spearmanr(values, block_mean.to_numpy(float))[0]
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(values)
        if abs(spearmanr(perm, block_mean.to_numpy(float))[0]) >= abs(rho) - 1e-12:
            count += 1
    return float(rho), float((count + 1) / (n_perm + 1))


def page_trend(values: np.ndarray, levels: np.ndarray, blocks: np.ndarray) -> tuple[float, float, float]:
    """Page trend test (ordered alternative, increasing) for a blocked design."""
    frame = pd.DataFrame({"v": np.asarray(values, float), "l": levels, "b": blocks})
    frame["rank"] = frame.groupby("b")["v"].rank()
    n_blocks = frame.b.nunique()
    k = frame.l.nunique()
    L = sum((level + 1) * frame.loc[frame.l == level, "rank"].sum() for level in sorted(frame.l.unique()))
    mean = n_blocks * k * (k + 1) ** 2 / 4.0
    var = n_blocks * k**2 * (k + 1) * (k**2 - 1) / 144.0
    z = (L - mean) / np.sqrt(var)
    return float(L), float(z), float(norm.sf(z))


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    rxy = spearmanr(x, y)[0]
    rxz = spearmanr(x, z)[0]
    ryz = spearmanr(y, z)[0]
    den = np.sqrt(max(1e-12, (1 - rxz**2) * (1 - ryz**2)))
    return float((rxy - rxz * ryz) / den)


def centered(series: pd.Series, groups: pd.Series) -> pd.Series:
    return series - series.groupby(groups).transform("mean")


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    ranked = p[order]
    m = len(p)
    adjusted = ranked * m / (np.arange(m) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adjusted, 0, 1)
    return list(out)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


# --------------------------------------------------------------------------- #
# master tables
# --------------------------------------------------------------------------- #
def build_geometry_master(batch70: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metrics = pd.read_csv(batch70 / "case_metrics.tsv", sep="\t")
    key = ["gradation", "mu", "stress_kpa"]
    structural = [
        "control_volume_porosity",
        "mean_coordination_number",
        "pore_radius_mean_m",
        "throat_radius_mean_m",
        "constriction_ratio_p50_segment_weighted",
        "estimated_pore_reynolds_p50",
        "estimated_capillary_number",
        "raw_pore_count",
        "raw_throat_count",
        "hydraulic_conductance_sum_m3_Pa_s",
    ]
    geo = metrics.groupby(key)[structural].mean().reset_index()
    for variant, tag in (("PNM-0.6", "06"), ("PNM-0.8", "08")):
        sub = metrics[metrics.variant == variant][
            key
            + [
                "pv_to_50pct_free_gas_loss",
                "auc_free_gas_retention_5pv",
                "gas_longitudinal_variance_normalized",
                "gas_center_x_normalized",
                "gas_entropy_normalized",
                "gas_gini",
                "initial_main_Sg",
                "initial_real_active_gas_units",
                "initial_bubble_radius_p50_m",
                "free_gas_retention_at_5pv",
                "final_injected_pore_volumes",
            ]
        ].copy()
        sub.columns = key + [
            f"pv50_{tag}",
            f"auc_{tag}",
            f"gasvar_{tag}",
            f"gascenter_{tag}",
            f"entropy_{tag}",
            f"gini_{tag}",
            f"sg_{tag}",
            f"units_{tag}",
            f"rbub_{tag}",
            f"ret5_{tag}",
            f"finalpv_{tag}",
        ]
        geo = geo.merge(sub, on=key)
    geo["pv50"] = (geo.pv50_06 + geo.pv50_08) / 2
    geo["auc"] = (geo.auc_06 + geo.auc_08) / 2
    geo["sg"] = (geo.sg_06 + geo.sg_08) / 2
    geo["units"] = (geo.units_06 + geo.units_08) / 2
    geo["gasvar"] = (geo.gasvar_06 + geo.gasvar_08) / 2
    geo["gascenter"] = (geo.gascenter_06 + geo.gascenter_08) / 2
    geo["entropy"] = (geo.entropy_06 + geo.entropy_08) / 2
    geo["gini"] = (geo.gini_06 + geo.gini_08) / 2
    geo["rbub"] = (geo.rbub_06 + geo.rbub_08) / 2
    geo = geo.sort_values(key).reset_index(drop=True)
    gradation_master = geo.groupby("gradation").mean(numeric_only=True)
    return geo, gradation_master


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def figure_psd(curves: dict[str, list[tuple[float, float]]], params: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    colors = plt.cm.tab10(np.linspace(0, 1, len(BATCH_GRADATIONS)))
    for color, gradation in zip(colors, BATCH_GRADATIONS):
        curve = curves[gradation]
        d = [p[0] for p in curve]
        q = [p[1] for p in curve]
        row = params.loc[gradation]
        ax.plot(d, q, marker="o", ms=3, lw=1.6, color=color,
                label=f"{gradation}: Cu={row.Cu:.2f}, Cc={row.Cc:.2f}, d50={row.d50_mm:.2f} mm")
        for pct, marker in ((10, "v"), (30, "s"), (60, "^")):
            dd = d_at_percent(curve, pct)
            ax.plot([dd], [pct], marker=marker, ms=6, color=color, mec="k", mew=0.4, zorder=5)
    ax.plot([], [], "kv", label="d10")
    ax.plot([], [], "ks", label="d30")
    ax.plot([], [], "k^", label="d60")
    ax.set_xscale("log")
    ax.set_xlabel("particle diameter d (mm, log scale)")
    ax.set_ylabel("percentage finer /%")
    ax.set_ylim(-3, 103)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.6, loc="upper left", framealpha=0.9)
    ax.set_title("DEM/test gradation curves of the 7 simulated gradations")
    fig.tight_layout()
    fig.savefig(out / "01_psd_curves.png", dpi=200)
    plt.close(fig)


def figure_gradation_links(gradation_master: pd.DataFrame, stats: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.0))
    panels = [
        ("Cu", "control_volume_porosity", "porosity /-", axes[0, 0]),
        ("Cu", "auc", "0-5 PV AUC /-", axes[0, 1]),
        ("span_d90_d10", "control_volume_porosity", "porosity /-", axes[1, 0]),
        ("span_d90_d10", "auc", "0-5 PV AUC /-", axes[1, 1]),
    ]
    for xcol, ycol, ylabel, ax in panels:
        x = gradation_master[xcol]
        y = gradation_master[ycol]
        ax.scatter(x, y, s=48, color="#1f77b4", zorder=3)
        for gradation in gradation_master.index:
            ax.annotate(gradation, (x[gradation], y[gradation]), fontsize=7.5,
                        xytext=(4, 3), textcoords="offset points")
        row = stats[(stats.param == xcol) & (stats.outcome == ycol)].iloc[0]
        ax.set_title(f"{xcol} vs {ycol}\nrho={row.rho_gradation_means:+.3f}  "
                     f"blocked perm p={row.perm_p_blocked:.4f}", fontsize=9)
        ax.set_xlabel(xcol)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    fig.suptitle("Gradation parameters vs packing porosity and early retention (7 gradations, exploratory)")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out / "02_gradation_params_vs_structure_durability.png", dpi=200)
    plt.close(fig)


def figure_compaction_stress(geo: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.6))
    colors = plt.cm.tab10(np.linspace(0, 1, len(BATCH_GRADATIONS)))
    cmap = dict(zip(BATCH_GRADATIONS, colors))
    mu_sub = geo[geo.stress_kpa == 100]
    stress_sub = geo[geo.mu == 0.5]
    specs = [
        (axes[0, 0], mu_sub, "mu", "control_volume_porosity", "porosity /-", "mu (DEM friction) vs porosity @100 kPa"),
        (axes[0, 1], mu_sub, "mu", "pv50", "PV50 /PV", "mu vs PV50 @100 kPa"),
        (axes[0, 2], mu_sub, "mu", "auc", "0-5 PV AUC /-", "mu vs AUC @100 kPa"),
        (axes[1, 0], stress_sub, "stress_kpa", "control_volume_porosity", "porosity /-", "overburden vs porosity @mu=0.5"),
        (axes[1, 1], stress_sub, "stress_kpa", "pv50", "PV50 /PV", "overburden vs PV50 @mu=0.5"),
        (axes[1, 2], stress_sub, "stress_kpa", "auc", "0-5 PV AUC /-", "overburden vs AUC @mu=0.5"),
    ]
    for ax, sub, xcol, ycol, ylabel, title in specs:
        for gradation, part in sub.groupby("gradation"):
            part = part.sort_values(xcol)
            ax.plot(part[xcol], part[ycol], marker="o", ms=4, lw=1.3,
                    color=cmap[gradation], label=gradation)
        ax.set_xlabel(xcol)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=9)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle("Compaction (mu) and overburden pressure within the same gradation")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out / "03_mu_stress_panels.png", dpi=200)
    plt.close(fig)


def figure_between_within(geo: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 8.2))
    cmap = dict(zip(BATCH_GRADATIONS, plt.cm.tab10(np.linspace(0, 1, len(BATCH_GRADATIONS)))))
    for col, (ycol, ylabel) in enumerate([("auc", "0-5 PV AUC /-"), ("pv50", "PV50 /PV")]):
        ax = axes[0, col]
        for gradation, part in geo.groupby("gradation"):
            ax.scatter(part.control_volume_porosity, part[ycol], s=34, color=cmap[gradation], label=gradation)
        means = geo.groupby("gradation")[["control_volume_porosity", ycol]].mean()
        ax.scatter(means.control_volume_porosity, means[ycol], s=90, facecolor="none",
                   edgecolor="k", linewidths=1.2, zorder=5)
        rho, p, lo, hi = spearman_ci(means.control_volume_porosity.to_numpy(), means[ycol].to_numpy())
        ax.set_title(f"between gradations (n=7): rho={rho:+.3f}, p={p:.3f}", fontsize=9)
        ax.set_xlabel("porosity /-")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        if col == 0:
            ax.legend(fontsize=6.6, ncol=2)
        ax = axes[1, col]
        wp = centered(geo.control_volume_porosity, geo.gradation)
        wy = centered(geo[ycol], geo.gradation)
        for gradation, part in geo.groupby("gradation"):
            idx = part.index
            ax.scatter(wp[idx], wy[idx], s=34, color=cmap[gradation])
        rho, p, lo, hi = spearman_ci(wp.to_numpy(), wy.to_numpy())
        ax.set_title(f"within gradations, centred (n=35): rho={rho:+.3f}, p={p:.4f}", fontsize=9)
        ax.set_xlabel("porosity - gradation mean")
        ax.set_ylabel(f"{ylabel} - gradation mean")
        ax.axhline(0, color="k", lw=0.6)
        ax.axvline(0, color="k", lw=0.6)
        ax.grid(alpha=0.3)
    fig.suptitle("Porosity-durability link: between- vs within-gradation components")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out / "04_porosity_durability_between_within.png", dpi=200)
    plt.close(fig)


def figure_gas_chain(geo: pd.DataFrame, gradation_master: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 8.2))
    cmap = dict(zip(BATCH_GRADATIONS, plt.cm.tab10(np.linspace(0, 1, len(BATCH_GRADATIONS)))))
    panels = [
        (axes[0, 0], "control_volume_porosity", "gasvar", "porosity /-", "gas longitudinal variance /-",
         "porosity vs initial gas spread (35 geometries)"),
        (axes[0, 1], "gasvar", "pv50", "gas longitudinal variance /-", "PV50 /PV",
         "initial gas spread vs PV50 (35 geometries)"),
        (axes[1, 0], "control_volume_porosity", "gascenter", "porosity /-", "gas centroid (0=inlet,1=outlet)",
         "porosity vs initial gas centroid"),
        (axes[1, 1], "gascenter", "auc", "gas centroid (0=inlet,1=outlet)", "0-5 PV AUC /-",
         "gas centroid vs AUC"),
    ]
    for ax, xcol, ycol, xlabel, ylabel, title in panels:
        for gradation, part in geo.groupby("gradation"):
            ax.scatter(part[xcol], part[ycol], s=30, color=cmap[gradation], label=gradation)
        rho, p, lo, hi = spearman_ci(geo[xcol].to_numpy(), geo[ycol].to_numpy())
        ax.set_title(f"{title}\nrho={rho:+.3f}, p={p:.4f}", fontsize=9)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=6.6, ncol=2)
    fig.suptitle("Structure -> initial trapped-gas pattern -> durability")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out / "05_structure_gas_durability_chain.png", dpi=200)
    plt.close(fig)


def figure_summary_heatmap(gradation_master: pd.DataFrame, out: Path) -> None:
    cols = ["d50_mm", "Cu", "Cc", "span_d90_d10", "control_volume_porosity",
            "mean_coordination_number", "sg", "gasvar", "gascenter", "entropy",
            "gini", "pv50", "auc"]
    frame = gradation_master[cols].copy()
    z = (frame - frame.mean()) / frame.std(ddof=0)
    fig, ax = plt.subplots(figsize=(11.6, 4.6))
    im = ax.imshow(z.to_numpy(), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(z.index)))
    ax.set_yticklabels(z.index, fontsize=8)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            ax.text(j, i, f"{frame.iloc[i, j]:.2f}", ha="center", va="center", fontsize=6.8)
    fig.colorbar(im, ax=ax, label="z-score across gradations")
    ax.set_title("Gradation-level summary (values shown; colours are z-scores)")
    fig.tight_layout()
    fig.savefig(out / "06_gradation_summary_heatmap.png", dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psd", type=Path, default=DEFAULT_PSD)
    parser.add_argument("--batch70", type=Path, default=DEFAULT_BATCH70)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--n-boot", type=int, default=5000)
    args = parser.parse_args()

    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    figures_dir = out / "figures"

    curves = parse_psd_workbook(args.psd)
    params = gradation_parameters(curves)
    params.to_csv(out / "gradation_parameters.csv")
    pd.DataFrame(
        [(g, d, p) for g, pts in curves.items() for d, p in pts],
        columns=["gradation", "d_mm", "percent_finer"],
    ).to_csv(out / "psd_curves_digitized.csv", index=False)

    geo, gradation_master = build_geometry_master(args.batch70)
    merge_params = params[["d10_mm", "d30_mm", "d50_mm", "d60_mm", "d90_mm", "Cu", "Cc",
                           "span_d90_d10", "span_ratio", "gap_fraction_pct", "gap_size_ratio"]]
    geo = geo.merge(merge_params, left_on="gradation", right_index=True, how="left")
    gradation_master = gradation_master.join(merge_params, how="left")
    geo.to_csv(out / "geometry_master.csv", index=False)
    gradation_master.to_csv(out / "gradation_master.csv")

    # ---- design checks ---------------------------------------------------- #
    design = geo[["gradation", "mu", "stress_kpa", "control_volume_porosity"]].copy()
    design.to_csv(out / "design_factors_and_porosity.csv", index=False)

    # ---- gradation-level statistics (n=7, blocked permutation) ------------ #
    # Condition block identity (mu, overburden) is the correct blocking unit:
    # the 7 gradations are compared inside each of the 5 identical conditions.
    block_id = (geo.mu.astype(str) + "kPa|" + geo.stress_kpa.astype(str)).to_numpy()
    gradation_rows = []
    for xcol in ["d50_mm", "Cu", "Cc", "span_d90_d10", "gap_fraction_pct"]:
        for ycol in ["control_volume_porosity", "mean_coordination_number", "auc", "pv50",
                     "sg", "gasvar", "gascenter", "entropy", "gini"]:
            rho_blocked, perm_p = blocked_permutation_pvalue(
                params[xcol], geo[ycol].to_numpy(), block_id, geo.gradation.to_numpy()
            )
            rho7 = spearmanr(params[xcol].reindex(gradation_master.index), gradation_master[ycol])[0]
            gradation_rows.append(dict(param=xcol, outcome=ycol,
                                       rho_gradation_means=rho7,
                                       rho_block_centered=rho_blocked,
                                       perm_p_blocked=perm_p))
    gradation_stats = pd.DataFrame(gradation_rows)
    gradation_stats["fdr_bh"] = bh_fdr(list(gradation_stats.perm_p_blocked))
    gradation_stats.to_csv(out / "gradation_vs_structure_durability.tsv", sep="\t", index=False)

    # ---- robustness: the gradation ranking inside every condition block ---- #
    consistency_rows = []
    for xcol in ["d50_mm", "Cu", "Cc", "span_d90_d10", "gap_fraction_pct"]:
        for ycol in ["control_volume_porosity", "auc", "pv50", "sg", "gascenter"]:
            block_rhos = {}
            for block, part in geo.groupby(block_id):
                part = part.sort_values("gradation")
                block_rhos[block] = spearmanr(part[xcol], part[ycol])[0]
            values = np.array(list(block_rhos.values()), float)
            consistency_rows.append(dict(param=xcol, outcome=ycol,
                                         n_blocks=len(values),
                                         rho_min=values.min(), rho_mean=values.mean(),
                                         rho_max=values.max(),
                                         blocks=";".join(f"{k}:{v:+.2f}" for k, v in block_rhos.items())))
    pd.DataFrame(consistency_rows).to_csv(out / "block_consistency.tsv", sep="\t", index=False)

    # ---- mu / stress ordered effects -------------------------------------- #
    trend_rows = []
    mu_sub = geo[geo.stress_kpa == 100].copy()
    stress_sub = geo[geo.mu == 0.5].copy()
    for label, frame, level_col, levels in (
        ("mu@100kPa", mu_sub, "mu", {0.1: 0, 0.3: 1, 0.5: 2}),
        ("overburden@mu=0.5", stress_sub, "stress_kpa", {20: 0, 100: 1, 200: 2}),
    ):
        codes = frame[level_col].map(levels).to_numpy()
        for ycol in ["control_volume_porosity", "pv50", "pv50_06", "pv50_08", "auc", "auc_06", "auc_08",
                     "ret5_06", "ret5_08", "sg", "gasvar", "gascenter"]:
            L, z, p_up = page_trend(frame[ycol].to_numpy(), codes, frame.gradation.to_numpy())
            _, z_dn, p_dn = page_trend(-frame[ycol].to_numpy(), codes, frame.gradation.to_numpy())
            trend_rows.append(dict(design=label, outcome=ycol, page_L=L, z=z,
                                   p_increasing=p_up, p_decreasing=p_dn))
    trends = pd.DataFrame(trend_rows)
    trends.to_csv(out / "ordered_trend_tests.tsv", sep="\t", index=False)

    # ---- porosity-durability decomposition -------------------------------- #
    decomp_rows = []
    for ycol in ["pv50", "pv50_06", "pv50_08", "auc", "auc_06", "auc_08"]:
        rho_t, p_t, lo_t, hi_t = spearman_ci(geo.control_volume_porosity.to_numpy(), geo[ycol].to_numpy(), args.n_boot)
        means = geo.groupby("gradation")[["control_volume_porosity", ycol]].mean()
        rho_b, p_b, lo_b, hi_b = spearman_ci(means.control_volume_porosity.to_numpy(), means[ycol].to_numpy(), args.n_boot)
        wp = centered(geo.control_volume_porosity, geo.gradation).to_numpy()
        wy = centered(geo[ycol], geo.gradation).to_numpy()
        rho_w, p_w, lo_w, hi_w = spearman_ci(wp, wy, args.n_boot)
        decomp_rows.append(dict(outcome=ycol, level="total", n=35, rho=rho_t, p=p_t, ci_low=lo_t, ci_high=hi_t))
        decomp_rows.append(dict(outcome=ycol, level="between_gradation", n=7, rho=rho_b, p=p_b, ci_low=lo_b, ci_high=hi_b))
        decomp_rows.append(dict(outcome=ycol, level="within_gradation", n=35, rho=rho_w, p=p_w, ci_low=lo_w, ci_high=hi_w))
    decomp = pd.DataFrame(decomp_rows)
    decomp.to_csv(out / "porosity_durability_decomposition.tsv", sep="\t", index=False)

    # ---- structure -> gas pattern -> durability --------------------------- #
    gas_rows = []
    pairs = [("control_volume_porosity", "sg"), ("control_volume_porosity", "gasvar"),
             ("control_volume_porosity", "gascenter"), ("control_volume_porosity", "entropy"),
             ("mean_coordination_number", "sg"), ("pore_radius_mean_m", "units"),
             ("throat_radius_mean_m", "units"), ("sg", "pv50"), ("sg", "auc"),
             ("gasvar", "pv50"), ("gasvar", "auc"), ("gascenter", "pv50"), ("gascenter", "auc"),
             ("entropy", "auc"), ("gini", "auc")]
    for xcol, ycol in pairs:
        rho_t, p_t, lo_t, hi_t = spearman_ci(geo[xcol].to_numpy(), geo[ycol].to_numpy(), args.n_boot)
        wxa = centered(geo[xcol], geo.gradation).to_numpy()
        wya = centered(geo[ycol], geo.gradation).to_numpy()
        rho_w, p_w, lo_w, hi_w = spearman_ci(wxa, wya, args.n_boot)
        gas_rows.append(dict(x=xcol, y=ycol, level="total", rho=rho_t, p=p_t, ci_low=lo_t, ci_high=hi_t))
        gas_rows.append(dict(x=xcol, y=ycol, level="within_gradation", rho=rho_w, p=p_w, ci_low=lo_w, ci_high=hi_w))
        gas_rows.append(dict(x=xcol, y=ycol, level="partial_porosity",
                             rho=partial_spearman(geo[xcol].to_numpy(), geo[ycol].to_numpy(),
                                                  geo.control_volume_porosity.to_numpy()),
                             p=np.nan, ci_low=np.nan, ci_high=np.nan))
    gas_stats = pd.DataFrame(gas_rows)
    gas_stats.to_csv(out / "structure_gas_durability_links.tsv", sep="\t", index=False)

    # ---- mediation / partial correlation chains --------------------------- #
    mediator_rows = []
    for xcol, mcol, ycol in [("Cu", "control_volume_porosity", "auc"),
                             ("Cu", "control_volume_porosity", "pv50"),
                             ("span_d90_d10", "control_volume_porosity", "auc"),
                             ("span_d90_d10", "control_volume_porosity", "pv50")]:
        a = spearmanr(gradation_master[xcol], gradation_master[mcol])[0]
        b = spearmanr(gradation_master[mcol], gradation_master[ycol])[0]
        c = spearmanr(gradation_master[xcol], gradation_master[ycol])[0]
        cp = partial_spearman(gradation_master[xcol].to_numpy(), gradation_master[ycol].to_numpy(),
                              gradation_master[mcol].to_numpy())
        bp = partial_spearman(gradation_master[mcol].to_numpy(), gradation_master[ycol].to_numpy(),
                              gradation_master[xcol].to_numpy())
        mediator_rows.append(dict(exposure=xcol, mediator=mcol, outcome=ycol,
                                  path_a=a, path_b=b, total_c=c,
                                  direct_c_given_m=cp, b_given_exposure=bp, indirect_a_times_b=a * b))
    # within-gradation stress mediation (35 geometries, centred)
    stress_code = geo.stress_kpa.map({20: 0, 100: 1, 200: 2})
    for ycol in ["auc", "pv50"]:
        s_c = centered(stress_code, geo.gradation).to_numpy()
        p_c = centered(geo.control_volume_porosity, geo.gradation).to_numpy()
        y_c = centered(geo[ycol], geo.gradation).to_numpy()
        mediator_rows.append(dict(exposure="overburden_rank", mediator="porosity",
                                  outcome=ycol + "_within_gradation",
                                  path_a=spearmanr(s_c, p_c)[0], path_b=spearmanr(p_c, y_c)[0],
                                  total_c=spearmanr(s_c, y_c)[0],
                                  direct_c_given_m=partial_spearman(s_c, y_c, p_c),
                                  b_given_exposure=partial_spearman(p_c, y_c, s_c),
                                  indirect_a_times_b=spearmanr(s_c, p_c)[0] * spearmanr(p_c, y_c)[0]))
    mediator = pd.DataFrame(mediator_rows)
    mediator.to_csv(out / "mediation_chains.tsv", sep="\t", index=False)

    # ---- unconditional within-variant gas -> durability (from batch70) ---- #
    vari_rows = []
    for tag in ("06", "08"):
        for xcol in ["sg", "gasvar", "gascenter", "entropy", "gini", "units", "rbub"]:
            x = geo[f"{xcol}_{tag}"].to_numpy()
            for ycol in [f"pv50_{tag}", f"auc_{tag}"]:
                rho, p, lo, hi = spearman_ci(x, geo[ycol].to_numpy(), args.n_boot)
                vari_rows.append(dict(variant=f"PNM-0.{tag[1]}", x=xcol, y=ycol, n=35,
                                      rho=rho, p=p, ci_low=lo, ci_high=hi))
    pd.DataFrame(vari_rows).to_csv(out / "variant_gas_durability_correlations.tsv", sep="\t", index=False)

    # ---- figures ---------------------------------------------------------- #
    figure_psd(curves, params, figures_dir)
    figure_gradation_links(gradation_master, gradation_stats, figures_dir)
    figure_compaction_stress(geo, figures_dir)
    figure_between_within(geo, figures_dir)
    figure_gas_chain(geo, gradation_master, figures_dir)
    figure_summary_heatmap(gradation_master, figures_dir)

    # ---- manifest --------------------------------------------------------- #
    manifest = {
        "format": "bubble_macro_linkage_analysis_v1",
        "status": "passed",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_is_read_only": True,
        "analysis_script": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "inputs": {
            "psd_workbook": {"path": str(args.psd), "sha256": sha256(args.psd)},
            "batch70_case_metrics": {
                "path": str(args.batch70 / "case_metrics.tsv"),
                "sha256": sha256(args.batch70 / "case_metrics.tsv"),
            },
            "batch70_manifest": {
                "path": str(args.batch70 / "analysis_manifest.json"),
                "sha256": sha256(args.batch70 / "analysis_manifest.json"),
            },
        },
        "design": {
            "gradations": BATCH_GRADATIONS,
            "conditions": [{"mu": mu, "stress_kpa": s} for mu, s in BLOCKS],
            "geometries": int(len(geo)),
            "durability_runs_per_geometry": 2,
        },
        "methods": {
            "gradation_parameters": "d10/d30/d50/d60/d90 by log-linear interpolation of the measured test curve; Cu=d60/d10; Cc=d30^2/(d10*d60); span=d90/d10; gap fraction from the longest cumulative plateau (<10 pp/decade)",
            "gradation_inference": "Spearman rho of the 7 gradation values against block-centred outcome means; exact permutation over 5040 gradation relabellings (blocked design)",
            "design_factor_inference": "Page trend test for ordered levels within gradation blocks",
            "decomposition": "total / between-gradation / within-gradation (centred) Spearman correlations",
            "mediation": "rank partial correlations (exploratory, n=7 for gradation-level chains; n=35 centred for stress chains)",
        },
        "outputs": {
            p.name: sha256(p) for p in sorted(out.glob("*")) if p.is_file()
        },
        "figure_outputs": {
            p.name: sha256(p) for p in sorted(figures_dir.glob("*.png"))
        },
    }
    (out / "analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(f"wrote {out}")
    print(gradation_stats.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
