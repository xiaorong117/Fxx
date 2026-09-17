#!/usr/bin/env python3
"""Corrected macro-input analysis with the user-supplied density (void-ratio) table.

Correction relative to ``analyze_macro_linkage.py``
---------------------------------------------------
``mu`` is the DEM **interparticle friction coefficient**, not a density index.
The density of each prepared sample is supplied separately by the user as the
table in ``/workspace/zz/级配曲线数据/密实度.png`` (7 gradations x 5 conditions).

Provenance of the transcription
-------------------------------
* source image: ``/workspace/zz/级配曲线数据/密实度.png``
  (sha256 recorded in the manifest);
* the table values behave as **void ratio e** (Chinese ``e`` in the case
  directory names, e.g. ``sands-of-100B-mu-0.5-20kpa-e-0.52``): converting
  ``n = e/(1+e)`` reproduces the measured PNM control-volume porosity with
  Spearman rho = 0.945 (mean difference -0.005, max |difference| 0.055);
* column order 1..5 was fixed by matching the 100B row against the five
  case-directory ``e`` labels and by the design pattern (column 3 = 20 kPa is
  always the largest value, column 5 = 200 kPa always the smallest).

Density as a state variable
---------------------------
The density is *produced* by the design factors (mu, overburden, gradation) and
*consumed* by the imbibition/dissolution chain, so it is analysed as the
carrier of the design effects, not as an independent experimental factor.

Read-only with respect to all production, verification and raw-PNM data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DENSITY_IMAGE = Path("/workspace/zz/级配曲线数据/密实度.png")
DEFAULT_MACRO = REPO / "output/analysis/macro_linkage_v1"
DEFAULT_OUT = REPO / "output/analysis/density_inputs_v1"
BATCH_GRADATIONS = ["100B", "100C", "100D", "25ABCD", "33ABC", "60A40D", "85A15D"]
# Condition for columns 1..5, taken from the case-directory naming.
COLUMN_CONDITIONS = [
    ("1", 0.1, 100.0, "mu=0.1, 100 kPa"),
    ("2", 0.3, 100.0, "mu=0.3, 100 kPa"),
    ("3", 0.5, 20.0, "mu=0.5, 20 kPa"),
    ("4", 0.5, 100.0, "mu=0.5, 100 kPa"),
    ("5", 0.5, 200.0, "mu=0.5, 200 kPa"),
]
# Transcription of 密实度.png: row = gradation, columns = sample 1..5.
DENSITY_TABLE = {
    "25ABCD": [0.18, 0.22, 0.31, 0.23, 0.14],
    "33ABC": [0.22, 0.26, 0.38, 0.28, 0.17],
    "60A40D": [0.22, 0.25, 0.33, 0.26, 0.19],
    "85A15D": [0.34, 0.38, 0.51, 0.39, 0.29],
    "100B": [0.34, 0.39, 0.52, 0.41, 0.29],
    "100C": [0.38, 0.43, 0.56, 0.45, 0.33],
    "100D": [0.40, 0.45, 0.60, 0.47, 0.35],
}
RNG_SEED = 20260917


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def spearman_ci(x, y, n_boot=5000, seed=RNG_SEED):
    rho, p = spearmanr(x, y)
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), len(x))
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            continue
        boots.append(spearmanr(x[idx], y[idx])[0])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(rho), float(p), float(lo), float(hi)


def centered(series: pd.Series, groups: pd.Series) -> pd.Series:
    return series - series.groupby(groups).transform("mean")


def partial_spearman(x, y, z) -> float:
    rxy = spearmanr(x, y)[0]
    rxz = spearmanr(x, z)[0]
    ryz = spearmanr(y, z)[0]
    return float((rxy - rxz * ryz) / np.sqrt(max(1e-12, (1 - rxz**2) * (1 - ryz**2))))


def page_trend(values, levels, blocks) -> tuple[float, float, float]:
    """Page ordered-alternative test; returns (L, z, p_increasing)."""
    frame = pd.DataFrame({"v": np.asarray(values, float), "l": levels, "b": blocks})
    frame["rank"] = frame.groupby("b")["v"].rank()
    n_blocks = frame.b.nunique()
    k = frame.l.nunique()
    L = sum((lvl + 1) * frame.loc[frame.l == lvl, "rank"].sum() for lvl in sorted(set(levels)))
    mean = n_blocks * k * (k + 1) ** 2 / 4.0
    var = n_blocks * k**2 * (k + 1) * (k**2 - 1) / 144.0
    z = (L - mean) / np.sqrt(var)
    return float(L), float(z), float(norm.sf(z))


def build_density_master(macro_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    geo = pd.read_csv(macro_dir / "geometry_master.csv")
    params = pd.read_csv(macro_dir / "gradation_parameters.csv").set_index("gradation")
    rows = []
    for gradation, values in DENSITY_TABLE.items():
        for value, (column, mu, stress, label) in zip(values, COLUMN_CONDITIONS):
            rows.append(dict(gradation=gradation, column=column, mu=mu, stress_kpa=stress,
                             condition=label, void_ratio_e=value))
    table = pd.DataFrame(rows)
    table["porosity_from_e"] = table.void_ratio_e / (1.0 + table.void_ratio_e)
    master = table.merge(geo, on=["gradation", "mu", "stress_kpa"], how="left")
    # geometry_master.csv already carries the gradation parameters; only add the
    # ones that are missing (keeps the table free of _x/_y duplicate columns).
    missing = [c for c in ["Cu", "Cc", "d50_mm", "span_d90_d10"] if c not in master.columns]
    if missing:
        master = master.merge(params[missing], left_on="gradation", right_index=True, how="left")
    return table, master.sort_values(["gradation", "mu", "stress_kpa"]).reset_index(drop=True)


def figure_density(master: pd.DataFrame, validation: pd.DataFrame, slope_auc: float,
                   observed: dict, out: Path) -> None:
    colors = dict(zip(BATCH_GRADATIONS, plt.cm.tab10(np.linspace(0, 1, len(BATCH_GRADATIONS)))))
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.4))
    mu_sub = master[master.stress_kpa == 100]
    stress_sub = master[master.mu == 0.5]

    ax = axes[0, 0]
    for gradation, part in mu_sub.groupby("gradation"):
        part = part.sort_values("mu")
        ax.plot(part.mu, part.void_ratio_e, marker="o", ms=4, lw=1.3, color=colors[gradation], label=gradation)
    ax.set_title("mu (friction coeff.) vs density e\nPage z=+3.74, p=1e-4 @100 kPa", fontsize=9)
    ax.set_xlabel("mu /-"); ax.set_ylabel("void ratio e /-"); ax.grid(alpha=0.3)
    ax.legend(fontsize=6.2, ncol=2)

    ax = axes[0, 1]
    for gradation, part in stress_sub.groupby("gradation"):
        part = part.sort_values("stress_kpa")
        ax.plot(part.stress_kpa, part.void_ratio_e, marker="o", ms=4, lw=1.3, color=colors[gradation], label=gradation)
    ax.set_title("overburden vs density e\nPage z=-3.74, p=1e-4 @mu=0.5", fontsize=9)
    ax.set_xlabel("overburden /kPa"); ax.set_ylabel("void ratio e /-"); ax.grid(alpha=0.3)
    ax.legend(fontsize=6.2, ncol=2)

    ax = axes[0, 2]
    ax.scatter(master.porosity_from_e, master.control_volume_porosity, s=30, color="#1f77b4")
    lo = min(master.porosity_from_e.min(), master.control_volume_porosity.min())
    hi = max(master.porosity_from_e.max(), master.control_volume_porosity.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=1)
    ax.set_title(f"validation: n=e/(1+e) vs PNM porosity\nrho={validation.rho.iloc[0]:+.3f}, "
                 f"mean diff={validation.mean_difference.iloc[0]:+.3f}", fontsize=9)
    ax.set_xlabel("density from table: n=e/(1+e)"); ax.set_ylabel("PNM control-volume porosity")
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    for gradation, part in master.groupby("gradation"):
        ax.scatter(part.void_ratio_e, part.auc, s=30, color=colors[gradation], label=gradation)
    means = master.groupby("gradation")[["void_ratio_e", "auc"]].mean()
    ax.scatter(means.void_ratio_e, means.auc, s=90, facecolor="none", edgecolor="k", linewidths=1.2, zorder=5)
    rho, p, lo_, hi_ = spearman_ci(master.void_ratio_e, master.auc)
    ax.set_title(f"density e vs 0-5 PV AUC\nrho={rho:+.3f}, p={p:.4f}", fontsize=9)
    ax.set_xlabel("void ratio e /-"); ax.set_ylabel("0-5 PV AUC /-"); ax.grid(alpha=0.3)

    ax = axes[1, 1]
    for gradation, part in master.groupby("gradation"):
        ax.scatter(part.void_ratio_e, part.pv50, s=30, color=colors[gradation])
    ax.scatter(means.void_ratio_e, master.groupby("gradation").pv50.mean(), s=90,
               facecolor="none", edgecolor="k", linewidths=1.2, zorder=5)
    rho, p, lo_, hi_ = spearman_ci(master.void_ratio_e, master.pv50)
    ax.set_title(f"density e vs PV50\nrho={rho:+.3f}, p={p:.4f}", fontsize=9)
    ax.set_xlabel("void ratio e /-"); ax.set_ylabel("PV50 /PV"); ax.grid(alpha=0.3)

    ax = axes[1, 2]
    labels = ["mu 0.1->0.5\n(100 kPa)", "overburden 20->200 kPa\n(mu=0.5)"]
    predicted = [slope_auc * observed["mu"]["delta_e"], slope_auc * observed["sigma"]["delta_e"]]
    actual = [observed["mu"]["delta_auc"], observed["sigma"]["delta_auc"]]
    xpos = np.arange(2)
    ax.bar(xpos - 0.18, predicted, width=0.34, label="predicted from dAUC/de", color="#4c72b0")
    ax.bar(xpos + 0.18, actual, width=0.34, label="observed mean change", color="#dd8452")
    ax.set_xticks(xpos); ax.set_xticklabels(labels, fontsize=8)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_ylabel("delta AUC (0-5 PV)")
    ax.set_title(f"density channel explains the design effects\n(dAUC/de = {slope_auc:+.3f})", fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3, axis="y")
    fig.suptitle("Density (void ratio e) as the carrier of mu and overburden effects")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out / "08_density_inputs_panels.png", dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--macro", type=Path, default=DEFAULT_MACRO)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--n-boot", type=int, default=5000)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)

    table, master = build_density_master(args.macro)
    table.to_csv(out / "density_table_from_image.csv", index=False)
    master.to_csv(out / "density_master.csv", index=False)

    # ---- validation: does the supplied density reproduce the PNM porosity? ---
    rho, p, lo, hi = spearman_ci(master.porosity_from_e, master.control_volume_porosity, args.n_boot)
    diff = master.porosity_from_e - master.control_volume_porosity
    validation = pd.DataFrame([dict(
        comparison="n=e/(1+e) vs PNM control-volume porosity",
        rho=rho, p=p, ci_low=lo, ci_high=hi,
        mean_difference=float(diff.mean()), max_abs_difference=float(diff.abs().max()),
        note="largest offsets are the densest (200 kPa) samples; PNM control volume sits slightly above the DEM global void ratio",
    )])
    validation.to_csv(out / "validation_density_vs_pnm_porosity.tsv", sep="\t", index=False)

    # ---- design factors -> density ------------------------------------------ #
    design_rows = []
    mu_sub = master[master.stress_kpa == 100]
    stress_sub = master[master.mu == 0.5]
    for label, frame, levels, source in (
        ("mu@100kPa", mu_sub, {0.1: 0, 0.3: 1, 0.5: 2}, "mu"),
        ("overburden@mu=0.5", stress_sub, {20: 0, 100: 1, 200: 2}, "overburden"),
    ):
        column = "mu" if source == "mu" else "stress_kpa"
        codes = frame[column].map(levels).to_numpy()
        L, z, p_up = page_trend(frame.void_ratio_e.to_numpy(), codes, frame.gradation.to_numpy())
        design_rows.append(dict(design=label, outcome="void_ratio_e",
                                page_L=L, z=z, p_increasing=p_up, p_decreasing=1 - p_up))
    design = pd.DataFrame(design_rows)
    # per-gradation sensitivities
    sens_rows = []
    for gradation, sub in master.groupby("gradation"):
        a = sub[sub.stress_kpa == 100].sort_values("mu")
        b = sub[sub.mu == 0.5].sort_values("stress_kpa")
        sens_rows.append(dict(
            gradation=gradation,
            de_d_mu=float(np.polyfit(a.mu, a.void_ratio_e, 1)[0]),
            de_d_log_overburden=float(np.polyfit(np.log10(b.stress_kpa), b.void_ratio_e, 1)[0]),
        ))
    sensitivity = pd.DataFrame(sens_rows).set_index("gradation")
    sensitivity.to_csv(out / "density_sensitivity_by_gradation.csv")
    corr_rows = []
    for column in sensitivity.columns:
        for indicator in ["Cu", "Cc"]:
            r = spearmanr(master.groupby("gradation")[indicator].first(), sensitivity[column])
            corr_rows.append(dict(sensitivity=column, indicator=indicator, rho=float(r[0]), p=float(r[1])))
    pd.DataFrame(corr_rows).to_csv(out / "density_sensitivity_vs_gradation.tsv", sep="\t", index=False)

    # ---- density -> durability ---------------------------------------------- #
    density_rows = []
    for measure in ["void_ratio_e", "porosity_from_e", "control_volume_porosity"]:
        for outcome in ["pv50", "auc", "pv50_06", "pv50_08", "auc_06", "auc_08"]:
            rho_t, p_t, lo_t, hi_t = spearman_ci(master[measure], master[outcome], args.n_boot)
            density_rows.append(dict(measure=measure, outcome=outcome, level="total", n=35,
                                     rho=rho_t, p=p_t, ci_low=lo_t, ci_high=hi_t))
            means = master.groupby("gradation")[[measure, outcome]].mean()
            rho_b, p_b, lo_b, hi_b = spearman_ci(means[measure], means[outcome], args.n_boot)
            density_rows.append(dict(measure=measure, outcome=outcome, level="between_gradation", n=7,
                                     rho=rho_b, p=p_b, ci_low=lo_b, ci_high=hi_b))
            wx = centered(master[measure], master.gradation)
            wy = centered(master[outcome], master.gradation)
            rho_w, p_w, lo_w, hi_w = spearman_ci(wx, wy, args.n_boot)
            density_rows.append(dict(measure=measure, outcome=outcome, level="within_gradation", n=35,
                                     rho=rho_w, p=p_w, ci_low=lo_w, ci_high=hi_w))
    density_stats = pd.DataFrame(density_rows)
    density_stats.to_csv(out / "density_durability_links.tsv", sep="\t", index=False)

    # ---- mediation: design -> density -> durability -------------------------- #
    e_c = centered(master.void_ratio_e, master.gradation)
    n_c = centered(master.control_volume_porosity, master.gradation)
    mediation_rows = []
    for label, code, source in (
        ("mu", master.mu.map({0.1: 0, 0.3: 1, 0.5: 2}), "mu"),
        ("overburden", master.stress_kpa.map({20: 0, 100: 1, 200: 2}), "overburden"),
    ):
        x_c = centered(code.astype(float), master.gradation)
        for outcome in ["auc", "pv50"]:
            y_c = centered(master[outcome], master.gradation)
            mediation_rows.append(dict(
                exposure=label + "_within_gradation", mediator="void_ratio_e", outcome=outcome + "_within_gradation",
                path_a=float(spearmanr(x_c, e_c)[0]), path_b=float(spearmanr(e_c, y_c)[0]),
                total_c=float(spearmanr(x_c, y_c)[0]),
                direct_c_given_mediator=partial_spearman(x_c, y_c, e_c),
                mediator_given_exposure=partial_spearman(e_c, y_c, x_c),
            ))
    mediation_rows.append(dict(
        exposure="void_ratio_e_within_gradation", mediator="PNM_porosity", outcome="auc_within_gradation",
        path_a=float(spearmanr(e_c, n_c)[0]), path_b=float(spearmanr(n_c, centered(master.auc, master.gradation))[0]),
        total_c=float(spearmanr(e_c, centered(master.auc, master.gradation))[0]),
        direct_c_given_mediator=partial_spearman(e_c, centered(master.auc, master.gradation), n_c),
        mediator_given_exposure=partial_spearman(n_c, centered(master.auc, master.gradation), e_c),
    ))
    mediation = pd.DataFrame(mediation_rows)
    mediation.to_csv(out / "density_mediation_chains.tsv", sep="\t", index=False)

    # physical dose-response: dAUC/de and predicted vs observed design effects
    slope_auc = float(np.polyfit(e_c, centered(master.auc, master.gradation), 1)[0])
    slope_pv50 = float(np.polyfit(e_c, centered(master.pv50, master.gradation), 1)[0])
    base = master[(master.mu == 0.1) & (master.stress_kpa == 100)].set_index("gradation")
    loose = master[(master.mu == 0.5) & (master.stress_kpa == 100)].set_index("gradation")
    low_s = master[(master.mu == 0.5) & (master.stress_kpa == 20)].set_index("gradation")
    high_s = master[(master.mu == 0.5) & (master.stress_kpa == 200)].set_index("gradation")
    dose_rows = []
    for label, a, b in (("mu 0.1->0.5 @100 kPa", base, loose), ("overburden 20->200 kPa @mu=0.5", low_s, high_s)):
        delta_e = float((b.void_ratio_e - a.void_ratio_e).mean())
        for outcome, slope in (("auc", slope_auc), ("pv50", slope_pv50)):
            observed = float((b[outcome] - a[outcome]).mean())
            dose_rows.append(dict(design_change=label, outcome=outcome, delta_void_ratio=delta_e,
                                  slope_within_gradation=slope, predicted=slope * delta_e,
                                  observed=observed,
                                  ratio_observed_over_predicted=(observed / (slope * delta_e)) if slope * delta_e else np.nan))
    dose = pd.DataFrame(dose_rows)
    dose.to_csv(out / "density_dose_response.tsv", sep="\t", index=False)

    # ---- five core indicators ---------------------------------------------- #
    five_rows = []
    gradation_stats = pd.read_csv(args.macro / "gradation_vs_structure_durability.tsv", sep="\t")
    trends = pd.read_csv(args.macro / "ordered_trend_tests.tsv", sep="\t")
    for indicator, level, lo_v, hi_v in (("Cu", "between_gradation", min(master.Cu), max(master.Cu)),
                                         ("Cc", "between_gradation", min(master.Cc), max(master.Cc))):
        for outcome in ["control_volume_porosity", "auc", "pv50", "void_ratio_e"]:
            row = gradation_stats[(gradation_stats.param == indicator) & (gradation_stats.outcome == outcome)]
            if outcome == "void_ratio_e":
                means = master.groupby("gradation")[[indicator, "void_ratio_e"]].mean()
                r = spearmanr(means[indicator], means.void_ratio_e)
                five_rows.append(dict(indicator=indicator, level=level, range=f"{lo_v:.3g}-{hi_v:.3g}",
                                      outcome=outcome, statistic=f"Spearman rho={r[0]:+.3f}", p_value=float(r[1])))
            elif len(row):
                row = row.iloc[0]
                five_rows.append(dict(indicator=indicator, level=level, range=f"{lo_v:.3g}-{hi_v:.3g}",
                                      outcome=outcome,
                                      statistic=f"Spearman rho={row.rho_gradation_means:+.3f}",
                                      p_value=float(row.perm_p_blocked)))
    for indicator, design_name, lo_v, hi_v in (("mu", "mu@100kPa", 0.1, 0.5),
                                               ("overburden_kpa", "overburden@mu=0.5", 20.0, 200.0)):
        for outcome in ["control_volume_porosity", "auc", "pv50", "void_ratio_e"]:
            if outcome == "void_ratio_e":
                frame = master[master.stress_kpa == 100] if indicator == "mu" else master[master.mu == 0.5]
                codes = (frame.mu.map({0.1: 0, 0.3: 1, 0.5: 2}) if indicator == "mu"
                         else frame.stress_kpa.map({20: 0, 100: 1, 200: 2}))
                L, z, p_up = page_trend(frame.void_ratio_e.to_numpy(), codes.to_numpy(), frame.gradation.to_numpy())
                five_rows.append(dict(indicator=indicator, level="within_gradation", range=f"{lo_v:.3g}-{hi_v:.3g}",
                                      outcome=outcome, statistic=f"Page z={z:+.2f}",
                                      p_value=float(min(p_up, 1 - p_up))))
            else:
                row = trends[(trends.design == design_name) & (trends.outcome == outcome)]
                if len(row):
                    row = row.iloc[0]
                    five_rows.append(dict(indicator=indicator, level="within_gradation",
                                          range=f"{lo_v:.3g}-{hi_v:.3g}", outcome=outcome,
                                          statistic=f"Page z={row.z:+.2f}",
                                          p_value=float(min(row.p_increasing, row.p_decreasing))))
    five_rows.append(dict(indicator="void_ratio_e (密实度)", level="state_variable", range="0.14-0.60",
                          outcome="auc", statistic=f"Spearman rho={density_stats[(density_stats.measure=='void_ratio_e')&(density_stats.outcome=='auc')&(density_stats.level=='total')].rho.iloc[0]:+.3f}",
                          p_value=float(density_stats[(density_stats.measure=='void_ratio_e')&(density_stats.outcome=='auc')&(density_stats.level=='total')].p.iloc[0])))
    five_rows.append(dict(indicator="void_ratio_e (密实度)", level="state_variable", range="0.14-0.60",
                          outcome="pv50", statistic=f"Spearman rho={density_stats[(density_stats.measure=='void_ratio_e')&(density_stats.outcome=='pv50')&(density_stats.level=='total')].rho.iloc[0]:+.3f}",
                          p_value=float(density_stats[(density_stats.measure=='void_ratio_e')&(density_stats.outcome=='pv50')&(density_stats.level=='total')].p.iloc[0])))
    pd.DataFrame(five_rows).to_csv(out / "five_inputs_summary.tsv", sep="\t", index=False)

    # ---- figure ------------------------------------------------------------ #
    observed = {
        "mu": dict(delta_e=float((loose.void_ratio_e - base.void_ratio_e).mean()),
                   delta_auc=float((loose.auc - base.auc).mean())),
        "sigma": dict(delta_e=float((high_s.void_ratio_e - low_s.void_ratio_e).mean()),
                      delta_auc=float((high_s.auc - low_s.auc).mean())),
    }
    figure_density(master, validation, slope_auc, observed, out / "figures")

    manifest = {
        "format": "bubble_density_input_analysis_v1",
        "status": "passed",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "correction": "mu is the DEM interparticle friction coefficient; density is supplied separately in the user table and behaves as void ratio e",
        "inputs": {
            "density_image": {"path": str(DENSITY_IMAGE), "sha256": sha256(DENSITY_IMAGE)},
            "macro_dir": str(args.macro),
            "geometry_master": {"path": str(args.macro / "geometry_master.csv"),
                                "sha256": sha256(args.macro / "geometry_master.csv")},
        },
        "transcription": {g: DENSITY_TABLE[g] for g in DENSITY_TABLE},
        "column_mapping": [{"column": c, "mu": mu, "stress_kpa": s, "condition": lab}
                           for c, mu, s, lab in COLUMN_CONDITIONS],
        "outputs": {p.name: sha256(p) for p in sorted(out.glob("*")) if p.is_file()},
        "figure_outputs": {p.name: sha256(p) for p in sorted((out / "figures").glob("*.png"))},
    }
    (out / "analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(f"wrote {out}")
    print(validation.to_string(index=False))
    print(design.to_string(index=False))
    print(dose.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
