#!/usr/bin/env python3
"""
Stage-1 peer-review revision analyses
=====================================
Implements minimum corrections with existing data:

- Primary cohort: 18 parent-confirmed wards (unstarred; n=126)
- Sensitivity: all 29 investigated wards (n=162)
- First-detection order taxonomy (NOT "dual seeding")
- Separate inpatient vs staff/outpatient numerators
- Unrepresented HA clusters coded missing (not zero rates)
- Pooled NT/urban ratios + means of cluster rates
- Exact Fisher OR with CI; underpowered regional comparison
- Leave-one-hospital / leave-one-ward-out; exclude BA.5.6 ward
- Tie windows: same day, ±1 day, ±2 days
- Cross-tab vs parent confirmation status
- Contribution matrix Gu et al. vs this secondary note

Does NOT claim: dual seeding mechanism, true incidence excess,
density paradox as causal, visitation effect, or psychiatric causation.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
HAI_PATH = ROOT / "data/hospital_data/metadata_cases_plot.csv"
CAPACITY_PATH = ROOT / "analysis/reference/ha_cluster_capacity_2022.csv"
THROUGHPUT_PATH = ROOT / "analysis/reference/ha_cluster_throughput_2021_22.csv"
ASCERTAIN_PATH = ROOT / "analysis/reference/ascertainment_presence_by_cluster.csv"
OUT = ROOT / "analysis/outputs"
FIG = ROOT / "manuscript/figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

NT = {"NTWC", "NTEC"}
URBAN_IN_SAMPLE = {"KWC", "KEC", "HKEC"}
ALL_CLUSTERS = ["NTWC", "NTEC", "KWC", "HKEC", "KEC", "KCC", "HKWC"]


def load_hai() -> pd.DataFrame:
    hai = pd.read_csv(HAI_PATH, parse_dates=["Date"])
    hai["starred"] = hai["Wards"].astype(str).str.endswith("*")
    hai["confirmed"] = ~hai["starred"]
    hai["case_group"] = np.where(
        hai["Case type"].eq("Outpatient/Staff"), "staff_outpatient", "inpatient"
    )
    hai["region"] = np.where(
        hai["Hospital clusters"].isin(NT), "NT", "Urban_in_sample"
    )
    return hai


def expand_cases(hai: pd.DataFrame) -> pd.DataFrame:
    """Expand aggregated N rows to one row per sequenced person (for dates)."""
    rows = []
    for _, r in hai.iterrows():
        for _ in range(int(r["N"])):
            rows.append(r.to_dict())
    return pd.DataFrame(rows)


def classify_pathways(hai: pd.DataFrame, tie_days: int = 0) -> pd.DataFrame:
    """
    First-detection order by ward.
    tie_days=0: same calendar day = same_day
    tie_days=1/2: |staff_date - inp_date| <= tie_days treated as same_window
    """
    rows = []
    for (h, w), g in hai.groupby(["Hospital_anonymized", "Wards"]):
        staff_dates = g.loc[g["case_group"] == "staff_outpatient", "Date"]
        inp_dates = g.loc[g["case_group"] == "inpatient", "Date"]
        staff = staff_dates.min() if len(staff_dates) else pd.NaT
        inp = inp_dates.min() if len(inp_dates) else pd.NaT
        if pd.isna(staff) and pd.isna(inp):
            path = "unknown"
        elif pd.isna(staff):
            path = "inpatient_only"
        elif pd.isna(inp):
            path = "staff_outpatient_only"
        else:
            delta = abs((staff - inp).days)
            if delta <= tie_days:
                path = "same_window" if tie_days > 0 else "same_day"
            elif staff < inp:
                path = "staff_outpatient_before_inpatient"
            else:
                path = "inpatient_before_staff_outpatient"
        rows.append(
            {
                "Hospital_anonymized": h,
                "Wards": w,
                "ha_cluster": g["Hospital clusters"].iloc[0],
                "region": g["region"].iloc[0],
                "confirmed": bool(g["confirmed"].iloc[0]),
                "path": path,
                "cases": int(g["N"].sum()),
                "n_inpatient": int(g.loc[g["case_group"] == "inpatient", "N"].sum()),
                "n_staff_outpatient": int(
                    g.loc[g["case_group"] == "staff_outpatient", "N"].sum()
                ),
                "lineages": ",".join(sorted(g["Lineages"].unique())),
                "is_ba56_ward": "BA.5.6" in set(g["Lineages"]),
            }
        )
    return pd.DataFrame(rows)


def fisher_with_ci(a: int, b: int, c: int, d: int):
    """2x2 Fisher exact; OR CI via conditional MLE approx (statsmodels if available)."""
    table = np.array([[a, b], [c, d]])
    oddsratio, p = stats.fisher_exact(table)
    # Exact CI via scipy hypergeom / conditional; use statsmodels if present
    try:
        from statsmodels.stats.contingency_tables import Table2x2

        t = Table2x2(table)
        ci_low, ci_high = t.oddsratio_confint(alpha=0.05, method="exact")
        or_mle = float(t.oddsratio)
    except Exception:
        # Woolf logit CI fallback (may be infinite with zeros)
        if min(a, b, c, d) == 0:
            ci_low, ci_high = np.nan, np.nan
            or_mle = oddsratio
        else:
            log_or = np.log((a * d) / (b * c))
            se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
            ci_low, ci_high = np.exp(log_or - 1.96 * se), np.exp(log_or + 1.96 * se)
            or_mle = float(np.exp(log_or))
    return {
        "a_nt_staff_before": a,
        "b_nt_inp_before": b,
        "c_urban_staff_before": c,
        "d_urban_inp_before": d,
        "odds_ratio": float(or_mle) if np.isfinite(or_mle) else float(oddsratio),
        "p_exact": float(p),
        "or_ci95_low": float(ci_low) if ci_low == ci_low else None,
        "or_ci95_high": float(ci_high) if ci_high == ci_high else None,
        "n_wards": int(a + b + c + d),
        "interpretation": (
            "Severely imprecise; failure to detect a difference is not evidence of "
            "equivalence or a shared regional mechanism."
        ),
    }


def regional_fisher(pathways: pd.DataFrame) -> dict:
    ordered = pathways[
        pathways["path"].isin(
            [
                "staff_outpatient_before_inpatient",
                "inpatient_before_staff_outpatient",
            ]
        )
    ]
    # NT staff_before, NT inp_before, Urban staff_before, Urban inp_before
    a = int(
        (
            (ordered["region"] == "NT")
            & (ordered["path"] == "staff_outpatient_before_inpatient")
        ).sum()
    )
    b = int(
        (
            (ordered["region"] == "NT")
            & (ordered["path"] == "inpatient_before_staff_outpatient")
        ).sum()
    )
    c = int(
        (
            (ordered["region"] == "Urban_in_sample")
            & (ordered["path"] == "staff_outpatient_before_inpatient")
        ).sum()
    )
    d = int(
        (
            (ordered["region"] == "Urban_in_sample")
            & (ordered["path"] == "inpatient_before_staff_outpatient")
        ).sum()
    )
    return fisher_with_ci(a, b, c, d)


def cluster_case_counts(hai: pd.DataFrame) -> pd.DataFrame:
    g = (
        hai.groupby("Hospital clusters")
        .agg(
            sequenced_all=("N", "sum"),
            sequenced_inpatient=("N", lambda s: s[hai.loc[s.index, "case_group"] == "inpatient"].sum()),
            sequenced_staff_outpatient=(
                "N",
                lambda s: s[hai.loc[s.index, "case_group"] == "staff_outpatient"].sum(),
            ),
            n_hospitals=("Hospital_anonymized", "nunique"),
            n_wards=("Wards", "nunique"),
        )
        .reset_index()
        .rename(columns={"Hospital clusters": "ha_cluster"})
    )
    # Fix inpatient/staff with explicit groupby
    rows = []
    for cl, sub in hai.groupby("Hospital clusters"):
        rows.append(
            {
                "ha_cluster": cl,
                "sequenced_all": int(sub["N"].sum()),
                "sequenced_inpatient": int(sub.loc[sub["case_group"] == "inpatient", "N"].sum()),
                "sequenced_staff_outpatient": int(
                    sub.loc[sub["case_group"] == "staff_outpatient", "N"].sum()
                ),
                "n_hospitals": int(sub["Hospital_anonymized"].nunique()),
                "n_wards": int(sub["Wards"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def build_denominator_table(
    counts: pd.DataFrame, capacity: pd.DataFrame, throughput: pd.DataFrame, ascertain: pd.DataFrame
) -> pd.DataFrame:
    base = capacity[["ha_cluster", "beds_total", "catchment_population_2017", "land_area_km2_core",
                     "n_hospitals_institutions", "region_group"]].merge(
        throughput, on="ha_cluster", how="left"
    )
    if ascertain is not None and len(ascertain):
        base = base.merge(
            ascertain[["ha_cluster", "in_official_bulletins", "in_sequenced_table", "mention_days"]],
            on="ha_cluster",
            how="left",
        )
    out = base.merge(counts, on="ha_cluster", how="left")
    for col in [
        "sequenced_all",
        "sequenced_inpatient",
        "sequenced_staff_outpatient",
        "n_hospitals",
        "n_wards",
    ]:
        out[col] = out[col].fillna(0).astype(int)

    # Representation flag: missing genomic sample vs true zero
    out["genomic_representation"] = np.where(
        out["sequenced_all"] > 0,
        "represented",
        "not_represented_in_genomic_sample",
    )

    # Rates only for represented clusters; NA otherwise
    def rate(num, den, scale):
        return np.where(
            out["genomic_representation"] == "represented",
            num / den * scale,
            np.nan,
        )

    out["inpatient_per_1000_beds"] = rate(out["sequenced_inpatient"], out["beds_total"], 1000)
    out["inpatient_per_100k_discharges"] = rate(
        out["sequenced_inpatient"], out["ip_discharges_2021_22"], 1e5
    )
    out["inpatient_per_100k_patient_days"] = rate(
        out["sequenced_inpatient"], out["patient_days_2021_22"], 1e5
    )
    # Descriptive only — mixed numerator; not an HAI incidence rate
    out["all_linked_per_1000_beds_DESCRIPTIVE"] = rate(
        out["sequenced_all"], out["beds_total"], 1000
    )
    out["pop_density"] = out["catchment_population_2017"] / out["land_area_km2_core"]
    out["hospitals_per_100km2"] = (
        out["n_hospitals_institutions"] / out["land_area_km2_core"] * 100
    )
    out["psych_patient_day_share"] = (
        out["patient_days_mentally_ill"] / out["patient_days_2021_22"]
    )
    return out


def nt_urban_ratios(den: pd.DataFrame, num_col: str, rate_col: str) -> dict:
    """Pooled ratio and mean-of-rates ratio for represented NT vs urban-in-sample."""
    d = den[den["genomic_representation"] == "represented"].copy()
    nt = d[d["ha_cluster"].isin(NT)]
    urb = d[d["ha_cluster"].isin(URBAN_IN_SAMPLE)]
    # Pooled: sum(num)/sum(den) style via rate reconstruction
    # Use beds for inpatient_per_1000_beds etc.
    if rate_col == "inpatient_per_1000_beds":
        den_col = "beds_total"
        scale = 1000
    elif rate_col == "inpatient_per_100k_discharges":
        den_col = "ip_discharges_2021_22"
        scale = 1e5
    elif rate_col == "inpatient_per_100k_patient_days":
        den_col = "patient_days_2021_22"
        scale = 1e5
    else:
        den_col = "beds_total"
        scale = 1000

    pooled_nt = nt[num_col].sum() / nt[den_col].sum() * scale
    pooled_urb = urb[num_col].sum() / urb[den_col].sum() * scale
    mean_nt = nt[rate_col].mean()
    mean_urb = urb[rate_col].mean()
    return {
        "numerator": num_col,
        "rate": rate_col,
        "pooled_nt": float(pooled_nt),
        "pooled_urban": float(pooled_urb),
        "pooled_ratio_nt_urban": float(pooled_nt / pooled_urb) if pooled_urb else None,
        "mean_rate_nt": float(mean_nt),
        "mean_rate_urban": float(mean_urb),
        "mean_ratio_nt_urban": float(mean_nt / mean_urb) if mean_urb else None,
        "caveat": (
            "Ratios describe concentration within the selected genomic sample after "
            "crude patient-throughput scaling; they are not selection-adjusted "
            "incidence ratios and may be explained by differential ascertainment."
        ),
    }


def leave_one_out(hai: pd.DataFrame, capacity, throughput, ascertain, by: str) -> pd.DataFrame:
    """by in {'hospital','ward'}."""
    keys = (
        hai["Hospital_anonymized"].unique()
        if by == "hospital"
        else list(hai[["Hospital_anonymized", "Wards"]].drop_duplicates().itertuples(index=False))
    )
    rows = []
    for key in keys:
        if by == "hospital":
            sub = hai[hai["Hospital_anonymized"] != key]
            label = str(key)
        else:
            h, w = key
            sub = hai[~((hai["Hospital_anonymized"] == h) & (hai["Wards"] == w))]
            label = f"{h}:{w}"
        counts = cluster_case_counts(sub)
        den = build_denominator_table(counts, capacity, throughput, ascertain)
        ratios = nt_urban_ratios(den, "sequenced_inpatient", "inpatient_per_1000_beds")
        rows.append(
            {
                "excluded": label,
                "remaining_cases": int(sub["N"].sum()),
                "pooled_ratio_nt_urban_inpatient_per_1000_beds": ratios["pooled_ratio_nt_urban"],
                "mean_ratio_nt_urban_inpatient_per_1000_beds": ratios["mean_ratio_nt_urban"],
            }
        )
    return pd.DataFrame(rows)


def pathway_summary(pathways: pd.DataFrame) -> pd.DataFrame:
    return (
        pathways.groupby("path")
        .agg(wards=("Wards", "count"), cases=("cases", "sum"),
             inpatient=("n_inpatient", "sum"), staff_outpatient=("n_staff_outpatient", "sum"))
        .reset_index()
        .sort_values("cases", ascending=False)
    )


def make_figures(den_primary: pd.DataFrame, pathways_primary: pd.DataFrame, hai_primary: pd.DataFrame):
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 150})

    # Fig 1: first-detection composition
    order = [
        "staff_outpatient_before_inpatient",
        "inpatient_before_staff_outpatient",
        "same_day",
        "staff_outpatient_only",
        "inpatient_only",
    ]
    labels = {
        "staff_outpatient_before_inpatient": "Staff/outpatient before inpatient",
        "inpatient_before_staff_outpatient": "Inpatient before staff/outpatient",
        "same_day": "Same-day first detection",
        "staff_outpatient_only": "Staff/outpatient only",
        "inpatient_only": "Inpatient only",
    }
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    regions = ["NT", "Urban_in_sample"]
    x = np.arange(len(order))
    width = 0.38
    for i, reg in enumerate(regions):
        vals = [
            pathways_primary.loc[
                (pathways_primary["region"] == reg) & (pathways_primary["path"] == p), "cases"
            ].sum()
            for p in order
        ]
        ax.bar(x + (i - 0.5) * width, vals, width, label=reg.replace("_", " "))
    ax.set_xticks(x)
    ax.set_xticklabels([labels[p] for p in order], rotation=25, ha="right")
    ax.set_ylabel("Sequenced persons (confirmed wards)")
    ax.set_title("First-detection order by region (18 confirmed wards)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "fig_first_detection_by_region.pdf")
    fig.savefig(OUT / "fig_first_detection_by_region.pdf")
    plt.close(fig)

    # Fig 2: inpatient rates — represented only
    d = den_primary[den_primary["genomic_representation"] == "represented"].copy()
    d = d.sort_values("inpatient_per_1000_beds", ascending=False)
    colors = ["#2a6f4e" if c in NT else "#4a6fa5" for c in d["ha_cluster"]]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    metrics = [
        ("inpatient_per_1000_beds", "Per 1,000 beds"),
        ("inpatient_per_100k_discharges", "Per 100k discharges"),
        ("inpatient_per_100k_patient_days", "Per 100k patient-days"),
    ]
    for ax, (col, title) in zip(axes, metrics):
        ax.bar(d["ha_cluster"], d[col], color=colors)
        ax.set_title(title)
        ax.set_ylabel("Sequenced inpatients (confirmed)")
        ax.tick_params(axis="x", rotation=0)
    fig.suptitle(
        "Sampled confirmed-ward inpatients after crude patient-throughput scaling\n"
        "(unrepresented clusters omitted — not coded as zero)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(FIG / "fig_inpatient_throughput_scaling.pdf")
    fig.savefig(OUT / "fig_inpatient_throughput_scaling.pdf")
    plt.close(fig)

    # Fig 3: density discordance (descriptive)
    d2 = den_primary[den_primary["genomic_representation"] == "represented"].copy()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for _, r in d2.iterrows():
        color = "#2a6f4e" if r["ha_cluster"] in NT else "#4a6fa5"
        ax.scatter(
            r["pop_density"],
            r["inpatient_per_1000_beds"],
            s=80 + 4 * r["sequenced_inpatient"],
            c=color,
            edgecolors="k",
            linewidths=0.4,
            zorder=3,
        )
        ax.annotate(r["ha_cluster"], (r["pop_density"], r["inpatient_per_1000_beds"]),
                    textcoords="offset points", xytext=(5, 4), fontsize=9)
    ax.set_xlabel("Catchment population density (persons/km², planning catchment)")
    ax.set_ylabel("Sequenced confirmed-ward inpatients / 1,000 beds")
    ax.set_title("Discordance: sampled inpatient concentration vs coarse density")
    fig.tight_layout()
    fig.savefig(FIG / "fig_density_discordance.pdf")
    fig.savefig(OUT / "fig_density_discordance.pdf")
    plt.close(fig)

    # Fig 4: ascertainment presence
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for _, r in den_primary.iterrows():
        x = r.get("mention_days", 0) or 0
        y = r["sequenced_all"]
        color = "#2a6f4e" if r["ha_cluster"] in NT else "#4a6fa5"
        marker = "o" if r["genomic_representation"] == "represented" else "x"
        ax.scatter(x, y, c=color, marker=marker, s=70, zorder=3)
        ax.annotate(r["ha_cluster"], (x, y), textcoords="offset points", xytext=(5, 3), fontsize=8)
    ax.set_xlabel("Official nosocomial bulletin mention-days (presence proxy)")
    ax.set_ylabel("Sequenced persons in confirmed wards (0 = not in sample)")
    ax.set_title("Geographic genomic representation vs official bulletin footprint")
    fig.tight_layout()
    fig.savefig(FIG / "fig_ascertainment_presence.pdf")
    fig.savefig(OUT / "fig_ascertainment_presence.pdf")
    plt.close(fig)


def contribution_matrix() -> list[dict]:
    return [
        {
            "item": "162 genomes / 29 suspected / 18 confirmed / 126 confirmed genomes",
            "source": "Gu et al. 2026",
            "this_note": "Reproduced as cohort definition",
        },
        {
            "item": "51.9% coverage of official nosocomial cases",
            "source": "Gu et al. 2026",
            "this_note": "Cited; not re-estimated",
        },
        {
            "item": "2/18 phylogenetically supported staff-before-inpatient clusters",
            "source": "Gu et al. 2026",
            "this_note": "Cited; temporal taxonomy is complementary, not equivalent",
        },
        {
            "item": "Mobility GAM null after time smooth",
            "source": "Gu et al. 2026",
            "this_note": "Cited as background",
        },
        {
            "item": "API / BA.5.6 single-ward deep outbreak description",
            "source": "Gu et al. 2026 (genomic)",
            "this_note": "Sensitivity exclusion of BA.5.6 ward only",
        },
        {
            "item": "Ward-level first-detection order taxonomy (all confirmed wards)",
            "source": "—",
            "this_note": "New descriptive analysis",
        },
        {
            "item": "HA-cluster sample concentration with inpatient-only throughput scaling",
            "source": "—",
            "this_note": "New; rates are sample concentration, not incidence",
        },
        {
            "item": "Official bulletin presence vs genomic representation audit",
            "source": "—",
            "this_note": "New ascertainment-aware geographic audit",
        },
        {
            "item": "Leave-one-out / BA.5.6 exclusion / tie-window sensitivity",
            "source": "—",
            "this_note": "New robustness checks",
        },
    ]


def main():
    hai = load_hai()
    capacity = pd.read_csv(CAPACITY_PATH)
    throughput = pd.read_csv(THROUGHPUT_PATH)
    ascertain = pd.read_csv(ASCERTAIN_PATH) if ASCERTAIN_PATH.exists() else pd.DataFrame()

    hai_primary = hai[hai["confirmed"]].copy()
    hai_all = hai.copy()

    # Pathways
    pw_primary = classify_pathways(hai_primary, tie_days=0)
    pw_all = classify_pathways(hai_all, tie_days=0)
    pw_tie1 = classify_pathways(hai_primary, tie_days=1)
    pw_tie2 = classify_pathways(hai_primary, tie_days=2)

    # Cross-tab confirmation × path (all 29)
    xtab = (
        pw_all.groupby(["confirmed", "path"])
        .agg(wards=("Wards", "count"), cases=("cases", "sum"))
        .reset_index()
    )

    fisher_primary = regional_fisher(pw_primary)
    fisher_all = regional_fisher(pw_all)

    counts_primary = cluster_case_counts(hai_primary)
    counts_all = cluster_case_counts(hai_all)
    den_primary = build_denominator_table(counts_primary, capacity, throughput, ascertain)
    den_all = build_denominator_table(counts_all, capacity, throughput, ascertain)

    ratios_primary = [
        nt_urban_ratios(den_primary, "sequenced_inpatient", "inpatient_per_1000_beds"),
        nt_urban_ratios(den_primary, "sequenced_inpatient", "inpatient_per_100k_discharges"),
        nt_urban_ratios(den_primary, "sequenced_inpatient", "inpatient_per_100k_patient_days"),
    ]

    # Exclude BA.5.6 ward (hospital I, A102)
    hai_no_ba56 = hai_primary[
        ~((hai_primary["Hospital_anonymized"] == "I") & (hai_primary["Wards"] == "A102"))
    ]
    den_no_ba56 = build_denominator_table(
        cluster_case_counts(hai_no_ba56), capacity, throughput, ascertain
    )
    ratios_no_ba56 = nt_urban_ratios(
        den_no_ba56, "sequenced_inpatient", "inpatient_per_1000_beds"
    )

    loo_hosp = leave_one_out(hai_primary, capacity, throughput, ascertain, "hospital")
    loo_ward = leave_one_out(hai_primary, capacity, throughput, ascertain, "ward")

    # Case-type totals
    case_totals = {
        "primary_confirmed_n": int(hai_primary["N"].sum()),
        "primary_wards": int(hai_primary[["Hospital_anonymized", "Wards"]].drop_duplicates().shape[0]),
        "sensitivity_all_n": int(hai_all["N"].sum()),
        "sensitivity_wards": int(hai_all[["Hospital_anonymized", "Wards"]].drop_duplicates().shape[0]),
        "primary_inpatient": int(hai_primary.loc[hai_primary["case_group"] == "inpatient", "N"].sum()),
        "primary_staff_outpatient": int(
            hai_primary.loc[hai_primary["case_group"] == "staff_outpatient", "N"].sum()
        ),
        "primary_api_ge3": int(
            hai_primary.loc[hai_primary["Case type"] == "Inpatient (API >= 3)", "N"].sum()
        ),
        "nt_share_primary_all": float(
            hai_primary.loc[hai_primary["Hospital clusters"].isin(NT), "N"].sum()
            / hai_primary["N"].sum()
        ),
        "nt_share_primary_inpatient": float(
            hai_primary.loc[
                (hai_primary["Hospital clusters"].isin(NT))
                & (hai_primary["case_group"] == "inpatient"),
                "N",
            ].sum()
            / max(hai_primary.loc[hai_primary["case_group"] == "inpatient", "N"].sum(), 1)
        ),
    }

    # Density discordance summary (descriptive)
    drep = den_primary[den_primary["genomic_representation"] == "represented"]
    dens = {
        "mean_pop_density_nt_over_urban": float(
            drep.loc[drep["ha_cluster"].isin(NT), "pop_density"].mean()
            / drep.loc[drep["ha_cluster"].isin(URBAN_IN_SAMPLE), "pop_density"].mean()
        ),
        "mean_hosp_per_100km2_nt_over_urban": float(
            drep.loc[drep["ha_cluster"].isin(NT), "hospitals_per_100km2"].mean()
            / drep.loc[drep["ha_cluster"].isin(URBAN_IN_SAMPLE), "hospitals_per_100km2"].mean()
        ),
        "note": (
            "Coarse catchment density measures; administrative catchments are not closed "
            "cohorts; Islands/Lantau handling may inflate urban density. Discordance is "
            "descriptive of the sample, not a proven paradox of incidence."
        ),
    }

    summary = {
        "case_totals": case_totals,
        "fisher_primary_confirmed": fisher_primary,
        "fisher_all_investigated": fisher_all,
        "nt_urban_ratios_primary_inpatient": ratios_primary,
        "nt_urban_ratio_excluding_ba56_ward": ratios_no_ba56,
        "density_discordance": dens,
        "loo_hospital_ratio_range": {
            "min": float(loo_hosp["pooled_ratio_nt_urban_inpatient_per_1000_beds"].min()),
            "max": float(loo_hosp["pooled_ratio_nt_urban_inpatient_per_1000_beds"].max()),
            "median": float(loo_hosp["pooled_ratio_nt_urban_inpatient_per_1000_beds"].median()),
        },
        "loo_ward_ratio_range": {
            "min": float(loo_ward["pooled_ratio_nt_urban_inpatient_per_1000_beds"].min()),
            "max": float(loo_ward["pooled_ratio_nt_urban_inpatient_per_1000_beds"].max()),
            "median": float(loo_ward["pooled_ratio_nt_urban_inpatient_per_1000_beds"].median()),
        },
        "contribution_matrix": contribution_matrix(),
        "claims_forbidden": [
            "dual seeding as identified mechanism",
            "true NT nosocomial incidence excess",
            "density paradox as causal finding",
            "geographic equivalence of routes",
            "causal role of psychiatric case mix",
            "effect of visitation relaxation",
            "KEC as negative control",
            "KCC/HKWC as zero HAI burden",
        ],
    }

    # Write tables
    pw_primary.to_csv(OUT / "stage1_pathways_confirmed.csv", index=False)
    pw_all.to_csv(OUT / "stage1_pathways_all29.csv", index=False)
    pathway_summary(pw_primary).to_csv(OUT / "stage1_pathway_summary_confirmed.csv", index=False)
    pathway_summary(pw_tie1).to_csv(OUT / "stage1_pathway_summary_tie1.csv", index=False)
    pathway_summary(pw_tie2).to_csv(OUT / "stage1_pathway_summary_tie2.csv", index=False)
    xtab.to_csv(OUT / "stage1_path_by_confirmation.csv", index=False)
    den_primary.to_csv(OUT / "stage1_denominator_confirmed.csv", index=False)
    den_all.to_csv(OUT / "stage1_denominator_all29.csv", index=False)
    loo_hosp.to_csv(OUT / "stage1_loo_hospital.csv", index=False)
    loo_ward.to_csv(OUT / "stage1_loo_ward.csv", index=False)
    pd.DataFrame(contribution_matrix()).to_csv(OUT / "stage1_contribution_matrix.csv", index=False)
    with open(OUT / "stage1_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    make_figures(den_primary, pw_primary, hai_primary)

    # Decision log (outside manuscript)
    decision = OUT / "DECISION_LOG_stage1.md"
    decision.write_text(
        """# Decision log (Stage-1 peer-review revision)

Moved out of the manuscript body per peer review.

1. Parent paper attributed to **Gu et al. (2026)**, DOI 10.1111/irv.70249 (PMC13097321).
2. Primary cohort = 18 confirmed (unstarred) wards, n=126; sensitivity = all 29, n=162.
3. Language: first-detection order, not dual seeding; sample concentration, not incidence excess.
4. Inpatient numerators for patient-throughput rates; staff/outpatient reported separately.
5. Unrepresented clusters (KCC, HKWC) marked not represented — never zero rates.
6. KEC = descriptive comparator, not negative control.
7. Psychiatric case-mix = candidate explanation only; BA.5.6 leave-one-ward sensitivity run.
8. Visitation appears only as study-period context, not a causal claim.
9. Named specialty hospitals omitted from manuscript text to reduce re-identification risk.
10. Full official sampling fractions remain unavailable; selection-adjusted inference deferred.
""",
        encoding="utf-8",
    )

    print(json.dumps({
        "primary_n": case_totals["primary_confirmed_n"],
        "fisher": fisher_primary,
        "pooled_ratio_beds": ratios_primary[0]["pooled_ratio_nt_urban"],
        "mean_ratio_beds": ratios_primary[0]["mean_ratio_nt_urban"],
        "excl_ba56_pooled": ratios_no_ba56["pooled_ratio_nt_urban"],
        "pathway_summary": pathway_summary(pw_primary).to_dict(orient="records"),
    }, indent=2, default=str))


if __name__ == "__main__":
    main()
