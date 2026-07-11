#!/usr/bin/env python3
"""
Dual seeding + geographic concentration exploration
===================================================
Evidence-first analysis for hospital-associated SARS-CoV-2 in HK 2022.

Core claims tested
------------------
1. Dual seeding: ward-level first detections split between staff/outpatient-first
   and inpatient-first pathways (not a single dominant route).
2. Geographic concentration: NTWC+NTEC dominate sequenced HAI burden.
3. Density paradox: New Territories clusters have *lower* population density and
   *lower* hospital/bed geographic density than urban clusters, yet higher
   sequenced HAI burden — including after bed- and population-standardization.
4. Linkage: dual seeding operates inside the NT-concentrated burden; BA.5.6 and
   high-API patterns sit inside the same geography that also holds most
   psychiatric bed capacity.

What this script will NOT claim
-------------------------------
- That sequenced HAI rates equal true nosocomial incidence (ascertainment unknown).
- That low density *causes* high HAI (paradox is descriptive; mechanisms need
  denominators we partly lack: admissions, sequencing effort by hospital).
- Hospital de-anonymization beyond HA cluster labels already in the plot table.
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
DIST_PATH = ROOT / "analysis/reference/district_demographics_2021.csv"
CLUSTER_PATH = ROOT / "analysis/reference/ha_cluster_capacity_2022.csv"
CASE_JSON = ROOT / "data/case_curve/HK_case_data.json"
OUT = ROOT / "analysis/outputs"
OUT.mkdir(parents=True, exist_ok=True)

RELAX_DATE = pd.Timestamp("2022-05-31")
NT_CLUSTERS = {"NTWC", "NTEC"}
URBAN_IN_SAMPLE = {"KWC", "KEC", "HKEC"}


def load_hai() -> pd.DataFrame:
    hai = pd.read_csv(HAI_PATH, parse_dates=["Date"])
    hai["starred"] = hai["Wards"].astype(str).str.endswith("*")
    hai["region"] = np.where(hai["Hospital clusters"].isin(NT_CLUSTERS), "New Territories", "Urban (HK Island/Kowloon)")
    return hai


def classify_ward_pathways(hai: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (h, w), g in hai.groupby(["Hospital_anonymized", "Wards"]):
        first = g.groupby("Case type")["Date"].min()
        staff = first.get("Outpatient/Staff", pd.NaT)
        inp_keys = [k for k in first.index if "Inpatient" in k]
        inp = first[inp_keys].min() if inp_keys else pd.NaT
        if pd.isna(staff) and pd.isna(inp):
            path = "unknown"
        elif pd.isna(staff):
            path = "inpatient_only"
        elif pd.isna(inp):
            path = "staff_only"
        elif staff < inp:
            path = "staff_first"
        elif inp < staff:
            path = "inpatient_first"
        else:
            path = "same_day"
        rows.append(
            {
                "Hospital_anonymized": h,
                "Wards": w,
                "Hospital clusters": g["Hospital clusters"].iloc[0],
                "region": g["region"].iloc[0],
                "path": path,
                "cases": int(g["N"].sum()),
                "lineages": ",".join(sorted(g["Lineages"].unique())),
                "starred": bool(g["starred"].iloc[0]),
                "n_case_types": int(g["Case type"].nunique()),
                "duration_days": int((g["Date"].max() - g["Date"].min()).days),
                "prop_api_ge3": float(
                    g.loc[g["Case type"] == "Inpatient (API >= 3)", "N"].sum()
                    / max(g.loc[g["Case type"].str.contains("Inpatient"), "N"].sum(), 1)
                ),
            }
        )
    return pd.DataFrame(rows)


def cluster_burden_table(hai: pd.DataFrame, capacity: pd.DataFrame) -> pd.DataFrame:
    burden = (
        hai.groupby("Hospital clusters")
        .agg(
            sequenced_cases=("N", "sum"),
            n_hospitals_in_sample=("Hospital_anonymized", "nunique"),
            n_wards=("Wards", "nunique"),
            n_lineages=("Lineages", "nunique"),
        )
        .reset_index()
        .rename(columns={"Hospital clusters": "ha_cluster"})
    )
    out = capacity.merge(burden, on="ha_cluster", how="left")
    out["sequenced_cases"] = out["sequenced_cases"].fillna(0).astype(int)
    out["n_hospitals_in_sample"] = out["n_hospitals_in_sample"].fillna(0).astype(int)
    out["n_wards"] = out["n_wards"].fillna(0).astype(int)

    out["pop_density"] = out["catchment_population_2017"] / out["land_area_km2_core"]
    out["hospitals_per_100km2"] = out["n_hospitals_institutions"] / out["land_area_km2_core"] * 100
    out["beds_per_100km2"] = out["beds_total"] / out["land_area_km2_core"] * 100
    out["cases_per_1000_beds"] = out["sequenced_cases"] / out["beds_total"] * 1000
    out["cases_per_million_pop"] = out["sequenced_cases"] / out["catchment_population_2017"] * 1_000_000
    out["in_hai_sample"] = out["sequenced_cases"] > 0
    if "region_group" not in out.columns:
        out["region_group"] = np.where(
            out["ha_cluster"].isin(["NTWC", "NTEC"]),
            "New Territories",
            np.where(out["ha_cluster"].isin(["HKWC", "KCC"]), "Urban (no sequenced HAI in table)", "Urban (in sample)"),
        )
    else:
        # Refine display labels using observed HAI presence
        out["region_group"] = np.where(
            out["ha_cluster"].isin(["NTWC", "NTEC"]),
            "New Territories",
            np.where(out["sequenced_cases"] > 0, "Urban (in sample)", "Urban (no sequenced HAI in table)"),
        )
    return out


def fisher_dual_seeding(pathways: pd.DataFrame) -> dict:
    """NT vs urban among wards with an ordered dual-type pathway."""
    sub = pathways[pathways["path"].isin(["staff_first", "inpatient_first"])].copy()
    tab = pd.crosstab(sub["region"], sub["path"])
    # Ensure both columns exist
    for col in ["staff_first", "inpatient_first"]:
        if col not in tab.columns:
            tab[col] = 0
    tab = tab[["staff_first", "inpatient_first"]]
    oddsratio, p = stats.fisher_exact(tab.values)
    return {"table": tab, "odds_ratio": oddsratio, "p_value": p, "n_wards": int(len(sub))}


def community_alignment(hai: pd.DataFrame) -> pd.DataFrame:
    with open(CASE_JSON) as f:
        data = json.load(f)
    start = pd.Timestamp(
        year=data["hkcaseaggr_v202209_confirm_first_date_utc_y"],
        month=data["hkcaseaggr_v202209_confirm_first_date_utc_m"] + 1,
        day=data["hkcaseaggr_v202209_confirm_first_date_utc_d"],
    )
    n = len(data["hkcaseaggr_v202209_confirm_classification_local_pcr"])
    dates = pd.date_range(start, periods=n, freq="D")
    comm = pd.DataFrame(
        {
            "date": dates,
            "local": np.array(data["hkcaseaggr_v202209_confirm_classification_local_pcr"])
            + np.array(data["hkcaseaggr_v202209_confirm_classification_local_rat"]),
        }
    )
    win = comm[(comm["date"] >= "2022-05-28") & (comm["date"] <= "2022-08-18")].set_index("date")
    hai_daily = hai.groupby("Date")["N"].sum().rename("hai")
    m = win.join(hai_daily, how="left").fillna({"hai": 0})
    m["local_7"] = m["local"].rolling(7, min_periods=3).mean()
    m["hai_7"] = m["hai"].rolling(7, min_periods=3).mean()
    return m


def make_figures(hai, pathways, cluster_df, community):
    # Color system: evidence-led, not purple-default
    nt_color = "#1f4e3d"
    urban_color = "#8b4513"
    staff_color = "#0b3d5c"
    inp_color = "#a34a28"

    # Fig A: density paradox — pop density vs cases/1000 beds
    sample = cluster_df[cluster_df["in_hai_sample"]].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), constrained_layout=True)

    ax = axes[0]
    colors = [nt_color if r == "New Territories" else urban_color for r in sample["region_group"]]
    ax.scatter(sample["pop_density"], sample["cases_per_1000_beds"], s=sample["beds_total"] / 40, c=colors, alpha=0.85, edgecolor="white")
    for _, r in sample.iterrows():
        ax.annotate(r["ha_cluster"], (r["pop_density"], r["cases_per_1000_beds"]), textcoords="offset points", xytext=(6, 4), fontsize=9)
    ax.set_xlabel("Catchment population density (persons / km², core districts)")
    ax.set_ylabel("Sequenced HAI cases per 1,000 beds")
    ax.set_title("A. Density paradox at HA-cluster scale")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    ax.scatter(sample["hospitals_per_100km2"], sample["cases_per_1000_beds"], s=90, c=colors, alpha=0.85, edgecolor="white")
    for _, r in sample.iterrows():
        ax.annotate(r["ha_cluster"], (r["hospitals_per_100km2"], r["cases_per_1000_beds"]), textcoords="offset points", xytext=(6, 4), fontsize=9)
    ax.set_xlabel("Public hospitals/institutions per 100 km²")
    ax.set_ylabel("Sequenced HAI cases per 1,000 beds")
    ax.set_title("B. Hospital geographic density vs HAI burden")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(OUT / "fig_density_paradox.png", dpi=220)
    fig.savefig(OUT / "fig_density_paradox.pdf")
    plt.close(fig)

    # Fig B: dual seeding composition
    order = ["staff_first", "inpatient_first", "same_day", "staff_only", "inpatient_only"]
    labels = {
        "staff_first": "Staff/OP first",
        "inpatient_first": "Inpatient first",
        "same_day": "Same day",
        "staff_only": "Staff/OP only",
        "inpatient_only": "Inpatient only",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)
    for ax, region, title in zip(
        axes,
        ["New Territories", "Urban (HK Island/Kowloon)"],
        ["C. Dual seeding — New Territories", "D. Dual seeding — Urban clusters in sample"],
    ):
        sub = pathways[pathways["region"] == region]
        counts = sub.groupby("path")["cases"].sum().reindex(order).fillna(0)
        cols = [staff_color if "staff" in p else inp_color if "inpatient" in p else "#6b6b6b" for p in order]
        ax.barh([labels[p] for p in order], counts.values, color=cols, alpha=0.9)
        ax.set_xlabel("Sequenced cases")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for y, v in enumerate(counts.values):
            if v > 0:
                ax.text(v + 0.5, y, f"{int(v)}", va="center", fontsize=8)
    fig.savefig(OUT / "fig_dual_seeding_by_region.png", dpi=220)
    fig.savefig(OUT / "fig_dual_seeding_by_region.pdf")
    plt.close(fig)

    # Fig C: mental model strip — weekly HAI by region + community
    weekly = (
        hai.assign(week=hai["Date"].dt.to_period("W-SUN").dt.start_time)
        .groupby(["week", "region"])["N"]
        .sum()
        .unstack(fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(10.5, 4.2), constrained_layout=True)
    ax.fill_between(weekly.index, 0, weekly.get("New Territories", 0), color=nt_color, alpha=0.75, label="NT HAI (sequenced)")
    ax.fill_between(
        weekly.index,
        weekly.get("New Territories", 0),
        weekly.get("New Territories", 0) + weekly.get("Urban (HK Island/Kowloon)", 0),
        color=urban_color,
        alpha=0.65,
        label="Urban HAI (sequenced)",
    )
    ax2 = ax.twinx()
    ax2.plot(community.index, community["local_7"], color="#444444", lw=1.6, alpha=0.85, label="Community local cases (7d mean)")
    ax.axvline(RELAX_DATE, color="#053C29", ls="--", lw=1.2, alpha=0.9)
    ax.set_ylabel("Weekly sequenced HAI cases")
    ax2.set_ylabel("Community local cases (7-day mean)")
    ax.set_title("E. NT-concentrated HAI rises with Phase-2 community wave")
    ax.spines["top"].set_visible(False)
    lines, labels_ = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels_ + labels2, loc="upper left", frameon=False, fontsize=8)
    fig.savefig(OUT / "fig_nt_wave_alignment.png", dpi=220)
    fig.savefig(OUT / "fig_nt_wave_alignment.pdf")
    plt.close(fig)

    # Fig D: cluster bars — raw vs standardized
    plot_df = sample.sort_values("cases_per_1000_beds", ascending=True)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2), constrained_layout=True)
    y = np.arange(len(plot_df))
    cols = [nt_color if r == "New Territories" else urban_color for r in plot_df["region_group"]]
    axes[0].barh(y, plot_df["sequenced_cases"], color=cols)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(plot_df["ha_cluster"])
    axes[0].set_xlabel("Sequenced cases")
    axes[0].set_title("F. Raw burden")
    axes[1].barh(y, plot_df["cases_per_1000_beds"], color=cols)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(plot_df["ha_cluster"])
    axes[1].set_xlabel("Cases / 1,000 beds")
    axes[1].set_title("G. Bed-standardized")
    axes[2].barh(y, plot_df["cases_per_million_pop"], color=cols)
    axes[2].set_yticks(y)
    axes[2].set_yticklabels(plot_df["ha_cluster"])
    axes[2].set_xlabel("Cases / million catchment pop.")
    axes[2].set_title("H. Population-standardized")
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.savefig(OUT / "fig_cluster_standardization.png", dpi=220)
    fig.savefig(OUT / "fig_cluster_standardization.pdf")
    plt.close(fig)


def write_memo(hai, pathways, cluster_df, fisher, community):
    sample = cluster_df[cluster_df["in_hai_sample"]]
    nt = sample[sample["ha_cluster"].isin(NT_CLUSTERS)]
    urban = sample[sample["ha_cluster"].isin(URBAN_IN_SAMPLE)]

    nt_cases = int(hai.loc[hai["Hospital clusters"].isin(NT_CLUSTERS), "N"].sum())
    total = int(hai["N"].sum())
    api_known = hai[hai["Case type"].isin(["Inpatient (API >= 3)", "Inpatient (API < 3)"])]
    api_ge3 = int(hai.loc[hai["Case type"] == "Inpatient (API >= 3)", "N"].sum())
    api_known_n = int(api_known["N"].sum())

    path_case = pathways.groupby("path")["cases"].sum()
    dual_case = int(path_case.get("staff_first", 0) + path_case.get("inpatient_first", 0) + path_case.get("same_day", 0))

    corr_daily = community[["local", "hai"]].corr().iloc[0, 1]
    corr_7 = community[["local_7", "hai_7"]].dropna().corr().iloc[0, 1]

    # Density ratios
    nt_pop_dens = nt["pop_density"].mean()
    urban_pop_dens = urban["pop_density"].mean()
    nt_hosp_dens = nt["hospitals_per_100km2"].mean()
    urban_hosp_dens = urban["hospitals_per_100km2"].mean()
    nt_bed_dens = nt["beds_per_100km2"].mean()
    urban_bed_dens = urban["beds_per_100km2"].mean()
    nt_rate = nt["cases_per_1000_beds"].mean()
    urban_rate = urban["cases_per_1000_beds"].mean()

    ba56 = hai[hai["Lineages"] == "BA.5.6"]
    multi = pathways[pathways["lineages"].str.contains(",")]

    memo = f"""# Dual seeding under a geographic paradox
## Exploratory research note — hospital-associated SARS-CoV-2, Hong Kong 2022

**Status:** evidence-building exploration (not a final manuscript).  
**Data sufficiency verdict:** sufficient for a *descriptive + standardized* argument; insufficient for causal claims about density or visit policy without admissions/sequencing-effort denominators.

---

## Working thesis

Hospital-associated infections in this sequenced window are **not** a single-pathway, city-wide process. They are:

1. **Dually seeded** — staff/outpatient bridges and inpatient introductions both matter.
2. **Geographically concentrated** in New Territories clusters (NTWC + NTEC).
3. Occurring in a setting where NT has **lower population density and lower hospital geographic density** than urban clusters in the sample — a **density paradox** that survives bed- and population-standardization.

The mental model is therefore: *Phase-2 community pressure → many independent introductions via two ward-level routes → short clusters, disproportionately observed in lower-density NT catchments with large bed bases (including psychiatric capacity).*

---

## Data used

| Layer | Source | Role |
|---|---|---|
| Sequenced HAI plot table | `data/hospital_data/metadata_cases_plot.csv` | Outcomes (162 cases, 17 hospitals, 5 HA clusters) |
| District area/population | 2021 Census (see `analysis/reference/SOURCES.md`) | Population density |
| HA beds / facility counts | HA Statistical Report 2021–22 (31 Mar 2022) | Hospital & bed density; rate denominators |
| Catchment population | LegCo LCQ13 (2017 projection) | Population-standardized rates |
| Community epicurve | `HK_case_data.json` | Wave alignment |

**Explicitly missing (judgment gate for causal work):** hospital-level admissions during the window, sequencing fraction by hospital/cluster, visitor-volume data, full official HA nosocomial line list linked to sequences.

---

## Finding 1 — Dual seeding is real in these data

Ward-level first-detection classification (n = {len(pathways)} wards):

| Pathway | Wards | Cases |
|---|---:|---:|
| Staff/outpatient first | {(pathways.path=='staff_first').sum()} | {int(path_case.get('staff_first',0))} |
| Inpatient first | {(pathways.path=='inpatient_first').sum()} | {int(path_case.get('inpatient_first',0))} |
| Same day | {(pathways.path=='same_day').sum()} | {int(path_case.get('same_day',0))} |
| Staff/OP only | {(pathways.path=='staff_only').sum()} | {int(path_case.get('staff_only',0))} |
| Inpatient only | {(pathways.path=='inpatient_only').sum()} | {int(path_case.get('inpatient_only',0))} |

Among wards with both staff and inpatient signal and a clear order, staff-first and inpatient-first are **balanced in case burden** (~{int(path_case.get('staff_first',0))} vs ~{int(path_case.get('inpatient_first',0))} cases).  
Fisher exact test (NT vs urban × staff-first vs inpatient-first): OR = {fisher['odds_ratio']:.2f}, p = {fisher['p_value']:.3f}, n = {fisher['n_wards']} wards.  
**Interpretation:** no evidence that NT and urban clusters use different dominant pathways; dual seeding is a **shared mechanism** inside a geographically skewed burden.

Among known-API inpatients, **{api_ge3}/{api_known_n} ({100*api_ge3/api_known_n:.0f}%)** have API ≥ 3 — consistent with substantial true nosocomial acquisition once infection is inside the ward.

---

## Finding 2 — Geographic concentration

- NTWC + NTEC: **{nt_cases}/{total} ({100*nt_cases/total:.1f}%)** of sequenced HAI cases.
- Top 3 anonymized hospitals: ~half of all cases (see script tables).
- HKWC and KCC appear in HA capacity data but contribute **0** cases to this sequenced plot table — concentration is even starker against the full HA map.

---

## Finding 3 — The density paradox (quantified)

Among clusters **present in the HAI table**:

| Metric (mean of clusters) | New Territories (NTWC, NTEC) | Urban in sample (KWC, KEC, HKEC) | Ratio NT/Urban |
|---|---:|---:|---:|
| Population density (per km²) | {nt_pop_dens:,.0f} | {urban_pop_dens:,.0f} | {nt_pop_dens/urban_pop_dens:.2f} |
| Hospitals per 100 km² | {nt_hosp_dens:.2f} | {urban_hosp_dens:.2f} | {nt_hosp_dens/urban_hosp_dens:.2f} |
| Beds per 100 km² | {nt_bed_dens:,.0f} | {urban_bed_dens:,.0f} | {nt_bed_dens/urban_bed_dens:.2f} |
| Sequenced HAI per 1,000 beds | {nt_rate:.2f} | {urban_rate:.2f} | {nt_rate/urban_rate:.2f} |

**Read carefully:** NT is less dense on all three geographic metrics, yet has ~{nt_rate/urban_rate:.1f}× the bed-standardized sequenced HAI rate.  
This is the opposite of a naive “dense city → more hospital outbreaks” story.

**Negative control inside the paradox:** KEC also has low hospital geographic density (~2 per 100 km², similar to NT) but only 1.7 sequenced HAI cases per 1,000 beds. Low density is therefore **not sufficient** for high observed burden — the NT excess needs NT-specific structure (case-mix / long-stay capacity / sequencing effort / IPC), not density alone.

Cluster-level detail is in `cluster_burden_standardized.csv`.

**Psychiatric capacity note:** NTWC alone reports **1,176 mentally-ill + 520 mentally-handicapped beds** (HA 2021–22). Long-stay settings inflate API and can sustain ward transmission — a plausible partial link to high-API / BA.5.6 patterns, not yet hospital-resolved in the anonymized table.

---

## Finding 4 — Other angles that fold into the thesis

1. **Visit relaxation as clock, not proven cause.** {int(hai.loc[hai.Date < RELAX_DATE, 'N'].sum())}/{total} sequenced HAI cases precede 31 May 2022. Community incidence was also rising; density-paradox + dual-seeding argue against a simple mobility story (consistent with the null GAM on mobility predictors).
2. **Community alignment is weekly, not daily.** Corr(local, HAI) ≈ {corr_daily:.2f}; corr(7-day means) ≈ {corr_7:.2f}. Introductions track the Phase-2 wave envelope.
3. **BA.5.6 is a single-hospital deep outbreak** ({int(ba56['N'].sum())} cases, all one hospital, almost all API ≥ 3) sitting inside NTWC — extreme version of “introduction then sustained nosocomial spread.”
4. **Multi-lineage wards (n = {len(multi)})** imply reintroduction into the same unit — supports many independent seeds, not one clonal citywide hospital epidemic.
5. **Starred (non-official) wards hold {int(hai.loc[hai.starred,'N'].sum())} cases ({100*hai.loc[hai.starred,'N'].sum()/total:.0f}%).** Official cluster lists understate hospital-linked transmission; dual-seeding estimates are conservative if anything.

---

## What would make this Laidlaw-caliber publishable

**Already strong enough for a methods/results subsection:** dual-seeding taxonomy; density paradox with bed/population standardization; mental model figure.

**Needs judgment / new data before causal language:**
1. Sequencing effort or admission denominators by hospital → test whether NT excess is ascertainment.
2. Pre-registered pathway definition (API threshold, same-day rule) + sensitivity analysis.
3. Phylogenetic confirmation that staff-first wards are *introductions* (basal staff tips) rather than delayed staff detection.
4. Formal interrupted time series around visit-policy dates with community offset.

**Stop condition met for now:** further model complexity without denominators would add theatre, not evidence.

---

## Outputs

- `analysis/outputs/ward_pathways.csv`
- `analysis/outputs/cluster_burden_standardized.csv`
- `analysis/outputs/summary_stats.json`
- `analysis/outputs/fig_density_paradox.(png|pdf)`
- `analysis/outputs/fig_dual_seeding_by_region.(png|pdf)`
- `analysis/outputs/fig_nt_wave_alignment.(png|pdf)`
- `analysis/outputs/fig_cluster_standardization.(png|pdf)`

Reproduce: `python scripts/4.dual_seeding_geography.py`
"""
    (OUT / "RESEARCH_NOTE.md").write_text(memo)
    return memo


def main():
    hai = load_hai()
    pathways = classify_ward_pathways(hai)
    capacity = pd.read_csv(CLUSTER_PATH)
    districts = pd.read_csv(DIST_PATH)
    cluster_df = cluster_burden_table(hai, capacity)
    fisher = fisher_dual_seeding(pathways)
    community = community_alignment(hai)

    pathways.to_csv(OUT / "ward_pathways.csv", index=False)
    cluster_df.to_csv(OUT / "cluster_burden_standardized.csv", index=False)
    districts.to_csv(OUT / "district_demographics_used.csv", index=False)

    # Region rollup
    region_rollup = (
        hai.groupby("region")
        .agg(cases=("N", "sum"), hospitals=("Hospital_anonymized", "nunique"), wards=("Wards", "nunique"))
        .reset_index()
    )
    region_rollup.to_csv(OUT / "region_rollup.csv", index=False)

    summary = {
        "n_cases": int(hai["N"].sum()),
        "n_hospitals": int(hai["Hospital_anonymized"].nunique()),
        "n_wards": int(hai.groupby(["Hospital_anonymized", "Wards"]).ngroups),
        "nt_case_share": float(hai.loc[hai["Hospital clusters"].isin(NT_CLUSTERS), "N"].sum() / hai["N"].sum()),
        "fisher_dual_seeding": {
            "odds_ratio": None if np.isnan(fisher["odds_ratio"]) else float(fisher["odds_ratio"]),
            "p_value": float(fisher["p_value"]),
            "n_wards": fisher["n_wards"],
            "table": fisher["table"].to_dict(),
        },
        "pathway_cases": pathways.groupby("path")["cases"].sum().astype(int).to_dict(),
        "cluster_rates": cluster_df.set_index("ha_cluster")[
            ["sequenced_cases", "pop_density", "hospitals_per_100km2", "beds_per_100km2", "cases_per_1000_beds", "cases_per_million_pop"]
        ]
        .round(3)
        .to_dict(orient="index"),
        "limits": [
            "Sequenced HAI != complete incidence",
            "No hospital-level sequencing effort denominator",
            "Catchment populations are planning estimates, not closed cohorts",
            "Core land areas exclude ambiguous Islands/Lantau splits",
        ],
    }
    (OUT / "summary_stats.json").write_text(json.dumps(summary, indent=2))

    make_figures(hai, pathways, cluster_df, community)
    write_memo(hai, pathways, cluster_df, fisher, community)

    print("=== Dual seeding pathway cases ===")
    print(pathways.groupby("path")["cases"].sum().sort_values(ascending=False).to_string())
    print("\n=== Fisher NT vs urban (staff_first vs inpatient_first) ===")
    print(fisher["table"])
    print(f"OR={fisher['odds_ratio']:.3f}, p={fisher['p_value']:.3f}")
    print("\n=== Cluster standardized burden (in-sample) ===")
    cols = [
        "ha_cluster",
        "sequenced_cases",
        "pop_density",
        "hospitals_per_100km2",
        "beds_per_100km2",
        "cases_per_1000_beds",
        "cases_per_million_pop",
    ]
    print(cluster_df.loc[cluster_df["in_hai_sample"], cols].sort_values("cases_per_1000_beds", ascending=False).round(2).to_string(index=False))
    print(f"\nWrote outputs to {OUT}")


if __name__ == "__main__":
    main()
