#!/usr/bin/env python3
"""
Geographic representation note (peer-review reconceptualisation)
================================================================
Primary estimand: geographic representation of sequenced confirmed-ward
inpatients relative to HA cluster throughput (FY2022-23 patient-days primary).

Adds:
- all-system (zeros) vs represented-cluster-only sample-concentration indices
- selection-bias threshold analysis
- ward/hospital influence metrics and leave-one-out plot
- binary official-presence audit
- secondary first-sampling chronology (supplement outputs only)

Does NOT claim incidence, density paradox, dual seeding, or psychiatric causation.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "analysis/reference"
OUT = ROOT / "analysis/outputs"
FIG = ROOT / "manuscript/figures"
SUPP = ROOT / "manuscript/supplement"
HAI = ROOT / "data/hospital_data/metadata_cases_plot.csv"
CAPACITY = REF / "ha_cluster_capacity_2022.csv"
THR_OPEN_22 = REF / "ha_cluster_throughput_open_2022_23.csv"
THR_OPEN_21 = REF / "ha_cluster_throughput_open_2021_22.csv"
ASCERTAIN = REF / "ascertainment_presence_by_cluster.csv"

NT = {"NTWC", "NTEC"}
for p in (OUT, FIG, SUPP):
    p.mkdir(parents=True, exist_ok=True)


def load_confirmed() -> pd.DataFrame:
    hai = pd.read_csv(HAI, parse_dates=["Date"])
    hai["confirmed"] = ~hai["Wards"].astype(str).str.endswith("*")
    hai["non_inpatient"] = hai["Case type"].eq("Outpatient/Staff")
    return hai[hai["confirmed"]].copy()


def cluster_counts(hai: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cl, g in hai.groupby("Hospital clusters"):
        rows.append(
            {
                "ha_cluster": cl,
                "sequenced_all": int(g["N"].sum()),
                "sequenced_inpatient": int(g.loc[~g["non_inpatient"], "N"].sum()),
                "sequenced_non_inpatient": int(g.loc[g["non_inpatient"], "N"].sum()),
                "n_hospitals": int(g["Hospital_anonymized"].nunique()),
                "n_wards": int(g["Wards"].nunique()),
            }
        )
    # ensure all 7 clusters present
    all_c = ["HKEC", "HKWC", "KCC", "KEC", "KWC", "NTEC", "NTWC"]
    out = pd.DataFrame(rows).set_index("ha_cluster").reindex(all_c).fillna(0).reset_index()
    for c in ["sequenced_all", "sequenced_inpatient", "sequenced_non_inpatient", "n_hospitals", "n_wards"]:
        out[c] = out[c].astype(int)
    return out


def build_denom_table(counts: pd.DataFrame) -> pd.DataFrame:
    cap = pd.read_csv(CAPACITY)[["ha_cluster", "beds_total", "region_group"]]
    thr = pd.read_csv(THR_OPEN_22)[["ha_cluster", "ip_discharges", "patient_days"]]
    thr = thr.rename(columns={"ip_discharges": "ip_discharges_fy2022_23", "patient_days": "patient_days_fy2022_23"})
    thr21 = pd.read_csv(THR_OPEN_21)[["ha_cluster", "ip_discharges", "patient_days"]].rename(
        columns={"ip_discharges": "ip_discharges_fy2021_22", "patient_days": "patient_days_fy2021_22"}
    )
    asc = pd.read_csv(ASCERTAIN) if ASCERTAIN.exists() else pd.DataFrame()
    d = counts.merge(cap, on="ha_cluster", how="left").merge(thr, on="ha_cluster", how="left").merge(
        thr21, on="ha_cluster", how="left"
    )
    if len(asc):
        d = d.merge(asc[["ha_cluster", "in_official_bulletins", "in_sequenced_table", "mention_days"]], on="ha_cluster", how="left")
    d["in_genomic_sample"] = d["sequenced_all"] > 0
    d["region"] = np.where(d["ha_cluster"].isin(NT), "New Territories", "non-New Territories")
    # sample-concentration indices (NOT incidence rates)
    d["sci_per_100k_patient_days"] = d["sequenced_inpatient"] / d["patient_days_fy2022_23"] * 1e5
    d["sci_per_1000_beds"] = d["sequenced_inpatient"] / d["beds_total"] * 1000
    d["sci_per_100k_discharges"] = d["sequenced_inpatient"] / d["ip_discharges_fy2022_23"] * 1e5
    # shares
    d["share_sequenced_inpatients"] = d["sequenced_inpatient"] / d["sequenced_inpatient"].sum()
    d["share_patient_days"] = d["patient_days_fy2022_23"] / d["patient_days_fy2022_23"].sum()
    d["representation_ratio"] = d["share_sequenced_inpatients"] / d["share_patient_days"]
    return d


def pooled_ratio(d: pd.DataFrame, num: str, den: str, scale: float, scope: str) -> dict:
    if scope == "represented":
        d = d[d["in_genomic_sample"]].copy()
    nt = d[d["ha_cluster"].isin(NT)]
    oth = d[~d["ha_cluster"].isin(NT)]
    r_nt = nt[num].sum() / nt[den].sum() * scale
    r_oth = oth[num].sum() / oth[den].sum() * scale
    return {
        "scope": scope,
        "denominator": den,
        "sci_nt": float(r_nt),
        "sci_other": float(r_oth),
        "ratio_nt_other": float(r_nt / r_oth) if r_oth else None,
    }


def selection_threshold(R_obs: float) -> dict:
    """
    Y = s * I  =>  R_obs = R_true * (s_NT / s_other).
    If R_true = 1, required inclusion-probability ratio = R_obs.
    Also tabulate R_true implied for a grid of s ratios.
    """
    grid = []
    for s_ratio in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.2, 5.0, 6.0, 7.0, 8.0]:
        grid.append(
            {
                "inclusion_prob_ratio_nt_over_other": s_ratio,
                "implied_true_rate_ratio_if_Robs": float(R_obs / s_ratio),
            }
        )
    return {
        "R_obs_patient_days_represented": R_obs,
        "inclusion_ratio_needed_for_equal_true_rates": R_obs,
        "grid": grid,
        "equation": "R_obs = R_true * (s_NT / s_other)",
    }


def ward_influence(hai: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    rows = []
    for (h, w, cl), g in hai.groupby(["Hospital_anonymized", "Wards", "Hospital clusters"]):
        rows.append(
            {
                "Hospital_anonymized": h,
                "Wards": w,
                "ha_cluster": cl,
                "region": "New Territories" if cl in NT else "non-New Territories",
                "sequenced_all": int(g["N"].sum()),
                "sequenced_inpatient": int(g.loc[~g["non_inpatient"], "N"].sum()),
            }
        )
    W = pd.DataFrame(rows).sort_values("sequenced_inpatient", ascending=False)
    total_ip = W["sequenced_inpatient"].sum()
    summary = {
        "n_wards": int(len(W)),
        "n_wards_nt": int((W["region"] == "New Territories").sum()),
        "n_wards_other": int((W["region"] == "non-New Territories").sum()),
        "n_hospitals": int(W["Hospital_anonymized"].nunique()),
        "median_inpatients_per_ward": float(W["sequenced_inpatient"].median()),
        "min_inpatients_per_ward": int(W["sequenced_inpatient"].min()),
        "max_inpatients_per_ward": int(W["sequenced_inpatient"].max()),
        "largest_ward_share_of_inpatients": float(W["sequenced_inpatient"].iloc[0] / total_ip),
        "top3_wards_share_of_inpatients": float(W["sequenced_inpatient"].head(3).sum() / total_ip),
        "largest_ward_label": f"{W['Hospital_anonymized'].iloc[0]}:{W['Wards'].iloc[0]}",
    }
    return W, summary


def leave_one_out(hai: pd.DataFrame, den: pd.DataFrame) -> pd.DataFrame:
    """Exclude one ward at a time; recompute represented-cluster patient-day SCI ratio."""
    wards = list(hai[["Hospital_anonymized", "Wards"]].drop_duplicates().itertuples(index=False))
    rows = []
    base_pd = den.set_index("ha_cluster")["patient_days_fy2022_23"]
    base_beds = den.set_index("ha_cluster")["beds_total"]
    for h, w in wards:
        sub = hai[~((hai["Hospital_anonymized"] == h) & (hai["Wards"] == w))]
        cnt = cluster_counts(sub)
        # merge dens
        m = cnt.set_index("ha_cluster")
        nt_ip = m.loc[list(NT), "sequenced_inpatient"].sum()
        oth_clusters = [c for c in m.index if c not in NT and m.loc[c, "sequenced_all"] > 0]
        # represented other among those still with sample
        oth_ip = m.loc[[c for c in m.index if c not in NT], "sequenced_inpatient"].sum()
        # use all non-NT with any remaining sample for represented scope
        rep_oth = [c for c in m.index if c not in NT and m.loc[c, "sequenced_all"] > 0]
        nt_pd = base_pd.loc[list(NT)].sum()
        oth_pd = base_pd.loc[rep_oth].sum() if rep_oth else np.nan
        nt_beds = base_beds.loc[list(NT)].sum()
        oth_beds = base_beds.loc[rep_oth].sum() if rep_oth else np.nan
        ratio_pd = (nt_ip / nt_pd) / (oth_ip / oth_pd) if oth_pd and oth_ip else np.nan
        ratio_beds = (nt_ip / nt_beds) / (oth_ip / oth_beds) if oth_beds and oth_ip else np.nan
        rows.append(
            {
                "excluded": f"{h}:{w}",
                "remaining_inpatients": int(sub.loc[~sub["non_inpatient"], "N"].sum()),
                "ratio_patient_days_represented": float(ratio_pd) if ratio_pd == ratio_pd else None,
                "ratio_beds_represented": float(ratio_beds) if ratio_beds == ratio_beds else None,
            }
        )
    return pd.DataFrame(rows).sort_values("ratio_patient_days_represented")


def first_sampling_secondary(hai_all_path: Path) -> pd.DataFrame:
    """Secondary descriptive chronology; not a transmission analysis."""
    hai = pd.read_csv(hai_all_path, parse_dates=["Date"])
    hai["confirmed"] = ~hai["Wards"].astype(str).str.endswith("*")
    hai["non_inpatient"] = hai["Case type"].eq("Outpatient/Staff")
    rows = []
    for (h, w), g in hai.groupby(["Hospital_anonymized", "Wards"]):
        ni = g.loc[g["non_inpatient"], "Date"]
        ip = g.loc[~g["non_inpatient"], "Date"]
        t_ni = ni.min() if len(ni) else pd.NaT
        t_ip = ip.min() if len(ip) else pd.NaT
        if pd.isna(t_ni) and pd.isna(t_ip):
            path = "unknown"
        elif pd.isna(t_ni):
            path = "inpatient_only"
        elif pd.isna(t_ip):
            path = "non_inpatient_only"
        elif t_ni < t_ip:
            path = "non_inpatient_before_inpatient"
        elif t_ip < t_ni:
            path = "inpatient_before_non_inpatient"
        else:
            path = "same_day_first_sample"
        rows.append(
            {
                "Hospital_anonymized": h,
                "Wards": w,
                "confirmed": bool(g["confirmed"].iloc[0]),
                "ha_cluster": g["Hospital clusters"].iloc[0],
                "first_sampling_order": path,
                "cases": int(g["N"].sum()),
            }
        )
    return pd.DataFrame(rows)


def make_figures(den: pd.DataFrame, loo: pd.DataFrame, sel: dict, ward_sum: dict):
    plt.rcParams.update({"font.size": 10, "figure.dpi": 150})
    # colour-blind-ish
    c_nt, c_oth = "#0072B2", "#E69F00"

    # Fig 1: sample share vs patient-day share
    d = den.sort_values("share_patient_days", ascending=False)
    x = np.arange(len(d))
    w = 0.38
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    ax.bar(x - w / 2, d["share_patient_days"], w, label="Share of FY2022–23 patient-days", color="#999999")
    ax.bar(x + w / 2, d["share_sequenced_inpatients"], w, label="Share of sequenced inpatients", color="#0072B2")
    ax.set_xticks(x)
    ax.set_xticklabels(d["ha_cluster"])
    ax.set_ylabel("Share")
    ax.set_title("Geographic representation: sequenced inpatients vs inpatient throughput")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_share_vs_throughput.pdf")
    fig.savefig(OUT / "fig_share_vs_throughput.pdf")
    plt.close(fig)

    # Fig 2: SCI by cluster (patient-days); zeros shown
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    colors = [c_nt if c in NT else c_oth for c in d["ha_cluster"]]
    ax.bar(d["ha_cluster"], d["sci_per_100k_patient_days"], color=colors)
    ax.set_ylabel("Sample-concentration index\n(sequenced inpatients / 100k patient-days)")
    ax.set_title("Observed genomic representation by HA cluster (0 = absent from sample)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_sci_patient_days.pdf")
    fig.savefig(OUT / "fig_sci_patient_days.pdf")
    plt.close(fig)

    # Fig 3: leave-one-ward-out
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    y = loo["ratio_patient_days_represented"].values
    ax.plot(range(len(y)), y, "o-", color="#0072B2", ms=4)
    ax.axhline(sel["R_obs_patient_days_represented"], color="#D55E00", ls="--", label="Full-sample ratio")
    ax.set_xlabel("Leave-one-ward-out analyses (sorted)")
    ax.set_ylabel("NT / non-NT SCI ratio (patient-days)")
    ax.set_title("Influence of individual wards on the represented-cluster concentration ratio")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_loo_ward_ratio.pdf")
    fig.savefig(OUT / "fig_loo_ward_ratio.pdf")
    plt.close(fig)

    # Fig 4: selection threshold
    grid = pd.DataFrame(sel["grid"])
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.plot(grid["inclusion_prob_ratio_nt_over_other"], grid["implied_true_rate_ratio_if_Robs"], "o-", color="#009E73")
    ax.axhline(1.0, color="gray", ls="--")
    ax.axvline(sel["inclusion_ratio_needed_for_equal_true_rates"], color="#D55E00", ls=":",
               label=f"s ratio for R_true=1 (≈{sel['inclusion_ratio_needed_for_equal_true_rates']:.1f})")
    ax.set_xlabel(r"Inclusion-probability ratio $s_{\mathrm{NT}}/s_{\mathrm{other}}$")
    ax.set_ylabel(r"Implied true rate ratio $R_{\mathrm{true}}$")
    ax.set_title(r"Selection-bias threshold: $R_{\mathrm{obs}}=R_{\mathrm{true}}\cdot(s_{\mathrm{NT}}/s_{\mathrm{other}})$")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_selection_threshold.pdf")
    fig.savefig(OUT / "fig_selection_threshold.pdf")
    plt.close(fig)

    # Fig 5: binary presence
    if "in_official_bulletins" in den.columns:
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        # categorical plot
        for i, r in den.iterrows():
            x = 1 if r.get("in_official_bulletins") in (True, "True", "yes", 1) or str(r.get("in_official_bulletins")).lower() == "true" else 0
            # ascertainment file may use bool-like
            off = r.get("in_official_bulletins")
            if isinstance(off, str):
                x = 1 if off.lower() in ("true", "yes", "1") else 0
            else:
                x = 1 if bool(off) else 0
            y = 1 if r["in_genomic_sample"] else 0
            ax.scatter(x + np.random.uniform(-0.03, 0.03), y + np.random.uniform(-0.03, 0.03),
                       s=80, c=c_nt if r["ha_cluster"] in NT else c_oth)
            ax.annotate(r["ha_cluster"], (x, y), textcoords="offset points", xytext=(6, 4), fontsize=8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["No official bulletin\npresence in window", "Official bulletin\npresence in window"])
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["No confirmed-ward\ngenomes in dataset", "Confirmed-ward\ngenomes present"])
        ax.set_xlim(-0.3, 1.3)
        ax.set_ylim(-0.3, 1.3)
        ax.set_title("Official reporting presence vs genomic-sample presence")
        fig.tight_layout()
        fig.savefig(FIG / "fig_binary_presence_audit.pdf")
        fig.savefig(OUT / "fig_binary_presence_audit.pdf")
        plt.close(fig)


def main():
    hai = load_confirmed()
    counts = cluster_counts(hai)
    den = build_denom_table(counts)

    ratios = []
    for scope in ("all_system", "represented"):
        for den_col, scale, label in [
            ("patient_days_fy2022_23", 1e5, "patient_days"),
            ("beds_total", 1000, "beds"),
            ("ip_discharges_fy2022_23", 1e5, "discharges"),
        ]:
            r = pooled_ratio(den, "sequenced_inpatient", den_col, scale, scope)
            r["label"] = label
            ratios.append(r)

    R_obs_pd = next(r["ratio_nt_other"] for r in ratios if r["scope"] == "represented" and r["label"] == "patient_days")
    R_obs_beds = next(r["ratio_nt_other"] for r in ratios if r["scope"] == "represented" and r["label"] == "beds")
    R_all_pd = next(r["ratio_nt_other"] for r in ratios if r["scope"] == "all_system" and r["label"] == "patient_days")
    R_all_beds = next(r["ratio_nt_other"] for r in ratios if r["scope"] == "all_system" and r["label"] == "beds")

    sel = selection_threshold(R_obs_pd)
    W, ward_sum = ward_influence(hai)
    loo = leave_one_out(hai, den)

    # exclude largest ward
    largest = ward_sum["largest_ward_label"]
    h, w = largest.split(":")
    hai_ex = hai[~((hai["Hospital_anonymized"] == h) & (hai["Wards"] == w))]
    den_ex = build_denom_table(cluster_counts(hai_ex))
    r_ex = pooled_ratio(den_ex, "sequenced_inpatient", "patient_days_fy2022_23", 1e5, "represented")
    r_ex_beds = pooled_ratio(den_ex, "sequenced_inpatient", "beds_total", 1000, "represented")

    # binary presence table
    presence = den[
        ["ha_cluster", "region", "in_official_bulletins", "in_genomic_sample", "sequenced_inpatient", "sequenced_all"]
    ].copy()
    presence["observed_genomes"] = presence["sequenced_all"]
    presence["underlying_burden"] = np.where(presence["in_genomic_sample"], "unknown (sample present)", "unknown (not identifiable from this sample)")

    # secondary first-sampling
    fs = first_sampling_secondary(HAI)
    fs_confirmed = fs[fs["confirmed"]]
    fs_sum = (
        fs_confirmed.groupby("first_sampling_order")
        .agg(wards=("Wards", "count"), cases=("cases", "sum"))
        .reset_index()
    )

    summary = {
        "cohort": {
            "confirmed_genomes": int(hai["N"].sum()),
            "confirmed_inpatients": int(hai.loc[~hai["non_inpatient"], "N"].sum()),
            "confirmed_non_inpatients": int(hai.loc[hai["non_inpatient"], "N"].sum()),
            "nt_share_inpatients": float(
                hai.loc[hai["Hospital clusters"].isin(NT) & ~hai["non_inpatient"], "N"].sum()
                / hai.loc[~hai["non_inpatient"], "N"].sum()
            ),
        },
        "ratios": ratios,
        "key_ratios": {
            "represented_patient_days": R_obs_pd,
            "represented_beds": R_obs_beds,
            "all_system_patient_days": R_all_pd,
            "all_system_beds": R_all_beds,
            "excl_largest_ward_patient_days": r_ex["ratio_nt_other"],
            "excl_largest_ward_beds": r_ex_beds["ratio_nt_other"],
        },
        "selection_threshold": sel,
        "ward_influence": ward_sum,
        "loo_patient_days_range": {
            "min": float(loo["ratio_patient_days_represented"].min()),
            "max": float(loo["ratio_patient_days_represented"].max()),
            "median": float(loo["ratio_patient_days_represented"].median()),
        },
        "table_s3_status": (
            "Supporting Table S3 cited in Gu et al. 2026 but not retrieved here "
            "(Wiley supplement HTTP 403; not in PMC full text). Routes: journal SI, "
            "corresponding authors, publisher support — not GitHub/Kaggle."
        ),
    }

    den.to_csv(OUT / "geo_denominator_table.csv", index=False)
    W.to_csv(OUT / "geo_ward_influence.csv", index=False)
    loo.to_csv(OUT / "geo_loo_ward.csv", index=False)
    presence.to_csv(OUT / "geo_binary_presence_audit.csv", index=False)
    pd.DataFrame(ratios).to_csv(OUT / "geo_pooled_ratios.csv", index=False)
    pd.DataFrame(sel["grid"]).to_csv(OUT / "geo_selection_threshold_grid.csv", index=False)
    fs.to_csv(SUPP / "supp_first_sampling_order.csv", index=False)
    fs_sum.to_csv(SUPP / "supp_first_sampling_summary_confirmed.csv", index=False)
    with open(OUT / "geo_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    make_figures(den, loo, sel, ward_sum)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
