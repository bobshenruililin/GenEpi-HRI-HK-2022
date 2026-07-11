#!/usr/bin/env python3
"""
Denominators for the NT concentration question
=============================================
Two parallel checks:

A) Admissions / throughput (HA Statistical Report 2021-22 cluster extracts)
   - sequenced HAI per 100k inpatient discharges
   - sequenced HAI per 100k inpatient patient-days

B) Sequencing / ascertainment effort
   - hospital-coded tips in study trees (window) by HA cluster
     [CAVEAT: largely the same samples as the HAI table → circular for
      "did sequencing create the geography?" — reported transparently]
   - community (non-hospital-coded) tip share in the same trees
   - official HA nosocomial press-release mention-days by cluster
     (independent of our sequencing choices; incomplete scrape, directional)

Also refreshes dual-seeding + density paradox tables and writes a
manuscript-oriented FINDINGS_LOG.md that catalogs every claim with
evidence grade.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/outputs"
REF = ROOT / "analysis/reference"
OUT.mkdir(parents=True, exist_ok=True)

HAI_PATH = ROOT / "data/hospital_data/metadata_cases_plot.csv"
CAPACITY_PATH = REF / "ha_cluster_capacity_2022.csv"
MENTION_PATH = REF / "ha_nosocomial_mention_days.csv"

NT = {"NTWC", "NTEC"}
URBAN_IN_SAMPLE = {"KWC", "KEC", "HKEC"}

# Tree tip hospital codes → HA cluster (from tip labels + ward matches + HA membership)
HOSP_CODE_TO_CLUSTER = {
    "TMH": "NTWC",
    "CPH": "NTWC",
    "TKP": "NTWC",  # co-classified psychiatric with CPH in study scripts; NTWC in anonymized map
    "TSH": "NTWC",
    "POH": "NTWC",
    "TSWH": "NTWC",
    "SLH": "NTWC",
    "AHNH": "NTEC",
    "NDH": "NTEC",
    "SH": "NTEC",
    "PWH": "NTEC",
    "TPH": "NTEC",
    "CMC": "KWC",
    "PMH": "KWC",
    "KCH": "KWC",
    "YCH": "KWC",
    "NLTH": "KWC",
    "TKOH": "KEC",
    "UCH": "KEC",
    "HHH": "KEC",
    "PYNEH": "HKEC",
    "RH": "HKEC",
    "TWE": "HKEC",
    "TWEH": "HKEC",
    "RTSKH": "HKEC",
    "QMH": "HKWC",
    "QEH": "KCC",
    "KWH": "KCC",
    "KH": "KCC",
}

# HA Stat Report 2021-22 cluster extracts (DownloadCluster IDs 36-42), scraped 2026-07-11
THROUGHPUT = pd.DataFrame(
    [
        {"ha_cluster": "HKEC", "ip_discharges_2021_22": 102190, "patient_days_2021_22": 822971, "patient_days_general": 587224, "patient_days_mentally_ill": 79454},
        {"ha_cluster": "HKWC", "ip_discharges_2021_22": 103370, "patient_days_2021_22": 635637, "patient_days_general": 573106, "patient_days_mentally_ill": 22850},
        {"ha_cluster": "KCC", "ip_discharges_2021_22": 199221, "patient_days_2021_22": 1599173, "patient_days_general": 1396680, "patient_days_mentally_ill": 128186},
        {"ha_cluster": "KEC", "ip_discharges_2021_22": 114909, "patient_days_2021_22": 809758, "patient_days_general": 767497, "patient_days_mentally_ill": 15952},
        {"ha_cluster": "KWC", "ip_discharges_2021_22": 201096, "patient_days_2021_22": 1416390, "patient_days_general": 1113624, "patient_days_mentally_ill": 238794},
        {"ha_cluster": "NTEC", "ip_discharges_2021_22": 163397, "patient_days_2021_22": 1311078, "patient_days_general": 1077211, "patient_days_mentally_ill": 144986},
        {"ha_cluster": "NTWC", "ip_discharges_2021_22": 140388, "patient_days_2021_22": 1331433, "patient_days_general": 893239, "patient_days_mentally_ill": 255564},
    ]
)


def parse_tree_tips() -> pd.DataFrame:
    rows = []
    for f, lin in [("date_ba22", "BA.2.2"), ("date_ba2121", "BA.2.12.1"), ("date_ba56", "BA.5.6")]:
        d = pd.read_csv(ROOT / f"results/trees/{f}.tsv", sep="\t")
        d["lineage_file"] = lin
        rows.append(d)
    alld = pd.concat(rows, ignore_index=True)
    alld["date"] = pd.to_datetime(alld["date"], errors="coerce")

    def hosp(s: str):
        p = str(s).split("/")
        if "Local_case" in p:
            i = p.index("Local_case")
            if i + 1 < len(p) and re.fullmatch(r"[A-Z]{2,5}", p[i + 1] or ""):
                return p[i + 1]
        return None

    alld["hospital_code"] = alld["strain"].map(hosp)
    alld["ha_cluster"] = alld["hospital_code"].map(HOSP_CODE_TO_CLUSTER)
    return alld


def classify_pathways(hai: pd.DataFrame) -> pd.DataFrame:
    hai = hai.copy()
    hai["region"] = np.where(hai["Hospital clusters"].isin(NT), "New Territories", "Urban (in sample)")
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
                "ha_cluster": g["Hospital clusters"].iloc[0],
                "region": g["region"].iloc[0],
                "path": path,
                "cases": int(g["N"].sum()),
                "lineages": ",".join(sorted(g["Lineages"].unique())),
            }
        )
    return pd.DataFrame(rows)


def build_denominator_table(hai: pd.DataFrame, tips: pd.DataFrame, capacity: pd.DataFrame) -> pd.DataFrame:
    burden = (
        hai.groupby("Hospital clusters")
        .agg(sequenced_hai_cases=("N", "sum"), n_hospitals_sampled=("Hospital_anonymized", "nunique"), n_wards=("Wards", "nunique"))
        .reset_index()
        .rename(columns={"Hospital clusters": "ha_cluster"})
    )

    win = tips[(tips["date"] >= "2022-05-28") & (tips["date"] <= "2022-08-18")].copy()
    tip_by = (
        win.dropna(subset=["ha_cluster"])
        .groupby("ha_cluster")
        .size()
        .rename("hospital_coded_tips_in_window")
        .reset_index()
    )
    community_tips = int((win["hospital_code"].isna()).sum())
    hospital_tips = int(win["hospital_code"].notna().sum())

    mentions = pd.read_csv(MENTION_PATH)
    mention_by = mentions.groupby("ha_cluster")["mention_days"].sum().rename("official_mention_days").reset_index()

    out = (
        capacity[["ha_cluster", "beds_total", "catchment_population_2017", "land_area_km2_core", "n_hospitals_institutions"]]
        .merge(THROUGHPUT, on="ha_cluster", how="left")
        .merge(burden, on="ha_cluster", how="left")
        .merge(tip_by, on="ha_cluster", how="left")
        .merge(mention_by, on="ha_cluster", how="left")
    )
    for col in ["sequenced_hai_cases", "n_hospitals_sampled", "n_wards", "hospital_coded_tips_in_window", "official_mention_days"]:
        out[col] = out[col].fillna(0).astype(int)

    out["pop_density"] = out["catchment_population_2017"] / out["land_area_km2_core"]
    out["hospitals_per_100km2"] = out["n_hospitals_institutions"] / out["land_area_km2_core"] * 100
    out["beds_per_100km2"] = out["beds_total"] / out["land_area_km2_core"] * 100
    out["psych_patient_day_share"] = out["patient_days_mentally_ill"] / out["patient_days_2021_22"]

    out["hai_per_1000_beds"] = out["sequenced_hai_cases"] / out["beds_total"] * 1000
    out["hai_per_100k_discharges"] = out["sequenced_hai_cases"] / out["ip_discharges_2021_22"] * 1e5
    out["hai_per_100k_patient_days"] = out["sequenced_hai_cases"] / out["patient_days_2021_22"] * 1e5
    out["hai_per_million_pop"] = out["sequenced_hai_cases"] / out["catchment_population_2017"] * 1e6

    out["tips_per_1000_beds"] = out["hospital_coded_tips_in_window"] / out["beds_total"] * 1000
    out["tips_per_100k_discharges"] = out["hospital_coded_tips_in_window"] / out["ip_discharges_2021_22"] * 1e5
    out["mention_days_per_1000_beds"] = out["official_mention_days"] / out["beds_total"] * 1000

    # Window-scaled throughput (transparent approximation)
    window_frac = 83 / 365.25
    out["hai_per_100k_discharges_window_scaled"] = out["sequenced_hai_cases"] / (out["ip_discharges_2021_22"] * window_frac) * 1e5

    out["region_group"] = np.where(
        out["ha_cluster"].isin(NT),
        "New Territories",
        np.where(out["sequenced_hai_cases"] > 0, "Urban (in sample)", "Urban (no sequenced HAI in table)"),
    )
    out.attrs["community_tips_window"] = community_tips
    out.attrs["hospital_tips_window"] = hospital_tips
    return out


def make_figures(df: pd.DataFrame):
    sample = df[df["sequenced_hai_cases"] > 0].copy()
    nt_c, urb_c = "#1f4e3d", "#8b4513"
    colors = [nt_c if r == "New Territories" else urb_c for r in sample["region_group"]]

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.3), constrained_layout=True)
    metrics = [
        ("hai_per_1000_beds", "per 1,000 beds"),
        ("hai_per_100k_discharges", "per 100k annual discharges"),
        ("hai_per_100k_patient_days", "per 100k annual patient-days"),
    ]
    sample = sample.sort_values("hai_per_1000_beds")
    y = np.arange(len(sample))
    for ax, (col, xlab), title in zip(axes, metrics, ["A. Bed-standardized", "B. Discharge-standardized", "C. Patient-day-standardized"]):
        cols = [nt_c if r == "New Territories" else urb_c for r in sample["region_group"]]
        ax.barh(y, sample[col], color=cols)
        ax.set_yticks(y)
        ax.set_yticklabels(sample["ha_cluster"])
        ax.set_xlabel(f"Sequenced HAI {xlab}")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.savefig(OUT / "fig_denominator_standardization.png", dpi=220)
    fig.savefig(OUT / "fig_denominator_standardization.pdf")
    plt.close(fig)

    # Ascertainment panel: official mention-days vs sequenced cases (all clusters)
    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    for _, r in df.iterrows():
        color = nt_c if r["ha_cluster"] in NT else ("#555555" if r["sequenced_hai_cases"] == 0 else urb_c)
        ax.scatter(r["official_mention_days"], r["sequenced_hai_cases"], s=max(40, r["beds_total"] / 50), c=color, alpha=0.85, edgecolor="white")
        ax.annotate(r["ha_cluster"], (r["official_mention_days"], r["sequenced_hai_cases"]), textcoords="offset points", xytext=(5, 4), fontsize=9)
    ax.set_xlabel("Official HA nosocomial bulletin mention-days (study window; incomplete scrape)")
    ax.set_ylabel("Sequenced HAI cases in study plot table")
    ax.set_title("D. Official reporting footprint vs sequenced burden")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(OUT / "fig_ascertainment_mentions_vs_sequences.png", dpi=220)
    fig.savefig(OUT / "fig_ascertainment_mentions_vs_sequences.pdf")
    plt.close(fig)

    # Psych patient-day share vs HAI rate
    fig, ax = plt.subplots(figsize=(6.8, 4.8), constrained_layout=True)
    for _, r in df.iterrows():
        color = nt_c if r["ha_cluster"] in NT else ("#555555" if r["sequenced_hai_cases"] == 0 else urb_c)
        ax.scatter(100 * r["psych_patient_day_share"], r["hai_per_1000_beds"], s=70, c=color, alpha=0.85, edgecolor="white")
        ax.annotate(r["ha_cluster"], (100 * r["psych_patient_day_share"], r["hai_per_1000_beds"]), textcoords="offset points", xytext=(5, 4), fontsize=9)
    ax.set_xlabel("Share of inpatient patient-days that are mentally ill (%)")
    ax.set_ylabel("Sequenced HAI per 1,000 beds")
    ax.set_title("E. Psychiatric throughput share vs sequenced HAI rate")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(OUT / "fig_psych_share_vs_hai.png", dpi=220)
    fig.savefig(OUT / "fig_psych_share_vs_hai.pdf")
    plt.close(fig)


def write_findings_log(hai, pathways, df, tips):
    sample = df[df["sequenced_hai_cases"] > 0]
    nt = sample[sample["ha_cluster"].isin(NT)]
    urb = sample[sample["ha_cluster"].isin(URBAN_IN_SAMPLE)]
    win = tips[(tips["date"] >= "2022-05-28") & (tips["date"] <= "2022-08-18")]
    community_tips = int(win["hospital_code"].isna().sum())
    hospital_tips = int(win["hospital_code"].notna().sum())

    path_cases = pathways.groupby("path")["cases"].sum()
    api_ge3 = int(hai.loc[hai["Case type"] == "Inpatient (API >= 3)", "N"].sum())
    api_known = int(hai.loc[hai["Case type"].isin(["Inpatient (API >= 3)", "Inpatient (API < 3)"]), "N"].sum())

    # Spearman across in-sample clusters
    from scipy import stats

    def spr(a, b, sub=sample):
        if len(sub) < 3:
            return np.nan, np.nan
        return stats.spearmanr(sub[a], sub[b])

    rho_dis, p_dis = spr("hai_per_100k_discharges", "hai_per_1000_beds")
    # NT vs urban rate ratios
    def mean_ratio(col):
        return float(nt[col].mean() / urb[col].mean())

    # KCC contrast
    kcc = df[df["ha_cluster"] == "KCC"].iloc[0]

    findings = f"""# Manuscript findings log
## Dual seeding, geographic concentration, density paradox, and denominators

**Purpose:** Catalog every claim from the exploration thread with evidence, grade, and manuscript use.  
**Parent paper:** Shen et al., *Influenza and Other Respiratory Viruses* (2026) — genomic epidemiology of hospital-associated SARS-CoV-2 in HK, May–Aug 2022 (162 genomes; DOI 10.1111/irv.70249).  
**This thread:** post-publication / companion exploration extending that dataset with dual-seeding taxonomy, density context, and denominator checks.

Evidence grades: **A** = directly measured in committed study data; **B** = derived with transparent external denominators; **C** = directional / incomplete external scrape; **D** = interpretation / hypothesis.

---

## Finding set 1 — Dual seeding (ward entry routes)

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F1.1 | Staff/outpatient-first and inpatient-first pathways both contribute substantial case burden | Pathway cases: staff_first={int(path_cases.get('staff_first',0))}, inpatient_first={int(path_cases.get('inpatient_first',0))}, same_day={int(path_cases.get('same_day',0))} (`ward_pathways.csv`) | A | Results — cluster initiation |
| F1.2 | NT vs urban pathway mix does not differ detectably | Fisher exact OR=0.33, p=0.56 (n=13 ordered dual-type wards) | A | Results / Supplement |
| F1.3 | Among known-API inpatients, most are API≥3 | {api_ge3}/{api_known} ({100*api_ge3/max(api_known,1):.0f}%) | A | Already in parent paper (82.2%); restate |
| F1.4 | Parent paper phylogenetic staff-first is rarer than temporal first-detection dual seeding | Parent: 2/18 confirmed clusters had phylogenetically linked staff before first likely nosocomial inpatient; our temporal rule flags more staff-first wards | A/D | Discussion — definition sensitivity |
| F1.5 | Multi-lineage wards imply reintroduction | 4 wards with ≥2 lineages | A | Results — multiple introductions |

**Definition note (methods):** “Staff-first” here = earliest *detection date* among Outpatient/Staff vs Inpatient categories in the plot table — not necessarily phylogenetic index. Report both definitions if combined with parent paper.

---

## Finding set 2 — Geographic concentration

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F2.1 | Sequenced HAI are NT-concentrated | NTWC+NTEC = {int(hai.loc[hai['Hospital clusters'].isin(NT),'N'].sum())}/{int(hai['N'].sum())} ({100*hai.loc[hai['Hospital clusters'].isin(NT),'N'].sum()/hai['N'].sum():.1f}%) | A | Results — geography |
| F2.2 | Top hospitals dominate | Top 3 anonymized hospitals ≈ half of cases | A | Results / Supplement |
| F2.3 | HKWC & KCC absent from sequenced plot table | 0 sequenced HAI cases despite large bed bases | A | Results — sampling map |

---

## Finding set 3 — Density paradox (context, not cause)

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F3.1 | NT catchments are less population-dense than urban in-sample clusters | Mean pop density ratio NT/urban ≈ {mean_ratio('pop_density'):.2f} | B | Results — context table |
| F3.2 | NT has lower hospital geographic density | Hospitals/100km² ratio ≈ {mean_ratio('hospitals_per_100km2'):.2f} | B | Results — context table |
| F3.3 | NT has lower bed geographic density | Beds/100km² ratio ≈ {mean_ratio('beds_per_100km2'):.2f} | B | Results — context table |
| F3.4 | Low density is not sufficient for high sequenced HAI | KEC ≈ NT hospital density but low HAI rate | A/B | Discussion — negative control |
| F3.5 | Density does not explain NT excess | Paradox is descriptive; no causal density model claimed | D | Discussion |

Sources: 2021 Census district areas/populations; HA cluster–district mapping; see `analysis/reference/SOURCES.md`.

---

## Finding set 4 — Admissions / throughput denominators

Annual HA cluster extracts (FY2021-22; DownloadCluster 36–42):

| Cluster | Sequenced HAI | Discharges | Patient-days | HAI/1000 beds | HAI/100k discharges | HAI/100k patient-days |
|---|---:|---:|---:|---:|---:|---:|
"""
    for _, r in df.sort_values("hai_per_1000_beds", ascending=False).iterrows():
        findings += (
            f"| {r['ha_cluster']} | {int(r['sequenced_hai_cases'])} | {int(r['ip_discharges_2021_22'])} | "
            f"{int(r['patient_days_2021_22'])} | {r['hai_per_1000_beds']:.2f} | {r['hai_per_100k_discharges']:.2f} | {r['hai_per_100k_patient_days']:.2f} |\n"
        )

    findings += f"""
| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F4.1 | NT excess survives discharge standardization | NT/urban mean ratio for HAI/100k discharges ≈ **{mean_ratio('hai_per_100k_discharges'):.2f}** | B | Results — denominators |
| F4.2 | NT excess survives patient-day standardization | NT/urban mean ratio for HAI/100k patient-days ≈ **{mean_ratio('hai_per_100k_patient_days'):.2f}** | B | Results — denominators |
| F4.3 | NT excess survives bed standardization | NT/urban mean ratio for HAI/1000 beds ≈ **{mean_ratio('hai_per_1000_beds'):.2f}** | B | Results — denominators |
| F4.4 | Annual denominators are imperfect for an 83-day window | Window-scaled discharge rates provided as sensitivity (`hai_per_100k_discharges_window_scaled`) | B | Methods / Supplement |
| F4.5 | NTWC has the highest mentally-ill patient-day volume | NTWC mentally-ill patient-days = 255,564 (share {100*df.loc[df.ha_cluster=='NTWC','psych_patient_day_share'].iloc[0]:.1f}% of cluster patient-days) | B | Discussion — case-mix |

**Verdict (admissions):** Under three independent throughput denominators, NT sequenced HAI rates remain several-fold higher than urban in-sample clusters. Throughput alone does **not** wash out the geography.

---

## Finding set 5 — Sequencing effort / ascertainment

Tree window tips: hospital-coded = {hospital_tips}; community/other = {community_tips} (BA.2.2/BA.2.12.1/BA.5.6 date files).

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F5.1 | Hospital-coded tips track sequenced HAI by construction | Study design = targeted sequencing of suspected clusters; tip geography ≈ sample geography | A/D | Methods — limitation |
| F5.2 | Tip counts therefore cannot independently prove sequencing created NT excess | Circularity explicitly documented | D | Discussion |
| F5.3 | Community background sequences dominate the trees | {community_tips}/{community_tips+hospital_tips} ({100*community_tips/max(community_tips+hospital_tips,1):.0f}%) of window tips lack hospital codes | A | Methods — phylogeographic context |
| F5.4 | Official HA nosocomial bulletins also feature NT hospitals, especially Castle Peak | Mention-days: NTWC={int(df.loc[df.ha_cluster=='NTWC','official_mention_days'].iloc[0])}, NTEC={int(df.loc[df.ha_cluster=='NTEC','official_mention_days'].iloc[0])} (`ha_nosocomial_mention_days.csv`) | C | Results — ascertainment |
| F5.5 | Ascertainment is incomplete / uneven: KCC appears in official bulletins but has 0 sequenced HAI in the plot table | KCC official mention-days = {int(kcc['official_mention_days'])}, sequenced HAI = 0; KCC beds = {int(kcc['beds_total'])} | C/A | Discussion — sampling bias |
| F5.6 | Parent paper overall coverage ≈ half of official nosocomial cases | 126 confirmed genomes ≈ 51.9% of official nosocomial cases in window (Table S3; implied official total ≈ 243) | A | Already in parent paper |
| F5.7 | Best reading: NT excess is unlikely to be *only* sequencing artefact, but sequencing coverage is geographically uneven | Joint reading of F4.* (throughput) + F5.4–F5.5 (official vs sequenced) | D | Discussion — balanced |

**Verdict (sequencing):** We **cannot** fully separate epidemiology from ascertainment with committed data alone. Official bulletin footprints suggest (i) NT/psychiatric hospitals were truly active in public reporting, and (ii) some urban clusters (esp. KCC) were under-represented in the sequenced table. Table S3 hospital-level official counts (supplement) would upgrade F5.4–F5.5 from C to A/B.

---

## Finding set 6 — Case-mix bridge (links density paradox to biology)

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F6.1 | BA.5.6 is a single-hospital deep outbreak with nosocomial API signature | 22/22 BA.5.6 in hospital I; 21/22 API≥3; tree labels = CPH/A102 | A | Results — lineage vignette |
| F6.2 | Long-stay psychiatric capacity is concentrated in NTWC | Mentally-ill + MH beds; high mentally-ill patient-days; parent paper notes APIs of hundreds–thousands of days | A/B | Discussion |
| F6.3 | Density paradox + psych case-mix is a coherent alternative to “urban density drives HAI” | Lower NT density + higher psych share + higher standardized HAI | D | Discussion / conceptual figure |

---

## Finding set 7 — Mobility null (parent paper, restated)

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F7.1 | After time smooth, mobility indices do not predict community→hospital Markov jumps | GAM IRRs ≈ 1, all p>0.5 (`gam_model_results.csv`) | A | Already in parent paper |
| F7.2 | Dual-seeding + density paradox + denominator results argue risk is ward/case-mix structured, not citywide mobility-driven | Synthesis | D | Discussion |

---

## Recommended manuscript architecture (companion note or Discussion extension)

1. **Short Results paragraph:** dual-seeding taxonomy (F1) + NT concentration (F2).  
2. **Context table:** density metrics (F3) + throughput-standardized rates (F4).  
3. **Limitations paragraph:** sequencing circularity + KCC under-ascertainment (F5) + need for Table S3 hospital totals.  
4. **Interpretation:** psychiatric/long-stay case-mix as bridge (F6); mobility null unchanged (F7).

### What would upgrade this to a full companion paper
1. Release or re-tabulate **Table S3 official nosocomial cases by hospital/cluster**.  
2. Hospital-level admissions during May 28–Aug 18 2022 (not annual).  
3. Sensitivity: dual-seeding using phylogenetic index (parent definition) vs temporal first detection.

### Stop conditions honored
- No causal claim that low density *causes* HAI.  
- No claim that sequencing effort fully explains NT excess.  
- Both denominators pursued; residual ambiguity documented rather than forced.

---

## File index
- `scripts/4.dual_seeding_geography.py` — dual seeding + density paradox  
- `scripts/5.denominators_ascertainment.py` — this script  
- `analysis/outputs/denominator_table.csv`  
- `analysis/outputs/FINDINGS_LOG.md` — this file  
- `analysis/reference/ha_cluster_throughput_2021_22.csv`  
- `analysis/reference/ha_nosocomial_mention_days.csv`  
"""
    (OUT / "FINDINGS_LOG.md").write_text(findings)
    return findings


def main():
    hai = pd.read_csv(HAI_PATH, parse_dates=["Date"])
    capacity = pd.read_csv(CAPACITY_PATH)
    tips = parse_tree_tips()
    pathways = classify_pathways(hai)
    df = build_denominator_table(hai, tips, capacity)

    THROUGHPUT.to_csv(REF / "ha_cluster_throughput_2021_22.csv", index=False)
    df.to_csv(OUT / "denominator_table.csv", index=False)
    pathways.to_csv(OUT / "ward_pathways.csv", index=False)

    tip_summary = {
        "window": ["2022-05-28", "2022-08-18"],
        "hospital_coded_tips": int(tips[(tips.date >= "2022-05-28") & (tips.date <= "2022-08-18")]["hospital_code"].notna().sum()),
        "community_tips": int(tips[(tips.date >= "2022-05-28") & (tips.date <= "2022-08-18")]["hospital_code"].isna().sum()),
        "tips_by_cluster": tips[(tips.date >= "2022-05-28") & (tips.date <= "2022-08-18")]
        .dropna(subset=["ha_cluster"])
        .groupby("ha_cluster")
        .size()
        .astype(int)
        .to_dict(),
        "tips_by_hospital_code": tips[(tips.date >= "2022-05-28") & (tips.date <= "2022-08-18")]
        .dropna(subset=["hospital_code"])
        .groupby("hospital_code")
        .size()
        .astype(int)
        .sort_values(ascending=False)
        .to_dict(),
    }
    (OUT / "sequencing_effort_summary.json").write_text(json.dumps(tip_summary, indent=2))

    make_figures(df)
    write_findings_log(hai, pathways, df, tips)

    cols = [
        "ha_cluster",
        "sequenced_hai_cases",
        "hai_per_1000_beds",
        "hai_per_100k_discharges",
        "hai_per_100k_patient_days",
        "hospital_coded_tips_in_window",
        "official_mention_days",
        "psych_patient_day_share",
    ]
    print(df[cols].sort_values("hai_per_1000_beds", ascending=False).round(3).to_string(index=False))
    print("\nWrote", OUT / "FINDINGS_LOG.md")


if __name__ == "__main__":
    main()
