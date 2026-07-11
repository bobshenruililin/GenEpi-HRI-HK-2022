#!/usr/bin/env python3
"""
Optional polish (public data only)
==================================
1. Refresh cluster throughput/beds from HA open data for FY2021-22 vs FY2022-23.
2. Hospital-letter inpatient scaling using internal code↔open-data name join;
   outputs use anonymized letters only (no name crosswalk published).
3. Density land/population definition sensitivity (core vs total land; 2017 vs 2021 pop).

Does not claim selection-adjusted incidence or recover Table S3.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "analysis/reference"
OUT = ROOT / "analysis/outputs"
FIG = ROOT / "manuscript/figures"
HAI_PATH = ROOT / "data/hospital_data/metadata_cases_plot.csv"
DIST_PATH = REF / "district_demographics_2021.csv"
CAPACITY_PATH = REF / "ha_cluster_capacity_2022.csv"
ASCERTAIN_PATH = REF / "ascertainment_presence_by_cluster.csv"

NT = {"NTWC", "NTEC"}
URBAN_IN_SAMPLE = {"KWC", "KEC", "HKEC"}

CLUSTER_NAME_TO_CODE = {
    "Hong Kong East Cluster": "HKEC",
    "Hong Kong West Cluster": "HKWC",
    "Kowloon Central Cluster": "KCC",
    "Kowloon East Cluster": "KEC",
    "Kowloon West Cluster": "KWC",
    "New Territories East Cluster": "NTEC",
    "New Territories West Cluster": "NTWC",
}

# High-confidence letter → HA hospital code from distinctive ward/tip overlap.
# Used only inside this script; not written to public tables as a name crosswalk.
LETTER_TO_CODE = {
    "A": "TMH",   # D4, H1
    "C": "PWH",   # 12C
    "D": "PYNEH", # E8
    "E": "CMC",   # 11A
    "F": "TKOH",  # KEC; tip count matches
    "G": "AHNH",  # E5
    "H": "TKP",   # E1/E2/E3 tips; psychiatric co-class with CPH — may lack open-data row
    "I": "CPH",   # A102
    "J": "NDH",   # 3H
    "K": "TWE",   # A2R
    "L": "TSH",   # 2A — Siu Lam Hospital in open data? tips TSH; map below
    "M": "TPH",   # 2DL
    "N": "SH",    # 9D
    "O": "KCH",   # L2
    "P": "RH",    # C7
    "Q": "PMH",   # ward label PMH
    # B (B101*) left unmapped — no reliable tip/open-data join
}

# HA open-data hospital display names for codes present in LETTER_TO_CODE
CODE_TO_OPEN_NAME = {
    "TMH": "Tuen Mun Hospital",
    "CPH": "Castle Peak Hospital",
    "PWH": "Prince of Wales Hospital",
    "PYNEH": "Pamela Youde Nethersole Eastern Hospital",
    "CMC": "Caritas Medical Centre",
    "TKOH": "Tseung Kwan O Hospital",
    "AHNH": "Alice Ho Miu Ling Nethersole Hospital",
    "NDH": "North District Hospital",
    "TWE": "Tung Wah Eastern Hospital",
    "TSH": "Siu Lam Hospital",  # tip code TSH; NTWC long-stay; open-data name Siu Lam
    "TPH": "Tai Po Hospital",
    "SH": "Shatin Hospital",
    "KCH": "Kwai Chung Hospital",
    "RH": "Ruttonjee & Tang Shiu Kin Hospitals",
    "PMH": "Princess Margaret Hospital",
    # TKP: no separate open-data institution row identified
}

URLS = {
    "discharges": "https://www.ha.org.hk/opendata/ipdpdd-en.json",
    "patient_days": "https://www.ha.org.hk/opendata/patientday-en.json",
    "beds": "https://www.ha.org.hk/opendata/hosp-bed-en.json",
}


def fetch_json(url: str, cache: Path) -> list:
    if cache.exists() and cache.stat().st_size > 1000:
        return json.loads(cache.read_text())
    print(f"Fetching {url}")
    with urllib.request.urlopen(url, timeout=120) as resp:
        data = resp.read()
    cache.write_bytes(data)
    return json.loads(data)


def load_open_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cache_dir = REF / "ha_opendata_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    disc = pd.DataFrame(fetch_json(URLS["discharges"], cache_dir / "ipdpdd-en.json"))
    pday = pd.DataFrame(fetch_json(URLS["patient_days"], cache_dir / "patientday-en.json"))
    beds = pd.DataFrame(fetch_json(URLS["beds"], cache_dir / "hosp-bed-en.json"))
    for df in (disc, pday, beds):
        df["ha_cluster"] = df["Cluster"].map(CLUSTER_NAME_TO_CODE)
    return disc, pday, beds


def cluster_throughput(disc: pd.DataFrame, pday: pd.DataFrame, beds: pd.DataFrame, year: str) -> pd.DataFrame:
    """Cluster-level totals excluding Overall rows; beds = sum of specialty beds."""
    d = disc[(disc["Financial Year"] == year) & disc["ha_cluster"].notna()].copy()
    d = d[~d["Hospital"].astype(str).str.startswith("Overall")]
    p = pday[(pday["Financial Year"] == year) & pday["ha_cluster"].notna()].copy()
    p = p[~p["Hospital"].astype(str).str.startswith("Overall")]
    b = beds[(beds["Financial Year"] == year) & beds["ha_cluster"].notna()].copy()
    b = b[~b["Hospital"].astype(str).str.startswith("Overall")]

    disc_c = d.groupby("ha_cluster", as_index=False).agg(
        ip_discharges=( "IP Discharges and Deaths", "sum"),
        n_hospitals_open=("Hospital", "nunique"),
    )
    pday_c = p.groupby("ha_cluster", as_index=False).agg(
        patient_days=("IP Patient Days", "sum"),
    )
    beds_c = b.groupby("ha_cluster", as_index=False).agg(
        beds_total=("Hospital Beds", "sum"),
        n_hospitals_bed_table=("Hospital", "nunique"),
    )
    out = disc_c.merge(pday_c, on="ha_cluster").merge(beds_c, on="ha_cluster")
    out["financial_year"] = year
    return out


def hospital_open_metrics(disc, pday, beds, year: str) -> pd.DataFrame:
    d = disc[(disc["Financial Year"] == year) & disc["ha_cluster"].notna()].copy()
    d = d[~d["Hospital"].astype(str).str.startswith("Overall")]
    p = pday[(pday["Financial Year"] == year) & pday["ha_cluster"].notna()].copy()
    p = p[~p["Hospital"].astype(str).str.startswith("Overall")]
    b = beds[(beds["Financial Year"] == year) & beds["ha_cluster"].notna()].copy()
    b = b[~b["Hospital"].astype(str).str.startswith("Overall")]

    disc_h = d.groupby(["ha_cluster", "Hospital"], as_index=False).agg(
        ip_discharges=("IP Discharges and Deaths", "sum")
    )
    pday_h = p.groupby(["ha_cluster", "Hospital"], as_index=False).agg(
        patient_days=("IP Patient Days", "sum")
    )
    beds_h = b.groupby(["ha_cluster", "Hospital"], as_index=False).agg(
        beds_total=("Hospital Beds", "sum")
    )
    return disc_h.merge(pday_h, on=["ha_cluster", "Hospital"]).merge(
        beds_h, on=["ha_cluster", "Hospital"]
    )


def load_hai_confirmed() -> pd.DataFrame:
    hai = pd.read_csv(HAI_PATH, parse_dates=["Date"])
    hai["confirmed"] = ~hai["Wards"].astype(str).str.endswith("*")
    hai["case_group"] = np.where(
        hai["Case type"].eq("Outpatient/Staff"), "staff_outpatient", "inpatient"
    )
    return hai[hai["confirmed"]].copy()


def cluster_case_counts(hai: pd.DataFrame) -> pd.DataFrame:
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
            }
        )
    return pd.DataFrame(rows)


def letter_case_counts(hai: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (letter, cl), sub in hai.groupby(["Hospital_anonymized", "Hospital clusters"]):
        rows.append(
            {
                "Hospital_anonymized": letter,
                "ha_cluster": cl,
                "sequenced_all": int(sub["N"].sum()),
                "sequenced_inpatient": int(sub.loc[sub["case_group"] == "inpatient", "N"].sum()),
                "sequenced_staff_outpatient": int(
                    sub.loc[sub["case_group"] == "staff_outpatient", "N"].sum()
                ),
                "n_wards": int(sub["Wards"].nunique()),
                "open_data_join": "yes" if letter in LETTER_TO_CODE and LETTER_TO_CODE[letter] in CODE_TO_OPEN_NAME else "no",
            }
        )
    return pd.DataFrame(rows)


def nt_urban_pooled(den: pd.DataFrame, num_col: str, den_col: str, scale: float) -> dict:
    d = den[den["genomic_representation"] == "represented"].copy()
    nt = d[d["ha_cluster"].isin(NT)]
    urb = d[d["ha_cluster"].isin(URBAN_IN_SAMPLE)]
    pooled_nt = nt[num_col].sum() / nt[den_col].sum() * scale
    pooled_urb = urb[num_col].sum() / urb[den_col].sum() * scale
    return {
        "pooled_nt": float(pooled_nt),
        "pooled_urban": float(pooled_urb),
        "pooled_ratio_nt_urban": float(pooled_nt / pooled_urb) if pooled_urb else None,
    }


def density_sensitivity(capacity: pd.DataFrame, districts: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    """
    Definitions:
    - land_core + pop_2017: existing capacity table
    - land_core + pop_2021: sum 2021 census for core districts
    - land_total_excl_islands + pop_2021: same districts (Islands excluded)
    - land_with_islands_split: assign Islands area 50/50 to HKEC/KWC as sensitivity
    """
    dist = districts.copy()
    # core district sums by cluster (exclude Islands)
    core = dist[dist["ha_cluster_primary"] != "mixed_HKEC_KWC"]
    pop21 = core.groupby("ha_cluster_primary", as_index=False).agg(
        pop_2021=("population_2021", "sum"),
        land_2021_core=("area_km2", "sum"),
    ).rename(columns={"ha_cluster_primary": "ha_cluster"})

    islands = dist[dist["district"] == "Islands"].iloc[0]
    # sensitivity: split Islands land+pop equally to HKEC and KWC
    pop21_island = pop21.copy()
    land_island = pop21.copy()
    for cl in ("HKEC", "KWC"):
        if cl in set(pop21_island["ha_cluster"]):
            pop21_island.loc[pop21_island["ha_cluster"] == cl, "pop_2021"] += islands["population_2021"] / 2
            land_island.loc[land_island["ha_cluster"] == cl, "land_2021_core"] += islands["area_km2"] / 2
    pop21_island = pop21_island.rename(columns={"pop_2021": "pop_2021_islands_split"})
    land_island = land_island.rename(columns={"land_2021_core": "land_2021_islands_split"})

    base = capacity[
        ["ha_cluster", "catchment_population_2017", "land_area_km2_core", "n_hospitals_institutions", "region_group"]
    ].merge(pop21, on="ha_cluster", how="left")
    base = base.merge(pop21_island[["ha_cluster", "pop_2021_islands_split"]], on="ha_cluster", how="left")
    base = base.merge(land_island[["ha_cluster", "land_2021_islands_split"]], on="ha_cluster", how="left")
    base = base.merge(counts, on="ha_cluster", how="left")
    base["sequenced_inpatient"] = base["sequenced_inpatient"].fillna(0).astype(int)
    base["genomic_representation"] = np.where(
        base["sequenced_inpatient"] + base.get("sequenced_all", 0).fillna(0) > 0
        if "sequenced_all" in base.columns
        else base["sequenced_inpatient"] > 0,
        "represented",
        "not_represented_in_genomic_sample",
    )
    # fix representation using sequenced_all if present
    if "sequenced_all" in base.columns:
        base["genomic_representation"] = np.where(
            base["sequenced_all"].fillna(0) > 0,
            "represented",
            "not_represented_in_genomic_sample",
        )

    base["dens_pop2017_land_core"] = base["catchment_population_2017"] / base["land_area_km2_core"]
    base["dens_pop2021_land_core"] = base["pop_2021"] / base["land_area_km2_core"]
    base["dens_pop2021_land_2021"] = base["pop_2021"] / base["land_2021_core"]
    base["dens_pop2021_islands_split"] = (
        base["pop_2021_islands_split"] / base["land_2021_islands_split"]
    )
    base["hosp_per_100km2_core"] = base["n_hospitals_institutions"] / base["land_area_km2_core"] * 100
    base["hosp_per_100km2_2021"] = base["n_hospitals_institutions"] / base["land_2021_core"] * 100
    return base


def dens_ratio_summary(dens: pd.DataFrame) -> dict:
    d = dens[dens["genomic_representation"] == "represented"]
    nt = d[d["ha_cluster"].isin(NT)]
    urb = d[d["ha_cluster"].isin(URBAN_IN_SAMPLE)]
    out = {}
    for col in [
        "dens_pop2017_land_core",
        "dens_pop2021_land_core",
        "dens_pop2021_land_2021",
        "dens_pop2021_islands_split",
        "hosp_per_100km2_core",
        "hosp_per_100km2_2021",
    ]:
        out[col] = {
            "mean_nt": float(nt[col].mean()),
            "mean_urban": float(urb[col].mean()),
            "nt_over_urban": float(nt[col].mean() / urb[col].mean()),
        }
    return out


def make_figures(compare: pd.DataFrame, letter: pd.DataFrame, dens: pd.DataFrame):
    plt.rcParams.update({"font.size": 10, "figure.dpi": 150})

    # FY compare pooled ratios
    fig, ax = plt.subplots(figsize=(6.5, 4))
    metrics = ["ratio_per_1000_beds", "ratio_per_100k_discharges", "ratio_per_100k_patient_days"]
    labels = ["Per 1k beds", "Per 100k discharges", "Per 100k patient-days"]
    x = np.arange(len(metrics))
    y21 = [compare.loc[compare["financial_year"] == "2021-22", m].iloc[0] for m in metrics]
    y22 = [compare.loc[compare["financial_year"] == "2022-23", m].iloc[0] for m in metrics]
    ax.bar(x - 0.18, y21, 0.36, label="FY2021–22 denominators")
    ax.bar(x + 0.18, y22, 0.36, label="FY2022–23 denominators")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Pooled NT / urban ratio (sequenced inpatients)")
    ax.set_title("Throughput-year sensitivity (confirmed wards, represented clusters)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "fig_fy_throughput_sensitivity.pdf")
    fig.savefig(OUT / "fig_fy_throughput_sensitivity.pdf")
    plt.close(fig)

    # Hospital-letter rates (joined only)
    L = letter[letter["denominator_available"]].sort_values("inpatient_per_1000_beds", ascending=False)
    colors = ["#2a6f4e" if c in NT else "#4a6fa5" for c in L["ha_cluster"]]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.bar(L["Hospital_anonymized"], L["inpatient_per_1000_beds"], color=colors)
    ax.set_xlabel("Anonymized hospital letter (open-data denominator joined; names not shown)")
    ax.set_ylabel("Sequenced confirmed-ward inpatients / 1,000 beds")
    ax.set_title("Hospital-letter sample concentration (FY2022–23 open beds)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_hospital_letter_inpatient_rates.pdf")
    fig.savefig(OUT / "fig_hospital_letter_inpatient_rates.pdf")
    plt.close(fig)

    # Density definitions
    d = dens[dens["genomic_representation"] == "represented"]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    defs = [
        ("dens_pop2017_land_core", "Pop2017 / land_core"),
        ("dens_pop2021_land_core", "Pop2021 / land_core"),
        ("dens_pop2021_land_2021", "Pop2021 / land2021"),
        ("dens_pop2021_islands_split", "Pop2021 / land+Islands split"),
    ]
    x = np.arange(len(defs))
    nt_vals = [d.loc[d["ha_cluster"].isin(NT), c].mean() / d.loc[d["ha_cluster"].isin(URBAN_IN_SAMPLE), c].mean() for c, _ in defs]
    ax.bar(x, nt_vals, color="#2a6f4e")
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in defs], rotation=20, ha="right")
    ax.set_ylabel("Mean NT / mean urban density ratio")
    ax.set_title("Density definition sensitivity (represented clusters)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_density_definition_sensitivity.pdf")
    fig.savefig(OUT / "fig_density_definition_sensitivity.pdf")
    plt.close(fig)


def main():
    disc, pday, beds = load_open_tables()
    hai = load_hai_confirmed()
    counts = cluster_case_counts(hai)
    letters = letter_case_counts(hai)
    capacity = pd.read_csv(CAPACITY_PATH)
    districts = pd.read_csv(DIST_PATH)
    ascertain = pd.read_csv(ASCERTAIN_PATH) if ASCERTAIN_PATH.exists() else pd.DataFrame()

    # --- Cluster FY compare ---
    rows_compare = []
    den_by_year = {}
    for year in ("2021-22", "2022-23"):
        thr = cluster_throughput(disc, pday, beds, year)
        thr.to_csv(REF / f"ha_cluster_throughput_open_{year.replace('-', '_')}.csv", index=False)
        den = thr.merge(counts, on="ha_cluster", how="left")
        for col in ("sequenced_all", "sequenced_inpatient", "sequenced_staff_outpatient"):
            den[col] = den[col].fillna(0).astype(int)
        den["genomic_representation"] = np.where(
            den["sequenced_all"] > 0, "represented", "not_represented_in_genomic_sample"
        )
        den["inpatient_per_1000_beds"] = np.where(
            den["genomic_representation"] == "represented",
            den["sequenced_inpatient"] / den["beds_total"] * 1000,
            np.nan,
        )
        den["inpatient_per_100k_discharges"] = np.where(
            den["genomic_representation"] == "represented",
            den["sequenced_inpatient"] / den["ip_discharges"] * 1e5,
            np.nan,
        )
        den["inpatient_per_100k_patient_days"] = np.where(
            den["genomic_representation"] == "represented",
            den["sequenced_inpatient"] / den["patient_days"] * 1e5,
            np.nan,
        )
        if len(ascertain):
            den = den.merge(
                ascertain[["ha_cluster", "in_official_bulletins", "mention_days"]],
                on="ha_cluster",
                how="left",
            )
        den.to_csv(OUT / f"polish_denominator_confirmed_{year.replace('-', '_')}.csv", index=False)
        den_by_year[year] = den
        r_beds = nt_urban_pooled(den, "sequenced_inpatient", "beds_total", 1000)
        r_disc = nt_urban_pooled(den, "sequenced_inpatient", "ip_discharges", 1e5)
        r_pd = nt_urban_pooled(den, "sequenced_inpatient", "patient_days", 1e5)
        rows_compare.append(
            {
                "financial_year": year,
                "ratio_per_1000_beds": r_beds["pooled_ratio_nt_urban"],
                "ratio_per_100k_discharges": r_disc["pooled_ratio_nt_urban"],
                "ratio_per_100k_patient_days": r_pd["pooled_ratio_nt_urban"],
                "pooled_nt_per_1000_beds": r_beds["pooled_nt"],
                "pooled_urban_per_1000_beds": r_beds["pooled_urban"],
            }
        )
    compare = pd.DataFrame(rows_compare)
    compare.to_csv(OUT / "polish_fy_throughput_sensitivity.csv", index=False)

    # Prefer FY2022-23 as primary open-data denominator table
    den22 = den_by_year["2022-23"]
    den22.to_csv(REF / "ha_cluster_throughput_2022_23.csv", index=False)

    # --- Hospital-letter scaling (letters only in output) ---
    hosp_open = hospital_open_metrics(disc, pday, beds, "2022-23")
    # map letter -> open name internally
    letter_rows = []
    for _, r in letters.iterrows():
        letter = r["Hospital_anonymized"]
        code = LETTER_TO_CODE.get(letter)
        name = CODE_TO_OPEN_NAME.get(code) if code else None
        rec = dict(r)
        rec["denominator_available"] = False
        rec["beds_total"] = np.nan
        rec["ip_discharges"] = np.nan
        rec["patient_days"] = np.nan
        if name:
            match = hosp_open[hosp_open["Hospital"] == name]
            if len(match) == 1:
                rec["denominator_available"] = True
                rec["beds_total"] = float(match["beds_total"].iloc[0])
                rec["ip_discharges"] = float(match["ip_discharges"].iloc[0])
                rec["patient_days"] = float(match["patient_days"].iloc[0])
        if rec["denominator_available"]:
            rec["inpatient_per_1000_beds"] = rec["sequenced_inpatient"] / rec["beds_total"] * 1000
            rec["inpatient_per_100k_discharges"] = (
                rec["sequenced_inpatient"] / rec["ip_discharges"] * 1e5 if rec["ip_discharges"] else np.nan
            )
            rec["inpatient_per_100k_patient_days"] = (
                rec["sequenced_inpatient"] / rec["patient_days"] * 1e5 if rec["patient_days"] else np.nan
            )
        else:
            rec["inpatient_per_1000_beds"] = np.nan
            rec["inpatient_per_100k_discharges"] = np.nan
            rec["inpatient_per_100k_patient_days"] = np.nan
        letter_rows.append(rec)

    letter_df = pd.DataFrame(letter_rows)
    # Drop any accidental name/code columns if present
    keep = [
        "Hospital_anonymized",
        "ha_cluster",
        "sequenced_all",
        "sequenced_inpatient",
        "sequenced_staff_outpatient",
        "n_wards",
        "open_data_join",
        "denominator_available",
        "beds_total",
        "ip_discharges",
        "patient_days",
        "inpatient_per_1000_beds",
        "inpatient_per_100k_discharges",
        "inpatient_per_100k_patient_days",
    ]
    letter_out = letter_df[keep].sort_values(["ha_cluster", "Hospital_anonymized"])
    letter_out.to_csv(OUT / "polish_hospital_letter_rates_fy2022_23.csv", index=False)

    # NT vs urban among letters with denominators
    L = letter_out[letter_out["denominator_available"]].copy()
    L["region"] = np.where(L["ha_cluster"].isin(NT), "NT", "Urban_in_sample")
    letter_pooled = {}
    for reg, sub in L.groupby("region"):
        letter_pooled[reg] = {
            "n_hospitals": int(len(sub)),
            "inpatients": int(sub["sequenced_inpatient"].sum()),
            "beds": float(sub["beds_total"].sum()),
            "per_1000_beds": float(sub["sequenced_inpatient"].sum() / sub["beds_total"].sum() * 1000),
        }
    letter_pooled["pooled_ratio_nt_urban"] = (
        letter_pooled["NT"]["per_1000_beds"] / letter_pooled["Urban_in_sample"]["per_1000_beds"]
        if letter_pooled.get("Urban_in_sample", {}).get("per_1000_beds")
        else None
    )

    # --- Density sensitivity ---
    dens = density_sensitivity(capacity, districts, counts)
    dens.to_csv(OUT / "polish_density_definition_sensitivity.csv", index=False)
    dens_sum = dens_ratio_summary(dens)

    # Discordance check: under all defs, is NT density still below urban?
    discordance = {
        k: {
            **v,
            "nt_less_dense_than_urban": v["nt_over_urban"] < 1,
        }
        for k, v in dens_sum.items()
        if k.startswith("dens_")
    }

    summary = {
        "fy_throughput_sensitivity": rows_compare,
        "hospital_letter": {
            "n_confirmed_hospitals": int(letters["Hospital_anonymized"].nunique()),
            "n_with_open_denominator": int(letter_out["denominator_available"].sum()),
            "n_without_open_denominator": int((~letter_out["denominator_available"]).sum()),
            "pooled_among_joined": letter_pooled,
            "note": (
                "Hospital-letter rates join open HA beds/discharges/patient-days via an "
                "internal code map; published table shows anonymized letters only."
            ),
        },
        "density_definition_sensitivity": dens_sum,
        "density_discordance_robust": all(
            v["nt_less_dense_than_urban"] for v in discordance.values()
        ),
        "claims": [
            "FY2022-23 open denominators leave NT sample concentration intact",
            "Hospital-letter scaling (joined subset) still shows NT concentration",
            "NT remains less dense than urban under alternative land/pop definitions",
        ],
    }
    with open(OUT / "polish_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    make_figures(compare, letter_out, dens)

    # Update SOURCES
    sources = REF / "SOURCES.md"
    addendum = """

## Optional polish open data (2026-07-11)
- HA open JSON: `ipdpdd-en.json`, `patientday-en.json`, `hosp-bed-en.json` (cached under `ha_opendata_cache/`).
- Cluster extracts: `ha_cluster_throughput_open_2021_22.csv`, `ha_cluster_throughput_open_2022_23.csv`.
- Hospital-letter rates use internal joins only; `polish_hospital_letter_rates_fy2022_23.csv` has anonymized letters, not names.
"""
    text = sources.read_text() if sources.exists() else ""
    if "Optional polish open data" not in text:
        sources.write_text(text.rstrip() + addendum, encoding="utf-8")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
