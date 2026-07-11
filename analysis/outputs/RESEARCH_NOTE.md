# Dual seeding under a geographic paradox
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

Ward-level first-detection classification (n = 29 wards):

| Pathway | Wards | Cases |
|---|---:|---:|
| Staff/outpatient first | 6 | 48 |
| Inpatient first | 7 | 54 |
| Same day | 3 | 24 |
| Staff/OP only | 5 | 15 |
| Inpatient only | 8 | 21 |

Among wards with both staff and inpatient signal and a clear order, staff-first and inpatient-first are **balanced in case burden** (~48 vs ~54 cases).  
Fisher exact test (NT vs urban × staff-first vs inpatient-first): OR = 0.33, p = 0.559, n = 13 wards.  
**Interpretation:** no evidence that NT and urban clusters use different dominant pathways; dual seeding is a **shared mechanism** inside a geographically skewed burden.

Among known-API inpatients, **74/90 (82%)** have API ≥ 3 — consistent with substantial true nosocomial acquisition once infection is inside the ward.

---

## Finding 2 — Geographic concentration

- NTWC + NTEC: **134/162 (82.7%)** of sequenced HAI cases.
- Top 3 anonymized hospitals: ~half of all cases (see script tables).
- HKWC and KCC appear in HA capacity data but contribute **0** cases to this sequenced plot table — concentration is even starker against the full HA map.

---

## Finding 3 — The density paradox (quantified)

Among clusters **present in the HAI table**:

| Metric (mean of clusters) | New Territories (NTWC, NTEC) | Urban in sample (KWC, KEC, HKEC) | Ratio NT/Urban |
|---|---:|---:|---:|
| Population density (per km²) | 4,420 | 16,131 | 0.27 |
| Hospitals per 100 km² | 2.09 | 10.54 | 0.20 |
| Beds per 100 km² | 1,780 | 6,223 | 0.29 |
| Sequenced HAI per 1,000 beds | 13.63 | 2.39 | 5.71 |

**Read carefully:** NT is less dense on all three geographic metrics, yet has ~5.7× the bed-standardized sequenced HAI rate.  
This is the opposite of a naive “dense city → more hospital outbreaks” story.

**Negative control inside the paradox:** KEC also has low hospital geographic density (~2 per 100 km², similar to NT) but only 1.7 sequenced HAI cases per 1,000 beds. Low density is therefore **not sufficient** for high observed burden — the NT excess needs NT-specific structure (case-mix / long-stay capacity / sequencing effort / IPC), not density alone.

Cluster-level detail is in `cluster_burden_standardized.csv`.

**Psychiatric capacity note:** NTWC alone reports **1,176 mentally-ill + 520 mentally-handicapped beds** (HA 2021–22). Long-stay settings inflate API and can sustain ward transmission — a plausible partial link to high-API / BA.5.6 patterns, not yet hospital-resolved in the anonymized table.

---

## Finding 4 — Other angles that fold into the thesis

1. **Visit relaxation as clock, not proven cause.** 2/162 sequenced HAI cases precede 31 May 2022. Community incidence was also rising; density-paradox + dual-seeding argue against a simple mobility story (consistent with the null GAM on mobility predictors).
2. **Community alignment is weekly, not daily.** Corr(local, HAI) ≈ 0.27; corr(7-day means) ≈ 0.51. Introductions track the Phase-2 wave envelope.
3. **BA.5.6 is a single-hospital deep outbreak** (22 cases, all one hospital, almost all API ≥ 3) sitting inside NTWC — extreme version of “introduction then sustained nosocomial spread.”
4. **Multi-lineage wards (n = 4)** imply reintroduction into the same unit — supports many independent seeds, not one clonal citywide hospital epidemic.
5. **Starred (non-official) wards hold 36 cases (22%).** Official cluster lists understate hospital-linked transmission; dual-seeding estimates are conservative if anything.

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
