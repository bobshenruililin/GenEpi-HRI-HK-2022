# Data strategy after peer review (no Stage 2+)

Author: Shen Ruililin — Laidlaw Scholars Programme, HKU  
Date: 2026-07-11

## Short answers

| Question | Answer |
|---|---|
| Will GitHub / Kaggle help? | **No** for the missing Stage-2 denominators. |
| Apply to Hospital Authority? | **Only if** you want a longer project; not needed to finish this note. |
| Work with what we have? | **Yes** — that is the right call for a Laidlaw secondary note. |

## What is already worthy (keep and defend)

1. **Ascertainment-aware geographic reanalysis** of Gu et al.’s sequenced sample  
   NT concentration after *inpatient-compatible* throughput scaling + explicit audit that some HA clusters are in official bulletins but absent from the genomic table.
2. **Descriptive first-detection-order taxonomy** on the 18 confirmed wards, clearly *not* equated to seeding/direction.
3. **Contribution matrix** separating Gu et al. findings from new secondary analyses.
4. **Honest estimands**: inpatient rates ≠ mixed HAI incidence; unrepresented clusters = NA; KEC = comparator; psych case-mix = candidate only.

That package is suitable as a **focused secondary research note / Laidlaw companion**, not a top-journal mechanistic paper.

## What peer-review items we can still do without restricted data

| Item | Status / action |
|---|---|
| Gu et al. attribution; first-detection language; confirmed primary cohort | Done (PR #3) |
| Separate inpatient vs staff numerators; NA not zero | Done |
| Pooled + mean ratios; Fisher CI; LOO; BA.5.6 exclusion; tie windows | Done |
| Contribution matrix; decision log outside MS | Done |
| Update FY2022–23 public HA throughput | Feasible polish (data.gov.hk / HA portal) |
| Hospital-level open throughput without name crosswalk | Feasible polish if open JSON used carefully |
| Density land-definition sensitivity | Feasible polish |
| Formal selection model / IPW / sampling fractions | **Blocked** without HA/parent totals |
| Visitation causal design | **Blocked** (and out of scope) |
| Staff FTE incidence | **Blocked** |

## External search: what you would find

- **GitHub:** Gu/Poon lab repos and other HK genomic projects → sequences/community epi, **not** HA nosocomial registries or sampling fractions.
- **Kaggle / open ML hubs:** no substitute for Table S3 or sequencing inclusion.
- **data.gov.hk / HA open data:** beds, discharges, patient-days by cluster/hospital → **better denominators only**.
- **info.gov.hk press:** presence/mention proxy (already used); cannot rebuild cumulative official counts.
- **GISAID:** community backdrop already used by Gu et al.; circular if used as “sequencing effort” for the same hospital tips.

## HA research-data application

Worth considering **only** if you want a follow-on paper months later. Typical ask:

- official nosocomial counts by hospital/cluster for 28 May–18 Aug 2022;
- specimens eligible / submitted / sequenced / included;
- optional staff headcount or FTE by cluster.

For the current Laidlaw note: **do not wait** on this.

## Recommended “further work” sentence (manuscript)

Already added to `manuscript/main.tex`: position the note as ascertainment-aware sample geography + descriptive first-detection patterns; list public-denominator polish as feasible; state that selection-adjusted incidence needs HA/parent access.
