# Manuscript findings log
## Dual seeding, geographic concentration, density paradox, and denominators

**Purpose:** Catalog every claim from the exploration thread with evidence, grade, and manuscript use.  
**Parent paper:** Shen et al., *Influenza and Other Respiratory Viruses* (2026) — genomic epidemiology of hospital-associated SARS-CoV-2 in HK, May–Aug 2022 (162 genomes; DOI 10.1111/irv.70249).  
**This thread:** post-publication / companion exploration extending that dataset with dual-seeding taxonomy, density context, and denominator checks.

Evidence grades: **A** = directly measured in committed study data; **B** = derived with transparent external denominators; **C** = directional / incomplete external scrape; **D** = interpretation / hypothesis.

---

## Finding set 1 — Dual seeding (ward entry routes)

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F1.1 | Staff/outpatient-first and inpatient-first pathways both contribute substantial case burden | Pathway cases: staff_first=48, inpatient_first=54, same_day=24 (`ward_pathways.csv`) | A | Results — cluster initiation |
| F1.2 | NT vs urban pathway mix does not differ detectably | Fisher exact OR=0.33, p=0.56 (n=13 ordered dual-type wards) | A | Results / Supplement |
| F1.3 | Among known-API inpatients, most are API≥3 | 74/90 (82%) | A | Already in parent paper (82.2%); restate |
| F1.4 | Parent paper phylogenetic staff-first is rarer than temporal first-detection dual seeding | Parent: 2/18 confirmed clusters had phylogenetically linked staff before first likely nosocomial inpatient; our temporal rule flags more staff-first wards | A/D | Discussion — definition sensitivity |
| F1.5 | Multi-lineage wards imply reintroduction | 4 wards with ≥2 lineages | A | Results — multiple introductions |

**Definition note (methods):** “Staff-first” here = earliest *detection date* among Outpatient/Staff vs Inpatient categories in the plot table — not necessarily phylogenetic index. Report both definitions if combined with parent paper.

---

## Finding set 2 — Geographic concentration

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F2.1 | Sequenced HAI are NT-concentrated | NTWC+NTEC = 134/162 (82.7%) | A | Results — geography |
| F2.2 | Top hospitals dominate | Top 3 anonymized hospitals ≈ half of cases | A | Results / Supplement |
| F2.3 | HKWC & KCC absent from sequenced plot table | 0 sequenced HAI cases despite large bed bases | A | Results — sampling map |

---

## Finding set 3 — Density paradox (context, not cause)

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F3.1 | NT catchments are less population-dense than urban in-sample clusters | Mean pop density ratio NT/urban ≈ 0.27 | B | Results — context table |
| F3.2 | NT has lower hospital geographic density | Hospitals/100km² ratio ≈ 0.20 | B | Results — context table |
| F3.3 | NT has lower bed geographic density | Beds/100km² ratio ≈ 0.29 | B | Results — context table |
| F3.4 | Low density is not sufficient for high sequenced HAI | KEC ≈ NT hospital density but low HAI rate | A/B | Discussion — negative control |
| F3.5 | Density does not explain NT excess | Paradox is descriptive; no causal density model claimed | D | Discussion |

Sources: 2021 Census district areas/populations; HA cluster–district mapping; see `analysis/reference/SOURCES.md`.

---

## Finding set 4 — Admissions / throughput denominators

Annual HA cluster extracts (FY2021-22; DownloadCluster 36–42):

| Cluster | Sequenced HAI | Discharges | Patient-days | HAI/1000 beds | HAI/100k discharges | HAI/100k patient-days |
|---|---:|---:|---:|---:|---:|---:|
| NTWC | 78 | 140388 | 1331433 | 16.38 | 55.56 | 5.86 |
| NTEC | 56 | 163397 | 1311078 | 10.89 | 34.27 | 4.27 |
| KWC | 15 | 201096 | 1416390 | 3.03 | 7.46 | 1.06 |
| HKEC | 8 | 102190 | 822971 | 2.42 | 7.83 | 0.97 |
| KEC | 5 | 114909 | 809758 | 1.71 | 4.35 | 0.62 |
| KCC | 0 | 199221 | 1599173 | 0.00 | 0.00 | 0.00 |
| HKWC | 0 | 103370 | 635637 | 0.00 | 0.00 | 0.00 |

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F4.1 | NT excess survives discharge standardization | NT/urban mean ratio for HAI/100k discharges ≈ **6.86** | B | Results — denominators |
| F4.2 | NT excess survives patient-day standardization | NT/urban mean ratio for HAI/100k patient-days ≈ **5.74** | B | Results — denominators |
| F4.3 | NT excess survives bed standardization | NT/urban mean ratio for HAI/1000 beds ≈ **5.71** | B | Results — denominators |
| F4.4 | Annual denominators are imperfect for an 83-day window | Window-scaled discharge rates provided as sensitivity (`hai_per_100k_discharges_window_scaled`) | B | Methods / Supplement |
| F4.5 | NTWC has the highest mentally-ill patient-day volume | NTWC mentally-ill patient-days = 255,564 (share 19.2% of cluster patient-days) | B | Discussion — case-mix |

**Verdict (admissions):** Under three independent throughput denominators, NT sequenced HAI rates remain several-fold higher than urban in-sample clusters. Throughput alone does **not** wash out the geography.

---

## Finding set 5 — Sequencing effort / ascertainment

Tree window tips: hospital-coded = 148; community/other = 764 (BA.2.2/BA.2.12.1/BA.5.6 date files).

| ID | Claim | Evidence | Grade | Manuscript home |
|---|---|---|---|---|
| F5.1 | Hospital-coded tips track sequenced HAI by construction | Study design = targeted sequencing of suspected clusters; tip geography ≈ sample geography | A/D | Methods — limitation |
| F5.2 | Tip counts therefore cannot independently prove sequencing created NT excess | Circularity explicitly documented | D | Discussion |
| F5.3 | Community background sequences dominate the trees | 764/912 (84%) of window tips lack hospital codes | A | Methods — phylogeographic context |
| F5.4 | Official HA nosocomial bulletins also feature NT hospitals, especially Castle Peak | Mention-days: NTWC=12, NTEC=9 (`ha_nosocomial_mention_days.csv`) | C | Results — ascertainment |
| F5.5 | Ascertainment is incomplete / uneven: KCC appears in official bulletins but has 0 sequenced HAI in the plot table | KCC official mention-days = 10, sequenced HAI = 0; KCC beds = 6005 | C/A | Discussion — sampling bias |
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
