# Reference data provenance
# Dual-seeding / geographic concentration exploration

## district_demographics_2021.csv
- Population, land area, and density by District Council district from the
  2021 Population Census Summary Results (Census and Statistics Department),
  with land areas consistent with Lands Department geographic data
  (as tabulated for districts of Hong Kong).
- `ha_cluster_primary`: conventional HA planning mapping of districts to clusters
  (HA refined population-based model / LegCo answers). Islands is split in practice
  (HKEC vs KWC/Lantau) and is flagged `mixed_HKEC_KWC`.

## ha_cluster_capacity_2022.csv
- Hospital/institution, SOPC, GOPC counts and bed totals: HA Statistical Report
  2021-2022 (service capacity as at 31 March 2022).
  https://www3.ha.org.hk/data/HAStatistics/StatisticalReport/2021-2022
- Catchment population (2017 projection): LegCo LCQ13, 13 December 2017.
  https://www.info.gov.hk/gia/general/201712/13/P2017121300666.htm
- `land_area_km2_core`: sum of core district areas only (excludes ambiguous
  Islands/Lantau splits for HKEC/KWC). Intentional conservatism for density contrasts.
- Psychiatric bed fields filled where the all-cluster bed table reports them clearly
  (notably NTWC mentally-ill + mentally-handicapped).

## ha_cluster_throughput_2021_22.csv
- Inpatient discharges and patient-days by HA cluster from HA Statistical Report
  2021-2022 cluster extracts (DownloadCluster IDs 36–42 on
  https://www3.ha.org.hk/data/HAStatistics/), scraped 2026-07-11.
- Mentally-ill patient-day component included where extractable.

## ha_nosocomial_mention_days.csv
- Directional ascertainment proxy: number of distinct calendar days a hospital
  appears in the “clusters of nosocomial infection” table (or clear narrative
  ward-cluster text) of Hospital Authority press releases on info.gov.hk
  between 2022-05-28 and 2022-08-18.
- Incomplete by construction (HTML variation; not all bulletins use the same
  table). Use as directional evidence only (FINDINGS_LOG F5.4–F5.5).
- Raw table rows: `ha_nosocomial_press_table_rows.json`.

## Important limits
- Sequenced hospital-associated cases ≠ complete nosocomial incidence.
- Cluster catchment ≠ closed population (cross-cluster care is common).
- Land-area density for KEC/NTEC is diluted by large rural/country-park districts;
  that is a real geographic feature — interpret hospital-per-km² carefully.
- KEC is a useful negative control: low hospital density like NT, but low sequenced
  HAI burden in this table.
- Hospital-coded phylogenetic tips are largely the study’s own cluster sequences;
  they are not an independent sequencing-effort denominator.
- Annual FY2021-22 throughput is only an approximation for the 83-day study window.

