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

## Important limits
- Sequenced hospital-associated cases ≠ complete nosocomial incidence.
- Cluster catchment ≠ closed population (cross-cluster care is common).
- Land-area density for KEC/NTEC is diluted by large rural/country-park districts;
  that is a real geographic feature — interpret hospital-per-km² carefully.
- KEC is a useful negative control: low hospital density like NT, but low sequenced
  HAI burden in this table.
