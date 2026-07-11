# Decision: Table S3 upgrade status

## Context
Parent paper Table S3 summarizes official HA nosocomial infections (May 28–Aug 18 2022)
and states confirmed genomic clusters cover 51.9% of official nosocomial cases
(126/≈243). The curated Excel used in `1.metadata_visualization.R` is gitignored
and the Wiley supplement returns HTTP 403 in this environment.

## What we tried
1. HA Statistical Report cluster extracts → strong **admissions/patient-day** denominators (success).
2. Public info.gov.hk HA press releases → nosocomial **tables are incremental/daily**
   (often 0–3 cases per hospital per day), so summing maxima reconstructs only ~15
   cases — not Table S3.
3. Narrative paragraphs sometimes describe larger ward events but are not a complete
   cumulative registry.

## Decision for the first draft
- Do **not** pretend we recovered Table S3 case counts.
- Use a **presence/ascertainment contrast**: which HA clusters appear in official
  nosocomial bulletins vs which appear in the sequenced plot table (KCC present in
  bulletins, absent in sequences).
- Anchor overall coverage with the parent paper’s published 51.9% statistic.
- Keep admissions standardization as the primary quantitative denominator result.
- Flag recovering author-curated Table S3 as the single highest-value revision item.

This is Laidlaw-caliber honesty: ambitious analysis, transparent limits, still a
coherent manuscript argument.
