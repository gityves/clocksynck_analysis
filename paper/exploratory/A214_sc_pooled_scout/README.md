# A214 scout: pooled single-cell round, step 1 (inventory and power; no outcome read)

Prompt: `paper/docs/prompts/sc_pooled_round.md`. Nothing here is registered. No malignant-vs-normal D difference, no
pseudobulk amplitude and no clone phase was computed.

## Inventory
- On disk: `inventory/INVENTORY_DISK.md`, `inventory_studies.csv`, `inventory_patients.csv` (`scripts/inventory_disk.py`).
  Eligible patients (>= 100 malignant and >= 100 normal epithelial cells): (a) adjacent-normal epithelium 108 (24
  studies), (b) in-tumour non-malignant epithelium 115, either 168. A149's 17 extra cohorts have pseudobulks only (no
  per-cell table): 48 more under (a) if their raw files were re-extracted.
- Clock UMI per cell over the 168 eligible patients (median of per-patient medians): malignant 2.0, adjacent normal 1.0,
  in-tumour non-malignant 0.5. The normal comparator fails the METHOD.md rule (median >= 2) in most patients: with the
  rule applied to both groups, 38 patients remain.
- Public sources (`inventory/INVENTORY_PUBLIC.md`): the 3CA portal was unreachable (connection reset on every try), so
  its v2 study table is unchecked; the old Zenodo export (5834545, 10 of 89 files read) has no patient, malignant or gene
  labels: unusable. CELLxGENE Census 2025-11-08: one new lineage-matched droplet study (Moorman 2025 colorectal
  metastasis, 14 eligible donors, about 14.8k UMI per cell, 0.75 GB) and LuCA (41 donors, mixed assays, partly on disk,
  17.6 GB). Numbers are the agent's metadata query, not re-read.
- No clone labels on disk: the copy-number calls keep a per-cell score and clonality, not per-cell profiles.

## Effect sizes (exact von Mises model, `output/power/model_*.csv`)
- Rewiring (independent N(0, sigma^2) acrophase shifts per gene and tumour): D falls to 0.92 / 0.74 / 0.51 / 0.30 of the
  intact value at sigma 1 / 2 / 3 / 4 h, 0.04 for uniform phases; the same over lam 1.5-2.5 and R 0.25-0.7. **The 2 h
  equivalence bound is a D difference of -0.26 D_intact**, about -1.0 at the normal-epithelium level D 3.90 (WP0 item 2).
- Damping vs dephasing at a matched pseudobulk amplitude ratio a: D ratio damping 0.26-0.27 / 0.50-0.53 / 0.73-0.75
  (a 0.5 / 0.7 / 0.85), dephasing 0.78-1.09 / 0.89-1.10 / 0.95-1.07. Gap between the two 0.52-0.83 / 0.39-0.60 /
  0.22-0.34 D_intact. Separation also needs a: the ratio of pseudobulk amplitudes is not identified without the malignant
  host coupling (H138.1 not estimable) or an assumption that malignant and normal pseudobulk phases have the same
  between-patient spread (Assumption W applied to both); `paired.log_var_ratio` also carries patient-specific variance.

## Measurement noise and power (`output/power/{noise_groups,power}.csv`, `scripts/power_scout.py`)
Split-half D per group (patient x compartment): half difference only, mean 0 by construction. MDE = 2.49 x SE of the
pooled mean of the per-patient differences; between-patient heterogeneity excluded, so MDEs are lower bounds.
- Per-group SE (rms / robust): malignant 17.8 / 7.6, adjacent normal 25.2 / 12.8, in-tumour non-malignant 23.7 / 13.6.
  Strongly depth-limited: robust SE falls about 4-fold from the lowest to the highest tertile of UMI x sqrt(cells).
- MDE of the paired D difference: (a) 108 patients 5.1-5.5; (b) 115, 7.3-7.9; either 168, 5.2-5.4; with the UMI rule
  (38 patients) 4.7-5.9. **Against the 2 h bound (about 1.0) every design is 5-fold short**; N for MDE 1.0 at the
  observed per-patient SD (14.7-21.5): about 1,300-2,900 patients. Even a fully scrambled malignant clock (D difference
  about -3.8) is below the MDE of the paired contrast.
- Damping vs dephasing (gap 0.9-3.2 D units at a 0.5-0.85): below the MDE at every a; and a itself not identified.
- The one-sample malignant D (228 groups) has SE 1.18, so D above 0 is detectable (the WP0 "scrambling rejected" rests on
  the malignant compartment alone); the comparator is what fails.

## Reading
Q1 (rewired vs intact, 2 h bound), Q2 (damping vs dephasing) and Q3 (the same by proliferation, which halves n) are not
estimable with the cohorts on disk; the shortfall is depth of the normal epithelium (about 1 clock UMI per cell), not the
number of studies. Q4 (clones) needs clone calls not on disk and a pseudobulk clock signal in malignant cells that A187 /
A195 did not find. See the proposal in the lead's report.

## Depth needed (`output/power/depth_scaling.csv`)
Median regression over the 543 groups: per-group D SE ~ (clock UMI per cell)^-0.61 (SE 0.06) x cells^-0.41 (SE 0.06).
At 500 cells per group and the same depth in both compartments, patients for MDE 1.0 (about the 2 h bound): 1,264 at 2
clock UMI per cell, 543 at 4, 233 at 8, 100 at 16 (x 4 if the tails stay as heavy as now: rms SE about 2x the robust
one). Clock UMI is about 3.5-4 per 10,000 total UMI, so 8-16 clock UMI means roughly 20,000-40,000 UMI per cell. The fit
rests mostly on shallow groups (90 % below 6.8 clock UMI per cell, 99 % below 15, maximum 26): 16 and above are
near the edge of the data. Scout of deep studies: `inventory/INVENTORY_DEEP.md`.
