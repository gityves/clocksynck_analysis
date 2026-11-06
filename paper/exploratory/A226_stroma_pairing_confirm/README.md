# A226: tumour vs normal stroma clock pairing in independent cohorts (registered; not scored)

- `scripts/scout_*.py`, `scripts/fetch_census.py`: dated scout and downloads (`scout/SCOUT.md`, `scout_candidates.csv`;
  data in `$SCCC_DATA/a226`, 9.9 GB, PROVENANCE per source).
- `scripts/build_a226.py`: per patient x compartment pseudobulks (A149 conventions; marker-rule cell typing; GSE308103
  invasive LUAD only; GSE183219 Cell Ranger v2-style cell calling; Qian 2020 donors from Census atlas b617ee1b).
- `scripts/compute_a226.py --stage power|tests`: power run (simulated independent states at the real cohort sizes,
  calibrated on the A151 pseudo-set spread); tests implemented, not run (both compartments not scored: 11 / 13 patients).
- Note: builder logs share `paper/results/logs/build_a226.log`; GSE308103 build peaked at 29.4 GB (cap 30).
