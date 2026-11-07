# A228: does the tumour immune clock state carry body time; is the T/NK pairing composition? (registered; power stage only)

- `scripts/compute_a228.py --stage power`: patients / cohorts per test and the within-cohort permutation null of the uncorrected
  |complex correlation|, shuffled pairings only (identity pairing never evaluated, no observed statistic computed or stored), pooled with
  patient-count weights; MDE = 2.49 x null SD; scored if >= 40 patients and MDE <= 0.25. Also builds the T_NK subset pseudobulks and
  the subset-fraction descriptives.
  Outputs `output/{eligibility,power,subset_status,subset_fractions,subset_fraction_summary}.csv`.
- `--stage tests`: implemented, NOT run (observed statistic, permutation p, 200 pseudo-clock sets, Holm over scored Part 1 primaries and,
  separately, Part 2 primaries, descriptive per cancer type with >= 12 patients). Not-scored primaries are skipped unless `--force`.
- Cohorts, tables, statistic, null and pseudo sets: the 34 A227 cohorts (`compute_a227`, `compute_a226`, `a149_lib`).
- Part 1: tumour T_NK / myeloid vs adjacent-normal fibroblast (H228.1 / H228.2); descriptive: vs adjacent-normal endothelium, tumour
  endothelium vs normal fibroblast (control), adjacent-normal T_NK / myeloid vs adjacent-normal fibroblast.
- Part 2: per-cell marker rule inside T_NK (raw counts; CD8 T, CD8-negative T, NK; others dropped), genome-wide subset pseudobulks per
  patient x sample kind (>= 50 cells; sample median UMI >= 1,500 as the A227 units), cached in `$SCCC_DATA/a228/<cohort>.parquet`
  (+ per-sample `<cohort>_fractions.csv`). The A149 compartment-presence rule is not applied to subsets. Marker genes absent from a
  cohort count as zero and are listed in `output/subset_status.csv`; a cohort without any CD3, CD8 or CD4 gene is excluded there.
  Pseudo sets are rebuilt on the subset pseudobulks (reference log-mean from the normal-side subset).
- Run: `SCCC_DATA=/mnt/max8tb/data/CoPhaser_data/cancer_circadian RUN_CAP_GB=30 paper/common/run.sh
  paper/exploratory/A228_immune_body_time/scripts/compute_a228.py --stage power|tests --jobs 2`.
- Patients are not de-duplicated across cohorts (as A227).
