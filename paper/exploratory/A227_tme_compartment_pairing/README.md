# A227: tumour vs normal pairing of TME compartments (registered; power stage only)

- `scripts/compute_a227.py --stage power`: patients / cohorts per compartment and the within-cohort permutation null of the
  uncorrected |complex correlation| on the real states, shuffled pairings only (identity pairing never evaluated, no observed
  value computed or stored), pooled with patient-count weights; MDE = 2.49 x null SD; scored if >= 40 patients and MDE <= 0.25.
  Outputs `output/{eligibility,eligibility_by_cancer,power}.csv`.
- `--stage tests`: implemented, NOT run (observed statistic, permutation p, 200 pseudo-clock sets, Holm per tier, descriptive).
- Cohorts: the 27 A149 cohorts with estimable S1 + the 7 A226 cohorts. Existing-cohort stored `cancer_sc` tables already carry
  myeloid / T_NK / B_plasma / mural / mast, so no rebuild was needed (A149 new cohorts: cached `A149/output/pb`; A226: `A226/output/pb`).
- Eligibility: >= 4 patients with tumour and normal pseudobulks (A149 table rules already applied); null = 10,000 permutations
  per cohort, independent streams `[227, kind index, cohort hash]`.
- Anchors (fibroblast, endothelial) are descriptive; power is reported for them but they are not hypotheses.
- Note: patients are not de-duplicated across cohorts (A149 / A226 overlap, if any, is as in the registered cohort sets).
