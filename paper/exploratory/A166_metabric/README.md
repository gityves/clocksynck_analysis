# A166: METABRIC replication of the ER and proliferation clock-coherence losses

Registered: `paper/docs/registered_experiments.md` (A166). Data: `$SCCC_DATA/metabric/brca_metabric/` (cBioPortal, ODbL 1.0), TCGA BRCA from
recount3 + PanCanAtlas subtypes. Scout (metadata only): `../scouts_2026-10-30/metabric/SCOUT.md`.

Files (run with `python -I`, in order)
- `scripts/prep_a166.py`: parses METABRIC and TCGA BRCA into `output/{metabric,tcga}_meta.parquet` (tumour metadata, proliferation and purity scores)
  and `output/{metabric,tcga}_expr.parquet` (clock, G2M, purity and pseudo-pool genes, float32). Caches: do not commit.
- `scripts/a166_lib.py`: loaders, pseudo-set construction, group builders and the two statistics (`coh`, `wf_for`, `tertile` imported from A127's `a127_lib`).
- `scripts/run_a166.py [--nperm 2000 --npseudo 200 --jobs 4]`: power first (`power.csv` is written before any observed T), then tests, specificity,
  sensitivities, gene floor, TCGA secondary. Seeds: primary nulls `default_rng(166)`; pseudo sets 1660 (METABRIC), 1661 (TCGA); sensitivities
  `default_rng([166, job])`; TCGA nulls 1662.
- Outputs (`output/`, git-ignored): `power.csv`, `a166_tests.csv`, `cells.csv`, `strata.csv`, `sensitivity.csv`, `tcga_secondary.csv`, `tcga_cells.csv`,
  `gene_floor.csv`, null and pseudo draws, pseudo sets, `*_prep.csv`.

Choices not in the registration
- METABRIC expression: 218 duplicated symbols averaged over probes; the 16 NaN cells set to the gene mean; array floor reported as the registered 4.56 (observed minimum 4.558).
- Y = ln2 x log2 intensity (METABRIC) or ln2 x log2 CPM (TCGA), centred within each group, no Poisson term (as A127).
- Scores: gene z within COHORT, mean over genes (A97 `_zmean`), then z-scored within COHORT, over all 1,980 tumours with expression (not only the analysed ones). TCGA: within the 1,093 primaries.
  Tertile cuts are rank-based, so the scale does not matter.
- Tertiles: H166.1 purity tertile within COHORT (cohort-wide in the "without COHORT" sensitivity); H166.2 proliferation tertile within PAM50 (among tumours with ER and PAM50);
  H166.3 purity tertile within PAM50, proliferation tertile by rank within the stratum. Rank ties broken by order (`rank(method="first")`).
- H166.1 uses every tumour with ER_IHC recorded (including claudin-low and NC, which the PAM50 restriction does not apply to since PAM50 is not used); H166.3 uses every tumour with a PAM50 in the A127 five,
  ER not required.
- Pseudo sets: A127 construction (each clock gene replaced by one of its 20 nearest pool genes in standardised mean and SD); pool = all METABRIC genes except clock-like (`CLOCK_LIKE`), G2M, purity and clock genes
  (no detection filter: every array gene has a value); TCGA pool additionally requires counts > 0 in >= 90 % of tumours.
  The pseudo-set comparison uses the same strata and the observed labels.
- Verdict: "pass" = Holm p < 0.05 (one-sided, loss). Specificity label is given only for passes; for fails the column states "not applicable" and records whether T is nevertheless below the pseudo q05.
  `scanb_effect_detectable` = MDE <= |SCAN-B effect| (descriptive).
- Sensitivities: ER_STATUS applies to H166.1-2, claudin-low as its own stratum to H166.2-3 (PAM50 strata only), without COHORT to H166.1 (cohort-wide purity tertile), per COHORT to all three with tertiles recomputed
  within the cohort; every sensitivity has `nperm` permutations (p one-sided lower, not tested). Cohort-level groups with no cell that meets the arm/size rules are NaN with n_groups 0.
- Gene floor: clock genes only: mean, SD, minimum, 5th percentile, mean minus floor (also in SDs), fraction of tumours within 0.25 of the floor, mean by ER arm.
- TCGA secondary: H166.2 and H166.3 only (H166.1 is not in the registration's secondary); raw one-sided p (no Holm); MDE and 200-set specificity computed on the same rules. ER = XML IHC (Indeterminate and missing out);
  PAM50 = PanCanAtlas `SUBTYPE` restricted to the five; one primary per patient = lowest barcode (repo `load_tumours` convention).
