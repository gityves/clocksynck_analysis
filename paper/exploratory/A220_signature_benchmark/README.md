# A220: published-style circadian risk signatures vs random gene signatures

Registered in `paper/docs/registered_experiments.md` (A220). Status: code written, **power stage run; tests stage implemented, not run** (no
clock-list C-index or percentile has been computed).

## Files
- `scripts/a220_lib.py`: numba kernels (Breslow Cox likelihood, unpenalised Newton Cox with Wald p, L1 Cox by cyclic coordinate descent on a
  lambda path, CV-partial-likelihood lambda choice, Harrell C with ties, log-rank). `scripts/check_a220_lib.py` checks them against statsmodels
  `PHReg` (beta and SE agree to 1e-16, univariate p to 1e-9), the lasso KKT conditions, lambda_max, lambda -> 0 limit, a toy Harrell C (incl. ties).
- `scripts/compute_a220.py --stage power|tests [--nlists 200] [--nsplit 10] [--jobs 8] [--cohorts ..] [--smoke --out DIR]`.
  `--smoke` replaces expression by N(0,1) noise (real survival, covariates and splits), few lists and splits.
- `output/`: `eligibility.csv`, `power.csv`, `null_percentiles.csv`, `cache/rand_<cohort>.parquet` (every random list x split: C, incremental C,
  log-rank p, genes passing the filter, non-zero coefficients), `cache/cdr_os.csv` (public TCGA-CDR columns used).
  Tests stage writes `tests.csv`, `per_cohort.csv`, `selected_genes.csv`.

## What is computed (power stage)
Cohorts: 18 TCGA projects with >= 60 OS events among primaries with expression, METABRIC (1,756 with PAM50 and OS), SCAN-B GSE96058 (3,073).
Per cohort: GOBP_CIRCADIAN_RHYTHM (MSigDB 2024.1) restricted to expressed genes; 200 random lists of the same size, one gene per circadian gene
from the same median-expression decile; 10 stratified 70/30 splits (same splits for every list); in training univariate Cox p < 0.05 (else the 5
lowest p), L1 Cox with 5-fold CV on the cross-validated partial likelihood (lambda.min), risk score on the held-out 30 %; Harrell C and
incremental C (held-out C of Cox(covariates + training risk score) minus C of Cox(covariates), both fitted on the training split); mean over splits.
Percentile of a list = (#others below + 0.5 #ties)/#others. Null of P (H220.1: C; H220.2: incremental C; H220.3: incremental C, METABRIC and
SCAN-B): each random list ranked among the other 199 in every cohort, averaged over cohorts (cohort r uses the r-th list of each cohort; lists
are independent across cohorts). MDE = 2.49 x null SD. 200 random lists were run (no reduction); runtime 1,053 s with 8 workers (peak 4.4 GB).

## Choices left open by the registration
- TCGA projects: all recount3 projects except LAML; eligibility from CDR OS (OS.time > 0) among primaries with expression. One primary per patient
  as A210 (sorted barcode, first). READ and COAD are separate projects as in A210; READ, KIRP, PRAD, THCA etc. fall below 60 events.
- Proliferation metagene: A210's (mean z of the 200 Hallmark G2M genes on log2 CPM of genes with mean CPM >= 1, standardised, computed on all
  primaries of the project before restriction to OS). METABRIC / SCAN-B: the `prolif` column of the A166 / A127 tables.
- "Expressed": median log2 CPM above the 25th percentile of gene medians over the genes of the loaded matrix (TCGA: mean CPM >= 1 genes, as A210;
  medians over tumours with OS). Deciles for matching: deciles of the expressed genes' medians. The random pool excludes the whole GO list and the
  14 waveform genes. Microarray (METABRIC) and SCAN-B use log2 intensity / transformed FPKM in the same rule.
- Covariates TCGA: age, stage (I-IV from the CDR; IS, X, NOS = missing), proliferation, PanCanAtlas subtype (A217 table) if >= 50 % labelled.
  A numeric covariate observed in < 50 % of tumours is dropped (GBM, LGG, OV, MESO partly: no stage); partly missing ones are median-imputed
  with a missing indicator; unlabelled subtype = own level, levels < 10 pooled. Imputation medians use the whole cohort, coefficients the training split.
  METABRIC: age, PAM50 (+ claudin-low class not separate), proliferation, tumour size, log(1 + positive nodes). SCAN-B: age, PAM50,
  proliferation, tumour size, node status. Stage / size / nodes are "where recorded".
- Risk score entering the covariate model is the training-set linear predictor (in-sample); this makes the incremental C of a random list
  negative on average (mean about -0.02), a property of the published-style adjusted workflow.
- Standardisation and lambda path (25 values, ratio 0.01, glmnet-style objective -loglik/n + lambda|b|_1) use the training split; CV folds are
  stratified by event; ties by Breslow. If the lasso keeps no gene the score is constant (C = 0.5, incremental C = 0).
- "Validated" (descriptive): median over the 10 splits of the held-out median-split log-rank p < 0.05; also the per-list validation rate.
- Tests stage p: (1 + #{null P >= P observed})/(1 + 200) against the random-list null of P (clock percentile is among 200, null among 199 others);
  Holm over the hypotheses with `scored` = TRUE. Clock list also run for the 14 waveform genes (descriptive).
- Threads: `limit_threads(1)`, joblib 8 workers.

## Runtime estimate (tests stage)
Only the GO list and the 14-gene list run (2 lists x 10 splits per cohort) plus the per-split selection frequencies: about 2-4 minutes in total
(dominated by loading the 18 recount3 projects), run through `run.sh` anyway for the ledger.
