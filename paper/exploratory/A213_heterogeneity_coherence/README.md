# A213: heterogeneity / copy-number burden and tumour clock coherence at matched proliferation

Registration: `paper/docs/registered_experiments.md`, section A213 (power note appended at the end of the file). Script:
`scripts/compute_a213.py` (`python -I scripts/compute_a213.py --stage power|tests [--jobs 4]`; minutes, no `run.sh` needed).
Template: A210. Reuses the A210 per-project cache (`A210_myc_coherence/output/tcga_cache`, `myc_cna_tcga.csv`) and its stored METABRIC pseudo sets (A166).

- `--stage power`: FGA (cached `output/fga_{tcga,metabric}.csv`), cells, `eligibility.csv`, `counts.csv`, `power.csv`, `null_draws.csv`. No real-group coherence.
- `--stage tests`: `tests.csv`, `pseudo_T.csv`, `s213_4.csv`, `s213_4_cells.csv`, `descriptive*.csv`.

Choices not fixed by the registration:
- FGA = fraction of genes with GISTIC != 0 over genes with no missing value in any sample of the matrix (TCGA Xena thresholded GISTIC2:
  24,776 genes; METABRIC `data_cna.txt`: 22,079 of 22,544). Samples joined on the 15-character barcode (TCGA) or patient id (METABRIC).
- ABSOLUTE rows with `call status` == called, first row per `array` (15-character barcode).
- Tertiles use `rank(method="first")` as A210 (`a127_lib.tertile`): ties (many tumours have subclonal fraction 0) are split by table order
  (barcode order), arbitrary but deterministic. Proliferation tertile within stratum over the tumours that have the measure; measure tertile
  within the cell over the same tumours; a cell needs >= 10 per arm. The two outer arms only are permuted (middle tertile unused).
- TCGA: all 23 A210 projects (A97 + A172 lists), primaries as A210. METABRIC: PAM50 strata with an FGA value. T pools all eligible cells (mean over cells).
- Holm over the scored hypotheses (MDE <= 0.25). Null: 2,000 within-cell permutations in 20 fixed chunks, seeds from crc32("A213").
- Pseudo sets: TCGA 200 per project from all that project's primaries (so one family serves H213.1 and H213.2; seeds crc32("A213")), A210 construction; METABRIC the 200 stored A166 sets.
- S213.4: A210 H210.2 cells (project x proliferation tertile over tumours with a MYC call; >= 10 amplified and >= 10 neutral), tumours without FGA
  dropped, then 1:1 matching on FGA without replacement with caliper 0.05, all candidate pairs taken in increasing distance (ties by position);
  a cell needs >= 10 pairs to enter T. Null: labels swapped within pairs with prob. 1/2 (2,000). "FGA carries the result" =
  (T_matched - null mean) > 0.5 x A210 excess (-0.0848) and p >= 0.05.
- Descriptive A176 components: per cell, groups centred, `components(y, 0, 1, wf)`; dlog_plane, dlog_off, ds = dlog_plane - dlog_off, mean over cells.
  Residual purity imbalance uses ABSOLUTE purity (TCGA only).
