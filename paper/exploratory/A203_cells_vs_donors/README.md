# A203: cells vs donors as the limit of internal-time precision (HLCA healthy lung)

Registration: `paper/docs/registered_experiments.md`, section A203 (commit 849d5b9). Scripts: `scripts/compute_a203.py`
(regrouping, acceptance, split-half fits), `scripts/score_a203.py` (V203.1-2, descriptive tables),
`scripts/rebuild_hlca_cells.sh` (the per-cell file with ARNTL, see below). Outputs (git-ignored): `output/`.

## Choices not in the registration
- **Regrouping rule** (to reproduce A17): cells of Census groups (donor x tissue x cell_type x assay x study x sampling
  method x ann_level_2 x ann_level_3) with >= 10 cells only, as the A17 pseudobulk build; lung sampling methods, not
  nose; the donor's majority assay (from the A17 stratum); classes from `common/datasets/_atlas/lung.py`.
- **Per-cell input rebuilt (2026-11-04).** The first acceptance run failed on ARNTL only: the stored per-cell file had
  ARNTL = 0 in every cell (Census symbol BMAL1; integrity note in results_preliminary). The builder was fixed and HLCA
  rebuilt (`rebuild_hlca_cells.sh`): cell ids and the 13 other genes equal the old file; acceptance then reproduces all
  1,321 A17 rows exactly (counts, cells, depth).
- **Half fits:** `Dataset(table, anchors)` with the A17 learned anchors present in the half, `fit(gate=False, n_boot=0)`,
  default LOO anchor readout; `offset_h` and `coupling` (truncated at the boundary, as the default readout).
- **Donor subsampling:** within strata, round(N / 198 x stratum size), at least one donor per stratum, so the realised N
  can differ from N by a few donors.
- **Thinning arm:** each row's counts binomially thinned by min(1, m / n_cells); depth scaled by the same factor.
- Seeds `SeedSequence([203, arm, level (0 = all cells), split])`.
