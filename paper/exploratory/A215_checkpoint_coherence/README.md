# A215 / A216: checkpoint-deregulated proliferation, matched-state proliferation loss, and tumour clock coherence

Registration: `paper/docs/registered_experiments.md`, sections A215 and A216 (power note appended at the end of the file). Script:
`scripts/compute_a215.py` (`python -I scripts/compute_a215.py --stage power|tests [--jobs 4]`; the power stage takes about 5 min, run through `run.sh`).
Template: A213 (same coherence, A210 per-project cache `A210_myc_coherence/output/tcga_cache`, METABRIC loader and stored A166 pseudo sets; FGA read from `A213_heterogeneity_coherence/output/fga_{tcga,metabric}.csv`).

- `--stage power`: `status_{tcga,metabric}.csv` (per-sample TP53 / RB1 / CDKN2A status derived from MC3 and METABRIC mutation files plus GISTIC), `status_fractions.csv`,
  `eligibility.csv`, `counts.csv`, `power.csv`, `null_draws.csv`. Nulls come from permuted labels only; no coherence of a real group is evaluated as an observed statistic.
- `--stage tests` (not yet run): `tests.csv`, `pseudo_T.csv`, `attenuation.csv`, `descriptive*.csv`, `descriptive_tp53_class.csv`.
  The A216 precondition (P0 below its own pseudo 5th percentile) needs the observed P0 and is evaluated there only; a H216 hypothesis failing it is "not scored".

Data: TCGA MC3 public MAF (`/mnt/max8tb/.../tcga/mc3/`, PROVENANCE.md with URL, date, sha256); METABRIC `data_mutations.txt`, `case_lists/cases_sequenced.txt` from the stored tarball.

Choices left open by the registration:
- MC3 variants: FILTER kept if every flag is in {PASS, wga, native_wga_mix} (PanCanAtlas recommendation; 5,461 non-silent hits in the three genes before, 5,162 after). Sequenced = every 15-character barcode present in the MAF (10,295), irrespective of FILTER. Multiple aliquots of a sample are merged. METABRIC: no FILTER column; sequenced = `cases_sequenced.txt` (2,433 samples; 2,369 with a row in the MAF).
- Variant_Classification set: Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del/Ins, In_Frame_Del/Ins, Splice_Site, Translation_Start_Site, Nonstop_Mutation.
- Intact requires a GISTIC call at the gene(s) (TCGA Xena thresholded matrix, METABRIC `data_cna.txt`); a sequenced tumour with no mutation and no call is undefined and excluded. Altered needs a mutation or GISTIC -2 only. RB pathway: RB1 or CDKN2A.
- TCGA tumours: the 23 A210 projects, primaries; METABRIC: PAM50 strata over the A166 cohort.
- Proliferation tertile (H215.1-3, H215.4 pooled) within stratum over tumours with the status, rank method "first" (A210). Cells: all tumours of the cell, altered vs intact, >= 10 per arm.
- H215.4: tertiles within project x status (rank "first"); a project enters if all four arms have >= 10 tumours. Permutation stratum: status permuted within project x pooled proliferation tertile (tertile over all tumours with a status), the status-specific tertiles recomputed after each permutation (arm sizes are invariant).
- S215.5: within each eligible H215.1 cell, altered matched 1:1 to intact on the continuous proliferation score, caliper 0.1 x SD of the project's score over tumours with TP53 status, all pairs within the caliper taken in increasing distance (A213 `match`); a cell needs >= 10 pairs; swap null with probability 1/2. Its null is also simulated (power.csv) although only counts are registered.
- A216: FGA tertile within project / PAM50 over tumours with both TP53 status and FGA (rank "first"); strata = project x FGA tertile x TP53 status, proliferation tertile within stratum, >= 10 per arm. P0 = mean over projects (PAM50 classes) of coh(top) - coh(bottom) of the project proliferation tertiles over all tumours (all three tertiles) of the eligible P1 strata. Null of P1: top / bottom labels permuted within stratum (the middle tertile is unused, as A213).
- Null: 2,000 permutations in 20 chunks, seeds from crc32("A215") / crc32("A216"). Pseudo sets (tests stage): one family of 200 per TCGA project (all primaries) per amendment, METABRIC the 200 stored A166 sets.
- Attenuation bootstrap: tumours resampled with replacement within each eligible stratum (1,000), proliferation tertiles recomputed on the resample, A = mean P1 / mean P0, percentile 95% interval; reported with the fraction of draws with P0 >= 0 (A unstable then).
- "Residual proliferation imbalance within strata" (A216) is read as top - bottom mean of FGA, TP53 altered fraction and purity within the P1 strata (proliferation is the arm variable); `prolif_top_minus_bottom` also reported.
- Holm: H215.1-3 together, H216.1-2 together, over scored hypotheses; H215.4 uses its own one-sided p (secondary).
