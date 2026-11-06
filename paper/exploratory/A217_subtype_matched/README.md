# A217 / A218: genomic-state and proliferation losses at matched molecular subtype (TCGA)

Registration: `paper/docs/registered_experiments.md`, sections A217 and A218 (power note after A218). Script: `scripts/compute_a217.py`
(`python -I scripts/compute_a217.py --stage power|tests [--jobs 4]`; power takes about 2 min). Template: A215 (same coherence, A210 cache
`A210_myc_coherence/output/tcga_cache`, MYC GISTIC from A210, FGA from A213, TP53 status from A215).

- `--stage power`: `subtypes_tcga.csv` (project, patient, subtype; from the stored cBioPortal API JSON in `$SCCC_DATA/tcga/pancan_subtypes/api/`),
  `purity_tcga.csv` (A127 purity proxy), `projects.csv` (coverage, levels, entry), `eligibility.csv`, `counts.csv`, `arm_composition.csv`,
  `power.csv`, `null_draws.csv`. Nulls come from permuted labels only; no coherence of a real group is evaluated as an observed statistic.
- `--stage tests` (not yet run): `tests.csv` (T, p, Holm, pseudo companion, original-design excess and ratio), `attenuation.csv` and
  `attenuation_jackknife.csv` (A), `pseudo_T.csv`, `descriptive_{cells,per_project,per_subtype}.csv`, `brca_descriptive.csv`.
  The A218 precondition and the jackknife need observed values and are evaluated there only. The stage was only run on random noise
  expression with the real labels and cells (`--smoke --out <scratch>`) to test the code path.

Choices left open by the registration:
- Subtype: cBioPortal `SUBTYPE` patient attribute (PanCanAtlas), as returned by the API 2026-11-06. Coverage = primaries with a label / primaries in the A210 cache.
  A project enters if coverage >= 50 % AND its labels have >= 2 levels. In BLCA, KIRC, KIRP, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, THCA the SUBTYPE equals
  the project name (no molecular subtype; matching is impossible, including them would only re-add the original design); they do not enter (OV also has 42 % coverage).
  Subtype levels are used as given (COAD / READ keep their `COAD_` / `READ_` prefixes; GBM `GBM`, IDHwt, IDHmut-non-codel; small levels drop out through the arm minimum).
- Tertiles: proliferation tertile within project x subtype over tumours with the subtype and the measure (rank "first", A210); FGA top / bottom tertile within cell.
  MYC arms GISTIC 2 vs 0; TP53 altered vs intact as A215 (intact requires a copy-number call).
- Purity proxy: mean of per-gene z (ddof 0) of log2 CPM of the 11 genes (all present, mean CPM >= 1, in all 23 projects) over the primaries of the project, then standardised;
  tertile within project over ALL primaries with a purity (before dropping tumours without subtype), rank "first". Tertile function: groups of < 3 tumours form one middle group.
- A218 cells: top / bottom proliferation tertile within each cell (middle tertile unused; null permutes top / bottom labels within cell); P_proj cells = project x purity tertile on
  the tumours (all three proliferation tertiles) of the eligible P_sub cells, proliferation tertile within the cell.
- Original design on the same tumours (descriptive, tests stage): projects with at least one eligible matched cell, tumours with a subtype and the measure, cells project x proliferation
  tertile, own permutation null (2,000) for the excess; excess ratio = matched excess / original excess.
- Jackknife: delete-one-project A_(-j); interval A +/- 1.96 x jackknife SE.
- Pseudo-clock sets (tests stage): 200 per project and amendment (seeds crc32("A217") / crc32("A218")), one family shared by the hypotheses of the amendment; P_proj and P_sub use the A218 family.
- Holm over H217.1-3 (scored ones); H218.1 single hypothesis. Null: 2,000 permutations, 20 chunks, seeds crc32("A217") / crc32("A218").
- BRCA: excluded from the tests; cells and nulls recomputed alone in the tests stage (`brca_descriptive.csv`).
