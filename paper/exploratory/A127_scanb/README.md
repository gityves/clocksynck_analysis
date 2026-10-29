# A127: SCAN-B (GSE96058) clock coherence vs proliferation and ER

Registered: `paper/docs/registered_experiments.md` (A127). Data: `$SCCC_DATA/scanb_gse96058/` (PROVENANCE.md).
Scripts (run with `python -I`, in order): `prep_scanb.py` (matrix read once, genes kept as float32, sample table), `feasibility.py`,
`run_a127.py` (all tests; 2,000 permutations, 200 pseudo sets; about 4 min). `a127_lib.py`: loader and A97 definitions.
`output/scanb_genes_kept.parquet` (344 MB) is a cache: do not commit.

Choices not in the registration
- Scale: published log2(FPKM + 0.1) (floor -3.32), not counts. y = ln2 x log2 FPKM, centred in each group; the A97 `components` call is used
  unchanged with zero Poisson variance (coherence does not use it). Scores z-scored within the 3,073 analysed tumours.
- Purity tertile for the H127.1 strata: within PAM50 subtype (sensitivity: cohort-wide). H127.2 uses cohort-wide purity tertiles.
- Pseudo sets: variance-matched (own draw: 20 nearest pool genes in standardised mean and SD, pool = detected in >= 90 % of samples, clock-like, G2M and
  purity genes excluded), because `pseudo_waveforms` matches mean only; that API version is reported as a sensitivity.
- Off-plane survival score: no normals in SCAN-B, so the deviation is from the cohort mean, without Poisson correction.
