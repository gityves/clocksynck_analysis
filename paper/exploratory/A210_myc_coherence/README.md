# A210: MYC amplification and tumour clock coherence at matched proliferation

Registration: `paper/docs/registered_experiments.md`, section A210. Script: `scripts/compute_a210.py`
(`python -I scripts/compute_a210.py --stage power|tests`; both stages take < 5 min with 4 workers, no `run.sh` needed).

- `--stage power`: extracts the MYC copy-number rows (METABRIC `data_cna.txt`, TCGA Xena GISTIC2), builds the cells, writes
  `output/{eligibility,counts,power,null_draws}.csv` (permutation null of T, MDE = 2.49 x null SD). No coherence of a real
  group is computed.
- `--stage tests`: observed T, one-sided permutation p, Holm over the scored hypotheses, 200 pseudo-clock sets, descriptive rows
  (`tests.csv`, `pseudo_T.csv`, `descriptive_cells.csv`, `descriptive.csv`). Run only after the power note is dated in the docs.

Choices not fixed by the registration: tertiles are taken over all tumours of the stratum with a copy-number call (any state); Holm
is over the hypotheses with MDE <= 0.25; TCGA T pools all eligible project x tertile cells (mean over cells); METABRIC pseudo sets are
the 200 of A166 (`A166_metabric/output/pseudo_sets_metabric.csv`); TCGA pseudo sets as A172 with seeds from crc32("A210"); TCGA primary
tumours are those of A172 `primaries()` (one per patient, lowest barcode; OV, BRCA etc. as the A97 projects) joined on the 15-character
barcode. `output/tcga_cache/` holds per-project clock counts, proliferation and MYC expression (open-access TCGA).
