# A222: clock-gene alterations against background (registered in `paper/docs/registered_experiments.md`)

`scripts/compute_a222.py --stage power|tests` (about 1.5 min; gene table cached in `output/gene_cache.pkl`).

Inputs: MC3 v0.2.8 public MAF (FILTER in {PASS, wga, native_wga_mix}, as A215; coding classes only), Xena GISTIC2
thresholded by gene, GENCODE v29 basic (`$SCCC_DATA/gencode_v29`, downloaded for this amendment), TCGA recount3.

Choices the registration left open:
- Genes: GENCODE v29 protein-coding, chrY / chrM dropped; CDS length = longest basic-transcript CDS. MC3 / GISTIC symbols
  matched by name; CIART is C1orf51 in both files. Genes not found by name count as 0 mutations / NaN copy number.
- Expression: median log2 CPM (protein-coding library) per project over the first 30 primary tumours by barcode, then
  the median over the 33 projects.
- Hypermutators: > 1,000 coding calls (silent included) after FILTER.
- Matched mutation null: each clock gene replaced by a non-clock gene of its CDS-length decile x expression decile cell.
- Copy number: project from the recount3 patient barcode; projects with >= 100 GISTIC samples (26). Flank = 25
  protein-coding genes each side in GENCODE start order, any arm; genes without a full flank excluded from the null.
- Copy-number null: 1,000 random sets of 14 non-clock protein-coding genes with a full flank (not matched further).
