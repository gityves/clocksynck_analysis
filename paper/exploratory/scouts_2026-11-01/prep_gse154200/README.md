# prep_gse154200: GSE154200 lesional / non-lesional skin, prep + power only

Status 2026-11-01: PARTIAL. Person counts from sample metadata (`scripts/meta_gse154200.py`, `output/gse154200_counts.csv`).
The count matrix (17.1 MB) could not be downloaded (`ftp.ncbi.nlm.nih.gov` gz bodies fail in the sandbox), so clock-gene
coverage, `Dataset.gate` and the label-permutation SD are NOT run. No contrast was computed.
Next: build the two-compartment (lesional, non-lesional) design as in A115/A148 (`clocksync.Dataset` with counts) per diagnosis
and pooled; null SD from permuted disease labels with labels only.
