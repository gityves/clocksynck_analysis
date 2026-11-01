# prep_nhs: Nurses' Health Study breast tumour / adjacent (GSE93601, GSE115577), prep + power only

Status 2026-11-01: PARTIAL. Metadata, deposit overlap, pair counts, ER / PAM50 / proliferation fields and platform clock-gene
coverage done (`scripts/meta_nhs.py`, `output/nhs_*.csv`). Expression files could not be downloaded (see below), so the gates and
the permutation nulls (A177/A175 code by import) are NOT run. No outcome statistic was computed.

Findings: the platform is Affymetrix (GPL22920 custom hGlue 3.0 HTA, and GPL17585 HTA 2.0), not Illumina DASL. GSE115577 holds
1,577 samples = 1,110 re-analyses of the GSE93601 CELs (`Reanalysis of: GSM...`, all 1,110 resolve) + 467 new HTA 2.0 samples.
No patient has samples on both platforms, so pairs are analysed per platform (445 + 178).
Blocker: `ftp.ncbi.nlm.nih.gov` gz bodies fail in the sandbox (curl 18, 0 bytes); data go to `$SCCC_DATA/nhs_gse93601`, `nhs_gse115577`.
Next steps once the matrices are present: adapt `A177/scripts/prep_a177.py` (log2 RMA values; probes are already one per gene on
GPL22920) and `power_a177.py` with stratifiers G2M tertile (A97) and ER (IHC). Use the deposited `proliferationscore` only as a
companion; the A97 G2M score is the registered one (178 of 200 G2M genes, clock genes removed, are in both annotations).
