# prep_nhs: Nurses' Health Study breast tumour / adjacent (GSE93601, GSE115577), prep + power only

Status 2026-11-01 (update): PREP, GATES AND POWER DONE for GSE93601 (GPL22920, 445 pairs); no outcome statistic computed (observed D,
ER or proliferation contrast never evaluated). Not registered.
- `scripts/prep_nhs.py` (`python -I`): pairs (exactly 1 tumour + 1 adjacent per patient), RMA log2, probes -> symbol (mean), A97 G2M
  tertiles (149/148/148), ER from GSE93601 `estrogen receptor status` (= GSE115577 `er ihc`, 0 disagreements in 443; ER+ 357, ER- 86,
  2 missing). Reuses A177 `prep_a177.finish` by import. Outputs `nhs_gpl22920_{tumour,normal}.parquet` (~50 MB each, untracked),
  `nhs_gpl22920_pairs.csv`, `prep_log.csv`.
- `scripts/power_nhs.py` (3 s): `gates.csv`, `power.csv`, `power_arms.csv` (A177 `pair_gate`, A175 `Contrast`, MDE = q95 + 0.84 SD).
- GPL17585 (HTA 2.0, 178 pairs): OUT. No downloaded file holds their expression (series matrix table empty, the 1,110-sample csv holds
  only re-analyses of the GSE93601 CELs; values only in RAW.tar), see `hta_check_*` rows of `prep_log.csv`.
- Clock genes: all 14 of the project set (incl. RORA, CIART) present; G2M: 153 genes of the A97 set in the matrix (the 178 of the
  annotation include probes absent from the series matrix).
- Gates: reference gate on normals passes (p 0.032); tumour-normal pair gate FAILS (p_pairing 0.0005, p_acrophase 0.098): uninformative.
- Power (10,000 label permutations): H175.1 null SD 0.0799, MDE 0.198 > TCGA D 0.158 (not powered); ER contrast null SD 0.0833,
  MDE 0.206 (<= 0.3 stop rule; TCGA ER MDE was 0.408).

(Earlier status, kept for the record:)
Status 2026-11-01 (before expression files): PARTIAL. Metadata, deposit overlap, pair counts, ER / PAM50 / proliferation fields and platform clock-gene
coverage done (`scripts/meta_nhs.py`, `output/nhs_*.csv`). Expression files could not be downloaded (see below), so the gates and
the permutation nulls (A177/A175 code by import) are NOT run. No outcome statistic was computed.

Findings: the platform is Affymetrix (GPL22920 custom hGlue 3.0 HTA, and GPL17585 HTA 2.0), not Illumina DASL. GSE115577 holds
1,577 samples = 1,110 re-analyses of the GSE93601 CELs (`Reanalysis of: GSM...`, all 1,110 resolve) + 467 new HTA 2.0 samples.
No patient has samples on both platforms, so pairs are analysed per platform (445 + 178).
Blocker: `ftp.ncbi.nlm.nih.gov` gz bodies fail in the sandbox (curl 18, 0 bytes); data go to `$SCCC_DATA/nhs_gse93601`, `nhs_gse115577`.
Next steps once the matrices are present: adapt `A177/scripts/prep_a177.py` (log2 RMA values; probes are already one per gene on
GPL22920) and `power_a177.py` with stratifiers G2M tertile (A97) and ER (IHC). Use the deposited `proliferationscore` only as a
companion; the A97 G2M score is the registered one (178 of 200 G2M genes, clock genes removed, are in both annotations).
