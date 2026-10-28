# A124 stage 1 inventory (blind: metadata and annotation only; no expression value was read)

Sources (PROVENANCE.md in `$SCCC_DATA/geo/GSE106800/`): series matrix (13.8 MB, sha256 14595aa3...), GPL17692 platform
table (458 MB, text, parsed only; the 2.4 GB family soft was not downloaded). Parser: `scripts/parse_inputs.py`
(`python -I`); outputs `output/samples.csv`, `output/probe_map.csv.gz`.

## Arrays and subjects
- Series matrix: 47 arrays (GSM2850460-), **12 subjects (RAS1, 2, 6-10, 13-17), not 14**; 53,617 probe sets (transcript
  clusters, RMA log2, Oligo 1.42.0, already normalised by the authors; no per-sample QC flag in the matrix).
- **11 subjects have all 4 biopsies**; RAS9 has 3 (misaligned 7 PM missing). Subjects missing entirely: RAS3-5, RAS11-12
  (the registration's n = 14 is not available; n = 11 complete, 12 with RAS9).
- Arrays per condition x time (as recorded): aligned 7 AM 12, aligned 7 PM 12, misaligned 7 AM 12, misaligned 7 PM 11.
- Times as recorded in the matrix: "7 AM (7h)" and "7 PM (19h)" (Methods say 8:00/20:00; the GEO record says 7/19).
  Titles label the 7 PM sample "evening (baseline)" and the 7 AM sample "morning". No other clock time is recorded.
- Crossover order from the CEL file names (P1/P2): aligned first in 7 subjects, misaligned first in 5.
- Other fields: gender (all male), age, height, weight, BMI, fasting glucose, insulin, FFA, TG (per array). Condition labels:
  "circadian aligned" (control) and "circadian misaligned".

## Probe-to-gene mapping (GPL17692 `gene_assignment`)
- 53,617 probe sets: 44,629 category main, 3,575 normgene->intron, 3,515 rescue, 1,626 normgene->exon, 273 control/reporter.
  All matrix IDs are in the platform table.
- Symbols are from the 2013-2020 platform annotation (pre-2019 HGNC names; e.g. C1orf51, FAM46B, SEPT4).
- Proposed rule (for the lead's data-rule note): keep category `main` probe sets; a gene's probe sets are those whose
  `gene_assignment` lists the symbol; if some of them list the symbol alone, use those only, else all (this matters for
  NR1D1 only: both of its probe sets are THRA;NR1D1 read-through, none is NR1D1 alone); collapse to genes by the mean of
  log2 values. Probe sets with several symbols are not otherwise dropped.
- Symbol aliases: exact match to the set symbols misses genes that are present under their old names. Option (a) exact
  symbols only; option (b) a fixed alias table (below) applied before any readout. Each alias maps to exactly one
  main probe set.

## Genes present
- Core (14): **13 exact**: ARNTL, CLOCK, NPAS2 (2 probe sets), PER1, PER2, PER3, CRY1, CRY2 (2), NR1D1 (2, read-through),
  NR1D2, DBP, RORA, TEF. Missing exact: CIART (platform name C1orf51, 1 probe set; alias option gives 14).
- Cue set (A81b, 40): **35 exact**. Missing exact: TENT5B (FAM46B), MTERF2 (MTERFD3), H2BC21 (HIST2H2BE),
  TMEM263 (C12orf23), SEPTIN4 (SEPT4); with aliases 40.
- Locked set (A81b, 40): **36 exact**. Missing exact: ADIRF-AS1 (C10orf116), CENPC (CENPC1), RETREG1 (FAM134B),
  CLUHP3 (no match under any name tried); with aliases 39.
- Probe-set counts per gene are in `output/probe_map.csv.gz` (columns probe, category, symbols_all).

## Design notes for the registered permutation
- n = 11 complete subjects: 2^11 = 2,048 within-subject label swaps (exact). With RAS9 included as a 3-array subject:
  2^12 = 4,096, but the swap of his 3 arrays is not defined for the missing cell; proposal is the 11 complete subjects for
  the primary and RAS9 only in descriptives (lead to decide).
- With 2 times per condition and states centred within subject over all 4 arrays, the within-subject label swap permutes
  the condition labels of the (7 AM, 7 PM) pairs.
- Smallest attainable exact p (n = 11): 1/2048 = 4.9e-4 (two-sided in the swap-symmetric sense, see stage 2).

## Decisions needed from the lead before readout
1. n: 11 complete subjects (proposal) vs 12.
2. Aliases: exact symbols (13 / 35 / 36 genes) vs alias table (14 / 40 / 39).
3. NR1D1: keep the two read-through probe sets (proposal) or drop.
4. Time label used only descriptively: 7 AM / 7 PM (GEO record).
