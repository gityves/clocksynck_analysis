# A156 stage 1 inventory (blind: metadata, annotation and per-probe missingness only; no clock statistic)

Sources (PROVENANCE.md in `$SCCC_DATA/geo/GSE48113/`): series matrix (72.9 MB), GPL15331 platform table (25.4 MB), parsed
only (`scripts/parse_inputs.py`, `python -I`). Outputs: `output/samples.csv`, `probe_map.csv.gz`, `probe_missing.csv.gz`.

## Samples and subjects
- 287 arrays (not 22 x 14 = 308), 22 subjects (BB0012 ... BB0289). Day 1: 149 arrays; day 4: 138.
- Platform GPL15331: Agilent-026817 University_Surrey_HumanSleep_44k_v1 (Agilent 4x44k, one colour). Matrix values are
  authors' QC-filtered, quantile-normalised (limma) log2-scale intensities, with missing values (empty) for filtered probes.
- **Subjects with both days: 21** (BB0104 has day 1 only, 6 arrays, no day-4 array).
- Samples per subject-day (both-day subjects): 15 subjects have 7 + 7; the others (day 1 / day 4): 5/7, 6/5, 6/7, 7/4 (x2),
  7/6. Subjects with at least m arrays in both days: m = 7: 15; m = 6: 17; m = 5: 19; m = 4: 21.
- Per time point (1-7) arrays: day 1: 20, 22, 22, 21, 20, 22, 22; day 4: 21, 21, 21, 21, 17, 19, 18.
- Other covariates: PER3 genotype (4/4 144 arrays, 5/5 143; constant within subject).

## Time fields (as recorded)
- `time sample taken` (clock time, hh:mm), `day sample was taken in protocol` (1 or 4), `time point` (1-7), `sleep timing`
  ("In phase" on day 1, "Out of phase" on day 4, with respect to melatonin). No melatonin phase value and no hours since
  wake are given.
- Sampling is every 4 h (steps 3.9-4.2 h within subject-day): time point 1 at 18:00-22:15, point 7 at 18:00-22:15 (the 24-h span
  plus the first point of the next window); the clock time of point 1 varies by subject (individual wake time).
  Day 1 and day 4 of a subject start at the same clock time (e.g. BB0121: 20:25 on both).
- One missing clock time: BB0264 day 1 point 7 ("NA").
- Metadata inconsistencies: BB0121 day 1 point 7 is labelled "Out of phase" (it shares clock time 20:25 with his day-4
  point 1, so it may be the boundary sample); no duplicate (subject, day, point) rows.
- A cosinor against sample time would use these clock times (24 h period); the melatonin phase needed to express time
  relative to the central clock is not in the matrix.

## Probe mapping (GPL15331 `GENE_SYMBOL`)
- 43,758 probes in the matrix, all present in the platform table; 317 probes have no value in any array. Platform control
  types: FALSE (genes) 42,119, pos 1,486, ignore 1,462, neg 153; 30,505 probes carry a gene symbol.
- Proposal (as A124): non-control probes (CONTROL_TYPE FALSE) whose GENE_SYMBOL, mapped with `canonical`, equals the gene;
  mean of log2 over the gene's probes (missing values skipped).
- Probes per core gene under canonical: ARNTL 1, CLOCK 1, NPAS2 2, PER1 1, PER2 2, PER3 3, CRY1 2, CRY2 3, NR1D1 2, NR1D2 2,
  DBP 2, RORA 3, TEF 3, CIART 2.

## Gene coverage (probes present in the matrix)
| set | `canonical` only | + 8 A124 aliases (not used yet) |
|---|---|---|
| core (14) | 14 | 14 |
| cue (40) | 34 (missing TENT5B, MTERF2, H2BC21, TMEM263, PMS2CL, SEPTIN4) | 39 (missing PMS2CL) |
| locked (40) | 35 (missing ADIRF-AS1, AGAP11, CENPC, CLUHP3, RETREG1) | 38 (missing AGAP11, CLUHP3) |

## Decisions needed before readout
1. Minimum arrays per subject-day for the per-subject offset (proposal 5 in both days: 19 subjects; 7+7 only gives 15, all 21
   gives n = 21 with 4-array offsets; the exact swap is 2^n: 2^19 = 524,288 enumerable, or 10,000 random per the registration).
2. Missing values in the matrix: probes with missing values in some arrays are averaged over present values per gene;
   gene-level missing (all probes missing in an array) must be imputed or the array dropped (proposal: gene mean of the
   subject-day, reported).
3. Timed cosinor: clock time (24 h), as recorded; BB0264 point 7 has no time and is dropped from the cosinor only.
4. Aliases: canonical only (primary) with the alias table as a labelled post hoc sensitivity, as in A124.
