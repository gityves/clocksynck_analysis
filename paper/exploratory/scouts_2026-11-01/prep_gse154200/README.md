# prep_gse154200: GSE154200 lesional / non-lesional skin, prep + gate + power (no contrast)

Status 2026-11-01: prep, gate and power done. No between-group or lesional-minus-non-lesional statistic was computed,
printed or saved (only permutation nulls). Not registered.

| file | content |
|---|---|
| `scripts/meta_gse154200.py` | person counts from sample metadata (`output/gse154200_counts.csv`, `gse154200_samples.csv`) |
| `scripts/prep_gse154200.py` | `--stage prep` (counts -> `clock_table.csv`, `clock_gene_coverage.csv`), `--stage gate` (`gates.csv`), `--stage power` (`power.csv`, ~5 s, no run.sh needed) |

Data: `$SCCC_DATA/skin_gse154200/GSE154200_countMat.txt.gz` (28,515 gene symbols x 564 libraries, parsed only). 376 libraries
have metadata (the other 188 are not in the three-diagnosis metadata and are ignored). Libraries are columns `MUC####`,
mapped to GSM by sample title.

## Design
- donor = patient, compartment = `lesional` | `nonlesional`, stratum = batch id (7 batches), one library per person x state,
  depth = genome-wide count, 14 clock genes (all present, median counts 84-8030, zeros <= 1.4 %: `clock_gene_coverage.csv`).
- One psoriasis patient has the two biopsies in different batches (clocksync needs one stratum per donor): dropped.
  Persons: psoriasis 87, nummular eczema 52, atopic eczema 46, pooled 185.
- `Dataset` anchors = non-lesional (the same person's unaffected skin as host reference), pooled and per diagnosis.

## Gate (`output/gates.csv`, `Dataset.gate`, 1,000 permutations)
All pass: psoriasis p_pairing 0.001, p_acrophase 0.001 (clock share 0.82); nummular 0.001 / 0.001 (0.65); atopic 0.002 / 0.001
(0.47); pooled 0.001 / 0.001 (0.76).

## Power (`output/power.csv`)
**The d* contrast is not estimable on a two-compartment design.** `power.mde` and `Dataset.contrast` need a host fitted on
at least 2 anchor compartments that exclude the target ("host needs >= 2 anchor compartments, got ['nonlesional']"), so
the registered MDE on d* (A153 `mde_table` logic) cannot be computed here, nor the 0.15 target checked. The estimable
quantities are the within-person paired ones of A54/A115 (`clocksync.paired`: offset of lesional against non-lesional, coherence,
r = log variance ratio, one projector for all persons). Power is therefore given on those readouts, from nulls only:
- between diagnoses: diagnosis labels permuted within batch, group sizes kept, 2,000 permutations; null SD of the
  difference of the group statistic (offset difference wrapped to +-12 h);
- lesional vs non-lesional within a diagnosis: per-person swap of the two labels (`paired.swap_*` logic), 2,000 swaps;
- MDE = 2.49 x null SD (one-sided 0.05, 80 %, the A115 rule). A153's s* damping route is d*-specific and not used.
No baseline per-compartment quantity was needed (the power step of A153 needs c0, which only the d* route uses).

| contrast | readout | MDE |
|---|---|---|
| psoriasis vs nummular (87 / 52) | offset, coherence, r | 0.70 h, 0.20, 0.38 |
| psoriasis vs atopic (87 / 46) | offset, coherence, r | 0.93 h, 0.20, 0.47 |
| atopic vs nummular (46 / 52) | offset, coherence, r | 1.42 h, 0.22, 0.42 |
| psoriasis vs other two (87 / 98) | offset, coherence, r | 0.69 h, 0.18, 0.35 |
| nummular vs other two (52 / 133) | offset, coherence, r | 0.85 h, 0.21, 0.39 |
| atopic vs other two (46 / 139) | offset, coherence, r | 1.01 h, 0.20, 0.42 |
| lesional vs non-lesional, psoriasis (87) | offset, r | 0.81 h, 0.40 |
| lesional vs non-lesional, nummular (52) | offset, r | 1.10 h, 0.49 |
| lesional vs non-lesional, atopic (46) | offset, r | 1.80 h, 0.51 |

For reference, the archive psoriasis effects were -1.0 to -1.7 h (A115 / A78 offsets); psoriasis lesional vs non-lesional
(MDE 0.81 h) is well powered on offset, with no d* equivalent. Specificity companion: the clock-gene design has no
non-clock compartment on the same permutations; the pseudo-clock and responsive-gene companions of `paired` (needs the
genome-wide mean log expression per series) are not yet run.
