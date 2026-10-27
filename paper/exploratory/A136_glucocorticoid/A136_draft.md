# A136 (DRAFT for PI review, 2026-10-27; not registered, nothing computed): is the loss of organ phase relations in brain-dead organ donors a drug effect?

Follows A132 and A132b (post hoc).
- Ventilated post-mortem cases keep the fast-death order of organ phase relations while being pulled toward 0 / 12 h
  (Fisher-Lee rho vs fast 0.678).
- Brain-dead organ donors show the pull without the order (rho 0.161; difference 0.517, permutation p 0.008, post hoc).
- Organ donors routinely receive high-dose corticosteroids during donor management. PER1 (and PER2) are direct
  glucocorticoid-receptor (GR) targets, so a GR pulse can move clock-gene expression organ by organ without brain death
  acting on the clocks.
- GTEx records no donor-management drugs. The only steroid field is long-term use (`MHSTRDLT`, protected). Exposure is
  therefore read from GR-target expression in the same tissues.

Data: GTEx v8, the A132 groups (Hardy 0 organ donors n 298; Hardy 0 ventilated post-mortem n 79; fast deaths as the
reference). Code: `paper/exploratory/A136_glucocorticoid/`. Controlled access: aggregates only, never per-donor rows.

## Exposure score (fixed before any look)
- GR-response genes: FKBP5, TSC22D3 (GILZ), ZBTB16, KLF15, DDIT4. Clock genes and PER1 are excluded. DUSP1 and SGK1 are
  excluded as immediate-early / ischaemia-responsive (A103).
- Per tissue: mean log CPM of the five genes, residualised on total ischaemic time within tissue, then z-scored within
  tissue over all GTEx donors. Donor score = median over the donor's tissues.
- Sensitivity: the same with `MHSTRDLT` donors removed.

## Hypotheses
- **H136.1 (exposure differs).** Organ donors have a higher GR score than ventilated post-mortem cases.
  - Statistic: median difference.
  - Null: donor-label permutation, 10,000 draws, one-sided.
  - This is a gate: if H136.1 fails, a drug cannot explain the A132b difference. H136.2 is then reported as descriptive.
- **H136.2 (dose-response within organ donors).** Fisher-Lee rho vs fast is lower in the top GR tertile of organ donors
  than in the bottom tertile.
  - Pairs: as A132b, >= 20 shared donors.
  - Null: tertile-label permutation within organ donors, 1,000 draws, one-sided.
  - Power (to be computed before scoring): about 99 donors per tertile; the A132b null SD at 79 vs 79 was 0.204, so the
    expected null SD is about 0.18 and the MDE about 0.45; minimum attainable p 0.001.
- **H136.3 (route through PER1/PER2).** The A132b organ-donor vs ventilated rho difference, recomputed with the 12-gene
  waveform without PER1/PER2 (as A103.3).
  - Read: "not through direct GR induction of PER1/PER2" if the difference stays >= 0.26 (half of 0.517) with
    p < 0.05; "through PER1/PER2" if it falls below 0.26 with p >= 0.05.
- **Specificity companions** (same permutations):
  - the same tertile split in ventilated post-mortem cases (expected: no dose-response);
  - a non-GR stress score (hypoxia: VEGFA, ADM, BNIP3, PGK1, SLC2A1) split into tertiles in organ donors (expected:
    weaker than GR).
- Holm over H136.1-H136.3.

## Readings
- **H136.1 and H136.2 pass:** the loss of order in organ donors follows GR exposure. The Figure 3 "loss of relations in
  brain death" becomes "loss of relations with donor-management glucocorticoids".
- **H136.1 passes and H136.2 fails:** exposure differs, but no dose-response at this power. The drug explanation is not
  supported, and the MDE is stated.
- **H136.1 fails:** drug exposure does not differ. The brain-death reading stands (still post hoc, A132b).

## Open choices for the PI
1. The GR gene set. Other candidates: PER1 (a clock gene, so it would make the score circular), MT2A, CDKN1C.
2. Whether H136.2 is the primary, with H136.1 as its gate.
3. An external replication (a brain-dead / multi-organ cohort; Q9 e) as the next step: none is identified yet.
