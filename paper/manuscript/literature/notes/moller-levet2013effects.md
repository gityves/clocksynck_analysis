# moller-levet2013effects
**Ref:** Moller-Levet et al., PNAS 2013 (110:E1132). DOI 10.1073/pnas.1217154110. Type: clinical / data (GEO GSE39445). Lab: Dijk (Surrey).
**Source:** full text

## Question and design
Human whole-blood transcriptome (microarray) after 1 week of sleep restriction vs 1 week of sufficient sleep, n = 26 healthy adults, balanced cross-over. **Sampling:** timed, full cycle: 10 blood samples per person and condition, every 3 h, during 39-41 h of constant-routine sleep deprivation after each week. Clock hours are expressed on a melatonin-aligned scale (individual melatonin midpoint 0415 control vs 0501 after restriction), not on wall time.

## Method (if any)
Not a method paper. Mixed ANOVA for condition x circadian time bin; per-person cosinor-type time-course fits; a gene counts as "prevalent circadian" if rhythmic in enough participants to give FDR < 5%. Assumes a full 24 h+ series per person and individual melatonin alignment.

## Key results
- Sleep restriction was 5.70 h (SEM 0.03) vs 8.50 h (SEM 0.11) of sleep per 24 h. "reduced the number of genes with a circadian expression profile from 1,855 to 1,481" (abstract); 793 genes rhythmic in both, 688 only after restriction.
- 711 genes (3.1%) changed level (444 down, 267 up), including PER2, PER3, RORA, DEC1, CSNK1E. "reduced the circadian amplitude of these genes" (abstract); amplitude reduction P < 2.2e-16 (Fig. 5D), also in genes rhythmic in both conditions.
- Phases are bimodal (night, in phase with melatonin: 1,072 genes, 57.8%; day: 591 genes, 31.9%). After restriction the loss falls mainly on day-peaking genes; "the average expression peak time of cluster 2 (processes associated with regulation of gene expression) had moved later to 0250 hours" and immune clusters moved earlier (1600 to 1430 h; 1906 to 1710 h), a narrowing of the window of rhythmic expression.
- Genes responding to time awake rose from 122 to 856 after restriction (P < 2.2e-16).
- Melatonin phase delay of about 46 min after restriction shows the central pacemaker also moved a little; transcript phases were therefore aligned per person to melatonin.

## Not addressed / limits
- One tissue (whole blood, a mixture of leukocytes); no second compartment, so no between-tissue (internal) desynchrony readout. Change is dampening and a within-blood phase redistribution.
- Rhythmic-gene counts are thresholded, so the drop from 1,855 to 1,481 mixes amplitude loss and power loss; amplitude is compared directly but only for called genes.
- Rhythms are read under sleep deprivation, which itself abolishes most sleep-wake-driven rhythms, and are aligned to melatonin; the result is conditional on a constant routine. Blood composition is not modelled.
- Arnardottir 2014 finds peaks about 2 h earlier on the same kind of data and attributes it to the time-awake term and to melatonin alignment: absolute phases are method-dependent.

## Link to us
- **Intro claim(s):** refs_circadian_desynchrony.md: loss of rhythmicity within one compartment, not between organs for this key with archer2014mistimed. Supported (indeed one compartment; the 1,855 to 1,481 and amplitude numbers are valid). The refs row calls it one week of sleep restriction, correct.
- **Relation to clocksync:** background; data (GSE39445; cited as the array-input format reference in registered_experiments.md line 391) but blood only.

## Ideas
- Not proposed as new: blood-only forced-desynchrony style test is done (A156, GSE48113, Archer 2014: untimed core-cue change not detected, MDE 13.1 h). GSE39445 is a within-person sleep restriction series, 10 samples x 2 conditions: a check of the A35 "blood clock weak" finding (rho 0.31) on a full 24 h+ series, not yet run (no mention of GSE39445 in results_preliminary.md).
- Compare per-person amplitude of the 13 clock genes, control vs restriction, as a within-person dampening positive for blood.
