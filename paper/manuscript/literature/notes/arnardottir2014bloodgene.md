# arnardottir2014bloodgene
**Ref:** Arnardottir et al., Sleep 2014 (37:1589). DOI 10.5665/sleep.4064. Type: clinical / data (GEO GSE56931). Lab: Pack (UPenn) / Merck.
**Source:** full text

## Question and design
Human whole blood (PAXgene, custom Affymetrix array), n = 14 (5 men, 9 women, age 29.7 +/- 8.9 y), selected as 7 behaviourally resistant and 7 sensitive to sleep deprivation (PVT lapses in an earlier twin study). **Sampling:** timed, full cycle: blood every 4 h for 19 time points per person (08:00 day 1, through 38 h of wakefulness and recovery sleep to 08:00 day 4); 249 arrays passed quality control. Normal sleep-wake baseline = day 1 08:00 to day 2 20:00.

## Method (if any)
Mixed linear model with cos/sin terms, subject random effect, sex fixed (amplitude from beta1, beta2; phase from atan2), cross-checked with cosinor and JTK_cycle; linear time-awake term fitted separately. Needs a regular 24 h series per person.

## Key results
- "We found 4,481 unique genes with a significant 24-h diurnal rhythm during a normal sleep-wake cycle in blood" (FDR < 5%; 8,064 probe sets, 23.6% of genes); 6,465 probe sets (80%) found by all three methods.
- Most rhythmic probes peak 01:00-03:00 or 13:00-15:00. Clock genes NPAS2, ARNTL, NR1D2, DBP peak near these times; PER1-3 and CSNK1E peak 3-6 h earlier; CRY1/2, RORA, CLOCK not rhythmic at FDR 5%.
- Only two genes (SREBF1, CPT1A) show a linear trend with time awake; slopes -0.03 and 0.06 log2/h.
- "The main change with sleep deprivation was a reduction in the amplitude of the diurnal rhythm of expression of normally cycling probe sets." In 1,397 robust probe sets, resistant subjects lose amplitude (lights-on set: change -0.0143, paired t-test P = 1.7e-76; lights-off P = 1.7e-6) while sensitive subjects do not (P = 0.89 and 0.17).
- The resistant-subject effect is mainly carried by myeloid-lineage and cytotoxic T-cell modules (788 lights-on probes).
- Phase estimate differs from Moller-Levet 2013 by about 2 h; the authors attribute it to the time-awake regressor or to melatonin alignment.

## Not addressed / limits
- Whole blood only; no second tissue and no central-clock marker (no melatonin), so no internal desynchrony readout.
- n = 7 per group; a "resistant" amplitude loss is a statement on thousands of probes from 7 people (P values treat probes as independent), so P = 1.7e-76 is not a person-level test.
- Group baseline rhythms are equal; the difference only appears during sleep deprivation, where the amplitude drop partly reflects sleep-wake-driven rhythm loss.
- Cell composition (myeloid vs lymphocyte) is inferred from modules, not measured.

## Link to us
- **Intro claim(s):** refs row reduced circadian rhythmicity of transcripts in resistant individuals: supported, and the direction is the surprising one (resistant individuals dampen, sensitive do not), worth keeping in the sentence.
- **Relation to clocksync:** background for between-person variation in blood rhythm amplitude; data (GSE56931).

## Ideas
- Blood clock genes show different peaks (ARNTL about 14:00, PER about 8-11:00): the phase relations inside blood are the "host-referenced" structure clocksync would use; check A32/A35 blood clock readouts against these timings (grep shows no use of GSE56931 in the amendments).
- Person-level re-analysis of GSE56931 amplitude loss of the 13 clock genes (n = 7 vs 7) with a permutation over persons.
