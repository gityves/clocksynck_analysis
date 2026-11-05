# hughey2016differential
**Ref:** Hughey & Butte, J Biol Rhythms 31:588 (2016). DOI 10.1177/0748730416668049. Type: data meta-analysis.
Lab: Atul Butte (UCSF) / Jake Hughey.
**Source:** full text

## Question and design
Are human brain and peripheral clocks in phase? Meta-analysis of 15 datasets (5 human, 7 mouse, 3 rat) + 13
published studies. Human blood: 3 timed in-lab studies (~8 samples per participant, ~20 participants each); human
brain: 2 post-mortem datasets (55 and 146 donors) with CT from time/date of death and location.
**Sampling:** blood timed around the clock; brain via time of death (rapid deaths only). Different people in each organ.

## Method
Core clock genes only (7 + DBP); ComBat by subject (blood) or area (brain); periodic smoothing spline (ZeitZeiger,
three knots) per gene and dataset; peak/trough times compared across datasets.

## Key results
- Headline: "after controlling for age, clock gene expression in brain was phase-delayed by ~8.5 h relative to that of blood".
- Quantified: "the absolute phase difference between blood and brain across all seven genes was 8.6 ± 2.1 h".
- Internal clock structure conserved: relative to ARNTL, "the differences between the two organs largely disappeared".
- Age: brain clock phase-advanced in donors >40 y (cf. Chen 2016).
- Human brain phasing "resembled that of the SCN"; other tissues ~12 h apart between diurnal and nocturnal species.
- Post-mortem progression cannot explain it: "even 8 h of postmortem progression of the clock could not explain the phase difference".

## Not addressed / limits
- Cohort-level comparison between different people (young volunteers' blood vs older donors' brain), different time references (lab clock time vs time of death converted to sunrise-based CT); not within person.
- Brain phase relies on recorded time of death being the circadian state at sampling.
- No disease contrast; no uncertainty on the 8.5 h at the person level.

## Link to us
- **Intro claim:** human organs may not be in phase (brain vs blood ~8.5 h), contrast with Talamanca 2023 (tissues within a few hours) — open question that within-person, time-free methods can address. Our refs table entry (brain vs blood ~8.5 h lag, pooled cohorts, not within donor) is correct.
- **Relation:** motivating evidence; within-donor CNS vs periphery offset is a direct test (GTEx brain vs other tissues; host anchored on neurons in CNS).

## Ideas
- GTEx: arg z of brain regions relative to the peripheral host. Is the offset ~8.5 h (Hughey) or a few hours (CHIRAL)? A discriminating result between two published claims, cheap to produce.
- Age split of brain offset (>40 vs ≤40) as a replication of their age advance, with aggregates only.
