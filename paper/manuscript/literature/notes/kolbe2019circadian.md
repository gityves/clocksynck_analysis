# kolbe2019circadian
**Ref:** Kolbe et al., Molecular Metabolism 2019 (30:140-151). DOI 10.1016/j.molmet.2019.09.012. Type: biology (SCN clock knockout, metabolism). Lab: Oster (Luebeck).
**Source:** full text

## Question and design
Which role do the SCN pacemaker and the peripheral clock network play in energy metabolism? SCN-enriched Bmal1 deletion (Syt10-Cre x Bmal1 flx/-, SCN-BKO) vs control, 10-week-old males, C57BL/6J, in 12:12 LD or constant darkness (DD), then DD with 12 h restricted feeding (DD-RF, food 06:00-18:00). **Sampling:** timed, four time points per cycle (ZT1-based in LD, CT-based in DD after 1 week), liver, pancreas, epididymal WAT by qPCR, n = 3-5 per time point and genotype; DD-RF tissues at only two opposite times (07:00 and 19:00). Metabolic phenotyping n = 5-6 per group over 10 weeks.

## Method (if any)
qPCR of Bmal1, Per2, Dbp, Rev-erb-alpha and metabolic genes, CircWave rhythm test per tissue and genotype; indirect calorimetry, MRI adiposity, glucose and insulin tolerance. No phase inference beyond CircWave.

## Key results
- LD: "SCN pacemaker function is dispensable for rhythmic behavior and tissue clock gene expression in light-dark conditions" but temporal differences were blunted in SCN-BKO tissues, "indicating rhythm dampening or phase-shift".
- DD: "After one week in constant darkness, time-of-day effects in clock gene expression were conserved in control animals" (p < 0.05, all genes and tissues) while in SCN-BKO no clock transcript was rhythmic in any tissue (p > 0.05); metabolic gene rhythms also lost.
- Metabolic consequence in DD only: more weight gain at equal food intake, more subcutaneous (p = 0.0409) and visceral (p = 0.0252) fat, impaired glucose tolerance (p = 0.0252), higher glucose (p = 0.0395), leptin (p = 0.0437) and triglycerides (p = 0.0304); none of this in LD.
- Time-restricted feeding rescues: "Restoring peripheral clock rhythmicity and synchrony by time-restricted feeding normalizes body weight and glucose metabolism." Per2 and Bmal1 phasing was comparable between genotypes under DD-RF; SCN-BKO activity onset sat 8.5 h before food vs 4.0 h in controls (Fig S3), i.e. feeding aligns the periphery while behaviour is anticipatory and SCN-independent.

## Not addressed / limits
- Not a misaligned light and feeding experiment: LD and feeding are never set in conflict. The desynchrony is loss of rhythmicity in the SCN-deficient animal, and synchrony is not measured directly (only 2 time points per tissue under DD-RF; no phase in hours between tissues).
- Group-mean sampling cannot separate arrhythmia from inter-animal phase dispersion (cf. Sinturel 2021, Welz 2019).
- Syt10-Cre also deletes Bmal1 in scattered other CNS sites; mouse only; males only.

## Link to us
- **Intro claim(s):** the refs_*.md tag (tissue clock network desynchrony from misaligned LD and feeding promotes weight gain) is OVERSTATED/WRONG on the design: it is SCN clock deficiency in DD, rescued by feeding. The metabolic claim (loss of coordinated rhythms goes with weight gain and glucose intolerance in DD) is supported.
- **Relation to clocksync:** background; feeding can serve as the systemic driver of peripheral phase when the SCN is absent.

## Ideas
- Cite as the mouse precedent that restoring rhythm and synchrony of peripheral clocks by feeding normalises metabolism (health consequence of desynchrony), not for misalignment.
- Missing: Husse 2014 and Izumo 2014 (SCN-specific Bmal1 loss), Oster 2019-type reviews of the federated clock model; the same lab's Heyde 2022 is the real misalignment study (card heyde2022induction).
