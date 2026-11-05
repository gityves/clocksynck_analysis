# tahara2012in
**Ref:** Tahara et al., Current Biology 2012 (22:1029-1034). DOI 10.1016/j.cub.2012.04.009. Type: biology (in vivo reporter imaging). Lab: Shibata (Waseda).
**Source:** full text

## Question and design
Do peripheral clocks of the intact mouse keep oscillating without the SCN? PER2::LUC and Bmal1-ELuc mice imaged non-invasively (kidney, liver, submandibular gland) after luciferin injection. **Sampling:** timed, six imaging time points at 4 h intervals over 24 h (ZT11-ZT7) in LD; SCN-lesioned (SCNx, n = 10) vs sham (n = 7) imaged once ~1 month after surgery in DD (starting 19:00); RF cohort n = 8.

## Method (if any)
In vivo photon imaging of one reporter per tissue, cosinor with an index of goodness of fit (g < 0.1 and amplitude threshold = rhythmic). Gives peak phase, amplitude, period per organ per mouse. Not a population inference method; needs a time series per animal.

## Key results
- Intact mice, LD, PER2::LUC peak ZT 18.7 (kidney), 18.3 (liver), 18.7 (submandibular); Bmal1-ELuc peaks ZT 2.4, 3.1, 3.6: nearly antiphase to PER2, so the three tissues sit within ~0.4 h of each other in LD (n = 3 and n = 4). "The phase in Bmal1-ELuc mice was almost antiphase to that of PER2::LUC mice and also to that of locomotor activity" (Fig 2C).
- SCNx: about half of the organs still oscillate, amplitude lower. "Circadian oscillations were detected in the kidney, liver, and submandibular gland studied in about half of the SCN-lesioned, behaviorally arrhythmic mice." (Summary)
- Peak phases of rhythmic SCNx organs: 14.8 h (kidney), 16.4 h (liver), 15.6 h (submandibular) vs sham 15.3, 15.5, 15.8 h; no difference, but wider SEM. "The peak phase of each organ (under 0.1 of goodness of fit) was similar among sham and SCNx mice within about 8 hr ranges, but bigger SEM was seen in SCNx group" (Fig 3F). The arrhythmic organs had a broader phase range.
- Restricted feeding (12 h, 20:00-08:00, 12 days) entrained SCNx peripheral clocks; two weeks after release into free feeding, rhythms persisted but phases spread: "Among individual mice, the range of phase difference was narrow under the RF condition but was broader after 2 weeks of FF" (Fig 4E). Free-running periods 23.7 h in all three organs (no difference), so inter-organ phase stays similar in a given mouse.
- Daytime feeding moved peripheral peak phase ~12 h relative to night feeding (Fig S4C-E, stated in text).

## Not addressed / limits
- Only three organs, one reporter each, never measured simultaneously in the same cell type or at higher-than-4 h resolution; no human tissue.
- Phase differences between organs in SCNx are shown as peak-phase scatter at one time, not as a coupling estimate; arrhythmic organs (g > 0.1) are excluded from the phase comparison.
- Does not test jet lag or shifted feeding in SCN-intact mice (only LD night-feeding control, Fig S4A).

## Link to us
- **Intro claim(s):** inter-tissue relations measured by reporters; peripheral clocks can run without the SCN (supports). The refs_*.md tag (different tissues re-entrain at different rates after jet lag and feeding shifts) is NOT what this paper tests: CLAIM WRONG (it is SCN lesion, dampened but oscillating clocks, RF entrainment).
- **Relation to clocksync:** background. Precedent that organ clocks in SCNx mice have mouse-specific, widely dispersed relative phases, i.e. the quantity (offset between organs within one animal) clocksync reads.

## Ideas
- Cite for the contrast that organs stay near-identical in phase when entrained (kidney, liver, gland within ~0.4 h in LD) but spread when the SCN is removed: a benchmark for the expected size of between-tissue offsets (compare with A36 baboon organ offsets 0.45 h).
- Missing reference: Yoo 2004 (already in notes) and Husse 2014 (SCN-specific Bmal1 loss, liver clock in DD).
