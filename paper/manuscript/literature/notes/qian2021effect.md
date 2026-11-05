# qian2021effect
**Ref:** Qian DC et al. (MEMOIR; Buchwald), Lancet Oncol 22:1777-1786, 2021. DOI 10.1016/S1470-2045(21)00546-5. Type: clinical (retrospective cohort). Lab: Buchwald, Emory Winship Cancer Institute.
**Source:** full text

## Question and design
Does clock time of immune checkpoint inhibitor (ICI) infusion affect overall survival (OS) in stage IV melanoma? Single tertiary centre, 2012-2020, 481 treated, 299 with stage IV and at least four infusions; median follow-up 27 months; 66 % male, median age 61. **Sampling / timing:** infusion clock times recorded retrospectively from the clinic record; infusions occur in clinic hours only (clinic closes at 18:00); no randomisation, and patients chose their own slots. Exposure = proportion of infusions after 16:30 h.

## Method (if any)
Cox regression with exposure per 20 % of infusions after 16:30; propensity-score matching (age, ECOG, LDH, corticosteroids, radiotherapy) of patients with at least 20 % evening infusions vs fewer. Cutoff 16:30: "midpoints (1630 h) as a biologically informed cutoff to" denote evening, derived from earlier studies that used 13:00-16:00, 15:00-17:00, 18:00-18:30 and 14:00-20:00. No circadian phase measured.

## Key results
- Continuous: each additional 20 % of infusions after 16:30 gives OS hazard ratio "of 1·31 (95% CI 1·00 to 1·71; p=0·046)" (n = 299).
- Matched (n = 73 vs 73): evening group shorter OS, "HR 2·04 [1·04 to 4·00; p=0·038]"; median OS 4.8 years vs not reached; multivariable HR 1.80 (1.08-2.98, p = 0.023).
- Adverse events (colitis 18 %, hepatitis 9 %, hypophysitis 5 %) balanced between groups; no treatment-related deaths.
- Infusion time pattern is not random: "Patients who received at least 20% of infusions after" 16:30 were more likely to have had an evening first infusion (50 % vs 14 %).
- Sensitivity: associations stable for cutoffs near 16:30 (appendix); seasonal start not associated with OS.

## Not addressed / limits
- "Limitations of our study include its single institution" and retrospective nature, dispersion of infusion times, and confounding by retrospective partitioning (authors' words).
- "inconsistent across patients because patients have" control over their schedule (evening slots may select working or less sick patients); number of infusions excluded as covariate; exclusion of < 4 infusions (immortal-time and survivor concerns).
- Clock time is not internal time: no chronotype, melatonin, tumour or immune-cell clock data. No mechanism tested.
- One cancer, one centre; only 54 + 73 comparisons in matched analysis; hypothesis-generating.

## Link to us
- **Intro claim(s):** our current claim (refs_*.md): MEMOIR, 299 stage IV melanoma, OS HR 1.31 per additional 20 % infusions after 16:30; propensity-matched 73 vs 73 shorter OS with evening infusions. Supported (matches Summary). The refs table wording (infusion before 16:30 associated with better overall survival) is the same finding inverted and is fine, but should be worded as association in a retrospective single-centre cohort.
- **Relation to clocksync:** background (cancer interaction: clock time of treatment matters); evidence that clinical practice is organised at clock times.

## Ideas
- Our framing: such studies use clock time as a proxy for internal time; a per-patient internal phase / coupling readout (needs tissue or blood clock genes) could explain the scatter. Pair with Landre 2024 meta-analysis (already cited).
- A prospective check would need blood or biopsy pseudobulks per donor before treatment, tested against the evening-infusion benefit.
