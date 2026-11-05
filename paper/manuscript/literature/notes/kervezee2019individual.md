# kervezee2019individual
**Ref:** Kervezee L, Cermakian N, Boivin DB, PLoS Biol 2019;17(6):e3000303. DOI 10.1371/journal.pbio.3000303. Type: clinical (human misalignment protocol, metabolomics). Lab: Boivin / Cermakian (McGill).
**Source:** full text

## Question and design
Within-subject laboratory study, n = 9 healthy (8 men, 1 woman; 18-29 y). Targeted plasma metabolomics (132 metabolites detected, 130 analysed), 229 samples. Single compartment: plasma metabolites.
**Sampling:** timed, every 2 h for two full 24 h periods: baseline and day 4 of a 4-day simulated night shift (habitual sleep delayed by 10 h; hourly isocaloric snacks while awake, fasting in the 8 h sleep period; not a constant routine). Time expressed relative to the individual habitual bedtime.

## Method (if any)
None new: mixed-effects cosinor per condition, model selection to classify each rhythmic metabolite as circadian-influenced (stays at the unshifted endogenous phase) or behaviour-influenced (follows the shifted sleep-wake and fasting-feeding cycles), individual cosinor and Rayleigh tests of phase shifts. The endogenous clock reference is the melatonin phase from the same protocol in earlier papers (refs 21, 22, 26), not measured here.

## Key results
- Rhythmic metabolites: 51 at baseline and 53 in the night shift condition (FDR below 0.05), 32 in both.
- Internal desynchrony at group level: "Twenty-four metabolites were identified as behaviour-influenced, while 7 were identified as circadian-influenced", the behaviour-influenced ones with "an average phase delay of 8.8 ± 0.62 h" (versus 10 h shift); circadian-influenced metabolites showed no significant shift (Rayleigh p = 0.16) and unchanged amplitude.
- Different outputs shift by different amounts: the authors cite their own earlier findings that "the phase of central circadian clock markers, peripheral clock gene expression, and transcriptome-wide gene expression rhythms do not shift in response to the night shift condition" while most metabolites shifted, so the "phase angle between different physiological rhythms is disturbed".
- Individual heterogeneity: individual average metabolite phase shifts "ranged from a 0.2-h advance in one subject to a 12.0-h delay in another subject", with a significant mean direction in 8 of 9 subjects; mean across subjects a 6.8 h delay (p = 0.0508, not significant). The number of rhythmic metabolites per subject ranged 6 to 56, pairwise overlap only about 21-22%.
- Levels: 27 metabolites changed in mean level; lysophospholipids lost rhythm.

## Not addressed / limits
- Only plasma; no tissue-resolved read-out. The authors: "daily variations in metabolite levels have been shown to be highly tissue specific". So this is not evidence of between-tissue clock offsets, only of behaviour-driven versus clock-driven outputs.
- No clock-gene or phase of any peripheral clock in this paper; the central-clock reference is from other papers on the same protocol.
- Small n (9, mostly young men), 4 days only, laboratory protocol, no real shift workers; non-constant routine (food and posture).
- Individual phase shifts are computed from few rhythmic metabolites per subject (uncorrected p < 0.05).

## Link to us
- **Intro claim(s):** human misalignment (night shift) separates behaviour-driven and clock-driven rhythms, with strong individual variability. Our refs line 75 (plasma metabolomic signatures, individual-specific phase shifts) is supported. The grouped line 102 (human protocols reading out two compartments in the same volunteers) is not right for this paper: it has one compartment (plasma) and no second tissue.
- **Relation to clocksync:** background (motivation for per-donor offsets); not data we can use (no clock-gene tissue data).

## Ideas
- The 0.2 h to 12 h per-person spread of shifts is a ready benchmark for the size of a donor-specific offset an internal-desynchrony method should distinguish from noise.
- Paired data from the same protocol (Kervezee 2018 PBMC transcriptome, Boivin 2022 PNAS style data) give blood-transcriptome phases next to plasma metabolite phases in the same individuals: a within-person two-compartment test for clocksync.
- Missing reference: Skene 2018 PNAS (human night shift, metabolome, 95% shifted) is cited but not yet read here.
