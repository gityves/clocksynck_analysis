# czeisler1999stability
**Ref:** Czeisler, Duffy, Shanahan, Brown, Mitchell, Rimmer, Ronda, Silva, Allan, Emens, Dijk, Kronauer, Science 1999 (284:2177). DOI 10.1126/science.284.5423.2177. Type: clinical (human physiology). Lab: Czeisler (Harvard/BWH).
**Source:** full text

## Question and design
Is the intrinsic period of the human circadian pacemaker near 24 h and does it shorten with age (the proposed cause of the advanced phase of older people)? 11 healthy young men (mean age 23.7 y) and 13 healthy older subjects (9 men, 4 women; mean age 67.4 y), 29 to 38 days each in a time-cue-free laboratory with controlled dim light. Forced desynchrony: a 28-hour "day" (bedtime 4 h later each day for about 3.5 weeks) so that sleep-wake and light exposure are spread over all circadian phases. A few subjects also did 20-h, 42.85-h or classical free-running protocols. **Sampling:** continuous core body temperature, plasma melatonin and cortisol (blood every 20 to 60 min in segments); human in vivo, not tissue.

## Method (if any)
Period of each rhythm estimated by nonorthogonal spectral analysis (maximum-likelihood fit of the imposed 28-h and the sought endogenous period plus harmonics). Subject period = mean of the available temperature, melatonin and cortisol periods.

## Key results
- Period is near 24 h and age independent. "the intrinsic period of the human circadian pacemaker averages 24.18 hours in both age groups, with a tight distribution consistent with other species" (abstract). Young 24.18 +/- 0.04 h (PCV 0.54%), older 24.18 +/- 0.04 h (PCV 0.58%), P = 0.961 (text, Table 1); nearly 90% of the 24 estimates between 24.00 and 24.35 h.
- Power limit: detectable difference "greater than 9 min between the young and older subjects" (alpha 0.05, power 0.90, SD 0.15 h); the authors state a much larger series might detect a small difference.
- The three rhythms give the same period within a person (temperature vs melatonin r = 0.951, vs cortisol r = 0.982, melatonin vs cortisol r = 0.984).
- Older subjects were phase advanced anyway. "the older subjects in this study exhibited the characteristically earlier entrained circadian phase and earlier morning awakening typically found in this age group" relative to young subjects; so it is "unlikely that the systematic" advance with age comes from a shorter period (text, discussion). No phase difference in hours is given in the text.
- Classical free-running gave longer periods because of self-selected light: two subjects had 25.1 and 25.0 h on free-run versus about 24.3 h on forced desynchrony; Kronauer-model simulations show room-light feedback can lengthen the apparent period "by more than 0.7 hours".
- Alternatives named for the age advance: "age-related changes in the sleep-homeostatic process and its interaction with the circadian and entrainment processes".

## Not addressed / limits
- Healthy, screened older volunteers (13) and young men only (no young women); small n, power limit of 9 min. Abnormal period in some older individuals not excluded (reported elsewhere).
- Zeitgeber strength, central amplitude and light exposure in daily life are not measured; the phase angle of entrainment in older people is not modelled here.
- Phase advance size is not quantified in the text (Methods and ref 23 are elsewhere).

## Link to us
- **Intro claim(s):** the human pacemaker has tau about 24.2 h, so tau > T = 24 h and a lag of the internal clock behind the zeitgeber is expected; period does not change with age.
- **Relation to clocksync:** background for A234 (internal time runs ahead of clock time with age, +0.64 h per decade in body, +0.99 h in brain). Period shortening is not the explanation in a healthy cohort: tau is the same, so the advance must come from the other terms of the phase of entrainment (zeitgeber strength or relative light, oscillator amplitude) or from behaviour/sleep homeostasis.

## Ideas
- In the entrainment model, phase of entrainment = f(tau - T, zeitgeber strength / amplitude): with tau = 24.18 h fixed across ages, test whether a lower central amplitude (the damped older rhythms in Talamanca 2023) gives an advance of the right size per decade; check with Grabe 2022 (Fig. 3: large periods lead to earlier phases, stronger zeitgebers weaker slope).
- Quote the 9 min detection limit: a small age shift of tau cannot be excluded by this study; compute in the model what tau change would be needed for +0.64 h per decade at fixed amplitude and compare with it.
- The phase advance of the older group in hours is not in this text; look for it in the cited reference 23 (not in our export).
