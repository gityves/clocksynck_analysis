# woller2016mathematical
**Ref:** Woller, Duez, Staels, Lefranc, Cell Rep 2016 (17:1087-1097). DOI 10.1016/j.celrep.2016.09.060. Type: model. Lab: Lefranc (Lille) with Staels.
**Source:** full text

## Question and design
Mechanistic model of the mouse liver clock driven by the feeding-fasting cycle through two metabolic sensors, AMPK (activity pattern from AMP) and SIRT1 (NAD+ level). 16 ordinary differential equations (mRNA and protein of Bmal1, Per, Cry, Rev-Erb, Ror, Nampt; Dbp mRNA; NAD+), 96 kinetic constants fitted to published liver profiles (Hughes 2009, 1 point/h, mice in constant darkness fed ad libitum; NAD+ and AMP from Hatori 2012), then checked on Sirt1 KO (Bellet 2013) and Lkb1 KO (Lamia 2009). Simulations of normal feeding, fasting-like (AMPK constantly high) and high-fat-diet-like (AMPK constantly low) states, Cry1 KO rescue, and timed REV-ERB agonist pulses.
**Sampling:** none new; the model is calibrated on literature time series (mouse liver, 24 h).

## Method (if any)
Deterministic ODE model of one entrained liver oscillator; forcing = a parametric two-pulse AMPK activity profile (peaks near ZT5 and ZT17 consistent with AMP data) and PGC1a activity. Outputs: phase and amplitude of clock genes and NAD+. The SCN input is dropped: "we neglected the influence of systemic signals from the SCN" and "This is legitimate as long as only the entrained regime is studied, not the resetting dynamics."

## Key results
- Premise used for liver: "when these two cues are discordant, the entrainment phase is determined by food timing" (Damiola 2000, Saini 2013), so the model includes feeding only.
- Flattening the metabolic drive lowers clock amplitude. Two conditions with a constant AMPK activity (fasted-like high, fed-like low, mimicking high-fat feeding) reduce "the amplitude of oscillations" significantly in both; decreased AMPK activity cuts Nampt mRNA, protein and NAD+, and constitutive high AMPK also shifts the phase of NAD+ and Nampt profiles. This reproduces the dampened clock gene and NAD+ rhythms in ad libitum high-fat-fed mice (Hatori 2012, Eckel-Mahan 2013).
- Cry1 loss (Cry transcription halved) restores the amplitude of NAD+ and of most clock genes under reduced AMPK rhythms; a timed REV-ERB agonist pulse can restore amplitude and phase.
- The model predicts a smaller non-circadian NAD+ peak that vanishes with both constitutively high and low AMPK activity.
- No model result on inversion, conflicting cues or phase shift per tissue. Resetting is explicitly deferred: "The next step is to study how it resets when AMP rhythms are phase-shifted from day to day".

## Not addressed / limits
- Liver only, entrained regime only; no SCN or glucocorticoid input; no inversion or shift simulations; no tissue-specific weights.
- Amplitude loss here arises from a flat or low-amplitude feeding drive acting on the clock network with parameters tuned to reproduce it; the text does not quantify how much of the amplitude is self-sustained.
- Drive is a prescribed AMPK profile, not feeding itself; no numbers for amplitude ratios in the main text.

## Link to us
- **Intro claim(s):** background for liver clock entrained by feeding-fasting; partly contradicts the saturation reading: in this model removing the metabolic rhythm lowers clock amplitude significantly (a drive-dependent oscillator), whereas data on arrhythmic feeding (Greenwell 2019) keep clock-gene amplitude near 0.7 (our A297 number). The weight of drive versus self-sustained amplitude is a free quantity here.
- **Relation to clocksync:** background; no data. Its AMPK-amplitude to clock-amplitude map is the only mechanistic prediction found for amplitude under loss of the feeding rhythm.

## Ideas
- Simulate the model with the AMPK profile at graded amplitude (1, 0.5, 0) and read the clock amplitude: if amplitude falls only to a plateau, it supports self-sustained plus saturation; if it falls to near 0, it conflicts with Greenwell (0.724 in A297).
- Shift the AMPK profile by 12 h and add a fixed SCN term (as in woller2018modeling) to get amplitude and phase for inversion in a liver-specific parameter set.
- Reference to read: Hughes 2009 (liver, 1 h resolution, constant darkness AL; also light-restricted feeding data) as test data for amplitude.
