# buhr2010temperature
**Ref:** Buhr, Yoo, Takahashi, Science 330:379 (2010). DOI 10.1126/science.1195262. Type: biology. Lab: Takahashi (Northwestern / UT Southwestern).
**Source:** full text (main text only; supplementary figures S1-S8 and Methods not in the export; the extracted file also starts with the reference list of the preceding paper)

## Question and design
Can a small temperature change within the physiological range (36 to 38.5 C, the mouse body-temperature range) reset circadian oscillators,
and is the SCN resistant to it? Mouse Per2::Luciferase explants (SCN, pituitary, lung and other peripheral tissues) cultured ex vivo in
Lumicycle photomultiplier machines; single-cell SCN imaging with an intensified CCD; drugs (TTX, nimodipine, mibefradil, KNK437). n = 4 to 14
explants per point (e.g. 4 SCN and 5 pituitary or lung in Fig. 1B). **Sampling:** continuous bioluminescence recording for days, temperature
pulses at all circadian phases (phase transition curves); no tissue is sampled in vivo, no time series of endogenous genes.

## Method (if any)
Not a method paper. Phase transition curves (new vs old phase after a 1 h or 6 h pulse to 38.5 C from a 36 C baseline) and entrainment to
12 h:12 h 36/38.5 C square waves in opposite phases.

## Key results
- Peripheral explants show type 0 resetting to a common phase: "pulses and exhibited high-amplitude, type 0 resetting." New phases cluster around CT12-18 (Fig. 1A).
- Opposite-phase temperature cycles split lung and pituitary by about 180 degrees within 3 days; the SCN does not move: "the SCN is resistant to cyclic temperature changes" (Fig. 1B).
- Resistance is a network property: TTX makes the SCN reset like a peripheral tissue ("of the voltage-gated Na+ channel blocker, tetrodotoxin (TTX), temperature pulses strongly reset"); n = 14 vehicle, 11 TTX, 13 nimodipine, 11 mibefradil (Fig. 2F).
- Ventral or dorsal SCN halves cultured alone reset like TTX-treated SCN, whereas "The sagittally separated SCN were resistant" (Fig. 3).
- Heat-shock pathway: KNK437 (HSF inhibitor) mimics cool pulses; in DMSO Q10 = 1.04 (SCN) and 0.97 (pituitary), with KNK437 "KNK437, the Q10 values became 1.34 for the SCN" and 1.38 for pituitary.
- Conclusion: "potent cue for resetting in peripheral tissues." Body-temperature rhythms driven by the SCN can therefore synchronise peripheral oscillators.

## Not addressed / limits
- Everything is ex vivo explants of Per2::Luc mice; no in vivo temperature manipulation. The in vivo argument is inferred from earlier work and from the explants.
- The authors themselves note the confound: "it is extremely difficult, if not impossible, to separate the effects of activity, sleep, feeding, and" temperature in vivo.
- Mouse only; no human tissue, no disease, no desynchrony measured; no quantification of how much temperature contributes relative to feeding or glucocorticoids.
- Neonatal and juvenile rat SCN are reported more temperature-sensitive, so resistance is age and species dependent.

## Link to us
- **Intro claim(s):** refs_circadian_desynchrony.md B1/B2 row says small temperature cycles reset the phase of SCN-independent peripheral clocks in culture and in vivo. CLAIM OVERSTATED: culture (ex vivo explants) is supported, in vivo is not shown here.
  Supported: temperature is a plausible internal coupling cue between SCN and periphery (internal synchrony).
- **Relation to clocksync:** background. Mechanistic basis for a shared host clock state that periphery reads (coupling c); a temperature or fever perturbation would shift all compartments together.

## Ideas
- Fever or hypothermia in disease or ICU is a candidate shared offset of all compartments of one donor; test on GTEx donors with ventilator or fever phenotypes (aggregates only): a common shift is cancelled by clocksync, a compartment-specific one is not.
- Missing in vivo reference: Saini 2012 (body-temperature cycles drive Hsf1-dependent circadian gene expression). Brown 2002 (in vivo body-temperature entrainment) is already a card (brown2002rhythms); cite it for the in vivo statement instead of this paper.
