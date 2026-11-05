# saini2013realtime
**Ref:** Saini, Liani, Curie, Gos, Kreppel, Emmenegger, Bonacina, Wolf, Poget, Franken, Schibler, Genes Dev 2013. DOI 10.1101/gad.210773.112. Type: method / biology. Lab: Schibler.
**Source:** full text

## Question and design
How does the SCN affect the phase of hepatocyte clocks when feeding conflicts with it? New device, the RT-Biolumicorder (round cage, reflecting cone,
photomultiplier), records bioluminescence of liver reporters (adenoviral Adv-Bmal1-luc and Adv-Rev-erba-luc in hairless mice; PER2::LUC knock-in) in
freely moving mice for many days, with luciferin from osmotic minipumps. Feeding inverted from night to day. Sham-operated (n = 3 recordings
shown) vs SCN-lesioned mice (n = 3 recordings shown, plus a fourth ad libitum case).
**Sampling:** continuous (1-min photon counts) for about 1 week per mouse, individual animals (not pooled), liver only.

## Method (if any)
Reports hepatocyte phase continuously per mouse (peak times by Gaussian fit, phase shift against a 24 h grid); calibrated to expected Bmal1, Rev-erba and Per2 phases. Needs a transduced liver or reporter mouse.

## Key results
- Gene-specific resetting by feeding: 7 days after inversion "the phase of Bmal1-luciferase expression was shifted by only 6–7 h", whereas Rev-erba-luc and PER2::LUC "were almost completely inverted (by" 12 h; about 60 degrees for Bmal1-luc vs 120 degrees for Rev-erba-luc vs 150 degrees for PER2::LUC at day 7 (Fig. 4).
- "These results indicated that different circadian genes responded with different kinetics to food entrainment and that the phase relationships between core clock components were uncoupled during the transition period." (Results)
- SCN-lesioned mice invert quickly: the lesioned mice "reached the new steady-state phase of hepatic Bmal1-luciferase expression within 1–2 d after inversion of the feeding rhythm" vs 6-7 h in 7 d in sham mice (Fig. 5).
- Summary: "hepatocyte oscillators of SCN-lesioned mice synchronized more rapidly to feeding cycles than hepatocyte clocks of intact mice." (Abstract); the SCN counteracts the feeding signal when its phase is in conflict.
- Persistence after SCN lesion: reporter rhythms continue for several days, with lower Bmal1-luc amplitude in lesioned mice fed ad libitum (possible partial desynchrony).
- The authors cite earlier results for "phase coherence is lost between different organs" in SCN-lesioned PER2::LUC mice (Yoo 2004; Tahara 2012).

## Not addressed / limits
- One organ (liver) per animal; no comparison of liver with another tissue in the same mouse; phase referenced to the feeding cycle and to a pre-inversion peak.
- Small n (3 sham, 3 lesioned recordings); some recordings with a prototype device; viral transduction and anaesthesia may reset hepatocyte phase (the authors state this caveat).
- Does not measure the SCN phase itself; sham/lesion difference in the steady-state phase (lesioned mice slightly delayed) is noted without hours.

## Link to us
- **Intro claim(s):** supported with adjustment. refs_circadian_desynchrony.md reads liver clock phase set by feeding; hepatocyte clocks phase-set by systemic cues. The main point is that the SCN slows the feeding-driven shift of liver clocks (brake), and that genes in the same tissue move at different speeds (6-7 h vs 12 h of a 12 h inversion at day 7).
- **Relation to clocksync:** background; supports the host = what the SCN and other tissues share reading, and the gene-wise heterogeneity caution.

## Ideas
- Not in docs: a simulated tissue where Bmal1 and Rev-erba/Per2 gene phases differ by 5-6 h (as at day 7 of an inversion) to see how the coherence gate and the compartment phase behave; A111/A112 only tested feeding and clock loss.
- Reference we lack: Asher et al. 2010 (Cell, PARP1 and food phase entrainment), cited here for the kinetics paradigm.
