# yamazaki2000resetting
**Ref:** Yamazaki, Numano, Abe, Hida, Takahashi, Ueda, Block, Sakaki, Menaker, Tei, Science 2000. DOI 10.1126/science.288.5466.682. Type: biology. Lab: Menaker / Tei.
**Source:** full text (the extracted text begins with the tail of an unrelated Science report; the paper starts at "Resetting Central and Peripheral Circadian Oscillators")

## Question and design
Per1-luciferase transgenic rats (mouse Per1 promoter), LD 12:12. Explants of SCN, skeletal muscle, liver and lung cultured in constant
darkness and 36 C, bioluminescence recorded continuously (SCN n = 62). The LD cycle was then shifted by 6 h (advance or delay; also 9 h for
SCN and muscle only) and animals were killed on the first or sixth complete cycle after the shift; locomotion recorded in vivo (n = 10).
**Sampling:** one explant per tissue per rat, killed 30-60 min before lights-off (single in-vivo time per rat); the phase is the first peak
of the in-vitro rhythm, i.e. phase is read from cultures, not from fresh tissue. Timed by design (known shift, known cycle after shift).

## Key results
- Peripheral oscillators lag the SCN in steady state by 7-11 h, and damp in culture. "Liver, lung, and skeletal muscle all showed circadian rhythms of light output that phaselagged the SCN rhythm by 7 to 11 hours." (Results)
- Table 1 peak phases (h after last lights-on): control SCN 30.9, muscle 41.3, liver 42.1, lung 38.1, i.e. muscle 10.4 h, liver 11.2 h, lung 7.2 h after SCN. "Circadian oscillation in muscle and liver peaked about 10 to 11 hours after" the SCN, lung about 7 h (Fig. 4, Table 1)
- 6 h delay, first cycle: SCN -4.9 h, muscle -3.9, lung -3.6, liver +1.1 (not shifted). "Liver, on the other hand, was not delayed at all on the first complete cycle." (Results; 4 of 8 liver cultures arrhythmic)
- 6 h delay, sixth cycle: SCN -6.0, muscle -6.6, lung -6.1, liver -3.5; "By the sixth complete cycle, the liver cultures had only delayed by about 3.5 hours." Liver still ~4 h delayed at cycle 16 (n = 4).
- 6 h advance, first cycle: SCN +5.0, liver +2.0, lung -1.1 (i.e. wrong direction, n.s.), muscle arrhythmic in 6 of 8. "On the other hand, none of the six liver cultures became arrhythmic following the advance, although they only shifted by about 2 hours." Sixth cycle: SCN +6.0, muscle +5.4, liver +6.3, lung +5.5, all complete.
- 9 h shifts: SCN almost fully shifted at cycle 1 (delay -9.3, advance +7.8); muscle arrhythmic after advance, only ~3 h shifted after delay.
- Behaviour lags the SCN: re-entrainment took 6.2 +/- 0.5 cycles (advance) vs 1.8 +/- 0.3 (delay). "The SCN shifted almost completely after the first cycle, but locomotor rhythmicity shifted more slowly" (Results)
- Transient internal misalignment: "SCN, muscle, and lung (but not liver) resumed their normal phase relationships to each other and to the light cycle within 6 days of the shift in either direction." (Discussion)

## Not addressed / limits
- Reference is the SCN (and the light cycle); no peripheral-vs-peripheral contrast in the same animal beyond the table; tissues from different rats at each time.
- Ex vivo peak phase of cultures; muscle and lung data from advanced animals are sparse (n of rhythmic cultures small, arrhythmia common).
- One reporter gene (Per1-luc), no multi-gene phase structure; rat; only 4 tissues; no sex or age variation; no disease.
- The liver advance/delay asymmetry is reported as unexpected and not explained (feeding was not controlled).

## Link to us
- **Intro claim(s):** supported. refs_circadian_desynchrony.md says SCN re-entrains within days while peripheral tissues re-entrain slowly at different rates: correct (liver slowest, about 3.5 h of 6 h after 6 cycles).
- **Relation to clocksync:** background. It gives the size of transient inter-tissue misalignment (up to ~6 h at cycle 1 of a 6 h shift) and a steady-state SCN-to-tissue spread of 7-11 h. Fits the offset cancelled by constant host reference logic only if all tissues are read at one clock time.

## Ideas
- Use Table 1 as a magnitude benchmark for simulated tissue-offset changes (3-6 h): done in part, A124 (simulated night shift in muscle, MDE 7.3 h); a mouse jet-lag test is not in METHOD.md, registered_experiments.md or results_preliminary.md (no "jet lag" hit).
- Scout public multi-tissue mouse jet-lag time courses (GEO) as an out-of-sample test of tissue offsets relative to a neuronal/SCN host; Yamazaki-style, one tissue per animal would not pair tissues, so pairing matters (cf. Manella pairing issue in A111).
