# heyde2022induction
**Ref:** Heyde & Oster, Scientific Reports 2022 (12:1601). DOI 10.1038/s41598-022-05624-x. Type: biology (zeitgeber misalignment). Lab: Oster (Luebeck).
**Source:** full text

## Question and design
How does acute misalignment of light and feeding zeitgebers desynchronise clocks across and within tissues? Mice entrained to LD 12:12 (LD-24), then moved to a 28 h LD cycle (14:14) combined with 24 h feeding-fasting (FF-24, 12 h food in the previous dark phase) or 28 h FF (FF-28, food in the 14 h dark). Light 300 lx or 3 lx. **Sampling:** timed; liver, adrenal and epididymal WAT sampled across 28 h on day 1 and on day 4 (the day of ~12 h LD phase displacement), n = 3-5 per time point, phase shifts n = 6-7; 7 clock genes (Bmal1, Per2, Dbp, Cry1, Nr1d1, Per3, Nr1d2) per tissue. Temporal resolution not stated in the results (the authors call it relatively low).

## Method (if any)
CircWave fits (sine, 25.8 h wavelength) per gene; phase shift = peak on day 4 minus day 1 per gene, relative to the shift of locomotor activity onset (SCN proxy); within-tissue phase coherence = overlap area of gene peak distributions relative to Bmal1 peak, day 4 vs day 1.

## Key results
- Locomotor period lengthens but does not entrain: 26.64 h (300 lx) vs 25.86 h (3 lx) under FF-24; 25.69 h under FF-28 (insensitive to feeding).
- Peripheral clocks lag the SCN proxy: day 4 minus day 1 phase delay under FF-24 vs FF-28: liver 2.6 vs 5.8 h, adrenal 3.4 vs 5.8 h, eWAT 1.6 vs 7.1 h (reference: same tissue on day 1). "In LD-28/FF-24, peripheral tissue clocks phase-shifted significantly less compared to locomotor activity" so across-tissue (SCN to periphery) desynchrony forms within 4 days.
- "Phase coherence across tissues was disrupted under LD-28/FF-24 while it was maintained to a higher extent under LD-28/FF-28 conditions." Under FF-28 adrenal and eWAT (not liver) shift more than activity.
- Within-tissue coherence of the clock genes falls more under FF-28: liver 73.5% (FF-28) vs 91.5% (FF-24), adrenal 87.6% vs 94.7%, eWAT 59.5% vs 88.4%. "within-tissue phase coherence of clock gene rhythms was consistently more affected" under FF-28.
- Strongest feeding effect is in eWAT (7.1 h vs 1.6 h), not liver; adrenal responds to both light and food.

## Not addressed / limits
- Abstract says cross-tissue coordination fell under FF-28, while Results describe within-tissue (gene-to-gene) coherence falling; read the Results.
- Only three tissues, one reference (day 1 and activity onset), mice are not the same animals across days (population sampling).
- "A limitation of this study is that clock gene tissue profiles were obtained at relatively low temporal resolution" (Discussion); no human data, no hormones measured (glucocorticoids, leptin invoked only as candidates).

## Link to us
- **Intro claim(s):** misaligned zeitgebers produce internal desynchrony within days: supported (4 days, liver/adrenal/eWAT). Evidence that feeding is the stronger systemic driver of peripheral phase and that coordination across and within tissues dissociate.
- **Relation to clocksync:** background; the relevant quantity (offset of a tissue relative to the host/SCN proxy) is what Fig 3d reports in hours.

## Ideas
- Test hypothesis for our human data: tissues with strong feeding input (adipose, liver) versus tissue less feeding-driven would carry larger offsets under shift-work-like misalignment; related human forced-desynchrony blood test already done: A156 (GSE48113, core-cue change not detected, MDE 13.1 h).
- Check whether Heyde data are deposited (no accession in the text).
