# sinturel2021circadian
**Ref:** Sinturel et al., Genes & Development 2021 (35:329-334). DOI 10.1101/gad.346460.120. Type: biology (in vivo real-time bioluminescence). Lab: Schibler / Dibner (Geneva); not a Naef-lab paper (Naef is not an author).
**Source:** full text

## Question and design
Does the SCN, or any extrahepatic clock, keep peripheral clocks and hepatocytes synchronised? RT-Biolumicorder adapted to weeks-long recording in freely moving mice, luciferin in drinking water. Arms: hairless Per2::Luc (whole-body signal, sham vs SCN lesion, n >= 5 per group); adenoviral hepatic Rev-erb-alpha-luc (and Bmal1-luc) by tail vein, sham vs SCNx (n = 3 per group); "hepatocyte clock-only" mice (Bmal1 stop allele rescued by Alb-Cre; n >= 4) vs Bmal1+/- and Bmal1-/-. **Sampling:** continuous 1-min photon counts for ~2 weeks or more, mostly in DD; feeding regimens (ad libitum, 12 h/12 h, inverted) imposed in some groups.

## Method (if any)
FFT of 120-min binned traces; Gaussian fit around 24 h gives peak amplitude, period and oscillation width (a measure of phase coherence); phase shift after feeding inversion from Gaussian-fitted peaks relative to an ideal 24 h cycle. Whole-body Per2::Luc is a sum over organs, so organ-specific phases are not resolved.

## Key results
- SCN needed for between-organ synchrony (whole-body signal): "In contrast, these rhythms were strongly dampened in SCN-lesioned mice, suggesting that the SCN is required for the synchronization of oscillators among different tissues." (Fig 2A,B,E)
- Feeding is sufficient to resynchronise: "Imposing feeding rhythms synchronized PER2::luc expression in both sham-operated and SCN-lesioned mice" (Fig 2); period and FFT amplitude then similar. After feeding inversion the phase shifted faster in SCN-lesioned than sham mice (Fig S2G; no hours given in the text).
- Hepatocytes stay synchronised without SCN: hepatic Rev-erb-alpha-luc oscillates robustly in SCNx, FFT peak at 24 h not different from sham (n = 3 each); period variability and width slightly larger, not significant.
- Hepatocyte clock-only mice: "Hence, hepatocytes were partially synchronized in behaviorally arrhythmic mice that harbor circadian oscillators exclusively in hepatocytes." Rhythms less stable than in Bmal1+/- or SCNx; dominant 24 h peak in all periodograms.
- Direct conflict with Koronowski 2019: "In contrast to previous reports (Guo et al. 2006; Koronowski et al. 2019), our findings show that the SCN is not required for maintaining phase coherence between hepatocytes." And: "The data presented here thus exemplify the importance of real-time recording when studying cyclic gene expression in behaviorally arrhythmic animals."

## Not addressed / limits
- Small n (3 per group for liver); adenoviral transduction efficiency varies, so amplitudes are not comparable between mice.
- The hours of phase offset between organs in SCNx mice are not quantified; the claim that organs have different phases (Fig 5B) rests on the dampened whole-body sum.
- Mouse; DD; only liver resolved.

## Link to us
- **Intro claim(s):** hepatocyte clocks stay synchronised without SCN (supported); the phrase in refs (cell-cell synchrony inside one tissue) is right, but the paper also supports the SCN as the inter-organ synchroniser, which is the stronger intro point.
- **Relation to clocksync:** background. Key caution for pooled readouts: a group-level flat profile (Koronowski) can come from between-animal phase dispersion; within-animal relations (what clocksync reads) are the safer quantity. Matches the rule that a failed gate is uninformative.

## Ideas
- Cite together with Tahara 2012 for: peripheral clocks without SCN keep oscillating but with animal-specific phases, so a donor-wise readout is needed.
- Missing: Saini 2013 (RT-Biolumicorder, feeding inversion kinetics), Guan 2020 (non-parenchymal liver cells and hepatic rhythms).
