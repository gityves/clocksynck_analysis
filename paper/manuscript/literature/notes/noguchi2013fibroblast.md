# noguchi2013fibroblast
**Ref:** Noguchi, Wang, Welsh, J Biol Rhythms 2013 (28:183-192). DOI 10.1177/0748730413487494. Type: biology (cell density, single-cell imaging). Lab: Welsh (UCSD).
**Source:** full text

## Question and design
Are single fibroblasts independent oscillators, or do they need neighbours? Primary tail-tip fibroblasts from neonatal PER2::LUC mice, serum-free explant medium, single-cell bioluminescence imaging (30 min frames, 5-7 days, 36 C). Perturbations: cell density (low 4-27 cells/mm2 versus high 46-410 cells/mm2), mixing 1-5% PER2::LUC cells with 20-100 fold excess nonluminescent wild-type, Bmal1-/- (arrhythmic) or Cry2-/- (long period) fibroblasts, conditioned medium (CM, 50%) from high-density cultures, high K+ (21 mM) or Ca2+ (3.6 mM). **Sampling:** in vitro, free-running after a medium change; first 12 h excluded.

## Method (if any)
FFT in LumiCycle: rhythmic if peak in 20-36 h and at least 6% of spectral power; period and amplitude from sine fit for rhythmic cells; Mann-Whitney, Kruskal-Wallis tests.

## Key results
- Low density abolishes rhythm: "The percentage of rhythmic cells in a culture increased dramatically as cell density increased"; mean 18.4% rhythmic at low density (4 cultures) versus 75.4% at high density (8 cultures), P < 0.01. Brightness, and period and amplitude of the cells that stay rhythmic, are not significantly different. "Cells that remained rhythmic in low-density culture were not clustered with other rhythmic or nonrhythmic cells", i.e. no local coupling clusters.
- The signal need not be rhythmic: wild-type cells in a low density of rhythmic cells but high total density recover (WT feeders 72.9%, Bmal1-/- feeders 59.2%, Cry2-/- feeders 69.7%; no significant difference, P > 0.05).
- Neighbours do not set period: with Cry2-/- feeders the wild-type period was "about 1 h shorter" but not significantly different (ANOVA, P > 0.05), although Cry2-/- cells have long periods.
- CM from high-density cultures rescues: 66.6% rhythmic (n = 3 cultures), not different from high density; high K+ or Ca2+ do not rescue consistently (high Ca2+ 20-77.8% per culture).
- Authors: "fibroblasts require paracrine signals from adjacent cells for normal expression of rhythmicity, but that these signals do not have to be rhythmic". Previously "no evidence has been found for local circadian communication in fibroblast cultures".

## Not addressed / limits
- Primary mouse tail fibroblasts, in vitro, imaged 5-7 days; no tissue, no in vivo context, no other cell type as the partner. Identity of the paracrine factor not found (later TGF-beta, Finger 2021).
- Amplitude of the rescued cells and period are measured only on rhythmic cells (selection); single-cell phase coherence is not quantified; coupled phase information (phase pulling) is not tested with phase-shifted partners.

## Link to us
- **Intro claim(s):** supports cell autonomy of period plus a permissive, non-rhythmic paracrine requirement for rhythmicity; argues against rhythm-transmitting coupling between fibroblasts (contrast finger2021intercellular, where phase pulling up to 2 h was seen in U-2 OS).
- **Relation to clocksync:** background. For the question: each fibroblast is autonomous in period; neighbours enable (not pace) rhythm; no dominant population. Time scale days. Not related to brain death directly; the permissive paracrine requirement is an intra-tissue cue that donor death would not remove for hours.

## Ideas
- Compare to the finding that stromal-vascular cells share most state in atlases: a permissive tissue context (density, secreted factors) would raise amplitude without synchronising phases.
