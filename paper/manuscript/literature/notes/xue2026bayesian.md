# xue2026bayesian
**Ref:** Xue X, Wang H, Zhu Z, Yu G, Seney ML, Zhang L, Tseng GC, bioRxiv 2026 (10.64898/2026.01.21.700905), submitted to Annals of Applied Statistics. Type: method. Lab: George Tseng (Pittsburgh).
**Source:** full text (preprint)

## Question and design
BayCT / BayCT-M: jointly detect circadian genes and estimate each subject's molecular circadian time (MCT) given recorded time of death (TOD). Real data: human post-mortem NAc, caudate, putamen of the same N = 59 subjects (Ketchesin 2021), bulk RNA-seq, 14,540 transcripts; mouse 12-tissue atlas (Zhang 2014, every 2 h for 2 days) with synthetic TOD error. **Sampling:** TOD known for every subject (ZT scale, about -6 to 18 h in brain); mouse full 48 h.

## Method
- Unit: subject (donor), with region-specific rhythmic parameters. Semi-supervised: TOD is an input; the unknown is a deviation delta t_i = MCT - TOD with a von Mises prior (kappa = 6, about 1.5 h SD), so it corrects TOD rather than inferring time de novo. Needs TOD for all subjects, hence a coverage of the day determined by the data.
- Gene acrophases learned: cosinor with spike-and-slab rhythmicity (reversible-jump MCMC), Gaussian noise; in BayCT-M phases phi_gj have a shared mean theta_g across regions, amplitudes and variances region-specific. delta t_i updated only with pre-selected core clock genes (11: NPAS2, PER2, CRY1, CIART, PER1, NR1D1, NR1D2, TEF, CRY2, PER3, DBP).
- Within-person tissues: yes, in BayCT-M the subject's delta t_i is shared across regions, but region phase differences are modelled as gene-level offsets from theta_g, not as donor-level desynchrony.

## Key results
- Simulation (30 samples, 1000 genes): BayCT beats ZeitZeiger, TimeSignature, CHIRAL, CYCLOPS, Cyclum when TOD is wrong; at A = 0.8 with 40% of samples off by 4 h, AUC 0.909 vs 0.869 (BayC) and 0.852 (cosinor).
- Multi-region gain: SD ratio BayCT-M/BayCT approaches 1/sqrt(3), consistent with pooling 3 regions.
- Mouse with synthetic delta t: "the 95% credible intervals (CIs) covered 94.4% of the artificial".
- Human brain: all BayCT-M delta t have 90% CI width below 2.4 h; 54/59 agree with TimeSignature within the stated tolerance (text says "44 hours", evidently a typo); examples delta t = 3.54 h [2.87, 4.19] (subject 43), 6.68 h [5.84, 7.51] (subject 103).

## Not addressed / limits
- No accuracy against an independent ground truth in human data; the 3 regions are same-tissue-type (striatum), so shared delta t is expected; no desynchrony between organs.
- Needs recorded TOD (absent in many biopsy and single-cell cohorts); a recorded hospital-hours time could serve as the prior, untested here; delta t prior (1.5 h) limits larger shifts.
- 4.45 h for 3000 MCMC iterations on 59 subjects; R package only.
- No claim of window coverage; our refs line is right on this.

## Link to us
- **Intro claim(s):** supports that recorded TOD differs from molecular time (technical challenges) and that pooling regions of the same donor sharpens the donor clock state: the same logic as the host reference in clocksync.
- **Relation to clocksync:** closest published precedent for a donor-level shared clock state across compartments (competitor, compare). Differs: it estimates one shared delta t, whereas clocksync estimates per-compartment coupling, |z| and offset, and does not need TOD.

## Ideas
- Run BayCT-M on GTEx donors with TOD (several tissues) and compare its delta t to clocksync's host state; check whether residual per-tissue offsets match arg z.
- Ask the authors for the NAc/caudate/putamen data (Ketchesin 2021, GEO) as a within-person test of coupling across compartments.
