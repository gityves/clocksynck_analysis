# yoo2004periodluciferase
**Ref:** Yoo et al., PNAS 2004. DOI 10.1073/pnas.0308709101. Type: biology (reporter mouse). Lab: Takahashi.
**Source:** full text

## Question and design
Are peripheral clocks self-sustained or damped, and does the SCN drive them? New mPer2Luc knock-in mouse (PER2::LUC
fusion at the endogenous locus; behaviour normal, n = 7/10/7 per genotype for period). Tissue explants (SCN, retrochiasmatic
area, pituitary, cornea, liver, lung, kidney, tail) in static culture, luminescence by PMT in 1 min bins for >7 days
(up to >20 days). SCN-lesioned animals (n = 11) versus intact controls in LD (n = 8) and DD (n = 8); explants of 3-5 weeks
after lesion. **Sampling:** timed; explants taken around lights-off (about ZT12) and phase is read from the first
full peak 12-36 h in culture; recordings continuous over days. One explant per tissue and animal.

## Method (if any)
Phase = highest point of the 2 h smoothed trace; period = mean peak-to-peak; ANOVA on phase and period (tissue,
lesion). The explant phase is taken to reflect the phase of the tissue in vivo at removal.

## Key results
- Peripheral explants oscillate for >20 cycles: "peripheral tissues are capable of self-sustained circadian oscillations for >20 cycles in isolation." (Abstract)
- Tissue-specific phase: SCN peak at circadian time 12; all peripheral tissues delayed (ANOVA F(6,21) = 7.19, P = 0.000285); SCN-liver delay about 4 h (Fig. 3B).
- Tissue-specific period (F(6,21) = 3.02, P = 0.027): cornea 22.2 h shortest, kidney 24.8 h longest, SCN 23.5 h (Fig. 3C). "the period of the SCN (23.5 h) was the same as that seen for the free-running period of locomotor behavior".
- Phase is not explained by period (R = 0.002, P = 0.99), so entrainment signals or their responsiveness differ by tissue.
- SCN lesion does not abolish peripheral rhythms but destroys inter-tissue phase coordination (tissue effect F(4,115) = 6.21, P = 0.000147; lesion effect F(2,115) = 22.95).
  "SCN lesions disrupted phase synchrony even from tissue to tissue within individual animals, which is a form of 'internal desynchronization.'" (Results)
- Periods stay tissue-specific after lesion (F(4,113) = 10.96); cornea stays coordinated (light-entrainable via the eye).
- Averaging desynchronised animals lowers the ensemble amplitude: "population sampling experiments of SCN-lesioned mice would be expected to have reduced amplitude as a consequence of desynchrony rather than loss of rhythmicity at the level of individual animals/tissues."

## Not addressed / limits
- Mouse, ex vivo reporter of PER2 protein, not endogenous RNA in vivo; phase read after explanting.
- Lesion experiment: phases at 12/16/32 days in DD only; no direct measure of inter-tissue coupling or of its
  mechanism; "organ-specific synchronizers" are inferred, not identified.
- No human data; no statement on how large between-tissue offsets are under normal life with feeding and activity.

## Link to us
- **Intro claim(s):** peripheral tissues carry autonomous, tissue-specific clocks; losing the central synchroniser
  gives inter-organ desynchrony (Internal desynchrony). refs summary (tissue-specific damping) is only partly right:
  the paper's headline is persistence (little intrinsic damping), tissue-specific period and phase, and desynchrony
  after lesion. Suggest rewording to tissue-specific period and phase.
- **Relation to clocksync:** supports. Gives the earliest explicit definition and demonstration of internal
  desynchronization between organs of one animal, and shows that pooling desynchronised individuals mimics rhythm loss
  (relevant to cohort-level phase-ordering methods and to our within-donor referencing).

## Ideas
- Their tissue phase map (SCN 12, liver about +4 h) is an external check of arg z offsets between compartments
  where a donor-level host exists (GTEx tissues).
