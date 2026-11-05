# droin2019low
**Ref:** Droin, Paquet, Naef, Nat Phys 2019 (15:1086-1094). DOI 10.1038/s41567-019-0598-1. Type: method + biology. Lab: Naef (EPFL).
**Source:** full text

## Question and design
Fit a low-dimensional stochastic model of clock and cell cycle from single-cell time-lapse traces: NIH3T3 Rev-erbα-YFP (n = 2,504 traces reanalysed from Bieler 2014; coupling fit on dividing traces n = 154 / 271 / 302 at 34 / 37 / 40 C) and human U2OS-Dual (n = 551 cells). Plus population Bmal1-Luc luminescence under temperature cycles and a reanalysis of public mRNA data of 12 adult mouse tissues. **Sampling:** live imaging, 24-72 h; tissue data are public timed series.

## Method (if any)
Hidden Markov model on the signal gives circadian phase posteriors; the cell-cycle phase is interpolated between divisions; the coupling function F(theta, phi) is optimised. Model is one-way: the clock-to-cell-cycle term is set to zero because it was negligible in the Bieler 2014 data ("the influence of the clock on the cell cycle was negligible in NIH3T3 cells"). Outputs: attractors, Arnold tongues, noise effects.

## Key results
- 1:1 locking over a range of cell-cycle periods. "we find 1:1 phase-locking over a range of cell-cycle times varying from 19 h to" 27 h; circadian phase at division 0.7 to 0.9 x 2 pi. In the 1:1 state "the endogenous circadian period of 24 h is shortened by" about 2 h (acceleration after cytokinesis through G1).
- Predicted shifts validated with 12 perturbation experiments, no extra free parameters.
- Arnold tongues 1:2, 1:1, 2:1, 3:1; data show 1:1 and 1:2; "the deterministically predicted 2:1 state is replaced in the data and stochastic system by 1:1-like orbits"; "noise extended the range of the 1:1 tongue, but asymmetrically".
- Human U2OS: divisions 4.96 +/- 2.6 h before a circadian peak. "Thus, similarly to mouse NIH3T3 cells, the coupling directionality is predominantly from the cell cycle to the circadian clock."
- Dividing cells are not entrained by temperature cycles (Fig. 6 title: "Temperature cycles do not entrain circadian oscillators in dividing cells").
- In vivo: "which revealed that clock phases span 1.5 h between the earliest and latest tissues"; Mki67 and Myc correlate with phase offsets "the levels of known markers of cell proliferation such as Mki67 or Myc were strongly correlated with the phase offsets"; amplitudes did not.

## Not addressed / limits
- Tissue result is correlational, 12 tissues, mean expression; confounded by tissue type and cell composition.
- Model assumes no clock-to-cell-cycle term; deterministic 2:1, 3:1 states are rarely observed.
- Cancer statement is a hypothesis: "the damped rhythms may also reflect circadian desynchrony of otherwise functional circadian oscillators", following "in the presence of non-coherent cell-cycle progression".

## Link to us
- **Intro claim(s):** refs_cancer_phase_inference (low-dimensional Arnold-tongue description, predicts locking and cell-cycle slowdown) needs rewording: the model predicts circadian phase shifts and lock states as a function of cell-cycle period; the clock is accelerated (about 2 h in 1:1), the cell cycle is the input, not slowed. Naef-lab attribution correct.
- **Relation to clocksync:** supports (hypothesis that tumour clock damping is desynchrony across cells with non-coherent cell cycles; tissue phase offsets track proliferation).

## Ideas
- Proliferation versus organ clock phase offsets across tissues of one donor (GTEx, A36-style): not found in registered/results by grep (A97/A126/A127/A166 test coherence). Registrable: arg z per tissue against Mki67/MYC tissue means, with a specificity set of non-proliferation genes.
- Missing reference: Zhang 2014 (CircaDB tissue series; in notes as zhang2014circadian), used for the 12-tissue analysis.
