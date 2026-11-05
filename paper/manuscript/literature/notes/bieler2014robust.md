# bieler2014robust
**Ref:** Bieler, Cannavo, Gustafson, Gobet, Gatfield, Naef, Mol Syst Biol 2014 (10:739). DOI 10.15252/msb.20145218. Type: biology + model. Lab: Naef (EPFL).
**Source:** full text

## Question and design
NIH3T3 mouse fibroblasts carrying a Rev-erbα-YFP circadian reporter (division visible as a 30-60 min dip from nuclear envelope breakdown), unstimulated, free-running, 72 h at 30 min intervals, serum 2-13% pooled. Over 10 000 single-cell traces: "totaling 20,000 circadian peaks and 13,000 cell divisions". Perturbations: temperature 34/37/40 C (n = 1,139 / 4,207 / 1,374 traces), CDK1/CDK2 inhibitors, Cry2 shRNA, Longdaysin (CK1), dexamethasone. **Sampling:** live imaging, continuous, 72 h (not a tissue sampling design).

## Method (if any)
Stochastic model of two coupled noisy phase oscillators (circadian phase, cell-cycle phase) fitted by likelihood to traces; yields coupling functions in each direction. Model-free check: Granger causality of nuclear size (cell-cycle proxy) on the reporter and the reverse.

## Key results
- 1:1 mode locking. "both oscillators tick in a 1:1 mode-locked state, with cell divisions occurring tightly 5 h before the peak in circadian Rev-Erba-YFP reporter expression" (division-to-peak interval 5 +/- 2 h, about CT1).
- Robust to fluctuations: "were highly correlated, although the individual circadian and cell cycle intervals varied by more than 30%" (R2 0.52, n = 1,230).
- Direction is predominantly cell cycle to clock. "our data combined with stochastic modeling unambiguously show that the reverse coupling is predominant in NIH3T3 cells" (reverse = cell cycle acts on the clock).
- Dividing cells have shorter circadian intervals: 21.95 +/- 3.8 h (n = 1,926) versus 23.7 +/- 3.1 h without division (n = 2,748), P < 1e-16.
- Granger: "showed evidence (P < 0.001, Granger-Wald test) for a causal influence of cell cycle progression on the circadian signal" in up to 60% of cells versus under 20% for the reverse (55% versus 12% when stronger in one direction).
- Robust: Cry2 knockdown (period 26.3 h) "did not affect cell cycle duration", divisions shift earlier in phase as predicted; Longdaysin lengthens circadian period from 24 to nearly 32 h with coupling unchanged.
- Gating not excluded but subordinate. "Although this finding did not completely exclude that a circadian gating of the cell cycle occurred as well"

## Not addressed / limits
- One immortalised mouse line, one reporter (Rev-erbα), in vitro; "Conceivably, this property may also depend on the model organism or cell types."
- Free-running only; no entrainment by systemic cues. Direction of coupling is a model inference. CRY dilution as mechanism is speculation.
- No tissue or in vivo data; no disease.

## Link to us
- **Intro claim(s):** refs_cancer_phase_inference says bidirectional coupling gives 1:1 phase locking: OVERSTATED. The paper's conclusion is that coupling is mainly one way, cell cycle to clock; circadian gating is weak or absent in NIH3T3. Locking arises mainly because divisions advance the clock. Abstract implication stated: "These findings have implications for circadian function in proliferative tissues, including epidermis, immune cells, and cancer."
- **Relation to clocksync:** background (why a proliferating compartment may have a shifted clock phase and a lower coupling to the host).

## Ideas
- Prediction for tissue data: proliferating compartments should show clock-phase offsets (division-driven acceleration, Droin: 24 h period shortened by about 2 h) and not only lower coherence; A97/A127/A126 test coherence versus proliferation, offsets versus proliferation were not found by grep. Test arg z versus G2M score in sc atlases (cycling versus quiescent cells of the same donor and tissue).
