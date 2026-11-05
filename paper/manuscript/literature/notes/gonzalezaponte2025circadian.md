# gonzalezaponte2025circadian
**Ref:** Gonzalez-Aponte MF et al., J Neuro-Oncol 176:36, 2026 (online 22 Oct 2025). DOI 10.1007/s11060-025-05242-3. Type: clinical + biology (in vitro, model). Lab: Herzog / Walch (Washington University in St. Louis).
**Source:** full text (Supplementary figures S1-S12 not in the text)

## Question and design
Does MGMT promoter methylation, the GBM biomarker for temozolomide (TMZ) benefit, vary with time of day, so that the biopsy hour biases the diagnosis? (i) Human GBM lines (LN229, primary B165, murine Nf1-/-DNp53, GL261) synchronised in vitro, sampled at CT4 and CT16 (qMSP, n = 3) or every 4 h for 24 h (immunostaining of MGMT, BMAL1, n = 6); BMAL1 knockdown. (ii) Retrospective chart review at Barnes-Jewish Hospital: newly diagnosed IDH-wildtype GBM, Jan 2020 - Oct 2024, 302 patients. (iii) ODE model of TMZ and MGMT. **Sampling:** biopsies timed by end of anaesthesia, from 10:00 to 18:00 only (daytime surgery window); in vitro full 24 h at 4-h steps.

## Method (if any)
JTK_CYCLE on the hourly-binned probability of an MGMT-methylated call (p < 0.05); Rayleigh test on time-of-day distributions per class; analysis restricted to hours with 6-15 biopsies. "Time of surgery was determined by the recorded end time of anesthesia to approximate when the tissue samples were embedded for analysis."

## Key results
- n = 302 patients. "We evaluated MGMT promoter status from patient biopsies collected at different times of day, from 10am to 6pm, binned every two hours." (Results)
- Probability of a methylated call peaks around noon, minimum (within window) near 16:00, JTK p < 0.05, "varying by at least 20% over the course of the day". (Results, Fig. 1D)
- "the probability of detecting MGMT promoter methylation significantly decreases from around 50% to 30% when patients were biopsied later in the day compared to in the morning". (Discussion)
- Rayleigh: methylated tumours peak 15:00 (arrow length 0.839), unmethylated 17:00 (0.842), both p < 0.05 (Fig. 1C legend).
- In vitro: lower MGMT promoter methylation at CT4 than CT16 in all lines; none in BMAL1 KD. MGMT and BMAL1 protein peak at CT4 (JTK p < 0.05).
- 100 µM TMZ gave "an average of 60% cell death when delivered at CT4, compared to 30% at CT16"; no rhythm at 10 or 1000 µM; abolished by BMAL1 KD.
- Model: best dosing CT2-CT8, never at the MGMT minimum; 73% vs 49% death at CT4-6 vs CT16-20.

## Not addressed / limits
- Window 10:00-18:00 only: "we do not know what happens to MGMT promoter methylation from around 6pm until 10am" (Discussion). Shows only an 8-h span; time of surgery is clock time, not internal phase.
- Single institution, retrospective; no control for fixation time, tumour location, anaesthesia duration, age, sex (the authors list these as to be controlled in future trials).
- The effect is on a binary diagnostic call, not on a clock-gene expression or phase readout; no tissue clock-gene data, no test that tumours oscillate in vivo. Patients are not sampled repeatedly.
- Rayleigh arrows per class largely reflect when surgeries cluster; only the JTK on the probability is informative.
- Cell-line rhythms need synchronised cultures; relevance to heterogeneous patient tumours is inferred.

## Link to us
- **Intro claim(s):** surgical sampling time varies and matters for tumour biomarkers. Our refs row (peak around noon, minimum near 16:00, 302 patients 2020-2024) is SUPPORTED, but should add that the span is 10:00-18:00 only, which also supports "hospital hours". It is not evidence that surgery times are spread over 24 h.
- **Relation to clocksync:** background (hospital-hours sampling is a daytime window of ~8 h; any phase-ordering method would see only a third of the cycle). Data we could not use (no expression).

## Ideas
- Quote the 10:00-18:00 window as the empirical hospital-hours range when motivating donor-referenced inference; combine with lamia2024telling.
- GBM-specific check (tumour vs normal brain compartment offsets) needs a matched dataset; tumour phase offset analyses done across cancers (A182 and related, see registered_experiments.md); no GBM-specific MGMT analysis exists (no match for MGMT/GBM in METHOD.md or registered_experiments.md).
- Related references in its list: Gonzalez-Aponte 2024 (MGMT circadian regulation), Damato 2021/2022 (TMZ timing and survival), Slat 2017 (GBM clock): check that they are in our bib.
