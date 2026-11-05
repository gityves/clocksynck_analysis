# vasilyev2026robust
**Ref:** Vasilyev V, Vlachou D, Giacchetti S, Bjarnason GA, Martino TA, Levi F, Dallmann R, Rand DA, bioRxiv 2026 (10.64898/2026.03.15.711880). Type: method application + clinical. Lab: David Rand / Robert Dallmann (Warwick).
**Source:** full text (preprint)

## Question and design
TimeTeller (ref [38] in the paper, trained on healthy time series) applied to 1286 early-stage breast cancer biopsies from six cohorts (REMAGUS n = 226, Mainz 200, Rotterdam 286, Transbig 198, Unt 125, Upp 251), Affymetrix microarrays, 10 clock-gene probes. Endpoints OS (n = 442), DMFS (809), RFS (776). **Sampling:** diagnostic biopsies, clock time recorded only for 108 REMAGUS samples; no 24 h coverage (working day).

## Method
- Unit: single sample (one tumour per patient, independent of other samples). Supervised: probabilistic model of the 10-gene vector (means and gene-gene correlations over time) learned from healthy timed data; output = likelihood curve L(t), TT time (argmax, "internal circadian phase"), ML and a dysfunction score Theta (0 to 1, sharpness, ML, secondary peaks).
- Gene acrophases fixed by the training clock. Does not need full-cycle coverage of test samples (one sample at a time). Noise model: not detailed in the main text (SI).
- Across compartments in a person: no (bulk tumour only).

## Key results
- 93.2% of tumours have TT time inside the working day (7 am to 8 pm): "for all but 88 (6.8%) of the 1,286 tumour samples, the estimated TT time fell within the working day window".
- Lasso day-night classifier: "only 24 of the 1286 samples are classified during nighttime".
- No link to external time (REMAGUS, n = 108): TT time vs sampling time slope -0.027, p = 0.80; interpreted as "temporal uncoupling: tumour clocks are running but have lost synchronization with the external light-dark cycle".
- Within-TT-time precision: SD about 1 h of TT time given the gene state.
- Survival (REMAGUS OS): nonlinear Cox, risk falls about 3-fold from Theta 0 to 0.35 (LRT nonlinearity p = 0.0009; Theta p = 0.003); Theta > 0.35 not significant; dichotomised HR = 0.48, log-rank p = 0.013. "good tumour clock function is antagonistic to survival".
- Theta barely tracks proliferation (r = 0.07).

## Not addressed / limits
- Breast cancer only: "Our study focused exclusively on breast cancer, and generalizability to other tumour types remains unknown." Our refs line ("across cancers") is OVERSTATED.
- Retrospective, 6 cohorts with different treatments; bulk tissue (tumour vs stroma/immune unresolved); microarray only (no RNA-seq training set).
- The near-universal working-day TT times and the absent correlation with sampling time (n = 108) could also reflect the training window or a collapsed signal; the authors add a lasso day/night check but no negative control (for example healthy hospital-hours tissue).

## Link to us
- **Intro claim(s):** supports that tumour clock function relates to outcome in breast cancer and that tumour phase is decoupled from clock time. Needs rewording: single cancer, non-monotone relation, good clock = worse survival.
- **Relation to clocksync:** competitor in spirit (clock quality score per sample); evidence that hospital-hours tumour samples carry a phase unrelated to wall time, which supports the premise that constant offsets cancel. Cannot give coupling between compartments of a donor.

## Ideas
- Compare clocksync |z| / coupling to Theta on breast cohorts with matched normal tissue (TCGA-BRCA, GTEx breast); check the nonlinear survival shape with a specificity companion.
- Ask whether TT times cluster at midday because of the training window: apply TimeTeller to GTEx hospital-hours subsets as negative control.
