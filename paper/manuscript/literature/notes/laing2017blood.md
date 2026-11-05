# laing2017blood
**Ref:** Laing EE et al., eLife 2017;6:e20214. DOI 10.7554/eLife.20214. Type: method (PLSR blood predictor of melatonin phase) + data. Lab: Dijk / Moller-Levet, Surrey Sleep Research Centre.
**Source:** full text

## Question and design
Human whole-blood microarray transcriptomes; 678 mRNA profiles from 53 participants in two Surrey experiments (mistimed sleep; sleep restriction/deprivation), four sleep-timing conditions, split into a training set (329 samples, 26 participants) and a validation set (349 samples). Each sample is assigned a melatonin phase (DLMO-referenced).
**Sampling:** timed: "blood samples (every 3 or 4 hr) covering at least one circadian cycle" in each experiment; samples covered the entire melatonin cycle, in dim light, constant lab conditions, young healthy volunteers.

## Method (if any)
Supervised partial least squares regression of melatonin phase (cos and sin) on about 26 K features; final models use 100 biomarkers (glucocorticoid signalling, immune function). One-sample and two-sample (differential, 12 h apart) variants. Compared to molecular timetable (73 features) and ZeitZeiger. Leave-one-participant-out CV. Tissue: whole blood on one microarray platform (GPL15331); models usable on the same platform or features. Needs timed training data across the cycle; a test subject needs one sample (R2 0.74) or two samples 12 h apart (R2 0.90).
On other tissues: "one can imagine ... the development of specific PLSR-based circadian phase prediction models for the liver and/or other specific tissues and organs", requiring "the identification of specific phase markers for the tissue/organ of clinical interest."

## Key results
- One-sample PLSR: "Fifty-four percent of samples fell within 2 hr of the true circadian phase and the R2 of 0.74" (validation).
- Two-sample differential: "the two-sample differential model is able to correctly predict 82% of samples within 2 hr, with an R2 of 0.9"; mean error 20 min, SD 1 h 41 min.
- Robustness: one-sample error did not differ across the four conditions; the two-sample model was significantly worse under "sleep out of phase with melatonin".
- Rhythmicity of the transcriptome is "most prominent in the 'sleep in phase with melatonin' condition and smallest in the 'sleep out of phase with melatonin' condition"; only 18% and 12% of model genes stayed rhythmic when sleep was out of phase.

## Not addressed / limits
- Training and validation limited: "all samples were collected in darkness and dim light and from young healthy individuals"; no patients.
- Blood only, so no evidence on how other organs shift; the model reads the SCN-driven melatonin phase through blood, and sleep-wake desynchrony is handled only as altered blood rhythms.
- Microarray platform bound; no external dataset.
- Median absolute error in h not the headline metric (mean error, R2 and fraction within 2 h are).

## Link to us
- **Intro claim(s):** supervised single-sample blood phase predictor from timed data. Our refs line (100 genes, constant routine / sleep deprivation, single sample) is supported but incomplete: the conditions are mistimed sleep and sleep deprivation, and the model works with one or two samples.
- **Relation to clocksync:** background / data. The same Surrey cohorts (GSE39445, GSE48113) have melatonin phase plus blood phase and could validate blood-versus-melatonin offsets.

## Ideas
- Use the sleep-out-of-phase condition as a positive control for internal misalignment between blood and SCN (melatonin): does a blood clock-gene read-out lose amplitude or shift relative to melatonin?
- The paper notes that models trained only in normal conditions degrade under desynchrony, so training must include misalignment: relevant for our gate (a failed gate is uninformative).
- Missing: Archer 2014 PNAS (mistimed sleep disrupts blood transcriptome).
- Done: forced desynchrony in human blood (GSE48113, registered_experiments.md / results_preliminary.md).
