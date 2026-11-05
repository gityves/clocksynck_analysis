# ballesta2017systems
**Ref:** Ballesta A, Innominato PF, Dallmann R, Rand DA, Lévi FA, Pharmacol Rev 69:161-199 (2017). DOI 10.1124/pr.116.013441. Type: review. Lab: Lévi / Rand (Warwick).
**Source:** full text

## Question and design
Review of systems chronotherapeutics: multiscale models (clock, cell cycle, cell populations, PBPK-PD) to design and personalise circadian-timed drug schedules, with irinotecan, oxaliplatin and 5-FU as examples; also rheumatology, cardiovascular, metabolic disease and circadian disruption. No new data.
**Sampling:** not applicable. Discusses patient circadian monitoring: wrist actigraphy, skin temperature, cortisol, PBMC clock-gene mRNA over 24 h.

## Key results
- Phase III trials of chronotherapy vs non-timed schedules gave "up to fivefold better tolerability and a near doubling in efficacy" (Section III.B, citing Dallmann 2016); meta-analyses suggest survival benefit in men (Giacchetti 2012).
- Many randomised morning vs evening comparisons found no difference; the authors list possible reasons: missed optimal time, dose, or interpatient differences masking the timing effect.
- Between-patient variability: pooled actigraphy of 436 metastatic colorectal patients showed activity maxima spread over a 10 h span; acrophase of skin temperature differed by up to 12 h among 24 gastrointestinal cancer patients. "Wrist actigraphy has further uncovered large intersubject differences in circadian phase and amplitude in a pooled analysis involving 436 patients with metastatic colorectal cancer." (II.C)
- Circadian disruption: altered rest-activity measured by actigraphy "was associated with worse prognosis in patients with mCRC" (VII.C); loss of diurnal cortisol predicts early mortality in metastatic breast cancer.
- Domomedicine pilot (InCASA): 31 metastatic cancer patients transmitted weight, symptoms and actigraphy daily for on average 58 days.
- Model-based personalisation: patient-specific PK-PD models fed by continuous circadian recordings; no patient yet treated with a model-derived schedule in the text.

## Not addressed / limits
- Models and mouse-to-human scaling remain unvalidated prospectively; no clinical trial shows that model-personalised timing beats standard chronotherapy.
- Human molecular clock state is measured only in blood (PBMC, healthy men over 24 h) or skin; tissue clocks of patients not accessible, which is the stated bottleneck ("Longitudinal measurement over several days of molecular biomarkers is currently not available in the clinics", Fig. 2 legend).
- Clock-cell-cycle coupling is modelled in cell lines, not in patients.

## Link to us
- **Intro claim(s):** models to design circadian-timed schedules ("verified"). Supported. Also supplies the clinical need: patient-specific internal timing is unobserved.
- **Relation to clocksync:** background; states the missing measurement (tissue clock state in patients) that single-biopsy host-referenced readout targets. Their Fig. 2 uses PBMC Per2 sampled over 24 h, the opposite of our untimed, hospital-hours setting.

## Ideas
- Inter-patient variability of up to 10-12 h (actigraphy, skin temperature) gives a prior for the spread of the shared host state across donors; compare with the GTEx spread d* in our results (not an existing comparison in our docs).
- Cite Mormont 2000 and Innominato 2014 (pooled actigraphy, n = 436) for rhythm-survival link.
