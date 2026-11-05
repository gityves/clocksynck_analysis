# levi2024circadian
**Ref:** Lévi FA, Okyar A, Hadadi E, Innominato PF, Ballesta A, Annu Rev Pharmacol Toxicol 64:89-114 (2024). DOI 10.1146/annurev-pharmtox-051920-095416. Type: review. Lab: Lévi / Ballesta.
**Source:** full text

## Question and design
Review: circadian regulation of ADMET, chronotherapy of immunotherapy, chemotherapy, cardiovascular and other drugs, sex-specific dosing, personalisation by circadian metrics. No new data.
**Sampling:** not applicable. Its Figure 2 reuses GTEx liver and small-intestine donors with CHIRAL-estimated internal phases (men vs women) and trial-derived least-toxic timings.

## Key results
- ICIs: "Seven independent retrospective reports have revealed significant and consistent improvements in progression-free and/or overall survival through the administration of a larger proportion of ICI infusions in the morning" in 1,019 patients (melanoma, NSCLC, renal, urothelial); median OS up to fourfold different (Section 4.1.1). All retrospective.
- Chemotherapy: systematic review of 18 randomised trials, 2,547 patients: efficacy improved in 3 (17%), toxicity reduced in 11 (61%), increased in 1; 14 colorectal chronotherapy trials, nearly 1,700 patients; peak 5-FU-leucovorin near 4:00 AM and oxaliplatin near 4:00 PM.
- Sex: "The overall survival of metastatic colorectal cancer patients receiving a chronomodulated oxaliplatin and 5-FU-leucovorin combination was significantly prolonged in male patients and significantly reduced in female patients" (3 randomised trials, 842 patients). Least-toxic timing 4-6 h later in women (113 patients, oxaliplatin-5-FU-leucovorin; 193 patients with irinotecan). Rituximab-CHOP (n = 210): afternoon better than morning in women only.
- Platinum drugs: toxicity reduced up to threefold in rodents at ZT15-19; tolerability better afternoon/evening in several cancers.
- Personalisation: DLMO from 152 studies, 4,397 healthy subjects shows large between-subject variation; rest-activity peaks and skin-temperature acrophases differ by up to 12 h between subjects. Concept of chronofit (robust circadian function) vs non-chronofit patients; the abstract states that molecular clock function is predictable "from single-timepoint biopsies" (ref 13, a ToD-type predictor).

## Not addressed / limits
- Sex-specific timing differences come from retrospective analyses of prospectively collected data, not randomised sex-by-time designs; no CIs given in the review. ICI evidence has no randomised arm and no circadian marker.
- No test of whether tissue-level clock disruption modifies timing benefit; chronofit stays a proposal.

## Link to us
- **Intro claim(s):** sex-specific, personalised chronotherapy and the need for per-patient circadian state ("not confirmed" in refs_*.md, now verified for the sex and personalisation statements). Supported; note ICI evidence is retrospective.
- **Relation to clocksync:** supports the problem statement (need to identify patient internal state from a biopsy); names untimed single-biopsy predictors as a future tool, i.e. our setting. Competing idea: ToD-type predictors need a reference model and timed training data.

## Ideas
- Chronofit as a clocksync endpoint: compartment coupling c / |z| of host-referenced readouts could define chronofit without time; check against GTEx donor traits already done (T2D/obesity A164, age A94, sex phase A99).
- Cite Dallmann 2016 (Trends Mol Med, dosing time) and Giacchetti 2012 meta-analysis of the sex interaction.
