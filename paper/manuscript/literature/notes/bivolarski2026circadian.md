# bivolarski2026circadian
**Ref:** Bivolarski I, Dimitrov M, Dechev H, Cureus 18(7):e113586, 2026. DOI 10.7759/cureus.113586. Type: review (narrative). Lab: none (Integrated Oncology Centre Burgas; clinical authors).
**Source:** full text

## Question and design
Narrative review of circadian regulation in prostate cancer (clock genes, AR signalling, RNA-seq, alternative splicing, chronotherapy). Search in PubMed, Scopus, Web of Science and Google Scholar through July 2026; no systematic selection, no risk-of-bias assessment, no meta-analysis. No new data. **Sampling:** not applicable (review); it states that prostate tumour RNA-seq data are single-time-point.

## Method (if any)
None. Mentions CYCLOPS, MetaCycle, JTK_CYCLE and blood phase predictors (Hughey; refs 25-27) as tools; no benchmark.

## Key results
- Most clinical tumour transcriptomes lack collection time. "Most clinical transcriptomic datasets were generated without recording the time of tissue collection." (Methodological landscape)
- Collection follows hospital logistics, not circadian considerations. "In routine clinical practice, tissue samples are collected according to institutional schedules rather than circadian considerations." (Methodological landscape)
- Cross-sectional tumour vs normal differences cannot be read as altered oscillation. "Differences between tumor and normal tissue may reflect disease biology, sampling time, or an interaction between the two." (Deciphering transcriptome oscillations)
- Prostate RNA-seq cannot establish 24-h oscillation. "In prostate tumors, most available RNA-Seq datasets were generated from samples obtained at a single recorded or unrecorded time point." (same section)
- Phase caveat: "Clock time is also not equivalent to internal biological time because chronotype, sleep patterns, light exposure, medication, and systemic illness may shift circadian phase." (Methodological landscape)
- Inferred phase is an estimate only. "Computationally inferred phase should therefore be regarded as an estimate rather than a substitute for prospectively collected time-resolved data." (Chronobioinformatic approaches)
- CYCLOPS is described as ordering samples without times, with performance depending on "sampling density, cohort size, signal strength, and the preservation of rhythmicity within the examined tissue". (Introduction)
- Disease model (qualitative, no data): progression may involve "selective remodeling of individual clock components, altered phase relationships, and loss of coordinated rhythmic output". (Alterations in clock gene networks)

## Not addressed / limits
- Secondary source: no measured distribution of collection times for prostate or any other cancer; the "institutional schedules" statement is asserted, uncited in that paragraph.
- Published 29 July 2026 in Cureus (low-tier journal); written largely from preclinical and cross-sectional work.
- No mention of host-referenced or multi-compartment inference; does not discuss tumour vs normal phase offsets within a patient.
- Does not treat the Hughey/Anafi-type methods quantitatively (training population, coverage).

## Link to us
- **Intro claim(s):** Technical challenges (hospital hours): tumour tissue is collected on institutional schedules, time not recorded. Our refs_cancer_phase_inference.md row (secondary source, weak, not a measured time distribution) is accurate: SUPPORTED, correctly hedged. Also supports the statement that disease-vs-sampling-time cannot be separated in single-time-point data.
- **Relation to clocksync:** background; the phase altered or time of sampling (our paraphrase) point is exactly what donor-matched host referencing removes.

## Ideas
- Cite as secondary support only; pair with lamia2024telling and gonzalezaponte2025circadian for primary evidence on biopsy timing.
- Prostate cohorts already scored in clocksync (GSE116918, 248 patients): done, A182 (registered_experiments.md); no new analysis proposed.
