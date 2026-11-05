# wu2018populationlevel
**Ref:** Wu G et al., PNAS 2018;115(48):12313-12318. DOI 10.1073/pnas.1809442115. Type: method + data (population-scale unordered skin, CYCLOPS ordering, ZeitZeiger biomarker). Lab: Hogenesch (Cincinnati Children's) with Procter & Gamble.
**Source:** full text

## Question and design
Human forearm skin, laser-capture-microdissected epidermis (and dermis; the epidermis is analysed), Affymetrix microarray, GSE112660. Hybrid design: 20 healthy Caucasian men (21-49 y) housed in a sleep facility, punch biopsies every 6 h over 24 h ("6 AM, 12 PM, 6 PM, and 12 AM", four per subject, 79 time-stamped samples including 19 with full series), plus 219 females (152 Caucasian, 67 African-American; 20-74 y) sampled once, "irrespective of wall time" and without sampling time; 298 samples from 238 donors in total.
**Sampling:** time-stamped, full 24 h at 6 h resolution for the 20 men; unordered (time unknown to the analysis) for the 219 donors; saliva melatonin and cortisol every 3 h confirm entrainment in the men.

## Method (if any)
CYCLOPS (unsupervised circular ordering on a seed set of clock and clock-associated genes, with Oscope to filter eigengenes) assigns a phase to every one of the 298 samples without wall time; modified cosinor then finds circadian genes. ZeitZeiger is then trained on CYCLOPS-assigned phases (219 unordered samples) plus 43 time-stamped samples from 11 subjects, and tested on 36 time-stamped samples of nine held-out subjects. Needs a population large enough to cover the cycle (the unordered samples are assumed to span all phases); one tissue (epidermis); a 29-gene biomarker set.
Authors on other tissues: "We believe that the hybrid design could be advantageous for circadian studies in many different human tissues."

## Key results
- Clock intact at the population level in unordered data: "We found that clock gene correlations are strong and evident in unordered data from human epidermis." Human clock-gene coexpression resembles 12 mouse tissues (Mantel test P = 3e-06) but not blood (P = 0.897; whole blood n = 201 time series).
- CYCLOPS order vs known time: Fisher circular correlation 0.955 between mean CYCLOPS phase and sampling time (20 subjects); human vs mouse clock-gene phase order 0.715; 188 circadian genes at P < 0.01, rAMP > 0.1; 110 rhythmic genes in the 19-subject series (peaks 8-9 AM and 8-9 PM).
- Biomarker accuracy: single-sample phase "within 3 h"; "the mean absolute error was 2.95 h for the nine subjects", "we were able to accurately assign 30 of 36 samples to their correct circadian phase"; 6 PM and 12 AM samples hard to separate in six subjects.
- Interindividual shifts: subjects 116 and 119 "show delayed clock phases of genes, salivary melatonin, and cortisol", consistent between skin and saliva.
- Comparison with Wittenbrink: BodyTime reports "accurate assignment of DLMO phase from a single sample to within 3 h" in 40 individuals; skin tested over the full cycle in 239 subjects.

## Not addressed / limits
- Not a misalignment or disease study: healthy volunteers only; the authors say "the consistency of results should be confirmed across a range of disease states and pathologies."
- Test subjects were chosen because "they had CYCLOPS-predicted phases surrounding the mean phase at each time point", which flatters the 2.95 h error; the 19 other subjects were used in training.
- Error is against clock time, with training phases derived by CYCLOPS, itself with limited accuracy; the unordered samples are all females and the time series all males (sex, ancestry confounded with design).
- Epidermis and dermis are compartments of the same biopsy, yet the paper does not compare their phases within a person; no inter-organ offsets.
- CYCLOPS requires the unordered cohort to span 24 h, so hospital-hours-only cohorts would not satisfy it (the authors do not state the hour distribution of the 219 single biopsies in the main text).

## Link to us
- **Intro claim(s):** human skin rhythms ordered at the population level; skin is a better clock source than blood. Our refs line 28 (skin rhythms from biopsies ordered by time of sampling) is wrong: the 219 samples were not ordered by sampling time; they were reconstructed by CYCLOPS without wall time, and only 79 samples carry a time stamp. Better framing: population-level untimed reconstruction (hybrid design), the nearest published precursor to untimed inference.
- **Relation to clocksync:** competitor / background for untimed inference, but it still presumes phase coverage across 24 h in the cohort; clocksync needs none and uses within-donor references.

## Ideas
- GSE112660 (public): epidermis and dermis from the same biopsy plus time stamps for 79 samples, a ready test of clocksync's within-donor compartment offset (epidermis vs dermis) and of the gate.
- Use the 20-man time series to validate arg z against saliva melatonin phase (subjects 116 and 119 as phase-delayed positive controls).
- Compare CYCLOPS-assigned phases against clocksync host-referenced offsets on the 219 unordered samples for specificity.
- Reference check: Brown 2008 (cell-type clock strength in blood monocytes vs fibroblasts) cited as ref 24.
- Done: A37 (living human skin, epidermis vs dermis: untimed offset -1.6 vs timed -2.2 h; METHOD.md validation table).
