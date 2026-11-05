# wittenbrink2018high
**Ref:** Wittenbrink N et al., J Clin Invest 2018;128(9):3826-3839. DOI 10.1172/JCI120874. Type: method (supervised phase predictor, BodyTime). Lab: Kramer (Charite Berlin).
**Source:** full text

## Question and design
Healthy human volunteers, CD14+ blood monocytes (MACS-sorted), RNA-Seq (discovery) then NanoString 48-plex (assay). Two cohorts.
BOTI (training): 12 men (10 passed QC), 22-30 y, 14 blood samples per subject. VALI (validation): 28 extreme/moderate larks and owls, 2 samples each.
**Sampling:** timed. BOTI: "14 blood samples were taken at regular intervals of 3 hours over a period of 40 hours in a constant routine protocol" (full 24 h and beyond, sleep deprivation, dim light). VALI: ambulatory, normal sleep/meals/light; morning sample about 2 h after habitual wake, second 6 h later (day only). Internal time = hours relative to DLMO (saliva melatonin).

## Method (if any)
Supervised ZeitZeiger (sparse PCA + maximum likelihood) on timed training data from one tissue (monocytes), constant-routine protocol, full-cycle coverage. Predicts internal time (DLMO-referenced) or external clock time from a single sample (1-sample) or the ratio of two samples 6 h apart (2-sample). Final predictors use 2 genes (NR1D2, PER2) or 12-13 genes (clock genes plus LGALS3, ELMO2, FKBP4, CRISPLD2...).
Applicable only to monocytes (sorted CD14+) measured on the NanoString panel; not to other tissues. Needs no 24 h coverage of the test donor (single sample at any time), but needs a training set spanning the cycle. Authors on other settings: "Whether BodyTime works similarly well in patients or cohorts with known lower circadian amplitudes (e.g., the elderly, shift workers) and/or in the presence of internal circadian desynchronization remains to be elucidated."

## Key results
- Internal cross-validation (leave-one-subject-out, RNA-Seq, 136 samples): median absolute error (MdAE) 1.4-1.7 h for internal time, "With 53.7%-59.6% of the predictions showing an error of" at most 2 h (Table 1).
- Migration to NanoString improved accuracy: The MdAEs improved "by about 1 hour ranging between 0.6 and 0.9 hours for the different predictors" (Fig 2C); NanoString MdAE 0.7-0.9 h.
- External validation (VALI, n = 28): morning samples MdAE 0.54 h (2-gene, 1-sample) and 0.69 h (12-gene); 100% and 96.3% within 2 h; afternoon 0.99 / 0.80 h; 2-sample 0.74-0.75 h (Table 2).
- Agreement with DLMO (circular Pearson r >= 0.9, P < 0.0001): "Moreover, Bland-Altman analyses showed that there is no systematic difference between the 2 methods".
- Accuracy does not depend on chronotype: "The accuracy of the BodyTime assay is the same for all chronotypes (Supplemental Figure 7), even the extreme ones as revealed in our external validation study."
- Candidate genes are mostly immune-related (119 genes, GO immune enrichment); only 4 known clock genes in the RNA-Seq set.

## Not addressed / limits
- One tissue (monocytes), one platform, young healthy (mostly male in training), no patients; no cross-tissue prediction and no multi-tissue comparison in the same person.
- Predicts DLMO (pineal proxy), not the phase of peripheral clocks; the authors state that "neither BodyTime nor DLMO determined from saliva melatonin concentrations reports the circadian phase of the suprachiasmatic nuclei (SCN)".
- Assumes normal phase angle between sleep and melatonin: robust for chronotype, untested under shift work or internal desynchronisation.
- Validation samples were taken in the morning and afternoon only; no night samples in VALI.

## Link to us
- **Intro claim(s):** supervised, timed-training, single-sample, blood-only phase predictor. Our refs line (12 individuals constant routine, 28 chronotypes, needs timed full-cycle training) is supported.
- **Relation to clocksync:** background / competitor for single-compartment phase (blood). It gives a clock time per sample, which presupposes a training set with that tissue and protocol; clocksync reads inter-compartment offsets without timed training.

## Ideas
- Benchmark: BodyTime-like within-person phase would be a monocyte reference if matched blood data exist; compare our blood compartment offsets with its DLMO-referenced phase.
- Hughey 2016 (cited here) reports that clock mutants change the SPC1-SPC2 circle radius and axis; a radius readout is analogous to our relative amplitude |z| (compare).
- Missing reference: Archer 2014 PNAS (sleep out of phase, blood transcriptome rhythms) cited as ref 31.
