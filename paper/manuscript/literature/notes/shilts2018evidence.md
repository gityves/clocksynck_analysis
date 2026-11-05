# shilts2018evidence
**Ref:** Shilts, Chen, Hughey, PeerJ 6:e4327 (2018). DOI 10.7717/peerj.4327. Type: method (CCD) + data.
Lab: Jake Hughey (Vanderbilt).
**Source:** full text

## Question and design
Is the clock disrupted in human tumours? Clock correlation distance (CCD): Spearman co-expression of 12 clock genes
in a group of samples vs a mouse reference (8 mouse datasets); validated in human tissues; 20 tumour/non-tumour
datasets (TCGA, GEO); lung module preservation.
**Sampling:** tumours "Presumably, most to all of the samples were acquired during surgical working hours (6:00 am to 6:00 pm)."

## Method
Group-level: correlation matrix of clock genes, distance to reference; ΔCCD tumour vs non-tumour with permutation.
No time labels, no phase per sample.

## Key results
- Motivation, verbatim: "Most computational methods for quantifying circadian rhythmicity and inferring the status of the clock require that samples be acquired over the entire 24-h cycle."
- Partial coverage: clock co-expression preserved in samples from only the first 8 h of day or night: "our approach can detect the signature of clock progression in groups of samples without using time of day information, even if the samples’ coverage of the 24-h cycle is incomplete".
- Non-tumour similar to reference "(P ≤ 0.005 for CCD relative to the mouse reference in 14 of 20 datasets", tumours lacked similarity "(P > 0.05 for CCD in 10 of 20 datasets"; tumour > non-tumour CCD "(P = 1.4 · 10−8 by paired t -test)".
- Not explained by clock-gene knockout alone (mouse KO comparison); rhythmic lung modules less preserved in tumours.

## Not addressed / limits (stated by authors)
- Phase: "cannot detect a phase difference between conditions".
- Group only: "our approach does not immediately lend itself to assessing clock progression in single samples".
- Large CCD does not mean arrhythmic, only altered phase relationships.
- Outlook, verbatim: "it may be possible to complement the 1CCD metric and infer clock progression in a given individual by directly comparing matched healthy and disease samples" (1CCD = ΔCCD in the PDF text) — this is what clocksync does.

## Link to us
- **Intro claims:** (i) the cleanest statement that existing methods need the full 24 h cycle; (ii) surgical cohorts are working-hours samples (6 am-6 pm); (iii) tumour clock co-expression is perturbed in many cancers.
- **Relation:** closest precursor of our idea (time-free, hospital-hours compatible). Differences to state: clocksync works within person (compartment vs host), estimates coupling, amplitude and phase offset, per compartment, with noise model; CCD is a group-level co-expression score that cannot see phase. Our model-free gate is in the CCD spirit.

## Ideas
- Benchmark: compute ΔCCD alongside clocksync coupling on TCGA pairs (A138): are low-coupling tumours the high-CCD ones? CCD is cheap (deltaccd R package).
- Cite their partial-coverage test as precedent for our hospital-hours simulation design.
