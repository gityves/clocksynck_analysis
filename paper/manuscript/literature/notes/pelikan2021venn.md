# pelikan2021venn
**Ref:** Pelikan, Herzel, Kramer, Ananthasubramaniam, FEBS J 289:6605-6621 (2022; online 2021). DOI 10.1111/febs.16095. Type: method (viewpoint + benchmark). Lab: Ananthasubramaniam / Herzel (Berlin).
**Source:** full text

## Question and design
Do Venn-diagram analyses (VDA: separately detected rhythmic lists, then overlap) overestimate rhythm reprogramming? Two artificial scenarios from the 48 h hourly mouse liver series of Hughes 2009 ("Mouse liver transcripts were quantified every hour for 48 h"): odd versus even time points (no true difference) and a spike-in of 500 transcripts with changed amplitude/phase (450 detectable). Reanalysis of three liver diet/KO studies (high-fat diet, ketogenic diet, p66Shc KO) plus one validation set. **Sampling:** timed, regular, one or more periods per condition.

## Method (if any)
compareRhythms R package. (i) Hypothesis testing: joint harmonic-regression test of equal vs different amplitude/phase in transcripts rhythmic in either group (limma, DESeq2/edgeR/voom, or DODR/RAIN on normalised data). (ii) Model selection: nested cosinor models (AIC/BIC) with a margin over the second best. Both give four categories: gain, loss, change, same (VDA gives three).

## Key results
- No true difference: "called 2296 DiffR hits or 48% of all rhythmic transcripts" for VDA at FDR 0.05; more than 40% over all thresholds (Fig. 1B-D). Hypothesis testing: 0 false hits; model selection: 189.
- Spike-in: "Only a small fraction of hits (8.3%) were true DiffR transcripts." (Results); "VDA never recovered more than 50% of true DiffR transcripts, and no more than 20% of the DiffR hits were truly DiffR" (Fig. 1F); random choice would give about 11%.
- Amplitude thresholds do not rescue VDA (still about 40% false).
- Hypothesis testing "almost perfectly recalled 70% of the true DiffR transcripts or 80% of the true DiffR transcripts with precision above 80%"; model selection had about 50% precision at 75% recall (Fig. 3, Fig. S2).
- Reanalysis: hypothesis testing/model selection call 8-23%/34-51% of rhythmic transcripts DiffR versus 72-85% reported by VDA (Table 1); "Only 25 of 90 DiffR hits from hypothesis testing were also DiffR hits in the original VDA" for HFD (Fig. S3F); same-condition WT comparison gives 75% (VDA) vs 6% (testing).

## Not addressed / limits
- Needs two timed designs on the same clock; liver only; assumes 24 h harmonic shape.
- The authors' own best method (hypothesis testing) recovers only part of true DiffR; they assume 50% in conservative statements.
- The error logic is general, but it is shown only for rhythms; no single-cell or untimed data.

## Link to us
- **Intro claim(s):** our entry says the model-selection recommendation: partly overstated. They recommend both approaches and rank hypothesis testing higher on precision; model selection over-calls. Supports independent thresholding of two conditions inflates differences.
- **Relation to clocksync:** supports our rule for group comparisons (compare directly, at common spread, with a null on the same permutations, not by separately thresholded gates).

## Ideas
- Venn-type control in our own disease contrasts: show that comparing the numbers of gate-passing donors/compartments between disease and control is the VDA error; direct permutation contrasts (already used, e.g. A33, A151) avoid it. Cite in Methods.
- Run compareRhythms on timed liver tissues (Hughes LD/DD style) versus our coupling difference to give a timed counterpart of a differential-rhythmicity call; not registered.
