# braun2018universal
**Ref:** Braun R et al., Universal method for robust detection of circadian state from gene expression, PNAS 2018 (TimeSignature). DOI not stated in the extracted text. Type: method (supervised). Lab: Allada / Braun (Northwestern).
**Source:** full text

## Question and design
Human whole-blood transcriptomes. Training: random 24 time courses of the Moller-Levet sleep restriction/extension microarray study (GSE39445). Test: the 24 remaining time courses; validation: V1 (GSE48113, forced desynchrony), V2 (GSE56931, sleep deprivation), V3 (new RNA-seq, 11 healthy subjects).
**Sampling:** timed, serial blood draws "every 2-4 h over a period of one or more days" (about 7 to 18 samples per subject), all full-cycle, laboratory sleep-intervention protocols.

## Method (if any)
Supervised: log2 data, within-subject renormalisation (each gene centred on its own mean over time in that subject), then elastic-net (group lasso on cos/sin pair) regression of clock angle on about 40 genes (41-gene panel). Output: clock time of each sample (hours of day).
Needs timed training data across 24 h (one study) but then applies to new studies without retraining; needs at least two draws per test subject (ideally 8-16 h apart, 12 h best) to estimate the subject mean. Tissue: whole blood only; the authors say "the predictor was designed using gene expression profiles in whole blood and so may require retraining to be applied in other tissues."

## Key results
- Median error under 2 h in each dataset: "with a median absolute error across all datasets of 1:37" (h:min); validation sets 1:21 to 1:42.
- Cross-study and cross-platform without recalibration: "The previously trained predictor was applied directly to each of these new datasets, without any retraining of the machine." Validation errors not significantly different from test set (Wilcoxon p = 0.1, 1, 1).
- Beats ZeitZeiger (median absolute errors above 3 h) and PLSR (about 2.25 h) on the same data; fewer genes (41 vs 100).
- Disrupted conditions: forced desynchrony median error 2:05 h vs 1:26 h controls, "with the predicted times in the FD subjects lagging the true times by ∼2 h on average"; consistent with predicting the internal clock.
- Two draws: accuracy depends on spacing; 8 to 16 h often near-optimal.

## Not addressed / limits
- Healthy adults in sleep-manipulation protocols only. The authors: "to date no method (including ours) has been tested across a broad range of conditions and diseases, in particular those that are associated with circadian dysregulation."
- Whole blood only; no other tissue and no comparison across tissues in the same individual.
- Within-subject renormalisation removes the mean and so cannot be used where the mean carries the signal; needs 2 draws spanning a large part of the day, i.e. not hospital-hours single biopsies.
- Predicts the blood clock phase referenced to clock time, not DLMO.

## Link to us
- **Intro claim(s):** supervised, timed-training predictor of blood circadian time that generalises across studies. Our refs line (within ~2 h, platform-independent without retraining, timed training required) is supported (median 1:21-1:42 h; "most" samples within 2 h).
- **Relation to clocksync:** background / competitor for single-tissue phase; its within-subject centring plus 2 draws 12 h apart is a different way to cancel constants than our cross-compartment host reference.

## Ideas
- Apply the trained TimeSignature panel to GTEx whole blood (if clock times known, check error in hospital-hours samples) as a baseline against clocksync blood offsets.
- Forced-desynchrony V1 (GSE48113) gives a human internal-desynchrony test with a lag of about 2 h: usable as a positive control for a drift of blood relative to the sleep-wake rhythm.
- Missing reference check: Ueda 2004 molecular timetable (not in this batch).
