# winkler2022gene
**Ref:** Winkler E, Ananthasubramaniam B, Herzel H, bioRxiv 2022 (10.1101/2022.11.02.514919). Type: biology / method application (pathway scoring). Lab: Hanspeter Herzel / Bharath Ananthasubramaniam (Berlin).
**Source:** full text (preprint)

## Question and design
Do gene-expression signatures of cancer-related signalling pathways (PROGENy, 14 pathways) show 24 h rhythms in mouse liver? Datasets: Atger 2015 (ad libitum and restricted feeding), Hughes 2009 and Zhang 2014 (microarrays). **Sampling:** published circadian time courses, full 24 h (CT) with replicates. Mouse only; no human or tumour data.

## Method
- Not phase inference. PROGENy score = weighted sum of expression of about 100 signature genes per pathway for each sample (extra normalisation: each gene divided by its mean); then RAIN rhythmicity on the score time series (q < 0.05), plus gene shuffling controls and a jackknife (drop one signature gene).
- Unit: sample in a timed time series. Supervised by known timing. No relevance to untimed or single snapshots.

## Key results
- Atger liver: "Out of 14 pathways, 9 were considered significantly rhythmic according to our criteria"; shuffled data gave on average no rhythmic pathways.
- Androgen and estrogen scores peak near CT13; EGFR and PI3K peak in the night (CT14-20).
- EGFR, PI3K and p53 are rhythmic in all liver datasets tested; peaks fall in the active phase.
- Rhythmicity is not due to one gene: removing any single signature gene "did not abolish rhythmicity".
- Delay between pathway activation and score peak estimated at about 4 h (Akt, PI3K) to 6-8 h (Erk).

## Not addressed / limits
- Mouse liver only; no human, no tumours, no tissue comparison, no snapshot data. The authors state that expression signatures may not equal pathway activity (post-translational control).
- Score rhythms partly reflect clock genes in the signature (Nr1d2 in EGFR), so they are not independent of clock phase.

## Link to us
- **Intro claim(s):** our refs line (predicts rhythms in oncogenic pathways in snapshot, TCGA-like, data) is WRONG. The paper uses timed mouse-liver time series, not snapshot or TCGA data; it supports only that oncogenic signalling pathways (EGFR, PI3K, p53) are rhythmic in liver, i.e. why clock disruption may matter for cancer.
- **Relation to clocksync:** background. If cited as a snapshot method it must be dropped; fine as a biology reference.

## Ideas
- Possible secondary use: PROGENy scores per compartment as a readout of whether the clock-derived host state predicts pathway activity in tumours (post hoc, would need registration).
