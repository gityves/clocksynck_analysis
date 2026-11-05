# kervezee2018simulated
**Ref:** Kervezee et al., PNAS 2018 (115:5540). DOI 10.1073/pnas.1720719115. Type: clinical / data (GEO GSE107537). Lab: Cermakian / Boivin (McGill).
**Source:** full text

## Question and design
Human PBMC transcriptome (microarray), n = 8 healthy subjects (mostly male), within-subject: baseline vs day 4 of a simulated night-shift schedule (habitual sleep delayed by 10 h). **Sampling:** timed; blood every 4 h over two 24 h measurement periods (baseline and after the third shifted night), in a time-isolation laboratory; 16 h constant posture with hourly snacks, plus a sleep episode, in each period. Times are expressed relative to habitual bedtime.

## Method (if any)
Group-level mixed-effects cosinor per condition, per-subject cosinor with Fisher combination, then a BIC model-selection scheme comparing both conditions simultaneously (10 categories: rhythm in both / one / none, same or shifted phase, same or different mesor). Needs sampling through the full 24 h in both conditions.

## Key results
- At FDR 0.05 rhythmic probe sets fell from 444 (3.8%) at baseline to 62 (0.5%) in the night-shift condition; rhythmic transcripts across subjects (Fisher) fell from 11.8% to 6.5%; 7 of 8 subjects lose rhythmic probe sets.
- With model selection, 1,092 probe sets (9.7%) are rhythmic; "73% of transcripts rhythmically expressed at baseline remained rhythmic during the night shift condition with a similar phase relative to habitual bedtimes, but with lower amplitudes" (abstract). So the binary loss is largely dampening: mean amplitude 0.185 to 0.115 log2 (471 probe sets, P < 0.0001).
- Only 24 probe sets shift; 23 of 24 are phase-delayed by 5.4 h (circular SD 1.2); 7 are NK-cell genes (PRF1, GZMB, GNLY, KLRD1) that move together. Their shifted rhythm is more likely to follow the shifted rest-activity and feeding cycle than the central clock.
- 331 probe sets change mean level (154 down, 177 up), with JUN down and STAT1/2/5A up.
- Peaks are bimodal at about 2 h and 14 h after habitual bedtime; the central clock (melatonin) did not adapt to the shifted schedule (cited from earlier work on the same subjects).

## Not addressed / limits
- Sleep-time is the anchor, not the clock: a transcript with "similar phase relative to habitual bedtimes" stays with the (unshifted) central clock; the main message is blood stays on central time while behaviour shifts, and a small subset follows behaviour. Different blood genes shift by different amounts (internal desynchrony within PBMC), but no second tissue.
- n = 8; one subject who gained rhythmic transcripts was female; sex not analysable.
- PBMC composition change under shifted sleep is not modelled. About 38 min less sleep on the last night (partial sleep deprivation is confounded with misalignment).
- Four days only: no information on full re-entrainment.

## Link to us
- **Intro claim(s):** refs row: PBMC transcriptome becomes misaligned, with gene-specific phase shifts and individual differences. Supported for gene-specific shifts (24 genes, 5.4 h); the majority (73%) stay put, which the refs text underplays. The grouped refs line two compartments ... in the same volunteers is WRONG for this key: PBMC only (one compartment); n = 8 is correct.
- **Relation to clocksync:** data we could use (GSE107537 listed as secondary in provisory_plan_paper.md, A124 line); background.

## Ideas
- Not done: GSE107537 is named in provisory_plan_paper.md as a descriptive secondary (blood readout weak, A35); no amendment scores it. Within-person test: does the untimed core-cue offset (A156 estimator) change day 1 to day 4 in these 8 subjects? Power is lower than A156 (MDE 13.1 h, 19 subjects), so only as descriptive.
- Use the 23 delayed vs 676 unshifted transcripts as a gene-set pair: cue-like vs core-like genes in PBMC (compare with the A81b cue/locked states).
