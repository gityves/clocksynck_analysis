# archer2014mistimed
**Ref:** Archer et al., PNAS 2014 (111:E682-E691). DOI 10.1073/pnas.1316335111. Type: clinical (forced desynchrony, blood transcriptome). Lab: Dijk (Surrey). Data: GEO GSE48113.
**Source:** full text

## Question and design
Does sleep mistimed relative to the central clock disrupt the human circadian transcriptome? Twenty-two healthy volunteers in a laboratory forced-desynchrony protocol: baseline day, then sleep-wake, dark/dim light and meals on a 28 h day (sleep episode 4 h later each cycle). Sleep 2 in phase with the melatonin rhythm; sleep 5 began at noon, ~12 h out of phase. **Sampling:** timed; blood RNA (whole blood, microarray) at seven samples within each of the two conditions (in-phase around sleep 2/wake 3, out-of-phase around wake 5/sleep 5; spacing not stated in the text read), 19 participants with paired data; melatonin (DLMO) as central phase marker (in phase 21:59, out of phase 23:03 clock time, group means).

## Method (if any)
Per-subject time-domain test for one near-24 h oscillation per probe, prevalent circadian = overt in a minimum number of subjects (FDR 5%); mixed-model ANOVA (condition x time point); linear model of 24 h circadian drive plus sleep-wake drive per transcript.

## Key results
- Rhythmic transcripts fall from 6.4% to 1.0%: "Here, we show that mistimed sleep leads to a reduction of rhythmic transcripts in the human blood transcriptome from 6.4% at baseline to" 1.0%; 1,502 transcripts (1,396 genes) in phase vs 237 (228 genes) out of phase.
- "Of the transcripts that were classified as circadian while sleeping in phase, the majority (97%) became arrhythmic when sleeping out of phase with the melatonin rhythm." Only 39 rhythmic in both, 27 day-peaking; mean amplitude 0.5169 to 0.4264.
- Clock genes ARNTL, NPAS2, PER2, PER3, CRY2 rhythmic in phase; "The only clock genes that remained rhythmic when sleeping out of phase were NR1D2 (Fig. 1C) and CSNK1E"; NR1D1 became rhythmic.
- Central clock unaffected: melatonin amplitude and waveform unchanged; sleep lost only ~1.9 h in total, so the effect is mistiming, not sleep loss.
- ANOVA condition x time: 10,848 transcripts (34%); the linear model separates circadian-only and sleep-only transcripts; "Thus, robustly rhythmic transcripts whose expression profiles were not affected by forced desynchrony were largely those related to intrinsic blood-specific functions."

## Not addressed / limits
- One compartment (blood); loss of rhythmicity within a tissue, not between tissues. No peripheral tissue phase relative to the central clock is reported in hours; the summary reports no between-condition phase shift of the rhythmic transcripts.
- Rhythm loss may include dispersion between subjects (per-subject rhythmicity threshold, 7 samples, short window).
- Healthy young volunteers, 3 days of delayed sleep; no disease; clock gene rhythm lost in whole blood is not a clock-loss statement about blood cells.

## Link to us
- **Intro claim(s):** the refs_*.md summary (misaligned sleep strongly reduces rhythmic transcripts and shifts their phases) is OVERSTATED on phases: the paper reports loss of rhythmicity (6.4% to 1.0%) and amplitude reduction, not a phase shift; the phase message is that sleep-driven transcripts move with the 28 h schedule. Loss within one compartment is supported.
- **Relation to clocksync:** data we used. Done: A156 (GSE48113, within-subject day swaps; H156.1, H156.2 fail, MDE 13.1 h: core-cue change not detected).

## Ideas
- Done: the main re-analysis (A156). Possible companion: report the sleep-driven versus circadian-driven split (their Fig 6-7 classes, 1,792 transcripts with R2 > 0.6) as the human feeding/sleep versus clock decomposition, comparable with the GTEx cue/locked split.
- Missing references: Moller-Levet 2013 (sleep restriction, blood); a paired multi-tissue forced-desynchrony study does not exist in humans.
