# li2013circadian
**Ref:** Li et al., PNAS 2013 (110:9950). DOI 10.1073/pnas.1305814110. Type: data / clinical (GEO GSE45642). Lab: Akil / Bunney (Pritzker consortium).
**Source:** full text

## Question and design
Human postmortem brain, six regions (dorsolateral prefrontal cortex, anterior cingulate, hippocampus, amygdala, nucleus accumbens, cerebellum), microarray, about 12,000 transcripts. 55 controls (52 with DLPFC) and 34 patients with major depressive disorder (33 in DLPFC); 29-55 samples per region. **Sampling:** untimed by design: each donor is one point, ordered by time of death (TOD), converted to zeitgeber time with sunrise = ZT0; pseudo-time series covering the 24 h day, average 1.2-2.3 points per hour; death hours come from the donor, not the protocol.

## Method (if any)
Sinusoidal regression of expression on TOD per region, empirical p by permuted TOD, Fisher combination across regions; TOD predicted from the top 100 cyclic genes (trained on 60 randomly selected subjects) and deviation between predicted and recorded TOD used as a measure of desynchrony with the solar day. Assumes the donor's clock is entrained to the sun at death.

## Key results
- Cyclic transcripts (P < 0.05) in controls: 922 DLPFC, 417 amygdala, 444 cerebellum, 565 NAcc, 566 ACC, 659 hippocampus. ">100 transcripts exhibited consistent rhythmicity and phase synchrony across regions" (abstract); 445 transcripts rhythmic in at least two regions with similar peak time and amplitude.
- ARNTL is the top gene in all six regions; PER1 peaks 0-2 h after sunrise, PER3 midday, PER2 afternoon; human peaks are about 6.5 h earlier than mouse (slope 1.18, r = 0.88; 7 genes).
- MDD: of the top 16 control genes, 11 are significant in 4 or more regions in controls but "only 2 had P < 0.05 in more than one region in patients with MDD"; the rank of ARNTL in MDD is 171st, PER2 532nd, PER3 10,191st.
- Predicted minus recorded TOD is larger in MDD than controls (P = 0.012; Mann-Whitney); in the 10 untreated suicide, toxicology-negative cases it is 3.3 h vs 1.9 h in all controls (P = 0.068) and 2.1 h in daytime controls (n = 30, P = 0.038).
- Gene-gene correlations are partly kept in MDD (Mantel statistic 0.38, P < 0.001) with some pairs reversed; this supports both a phase shift in many patients and disrupted gene-gene relationships.

## Not addressed / limits
- Between-region phase offsets are described as consistent, but no per-donor between-region offset is estimated; regions are analysed separately and combined by p-values.
- No record of sleep or consciousness at death; deaths include suicide, drug and medication effects (antidepressants, toxicology) confounded with MDD; TOD is not an independent measure of the clock.
- MDD signal loss could be a flattened rhythm or shifted phases; the data cannot separate them cleanly (the authors argue both). The smaller MDD group is checked by subsampling controls.
- Cerebellum and amygdala (smallest n) are the weakest.

## Link to us
- **Intro claim(s):** refs row: six regions, 55 controls: rhythms by time of death; disrupted in MDD. Supported (55 controls and 34 MDD, six regions). Add: this is the closest published result to our ask: phase-synchrony across regions in controls and a loss of it in MDD, with a time-of-day-free correlation analysis that is "immune to any desynchronization" (gene-gene correlations).
- **Relation to clocksync:** data we could use (GSE45642) and a conceptual precedent for host-referenced readout; the region-consistent phases argue that within-brain offsets are small in controls.

## Ideas
- Planned, not scored: A117 (timed brain regions) in provisory_plan_paper.md lists Li 2013 as scout; H117.1 untimed between-region offsets vs timed region phases, descriptive MDD vs controls. Not in registered_experiments.md as an amendment (grep GSE45642: no hits).
- Within-subject region pairs: does clocksync coupling c (six regions per donor) fall in MDD, at matched spread (specificity companion: non-clock gene-pair coupling)? Needs per-donor region matching from GSE45642.
