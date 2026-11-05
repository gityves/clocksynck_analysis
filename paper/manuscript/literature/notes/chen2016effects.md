# chen2016effects
**Ref:** Chen et al., PNAS 2016 (113:206). DOI 10.1073/pnas.1508249112. Type: data / biology (GEO GSE71620). Lab: McClung (Pittsburgh).
**Source:** full text

## Question and design
Human postmortem prefrontal cortex, Brodmann areas 11 and 47 (BA11, BA47), Affymetrix Human Gene 1.1 ST arrays. 146 individuals without lifetime psychiatric disorder, mean age 50.7 (16-96) y, 78% male, postmortem interval 17.3 h (4.8-28 h); only witnessed, rapid deaths were kept (64 of 210 removed). **Sampling:** untimed by design, one sample per donor ordered by time of death converted to zeitgeber time from sunrise at the place of death (ZT -6 to 18); coverage of the 24 h day not given as a per-hour count in the text.

## Method (if any)
Nonlinear least-squares sinusoid per gene and region, empirical P from 1,000 TOD-permuted datasets, adaptive weighted Fisher meta-analysis across BA11 and BA47. Age effect: younger (< 40 y, n = 31) vs older (>= 60 y, n = 37) fitted separately; age-effect P from 1,000 age-shuffled matrices.

## Key results
- Rhythmic (P < 0.05): 2,475 genes in BA11 (about 12%) and 1,615 in BA47 (about 8%); "a core set of 235 genes (q < 0.05) with significant circadian rhythms of expression", with 92% concordant acrophase between areas (Pearson r = 0.95).
- Age-dependent rhythm change in 1,186 genes (BA11) and 1,591 (BA47): in BA11, 588 lose rhythmicity, 533 gain it, 67 shift phase, 30 change level, only 2 lose amplitude without loss (BA47: 1,063 lost, 434 gained, 94 phase shifts).
- Clock genes: PER1 acrophase moves from ZT5 (BA11) and ZT7 (BA47) to ZT1-2 with lower amplitude; PER2 from about ZT12 to ZT6; PER3 shows no obvious age effect; CRY1 loses rhythmicity in BA11 and ARNTL changes level in BA47.
- "a set of transcripts gained rhythmicity in older individuals, which may represent a compensatory mechanism due to a loss of canonical clock function" (abstract; interpretive).

## Not addressed / limits
- Two neighbouring prefrontal areas of the same donors; no other regions or tissues, so no internal desynchrony across organs.
- Old vs young is a between-person contrast of 31 vs 37 donors, with the intermediate 40-60 y excluded; rhythm calls at P < 0.05 per group, so "loss" and "gain" partly reflect power and false positives (the authors note the paper cannot detect rhythms without lifestyle, actigraphy or sleep data).
- Time-of-death ordering treats a donor's clock as sun-entrained; an age-related phase shift of the donor's own clock cannot be separated from a change in the gene's rhythm within the donor.
- Heterogeneity in death time and cause between age groups is not modelled.

## Link to us
- **Intro claim(s):** refs row n=146 by time of death: ageing reduces/shifts rhythmicity of gene expression. Mostly supported; precisely: ageing causes loss, gain and phase advance of PER1/PER2 rather than a uniform amplitude reduction (only 2 genes with amplitude reduction in BA11).
- **Relation to clocksync:** data we could use (GSE71620; the only timed cortex study named in A117); background for age effects on clock-gene phase.

## Ideas
- Planned, not scored: A117 (provisory_plan_paper.md) H117.2 (log2 amplitude, older vs younger, pooled over regions, replicating A3) and H117.1 on regions; not in registered_experiments.md as an amendment (grep GSE71620: no hits).
- Phase check for the A130 chronotype hypothesis (internal time advances with age): Chen reports PER1 advanced by about 4-6 h in older donors (ZT5/7 to ZT1-2); compare with the untimed per-donor internal-minus-wall time against age (A130 design) as a timed positive.
- Missing reference: Seney 2019 (prefrontal cortex time of death, schizophrenia; already scouted in A117).
