# zhang2014circadian
**Ref:** Zhang, Lahens, Ballance, Hughes, Hogenesch, PNAS 111:16219 (2014). DOI 10.1073/pnas.1408886111. Type: data. Lab: Hogenesch (Penn).
**Source:** full text (main text; SI and Methods not in the export)

## Question and design
Circadian transcriptome atlas of 12 mouse organs (adrenal gland, aorta, brainstem, brown fat, cerebellum, heart, hypothalamus, kidney, liver,
lung, skeletal muscle, white fat) by RNA-seq and Affymetrix arrays; GEO GSE54652. **Sampling:** timed, 24 h time courses in healthy wild-type
mice: "We sampled organs every 6 h by RNA-seq and every 2 h by arrays, to develop an atlas of mouse biological space and time." The
conditions (constant darkness, ad libitum feeding, 2 days) are not in the extracted text; Yeung 2018 states the Zhang set was DD, ad libitum, every 2 h.
Number of animals per time point is not given in the main text.

## Method (if any)
JTK_CYCLE on the array data (24 h period, FDR 5 %); RNA-seq used for splice forms and non-coding RNAs. Phase distance between organs from shared
rhythmic genes (Fig. 2B). No inference of time or of organ phase offsets from untimed data.

## Key results
- Breadth: "Transcript abundance for 43% of protein-coding genes oscillated in at least one organ" (19,788 background genes; Fig. 1B). Extrapolation: about 10,901 genes (55 %) with more organs (asymptote of y = a(1-exp(-bx)), R2 > 0.99).
- Per organ: "Following these analyses, we found liver had the most circadian genes (3,186), whereas hypothalamus had the fewest (642)".
- Organ specificity: only 10 genes cycle in every organ ("all organs: Arntl, Dbp, Nr1d1, Nr1d2, Per1, Per2, and Per3 (core"; plus Usp2, Tsc22d3, Tspan4).
- Clock genes are synchronous across organs: "of a given gene falling within 3 h of each other across all organs" (Fig. S5); 1 to 2 h advances in skeletal muscle and delays in cerebellum are called not resolvable at 2 h resolution.
- Output genes are not: "Remarkably, 1,400 genes were phase-shifted with respect to themselves by at least 6 h between two organs, with 131 genes" antiphased (e.g. Vegfa, brown fat vs heart). Authors attribute this to organ-specific timing between core clock and outputs.
- Phases cluster before dawn and dusk ("rush hours"), heart and lung diverge; tissue tree from mean phase differences recovers organ lineage (Fig. 2B).
- Medicine: most best-selling drugs and 119 WHO essential medicines target products of rhythmic genes.

## Not addressed / limits
- Healthy, entrained wild-type mice only; no disease, no human tissue, no inter-individual variation in organ phase (organ means per time point).
- Brain organs are mixtures: "regions, we assume that heterogeneous mixtures of cell types" may have different sets of genes or phases, which limits detection there.
- Resolution: "precisely resolve small (<2 h) phase differences from data with" 2 h sampling, so tissue offsets below 2 h cannot be claimed.
- Compares organs only through gene-level phases; never estimates a per-organ clock state or its coupling to a host.

## Link to us
- **Intro claim(s):** refs_circadian_desynchrony.md and refs_cancer_phase_inference.md claim: 43 % of protein-coding genes rhythmic in at least one of 12 mouse tissues, mostly tissue-specific. SUPPORTED (43 %, 12 organs, mostly organ-specific confirmed); the 43 % is for array JTK_CYCLE at FDR 5 %.
- **Relation to clocksync:** background and data we could use. It documents that core clocks stay within about 3 h across organs in a healthy mouse (the shared host state) while outputs diverge, which is the premise of reading each compartment against the host.

## Ideas
- Use GSE54652 as a ground-truth test of clocksync on mice: per-time-point organ pseudobulks of clock genes, check c, |z|, arg z against the known 3 h clock-gene spread and the muscle/cerebellum 1 to 2 h offsets (with only 24 samples per organ, one donor per time point).
- Simulate hospital hours by keeping only 8 to 12 time points of the day and verify that offsets survive (our central claim for untimed day-only sampling).
- Citable for tissue specificity of rhythmic genes together with Mure 2018 (baboon, 82 %) and Yeung 2018.
