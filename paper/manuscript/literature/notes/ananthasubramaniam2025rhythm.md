# ananthasubramaniam2025rhythm
**Ref:** Ananthasubramaniam & Venkataramanan, PLOS Biol (2025). DOI 10.1371/journal.pbio.3003196. Type: method
(COFE) + data. Lab: Bharath Ananthasubramaniam (Humboldt Berlin).
**Source:** full text

## Question and design
Reconstruct population rhythms from untimed single samples per person; apply to 11 TCGA adenocarcinomas (>250
samples each: BRCA, KIRC, KIRP, PRAD, OV, UCEC, LUAD, COAD, THCA, BLCA, LIHC), proteome (TCPA), matched healthy
tissue. **Sampling:** surgical biopsies (unrecorded times); authors expect daytime procedures.

## Method
Sparse cyclic PCA (scPCA): simultaneously assigns time labels and selects rhythmic features de novo (no seed genes),
cross-validation chooses sparsity and detects invalid orderings. Cohort-level, one time label per sample.

## Key results
- Aim: "reconstruct circadian rhythms across cohorts" from "single high-throughput omics samples (without time labels)".
- Predicted tumour times span the whole cycle despite daytime surgery: "assuming these tissue biopsies were gained from procedures conducted during the day (standard in most countries), we expect the predicted time labels to be restricted to a fraction of the circadian cycle. Surprisingly, predicted time labels lay throughout the circadian cycle in all ACs".
- Tumour vs matched healthy: "The difference between the predicted time labels for healthy and tumor tissue was distributed throughout the circadian cycle in most ACs"; interpreted as tumour clocks "uncoupled from the systemic clock in the subject".
- Clock structure in tumours resembles healthy, with delays in NR1D2, PER2, TEF, DEC1.
- Rhythmic gene fraction: ovarian and prostate ">40% rhythmic genes", breast and uterine "less than 15%".
- Cell cycle scores rhythmic in all ACs except KIRC "(cosinor, p < 0.003)"; interpreted as clock-cell cycle coupling.
- Rationale against two-step methods: "This two-step procedure using ML-predicted covariates is known to be statistically suboptimal".

## Not addressed / limits
- **Key weakness for us:** full-cycle spread of predicted times in surgical cohorts is taken as biology (tumour uncoupled), but an orderer fitted to a hospital-hours sample may spread samples over the circle by construction; no simulation of daytime-only sampling is shown in the main text. Their own CV argument ("samples were indeed collected around the clock ... is unlikely") flags the tension.
- Healthy-vs-tumour comparison uses separately reconstructed orderings (reference phase of each ordering arbitrary up to the alignment).
- Li 2024 finds TCGA normal breast phases concentrated in the working day, in contrast.
- Rhythmic "cell cycle" scores could reflect the ordering being driven by proliferation variation rather than time.

## Link to us
- **Intro claims:** (i) TCGA is a hospital-hours cohort analysed with orderers; (ii) prior claim that tumour clocks are decoupled from healthy tissue — rests on orderings that need full-cycle coverage. Our refs entry (COFE, tumour vs normal, time-series-free) is correct but should add this caveat.
- **Relation:** direct competitor on the tumour-vs-normal question; already in A145 (CoFE).

## Ideas
- Already partly done: A145 (results_preliminary.md) W arm 08-18 h, coupling MAE CoFE 0.206 vs clocksync 0.060 (clock-gene-only simulation); A138 finds a tumour-normal clock match in TCGA pairs (H138.4). Still missing: COFE's own outputs (full-cycle labels, random tumour-normal differences) under daytime sampling, to show the mechanism of the false uncoupling.
- Key simulation: hospital-hours sampling (08-18 h) of a coupled tumour-normal system; run COFE: does it return full-cycle labels and random tumour-normal differences (false "uncoupling")? clocksync on the same data: coupling recovered. Strong figure for the technical-challenges section.
- Run clocksync on TCGA tumour/adjacent-normal pairs in the same 11 ACs: do we find coupling where COFE claims none?
