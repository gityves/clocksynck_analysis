# duan2024taufisher
**Ref:** Duan, Ngo, Karri, Tsoi, Gudjonsson, Shahbaba, Lowengrub, Andersen, Nat Commun 15:3840 (2024). DOI 10.1038/s41467-024-48041-6. Type: method + data. Lab: Bogi Andersen / John Lowengrub (UC Irvine).
**Source:** full text

## Question and design
Predict the circadian time of one transcriptomic sample (bulk, or scRNA-seq pseudobulk) from within-sample gene-gene differences. Benchmarks: 12 timed datasets (GSE54650/54651 kidney, liver, brainstem, cerebellum; mouse skin GSE38622; SCN GSE157077; human blood GSE56931, GSE113883), 100 random 80/20 splits. New data: scRNA-seq of mouse dorsal dermis, "every four hours for 72 hours under 12:12 light-dark cycle" (18 samples, GSE223109); scRNA-seq SCN (GSE117295; CT14-58, 12 samples, constant darkness). **Sampling:** timed training series 24-72 h at 2-6 h resolution; test samples one time point each.

## Method
Training: rhythmic 24 h genes (JTK_Cycle or Lomb-Scargle), functional data analysis to fill missing time points, within-sample normalisation (rescaled pairwise gene differences), PCA, multinomial regression on PC1-2. Test: trim to predictor genes, same normalisation, project, predict. Supervised; "it does not require the training data to be a complete time series". Unit of inference: one sample or one pseudobulk (not single cells). Cell types are handled one at a time.

## Key results
- 2 h accuracy counted as correct. "For the two human blood datasets, TimeSignatR, using two test samples, outperformed ZeitZeiger and tauFisher". "tauFisher achieved lower but comparable accuracy (difference < 10%) when compared to TimeSignatR in two of the remaining four mouse datasets."
- Cross-platform skin (microarray to RNA-seq, 8 tests): accuracy 0.75, RMSE 4.704 h.
- Bulk to scRNA-seq SCN: "In ten out of the twelve tests, tauFisher predicted a time that is within 2-hour of the labeled time" (accuracy 0.833, RMSE 1.936; beat TimeSignatR and ZeitZeiger).
- Dermis pseudobulk: "In 14 out of the 18 tests, tauFisher predicted circadian time within 2 hours of the labeled time" (accuracy 0.778, RMSE 2.198).
- Disturbed systems: trained on control skin, predictions on time-restricted feeding skin move away from sampling time, but "can only tell whether the test system is disturbed or not". Trained within the disturbed schedule, LOO accuracy 0.875 and 1.
- Cell-type heterogeneity: bootstrapped pseudobulks of dermal immune cells give more variable predicted phases than fibroblasts. "Incorporating tauFisher with bootstrapping revealed that circadian phase heterogeneity contributes to the dampened collective rhythm as well as fewer rhythmic genes" (Rao and Wallraff tests, 500 pseudobulk replicates per type, within each ZT; p-values in Source Data).

## Not addressed / limits
- Human blood is poor: "its performance on human blood samples is subpar".
- Assumes pairwise gene differences depend only on circadian time: "the differences between predictor genes are solely dependent on circadian time"; no donor centring.
- Trained on healthy data, a disturbed tissue is only mapped to the healthy time scale; a disease time scale needs disease time series. The authors note dampened rhythms in cancer make times hard to separate.
- No single-cell prediction (dropout). Cell-type comparison is a hypothesis; they "did not capture enough cells" per finer cell type.
- Cell types are not compared against the same animal's other compartment as a reference; each is projected on the same bulk skin reference.

## Link to us
- **Intro claim(s):** supervised, timed training, single bulk or pseudobulk sample. Supported. Add: pseudobulk SCN and dermis are mouse, full-cycle test sets, human accuracy poor.
- **Relation to clocksync:** competitor (compare). It is the closest published pseudobulk phase tool, and its within-sample normalisation removes global scale but not donor effects. Its cell-type result (predicted-phase dispersion within a time point, immune vs fibroblast) overlaps with our relative amplitude and coupling by compartment.

## Ideas
- Timed truth for compartments: GSE223109 (dermis, 18 samples, fibroblasts and immune cells) and GSE117295 (SCN). No hit for these accessions or for tauFisher in METHOD.md, registered_experiments.md, results_preliminary.md. Test clocksync coupling and |z| of immune vs fibroblasts against their reported dampening. Related done: A31 / V13 (mouse colon and skin cell types).
- Run tauFisher per organ on the baboon (A36) and mouse multi-organ series and compare the organ offsets with clocksync offsets (also under 8 h windows, A40).
