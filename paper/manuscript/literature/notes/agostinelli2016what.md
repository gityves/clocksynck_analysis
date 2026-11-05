# agostinelli2016what
**Ref:** Agostinelli, Ceglia, Shahbaba, Sassone-Corsi, Baldi, Bioinformatics 32:i8-i17 (2016). DOI 10.1093/bioinformatics/btw243. Type: method. Lab: Pierre Baldi (UC Irvine).
**Source:** full text

## Question and design
Two deep-learning tools: BIO_CYCLE (is a time series periodic; period, phase, amplitude) and BIO_CLOCK (time of a single-time-point
expression sample). Data: CircadiOmics. BioCycleReal is built from 36 circadian transcriptome datasets, all mouse except one
Arabidopsis set. **Sampling:** training data are timed series. "BioCycleReal comprises experiments carried over a: 24-h period with a 4 h sampling rate"; also 48 h at 2 h and 1 h. BioClock is mouse single-time-point measurements labelled with their time, by tissue and condition (liver, kidney, heart, colon, glands, muscle, bone, fat, SCN, hippocampus, hypothalamus, cerebellum, fibroblasts, macrophages). Test samples: one time point each; no coverage needed.

## Method
BIO_CLOCK is a supervised DNN with two coupled outputs (cos, sin of the phase), converted into ZT. "BIO_CLOCK is trained using 16 core clock genes" (Arntl, Per1-3, Cry1-2, Nr1d1-2, Bhlhe40-41, Dbp, Npas2, Tef, Fmo2, Lonrf3, Tsc22d3). 70 % / 30 % train/test, 5-fold CV for architecture; one predictor per tissue or one on all tissues. Unit of inference: one expression sample. Compartments of one individual are not compared.

## Key results
- Overall: "it accurately predicts the time of the experiment with a mean absolute error of 1.22 h (less than 75 min)" (train and test on all mouse tissues, WT).
- Gene number: "Adding Nr1d1 and Nr1d2 to this set reduces the mean absolute error to 1.65 h" (from 3.72 h with six core genes).
- Cross-tissue: "The DNNs perform significantly worse when trained and tested on data with diverging origins." Liver, brain and two pooled-organ sets; the all-tissue DNN does about as well as the specialist (Table 6; numbers not usable from the extracted table).
- Cross-condition (diet, epilepsy, SIRT1/SIRT6): "This experiment yields a mean absolute error of 2.57 h."
- BIO_CYCLE: DNN AUC above JTK_Cycle, ARSER, Lomb-Scargle, MetaCycle on synthetic data (e.g. 0.92 vs 0.85-0.87, 24 h at 4 h sampling), and similar on real data (AUC 0.94-0.98, t-test P 0.3-0.7).
- All mouse GEO experiments annotated with an inferred time stamp.

## Not addressed / limits
- Mouse only. The conclusion names human as future work: "developing methods for annotating the time of all the human gene expression experiments" would be valuable. Our refs line says (mouse, human): human is not tested.
- Training on all tissues shares peripheral clock phase: the model implicitly assumes tissue-independent phase relations. No test of phase offsets between tissues of one animal.
- Disease, tumour and shifted-clock conditions only through mouse knockouts and diets (MAE 2.57 h, higher).
- Test and training time points come from the same ZT scale; no hospital-hours or partial-coverage test (not needed, but not examined either).
- Per-sample vs per-gene normalisation is not specified clearly ("training and testing examples are normalized").

## Link to us
- **Intro claim(s):** supervised phase inference from a single sample with timed training data (technical challenges). Supported for mouse; **CLAIM OVERSTATED** for "(mouse, human)": the paper is mouse only.
- **Relation to clocksync:** background / competitor on the single-sample time task. It estimates ZT in absolute time per sample, so it needs a timed training series from the same species; clocksync infers relative offsets between a donor's compartments.

## Ideas
- Comparator on a timed multi-organ dataset (mouse or baboon): BIO_CLOCK or tauFisher per organ, then the difference of the predicted time between organs of one animal vs the clocksync offsets. Not in METHOD.md, registered_experiments.md or results_preliminary.md; baboon timed truth exists (A36, A40).
- Use the 1.2-3.7 h error scale and the cross-tissue degradation as a reference for how much absolute-time inference loses compared with relative offsets.
