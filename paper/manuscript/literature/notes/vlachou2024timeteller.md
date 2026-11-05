# vlachou2024timeteller
**Ref:** Vlachou D et al., PLoS Comput Biol 2024;20(2):e1011779. DOI 10.1371/journal.pcbi.1011779. Type: method (TimeTeller, supervised phase plus clock-dysfunction score). Lab: Rand (Warwick) with Levi, Dallmann, Bjarnason.
**Source:** full text (main text; the S1 Appendix with dataset notes is not included)

## Question and design
Estimate clock function from a single transcriptome by modelling the multi-dimensional state of a small panel of rhythmic genes. Training datasets (four): Zhang et al. mouse microarray and RNA-seq (several tissues, 24 h series), Mure et al. baboon (33 tissues of 64, one animal sampled every 2 h over 24 h), and Bjarnason et al. human timecourse data (pooled across individuals; tissue and n not stated in the main text). Test: mouse Fang, Kinouchi (fasting), Acosta-Rodriguez (calorie restriction/time-restricted feeding), Koronowski (Arntl KO), and human oral mucosa/buccal biopsies (Boyle: 39 smokers vs 40 never smokers; Feng: oral squamous cell carcinoma, 167 cancer vs 45 normal mucosa; Cadenas: breast cancer).
**Sampling:** training data are timed, covering 24 h (human RNA-seq training only "at a few training time points around the day"); cancer test data are untimed-in-practice (collection times unknown) and are analysed without time stamps.

## Method (if any)
Supervised: from training timecourses it learns, per time, the mean and covariance of the normalised gene vector (REV), splines them over the day, and gives for a test sample a likelihood curve over time. Outputs: phase T (maximum of the curve), maximum likelihood ML and a dysfunction metric Theta (proportion of the time the likelihood stays above a reference curve, i.e. width/multiple peaks), plus "timing displacements" per individual, tissue or condition, and phase displacement plots.
Single sample at test time; needs timed training data across the day. Cross-tissue use: pooling tissues requires a gene panel with "good cross-tissue synchronicity" and normalisation; tissue-matched (timecourse) normalisation needs a training series for each test tissue, otherwise intergene normalisation (centre/scale of the gene vector) is used.
Authors on other settings: they explicitly aim at patients and cancer, but "TimeTeller's view of functionality is based on statistical analysis of gene expression and not on timing of physiological processes."

## Key results
- Accuracy: baboon leave-one-tissue-out MdAE around 1.23 h; correcting by tissue displacement gives "a mean absolute error of just over one hour"; human training MAE reduced about 50% after adjusting for individual displacements. Compared with TimeSignature 1.21-1.49 h and ZeitZeiger 0.6-1.1 h (mouse); "we do not wish to claim timing superiority".
- Between-tissue desynchrony in one animal (baboon): many tissues had "a substantial timing displacement (Fig 5H and 5J) ranging from approximately -3.5h to +2.5h compared to the time the samples were taken", the null that the m-th most advanced and m-th most retarded tissue share a displacement is rejected at p = 0.01 for all m < 7. A phase displacement plot gives low r2: tissue differences are "more complex than a simple phase shift in the core clock genes".
- Humans: individuals have consistent timing displacements (molecular chronotype).
- Cancer: oral mucosa, "a highly significant difference in median Θ values between the cancer group (167 individuals) and the the normal mucosa group (45 individuals) (p < 0.002)"; ML not different. Smokers vs never-smokers differ in Theta (Wilcoxon).
- Timing in cancer: "109 samples whose timing T is before 7am and 93 of these have a second peak"; after choosing the second peak, times fall in the early afternoon, within about 0.4 h SD (0.83 h for the breast set).
- Fasting mouse muscle (Kinouchi): MAE 0.44 h fed vs 4.04 h fasted.

## Not addressed / limits
- Not independent of the training panel: phase T is relative to the training clock structure; displacements are relative to external sampling time that is unknown for cancer data, so the absolute T is unvalidated there.
- Theta mixes true clock dysfunction with low amplitude, noise, tissue composition and sample type (tumour vs mucosa differ in cell mix); no specificity control against these (our Statistics rule).
- Human training set is small and heterogeneous; RNA-seq training time points are few.
- Between-tissue displacements shown only in a single baboon (not humans); no human multi-organ same-donor data.
- No demonstration that the predicted clock reflects physiology rather than expression differences.

## Link to us
- **Intro claim(s):** supervised, timed-training single-sample method giving clock phase and a dysfunction metric, applied to cancer. Our refs line is supported; the claim of cancer prognosis application is overstated: the cancer analyses here are Theta distributions in oral and breast tumours versus normal mucosa, with no survival or prognosis analysis in the text read.
- **Relation to clocksync:** competitor (compare). Same goal (dysfunction in tumours without time stamps) but a single-compartment, training-referenced measure. Its baboon result (-3.5 to +2.5 h tissue displacements) is the best existing support for inter-tissue offsets that clocksync estimates in humans.

## Ideas
- Run TimeTeller (R package on GitHub) on our tumour/normal bulk sets and compare Theta, ML and T with clocksync coupling c and arg z, using the same permutations (specificity companion).
- Check that tumours are mistimed by second-peak artefacts (93 of 109) in our data: a clock panel with a wrong peak equals a failed gate, not decoupling.
- Use Mure baboon 33-tissue displacements as an organ-offset benchmark for clocksync (mouse/baboon multi-organ truth for arg z).
- Check the S1 Appendix (Note A) for the human Bjarnason et al. training series (n, tissue, sampling times), which is not in the main text.
- Done: TimeTeller is in the GTEx hidden-times benchmark (V2, METHOD.md).
