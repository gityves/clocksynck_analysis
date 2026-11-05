# Verified references: Introduction, cancer interaction and phase inference

Exact CoPhaser citation: Paychere Y, Salati A, Gobet C, Naef F. CoPhaser: generic modeling of biological cycles in scRNA-seq with context-dependent periodic manifolds. bioRxiv 2025 (posted 26 Dec 2025; v2 exists). doi:10.64898/2025.12.24.696376. Key: paychere2025cophaser.
scRitmo: Salati A, Paychere Y, Hahaut V, Gobet C, Naef F. Inferring circadian phases and quantifying biological desynchrony across single-cell transcriptomes. bioRxiv 2026. doi:10.64898/2026.03.30.715278.

# Introduction references: cancer interaction and phase inference

Metadata via Crossref/PubMed; claims from abstracts read unless marked. OA column: PMC id exists (may be embargoed/author manuscript) or bioRxiv.

## Section 1. Cancer interaction

### A1 clock-cell cycle coupling

| key | claim | OA | status |
|---|---|---|---|
| matsuo2003control | Regenerating mouse liver: clock controls cell-cycle genes (Wee1 directly) and Cyclin B1-Cdc2 activity; in single cells the clock runs independently of the cell cycle, i.e. unidirectional clock-to-cell-cycle control. | not confirmed | verified |
| bieler2014robust | Single-cell time-lapse imaging of NIH3T3 with clock (Rev-erbα-YFP) and cell-cycle reporters: bidirectional coupling gives 1:1 phase locking of the two oscillators (mouse fibroblasts). | PMC PMC4299496 | verified |
| feillet2014phase | Single cells: clock and cell cycle phase-lock in multiple modes (1:1, 3:2, 2:1 attractors), tuned by period mismatch and coupling (mouse fibroblasts, imaging + model). | PMC PMC4103330 | verified |
| droin2019low | Theory/data: low-dimensional Arnold-tongue description of two coupled oscillators (clock, cell cycle), predicts locking and cell-cycle slowdown from single-cell time series (Naef lab). | PMC PMC7289635 | verified |
| gaucher2018molecular | Review of molecular links between clock and cell cycle (Wee1, p21, cyclins, checkpoints) and cancer relevance. | not confirmed | verified |
| hunt2007riding | Cell preview/review of coupling between clocks and cell cycle ('Riding tandem'). | not confirmed | verified |
| grechezcassiau2008circadian | BMAL1 regulates p21WAF1/CIP1 via REV-ERBα/ROR and controls hepatocyte proliferation (mouse liver, Bmal1-/- cells). | not confirmed | verified |
| nagoshi2004circadian | Single fibroblasts carry cell-autonomous, self-sustained circadian oscillators, passed to daughter cells across division; clock continues through the cell cycle (mouse fibroblasts). | not confirmed | verified |
| paychere2025cophaser | CoPhaser (Paychere, Salati, Gobet, Naef; bioRxiv 2025): VAE with context-dependent periodic manifolds; unsupervised, single-cell scRNA-seq, no timed training data; recovers cell-cycle phase and circadian phase (aorta), cell-type context absorbs constant shifts. Source of the shared-acrophase-with-fitted-amplitude idea. | bioRxiv (open) | verified |

### A2 circadian disruption and tumours in mice / host-tumour clock

| key | claim | OA | status |
|---|---|---|---|
| filipski2002host | Mice with chronic jet-lag/SCN lesion: disrupted rest-activity and circadian coordination accelerate growth of two transplanted tumour types (Glasgow osteosarcoma, pancreatic adenocarcinoma). | not confirmed | verified |
| papagiannakopoulos2016circadian | Kras;p53 lung adenocarcinoma mouse: jet lag and Per2/Bmal1 mutation decrease survival and promote tumour growth/progression; Per2 and Bmal1 loss linked to Myc and metabolic changes. | PMC PMC5367626 | verified |
| kettner2016circadian | Wild-type mice: chronic jet lag alone induces spontaneous hepatocellular carcinoma via NAFLD/steatohepatitis, nuclear-receptor (FXR/CAR) metabolic dysregulation. | PMC PMC5695235 | verified |
| fu2002circadian | mPer2-mutant mice are cancer-prone after gamma radiation; Per2 controls c-Myc, Cyclin D1, Mdm2, Gadd45a. | not confirmed | verified |
| masri2016lung | Lung adenocarcinoma in mice distally rewires hepatic circadian gene expression and metabolism (tumour-to-host clock effects). | PMC PMC5373476 | verified |
| hadadi2020chronic | Chronic circadian disruption (jet lag) in mice promotes breast cancer stemness, immune-microenvironment changes and lung metastasis. | PMC PMC7314789 | verified |
| aiello2020circadian | Murine melanoma: circadian disruption (shifted LD schedule) speeds tumour growth and remodels tumour-associated macrophages (loss/inversion of daily M1/M2 patterns). | PMC PMC7556830 | verified |

### A3 clock dysregulation in human tumours (tumour vs normal, microenvironment)

| key | claim | OA | status |
|---|---|---|---|
| shaashua2020stromal | Stromal (not tumour) Per2 is required for tumour initiation and metastatic colonisation in mice: host-tissue clock acts on tumours. | PMC PMC7573548 | verified |
| shilts2018evidence | Clock correlation distance (CCD): infers clock progression of a group of untimed samples from co-expression of 12 clock genes; applied to human tumours vs normal (TCGA and others): clock gene co-expression is disrupted in many tumour types (cohort-level, bulk, no per-sample time, no multi-compartment design). | PMC PMC5797448 | verified |
| ye2018genomic | Pan-cancer (TCGA, cell lines) landscape of clock-gene alterations (CNV, mutation, expression) and interaction with drug response; clock genes as chronotherapy targets. | PMC PMC6056007 | verified |
| cadenas2014loss | 766 node-negative breast cancer tumours: higher clock gene expression (PER3, RORC independent) associated with longer metastasis-free survival; loss of clock gene co-expression relationships with progression. | PMC PMC4613905 | verified |
| lellupitiyagedon2019circadian | Luciferase reporters in breast cancer lines: oscillations of BMAL1/PER2 persist in low-grade luminal A MCF7 cells but are lost or weak in high-grade lines. | PMC PMC6739049 | verified |
| wu2019pan | Pan-cancer TCGA: negative-arm clock genes (PER1-3, CRY1-2) consistently down in tumours, via copy-number loss; associated with T-cell exhaustion/immune features. | PMC PMC6821711 | verified |
| broadberry2018disrupted | Primary breast tumour vs matched normal tissue from the same patients: tumour cells show disrupted/abnormal molecular clocks in culture, linked to altered tissue mechanics (paired human samples). | PMC PMC6198506 | verified |
| relogio2014ras | Oncogenic Ras deregulates the circadian clock in cancer cells (cell lines, model). | PMC PMC4038477 | verified |

### A4 clock-dependent cancers

| key | claim | OA | status |
|---|---|---|---|
| puram2016core | In vivo RNAi screens in murine AML: Clock and Bmal1 required for AML/leukaemia stem cell growth in vitro and in vivo; clock disruption impairs proliferation and promotes differentiation (human AML lines also dependent). | PMC PMC4826477 | verified |
| dong2019targeting | Patient-derived glioblastoma stem cells (GSC), but not differentiated cells or non-malignant brain cultures, depend on BMAL1/CLOCK; CRY and REV-ERB agonists reduce GSC growth, synergistic in combination (human GSCs, xenografts). | PMC PMC6983300 | verified |
| sulli2018pharmacological | REV-ERB agonists (SR9009/SR9011) are selectively lethal to cancer cells and oncogene-induced senescent cells, independent of the clock (mouse, cell lines). | PMC PMC5924733 | verified |

### A5 clock strength and prognosis

| key | claim | OA | status |
|---|---|---|---|
| li2024tumor | Luminal A breast cancer: biopsies (Manchester, time-stamped) + public data; CYCLOPS-ordered samples show dampened but present rhythms; HIGH clock strength (CYCLOPS magnitude) tumours have reduced 5-y survival and EMT cycling; clock disruption reduces invasion in 3D culture. Note direction: strong clock = worse prognosis. | PMC PMC10873596 | verified |
| vlachou2024timeteller | TimeTeller: supervised, trained on timed samples; single-sample clock phase and a clock-function (dysfunction) metric from 10-ish genes; applied to cancer prognosis (human, bulk). | PMC PMC10931517 | verified |
| vlachou2019timeteller | TimeTeller preprint 2019: first version; tumour clock function and prognosis (bioRxiv). | bioRxiv (open) | verified |
| vasilyev2026robust | Tumour clock function estimated with TimeTeller negatively associated with survival across cancers (bioRxiv 2026; check final claim). | bioRxiv (open) | verified |
| ananthasubramaniam2025rhythm | COFE: unsupervised, time-series-free ordering of samples (matrix factorisation) used on TCGA/CPTAC to profile rhythms in human cancers vs normal; cohort-level, no per-sample time. | PMC PMC12136439 | verified |
| zhang2025inferring | Sample-time-free test of clock disruption from phase differences between clock gene pairs, tissue-specific reference matrices, ssDistance single-sample score; applied to TCGA tumour vs normal and survival (preprint, single author). | bioRxiv (open) | verified |

### A3 clock dysregulation in human tumours (tumour vs normal, microenvironment)

| key | claim | OA | status |
|---|---|---|---|
| winkler2022gene | Gene-expression signatures predict circadian rhythms in oncogenic pathways in snapshot (TCGA-like) data (bioRxiv). | bioRxiv (open) | verified |

### A6 time-of-day dosing, chronotherapy, TME clocks

| key | claim | OA | status |
|---|---|---|---|
| landre2024effect | Study-level meta-analysis, 13 studies, 1663 patients (NSCLC, RCC, melanoma...): early time-of-day ICI infusion associated with longer OS (HR 0.50) and PFS (HR 0.51); retrospective, hypothesis-generating. | PMC PMC10937202 | verified |
| qian2021effect | MEMOIR, single-centre retrospective (299 stage IV melanoma): each additional 20% of ICI infusions after 16:30 h gave OS HR 1.31; propensity-matched (n=73 v 73) evening infusions associated with shorter OS. | not confirmed | verified |
| wang2024circadian | Mouse and human cancers: CD8 T-cell infiltration and function oscillate through leukocyte clocks and endothelial-clock-dependent trafficking; ICB and CAR-T efficacy improve with time-of-day adjustment. | not confirmed | verified |
| fortin2024circadian | Colorectal cancer mouse model + scRNA-seq: epithelial clock disruption raises inflammatory signalling and PD-L1+ MDSC rhythms, and affects ICB efficacy (time-of-day). | PMC PMC11374317 | verified |
| diamantopoulou2022metastatic | Breast cancer patients (blood draws) and mice: circulating tumour cell intravasation peaks during sleep; rest-phase CTCs are far more metastatic. | not confirmed | verified |
| levi2010circadian | Review: circadian timing in cancer treatment, evidence for chronotolerance and chronoefficacy of anticancer drugs. | not confirmed | verified |
| levi2024circadian | Review: sex-specific, personalized chronotherapy; clinical chronotherapy trials and mathematical personalisation. | not confirmed | verified |
| levi1997randomised | RCT (186 patients, previously untreated colorectal metastases): chronomodulated vs constant-rate oxaliplatin/5-FU/folinic acid improved tolerability and objective response. | not confirmed | verified |
| giacchetti2006phase | EORTC 05963 RCT (564 patients, metastatic colorectal cancer): 4-day chronomodulated FLO vs conventional 2-day FOLFOX2: no survival gain overall (hypothesis of +10% survival not met); sex-dependent toxicity/efficacy. | not confirmed | verified |
| ballesta2017systems | Review: systems chronotherapeutics, using models to design circadian-timed drug schedules. | PMC PMC5394920 | verified |
| karaboue2022time | Pilot, metastatic NSCLC patients on nivolumab: efficacy differs by time of day of infusion (small observational cohort). | PMC PMC8870559 | verified |
| innominato2009circadian | Prospective cohort in a phase III trial (metastatic colorectal cancer): disrupted rest-activity rhythm predicts survival and quality of life. | PMC PMC2690615 | verified |

### A7 reviews

| key | claim | OA | status |
|---|---|---|---|
| sulli2019interplay | Review: circadian clock and cancer, targeting the clock for treatment. | PMC PMC7120250 | verified |
| sancar2021clocks | Review: clocks, cancer and chronochemotherapy (DNA repair, UV, drug timing). | not confirmed | verified |
| kinouchi2020metabolic | Review: circadian homeostasis, metabolism and tumorigenesis. | not confirmed | verified |
| mello2026rhythms | Review (2026): clocks, cancer and chronotherapy, immune-tumour clock links. | PMC PMC12867155 | verified |
| lee2021roles | Review: roles of circadian clocks in cancer pathogenesis and therapy. | PMC PMC8568965 | verified |

## Section 2. Technical challenges / phase inference

### B1 supervised single-sample time predictors

| key | claim | coverage / partial-day | OA | status |
|---|---|---|---|---|
| ueda2004molecular | Molecular timetable: timed mouse (liver) time-course defines marker genes; single-time-point expression profile gives body time, supervised, needs timed training, single sample, bulk. | Supervised, timed mouse training series (full-cycle) -> single-time-point profile; needs timed training data; test sample needs no coverage. | PMC PMC509173 | verified |
| hughey2016zeitzeiger | ZeitZeiger: supervised sparse penalised method predicting periodic variable (time) from high-dimensional expression; timed training data; single sample; bulk. | Supervised; needs timed training data spanning the cycle; single-sample prediction afterwards. | PMC PMC4856978 | verified |
| hughey2017machine | Machine-learning (ZeitZeiger) on human blood: compact gene set for single-sample clock phase; supervised, timed training. | Supervised (ZeitZeiger); timed human blood training data; single-sample prediction. | PMC PMC5329904 | verified |
| braun2018universal | TimeSignature: supervised, trained on timed blood series (one study), predicts circadian time of a new sample within ~2 h for most samples, platform-independent without retraining; blood, human. | Supervised; timed training study required; once trained applies to new studies without retraining. | PMC PMC6166804 | verified |
| wittenbrink2018high | BodyTime: supervised; 12 individuals constant routine blood monocytes to pick biomarkers, NanoString, validated on 28 chronotypes; single blood sample; needs timed training data with full-cycle coverage. | Supervised; timed training (constant routine, full cycle); single blood sample. | PMC PMC6118629 | verified |
| laing2017blood | Blood transcriptome biomarkers for human circadian phase: 100 genes, single-sample, supervised, human constant routine/sleep deprivation. | Supervised; timed training (constant routine/sleep deprivation); single blood sample. | PMC PMC5318160 | verified |
| wu2020population | BHTc/epidermal: population-based single-sample molecular clock phase signature from skin (human), supervised with timed samples. | Supervised; timed training; single epidermal sample. | PMC PMC7441562 | verified |
| agostinelli2016what | Deep-learning (supervised) prediction of circadian time from expression (mouse, human) from timed data. | Supervised; timed training. | PMC PMC4908327 | verified |
| duan2024taufisher | tauFisher: supervised, requires small gene set; predicts time from single bulk or pseudobulk sample; timed training data. | Supervised; timed training; single bulk or pseudobulk sample. | PMC PMC11076472 | verified |

### B2 unsupervised ordering (bulk and single-cell) and multi-tissue human inference

| key | claim | coverage / partial-day | OA | status |
|---|---|---|---|---|
| anafi2017cyclops | CYCLOPS: unsupervised cyclic ordering (autoencoder/ elliptical structure) of cohorts without time labels; applied to human liver, lung, GTEx tissues and cancer; needs cohort with full-cycle coverage; bulk, cohort. | UNSUPERVISED cohort ordering; cyclic ordering presupposes samples spread around the cycle. No explicit full-cycle requirement quoted (PMC text not retrievable). Related quotes: PNAS prefrontal paper notes reconstruction 'did not uniformly distribute samples across the circadian cycle' (search-result paraphrase); DCPR paper: CYCLOPS declines on non-uniform data and 'typically requires more than 250 samples'. | PMC PMC5441789 | verified |
| hammarlund2025cyclops | CYCLOPS 2.0 (Hammarlund thesis, Drexel 2025): covariate-aware CYCLOPS; breast cancer and GTEx; thesis, not peer reviewed; unverified journal version. | UNSUPERVISED covariate-aware CYCLOPS; coverage assumption unverified. | not confirmed | unverified |
| hollis2025reconstructed | CYCLOPS-reconstructed cell-type-specific rhythms in human brain linking Alzheimer's pathology and circadian dysfunction (cohort, bulk + deconvolution). |  | PMC PMC12352436 | verified |
| talamanca2023sex | CHIRAL: unsupervised, GTEx 914 donors, 46 tissues; assigns donor circadian phase by integrating tissues; clock transcripts tightly synchronised across body, rhythms depend on age and sex; per-tissue and donor phase, no coupling/amplitude comparison between compartments. | UNSUPERVISED; fits 24 h harmonic regression over clock genes jointly across samples. No explicit requirement stated in the text read; GTEx donor phases 'distributed fairly uniformly along the 24 hour cycle' (bioRxiv v1) so coverage holds in GTEx. DCPR: robust on non-uniform data, but 'relies on the core clock genes of organisms with an intact, healthy circadian system'. | not confirmed | verified |
| palumaa2026local | GTEx v11: CHIRAL local tissue phase (45 tissues, 927 donors) and donor global phase; closest to a within-donor comparison but donor phase = consensus, no coupling/amplitude model. | UNSUPERVISED (CHIRAL per tissue + donor consensus); coverage inherited from CHIRAL. | bioRxiv (open) | verified |
| grabe2024quantification | Snapshot lung data (human GTEx lung): quantification of circadian amplitude/phase from unsorted bulk samples (Herzel lab). |  | PMC PMC11247089 | verified |
| han2025dcpr | DCPR: unsupervised deep-learning circadian phase reconstruction, no time labels; benchmarks non-uniform and incomplete-cycle (<24 h span) sampling: 'When dealing with non-uniform datasets, DCPR and CHIRAL exhibit more robust performance ... than Cyclum and CYCLOPS'; all methods degrade on incomplete cycles; says CYCLOPS typically needs >250 samples; CHIRAL 'relies on the core clock genes of organisms with an intact, healthy circadian system'. Closest method paper on partial-day coverage (simulation + real sets). | UNSUPERVISED; explicitly tested on non-uniform and incomplete (<24 h) cycles in simulation/real data: all methods degrade; DCPR best. | PMC PMC12866578 | verified |
| xue2026bayesian | BayCT (Xue ... Tseng, bioRxiv 2026): Bayesian rhythmic model, von Mises prior, predicts molecular circadian time in post-mortem brain (59 subjects, time of death as input) and mouse atlas; frames mismatch between recorded time of death and true molecular time; no statement on window coverage. |  | bioRxiv (open) | verified |

### B3 human tissue rhythm atlases / time-of-death studies

| key | claim | OA | status |
|---|---|---|---|
| ruben2018database | CYCLOPS applied to 13 GTEx tissues (632 donors) without time-of-day: database of tissue-specific rhythmic human genes (unsupervised, cohort-level, bulk). | PMC PMC8961342 | verified |
| chen2015effects | Time-of-death analysis, human prefrontal cortex (n=146, microarray): ~10% transcripts rhythmic; aging shifts and reduces rhythms. | PMC PMC4711850 | verified |
| li2013circadian | Time-of-death analysis in six postmortem cortical/limbic regions (55 controls): region-specific rhythms, dysregulated in major depression. | PMC PMC3683716 | verified |
| hughey2016differential | Five public studies (>1000 samples): clock gene expression in human brain phase-delayed ~8.5 h relative to blood; compares phases between organs across cohorts, not within individuals. | PMC PMC5105327 | verified |
| zhang2014circadian | Mouse circadian atlas: 43% of protein-coding genes rhythmic in at least one of 12 tissues; tissue-specific rhythms. | PMC PMC4234565 | verified |
| mure2018diurnal | Baboon diurnal atlas of 64 tissues/ brain regions; time-resolved tissue clock phase relationships (Naef lab, with Mure). | PMC PMC5924732 | verified |
| wucher2023day | GTEx (932 donors): uses recorded time and season of death to identify day-night and seasonal variation across tissues; tissue clock genes vary day-night. | PMC PMC9934459 | verified |

### B4 single-cell phase inference and limits

| key | claim | coverage / partial-day | OA | status |
|---|---|---|---|---|
| auerbach2022tempo | Tempo: unsupervised Bayesian single-cell circadian phase inference (needs clock gene prior), single-cell. | UNSUPERVISED single-cell; cells assumed to span cycle (von Mises phase prior); no explicit statement verified. | PMC PMC9630322 | verified |
| leng2015oscope | Oscope: unsupervised identification of oscillatory genes and cell ordering in unsynchronised scRNA-seq; cell cycle. | UNSUPERVISED; ordering assumes cells cover the cycle of an oscillation; no explicit statement verified. | PMC PMC4589503 | verified |
| liang2020latent | Cyclum: autoencoder with circular latent for periodic process (cell cycle) inference from scRNA-seq; unsupervised. | UNSUPERVISED; circular latent; assumes cells cover the cycle; not verified. | PMC PMC7080821 | verified |
| riba2022cell | DeepCycle: deep learning with RNA velocity to estimate cell-cycle phase from scRNA-seq, unsupervised. | UNSUPERVISED (RNA velocity + DL); cells assumed to span cycle; not verified. | PMC PMC9126911 | verified |
| schwabe2020transcriptome | Cell cycle transcriptome dynamics of single cells (scRNA-seq), cell-cycle ordering from unspliced/spliced counts; Schwabe et al. |  | PMC PMC7672610 | verified |
| nikhat2026transcriptional | Information-theoretic limit: noise limits single-cell clock phase decoding from scRNA-seq; justifies pseudobulking (iScience 2026). |  | PMC PMC13084431 | verified |
| salati2026inferring | scRitmo: unsupervised Bayesian single-cell circadian phase from mRNA counts with posterior uncertainty and desynchrony quantification (Naef lab, preprint). | UNSUPERVISED single cell; not verified. | bioRxiv (open) | verified |
| phillips2020circadian | Naef lab preprint: smFISH of Bmal1, Cry1, Nr1d1 in single fibroblasts: 10-35 mRNA per cell oscillate; count noise and cell-size heterogeneity make single-cell phase from mRNA counts intrinsically imprecise. |  | bioRxiv (open) | verified |

### B5 rhythmicity detection

| key | claim | OA | status |
|---|---|---|---|
| hughes2010jtk | JTK_CYCLE: nonparametric rhythm detection in timed time series. | PMC PMC3119870 | verified |
| thaben2014detecting | RAIN: rank-based nonparametric detection of rhythms of arbitrary waveform in timed data. | PMC PMC4266694 | verified |
| wu2016metacycle | MetaCycle: integrates JTK, LS, ARS for periodicity detection in timed series. | PMC PMC5079475 | verified |
| thaben2016differential | Differential rhythmicity method (limma-based) for comparing rhythms between conditions in timed data. | not confirmed | verified |
| weger2021systematic | dryR: Bayesian model selection to classify rhythmic/differentially rhythmic genes (Naef lab), requires timed series. | PMC PMC7826335 | verified |
| pelikan2021venn | Venn-type comparison overestimates rhythm reprogramming; recommends model-selection (compareRhythms) (Herzel/Ananthasubramaniam). | not confirmed | verified |
| atger2015circadian | Mouse liver: circadian vs feeding rhythms in transcription and translation; dryR-type analysis precursor (Naef lab). | PMC PMC4664316 | verified |

### B6 human sampling constraints

| key | claim | OA | status |
|---|---|---|---|
| lamia2024telling | PNAS commentary on Li et al. 2024: 'Sample collection time is rarely recorded when biopsies are performed'; collection time inferred from high-amplitude genes. | PMC PMC10945801 | verified |
| bivolarski2026circadian | Review (prostate cancer): 'tissue samples are collected according to institutional schedules rather than circadian considerations'; most clinical transcriptomic datasets lack collection time (secondary source; not a measured time distribution). | PMC PMC13423029 | verified |
| gonzalezaponte2025circadian | Glioblastoma, 302 patients operated 2020-2024: probability that a biopsy is called MGMT-methylated depends on hour of biopsy (peak around noon, minimum ~16:00); example that time-of-day of surgical sampling confounds tumour biomarkers. | PMC PMC12546296 | verified |
| thegtexconsortium2020gtex | GTEx v8: 838 donors, 15,201 RNA-seq samples, 49 tissues; time of death recorded (controlled access). | PMC PMC7737656 | verified |
| lonsdale2013genotype | GTEx pilot description: post-mortem multi-tissue collection from donors, sampling after death. | PMC PMC4010069 | verified |

## Evidence on sampling-time clustering (requested by lead; read before using)

- GTEx time of death is NOT clustered in working hours. wucher2023day: "since this time is approximately uniformly distributed throughout the day and the year" (Discussion); talamanca2023sex (bioRxiv v1): inferred donor phases "are distributed fairly uniformly along the 24 hour cycle (Table ST1), as the TODs". So for GTEx, full-cycle coverage of the cohort holds; the clocksync argument for GTEx must rest on single-donor/within-donor inference and on small-n tissues, not on a working-hours bias. Hospital-hours clustering applies to surgery/biopsy cohorts (tumour vs matched normal).
- No paper was found that measures the time-of-day distribution of surgical resections or biopsies. Indirect citable support: lamia2024telling ("Sample collection time is rarely recorded when biopsies are performed"); bivolarski2026 review in B6 (samples "collected according to institutional schedules rather than circadian considerations"; secondary source, weak); gonzalezaponte2025 (302 GBM biopsies, MGMT call depends on hour of biopsy, so surgical timing varies and matters); li2024tumor (time-stamped Manchester breast biopsies, supplement likely has the hour distribution: not read, PMC blocked). Recommendation: cite the Li 2024 time-stamped biopsy cohort and, if needed, compute the hour histogram ourselves from GTEx-independent data with recorded times; flag as unverified until the Li 2024 supplement is checked.

## Close to our gap (cross-compartment, untimed, same individual / partial-day coverage)

1. talamanca2023sex (CHIRAL on GTEx, Science 2023): donor internal phase integrates temporal information from several tissues of one donor; clock transcripts "tight synchrony across the body". Unsupervised, untimed, multi-tissue per donor. Gap: yields one consensus phase per donor and per-tissue phase; no per-compartment coupling/relative amplitude/offset model, no disease contrast, assumes cohort with spread phases.
2. palumaa2026local (bioRxiv 2026, GTEx v11, 45 tissues, 927 donors): CHIRAL local phase per tissue sample vs donor global phase; concluded that donor phase captures tissue clocks only partly. Closest competitor; same gap as above; check for scooping before submission.
3. ruben2018database (CYCLOPS on 13 GTEx tissues) and hughey2016differential (brain vs blood ~8.5 h lag, pooled cohorts, not within donor).
4. zhang2025inferring (bioRxiv 2025, single author): sample-time-free clock disruption via clock-gene phase differences, ssDistance single-sample score, tumour vs normal TCGA. No time needed at test (like ours) but within-sample gene phase relations, not cross-compartment within person.
5. shilts2018evidence (CCD, no sample times, tumour vs normal), ananthasubramaniam2025rhythm (COFE, tumour vs normal, time-series-free), li2024tumor + hammarlund thesis (CYCLOPS in tumours): cohort-level ordering within tissue, not within person across tissues.
6. han2025dcpr (DCPR): the only method paper found that explicitly benchmarks non-uniform and <24 h-span (incomplete cycle) sampling; all methods degrade, CHIRAL second best on non-uniform. Quote available in the claim cell. No method found that handles restricted sampling windows by design; clocksync's "constant offsets cancel" is, to our search, unmatched.
7. xue2026bayesian (BayCT, bioRxiv 2026): Bayesian molecular circadian time in post-mortem brain; addresses recorded-vs-molecular time mismatch, not coverage or cross-tissue.
8. grabe2024quantification (snapshot lung data) and hughey2016zeitzeiger/BodyTime family: single-compartment.
Searches for "within-person multi-tissue untimed clock comparison" found nothing beyond 1-2; no coupling/amplitude-relative-to-host method found.

## Considered and excluded

- Huang et al. Nat Med 2026 (10.1038/s41591-025-04181-w, LungTIME-C01 RCT, early vs late immunochemotherapy): RETRACTED (retraction note 10.1038/s41591-026-04508-1). Do not cite as evidence; coverage in news and Cancer Discovery highlights is stale.
- Filipski 2004 Cancer Res (10.1158/0008-5472.CAN-04-0674) and Filipski 2005 JNCI: other Filipski papers; kept only the 2002 JNCI paper as asked (the Crossref record for 2002 lists the first author only; full author list from secondary records).
- Hammarlund 2023 bioRxiv (10.1101/2023.05.17.540386): superseded by Li 2024 PNAS.
- Qian "Lancet Oncol 2021" confirmed as MEMOIR (retrospective, not RCT). Candidate "Landre 2024 immunotherapy infusion time" is a study-level meta-analysis (13 studies, 1663 patients), not an original cohort.
- Candidate "Ye et al. 2018 Cell Syst" title corrected: "The Genomic Landscape and Pharmacogenomic Interactions of Clock Genes in Cancer Chronotherapy". Candidate "Lellupitiyage Don 2019 Cell Cycle" title corrected: "Circadian oscillations persist in low malignancy breast cancer cells". Candidate "Dong 2019" title: "Targeting Glioblastoma Stem Cells through Disruption of the Circadian Clock". Sancar & Van Gelder title: "Clocks, cancer, and chronochemotherapy" (not chronotherapy).
- Candidate "Wucher et al. ... Naef": the Wucher paper is Guigo/Irimia lab, not Naef lab (Wucher, Sodaei, Amador, Irimia, Guigo, PLoS Biol 2023). CYCLOPS 2.0: thesis only (Drexel 2025); no journal/preprint found. "Transcriptional noise ..." is Nikhat, Mandal, Veerasubramanian, Chakrabarti (iScience 2026, 10.1016/j.isci.2026.115394), not Naef lab. "PRECISE": no such paper found. "Hughey limorhyde/deltaCCD": limorhyde (Singer & Hughey 2019 JBR) not retrieved; deltaCCD is covered by Shilts 2018 (CCD); add limorhyde if wanted (unverified).
- Naef lab clock-cell cycle work beyond Bieler/Feillet/Droin: none found in 2022-2026 beyond CoPhaser (somite clock-cell cycle coupling) and phillips2020circadian (single-transcript). Gutu 2024 "Circadian Coupling Orchestrates Cell Growth" bioRxiv is not Naef lab (not included).
- CircaDB, Hanahan/Weinberg hallmarks, MetaCycle-like extras beyond list: not needed.
