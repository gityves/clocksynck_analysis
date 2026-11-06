# Literature synthesis

Built from the cards in `notes/` (full texts of the Zotero export, quotes checked by `check_quotes.py`). Numbers here
are taken from the cards' verbatim quotes; follow the key to the card for the quote and its location.
Status: all 148 full texts carded (144 + 4), 865+ quotes verified; claims matrix to be built with the introduction.

## 1. Hospital hours: what the literature states

| statement | source (card) | strength |
|---|---|---|
| Orderers need the whole cycle: "It requires data from the entire periodic cycle to form an ellipse. Biopsies are almost exclusively obtained during the day." | anafi2017cyclops | author's own limitation, verbatim |
| "Most computational methods for quantifying circadian rhythmicity and inferring the status of the clock require that samples be acquired over the entire 24-h cycle." | shilts2018evidence | verbatim |
| Surgical biopsies in working hours (~6 am-6 pm) | anafi2017cyclops, shilts2018evidence, grabe2024quantification ("clinical samples are collected during surgical working hours that are not evenly spread across the day") | stated, not measured |
| TCGA inferred phases "concentrated during the presumed working day"; GTEx more uniform but biased | li2024tumor | inferred phases (CYCLOPS 2.0), not recorded times |
| GTEx deaths cover the day: 351 (08-17 h), 315 (21-05 h), 266 twilight donors ≈ 38 per hour in each window | wucher2023day | published aggregate (recorded TOD) |
| Collection time usually unrecorded; most dosing-time trials compare two time points | ruben2019dosing | perspective |
| Recorded surgery times (end of anaesthesia) of 302 GBM patients analysed in 2 h bins "from 10am to 6pm"; authors: "we do not know what happens to MGMT promoter methylation from around 6pm until 10am" | gonzalezaponte2025circadian | only recorded hospital-hours distribution in our set |
| Recorded resection times (Manchester breast, REMAGUS breast) exist only for small subsets | li2024tumor, vasilyev2026robust | data |
| Orderers degrade on incomplete cycles (14-18 h synthetic spans) | han2025dcpr | only benchmark of partial coverage; numbers in SI |

**Consequence for the intro:** two published claims of tumour-clock "uncoupling" come from methods that need full
coverage applied to working-hours cohorts: COFE (TCGA tumour times span the whole cycle, tumour-normal differences
random; authors note daytime surgery makes this "unlikely") and TimeTeller (93.2 % of tumour times in the working
day only after reassigning off-window samples; r ≈ 0 with sampling time within a narrow window, n = 108). Li 2024
sees poor tumour-normal phase agreement except in strong-clock tumours, compatible with noise in weak clocks.
Our A138 finds a tumour-normal clock match in TCGA pairs; A145 W arm: COFE coupling MAE 0.206 vs clocksync 0.060.

## 2. Methods table

| method (card) | unit | supervision | needs 24 h coverage | gene phases | within-person compartments | noise model |
|---|---|---|---|---|---|---|
| CYCLOPS (anafi2017cyclops) | sample, per tissue cohort | none (seed genes) | yes (stated); >250 biopsies | learned (ellipse) | no | none |
| CYCLOPS 2.0 (li2024tumor) | sample | optional timed seeds + covariates | benchmarked under temporal bias (SI) | learned | no (tumour projected on normal) | none |
| CHIRAL (talamanca2023sexdimorphic) | sample → donor (DIP) | none (12 clock genes) | yes (DIPs ~uniform in GTEx) | fixed-ish seed set | assumes common tissue shift; no per-donor scatter | Gaussian |
| Palumaa LP (palumaa2026local) | tissue sample | CHIRAL TIP step | yes | CHIRAL | LP-GP r 0.20-0.70, noise not separated | — |
| COFE (ananthasubramaniam2025rhythm) | sample | none, de novo features | yes in practice (CV assumes spread) | learned | separately ordered tumour/normal | — |
| DCPR (han2025dcpr) | sample | none | degrades on 14-18 h spans | learned | no | — |
| Cyclum (liang2020latent) | cell (cell cycle) | none | — | learned | no | — |
| CCD (shilts2018evidence) | group | none | no (works on 8 h windows) | reference correlations | no; "cannot detect a phase difference" | — |
| TimeTeller (vlachou2024timeteller, vasilyev2026robust) | sample | trained on timed healthy series | no, but phase restricted by likelihood | fixed by training | no | Gaussian on gene vector |
| ZeitZeiger / TimeSignature / BodyTime / PLSR (hughey2016zeitzeiger, braun2018universal, wittenbrink2018high, laing2017blood) | sample | timed training, blood/monocytes (ZeitZeiger: mouse) | no (supervised) | trained | no | — |
| CYCLOPS→ZeitZeiger skin (wu2018populationlevel, wu2020population) | sample | hybrid | CYCLOPS step yes | learned | no (epidermis/dermis separate) | — |
| Tempo (auerbach2022tempo) | cell | acrophase prior | not stated | learned with prior | no | NB |
| scRitmo (salati2026inferring) | cell | none (fixed reference) | no (cells within sample) | fixed | no (cell scale σbio) | NB, technical σ simulated |
| CoPhaser (paychere2025cophaser) | cell | none | cells must span cycle | shared, context mean/amplitude | no donor level | NB |
| BayCT-M (xue2026bayesian) | donor shift shared by 3 brain regions | TOD prior | needs TOD | learned | one shared shift, no per-compartment coupling | — |
| CYCLOPS 2.0 on snRNA brain (hollis2025reconstructed) | donor (one phase from excitatory neurons, applied to all cell types) | none | yes | learned | no (cell types share the donor phase by construction) | — |
| ORPHEUS (hollis2026singlecell) | cells within a cell type: gene-wise phase dispersion | needs subject phases (time course or CYCLOPS) | yes | — | no | — |
| Grabe 2024 (grabe2024quantification) | sample; tumour ordered by its paired normal's inferred time | none on target (gene list from timed mouse/baboon) | assumed ("facilitated, if a dataset covers a full circadian cycle"); partial coverage untested | reference correlations to PER3 | tumour vs own normal (closest host-referencing precedent), no coupling estimate, no gate, no CI | — |
| Zhang 2025 (zhang2025inferring) | phase differences from co-expression | reference matrices | no | reference | no | — |
| **clocksync** | donor × compartment | none (host-referenced) | **no** (constant offsets cancel) | fixed acrophases shared (CoPhaser idea) | **yes: c, abs(z), arg z per compartment** | Poisson V_ik + compartment σ² |

## 3. Corrections to `refs_*.md` (to apply)
- anafi2017cyclops: full-cycle requirement IS explicit (quote above).
- palumaa2026local: not a competitor (CHIRAL TIP step + catalogue).
- stokkan2001entrainment: rats, not mice. buhr2010temperature: explants only. brown2002rhythms: liver only in vivo.
- kume1999mcry: Cry double KO not tested (van der Horst 1999, missing). preitner2002orphan: REV-ERBα dispensable for rhythm generation.
- mure2018diurnal: Panda/Cooper, not Naef. mermet2018clockdependent: Rev-erb not perturbed.
- kervezee2019individual: plasma only (not two compartments). wefers2018circadian: muscle clock did not shift. wu2018populationlevel: 219 biopsies untimed, CYCLOPS-ordered.
- vlachou2024timeteller: no prognosis analysis. vasilyev2026robust: breast only, nonlinear Θ-survival. winkler2022gene: PROGENy on timed mouse liver, not TCGA phase inference.
- yoo2004periodluciferase: headline is persistence + tissue-specific period/phase, not damping.
- hughey2016zeitzeiger: mouse only.
- reddy2002differential: no SCN shell/core; mPer vs mCry1 dissociation within SCN. kiessling2010adrenal: adrenal clock ablation accelerates re-entrainment (refs say GC needed). vollmers2009time: Cry1/Cry2 double KO, not Bmal1.
- tahara2012in: no jet lag (SCN lesion + restricted feeding). kolbe2019circadian: no misalignment design. archer2014mistimed: loss of rhythmicity (6.4 % → 1.0 % of transcripts), not phase shifts. sinturel2021circadian: Schibler/Dibner, not Naef.
- kervezee2018simulated (PBMC only) and skene2018separation (plasma, between subjects): not two compartments in the same people.
- feillet2014phase: not Naef lab (Delaunay, Rand, van der Horst); modes 1:1, ~3:2, ~5:4. bieler2014robust: coupling predominantly cell cycle → clock. droin2019low: cell cycle accelerates the clock (not a cell-cycle slowdown). matsuo2003control: gating shown, unidirectionality not.
- filipski2002host: SCN lesion only (no jet lag). masri2016lung: host core clock unchanged; output rewiring (46 % of 505 shared rhythmic liver genes shift ≥ 1 h).
- wu2019pan: loss of rhythm inferred from mean mRNA only (cite for copy-number down-regulation). shaashua2020stromal: global host Per knockouts, not stroma-specific. broadberry2018disrupted: rhythm in 3 patients, qPCR in 8.
- sulli2018pharmacological: clock-independence of SR9009/SR9011 lethality not shown. fortin2024circadian: epithelial clock loss abolishes (not raises) the MDSC rhythm; scRNA at one time point.
- schwabe2020transcriptome (Revelio): PCA rotation, not spliced/unspliced (that is riba2022cell). thaben2016differential (DODR): not limma-based. pelikan2021venn: hypothesis testing preferred over model selection.
- grabe2024quantification: human data are TCGA-LUAD paired normal/tumour, not GTEx lung. agostinelli2016what: mouse only. mello2026rhythms: immune-tumour content minor.
- thegtexconsortium2020gtex: no time of death, Hardy scale or ventilation in the paper (dbGaP only); cite for cohort sizes. GTEx donors are post-mortem/organ donors without exclusion of common diseases (lonsdale2013genotype): write "non-cancer reference", not "healthy".
- hollis2025reconstructed: snRNA-seq pseudobulk per cell type, not bulk deconvolution. galinde2023different: does not state a lack of human biomarkers (cite for the 4-level desynchrony taxonomy only).

## 4. Verification ideas (checked against METHOD.md / registered / results; "new" = not found there)
- done (A204, descriptive; A204b): COFE's own outputs under simulated 08-18 h sampling (full-cycle labels, random tumour-normal differences) to show the mechanism of false uncoupling (A145 has the coupling MAE, not the mechanism).
- new: Wu 2020 skin (154 women × arm/buttock/cheek × epidermis/dermis, untimed): within-person multi-site compartments.
- new: GTEx brain arg z vs periphery: discriminates Hughey 2016 (~8.5 h, between people) from CHIRAL (few hours).
- new: ΔCCD next to clocksync coupling on TCGA pairs (A138): do low-coupling tumours have high CCD?
- new: stress-test the gate with gene-specific partial shifts inside one tissue (Kiessling 2010, Saini 2013: Bmal1-luc 6-7 h vs Rev-erbα/PER2 ~12 h after feeding inversion); METHOD.md has scrambled phases only.
- new: Cry1/2 DKO under restricted feeding (Vollmers 2009, GSE13093): feeding-rhythmic clockless liver as a gate negative control.
- new: adipose vs liver offsets as a function of feeding-cue strength, prior from Heyde 2022 (day-4 delays eWAT 1.6 vs 7.1 h, liver 2.6 vs 5.8 h across feeding regimes).
- planned, not registered: A117 timed brain regions (provisory_plan_paper.md WP7): GSE45642 (Li 2013: cross-region consistency lost in MDD) and GSE71620 (Chen 2016: age). Strongest literature precedent for disease-related loss of cross-region synchrony; register before scoring.
- new: organ phase offset arg z vs proliferation (MKI67, MYC) across tissues of the same GTEx donor, with a specificity gene set (Droin 2019: 12 mouse tissues span 1.5 h, offsets correlate with Mki67/Myc; cell cycle accelerates the clock). A97/A126/A127/A166 test coherence only.
- new: Masri 2016 GSE73222 (timed liver, lung-tumour-bearing vs WT): core clock offset ≈ 0 under a distant tumour while outputs shift; a specificity positive control for the clock panel.
- planned, not registered: A125 mouse chronic jet lag (Kettner 2016, GSE75475; provisory_plan_paper.md).
- done: Cadenas 2014 replication (A192, mostly proliferation); Ye 2018 co-expression loss in TCGA pairs (A138 H138.4); copy number and malignant coherence (A134).
- new: ER-stratified malignant coupling at common d* with matched non-clock pairs (partly A134/A138.5); stromal |z| tumour-adjacent vs normal stroma of the same donors (Broadberry 2018 stiffness).
- new: sex-stratified coupling c and offset arg z at common d* (only the sex phase A99 exists), motivated by the sex × schedule interaction of chronotherapy (levi2010circadian, levi2024circadian).
- framing: map the three published tumour scenarios (arrhythmic / out of phase / in phase; sancar2021clocks) onto |z| ≈ 0, low c with offset, high c.
- done: chronotherapy benefit vs clock readout (A198, null), tumour clock vs survival (A191-A194), T2D/obesity (A164).
- new: glioblastoma single-cell cohorts (malignant vs non-malignant brain per patient; Dong 2019 GSC clock dependence); no glioma cohort in the docs (only paragangliomas).
- new: Fortin 2024 Apc/epithelial-Bmal1 mouse (epithelial clock off, spleen intact) as a controlled cross-compartment desynchrony positive, if timed data are public.
- motivation: Diamantopoulou 2022 (CTC timing driven by host hormones, tumour clock genes not differential) supports host-referenced readout; parallels A138/A149.
- methods support: pelikan2021venn (Venn comparisons of thresholded lists: 48 % false differential hits) justifies direct permutation contrasts over counting gate-passing units; phillips2020circadian and nikhat2026transcriptional justify pseudobulks (24 h explains ~all variance when pooling ~1,000 cells; 20-70 cells × 3 genes ≈ 1 h): compare with A203 (registered, unscored) when scoring.
- new: Grabe 2024 procedure (PER3-correlation shifts; tumour ordered by the paired normal's time) on TCGA-LUAD vs clocksync coupling (A171 covers LUAD coherence only).
- done: tauFisher timed dermis scRNA (GSE223109) is the V13 cell-type truth set (`common/datasets/skin_dermis_zt`); the earlier "new" label was wrong (checked 2026-11-06).
- new: gate vs technical cyclic artefacts (Oscope C1 chip position, leng2015oscope).
- done: ORPHEUS bias under subject-phase error (V14); RADC/brain AD contrasts (A135, A141, mostly not estimable).
- done: baboon (A36/A40/A41), living skin epidermis/dermis (A37, Del Olmo), forced desynchrony blood (GSE48113), TimeTeller/COFE/CYCLOPS 2.0/CHIRAL benchmark (V2), Manchester breast pairs (A25 V25.5, A146), orderers under working hours (A145, A123).

## 5. Missing literature
- Giacchetti 2012 (sex and chronotherapy, 842 patients), Mormont 2000, Innominato 2014 (rest-activity and survival), Dallmann 2016 (dosing time).
- Fuhr 2018/2019, Gutierrez-Monreal 2016 (tumour clocks).
- Filipski 2005, Lee 2010 (Cancer Cell), Hojo 2017 (tumour reprograms hepatic circadian transcriptome), Gery 2006, Altman 2015 (MYC and clock).
- Husse 2014, Izumo 2014 (SCN-specific Bmal1 loss); Leproult 2014, Morris 2015 (human misalignment and glucose); Janich 2011 (epidermal clock).
- van der Horst 1999 (Cry1/Cry2 double KO arrhythmia); Kornmann 2007 (liver-specific clock); Cho 2012 (REV-ERB double KO).
- Giacchetti 2012 (sex and chronotherapy survival); Iurisci 2006 (tumour clocks).
- hughey2017machine, hastings2018generation: in bib, PDFs missing from the export.
