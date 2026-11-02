# Literature positioning

_Status: references checked by web search on 2026-11-01 (PubMed connector unavailable); entries marked candidate or unverified must be checked before citing._

Compiled 2026-11-01 from the table "Literature positioning to do before writing" in `HANDOFF_2026-10-31.md`.
Every citation below was checked against a PubMed Central, publisher, bioRxiv or EPFL record in the search session
(one record per reference). PMIDs are given only where a record stated them. Status codes: **verified**,
**corrected** (the handoff description was wrong or vague; the correct citation is given), **candidate**
(closest match found, the handoff reference itself not pinned down), **not found**.

Our results referred to (source: `results_preliminary.md`, A174-A181 and STATE section 3):
cancer coherence loss in proliferating and ER-negative tumours is a swamped clock, not disordered relations
(A176, 7 of 7 contrasts, in-plane log-odds falls, off-plane variance rise not clock-specific); C1 host-pairing is
lower in the high-proliferation tertile in TCGA (A175, secondary because TCGA is already seen; 0.459 vs 0.301,
D +0.158, p 0.021, clock-specific against pseudo-clock and responsive sets); the coherence loss is not explained by
linear cell-cycle composition in bulk (A179, 0 of 5, resolved in breast and lung, unresolved in gastric); the fetal
gate failure is uninformative and design-limited (A181).

## 1. Untimed clock analysis in cancer

- Shilts J, Chen G, Hughey JJ. Evidence for widespread dysregulation of circadian clock progression in human cancer.
  *PeerJ* 2018;6:e4327. doi:10.7717/peerj.4327. PMID 29404219. **Verified.**
  The clock correlation distance (CCD) compares clock-gene co-expression of a sample group with a mouse reference.
  Tumour CCD exceeds non-tumour CCD in 12 TCGA types (paired t-test P = 1.4e-8) and in independent datasets; no clear
  trend with grade; the change in CCD correlates with median tumour purity (rho 0.67, P = 0.059, not significant).

How we differ. CCD is a group-level co-expression distance to a fixed mouse reference; clocksync reads each
compartment against the clock state shared with the same donor's other compartments, so tumour and its matched normal
tissue define the reference and no external reference is needed. We add what CCD does not report (to our reading of the
paper; proliferation was not checked in the full text): a specificity companion on the same permutations
(pseudo-clock and responsive gene sets), purity strata, proliferation and ER as predictors, external replication
cohorts and a separation of "less clock signal" from "disordered relations". A176 shows the loss is a swamped
clock (in-plane variance falls, off-plane rise is not clock-specific), and A175 shows lower pairing with the normal
tissue in the high-proliferation tertile in TCGA (secondary; primary replication pending). CCD, being a correlation
property, cannot make that distinction.

## 2. Cell cycle - clock coupling

- Bieler J, Cannavo R, Gustafson K, Gobet C, Gatfield D, Naef F. Robust synchronization of coupled circadian and cell
  cycle oscillators in single mammalian cells. *Mol Syst Biol* 2014;10:739. doi:10.15252/msb.20145218. PMID 25028488.
  **Verified.** NIH3T3 time-lapse: 1:1 mode locking, divisions about 5 h before the Rev-erbα-YFP peak, reverse
  coupling (cell cycle to clock) predominant; the authors note implications for proliferative tissues and cancer.
- Feillet C, Krusche P, Tamanini F, Janssens RCC, Downey MJ, Martin P, Teboul M, Saito S, Lévi FA, Bretschneider T,
  van der Horst GTJ, Delaunay F, Rand DA. Phase locking and multiple oscillating attractors for the coupled
  mammalian clock and cell cycle. *Proc Natl Acad Sci USA* 2014;111(27):9828-9833. doi:10.1073/pnas.1320474111.
  **Verified** (title, journal, pages, DOI confirmed; the author list is from memory and the sources disagreed on
  initials, check on the PNAS page). Fibroblasts show 1:1 phase locking and coexisting higher-period states.
- Droin C, Paquet ER, Naef F. Low-dimensional dynamics of two coupled biological oscillators. *Nat Phys*
  2019;15:1086-1094. doi:10.1038/s41567-019-0598-1. **Verified.** Inferred coupling from mouse and human single-cell
  traces predicts several locked states and explains why tissues with different proliferation rates have shifted
  circadian phases.

How we differ. These are single-cell, mostly fibroblast, imaging studies; none measures an inter-compartment
coherence in patient material. They make locking a prediction rather than a confound: proliferation can shift or
reshape the clock state in a cell-type-dependent way. In bulk, however, the loss of coherence in proliferating
tumours is not accounted for by a linear G1/S and G2M composition term (A179: 0 of 5; the loss grows in SCAN-B,
GSE72094 and GSE68465), and A175 shows the proliferation effect on host pairing is clock-specific. The C4
single-cell test of locking is untested, so locking stays a prediction, not a finding.

## 3. Breast subtype clocks

- Lellupitiyage Don SS, Lin HH, Furtado JJ, Qraitem M, Taylor SR, Farkas ME. Circadian oscillations persist in low
  malignancy breast cancer cells. *Cell Cycle* 2019;18(19):2447-2453. doi:10.1080/15384101.2019.1648957.
  PMID 31357909. **Corrected.** The handoff's "Lin & Farkas review around 2018" is this primary study (not a review,
  2019, first author Lellupitiyage Don; Lin is second author): luminometry shows BMAL1/PER2 oscillations in luminal A
  MCF7 but not in basal MDA-MB-231 cells, supporting a link between grade and clock disruption.
- Candidate for a broader basis: "Circadian clock features define novel subtypes among breast cancer cells and shape
  drug sensitivity" (bioRxiv 10.1101/2024.10.07.617001; PMC11965565; Berlin group, apparently Ector, Didier and
  colleagues; author list and journal version not checked). Screens 19 cell lines (14 breast models) and reports four
  clock phenotypes (functional, weak, unstable, dysfunctional). **Candidate; verify before citing.**

How we differ. Cell-line evidence compares one luminal and one basal line, or a screen of lines, and cannot separate
lineage from proliferation or from the tumour microenvironment. In SCAN-B and METABRIC the coherence loss in
ER-negative tumours is swamped-clock, with the in-plane variance fall not clearly clock-specific in the ER contrasts
(A176; SCAN-B -0.087, p 0.22) and weaker without PER1/NR1D1. BRCA ER pairing contrast is not estimable at this n
(A175, MDE 0.408). Acute ER blockade does not reproduce the loss (A168/A169), so we frame ER-negative loss as a
lineage-state property, not an acute estrogen-signalling effect.

## 4. Untimed organ ordering in GTEx

- Anafi RC, Francey LJ, Hogenesch JB, Kim J. CYCLOPS reveals human transcriptional rhythms in health and disease.
  *Proc Natl Acad Sci USA* 2017;114(20):5312-5317. doi:10.1073/pnas.1619320114. **Corrected attribution:** this, not
  Ruben 2018, is the CYCLOPS method paper (circular node autoencoder, human liver and lung).
- Ruben MD, Wu G, Smith DF, Schmidt RE, Francey LJ, Lee YY, Anafi RC, Hogenesch JB. A database of tissue-specific
  rhythmically expressed human genes has potential applications in circadian medicine. *Sci Transl Med*
  2018;10(458):eaat8806. doi:10.1126/scitranslmed.aat8806. **Corrected** (title and authors; the handoff gave no
  title). Applies CYCLOPS to GTEx (preprint: 632 donors, 13 tissues): about 7,500 cycling genes, 917 with a drug
  target, transporter or metabolising function, core clock shared and most rhythmic expression tissue-specific.
- Talamanca L, Gobet C, Naef F. Sex-dimorphic and age-dependent organization of 24-hour gene expression rhythms in
  humans. *Science* 2023;379:478-483. doi:10.1126/science.add0846. **Verified.** CHIRAL assigns a donor internal
  phase from multiple tissues in GTEx (914 donors, 46 tissues); rhythms are more sustained in females and dampen with
  age.

How we differ. These orderers need coverage of the full day, which GTEx lacks (deaths cluster in working hours).
In our test every untimed orderer fails on 08-18 h donors (CHIRAL -0.165, CYCLOPS -0.296 donor-order correlation,
A123/H123.4) while clocksync recovers the organ map (offsets r 0.945; against a donor-disjoint truth 0.679, a failed
0.7 bar, post hoc ceiling about 0.83, A123c). With full-day sampling the two-step pipelines recover coupling about as
well (A145), so our advantage is partial-day coverage, not exclusivity. The CHIRAL age dampening is largely a
death-type artefact (-0.73 pooled to -0.13 matched on Hardy; archive reanalysis, not re-verified here; D4).

## 5. Human skin clock

- Wu G, Ruben MD, Schmidt RE, Francey LJ, Smith DF, Anafi RC, Hughey JJ, Tasseff R, Sherrill JD, Oblong JE, Mills KJ,
  Hogenesch JB. Population-level rhythms in human skin with implications for circadian medicine. *Proc Natl Acad Sci
  USA* 2018;115(48):12313-12318. doi:10.1073/pnas.1809442115. PMID 30377266. **Verified** (Wu is first author;
  Hughey and Anafi are co-authors). Timed and single-sample skin cohorts ordered with CYCLOPS; a panel estimates
  phase within about 3 h from one sample.
- del Olmo M, Spörl F, Korge S, Jürchott K, Felten M, Grudziecki A, de Zeeuw J, Nowozin C, Reuter H, Blatt T,
  Herzel H, Kunz D, Kramer A, Ananthasubramaniam B. Inter-layer and inter-subject variability of diurnal gene
  expression in human skin. *NAR Genom Bioinform* 2022;4(4):lqac097. doi:10.1093/nargab/lqac097. PMID 36601580.
  **Verified** (added; the closest match to our A37 epidermis-dermis design). Core-clock expression is similar across
  layers at the population level but varies between subjects.

How we differ. Both studies place samples in time or estimate a phase per sample; the psoriasis trait component
(A115) needs a reference for between-subject phase variability, which these provide. We use the within-person
epidermis-dermis offset instead (A37: untimed -1.69 h vs timed -2.14 h, 11 subjects, moderate).

## 6. Surgery and treatment time of day

- Montaigne D, et al. (author list not retrieved). Daytime variation of perioperative myocardial injury in cardiac surgery
  and its prevention by Rev-Erbα antagonism: a single-centre propensity-matched cohort study and a randomised study.
  *Lancet* 2018;391:59-69. PMID 29107324. **Corrected.** It concerns aortic valve replacement (596 matched patients
  plus an 88-patient randomised study), not cardiac surgery in general; fewer major cardiac events in the afternoon
  (HR about 0.5) with a Rev-erbα mechanism. DOI and full author list not retrieved (the PMID is from a PubMed
  result).
- Fudulu DP, et al. Lancet Reg Health Eur 2021 (meta-analysis and UK database analysis; first author and full title
  seen only in a secondary summary). **Candidate, contradicts Montaigne:** no significant diurnal difference in
  troponin release or outcomes in aortic valve replacement. Citation details to verify.
- Landré T, Karaboué A, Buchwald ZS, Innominato PF, Qian DC, Assié JB, Chouaïd C, Lévi F, Duchemann B. Effect of
  immunotherapy-infusion time of day on survival of patients with advanced cancers: a study-level meta-analysis.
  *ESMO Open* 2024;9(2):102220. doi:10.1016/j.esmoop.2023.102220. PMID 38232612. **Verified** (the handoff's
  unspecified "infusion studies"): 13 retrospective studies, 1,663 patients, early infusion OS HR 0.50 (0.42-0.58).

How we differ. Time of day of an intervention is recorded; the physiological internal time of the patient is not,
and chronotype, shift work and illness move it by hours. Our untimed readout gives internal time on resected
material (V25.5: 42 breast patients, order Spearman 0.558, error 1.42 h without truth versus 2.10 h for "always 13:00";
V3 heart in the archive, to be ported, A174 confirms on current code). The literature is correlational and partly
contradictory, which strengthens the motivation but means the paper must not claim that internal time explains the
time-of-day effect.

## 7. Fetal clock ontogeny

- Li C, Yu S, Zhong X, Wu J, Li X. Circadian rhythms of fetal liver transcription persist in the absence of canonical
  circadian clock gene expression rhythms in vivo. *PLoS One* 2012;7(2):e30781. doi:10.1371/journal.pone.0030781.
  PMID 22383974. **Verified; correction to the description:** mouse, no rhythms of canonical clock genes in fetal
  liver, but JTK_CYCLE finds some low-amplitude rhythmic transcripts, so "no rhythms" is too strong.
- Stanton DL, Zhu L, Hansen PJ. The ontogeny of circadian clock gene expression during mouse fetal development.
  *Biochem Biophys Rep* 2025;44:102216. doi:10.1016/j.bbrep.2025.102216. PMID 40917721. **Corrected:** primary
  mouse study (not a review). Core clock genes in fetal brain, heart, liver and kidney rise late in gestation or after
  birth; more developed in brain and liver than in heart and kidney.
- Canaple L, Gréchez-Cassiau A, Delaunay F, Dkhissi-Benyahya O, Samarut J. Maternal eating behavior is a major
  synchronizer of fetal and postnatal peripheral clocks in mice. *Cell Mol Life Sci* 2018 (volume and pages not
  retrieved; journal inferred from the DOI prefix). doi:10.1007/s00018-018-2845-5. **Corrected/verified:** the
  "maternal feeding 2018" reference. Rev-Luc mouse imaging shows that in-utero peripheral clocks are mainly set by
  maternal feeding.
- Ohta H, et al. Maternal feeding controls fetal biological clock. *PLoS One* 2008 (PMC2432029; only the first author
  and title verified). Rat: restricted maternal feeding shifts fetal SCN and liver by 5 and 7 h independent of the
  maternal SCN. **Candidate** earlier basis.
- Torres-Farfan C, et al. Maternal melatonin effects on clock gene expression in a nonhuman primate fetus.
  *Endocrinology* 2006;147(10):4618-26. doi:10.1210/en.2006-0628. PMID 16840546. **Abstract verified (2026-11-02).**
  Capuchin fetus at ~90 % gestation: Bmal1, Per2 and MT1 oscillate in the fetal SCN and adrenal; maternal melatonin
  suppression changes the fetal SCN clock but not the adrenal one ("probably not" its Zeitgeber). Related: functional
  MT1 in the capuchin fetal adrenal (J Physiol 2004, PMID 14673186, listing only); a melatonin-responsive rat fetal
  adrenal clock (Endocrinology 2011, doi 10.1210/en.2010-1260, title unverified).
- Umemura Y, Koike N, et al., Yagita K. Cell and tissue-autonomous development of the circadian clock in mouse embryos.
  *PNAS* 2017. doi:10.1073/pnas.1703170114. PMID 28827343. Abstract level: no cell-autonomous oscillation in E10
  hearts, clear at E18; CLOCK translation suppressed early (Dicer/Dgcr8).

How we differ. Fetal peripheral clocks are weak, partly driven by the mother, and mature late, so the clock-gene
shared state across organs that clocksync needs may not exist or may be maternal-cue-driven (D1 prediction). A181
could not test this: adult GTEx at the Cao fetal design passes the gate in only 0.20 of subsamples, so the fetal
gate failure (whole organs p 0.746) is uninformative and design-limited, not depth-limited. Within one compartment
(A189) the fetal adrenal carries a clock-specific coherent clock at 13-17 weeks, the first human fetal clock evidence we
know of, matching the primate fetal adrenal oscillation above; its cross-class coupling (A190) is confounded by cortex
RNA. The D3 scout found no paediatric multi-compartment series.

## 8. Fetal atlases

- Cao J, O'Day DR, Pliner HA, et al. A human cell atlas of fetal gene expression. *Science* 2020;370:eaba7721.
  doi:10.1126/science.aba7721. GEO GSE156793. **Verified** (full author list from memory beyond "J. Cao et al.").
  Over 110 samples, 15 organs, about 4 million cells, 77 main cell types (second trimester).
- "Human Developmental Cell Atlas preprint 2026 (4-22 post-conception weeks)": **not found.** Closest verified
  candidates: the Human Developmental Cell Atlas framework paper (*Development* 2018, seen as a search hit, citation not opened); Cell 2022 fetal lung atlas (5-22 post-conception weeks, lung only, citation not pinned
  down); CS12-16 embryo atlas, bioRxiv 10.1101/2021.11.30.470583 (4-6 weeks). Ask the PI for the preprint.

How we differ. Cao covers a narrow window of mid-gestation with few donors sharing organs; it cannot show a shared
state even where adults have one (A181). Later gestation and postnatal series need a different design (D2).

## 9. Other untimed multi-tissue and single-cell phase work (2023-2026, searched once)

| reference | relevance | status |
|---|---|---|
| Duan J, Ngo MN, Karri SS, Tsoi LC, Gudjonsson JE, Shahbaba B, Lowengrub J, Andersen B. tauFisher predicts circadian time from a single sample of bulk and single-cell pseudobulk transcriptomic data. *Nat Commun* 2024;15:3840. doi:10.1038/s41467-024-48041-6 | timed supervised prediction per sample, cell-type phase differences in mouse skin; not host-referenced | verified |
| Salati A, Paychère Y, Hahaut V, Gobet C, Naef F. Inferring circadian phases and quantifying biological desynchrony across single-cell transcriptomes (scRitmo). bioRxiv 2026. doi:10.64898/2026.03.30.715278 (DOI prefix unusual for bioRxiv: verify against the lab record) | lab's own; single-cell phase and desynchrony within a tissue | verified |
| Tempo: an unsupervised Bayesian algorithm for circadian phase inference in single-cell transcriptomics. *Nat Commun* 2022 (PMC9630322) | single-cell untimed phase; authors and DOI not checked | candidate |
| Transcriptional noise sets fundamental limits to decoding circadian clock phase from single-cell RNA snapshots. bioRxiv 2024 (10.1101/2024.06.30.601408), iScience 2026 | limits of single-cell phase from core clock genes (mouse fibroblasts) | candidate (authors not checked) |
| DCPR: a deep learning framework for circadian phase reconstruction. *BMC Bioinformatics* 2025 (doi:10.1186/s12859-025-06363-2) | untimed bulk phase reconstruction, single tissue | candidate (authors not checked) |
| CYCLOPS 2.0: improving and applying circadian phase reconstruction in confounded datasets (Drexel thesis) and a 2025 bioRxiv on human brain cell-type rhythms (10.1101/2025.02.21.639499) | covariate-aware ordering; brain; one tissue | candidate (thesis, preprint) |
| BayCT, Bayesian rhythmic model for post-mortem brain phase (bioRxiv 10.64898/2026.01.21.700905) | post-mortem brain phase only | candidate |
| TimeTeller (multigene dynamical-systems tool, Vlachou et al.) | single-sample phase and clock functionality; no inter-organ readout | candidate (record not opened; authors from memory) |

No paper found in these two searches estimates untimed, host-referenced coupling, relative amplitude and phase offset
between compartments of the same person, or reports inter-organ clock coherence in post-mortem or biopsy data as a
disease readout. This is an absence in a limited search, not a proven gap; a PubMed search of "inter-organ" or
"inter-tissue" circadian coherence and a bioRxiv alert for 2026 should be repeated before submission.

## 10. Tumour clock state and survival (searched 2026-11-02, after A191 was scored)

- Cadenas C, van de Sandt L, Edlund K, et al. Loss of circadian clock gene expression is associated with tumor
  progression in breast cancer. *Cell Cycle* 2014;13(20):3282-3291. doi:10.4161/15384101.2014.954454. **Verified**
  (Europe PMC full text). 766 node-negative, systemically untreated patients (Mainz, Rotterdam, Transbig), metastasis-free
  survival. Univariate: higher CLOCK, PER1-3, CRY2, NPAS2, RORC protective (e.g. CRY2 HR 0.781, PER3 0.806), TIMELESS
  adverse (1.406). Multivariate (age, size, grade, ER, HER2; n 480): only PER3 (0.658), RORC (0.421) and TIMELESS (1.607)
  remain. Not adjusted for proliferation, although the genes correlate strongly with the proliferation metagene (CRY2
  R -0.612, PER2 -0.440, PER3 -0.409, TIMELESS +0.674; TIMELESS is a replication gene). Single-gene, per-gene tests.
- Li et al. (with Hammarlund and Clarke; initials and author order not checked). Tumor circadian clock strength influences metastatic potential and
  predicts patient prognosis in luminal A breast cancer. *PNAS* 2024;121(7):e2311854121. doi:10.1073/pnas.2311854121.
  PMC10873596. **Verified** (full text; middle authors not checked). CYCLOPS 2.0 magnitude in untimed TCGA luminal A
  tumours; the high tertile has higher 5-year mortality (ANOVA p 0.047; logistic with age and metastasis at diagnosis
  p < 0.05, about 1.5-fold risk; p < 0.01 with MKI67 added). Survival subset apparently n 70; no HR, no interval,
  no grade, stage or PAM50 adjustment in a model; the death-or-new-event endpoint is not significant. Direction is the
  opposite of "clock loss is bad": stronger rhythm, worse outcome (via EMT genes).
- Vlachou D, et al. TimeTeller (bioRxiv 10.1101/622050, 2019; methods paper *PLoS Comput Biol* 2024;20(2):e1011779).
  **Verified (preprint).** Single-sample clock dysfunction metric in pre-treatment biopsies of the REMAGUS neoadjuvant
  trial (224 tumours): good clock function is favourable, multivariate HR 0.24 (p 0.021) in grade 1-2 tumours (n 96),
  adjusted for ER, PR, HER2, triple-negative, grade, size, nodes, LVI. No proliferation (Ki67) term and no check of
  the metric against proliferation; subgroup result. The survival part is not in the 2024 journal version as far as
  the search shows.
- Ye Y, et al. The genomic landscape and pharmacogenomic interactions of clock genes in cancer chronotherapy. *Cell
  Syst* 2018. doi:10.1016/j.cels.2018.01.013. **Candidate** (abstract only): clock-gene expression alterations relate
  to survival, stage and subtype across TCGA types; per-gene, no proliferation adjustment checked.
- Chang C, Lai L-C (as indexed; authors not checked), *J Transl Med* 2019 (10.1186/s12967-019-1880-9). **Candidate.**
  "ClockLoss" gene down-regulation associated with mortality in five of 21 cohorts, "independent of clinicopathological
  features"; proliferation adjustment not seen.
- Clock-gene risk-score papers (LASSO signatures in lung, liver, kidney and other cancers) are numerous; those seen
  report enrichment of proliferation pathways in the high-risk group and do not adjust for it.
- Background for the confound: Venet D, Dumont JE, Detours V. Most random gene expression signatures are
  significantly associated with breast cancer outcome. *PLoS Comput Biol* 2011. **Candidate** (from memory; record not
  opened): the association is carried by proliferation, which is why A191 adjusts for G2M and uses pseudo-clock sets.

How we differ. Earlier survival associations use single clock genes or untimed orderings, without a proliferation
term and without matched random gene sets, mostly in single cohorts or subgroups. Several of the genes (TIMELESS,
CRY2, PER2/3) track the proliferation metagene closely, so "clock loss predicts poor outcome" is largely expected from
the proliferation-outcome association (Venet 2011) together with A176 (proliferation swamps the clock). A191 tests the
clock-plane share with age, size, grade, nodes, G2M and PAM50 in the model and against 200 pseudo-clock sets: HR 0.977
[0.881, 1.083] per SD in 1,020 ER+ endocrine-treated METABRIC tumours, unadjusted 0.923 (p 0.105), consistent with an
association that disappears once proliferation and stage are in (A127 SCAN-B off-plane HR 0.94 likewise). Li et al.
2024 is the one result not explained this way (MKI67 in the model, opposite direction); it is luminal A only, n about
70, and uses CYCLOPS magnitude, not a clock-specific share. A luminal A-only METABRIC check would be post hoc.

Tested since (registered A192-A194, `results_preliminary.md`):
- **Cadenas (A192):** the same three public cohorts (684 tumours, 215 metastases). The direction reproduces. A score of
  their 7 protective genes (HR 0.81 per SD) sits at the 7th percentile of 1,000 random gene sets matched for correlation
  with proliferation, so it is not clock-specific. It falls to 0.88 (p 0.12) with G2M, and to 0.96 with AURKA.
  TIMELESS is fully absorbed. RORC keeps an independent association (0.80, p 0.019). It is not immune-driven (post hoc), and the other six genes are null after adjustment (0.93).
- **TimeTeller (A194):** Theta is not prognostic in SCAN-B and TCGA (400 deaths; HR 0.97 per SD). The grade 1-2 claim
  goes the other way (0.83, p 0.03, nominal).
- **Li et al. (A193):** CYCLOPS magnitude in TCGA, METABRIC and SCAN-B luminal A; running.

## Status table

| # | reference | status | note |
|---|---|---|---|
| 1 | Shilts, Chen, Hughey, PeerJ 2018 | verified | PMID 29404219 |
| 2 | Bieler et al., Mol Syst Biol 2014 | verified | PMID 25028488 |
| 3 | Feillet et al., PNAS 2014 | verified | author list from memory |
| 4 | Droin, Paquet, Naef, Nat Phys 2019 | verified | |
| 5 | "Lin and Farkas review ~2018" | corrected | Lellupitiyage Don et al., Cell Cycle 2019, primary study |
| 5b | Berlin breast clock subtypes 2024-25 | candidate | authors/journal unverified |
| 6 | "Ruben 2018 (CYCLOPS)" | corrected | CYCLOPS is Anafi et al., PNAS 2017; Ruben 2018 applies it to GTEx |
| 7 | Anafi et al., PNAS 2017 | verified | |
| 8 | Ruben et al., Sci Transl Med 2018 | verified | |
| 9 | Talamanca, Gobet, Naef, Science 2023 | verified | |
| 10 | Wu et al., PNAS 2018 skin | verified | |
| 11 | del Olmo et al., NAR Genom Bioinform 2022 | verified | added |
| 12 | Montaigne et al., Lancet 2018 | corrected | aortic valve replacement; PMID only |
| 13 | Fudulu et al., 2021 | candidate | contradicts Montaigne; unverified details |
| 14 | Landré et al., ESMO Open 2024 | verified | |
| 15 | "PLOS One 2012 fetal liver" | verified | Li et al.; mouse; low-amplitude rhythms present |
| 16 | "Biochem Biophys Rep 2025" | corrected | Stanton et al.; primary mouse study |
| 17 | "maternal feeding 2018" | corrected | Canaple et al., Cell Mol Life Sci; volume and pages unverified |
| 18 | Ohta et al., PLoS One 2008 | candidate | only first author and title verified |
| 19 | Cao et al., Science 2020 | verified | |
| 20 | "HDCA preprint 2026" | not found | candidates listed in section 8 |
| 21 | tauFisher, Nat Commun 2024 | verified | |
| 22 | scRitmo, bioRxiv 2026 | verified | |
| 23 | Tempo, transcriptional-noise paper, DCPR, CYCLOPS 2.0, BayCT, TimeTeller | candidate | see section 9 |
| 24 | Cadenas et al., Cell Cycle 2014 | verified | section 10; not proliferation-adjusted |
| 25 | Li, Hammarlund, Clarke et al., PNAS 2024 | verified | luminal A, CYCLOPS magnitude, MKI67-adjusted |
| 26 | TimeTeller survival (Vlachou et al., bioRxiv 2019) | verified (preprint) | REMAGUS, grade 1-2 subgroup |
| 27 | Ye et al., Cell Syst 2018; Chang and Lai, J Transl Med 2019; Venet et al., PLoS Comput Biol 2011 | candidate | section 10 |
