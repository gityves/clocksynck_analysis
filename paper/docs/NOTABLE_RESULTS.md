# Notable results (unexpected, strong, or likely to need a follow-up)

One entry per result: date, amendment, the result in one or two lines with its numbers, why it is notable, the
follow-up I would suggest, and the source file. Newest first. Numbers come from output files. Full write-ups are in
`results_preliminary.md`.

## 2026-10-31

- **A173: the proliferation clock-coherence loss replicates in gastric cancer (pre-specified; H173.1 −0.187 Holm 0.032,
  H173.3 −0.533 Holm 0.0015, both clock-specific; H173.2 −0.134 fails).** Lung adenocarcinoma minus squamous −0.230
  (p 0.058, MDE 0.355): adenocarcinomas lose coherence in all three lung series (mean −0.303), squamous tumours do not
  (−0.073).
  - *Why notable:* third cancer type with a replicated, clock-specific loss (breast, lung adenocarcinoma, gastric). All
    three are adenocarcinomas; squamous lung and A172's non-epithelial tumours do not show it.
  - *Caveats:* GSE62254 is at the 5 % pseudo line; the FFPE DASL cohort GSE26253 is +0.18; H173.4 is not significant.
  - *Follow-up I would suggest:* a registered "adenocarcinoma vs squamous" test where both are large (TCGA ESCA/HNSC/CESC
    are seen; oesophageal GEO cohorts are too small), or single-cell malignant-cell pseudobulks in gastric cancer.
  - *Source:* `A173_gastric_lung_histology/output/tests.csv`, `descriptive.csv`.

- **A172: the proliferation clock-coherence loss is not pan-cancer (pre-specified; H172.1 −0.032, p 0.18, MDE 0.087,
  12 unseen TCGA projects).** No project passes Holm. Descriptively the six carcinomas lean negative (−0.114, p 0.011,
  3rd pseudo percentile; STAD −0.305), while LUSC (+0.21) and LGG (+0.20) go the other way.
  - *Why notable:* bounds the breast / LUAD claim: the loss is cancer-type specific. LUSC vs LUAD makes histology, not
    organ, the axis.
  - *Follow-up I would suggest:* a registered gastric replication (GEO STAD cohorts, e.g. ACRG GSE62254) and an
    adeno vs squamous contrast within lung and oesophagus.
  - *Source:* `A172_tcga_unseen/output/per_project.csv`, `pooled.csv`.

- **A171: the proliferation clock-coherence loss replicates in lung adenocarcinoma (pre-specified; H171.1-2 pass,
  Holm 0.002 and 0.0135).** GSE72094 −0.242 (442 tumours), GSE68465 −0.293 (443, clock-specific, 1st pseudo percentile);
  GSE31210 −0.232 descriptively. Hepatocellular carcinoma does not replicate (−0.12, −0.10; powered for the TCGA −0.51).
  - *Why notable:* the first cancer beyond breast with an independently replicated, pre-specified clock positive. It
    makes "proliferating tumours lose clock coherence" a multi-cancer claim.
  - *Caveats:* GSE72094 is at the 5.5th pseudo percentile (not specific by the 5 % rule). No LUAD driver mutation or
    smoking group carries the loss. HCC: tumours are less coherent than paired non-tumour liver (−0.327, descriptive).
  - *Source:* `A171_luad_lihc/output/tests.csv`, `descriptive.csv`.

- **A170: ER− malignant cells lose clock-gene coherence, but not clock-specifically (pre-specified, p 0.044, 12th pseudo
  percentile).** The whole-tumour pseudobulks of the same 38 patients show no ER effect.
  - *Source:* `A170_sc_breast_er/output/tests.csv`.

- **A168: two weeks of aromatase inhibition leave breast-tumour clock coherence unchanged (pre-specified; H168.1-3
  fail).** POETIC treated change −0.010 (MDE 0.196) while proliferation (−0.90 SD) and oestrogen response (−0.96 SD)
  fall; GSE20181 and Z1031B point the same way.
  - *Why notable:* the ER− coherence loss (A127, A166) is not reproduced by acute ER blockade. It reads as a stable
    property of the ER− tumour state, not ongoing ER signalling.
  - *Caveats:* the causal difference-in-differences is underpowered (MDE 0.43); the 56 untreated pairs fall in coherence
    (0.680 → 0.289), partly through PER1 / NR1D1, possibly a batch.
  - *Follow-up:* longer deprivation tested in A169 (pre-specified, H169.1-2 fail): no fall after 3 months of letrozole
    (GSE59515 −0.049, MDE 0.307) and no progressive trend in three series. Open: single-cell ER+ vs ER− malignant cells
    (is the loss cell-intrinsic?).
  - *Source:* `A168_endocrine/output/tests.csv`, `A169_endocrine_duration/output/tests.csv`.

## 2026-10-30

- **A166: the ER-negative breast-tumour clock-coherence loss replicates in METABRIC, clock-specifically (pre-specified).**
  H166.1 T −0.229, Holm 0.003, 1st pseudo percentile. Within PAM50 × proliferation, T −0.168, Holm 0.047, 2nd
  percentile. The effect beyond the null mean is about half of SCAN-B's. The proliferation loss does not replicate
  (p 0.15, underpowered).
  - *Why notable:* this is the first pre-specified, independently replicated cancer clock positive across two platforms
    (RNA-seq and microarray). It links ER signalling to tumour clock integrity, beyond subtype and proliferation.
  - *Caveats:* the within-subtype result is at the Holm line, and with ER_STATUS labels it is not detected.
  - *Follow-up (PI):* survival or endocrine response (not registered; confounded by indication); a mechanism for ER and
    clock coherence (ER-target signature in ER+ Basal tumours).
  - *Source:* `A166_metabric/output/a166_tests.csv`.

- **A164: T2D and obesity are equivalent to no coupling change in the metabolic organs (matched GTEx).** T2D pooled
  module d* +0.019 [−0.021, 0.060] and obesity +0.004 [−0.081, 0.089], both within ±0.15 and powered. No effect.
  - *Why notable:* the first registered, powered equivalence. "Metabolic disease leaves organ clocks coupled to the
    body clock" becomes a tested claim.
  - *Follow-up:* the A165 phase, amplitude and noise readouts could be applied to the same matched design (register
    first). Vent obese pancreas d* +0.36 is a descriptive lead only.
  - *Source:* `A164_gtex_metabolic_matched/output/a164_tests_hypotheses.csv`.

- **A153: the corrected disease contrasts have no primaries.** No A17 disease cohort reaches MDE ≤ 0.2 on d* in an
  eligible family. The smallest is 0.232 (dcm_healthy heart fibroblast), and 0.214 for endothelial in the descriptive
  scope. The A17 islets beta effect does not survive (p* 0.189).
  - *Why notable:* this is the main obstacle to the "pre-specified disease positive" Cell needs. Per-arm cohort sizes
    of current single-cell disease atlases are too small for a 0.2 coupling change.
  - *Follow-up:* within-person designs (6.6) and GTEx matched contrasts (A164); the PI may consider option (b) or larger cohorts.
  - *Source:* `A153_corrected_contrasts/output/readout_all5.csv`, `mde_table.csv`.

- **A127 post hoc: the ER− clock coherence loss holds within PAM50 subtype and proliferation tertile, clock-specifically
  (post hoc).** Basal −0.379, Her2 −0.305; within subtype × proliferation tertile T −0.268, p 0.0045, below all 200
  pseudo sets.
  - *Why notable:* it separates an ER (hormone-receptor) effect on clock coherence from intrinsic subtype and
    proliferation. It may link to the oestrogen regulation of clock genes.
  - *Follow-up:* replicate in TCGA BRCA or METABRIC Basal/Her2 (register first), and test whether ER+ Basal tumours
    have an ER-target signature.
  - *Source:* `A127_scanb/output/posthoc_er_within_pam50.csv`.

- **A162: the coupling estimate is attenuated by host error (method).** Coverage fails at c ≥ 0.6 even with anchors
  at 0.8 (V162.1 19/32 cells). The shrinkage (×0.87-0.92) is exactly the errors-in-variables term: the host's
  posterior error enters psi.
  - *Follow-up:* A163, a correction that removes |z|² tr S / 2 from psi, validated on fresh seeds.
  - *Source:* `A154_gate_calibration/output/a162/`.

## 2026-10-29

- **A154 V154.3: coupling estimates are biased low and their bootstrap intervals undercover when anchors are weakly
  coupled (method weakness).** Coverage is median 0.83 / 0.74 / 0.35 at c 0.3 / 0.6 / 0.9. In a noise-free check, ĉ for
  0.3 is 0.084 with weak anchors and 0.288 with precise ones.
  - The gate is close to calibrated: N0 0.040-0.085, the N1 factor map about 0.07. Contrast p is uniform.
  - *Follow-up:* A162 simulation at realistic anchor coupling (Q17).
  - *Source:* `A154_gate_calibration/output/v154_3_coverage.csv`, `diag_nominal_c.csv`.

- **A126: endometrial epithelial clock coherence collapses from the proliferative (0.653) to the secretory phase
  (0.126), p 0.0005 (unexpected, descriptive).** Meanwhile the registered tests fail:
  - H126.0, cycling vs non-cycling −0.141, p 0.10;
  - H126.1, breast luteal tertile +0.029, p 0.58.
  - *Why notable:* progesterone may reorganise the endometrial clock. This is a hormone, not proliferation, effect,
    since cycling cells show no loss.
  - *Follow-up:* a registered test on an independent endometrium series (Ulrich 2024, Census, 20 donors; phase labels
    to verify). Check whether the drop is amplitude or coherence (host-referenced coupling with stromal host).
  - *Source:* `A126_menstrual/output/h126_0_endometrium_phase_descriptive.csv`.
  - **Update (A161, same day):** not confirmed in bulk. In GSE51981 secretory coherence is higher (T +0.096, MDE
    0.117), and also in GSE6364 and GSE141549. The drop is epithelium- or GSE266265-specific at most.

- **A127: SCAN-B replicates the A97 proliferation-coherence loss at 3,073 tumours (strong, a pre-specified positive).**
  - H127.1: T −0.123, Holm 0.0045, within PAM50 × purity, clock-specific at the 2.5th pseudo percentile.
  - H127.2: ER− −0.343, Holm 0.001, 1.5th percentile.
  - Low-purity tumours are more coherent, so purity is not the route.
  - No survival association (HR 0.94, p 0.35).
  - LumB coherence 0.267 vs LumA 0.754.
  - *Follow-up:* a within-PAM50 ER contrast (ER is nearly all Basal/Her2); is the LumB vs LumA difference only
    proliferation?
  - *Source:* `A127_scanb/output/a127_tests.csv`, `a127_specificity.csv`.

- **Bounded ccorr removes much of the cancer pairing evidence (reframing for Figure 5).**
  - Kept: H149.1 0.490 (p 0.011) and H149.3 0.707 (p 0.0004).
  - Not estimable on the bounded statistic: H149.2 (3 cohorts), H151.2, the endothelial pseudo companion, and H138.5
    (2 cohorts).
  - The A149 pairing is still not clock-specific (pseudo q95 0.602).
  - *Follow-up:* the fresh A151 registration on new cohorts (7.2) must use the bounded statistic and enough patients
    per cohort (flags remove units of n 4-11).
  - *Source:* `A149_paired_tumour_normal/bounded_2026-10-29/output/`, `A151_fibroblast_pairing/bounded_2026-10-29/output/`.
- **A134: the permutation p (0.017) and the patient bootstrap disagree (share ≤ 0 is 0.30).** The null tests association
  against unrelated covariates, not precision. "Copy number explains part of the gap" needs its CI [−0.16, 0.27].
  - *Source:* `A134_cancer_coherence/bounded_2026-10-29/output/h134_1_bootstrap.csv`.

- **A159: breast plasma cells are clock-specific (c 0.81, Δc 0.50), the only lymphocyte in 13 that is (unexpected).**
  Clock-specific means above both the matched and the responsive pseudo-gene q95. In the 3 covered atlases this holds
  for 3 of 31 compartments, all in breast: endothelial (Δc 0.31), epithelial (0.40) and plasma (0.50). Lymphocyte
  family median Δc is 0.015.
  - *Why notable:* this cuts against "lymphocytes carry no clock-specific coupling". It may be a resident plasma
    population, or a breast-specific artefact.
  - *Follow-up:* check plasma in other atlases once genome tables exist (Q15), and its donor count and depth.
  - *Source:* `A159_specific_coupling/output/d159_compartments.csv`.
- **A159: the per-compartment clock-specific flag is underpowered by construction (design problem).** A single
  compartment's pseudo-gene q95 is 0.37-0.92 (mostly 0.6-0.9), so only compartments with c above about 0.75 can pass. Heart
  fibroblast (c 0.718 vs matched median 0.579, q95 0.912) fails. Colon endothelial clears the matched q95 (0.472 vs 0.467; median 0.023) but not the responsive q95 (0.736).
  Fibroblast_stromal family median Δc is 0.233, but share clock-specific is 0. In addition, only 3 atlases are covered,
  and no disease control arm has a genome-wide table. So A153 has no primaries as registered (Q15).
  - *Source:* `A159_specific_coupling/output/{d159_compartments,d159_families,row_check,h159}.csv`.

## 2026-10-28

- **H109.1 and H138.4 survive the review's main alternative (strong).** Against responsive (GC / immediate-early) gene
  sets: ventilated C 0.463 vs q95 0.257, and tumour-normal pairing 0.365 vs q95 0.310 (both Holm 0.0099).
  PER1/PER2-free and TSS/plate-stratified nulls change nothing (stratified null mean at most 0.129).
  - *Caveat:* the responsive pool is small (Jaccard 0.79-0.91).
  - *Follow-up:* the Figure 3 / Figure 5 wording can drop the shared-stress caveat.
  - *Source:* `A138_cancer_estimability/output/responsive/holm_responsive.csv`.
- **Psoriasis core delay is robust and specific (strong, reassuring after A156).** Leave-one-gene-out and 80 % subsets
  keep it 100 % negative. 0 of 600 responsive or matched pseudo sets reach the core offset (-1.43 h; pseudo
  2.5-97.5 % about [-0.6, 0.6]). NR1D1 alone carries +0.54 h of it.
  - *Follow-up:* state the NR1D1 dependence; a composition-controlled readout (roadmap 6.9) is still open.
  - *Source:* `psoriasis_robustness_2026-10-28/output/pseudo_sets_summary.csv`.
- **A156: the untimed core-cue readout fails in timed human blood under forced desynchrony (unexpected weakness).**
  The timed cosinor shows the biology: cue genes follow the displaced sleep (-6.85 h [-11.1, -2.8]) and the core clock
  stays (-1.36 h). The untimed per-subject core-cue offset does not track it: r 0.04 against the timed change, MDE
  13 h. It moves 4.4 h with the mapping of 5 cue genes.
  - *Why notable:* the core-cue estimator carries the psoriasis claim. It is fragile at low signal (blood, 7 samples).
  - *Follow-up:* a host-referenced or timed-aware readout for perturbation designs. Check the psoriasis core-cue
    estimate's sensitivity to the cue-gene set (leave-one-gene-out).
  - *Source:* `A156_forced_desynchrony/output/descriptive_timed.csv`.
- **A124: no human perturbation positive from the simulated night shift (MDE 7.3 h, n = 11).** H124.1 |Δ core-cue|
  3.35 h, p 0.21. Descriptive: the muscle core clock keeps its 07->19 h direction under a 12 h behavioural shift but
  loses half its amplitude (3.20 -> 1.67), consistent with damping rather than shifting. H124.2 flips with the mapping of
  8 old gene names (p 0.77 vs 0.054).
  - *Follow-up:* a better-powered perturbation (timed series: GSE48113 forced desynchrony, blood, 22 subjects x 7
    times) or GTEx steroid users (roadmap 6.2-6.3).
  - *Source:* `A124_night_shift/output/tests.csv`.
- **A123c diagnostic: truth (iii) is too noisy to test offsets at 0.7 (explains the fail).** Split-half truth
  reliability 0.54, so the ceiling for a perfect estimate is about 0.83. The estimate disattenuated is 0.73-0.81. Whole
  Blood alone (truth (iii) R2 0.028) moves r from 0.679 to 0.754. The 0.807 split was a 95th-percentile draw.
  - *Follow-up:* a Figure 1 legend with both truths and the ceiling. A future working-hours test needs a truth with
    reliability > 0.9 (all-day donors, or a second cohort).
  - *Source:* `A123_working_hours/output/a123c/diag/split_half.csv`.
- **A123c: the Figure 1 offsets miss the bar against a donor-disjoint truth (unexpected).** Offsets r 0.679
  [0.467, 0.839] vs truth (iii) (H123.5 fail, threshold 0.7); amplitudes 0.965 (H123.6 pass). Truth (i) vs truth (iii) is
  only 0.821. A truth built from half the donors gives 0.807 on the same tissues, which is better than the full truth.
  - *Why notable:* the method's main real-data advantage rests on A123. Every untimed orderer, and the timed oracle on
    the same donors (0.491), stays far below.
  - *Follow-up:* the post hoc noise-ceiling diagnostic is running (`output/a123c/diag/`).
  - *Source:* `paper/exploratory/A123_working_hours/output/a123c/tests.csv`.
- **Review: pseudo-clock genes already give a coupling of 0.249 (vs 0.404 for clock genes) in 5 atlases (A150)
  (strong, reframing).** Lymphocyte coupling (0.094) equals its pseudo-clock level (0.087).
  - *Why notable:* raw coupling is not clock-pure, and Figure 4 nulls on immune classes are uninformative.
  - *Follow-up:* roadmap 0.1 / 3.3 (Δc readout).
  - *Source:* `paper/exploratory/A150_atlas_specificity/output/h150_1.csv`.
- **Review: with the anchor-free host, the "stromal-vascular backbone" largely flattens (unexpected).** Fibroblast 0.709,
  parenchyma 0.670, epithelial 0.600; only immune cells and glia separate.
  - *Follow-up:* show the learned, fixed and free rankings side by side after the class fix (roadmap 3.2).
  - *Source:* `paper/results/figure_2/output/backbone/pre_classfix/family_medians.csv`.
- **Review: A111 H111.2 is a design failure, not a wrong hypothesis.** Under daytime feeding, WAT and lung are out of
  phase and their joint gate fails (p 0.148), so no coherent host exists.
  - *Follow-up:* the per-arm host gate in `contrast()` (roadmap 1.4).
  - *Source:* `paper/exploratory/A111_feeding_controls/output/a111_gates.csv`.
