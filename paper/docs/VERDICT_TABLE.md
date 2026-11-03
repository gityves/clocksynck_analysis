# Verdict table: every registered test from A109 onward

Written 2026-11-01; A174-A181 and the A150 rerun added the same day (counts below are those before the additions: A174-A181 add 6 passes, 6 failures, 1 not scored, 1 not informative, with A180 counted as two rows; A182-A185 add 2 passes (not clock-specific), 3 failures, 4 not estimable; A136 adds 1 pass, 2 failures; A186 adds 1 pass (not clock-specific); A187 adds 2 not scored; A188 adds 1 descriptive, 1 not scored). Procedure: analysis plans committed before scoring (each hypothesis is an amendment in
`paper/docs/registered_experiments.md`, committed before its data were scored); this is internal pre-specification with a git
history, not external pre-registration. Anything chosen after seeing results is labelled post hoc and is not counted below.
Verdicts are those of `paper/docs/results_preliminary.md` (RP), which governs over `state_2026-10-31/inventory.md`; numbers are copied
from RP. Checked against the output tables on 2026-11-01: A109, A111, A112, A113, A121-A124, A127, A130, A138, A142-A146, A149-A152, A154, A156,
A157, A161, A163-A173. Taken from RP and the inventory without reopening the table: A110, A126, A132-A135, A139, A140, A147, A153, A155, A159, A162.

## Counts (amendments A109 onward: registered hypotheses and validation criteria)

- **Tests listed: 145**
- **Passed: 45**  (a pass of an exclusion rule counts as a pass: H109c.1, H110.3, H110.5; a pass is not necessarily clock-specific, see the verdict text)
- **Failed: 68**  (includes 5 'no effect detected at the stated MDE' and the gate and simulation-criterion failures)
- **Not scored: 22**  (never computed: infeasible, withdrawn, gate failed, parked, dropped, not registered)
- **Scored, no verdict: 10**  (descriptive by registration, uninformative, or not scorable)
- Archive results cited by a figure: 16 rows, listed separately and not in the counts.

Rows that bundle several hypotheses (A116, A140 H140.1-3 and P140.1-3, A148, A159) count once per row. The inventory tally
(`state_2026-10-31/inventory.md` section 16) uses hypothesis-ID level with other conventions (it counts H110.5 as a failure and
lists A109b and the simulations separately); the counts here follow the verdict column below.

## Conventions

- **Level.** L1 gene order within a compartment (coherence, core-cue offset); L2 compartment-to-host coupling c, relative amplitude |z|, phase offset arg z; L3 organ phase relations (C between organs, offset preservation, r*); method (simulation, benchmark, validation); other.
- **P / S.** P = labelled primary in the registration, or the headline test of an amendment with no label; S = labelled secondary, specificity or other companion, sensitivity, or outside the confirmatory family.
- **p.** Raw p (permutation unless stated; floors are 1/(n+1)); Holm as registered, over the registered family, unscored members taken as not smaller (see the 2026-11-01 errata in `results_preliminary.md` for the A109 and A138 families). m = family size used.
- **Source paths** are relative to `paper/exploratory/` unless they start with `paper/` or `RP` / `RE` (`paper/docs/results_preliminary.md`, `registered_experiments.md`).
- **Not in the table (not tests):** registered descriptives (D-items, A141, A158, A144 part 2), post hoc analyses (A132b, A113 muscle with Fisher-Lee, A127 within PAM50, A149 / A151 bounded statistics, A142 bounded rescore, psoriasis robustness) and integrity re-readings (WP0). H114.1 was never registered (A114 was replaced by A138); A136 is a draft.
- **A111 H111.1** is one hypothesis registered on two datasets; it is listed once per dataset.

## Tests, A109 onward

| amendment | hypothesis | level | P/S | statistic and value | p (raw / Holm) | verdict as registered | MDE or bound | source |
|---|---|---|---|---|---|---|---|---|
| A109 | H109.1 | L3 | P | C_clock 0.463 vs pseudo-clock median 0.154; excess 0.310 [0.225, 0.348] | 0.005 (floor 1/201) / Holm 0.025 (m=5), 0.010 (m=2) | PASS | MDE 0.018 | A109_coupling_nulls/output/h109_1_null.csv; A109_coupling_nulls/output/excess_bootstrap.csv |
| A109 | H109.1 sec (a, e): fast, Hardy 3, Hardy 4 | L3 | S | excess 0.469 / 0.246 / 0.255 over pseudo-clock null | 0.005 each (floor) | PASS (each) | - | A109_coupling_nulls/output/h109_1_null.csv |
| A109 | H109.2 | L3 | P | Fisher-Lee rho 0.261 (66 pairs); folded distance 0.60 vs 1.02 h, p_closer 0.001 | 0.165 / Holm 0.660 (m=5), 0.165 (m=2) | FAIL (ventilated donors) | power 80 % at rho 0.29 | A109_coupling_nulls/output/h109_2_offsets.csv; A109_coupling_nulls/output/h109_2_power.csv |
| A109 | H109.2 sec (e): Hardy 4 | L3 | S | rho 0.784; folded 1.06 vs 1.02 h, p_closer 0.57 | 0.0002 | PASS | - | A109_coupling_nulls/output/h109_2_offsets.csv |
| A109b | C vs pseudo-clock null (Tabula Sapiens) | L3 | S | C_clock 0.683 vs null 0.675; excess +0.009 | 0.438 | FAIL (not replicated, underpowered) | MDE 0.104 | A109_coupling_nulls/output/a109b_c_null.csv; A109_coupling_nulls/output/a109b_power.csv |
| A109b | offset preservation (Tabula Sapiens) | L3 | S | rho 0.091 (6 pairs, 24 permutations) | 0.417 | uninformative (registered feasibility: power 0.06-0.13) | - | A109_coupling_nulls/output/a109b_offsets.csv; A109_coupling_nulls/output/a109b_power.csv |
| A109c | H109c.1 (Hardy 4 wall time; exclusion) | L3 | S | R 0.097 vs null median 0.101; null q95 0.145 | 0.56 | PASS of the exclusion reading (q95 < 0.15; the 80 %-power MDE 0.208 is not) | MDE 0.208 (80 % power) | A109_coupling_nulls/output/a109c_walltime.csv; A109_coupling_nulls/output/a109c_walltime_power.csv |
| A109c | H109c.2 (Hardy 4 offsets without wall time) | L3 | S | residual rho 0.750; p_closer 0.39 | 0.0002 | PASS | - | A109_coupling_nulls/output/a109c_offsets_residual.csv |
| A109 / A138 note | H109.1-R (responsive-gene companion) | L3 | S | C_clock 0.463 vs responsive median 0.214, q95 0.257 | 0.005 (floor) / Holm 0.0099 (family of 2) | PASS | - | A109_coupling_nulls/output/responsive/h109_1_responsive.csv |
| A110 | H110.1 (decay with time ill) | L3 | P | interval recorded for 16 of 103 Hardy 4 donors | - | not scored (infeasible: interval completeness < 60 %) | - | A110_mechanism/ |
| A110 | H110.2a (strong clocks carry coupling, Hardy 4) | L2 | S | Spearman(amplitude, coupling) 0.498 | 0.007 / Holm 0.014 | PASS | - | A110_mechanism/output/a110_h2.csv |
| A110 | H110.2b (retained coupling vs stromal fraction) | L2 | S | rho -0.857 over 7 tissues | 0.994 / Holm 0.994 | FAIL (not supported, opposite sign) | - | A110_mechanism/output/a110_h2b_tissues.csv |
| A110 | H110.2c (hubs) | L3 | S | mean drop 0.010 vs 0.003 | 0.32 | descriptive (registered as descriptive) | - | A110_mechanism/output/a110_h2c_hubs.csv |
| A110 | H110.3 (temperature; exclusion) | L3 | S | drop in C 0.017 [-0.000, 0.032] | random sets 0.44 | PASS of the exclusion rule (drop < 0.05) | bound 0.05 | A110_mechanism/output/a110_h3.csv |
| A110 | H110.4 (hypothalamus, organ donors vs Hardy 3-4) | L2 | S | fallback fast vs Hardy 3-4: +0.208 (hypothalamus +0.087 vs other regions -0.121) | 0.002 | descriptive fallback (organ-donor contrast infeasible: 0 of 301 with brain) | MDE 0.18 | A110_mechanism/output/a110_h4.csv |
| A110 | H110.5 (feeding; exclusion) | L3 | S | -0.0005 (2 within-digestive pairs) | 0.49 | PASS of the exclusion rule, weak (2 pairs); the inventory tallies it FAIL | MDE 0.20 | A110_mechanism/output/a110_h5.csv |
| A132 | H132.1 (organ donors vs ventilated post-mortem) | L3 | P | folded-distance difference +0.115 h (equal n), -0.057 h (full n) | 0.72, 0.68 | FAIL (not supported) | MDE 0.54 h / 0.44 h | A132_brain_death/output/a132_contrast.csv |
| A132 | H132.2 (ventilation-duration tertiles) | L3 | S | top minus bottom tertile -0.192 h | 0.155 | FAIL (not supported) | MDE 0.45 h | A132_brain_death/output/a132_contrast.csv |
| A132 | D132.3 (steroid users) | L3 | S | 8 and 5 users (< 15 per arm) | - | not scored (dropped) | - | RP A132 |
| A133 | V133.1 (circular-part estimator gate) | method | P | rho with truth at f 0.5 0.646 (< 0.8); agreement at f 0 0.941 (< 0.95) | - | FAIL (estimator stops) | - | A133_axial_decomposition/output/v133_gate.csv |
| A133 | H133.1 | L3 | S | - | - | not scored (V133.1 failed) | - | RP A133 |
| A147 | V147.1 (order separates loss from pull) | method | P | n 79 AUC 0.977; f 0.8 order above loss median 0.805 | false positive 0.050-0.065 | PASS | - | A147_pull_order/output/v147.csv |
| A147 | V147.2 (shrink ratio unbiased at no pull) | method | P | ratio at f 0 = 0.809 (n 79), limit [0.85, 1.15] | - | FAIL (overall A147 FAIL) | - | A147_pull_order/output/v147.csv |
| A111 | H111.1 Xin (offset change tracks timed change) | L3 | P | Fisher-Lee 0.988, 4 tissues | exact 0.042 (minimum attainable) | PASS | - | A111_feeding_controls/output/a111_h111_1.csv |
| A111 | H111.1 Manella | L3 | P | Fisher-Lee 0.986, 3 tissues | exact 0.167 (minimum attainable) | FAIL by design (registration flaw) | - | A111_feeding_controls/output/a111_h111_1.csv |
| A111 | H111.2 (liver coupling falls under daytime feeding) | L2 | P | Manella liver d -0.556, d* -0.370; timed R2 change +0.164 | 0.010 (d), 0.060 (d*) | FAIL (direction condition not met; WAT and lung drop as much) | - | A111_feeding_controls/output/a111_h111_2.csv |
| A111 / A129 | H111.3 (hepatocyte Bmal1 KO decouples liver) | L2 | S | d* -0.287 [-0.411, -0.118] (c 0.967 -> 0.676); WAT -0.009, lung +0.019 | 0.001 (floor) | PASS, specific | - | A111_feeding_controls/output/a111_h111_3.csv (Liver row) |
| A112 | H112.1 (cue vs locked genes lose lock under arrhythmic feeding) | L1 | P | median log lock ratio cue minus locked +0.230 [-0.096, 0.819] | 0.91 (less) | FAIL (wrong sign) | - | A111_feeding_controls/output/a112_h112_1.csv |
| A113 | H113.1 adipose GSE104674 | method | P | JS 0.888; Fisher-Lee 0.822; median error 1.11 h | 0.001 | PASS | threshold 0.3 | A113_living_time/output/a113_summary.csv |
| A113 | H113.1 muscle GSE108539 | method | P | registered JS -0.854 (Fisher-Lee 0.576, adopted after JS was seen) | 1.00 (JS); 0.001 (FL, post hoc) | FAIL as registered (post hoc pass with Fisher-Lee; A157 sends muscle to ED) | threshold 0.3 | A113_living_time/output/a113_summary.csv |
| A130 | H130.1 (internal minus wall time vs age) | method | S | -0.339 h per decade (54 subjects) | 0.054 | FAIL (not supported) | - | A113_living_time/output/a130_summary.csv |
| A131 | SMCENTER sun-time test | other | S | location codes absent | - | not scored (dropped, infeasible) | - | inventory section 3 (RE A113 / A123 / A130 / A131 block) |
| A123 | H123.1 (offsets vs truth (i)) | method | P | circular r 0.945 [0.907, 0.972]; median error 0.35 h | < 0.001 | PASS (truth (i) includes the test donors) | threshold in RE | A123_working_hours/output/tests.csv |
| A123 | H123.2 (relative amplitudes) | method | P | Pearson 0.976 [0.956, 0.987] | < 0.001 | PASS | - | A123_working_hours/output/tests.csv |
| A123 | H123.3 (donor order in window) | method | P | Spearman 0.613 [0.476, 0.723] | - | PASS | - | A123_working_hours/output/donor_time.csv |
| A123b | H123.4 (clocksync vs CoFE genome, CoFE clock, CYCLOPS) | method | S | offset r 0.945 vs -0.283 / -0.417 / -0.296 | 0.0005 each (floor) / Holm 0.0015 | PASS (all three) | - | A123_working_hours/output/a123b/h123_4.csv |
| A123c | H123.5 (offsets vs disjoint truth (iii)) | method | P | circular r 0.679 [0.467, 0.839]; median error 0.58 h | 0.0001 | FAIL (r > 0.7 required; p criterion alone passes) | threshold 0.7 | A123_working_hours/output/a123c/tests.csv |
| A123c | H123.6 (amplitudes vs truth (iii)) | method | P | Pearson 0.965 [0.948, 0.980] | 0.0001 | PASS | - | A123_working_hours/output/a123c/tests.csv |
| A134 | H134.1 (copy number) | L1 | P | malignant minus normal gain +0.077 (0.228 -> 0.291) | 0.017 | PASS, small effect (patient-bootstrap interval [-0.160, 0.266] includes 0) | MDE 0.097 | A134_cancer_coherence/output/a134_h1_h2.csv |
| A134 | H134.2 (MYC, hypoxia, cycling) | L1 | S | -0.092 | 0.96 | FAIL (not supported) | - | A134_cancer_coherence/output/a134_h1_h2.csv |
| A138 | H138.1 (malignant host coupling > 0) | L2 | P | pooled c 0.076 [0.014, 0.139], 225 patients; threshold 0.083 (rule changed in note 5) | 0.097 (3/31) / Holm 0.484 (m=6), 0.097 (m=2) | FAIL under the changed rule; the registered rule (lower bound > 0) would PASS (design change, see errata) | MDE about 0.2 | A138_cancer_estimability/output/score/tests.csv; A138_cancer_estimability/output/score/holm.csv; A138_cancer_estimability/output/power/h138_1_sim.csv |
| A138 | H138.2 (CoFE false validity) | other | P | - | - | not scored (Q2 simulator gate not passed; parked, ROADMAP 4.8) | - | A138_cancer_estimability/output/q2_gate_check/gate_verdict.csv |
| A138 | H138.3 (uncoupled inference) | other | P | - | - | not scored (as H138.2) | - | A138_cancer_estimability/output/q2_gate_check/gate_verdict.csv |
| A138 | H138.4 (tumour clock follows own normal) | L2 | P | pooled /corr/ 0.365 vs null 0.122 +/- 0.021 (q95 0.157), 517 pairs; pseudo-clock q95 0.261 | 1e-4 (floor) / Holm 6e-4 (m=6), 2e-4 (m=2) | PASS, clock-specific | MDE 0.174 | A138_cancer_estimability/output/score/tests.csv; A138_cancer_estimability/output/score/holm.csv |
| A138 | H138.5 (single-cell malignant vs normal epithelium) | L2 | S | 0.583 vs null 0.443 +/- 0.092 (52 patients) | 0.069 | descriptive (note 5 ii; out of the Holm family) | MDE 0.677 | A138_cancer_estimability/output/score/tests.csv |
| A138 | H138.6 (genome vs clock-gene CoFE) | other | P | - | - | not scored (as H138.2) | - | A138_cancer_estimability/output/q2_gate_check/gate_verdict.csv |
| A138 | H138.7 (healthy anchoring) | other | P | - | - | not scored (as H138.2) | - | A138_cancer_estimability/output/q2_gate_check/gate_verdict.csv |
| A138 note 6 | H138.8 (match does not depend on purity) | L2 | S | D (low minus high purity) +0.025 | 0.35 | FAIL | MDE 0.16 | A138_cancer_estimability/output/h138_8/tests.csv |
| A109 / A138 note | H138.4-R (responsive-gene companion) | L2 | S | 0.365 vs responsive median 0.252, q95 0.310 | 0.005 / Holm 0.0099 (family of 2) | PASS | - | A138_cancer_estimability/output/responsive/h138_4_responsive.csv |
| A139 | V139.1 (per-gene refit error) | method | P | S3 median error 1.88 / 1.12 h (limits 1.5 / 1.0) | - | FAIL | - | A139_waveform_refit/output/v139.csv |
| A139 | V139.2 (gene flags) | method | P | false flags up to 0.15 (E-a), 0.54 (E-b) | - | FAIL | - | A139_waveform_refit/output/v139.csv |
| A140 | Q1 simulator gate (three runs) | method | P | CYCLOPS median circular r 0.067 / 0.060 / 0.058; 0 of 50 at >= 0.7 | - | FAIL (gate; Q1 stops) | - | A140_dip_cyclops/output/q1/simulator_gate.csv; A140_dip_cyclops/output/q1rev/; A140_dip_cyclops/output/q1rev/simulator_gate_core.csv |
| A140 | H140.1-H140.3 (3 hypotheses) | method | P | - | - | not scored (simulator gate failed) | - | RP A140 Q1 |
| A140 | Q2 (i) gate, GSE266919 immune compartments | L2 | P | clock share 0.025 vs null 0.028 | pairing 0.001; acrophase 0.459 (0.540 all 78) | gate FAILS (uninformative, not 'decoupled') | - | A140_dip_cyclops/output/q2/part_i_summary.csv |
| A140 | P140.1-P140.3 (3 predictions) | L2 | S | - | - | not scored (not estimable: gate failed) | - | RP A140 Q2 |
| A140 | Q2 (ii) published pipeline | method | S | - | - | not scored (waits on Q1 gate, never run) | - | RP A140 Q2 |
| A115 | H115.1 (trait component, uninvolved vs healthy) | L1 | P | D115.1 core-cue offset -1.806 h [-3.811, -0.689] | 0.0084 | descriptive (PI 2026-10-27; left the Holm family) | MDE 1.91 h | A115_A116_psoriasis/output/a115_d115_1.csv |
| A116 | H116.1-H116.2 (2 hypotheses; treatment vs placebo) | L1 | P | - | - | not scored (withdrawn: data already used in A56) | MDE 2.15 h | A115_A116_psoriasis/output/a116_null.csv |
| A135 | reference rule (no H id): one family first-ranked in >= 4 of 5 control sets | L2 | P | excitatory neurons first in 2 of 5 (MSSM, HBCC); astrocytes 1, inhibitory 1, vascular 1 | - | FAIL (rule not met; no brain reference) | needs 4 of 5 | A135_brain_reference/output/a135_decision.csv |
| A142 | H142.1 (brain split brain-specific) | L3 | P | brain r* 0.639 [0.530, 0.728] vs random 13-tissue blocks median 0.978, q05 0.936 | 0.001 (floor) / Holm 0.002 | PASS | reference q05 0.936 | A142_organ_blocks/output/tests.csv |
| A142 | H142.2 (digestive departure) | L3 | P | r* 0.883 [0.819, 0.971] vs median 0.970, q05 0.900; uncorrected r 0.782 vs reference median 0.821 | 0.028 / Holm 0.028 (r*); 0.207 (uncorrected r) | PASS on the registered r*; NOT supported on the bounded r (erratum, 2026-11-01) | reference q05 0.900 | A142_organ_blocks/output/tests.csv; A142_organ_blocks/bounded_2026-10-29/output/targets_bounded.csv |
| A143 | H143.1 (not explained by quality) | L3 | P | adjusted r* 0.627 [0.550, 0.712] | rule: r* < 1 and upper < 1 | PASS | - | A143_brain_artefact/output/tests.csv |
| A143 | H143.2 (isolation batch) | L3 | S | same- minus different-batch +0.039 (30 pairs) | 0.19 / Holm 0.38 | FAIL (no brain-specific batch effect; excess above about 0.07 excluded) | null q95 0.070 | A143_brain_artefact/output/tests.csv |
| A143 | H143.3 (hypothalamus follows body) | L3 | S | alignment ratio 0.582 vs about 0.66; difference -0.082 [-0.248, 0.152] | 0.74 | FAIL | - | A143_brain_artefact/output/tests.csv |
| A144 | H144.1 (strong-clock association not shared noise) | L2 | P | cross-half rho 0.690 (5-95 %: 0.645-0.732) | 1e-4 (floor) | PASS | - | A144_strong_clocks/output/h144_1.csv |
| A145 | H145.1 (working hours) | method | P | MAE 0.060 vs best two-step (CYCLOPS) 0.074; wins 21 of 30 | 0.021 / Holm 0.064 | FAIL after Holm | - | A145_two_step/output/h145.csv |
| A145 | H145.2 (shared factor) | method | P | MAE 0.239 vs CHIRAL 0.234 | 0.82 | FAIL | - | A145_two_step/output/h145.csv |
| A145 | H145.3 (no clock) | method | P | gate pass 0.133 (4 of 30); false-claim 0.127 vs CoFE 0.06 | 0.99 | FAIL | - | A145_two_step/output/h145.csv |
| A146 | H146.1 (programs vs internal time, breast) | other | P | S -0.056 (42 patients) | 0.554 | FAIL (not detected) | MDE 1.09 | A146_tme_time/output/part1_tests.csv |
| A146 | H146.2 (programs vs recorded time) | other | P | S 1.339 | 0.0006 / Holm 0.0012 (H146.1-2), 0.0018 (H146.1-3) | PASS | MDE 0.94 | A146_tme_time/output/part1_tests.csv; A146_tme_time/output/part2_tests.csv |
| A146 | H146.3 (pan-cancer, internal time) | other | P | S 0.528 (213 patients, 16 cohorts) | 0.092 / Holm 0.184 | FAIL (not detected) | MDE 1.00 | A146_tme_time/output/part2_tests.csv |
| A146 note 2 | H146.4 (glucocorticoid control) | other | S | t slope -1.61 | 0.054 / Holm 0.054 | FAIL (in predicted direction) | - | A146_tme_time/output/checks_tests.csv |
| A146 note 2 | H146.5 (adjusted for age, subtype) | other | S | S adjusted 1.234 | 0.0005 / Holm 0.001 | PASS | - | A146_tme_time/output/checks_tests.csv |
| A121 | H121.1 (embryo clock onset in late window) | other | P | R late 0.594 - R early 0.392 = 0.202 (8 vs 8 embryos) | 0.226 / Holm 0.45 | FAIL | - | A121_embryo_timing/output/h121.csv |
| A121 | H121.2 (embryo phase) | other | P | frame offset 21.54 h (SE 0.30); datasets disagree 3.07 h (> 3 h) | - | uninformative (registered rule) | embryo offset SE 3.29 h | A121_embryo_timing/output/frame_offsets.csv |
| A122 | H122.1 (cue minus locked set) | L2 | P | +0.007 [-0.090, 0.098] (rerun on relearned anchors -0.017) | 0.44 / Holm 0.55 (rerun 0.645 / 0.645) | FAIL | - | A122_cue_controls/output/h122.csv (holds the rerun) |
| A122 | H122.2 (pseudo-cue sets) | L2 | P | observed 0.087 vs pseudo median 0.072, q95 0.114 (rerun 0.097) | 0.27 / Holm 0.55 (rerun 0.254 / 0.507) | FAIL | - | A122_cue_controls/output/h122.csv (holds the rerun) |
| A148 | H148.1-H148.2 (2 hypotheses; psoriasis decomposition) | L1 | S | acceptance check: cycling -0.020 [-0.050, 0.011] | - | not scored (stopped at the registered acceptance check) | - | A148_psoriasis_decomposition/output/acceptance.csv |
| A149 | H149.1 (malignant vs normal epithelium, single cell) | L2 | S | 0.534 vs null 0.446 +/- 0.042 (q95 0.517), 190 patients; pseudo-clock median 0.579 | 0.022 / Holm 0.044 | PASS nominal; not clock-specific | MDE (statistic level) 0.552 | A149_paired_tumour_normal/output/score/tests.csv |
| A149 | H149.2 (via stroma) | L2 | S | 0.550 vs 0.539 (61 patients) | 0.45 | FAIL | MDE 0.768 | A149_paired_tumour_normal/output/score/tests.csv |
| A149 | H149.3 (stroma pairing) | L2 | S | 0.607 vs 0.410 (51 patients); no specificity control registered | 0.0077 / Holm 0.023 | PASS (no specificity control) | MDE 0.608 | A149_paired_tumour_normal/output/score/tests.csv |
| A150 | H150.1 (atlas coupling clock-specific) | L2 | S | mean c 0.393 vs pseudo mean 0.242, q95 0.313 (53 compartments; rerun 2026-11-01 on class-fixed tables; was 0.404 vs 0.249, p 0.005) | 0.00995 | PASS | - | A150_atlas_specificity/output/h150_1.csv |
| A151 | H151.1 (tumour fibroblasts pair with normal) | L2 | S | 0.860 vs null 0.504 +/- 0.058 (98 patients) | 1e-4 / Holm 2e-4 | PASS | MDE 0.650 | A151_fibroblast_pairing/output/tests.csv |
| A151 | H151.2 (clock-specific) | L2 | S | pseudo median 1.013, q95 1.693 | 0.77 | FAIL (uninterpretable: statistic exceeds 1) | - | A151_fibroblast_pairing/output/tests.csv |
| A152 | H152.1 (residency carries lymphocyte coupling) | L2 | S | resident minus circulating -0.092 | 0.856 | FAIL | MDE 0.220 | A152_lymphocyte_residency/output/h1521.csv |
| A124 | H124.1 (core-cue offset, simulated night shift) | L1 | P | /d/ 3.35 h | 0.211 / Holm 0.422 | FAIL | MDE 7.3 h | A124_night_shift/output/tests.csv; A124_night_shift/output/mde.csv |
| A124 | H124.2 (core vs core-locked) | L1 | S | -0.43 h | 0.771 | FAIL | MDE 7.3 h | A124_night_shift/output/tests.csv |
| A156 | H156.1 (core-cue offset, forced desynchrony) | L1 | P | 4.80 h | 0.365 / Holm 0.731 | FAIL | MDE 13.1 h | A156_forced_desynchrony/output/tests.csv; A156_forced_desynchrony/output/mde.csv |
| A156 | H156.2 (core vs core-locked) | L1 | S | 3.04 h | 0.468 / Holm 0.731 | FAIL | MDE 13.1 h | A156_forced_desynchrony/output/tests.csv |
| A154 | V154.1 (N0 gate pass rate <= 0.075 in every design) | method | P | n100_K8_single_cell 17/200 = 0.085; the other 7 designs 0.040-0.070 | binomial 0.024 | FAIL, narrowly (literal reading; see errata) | 0.075 | A154_gate_calibration/output/v154_1_n0.csv |
| A154 | V154.2 (global-null contrast share < 0.05 <= 0.075) | method | P | 0.038-0.058 in all 8 designs | KS 0.17-0.94 | PASS | 0.075 | A154_gate_calibration/output/v154_2_null.csv |
| A154 | V154.3 (coverage of 95 % interval of c >= 0.90) | method | P | fails in 22 of 24 cells; median coverage 0.83 / 0.74 / 0.35 at c 0.3 / 0.6 / 0.9 | - | FAIL | 0.90 | A154_gate_calibration/output/v154_3_coverage.csv |
| A155 | decision rule (LOO anchor readout) | method | P | anchor amplitude bias 13.0-13.6 % (bulk), 18.2-19.6 % (single cell) at 3 anchors | - | met | - | A155_anchor_self_inclusion/output/decision.csv |
| A157 | V157.1 (JS usable on a 4-h grid) | method | P | muscle JS < 0 in 37.2 %, FL > 0.5 in 99.6 % (pass); adipose JS < 0 in 5.0 % (25/500) vs '< 5 %' | - | FAIL (literal: 5.0 % is not < 5 %; see errata) | 5 % | A157_js_grid/output/a157_summary.csv |
| A162 | V162.1 (coverage at c_anchor 0.8) | method | P | 19 of 32 cells below criterion; median 0.927, min 0.000 | - | FAIL | 0.90 | A154_gate_calibration/output/a162/summary.csv |
| A162 | V162.2 (coverage at c_anchor 0.6) | method | P | 16 of 32 cells below criterion; median 0.898 | - | FAIL | 0.90 | A154_gate_calibration/output/a162/summary.csv |
| A163 | V163.1 (coverage, c_anchor 0.8) | method | P | 30/32 cells (min 0.854) | - | FAIL | - | A163_reliability_correction/output/validate/summary.csv |
| A163 | V163.2 (coverage, c_anchor 0.6) | method | P | 27/32 cells (min 0.802) | - | FAIL | - | A163_reliability_correction/output/validate/summary.csv |
| A163 | V163.3 (bias) | method | P | 58/64 cells; max /bias/ 0.146 | - | FAIL | - | A163_reliability_correction/output/validate/summary.csv |
| A163 | V163.4 | method | P | 16/16 (min 0.940) | - | PASS | - | A163_reliability_correction/output/validate/summary.csv |
| A153 | H153.1 (effect, Holm over all primaries) | L2 | P | 0 primaries in either scope | - | not scorable (no eligible primary; every eligible MDE on d* > 0.2) | smallest MDE 0.214 | A153_corrected_contrasts/output/family_eligibility.csv; A153_corrected_contrasts/output/mde_table.csv |
| A153 | H153.2 (equivalence, TOST) | L2 | P | 0 primaries | - | not scorable | bound /d*/ 0.15 | A153_corrected_contrasts/output/family_eligibility.csv |
| A153 | H153.3 (pooled family equivalence) | L2 | P | < 3 targets | - | not scorable | bound /d*/ 0.15 | A153_corrected_contrasts/output/h1533_all5.csv |
| A159 | H159.1-H159.2 (2 hypotheses; replication of H150.1) | L2 | P | all covered atlases are A150 atlases | - | not scorable | - | A159_specific_coupling/output/h159.csv |
| A164 | H164.1 (T2D equivalence) | L2 | P | module d* +0.019, 90 % [-0.021, 0.060] | - | PASS (within +/-0.15) | bound 0.15; MDE 0.083 | A164_gtex_metabolic_matched/output/a164_tests.csv |
| A164 | H164.2 (obesity equivalence) | L2 | P | module d* +0.004, 90 % [-0.081, 0.089] | - | PASS (within +/-0.15) | bound 0.15; MDE 0.130 | A164_gtex_metabolic_matched/output/a164_tests.csv |
| A164 | H164.3 (effect) | L2 | P | T2D +0.019; obesity +0.004 | 0.568, 0.935 / Holm 1.0 | FAIL (no effect detected) | MDE 0.083 / 0.130 | A164_gtex_metabolic_matched/output/a164_tests.csv |
| A165 | H165.1 (phase / amplitude / noise per compartment, 89 tests) | L2 | P | 0 of 89; smallest Holm p 0.089 | Holm >= 0.089 | FAIL (no effect detected) | phase MDE median 2.20 h | A165_disease_readouts/output/h1651.csv |
| A165 | H165.2 (family damping / noise, 4 tests) | L2 | P | fibroblast_stromal log z2 +0.117; endothelial +0.409; 0 of 4 | Holm 1.0 | FAIL (no effect detected) | amplitude MDE 2.11-fold (median) | A165_disease_readouts/output/h1652.csv |
| A167 | H167.1 (stellate, endothelial, hepatocyte d*) | L2 | P | -0.060 / -0.136 / -0.209 (0 of 3) | 0.429 / 0.144 / 0.034; Holm 0.429 / 0.288 / 0.102 | FAIL (no effect detected) | MDE 0.29 / 0.43 / 0.35 | A167_liver_fibrosis/output/h1671.csv |
| A167 | H167.2 (phase, amplitude, noise, 9 tests) | L2 | P | 0 of 9 | Holm 1.0 (all) | FAIL (no effect detected) | phase 1.5-1.8 h | A167_liver_fibrosis/output/h1672.csv |
| A126 | H126.0 (cycling vs non-cycling coherence, gate) | L1 | P | T0 -0.141 (0.241 vs 0.382) | 0.103 | FAIL | MDE 0.276 | A126_menstrual/output/h126_0_tests.csv |
| A126 | H126.1 (luteal tertile, breast) | L1 | P | T +0.029 | 0.575 / Holm 1.0 | FAIL | MDE 0.335 | A126_menstrual/output/h126_main_tests.csv |
| A126 | H126.2 (core-cue offset vs luteal score) | L1 | S | no per-donor estimand defined | - | not scored (descriptive in A126, not run) | - | RP A126 Deviations |
| A126 | H126.3 (BMI >= 30) | L1 | S | T -0.052 | 0.766 / Holm 1.0 | FAIL | MDE 0.372 | A126_menstrual/output/h126_main_tests.csv |
| A127 | H127.1 (coherence vs proliferation, SCAN-B) | L1 | P | T -0.123; 2.5th pseudo percentile | 0.0045 / Holm 0.0045 | PASS, clock-specific (narrowly) | MDE 0.114 | A127_scanb/output/a127_tests.csv |
| A127 | H127.2 (ER- minus ER+) | L1 | P | T -0.343; 1.5th pseudo percentile | 0.0005 / Holm 0.0010 | PASS, clock-specific | MDE 0.158 | A127_scanb/output/a127_tests.csv |
| A127 | H127.3 (single-cell core-cue offset) | L1 | S | - | - | not scored (not registered: header says 'not registered here') | - | RE A127 header |
| A161 | H161.1 (secretory vs proliferative, bulk endometrium) | L1 | P | T +0.096 (0.568 vs 0.473) | 0.987 (T < 0) | FAIL (opposite sign) | MDE 0.117 | A161_endometrium_phase/output/a161_results.csv |
| A166 | H166.1 (ER- loss, METABRIC) | L1 | P | T -0.229; 1.0th pseudo percentile | 0.0010 / Holm 0.0030 | PASS, clock-specific | MDE 0.131 | A166_metabric/output/a166_tests.csv |
| A166 | H166.2 (ER- within PAM50 x proliferation) | L1 | P | T -0.168; 2.0th pseudo percentile | 0.0235 / Holm 0.047 | PASS, clock-specific | MDE 0.172 | A166_metabric/output/a166_tests.csv |
| A166 | H166.3 (proliferation loss, METABRIC) | L1 | P | T -0.057; 13.0th pseudo percentile | 0.148 / Holm 0.148 | FAIL (MDE 0.137 > SCAN-B effect 0.123) | MDE 0.137 | A166_metabric/output/a166_tests.csv |
| A168 | H168.1 (POETIC difference-in-differences) | L1 | P | T +0.381 | 0.993 / Holm 1.0 | FAIL | MDE 0.431 | A168_endocrine/output/tests.csv |
| A168 | H168.2 (AI-treated change) | L1 | S | T -0.010 | 0.452 / Holm 1.0 | FAIL | MDE 0.196 | A168_endocrine/output/tests.csv |
| A168 | H168.3 (GSE20181 letrozole) | L1 | S | T +0.191 | 0.823 / Holm 1.0 | FAIL | MDE 0.504 | A168_endocrine/output/tests.csv |
| A169 | H169.1 (GSE111563, duration) | L1 | P | T +0.370 | 0.917 / Holm 0.917 | FAIL | MDE 0.653 | A169_endocrine_duration/output/tests.csv |
| A169 | H169.2 (GSE59515) | L1 | S | T -0.049 | 0.340 / Holm 0.681 | FAIL | MDE 0.307 | A169_endocrine_duration/output/tests.csv |
| A170 | H170.1 (ER- malignant cells, single cell) | L1 | P | T -0.231 (0.163 vs 0.394); 12th pseudo percentile | 0.0435 | PASS by p; NOT clock-specific | MDE 0.385 | A170_sc_breast_er/output/tests.csv |
| A171 | H171.1 (LUAD GSE72094) | L1 | P | T -0.242; 5.5th pseudo percentile | 0.0005 / Holm 0.0020 | PASS (not clock-specific, narrowly) | MDE 0.167 | A171_luad_lihc/output/tests.csv |
| A171 | H171.2 (LUAD GSE68465) | L1 | P | T -0.293; 1.0th pseudo percentile | 0.0045 / Holm 0.0135 | PASS, clock-specific | MDE 0.304 | A171_luad_lihc/output/tests.csv |
| A171 | H171.3 (HCC GSE14520) | L1 | P | T -0.121 | 0.175 / Holm 0.351 | FAIL | MDE 0.318 | A171_luad_lihc/output/tests.csv |
| A171 | H171.4 (HCC GSE36376) | L1 | P | T -0.102 | 0.257 / Holm 0.351 | FAIL | MDE 0.363 | A171_luad_lihc/output/tests.csv |
| A172 | H172.1 (mean over 12 unseen TCGA projects) | L1 | P | mean T -0.032; 22.5th pseudo percentile | 0.183 | FAIL | MDE 0.087 | A172_tcga_unseen/output/pooled.csv |
| A172 | H172.2 (per project, Holm) | L1 | S | best STAD -0.305 (Holm 0.180); no project passes | Holm >= 0.180 | FAIL (no project passes) | MDE 0.357 (STAD) | A172_tcga_unseen/output/per_project.csv |
| A173 | H173.1 (gastric GSE62254) | L1 | P | T -0.187; 5.0th pseudo percentile (at the line) | 0.016 / Holm 0.032 | PASS, clock-specific | MDE 0.224 | A173_gastric_lung_histology/output/tests.csv |
| A173 | H173.2 (gastric GSE84437) | L1 | P | T -0.134; 48.0th pseudo percentile | 0.071 / Holm 0.071 | FAIL | MDE 0.229 | A173_gastric_lung_histology/output/tests.csv |
| A173 | H173.3 (gastric GSE15459) | L1 | P | T -0.533; 0.0th pseudo percentile | 0.0005 / Holm 0.0015 | PASS, clock-specific | MDE 0.332 | A173_gastric_lung_histology/output/tests.csv |
| A173 | H173.4 (lung ADC vs SQC) | L1 | P | D -0.230; 6.0th pseudo percentile | 0.058 / Holm 0.058 | FAIL (predicted direction) | MDE 0.355 | A173_gastric_lung_histology/output/tests.csv |
| A174 | V174.1 (V3 heart gate passes, current code) | method | P (confirmation) | gate on 16 patients | pairing 0.0005, acrophase 0.0010 | PASS | - | A174_v3_heart/output/v3_scores.csv |
| A174 | V174.2 (PM - AM host displacement) | method | P (confirmation) | length vs window-label permutation | 0.0004 | PASS | - | A174_v3_heart/output/v3_scores.csv |
| A174 | V174.3 (displacement angle within +/-3 h of 18 h) | method | P (confirmation) | angle 19.73 h, error +1.73 h | - | PASS | - | A174_v3_heart/output/v3_scores.csv |
| A175 | H175.1 (pairing lower in high-proliferation tumours) | L2 | S | D +0.158 (low 0.459, high 0.301); pseudo 98.0th, responsive 99.5th pct | 0.021 (Holm 0.021) | PASS (clock-specific; reading pending spread match, A176) | MDE 0.194 | A175_c1_pairing_strata/output/tests.csv |
| A175 | H175.2 (pairing lower in ER- BRCA; runs only if MDE <= 0.3) | L2 | S | - | - | not scored (MDE 0.408 > 0.3, stop rule) | MDE 0.408 | A175_c1_pairing_strata/output/power.csv |
| A176 | H176.1 (swamped clock: in-plane log-odds < 0 in >= 5 of 7 contrasts) | L1 | P | 7 of 7, Ds -0.25 to -0.86, all below pseudo q05 | Holm <= 0.006 | PASS | MDE 0.16-0.54 | A176_c3_plane_variance/output/tests.csv |
| A176 | H176.2 (disordered clock, per contrast) | L1 | P | 0 of 7 disordered; calibration degenerate | - | not informative (design flaw) | - | A176_c3_plane_variance/output/verdict.csv |
| A178 | V178.1 (leave-one-gene-out robust in 4 breast contrasts) | L1 | S | robust: H127.2, H166.1; not: H127.1 (ARNTL, NPAS2, PER3), H166.2 (ARNTL carries it) | - | FAIL | - | A178_c7_robustness/output/verdict.csv |
| A178 | V178.2 (GSE72094, GSE62254 clock-specific with 1,000 pseudo sets) | L1 | S | 2.8th [1.9, 4.0]; 5.7th [4.3, 7.3] | - | FAIL (GSE62254 at the line) | - | A178_c7_robustness/output/pseudo1000.csv |
| A179 | H179.1 (cell-cycle conditioning shrinks the proliferation loss) | L1 | P | 0 of 5; median shrinkage -0.13 | Holm >= 0.37 | FAIL (resolved negative in 4 of 5) | median MDE 0.103 | A179_c4_cellcycle/output/tests.csv |
| A180 | A19 age TOST (+/-0.0375 per decade) | L2 | S | 90 % CI [-0.147, 0.089] | - | FAIL (not equivalent) | bound 0.0375 | A180_equivalence/output/tost.csv |
| A180 | A165 TOST (phase +/-1.5 h, amplitude 1.5-fold; testable rows) | L3/L2 | S | phase 7 of 8, amplitude 7 of 10 equivalent | - | PASS per row (14 of 18) | 24 phase, 22 amplitude rows not testable | A180_equivalence/output/tost.csv |
| A181 | V181.1 (adult design-matched gate pass rate >= 0.8) | method | P | 0.20 [0.15, 0.26] thinned; 0.21 unthinned | - | FAIL (fetal gate uninformative) | - | A181_fetal_depth/output/v181_1_passrates.csv |
| A181 | H181.2 (fetal responsive-gene shared state) | other | P | whole organs p 0.461, stromal 0.366 | Holm 0.731 | FAIL (descriptive under V181.1) | - | A181_fetal_depth/output/h181_2_responsive.csv |
| A182 | H182.1 GSE39582 colon (blind) | L1 | P | T −0.220, pseudo pct 5.5 | 0.041 / Holm 0.123 | not estimable (MDE > 0.30); not equivalent | MDE 0.305 | C6_blind_types/output/tests.csv |
| A182 | H182.2 GSE17536 colon (blind) | L1 | P | T −0.090 | 0.240 / Holm 0.480 | FAIL (not resolved; not equivalent) | MDE 0.292 | C6_blind_types/output/tests.csv |
| A182 | H182.3 GSE116918 prostate (blind) | L1 | P | T −0.034 | 0.407 / Holm 0.480 | not estimable; not equivalent | MDE 0.358 | C6_blind_types/output/tests.csv |
| A182 | H182.4 GSE157009 lung squamous (blind) | L1 | P | T −0.089 | 0.249 / Holm 0.249 | not estimable; not equivalent | MDE 0.329 | C6_blind_types/output/tests.csv |
| A182 | H182.5 GSE157010 lung squamous (blind) | L1 | P | T −0.158 | 0.107 / Holm 0.214 | not estimable; not equivalent | MDE 0.318 | C6_blind_types/output/tests.csv |
| A183 | H183.1 psoriasis lesion delay (GSE154200) | L3 | P | offset −1.45 h; matched pct 5.5, responsive 0.5 | 0.0001 / Holm 0.0002 | PASS, not clock-specific | MDE 0.81 h | A183_skin_gse154200/output/tests.csv |
| A183 | H183.2 psoriasis minus eczema delay | L3 | P | −1.09 h; matched pct 20, responsive 5.5 | 0.037 / Holm 0.037 | PASS, not clock-specific (NR1D1-dependent) | MDE 0.69 h | A183_skin_gse154200/output/tests.csv |
| A184 | H184.1 cord × maternal blood coherence | L2 | S | 0.193 vs null q95 0.221 | 0.103 | FAIL (not resolved) | - | A184_mother_fetus/output/tests.csv |
| A185 | H185.1 tumour-amplitude-matched pairing D | L2 | S | D_match +0.116 (A175 D 0.158) | 0.064 | FAIL (not resolved) | MDE 0.187 | A185_c1_spread_matched/output/tests.csv |
| A136 | H136.1 GR exposure organ donors > ventilated (gate) | other | P | +0.103 (null SD 0.076) | 0.098 / Holm 0.197 | FAIL (gate; p 0.021 without long-term steroid users) | - | A136_glucocorticoid/output/tests.csv |
| A136 | H136.2 GR dose-response in organ donors | L3 | P | ρ top − bottom −0.194 | 0.185 / Holm 0.197 | FAIL (descriptive under failed gate) | MDE 0.515 | A136_glucocorticoid/output/tests.csv |
| A136 | H136.3 A132b difference without PER1/PER2 | L3 | P | 0.458 (≥ 0.26) | 0.012 / Holm 0.036 | PASS ("not through PER1/PER2") | - | A136_glucocorticoid/output/tests.csv |
| A186 | H186.1 CPTAC replication of H175.1 (pairing lower in high-proliferation tumours) | L2 | P | D +0.113 (low 0.506, high 0.392; TCGA 0.158); pseudo 90.0th, responsive 85.0th pct | 0.017 | PASS (not clock-specific) | MDE 0.135 | A186_cptac_c1/output/tests.csv |
| A187 | H187.1 G2/M − NC clock offset in malignant cells (cell-cycle locking) | L3 | P | - | - | not scored (no study passes the NC gate; MDE not computable) | - | C4_single_cell/A187/output/gates.csv |
| A187 | H187.2 locking larger in high-proliferation tumours | L3 | P | - | - | not scored (H187.1 not estimable) | - | C4_single_cell/A187/output/power.csv |
| A187b | H187.1 with the gate over all NC patients (POST HOC, not counted) | L3 | - | −0.069 h (SRP327758, 12 patients); pseudo 18th, responsive 21st pct | 0.83 | FAIL (no offset; post hoc) | MDE 0.48 h | C4_single_cell/A187/output/tests_b.csv |
| A188 | H188.1 malignant coupling to own stroma, high − low proliferation (d*) | L2 | S | d* −0.166 [−0.541, 0.208] (1 cohort, SRP329496) | - | descriptive only (MDE 0.476 > 0.3) | MDE 0.476 | A188_c5_within_tumour/output/descriptive.csv |
| A188 | H188.2 ER− vs ER+ | L2 | S | - | - | not estimable (no eligible cohort) | - | A188_c5_within_tumour/output/descriptive.csv |

## Archive results cited by a figure (A1-A108c; flagged archive, not rescored here)

Source for the rows with a figure panel: `paper/exploratory/WP0_integrity/output/figure_number_sources.csv` (columns `archive_verdict`, `value_in_file`), with the 2026-10-26/30 re-readings of RP where noted.

| archive amendment | test | statistic and value | p | verdict (archive, with later re-reading) | MDE | figure |
|---|---|---|---|---|---|---|
| A3 | primary: prefrontal amplitude above 60 y | log2 amplitude old vs young -0.532 | 0.022 | primary met | - | Fig 2 |
| A17 | V17.2 (anchor rule on colon injury) | lymphatic EC decouples after DSS | - | V17.2 FAILS (anchors then set by hand, A21 panel restored) | - | Fig 1 |
| A21 | R21.4 (epithelial shift blind) | untimed -8.70 h vs timed -9.40 h; without lymphatic anchor -8.95 vs -9.03 h | - | passes | - | Fig 1 |
| A25 | V25.1 / R25 T1 / T2 / T3 (M3 rerun) | axis error 0.34 / 0.27 / 0.40 h; T1 R 0.804; T2 14 of 17; T3 r 0.275 vs 0.332 | T1 0.0014; T2 0.006; T3 0.176 | V25.1 pass, T1 pass, T2 pass, T3 not supported | - | Fig 1 |
| A25 / V25.5 | breast internal time vs resection time | Spearman 0.56; median error 1.1 h | - | post hoc note, not an analysis plan | - | Fig 1 |
| A36 | baboon offsets | median error 0.455 h; R 0.987 | - | primary JS statistic degenerate (fails as written) | - | Fig 1 |
| A42 / A103 | brain vs periphery r* | brain 0.639 [0.59, 0.71]; sheep 0.625 [0.36, 0.87]; baboon 0.97 | - | H103.1, H103.2 pass; sheep passes, baboon readout 1 fails | - | Fig 2 |
| A85 (corrected) | immune cue-following | immune 0.182 vs structural 0.095 | 0.002 | PASS after the 2026-10-22 correction (A122 then shows it is not cue-specific) | - | Fig 2 |
| A104 | H104.1 / H104.2 (strong clocks decoupled) | rho +0.734 (predicted < 0) | 1.00 (rho < 0) | H104.1 not supported, H104.2 not replicated | - | Fig 2 |
| A108 | H108.1 (embryo late window) | R 0.52; window gate p_acrophase 0.08 | 0.016 | passes narrowly; 'from E15' not established | - | Fig 2 |
| A91 | rhythmic genes follow internal time | liver clock excess 2,691 vs time excess -14.5 | - | pass | - | Fig 3 |
| A94 | H94.1-H94.3 (age effects to 70 y) | no age effect | - | all FAIL (no age effect) | - | Fig 4 |
| A19 (rerun 2026-10-30) | H19.1 age amplitude | gamma -0.029 (SE 0.072) | 0.682 | not supported | MDE 0.201 per decade | Fig 4 |
| A19 (rerun 2026-10-30) | H19.2 immune vs stromal age | +0.000 (SE 0.041) | 0.994 | uninformative (immune arm not eligible) | MDE 0.114 | Fig 4 |
| A54 | psoriatic delay, discovery H54.1-3; held-out H54.1-3 | discovery -1.50 h; held-out -0.93 h | held-out H54.2 0.004; H54.1 0.11 | discovery pass; held-out H54.2 pass, H54.1 and H54.3 fail | - | Fig 5 |
| A5 / A17 / A98 | fibrotic-lung T cells | HLCA 0.14 -> 0.00; Natri 0.45 -> 0.00 (free fit; WP0 Q3 anchored: |z| 0.659 -> 0.000 and 0.344 -> 0.000) | 0.027; 0.009; Fisher 0.0023 | holds, conditional on procurement (A119 not registered) | - | Fig 5 |
| A189 | H189.1 L1 coherence, fetal adrenal (Cao) | L1 | P | R 0.734; matched 98.5th, responsive 99.5th pct | 0.0055 / Holm 0.0495 | PASS (clock-specific) | sens 0.705 | A189_fetal_within/output/h189_1_coherence.csv |
| A189 | H189.1 L1 coherence, 8 other fetal compartments | L1 | P | R 0.189-0.422 | Holm 1.0 | FAIL (uninformative) | sens 0.51-0.935 | A189_fetal_within/output/h189_1_coherence.csv |
| A189 | H189.2 within fetal kidney gate | L2 | P | p_acrophase 0.858 | p_pairing 0.920 | FAIL (uninformative) | sens 0.76 | A189_fetal_within/output/h189_2_3_gates.csv |
| A189 | H189.3 within fetal heart / lung gates (A107b re-read, not counted) | L2 | S | lung pseudo 22.5th pct | heart 0.095; lung 0.028 | heart FAIL; lung PASS not clock-specific | sens 0.94 / - | A189_fetal_within/output/h189_2_3_gates.csv |
| A190 | H190.1 within fetal adrenal gate (parenchyma, stromal, endothelial; Cao) | L2 | P | p_acrophase 0.0018; matched 4.5 % (q05 0.0060), responsive 0.5 % | p_pairing 0.0047 | PASS (clock-specific, narrow); confounded: shallow classes ~1/3 cortex RNA, corrected gate p 0.672 (post hoc) | proxy sens 0.875-0.90 | A190_fetal_adrenal/output/h190_1_gate.csv |
| A190 | L1 per adrenal class (descriptive) | L1 | D | R parenchyma 0.731, stromal 0.592, endothelial 0.555 | 0.0065 / 0.054 / 0.084 | descriptive | - | A190_fetal_adrenal/output/l1_classes.csv |
| A191 | H191.1 METABRIC ER+ endocrine tumour clock signal vs RFS (adjusted) | L1 | P | HR 0.977 [0.881, 1.083] per SD; pseudo 29th pct | 0.657 / Holm 1.0 | FAIL | MDE HR 1.15 | C8_internal_time_outcomes/output/tests.csv |
| A191 | H191.2 GSE14520 tumour-liver phase distance vs recurrence (adjusted) | L3 | P | HR 0.956 [0.790, 1.157] per SD; pseudo 43rd pct | 0.643 / Holm 1.0 | FAIL | MDE HR 1.31 | C8_internal_time_outcomes/output/tests.csv |
| A192 | H192.2 Cadenas clock-gene score vs DMFS: proliferation-matched specificity and G2M-adjusted HR (3 fresh untreated N0 cohorts) | L1 | P | unadj HR 0.809 [0.709, 0.923], 7.1th pct of prolif-matched sets; adj HR 0.883 [0.756, 1.031] | 0.115 / Holm 0.136 | carried by proliferation (prediction met) | MDE HR 1.25 | C8b_clock_prognosis_replication/output/a192_tests.csv |
| A192 | H192.3 in-plane log-odds vs DMFS, adjusted (fresh test of A191) | L1 | P | HR 1.145 [0.990, 1.324]; pseudo 90.5th pct; pooled with A191 1.030 [0.947, 1.121], TOST(1.15) p 0.005 | 0.068 / Holm 0.136 | FAIL | MDE HR 1.21 | C8b_clock_prognosis_replication/output/a192_tests.csv |
| A194 | H194.1 TimeTeller Theta vs OS, clinical model (SCAN-B + TCGA) | L1 | P | HR 0.966 [0.879, 1.062] per SD | 0.475 / Holm 0.732 | FAIL | MDE HR 1.15 | C8b_clock_prognosis_replication/output/a194_tests.csv |
| A194 | H194.2 TimeTeller Theta vs OS, clinical + G2M | L1 | P | HR 0.957 [0.870, 1.053] per SD; grade 1-2 secondary 0.833 (opposite to published) | 0.366 / Holm 0.732 | FAIL | MDE HR 1.15 | C8b_clock_prognosis_replication/output/a194_tests.csv |
| A195 | H195.1 cycling − NC clock-plane variance (malignant single cell) | L1 | P | not run | — | NOT ESTIMABLE (signal check: NC coherence R 0.243, 42nd pseudo pctile) | MDE 0.158 | C4_single_cell/A195/output/signal_check.csv |
| A195 | H195.2 NC clock-plane variance, high − low cycling-fraction tertile | L1 | P | not run | — | NOT ESTIMABLE (signal check fails; MDE 0.630 > 0.40) | MDE 0.630 | C4_single_cell/A195/output/power.csv |
| A196 | H196.1 NeoPalAna coherence C1D15 − C1D1 (palbociclib arrest) | L1 | P | T −0.139 (prolif −0.54 SD) | 0.485 | NOT ESTIMABLE (descriptive) | MDE 0.478 | A196_proliferation_arrest/output/tests.csv |
| A196 | H196.2 reversal calibration, 4 series (post hoc) | L1 | S | T − P −0.019; pooled change +0.019 (SD 0.062, 319 pairs) | 0.817 | UNINFORMATIVE (no pre-treatment slope, p 0.41) | — | A196_proliferation_arrest/output/calibration_posthoc.csv |
| A197 | H197.1 SCAN-B coherence, high − low ES-like within PAM50 × proliferation | L1 | P | T +0.092; pseudo 98.0th pct | 0.979 / Holm 0.979 | FAIL | MDE 0.114 | A197_dedifferentiation/output/tests.csv |
| A197 | H197.2 METABRIC, same | L1 | P | T −0.067; pseudo 13.5th pct | 0.100 / Holm 0.200 | FAIL | MDE 0.131 | A197_dedifferentiation/output/tests.csv |
| A197 | H197.3-4 LUAD GSE72094, GSE68465 | L1 | S | T −0.090, +0.008 | 0.100, 0.546 | NOT ESTIMABLE | MDE 0.169, 0.193 | A197_dedifferentiation/output/tests.csv |
