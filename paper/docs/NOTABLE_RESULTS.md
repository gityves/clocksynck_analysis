# Notable results (unexpected, strong, or likely to need a follow-up)

One entry per result: date, amendment, the result in one or two lines with its numbers, why it is notable, the
follow-up I would suggest, and the source file. Newest first. Numbers come from output files. Full write-ups are in
`results_preliminary.md`.

## 2026-10-28

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
