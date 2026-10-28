# Notable results (unexpected, strong, or likely to need a follow-up)

One entry per result: date, amendment, the result in one or two lines with its numbers, why it is notable, the
follow-up I would suggest, and the source file. Newest first. Numbers come from output files. Full write-ups are in
`results_preliminary.md`.

## 2026-10-28

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
