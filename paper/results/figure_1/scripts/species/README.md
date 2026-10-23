# species: known-time validations in mouse, sheep and baboon

Mouse colon injury with known ZT (A21), SheepGTEx (A34), baboon around the clock (A36, A40, A41). Brain vs periphery in sheep and
baboon (A42) is in `organs`; mouse ageing (A48) is `resilience/compute_mouse_ageing.py`; mouse skin (A72) is here. Everything runs on `clocksync`
(`Dataset.fit` / `contrast` / `timed` / `gate` / `subset` / `resample`, `timed_agreement`, `waveform.fit_waveform`) plus `common/` and the builders
`data/{colon_chrono,colon_chrono_raw,sheep,baboon}`; GTEx gains through `common/gtex.py` (aggregates only). No plotting code.
Outputs: `paper/results/figure_1/output/species/` (`.../subsample/` for checks).

```
run.sh [--subsample N]   all scripts (full run by default)
_shared.py               local helpers: known-host regression, Poisson-weighted timed fit, GTEx timed gains, per-tissue timed gain
```

## Analyses

| script | amendment | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_colon_injury.py` | A21: R21.1-R21.4, power | untimed vs known-time truths in C-DC / R-DC / R-PC: anchored coupling vs coupling to time; d*, log z2 and log psi ratios of R-DC minus C-DC (`Dataset.contrast`, permutation within ZT blocks); immune verdicts; epithelial phase shift; MDE; `--counts raw` sensitivity | `21_colon_injury/scripts/colon_injury.py` |
| `compute_sheep.py` | A34 | gate; all-tissue host offsets with bootstrap; circular correlation of sheep offsets with the human timed phase over 21 matched tissues (primary), secondary matches, GTEx positive control, amplitude Spearman | `farm_gtex/scripts/{fit_sheep,human_truth,match}.py` |
| `compute_baboon.py` | A36, A40 (readout 1), A41 | host phase vs ZT, offsets vs timed phase, baboon vs human gain phase over 25 tissue pairs; working-hours windows; baboon waveform vs human | `baboon/scripts/{score,a40_window,a41_reference,human_gain,common}.py` |
| `compute_baboon_stability.py` | A36 (stability, gate, leave-one-animal-out) | stability of the host over 200 donor bootstraps; descriptive gate of the 63 tissues; similarity of the host without each animal | `baboon/scripts/score.py` (stability, gate), `stability.csv` |
| `compute_sheep_study_gates.py` | A34 per-study gates | model-free gate per BioProject (>= 10 animals, >= 2 tissues) | `farm_gtex/scripts/fit_sheep.py` (`study_gates`) |
| `compute_mouse_skin_timed.py` | A72 (H72.1; post hoc dataset permutation) | epidermal stem cells vs whole epidermis in six timed mouse datasets: clock offset (h) of the clock-plane state against ZT, bootstrap 1,000 | `paired_damping/scripts/a72_mouse_timed.py` |

## Outputs
`a21_{r211_couplings,r211_summary,r212_contrast,r212_recovery,r213_immune,r214_phase,mde}[_raw].csv`; `a34_{gate,offsets,primary,human_gain,human_control_offsets,amplitude}.csv`
(`--tissue-gates`: `a34_gate_tissue.csv`); `a36_{readouts,tissues,pairs,stability}.csv`, `a40_{windows,summary}.csv`, `a41_{readouts,genes}.csv`.

## Checks
Every script runs with `--subsample` (A21: reduces bootstraps and permutations only; 24 mice). Full point estimates vs legacy outputs:

| readout | result here | legacy | note |
|---|---|---|---|
| A21 R21.1 anchored coupling vs truth, C-DC / R-DC / R-PC: Spearman, mean error | 0.945, -0.001 / 0.673, -0.040 / 0.958, 0.122 | identical | truth per compartment identical (13 classes) |
| A21 R21.2 Spearman d* / log z2 / log psi | 0.964 / 0.764 / 0.483 | 0.964 / 0.600 / 0.434 | log ratios differ in which compartments are estimable |
| A21 R21.4 epithelial shift untimed / timed (h) | -8.70 / -9.40 | -8.70 / -9.40 | per class identical; 20 bootstraps in the check |
| A21 R21.3 verdicts | all indeterminate | all indeterminate | |
| A34 primary r (n 21), p; secondary (n 25) | -0.199, 0.81; -0.188 | -0.194, 0.80; -0.186 | gate share 0.26 vs null 0.05; jackknife -0.27 to -0.12 |
| A34 positive control (matched 21; all 27) | 0.930; 0.915 | 0.931; 0.934 | p 1e-4; amplitude Spearman 0.575 (p 0.006), legacy 0.56 |
| A36 R1: R, offset, median error, Fisher-Lee r | 0.987, 7.31 h, 0.39 h, 0.980 | 0.987, 7.31 h, 0.39 h, 0.980 | peripheral 0.990 / neural 0.940 identical; JS r 0.68 vs -0.87 (host not centred) |
| A36 R2 untimed offset vs timed, Fisher-Lee r | 0.912 | 0.859 | offsets read on the full host (see deviations) |
| A36 R3 baboon vs human gain phase (25 pairs): r_fl, p | 0.175, 0.199 | 0.177, 0.195 | |
| A40 8 h / 10 h windows, median Fisher-Lee r | 0.552 / 0.557 | 0.516 / 0.55 | 4 h, 6 h infeasible in both |
| A36 host stability (200 bootstraps), median similarity; gate clock share vs null, p; leave-one-animal-out similarity (min, median) | 1.0, 0.9985; 0.522 vs 0.068, p 0.005; 0.999, 1.0 | 1.0, 0.9985; 0.522 vs 0.068, 0.005; `stability.csv` min 0.997, median 0.9986 | bootstrap and gate identical; leave-one-animal-out similarity 0.999 / 1.0 vs 0.997 / 0.9986 (cause not traced) |
| A34 per-study gates (10 BioProjects): clock share, p_pairing, p_acrophase, pass | identical in all 10 rows (PRJNA1304012 0.377, PRJNA1198671 0.253, PRJEB19199 0.104 fails) | `sheep_gate_study.csv` | 8 of 10 pass; 2 fail: PRJEB19199 (p_acrophase 0.08), PRJNA542078 (3 tissues, p_pairing 0.27); run note "9/10 large studies" counts large ones |
| A72 H72.1 primary: whole minus stem offset, 95 % interval; strict | +0.533 h [0.300, 0.738], PASS; +0.825 h [0.558, 1.074], PASS; dataset permutation p 0.10 (20 labellings) | identical (`a72_*_verdict.csv`) | post hoc permutation p: the smallest possible is 0.05 for 6 datasets |
| A41 baboon vs human waveform: shift, R, max residual | 17.03 h, 0.918, 3.71 h (NR1D2) | identical | host phase R with baboon waveform 0.993 |

Not rerun in full: A21 bootstraps (500; 20 in the check), A34 bootstrap (`--n-boot 100`; point estimates only), A36 gate (now `compute_baboon_stability.py`).

## Deviations
- **A21 timed truths** use a local known-host regression (`known_host`, the `on_host` formulas); the contrast uses `Dataset.contrast` with random permutations (1000) within ZT blocks instead of the exact 1296 labellings; untimed immune-family intervals use the mean of the members' bounds (draws not exposed).
- **A34** needs `mixed_strata=True` (124 animals span several BioProjects); per-tissue gates by `gate(targets, hosts)` (`--tissue-gates`), per-study gates not run.
- **A36 R2 and A34 offsets** are read on the full host (`Dataset.fit`), the legacy used each tissue's leave-one-out host; correlations are slightly higher.
- **A40** here is readout 1 only; readouts 2-3 and the 4 h / 6 h windows (in-fit loadings) are `windows/compute_baboon_windows.py` (A25 curvature estimator).

## Not ported
Nothing of the work package: A36 stability, A34 per-study gates and A72 are here; A40 readouts 2-3 (ordering, centre, all window sizes) are in `windows/compute_baboon_windows.py`, A41 readouts 4a-4c in `gate`, A48 in figure_5, A42 in organs.

## API gaps
Known-host `on_host` / `Dataset.timed_loadings`; `Dataset.contrast` draws; A25 window / curvature estimator; `Dataset.stability`.
