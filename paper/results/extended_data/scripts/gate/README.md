# gate: gate power, extended waveform, species reference

What a failed gate means, and whether the 14-gene waveform can be improved or replaced. Public `clocksync` (`Dataset.gate`, `power.gate_power`,
`waveform.fit_waveform`, `Dataset.timed / fit`, `timed_agreement`); A25 estimator and baboon loaders from `windows`. No plotting.
Outputs: `paper/results/extended_data/output/gate/` (`.../subsample/`). Helper module `_gate_shared.py`.

| script | amendment | content | legacy |
|---|---|---|---|
| `compute_gate_power.py` | A17 step 5 (gate power) | gate of the 34 healthy fits; power of the failing ones to a lung-sized shared clock | `common/run_anchored.py --step gatepower` |
| `compute_extended_clock.py` | A28 | W14 / W18 / Wcons / Wgen on the GTEx selection half; E1 single sample, E2 host, E3 ventilated vs CHIRAL DIP | `validation/extended_clock/scripts/evaluate.py` |
| `compute_baboon_reference.py` | A41 | baboon vs human waveform: genes, A40 windows paired, full day, sheep gate / couplings | `validation/baboon/scripts/a41_*.py` |
| `compute_waveform_sensitivity.py` | A17 step 6 (descriptive) | healthy fits rerun with each tissue's own GTEx waveform: gate verdicts, coupling correlation | `common/waveform_sens.py` |
| `compute_gene_screen.py` | A28 (screen on the selection half) | per-gene shared acrophase across the 15 GTEx tissues; strict / generous gene sets | `validation/extended_clock/scripts/screen.py` |
| `compute_extended_e4.py` | A28 E4 | W14 / W18 / Wcons / Wgen on Li 2024 breast with the V25.5 estimator | `validation/extended_clock/scripts/evaluate.py` (e4) |
| `compute_sheep_reference_primary.py` | A41 readout 4c | A34 primary (sheep offsets vs human tissue phase) with human and baboon references | `validation/baboon/scripts/a41_sheep.py` |

## Checks
| readout | here | legacy |
|---|---|---|
| A28 E1 (full): median error W14 / W18; change; CI; tissues better | 2.648 / 2.583 h; -2.46 %; [-0.148, 0.082]; 10 | identical |
| A28 E1 Wcons, Wgen change | +2.41 %, -1.43 % | identical |
| A28 E2 15 tissues, R (W14 ... gtex_orig) | 0.59-0.60 | 0.60 |
| A41 genes: shift, R, within 2 h, amp Spearman | 17.0 h, 0.918, 11/14, 0.43 | 17.0, 0.92, 11/14, 0.43 |
| A17 step 6 (20 fits, 161 compartments): gate verdict unchanged; median abs coupling difference | 20 of 20; 0.026 (r 0.927) (rerun on the relearned anchors 2026-10-29; old values in m3_snapshot_backup/2026-10-29) | 20 of 20 (gates, `a17_waveform_sens_gates.csv`); 0.026 | |
| A28 screen: 17,404 genes; amp, z, phase SD, tissues, n_sig | max difference 4e-15, 1e-14, 9e-13, 0, 0 | | I^2 differs for 123 genes expressed in one tissue (0 vs 1 convention); sets W18 4 genes (BHLHE41, HLF, PPL, ADIRF-AS1), Wcons 14, generous 20 |
| A28 E4 Spearman: W14, W18, Wcons, Wgen, original | 0.559, 0.582, 0.700, 0.611, 0.558 | 0.559, 0.582, 0.700, 0.578, 0.558 | MAE 1.20, 1.21, 1.01, 1.07 h vs 1.20, 1.21, 1.01, 1.80 h; Wgen differs (cause not traced); no-harm holds for all |
| A41 4c: JS r of sheep offsets vs human d_t (21 tissues), one-sided p; human / baboon reference | -0.199 (p 0.81) / -0.456 (p 0.99) | -0.194 (0.80) / -0.378 (0.96) | A34 primary fails with both references; baboon waveform refit here |
| gate power (subsample 60 donors, 5 healthy fits, 5 injections) | runs; not comparable | |

## M3 full run (2026-10-26)
`run.sh --jobs 4`, 64 min, peak 1.3 GB. Against the archive copies (`outputs/paper/figure_S03_gate`): every table identical
or within 2e-8, except 5 of 17,404 `gene_screen_genes` rows whose I^2 flips 0 / 1 (genes in one tissue, I^2 undefined).
First full runs:
- **Gate power** (100 injections; `gate_power.csv`) against the legacy `common_output/a17_gatepower.csv`: every reading
  unchanged (blood OneK1K and adipose indeterminate, 0.64 / 0.56 vs 0.72 / 0.64; AIDA, bone marrow, ileum, spleen,
  jejunum low power; thoracic lymph node shared, not clock-shaped, 0.81 vs 0.88). Immune |z| 0.929 here vs 0.991 (A17
  estimator). Target sets differ by one each: SSc control skin (legacy, 0.24, low power) is not in the healthy fits here;
  IA lung is new (0.15, low power).
- **A41 windows** (`baboon_reference_summary.csv`) against the legacy `validation/baboon/output/a41_summary.csv`: m = 5 /
  6 median LOO r_FL human 0.552 / 0.557 (legacy 0.516 / 0.546), baboon 0.475 / 0.577 (0.456 / 0.561); baboon better in
  6 of 12 at both m (as legacy). A41 readouts 1 and the truth check identical; offsets vs timed (human-waveform truth)
  r_FL 0.755 (legacy 0.719; one full host per window here, legacy leave-one-out).
<!-- src: extended_data/output/gate/{gate_power,baboon_reference_summary,baboon_reference_readouts,gene_screen_genes}.csv (M3, 2026-10-26) -->

## Deviations
- Selection half regenerated (seed 0 over all fast-death subjects); the gene screen is rerun in `compute_gene_screen.py` (identical to the legacy table).
- A17 step 6: couplings are point estimates (no bootstrap), gate with the free host and 500 permutations; ileum (no GTEx entry) not rerun.
- A41 4c: the human truth d_t is `a34_human_gain.csv` of `species` (run it first).
- Gate power: Sigma_f of the lung host is `host_cov` of the anchored HLCA fit; |z| = median sqrt(z2) of immune classes.
- A41 offsets are read on one full host per window (legacy leave-one-out).

## Not ported / not rerun
Nothing of the work package. A41 readouts 4a, 4b are in `compute_baboon_reference.py`.
Not rerun in full: gate power (100 injections), A41 windows and sheep.

## API gaps
Same as figure_S02 (Sigma_f; 2-anchor host: E4 uses the local two-compartment estimator); `Dataset.gate_power` takes the dataset's own host spread (here `power.gate_power` with the lung spread).
