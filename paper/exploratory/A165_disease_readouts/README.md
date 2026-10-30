# A165: phase, amplitude and noise readouts of the disease contrasts

Registration: `paper/docs/registered_experiments.md`, A165 (commit eaae08a). clocksync: `Dataset.contrast` columns
`p_log_z2`, `p_log_psi` (commit 9db4ed8 of the clocksync repo).

| file | role |
|---|---|
| `scripts/compute_a165.py` | A153 stage 1 rerun (same targets, seeds, anchors, gates), companion for every readout |
| `scripts/score_a165.py` | H165.1 (Holm over primary rows x readouts), H165.2 (HK pooled, one-sided, Holm over 4), MDEs |

Outputs (git-ignored, `output/`): `targets/`, `target_status.csv`, `readout.csv`, `h1651.csv`, `h1652.csv`,
`gate_share.csv`, `summary.csv`.

## Choices not in the registration
1. **Companion of the phase readout:** the circular difference between the compartment offset and the circular mean of
   the companion offsets, with the permutation p taken without a tie margin (as `p_offset`).
2. **MDE:** 2.80 x bootstrap SE (normal approximation). Natural units: hours (phase), exp(MDE/2) fold in amplitude,
   exp(MDE) fold in psi.
3. **A phase effect whose companion p ≥ 0.05** is classed "effect shared with the companion (not internal
   desynchrony)".
4. **H165.2 per-target mean:** inverse-variance weighted (1/SE²) over the primary rows of the family; SE = mean SE.
5. **Pipeline checks:** a `--subsample 40` run of skin_ad and dcm, and a synthetic check of the H165.2 pooling. These
   are not results; the subsample output is under `output/subsample/`.
