# A158: why the clocksync coupling drop in hepatocyte Bmal1 KO exceeds the timed oracle's

Registered in `paper/docs/registered_experiments.md` (A158). `scripts/compute_a158.py` (4 workers, 9 min, run directly: under the 20-min rule),
`scripts/lib158.py`. Outputs: `a158_raw.csv`, `a158_bias.csv`, `a158_reproduce.csv`, `run.log`.

Choices beyond the registration:
- **Design: Manella 2021, not Jouffe.** The registration names `jouffe2022_bmal1ko` (whole-body KO, all samples at ZT12). The observed values to reproduce
  (clocksync -0.293, oracle -0.186; D145.2) come from Manella 2021 (hepatocyte KO vs Alb-Cre, AL; liver on WAT + lung; 24 + 24 mice, ZT 0-22 every
  2 h), and a timed oracle needs several ZT. Simulated: the Manella design (real ZTs, genes after QC, depth 3.2e6).
  Needs confirmation by the author; Jouffe would need a different oracle.
- Control true coupling 0.97 (observed Alb-Cre 0.969), control liver z2 1.926 (observed), host amplitude sqrt(14.0) (observed V), anchors z = 1 and
  coupling 0.9 (A145 value, not estimated), simulator defaults for tau, shift_sd, mesor_sd. Local generator (`lib158.sim_table`, a copy of
  `cs.simulate` with the times given).
- Multiplier 1.0 added (no change, descriptive). Two amplitude arms: KO z2 = control, or the observed 0.215 (BLKO fit) so the z2 floor is exercised.
- Change estimate = coupling(KO) - coupling(control) of independent cohorts at the same ZT design (same spread, so equals d*); the timed oracle is
  A145 `two_step` on the true ZT. Anchor effect: liver added as third anchor (`Dataset.fit` reads an anchor on the leave-one-out host).
- True change = population coupling difference (x - 1) c0. Reproducing true change: linear fit of mean estimate on true change over the 5 multipliers.
