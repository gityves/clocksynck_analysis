# A19 acceptance check (PI Q10, post hoc diagnostic, 2026-10-26)

`check_null.py --n-perm 200 --n-sim 400` (1 worker, 105 s). At the six acceptance compartments, on the real tables with
ages permuted within stratum (no age effect), the permutation distribution of each A19 statistic is compared with the
simulated null that `Dataset.covariate` uses for its `_p`. Real ages are not read.

**Result.** The simulated null does not match the permutation null, in two ways:
- **Too narrow (PsychAD, all three compartments).** SD ratio permutation / simulation 1.10-1.28 for locked amplitude,
  gamma and coupling slope (gamma 1.24, 1.20, 1.26 in OPC, inhibitory, excitatory); 8-13 % of permutation statistics fall
  outside the simulated central 95 % (0.05 if calibrated). The real data are more variable than the fitted Gaussian
  model.
- **Off-centre (lung HLCA, all three compartments).** Under permuted ages the statistics are not centred where the
  simulation is: locked amplitude t median 1.84 vs 0.17 (AT1), 0.42 vs -0.02 (multiciliated), -0.53 vs 0.23 (NK/ILC);
  gamma t -1.20 vs 0.26 (multiciliated); 34 % of AT1 locked-amplitude statistics fall outside the simulated 95 %. The
  shift is present with and without boundary fits (AT1 median 1.95 interior, 1.80 at the boundary; multiciliated has no
  boundary fits), so it is not the boundary alone.
- **Degenerate fits (lung AT1, NK/ILC).** The age model puts own noise at its floor (coupling >= 0.999) in 78 % and 84 %
  of permutations (median coupling 1.0; NK/ILC's full-data coupling is 0.045) and gamma t is undefined in 7 % and 48 %.
  The age model's coupling column is not the A17 coupling and is unreliable for weakly coupled compartments.

**Reading.** The `_p` of `Dataset.covariate` is anti-conservative at real designs: mildly where the model fits (brain,
overdispersion) and strongly where it does not (lung atlas). A within-stratum permutation of the covariate, refitting
the statistic each time, is calibrated by construction and cheap (seconds per compartment).
<!-- src: paper/exploratory/A19_acceptance_check/output/{null_compare,perm_fits}.csv -->
