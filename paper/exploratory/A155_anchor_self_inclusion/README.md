# A155: bias from reading anchors on a host that contains them

Registered in `paper/docs/registered_experiments.md` (A155). Code: `scripts/compute_a155.py` (simulation and scoring),
`scripts/depth_values.py` (depth values). Outputs in `output/`: `raw/<cell>.csv.gz` (replicate x compartment), `summary_bias.csv`,
`summary_a123like.csv`, `a123like_per_rep.csv`, `decision.csv`, `depth_values.csv`.
Run: `RUN_CAP_GB=12 paper/common/run.sh paper/exploratory/A155_anchor_self_inclusion/scripts/compute_a155.py --jobs 6` (13 min, peak 0.6 GB).

## Choices not in the registration
- **Readouts (a) and (b).** (a) = `Dataset.fit`'s readout: z of every compartment on the full anchor host (`_core.host.host`, `on_host`),
  centred on the mean log amplitude and the circular mean offset over all compartments of the design. (b) = for each anchor, the
  host fitted on the other anchors (`host(..., exclude=k)`, read-only use of `_core`), z read on it; non-anchors on the full host in both.
  Coupling in (a) is read on the full host too (note that `Dataset.fit`'s own coupling column already uses the leave-one-out host).
- **Gauge of the leave-one-out host.** A host is defined up to the loading of its reference compartment (z_ref = 1), so z values
  read on hosts with different references are not comparable. The leave-one-out hosts use the full host's reference
  (`ref=`); when the reference compartment itself is left out, the host's reference is another anchor and z is multiplied by
  the complex scalar s = sum(conj(z_loo) z_full) / sum(|z_loo|^2) over the other anchors (their loadings on the two hosts).
- **Truth.** The 40 A123 tissues (truth (i) amplitude and offset, `A123_working_hours/output/truth_i.csv`; coupling = the A123
  per-tissue estimate of `organ_map.csv`, clipped to [0.15, 0.95]) are permuted once (seed 155): the first 4 are the targets, the
  next m the anchors (nested in m). psi per compartment from c with the sampling spread of the window 8-18 h, jitter 1 h
  (psi = |z|^2 V (1 - c) / (2 c)). The A123-like design uses all 40 tissues as anchors, n = 111.
- **Depth.** `depth_values.csv`: single cell = median clock-gene UMI per pseudobulk row over the 35 A17 single-cell pseudobulk tables
  present in `$SCCC_DATA/pseudobulk` (35,710 rows); bulk = median clock-gene count per GTEx fast-death sample (5,611 rows; aggregate
  only). The simulator's genome-wide `depth` is set so that the expected clock-gene UMI equals that median (`depth_for`).
- **Decision statistic.** Per cell, the anchors' per-compartment bias (mean over 200 replicates) of exp(log |z|) - 1 and of the wrapped
  offset error; the rule is applied to the mean over anchors of the absolute per-compartment bias (the signed mean is also in
  `decision.csv`). Typical designs: m 3-6 x n 30-100, both depths shown (A17 is single cell).
- Simulator defaults otherwise (tau, shift_sd, mesor_sd), `amplitude` 1, one pseudobulk per donor and compartment, all present.
