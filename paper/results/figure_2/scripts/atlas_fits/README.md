# atlas_fits: cohort table and the anchored fit of every healthy atlas

The A17 re-analysis of the healthy fits (A1-A16 and the later atlases, 34 datasets; source amendments in
`docs/PLAN_internal_desync.md (archive)`): gate, anchored coupling, relative amplitude and offset of every compartment under the two
frozen anchor sets (learned leave one study out, and the fixed stromal-vascular / neuronal rule), and the cohort table.
Host identity (agreement, stability) and the anchor learning are in `host_identity`; the family summary is in
`backbone`; the disease contrasts of A17 are in figure 5. No plotting code. Outputs:
`paper/results/figure_2/output/atlas_fits/` (`.../subsample/` for checks).

```
_shared.py          A17 target registry (data/_atlas/a17_targets.csv), `Target` (donors, compartments, labels, blocks of one
                    A17 fit), host-identity helpers (imported by figure_S04 and figure_4)
compute_atlas_fits.py
compute_step7_confirmations.py   A17 step 7: Tabula Sapiens pair excess coherence of anchors, GTEx bulk vs snRNA stromal-vascular fraction
run.sh [--subsample N]
```

## Analyses
| script | amendment, hypothesis | estimand | legacy script replaced |
|---|---|---|---|
| `compute_atlas_fits.py` | A17 steps 0, 3, 8 (descriptive; no hypothesis id) | per healthy dataset: gate (1,000 permutations; nothing estimated if it fails), coupling c_k, relative amplitude and offset of every compartment with donor-bootstrap intervals (host refitted, default 500), learned and fixed anchors from `common/anchors_a17.csv`, each anchor on the host of the others; cohort table (donors, compartments, median cells, depth) | `desync/common/run_anchored.py` (step healthy), `anchors.py` (gates) |

## Outputs
`a17_healthy_fits.csv` (legacy `a17_healthy.csv`: target, rule, gate p, status, compartment, coupling, interval, amplitude,
offset_h, psi, V, z2, n, anchor, boundary, host_n, n_boot), `a17_healthy_gates.csv` (legacy `a17_gates.csv`),
`cohort_table.csv`, `a17_healthy_hosts_summary.csv` (host spread per fit; no per-donor states). Column `amplitude` is
the package's reference-free relative amplitude.

## Checks
`run.sh --subsample 30` (17 min; 30 random donors per fit, 100 gate permutations, 20 bootstraps) runs every dataset. Full-data point
estimates (no bootstrap, 200 gate permutations) against `legacy/desync/common_output/a17_healthy.csv` and `a17_gates.csv`,
8 healthy fits (psychad_hbcc, seaad_dlpfc, heart, kidney_sc, organ_colon, islets, endometrium, ia_lung):

| readout | here | legacy | note |
|---|---|---|---|
| gate clock covariance, 8 fits | identical (e.g. psychad_hbcc 3.1784, heart 13.907) | | verdict identical (7 pass, ia_lung fails); p at the permutation floor |
| anchored coupling, learned and fixed, 75 compartments of 7 passing fits | max difference 0.0 | | also psi and V (0.0) |
| rule status | learned: 7 ok; fixed: 3 ok, 4 no anchor | same | fixed rule has < 3 anchors in kidney, islets, colon, endometrium |
| relative amplitude and offset | differ (log-amplitude by 0.04 to 2.2, offsets by up to 4.4 h) | | different reference frame, see deviations |

## Deviations
- **Amplitude and offset frame.** `Dataset.fit` reports amplitude and offset relative to the mean over compartments on one
  host (reference-free); the legacy values are relative to the reference compartment. Couplings, psi and V do not depend on
  this and are identical. Shape of the amplitude and offset differences is not compared further.
- **Compartments of a fit** are those of the frozen anchor table (`common/anchors_a17.csv` lists every compartment that forms
  a pair, i.e. the legacy fit): 33 of 34 healthy tables equal the legacy count; tabula_sapiens has 62 of 63.
- **Not ported:** the simulation lower bound of boundary couplings (`c_lo_bound`, A17 step 4; `lower_bound` is not public),
  MDE and gate power (`power.mde`, `power.gate_power` are in the package but belong to figure 5 and figure_S03), waveform
  sensitivity (A17 step 6) is in `gate/compute_waveform_sensitivity.py`.
- **A17 step 7** (`compute_step7_confirmations.py`): Tabula Sapiens anchor x anchor pairs: excess 0.189 vs 0.042 (diff 0.147 [-0.012, 0.356], 513 pairs, 15 anchor pairs) against legacy 0.184 / 0.041 (0.144 [-0.011, 0.344]); GTEx bulk vs Eraslan snRNA stromal-vascular fraction over 8 tissues: Spearman with `phase_dep_sd` -0.952 (p 0.000), with coupling +0.762 (p 0.028), partial given amplitude -0.940 / +0.701; identical to the legacy tables. Bootstrap over compartments (pair table, not donors).

**M3 full run (2026-10-26, `run.sh --jobs 4`, 2.8 h, peak 2.2 GB; comparison `paper/results/logs/atlas_fits_m3_compare.md`):**
- Gates: 34 of 34 targets as legacy (24 pass), no flips.
- Healthy fits: 478 estimated rows, coupling within 1.6e-6 of legacy, status, boundary and anchor identical; only bootstrap
  intervals differ (new stream; lower bound > 0 flips in 9 boundary rows, none in a registered verdict).
- The consumers compared earlier (host_identity, resilience, spread_multiplicity) are therefore unaffected.
- Fix: the organs table now holds both GTEx cohorts, and step 7 had joined both (16 rows; sv vs coupling +0.527).
  It now reads the fast cohort (rerun 15:35): the GTEx table equals legacy exactly.
- Tabula step 7 unchanged: 0.189 vs 0.042, diff 0.147 [-0.012, 0.356].
<!-- src: paper/results/figure_2/output/atlas_fits/{a17_healthy_gates,a17_healthy_fits,a17_step7_gtex,a17_step7_tabula}.csv -->
- **Seeds.** The gate uses seed = position of the target in the registry (as the legacy coherence step), the bootstrap
  seed 0; the legacy bootstrap seeds were job indices, so intervals agree up to Monte Carlo error only.
- **Full run.** Not rerun here (500 bootstraps on 34 datasets takes hours: about 0.8 s per bootstrap and compartment).

## API gaps (for clocksync)
- `Dataset.fit` could return the legacy-style reference-compartment amplitude and offset as an option.
- Simulation lower bound for boundary couplings (`lower_bound`) and `stability` / `similarity` / free-vs-anchored
  `agreement` of the host are in `_core` only (see figure_S04 for the local helpers).
