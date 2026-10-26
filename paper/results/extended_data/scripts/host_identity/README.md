# host_identity: coherence, anchors, host identity, A12 and A22

Which compartments carry the host clock state, and is that state the same whatever the anchors? Model-free pair coherence of
every healthy fit and contrast table (A17 steps 0-1), the leave-one-study-out anchor rule (A17 step 2; re-derived here and
compared with the frozen table `common/anchors_a17.csv`), the identity diagnostics of the host (A17 step 3), the common-axis
check of coupling contrasts (A12) and the group-wise anchor confirmation with its revised detector (A22). Everything runs on
`clocksync.Dataset` (`gate`, `coherence`, `fit`, `contrast`, `resample`) plus `power.damped` and `stats.meta`. No plotting
code. Outputs: `paper/results/extended_data/output/host_identity/` (`.../subsample/` for checks).

```
anchor_rule.py                  A17 families and the leave-one-study-out rule (local: not in clocksync)
compute_coherence.py            gate and coherence of the 36 healthy and contrast-table targets (and 34 healthy fits)
compute_anchors.py              anchors learned from those tables vs the frozen table
compute_host_identity.py        agreement free vs anchored, stability, learned vs fixed host
compute_common_axis.py          A12
compute_axis_reading.py         A12 reading rule: separate-fit differences (A5-A9 style) against the common-axis ones
compute_cortex_modules.py       A18 part 3.5: neuron-host vs non-neuronal-host similarity per cortex dataset
compute_group_anchors.py        A22 (V17.1 null, Part 1 diagnostic, revised detector, V17.3 simulation)
run.sh [--subsample N]          in this order (coherence first: later scripts read its tables from the same folder)
```
Helpers (target tables, similarity, stability) are in `../figure_S05_atlas_fits/_shared.py`.

## Analyses
| script | amendment, hypothesis | estimand | legacy script replaced |
|---|---|---|---|
| `compute_coherence.py` | A17 steps 0-1 | gate (1,000 permutations, healthy fits); pair coherence rho_kl (|complex correlation| over shared donors, within strata), excess over 200 shuffled pairings, class excess (mean over partners) with a donor-bootstrap SE (200); contrast tables: pairs on the pooled donors, no bootstrap | `desync/common/anchors.py` (step coherence) |
| `compute_anchors.py` | A17 step 2 | anchors per target: family excess pooled over gate-passing healthy fits of other studies (same tissue group if available) by DerSimonian-Laird, qualifying families (lower 95 % limit > 0) in the top half plus lower ones until >= 3 compartments, confirmation (an anchor with pair excess <= 0 with every other anchor is dropped), < 3: no anchor; fixed rule as sensitivity | `anchors.py` (step select), `src/sccc/atlas/anchor_rule.py` |
| `compute_host_identity.py` | A17 step 3 identity diagnostics | similarity = |complex correlation| of centred host states; agreement (free host vs anchored host), stability (share of 200 bootstrap refits with similarity > 0.7; < 0.8 "unstable"), learned vs fixed host | `run_anchored.py` (step healthy: `agreement_free`, `stability`) |
| `compute_common_axis.py` | A12 (reading rule stated there) | the two groups pooled, one host fitted without the tested compartment k; c_gk per group on that common state; difference case - control, permutation p (within stratum or free, 1,000), conditional donor-bootstrap interval (500); 16 contrasts (A5-A9, A11, A13, A14, A11 part 3) | `run_common_axis.py`, `common_axis.py` |
| `compute_group_anchors.py` | A22 (V17.1, Part 1, note of 2026-10-13, V17.3) | group confirmation by anchor-only pair coherence in each group (delta_k = e_case - e_control, drop if delta < -1.96 SE and e_case <= 0, or the mirror); V17.1 null on 20 random splits of the controls of the 22 A20 contrasts (pass <= 5 % dropped); revised detector (own d* with p* < 0.05 and |d*| >= 0.2 on the leave-one-out pooled host); V17.3 (two anchors per contrast damped to s = 0 and 0.5 in cases by thinning, 3 replicates, detection = p* < 0.05 and d* <= -0.2) | `a22_group_anchors.py` (null, diagnostic, flags, sim) |

## Outputs
`a17_gates.csv`, `a17_coherence_classes.csv`, `a17_coherence_pairs.csv` (legacy columns), `anchors_a17_relearned.csv` (columns of
`common/anchors_a17.csv`), `a17_anchor_families.csv`, `anchors_vs_frozen.csv`, `a17_identity.csv`, `common_axis.csv` (legacy
columns plus `d_star`, `p_star`, `status`), `a22_{null,diagnostic,flags,sim}.csv` (legacy columns).

## Checks
`run.sh --subsample 30` runs every script (22 min; 30 donors, 20 bootstraps, 50 permutations). Full data against
`legacy/desync/common_output/`:

| readout | here | legacy | note |
|---|---|---|---|
| gate, 34 healthy fits: clock covariance, p_pairing, p_acrophase, verdict | identical (max difference 7e-15, 0, 2e-3 (one p at 1,000 perms)); 24 pass | 24 pass | |
| class coherence excess and SE, 361 classes | identical for 33 of 34 fits (max 0) | | tabula_sapiens differs by up to 0.014 / 0.038 (one compartment absent from the frozen table shifts the shuffle stream) |
| pair excess, 4,334 pairs | identical for 66 of 70 tables | | pd, skin_ssc, skin_ad_nonlesional, tabula_sapiens differ by up to 0.07 (same reason; rho identical) |
| anchors re-learned from these tables | learned and fixed rule identical to the frozen table for 70 / 70 targets (34 healthy, 36 contrast tables) | frozen | the cancer targets (21) are not re-learned |
| host identity, 5 fits (hbcc, seaad_dlpfc, heart, kidney_sc, islets): agreement free vs learned | 0.603, 0.796, 0.969, 0.947, 0.994 | identical | |
| stability (60 refits, not 200) | 1.0 in all fits, both rules | 1.0 (learned: hbcc, seaad_dlpfc, heart; fixed: hbcc, heart) | agreement of the fixed rule 0.994, 0.997, 0.991 (hbcc, seaad_dlpfc, heart) |
| A12 common axis, 53 compartment rows | correlation of the differences 0.76, same sign in 85 %; cortex contrasts (DLBD, DM) within 0.01-0.03 | `common_axis.csv` | see deviations: compartment sets differ in the other contrasts |
| A22 revised detector, 4 contrasts (21 anchors; 1,000 permutations) | d* identical (0.0), p* within Monte Carlo (max 0.064); same drop flags; ild_natri capillary EC +0.614 dropped | +0.61 (A22 note) | |
| A22 Part 1 diagnostic, 4 contrasts, both rules (34 rows) | e_case, delta identical (0.0); same flags | | same seeds |
| A22 V17.1 null (22 contrasts x 20 splits), V17.3 simulation | not rerun (checked on the subsample only) | 1.6 % dropped of 2,180 | long run: about 440 coherence jobs |

## M3 full run (2026-10-26)
`run.sh --jobs 4`, 55 min, peak 0.7 GB. Against the legacy full outputs (`legacy/desync/common_output/`) and the archive
A22 section: `a22_diagnostic` (189 rows) identical (max difference 0.0; no anchor dropped under either rule), `a22_flags`
(109 rows) identical except p* (Monte Carlo, max 0.064): the same five anchors are dropped (ild_natri capillary EC +0.614,
endometriosis endothelium +0.405, pd DMV neurons -0.467, SCZ_hbcc OPC -0.371, BD_hbcc oligodendrocytes -0.284; only
endometriosis falls to 2 anchors), as in the archive. `a17_identity` grows from 10 (subsample) to 68 rows; the 22 learned
and 13 fixed fits shared with legacy agree in agreement-free-vs-anchored to 2.2e-6 (tabula_sapiens: 0.404 vs 0.390), stability
differs by up to 0.04 apart from tabula_sapiens (0.400 vs 0.745; cold starts): organ_colon (learned) is now 0.795, below the
0.8 limit (legacy 0.835), so 2 learned fits (organ_colon, tabula_sapiens) are flagged unstable instead of 1. `a22_null`:
35 of 2,180 split x anchor dropped (1.6 %, as in the archive; limit 5 %); endometriosis and ia_crohn draw other pseudo-splits
(the legacy label order was unsorted; 9 vs 6 and 1 vs 4 dropped), the other 20 contrasts are identical. `a22_sim` (V17.3):
detected 114 / 132 at s = 0 (86 %, archive 86 %) and 51 / 132 at s = 0.5 (39 %, archive 42 %); legacy seeds came from
`hash(name)`, here from the contrast index. The 2026-10-23 empty-stratum erratum does not apply: no A22 table has an empty
stratum. `a18_modules` within 8e-10; the other tables identical to the pre-M3 copies. No verdict changes.
<!-- src: extended_data/output/host_identity/{a17_identity,a22_diagnostic,a22_flags,a22_null,a22_sim}.csv vs
legacy/desync/common_output/{a17_unified_healthy,a22_diagnostic,a22_flags,a22_null,a22_sim}.csv (M3, 2026-10-26) -->
The Checks rows "A22 V17.1 null not rerun" and the stability row (60 refits) are superseded by this section.

## Deviations
- **Contrast tables** are the standard tables restricted to the labelled donors and to the compartments of the frozen anchor
  table. Seven Parkinson compartments and one each in skin_ssc, skin_ad_nonlesional and Tabula Sapiens have no pair and are
  not in the frozen table; their absence changes only the shuffle stream.
- **A12 tables.** The legacy A12 loaders used their own compartment sets (classes with >= 20 donors per contrast); here the
  A17 target tables are used (all compartments of the A17 fit are host compartments). The pooled host differs, so
  contrasts with a different compartment set move (natri T: common-axis difference -0.46 here, -0.90 there; the host-identity
  sensitivity A12 describes). The reading rule needs the separate-fit differences of A5-A9 and is not applied here.
  `contrast` returns the same estimand for A17 anchored contrasts; `d_star` and `p_star` are added.
- **A22 V17.2** (A21 mouse colon: lymphatic EC) and Part 2 (A20 inputs on the group-confirmed host, A20 model refit) belong
  to `species` and `resilience`; not here.
- **Stability** is run without the legacy warm start, so the refits are cold starts.
- **Seeds.** Coherence and gate seeds are the registry positions (as legacy); bootstrap intervals elsewhere differ by Monte
  Carlo error.

## Checks of the added items
| readout | here | legacy | note |
|---|---|---|---|
| A18-3.5 similarity neuron vs non-neuronal host: PsychAD MSSM, HBCC, SEA-AD DLPFC, MTG, Parkinson controls | 0.395, 0.573, 0.752, 0.644, 0.242 | identical (`a18_modules.csv`) | all 85 rows identical to 3 decimals; within-module reference 0.70-0.90 (neurons), 0.28-0.93 (non-neuronal); donors per row 257 / 145 / 44 / 47 / 24 here (legacy: donors with all compartments of the half) |
| A12 reading (53 contrast x compartment rows, 100 bootstraps per group) | robust 5, robust null 44, identity artefact 1, non-estimable 3 | A12: lung (Natri) T decoupling robust, kidney T "increase" an artefact | the reading follows the registered rule; values differ from the A5-A9 originals (separate fits on A17 learned anchors: Natri lung T 0.42 vs 0.01, here p-less, legacy 1.00 vs 0.10) |

Deviation: the separate-fit differences of A12 are recomputed with `Dataset.subset(donors).fit` and the table's own host, not the A5-A9 reference classes; the reading rule is applied to those and to `common_axis.csv`.

## Local helpers (not in clocksync)
- `anchor_rule.py`: families, `select`, `fixed` (the learned rule; the package ships only `suggest_anchors`).
- `_shared.similarity`, `stability`, `agreement_free`, `host_of`: built on `Dataset.fit(...).host` and `Dataset.resample`.

## API gaps (for clocksync)
- Learned anchors (`select`) and the family map; `stability`, `similarity`, free-vs-anchored `agreement` as public calls.
- `contrast(..., targets=)` for anchors uses the leave-one-out pooled host, which is what A22 needs; a per-group
  `z2`, `psi`, `V` table would let the per-cell analysis (figure_4) use the pooled host.
