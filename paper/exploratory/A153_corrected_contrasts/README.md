# A153: corrected disease contrasts

Registration: `paper/docs/registered_experiments.md`, A153 and its note of 2026-10-30 (PI Q15 a).

| file | role |
|---|---|
| `scripts/family_eligibility.py` | family screen (output `output/family_eligibility.csv`, scopes `all5`, `rowcheck3`) |
| `scripts/compute_a153.py` | per-target gate and contrast, MDE on d* (resumable: `output/targets/`, `output/mde/`) |
| `scripts/score_a153.py` | primaries, H153.1-3, readout, comparison with A17 |

Outputs (git-ignored, `output/`): `readout_<scope>.csv`, `summary_<scope>.csv`, `h1531_<scope>.csv`, `h1533_<scope>.csv`,
`a17_comparison_<scope>.csv`, `mde_table.csv`, `target_summary.csv`. Scope `all5` is registered; `rowcheck3` is a
descriptive sensitivity.

## Choices not in the registration

1. **Brain targets** (18 of the 36: all PsychAD, SEA-AD, PD) are not run: not in the family, Extended Data (A141).
   They appear as one "not eligible" row each (`brain_targets.csv`).
2. **Anchors:** the learned rule of `anchors_a17.csv`; the three IBD contrasts have no learned anchor set
   ("no anchor"), so they are not estimable.
3. **Seeds:** crc32 of the target name (plus `--seed`, default 0), for gate, contrast and MDE.
4. **Resampling counts** as in `compute_contrasts_a17.py` (1,000 permutations, 200 bootstrap draws, 1,000 for the pooled
   and per-arm gates); clocksync warns that the paper setting for the bootstrap is 500.
5. **Companion:** computed here from the saved permutation draws (`Contrast.null`) as d* (and d) of the compartment minus
   the mean of the fixed-rule stromal-vascular anchors, for every compartment. A compartment that is itself a member
   of that set is read against the other members (leave-self-out; clocksync's companion table leaves the members out).
   The specificity label uses the companion p on d* (`p_companion_d_star`; the p on d is also stored).
6. **MDE:** `power.mde` is called directly, with the permutation blocks of the real contrast passed (`Dataset.contrast(mde=True)`
   does not pass donor-table blocks). The conversion to d* follows the registration: s* by linear interpolation of
   power on the damping grid (first crossing from s = 1 down), MDE = c0 - s*^2 c0 / (s*^2 c0 + 1 - c0). c0 is
   `c_pooled` of `power.mde`: the pooled coupling is read at the pooled V, which is the common spread V0 of the
   contrast, so c*0 = c_pooled. Power never reaching 0.8 gives MDE = inf ("> max"), which is not eligible.
   Replicates with a rank-deficient case group are out of the power denominator (clocksync).
7. **Estimable compartment:** contrast status "ok" with finite p* and se_d_star (a row with both arms at the boundary is
   "not estimable"; a one-arm boundary stays, with `boundary_type_*` reported).
8. **Class of a row:** not estimable (target: no anchor, pooled gate or per-arm gate failed; compartment: both arms at
   the boundary or non-finite estimate) > not eligible (family, or MDE > 0.2, or MDE not computed) > primary, then
   effect (Holm p* < 0.05 over all primaries of all targets, scope-wise), equivalent (90 % interval inside +-0.15),
   inconclusive.
9. **H153.3:** the per-target inverse-variance mean uses the d* of the primaries; the 90 % HK interval is
   mean +- t(0.95, k-1) se_random; effect if the 95 % HK interval excludes 0 (checked first). Targets are not independent
   (shared control arms: dcm / dcm_healthy, skin_ad / skin_ad_nonlesional, endometriosis / severe, lung contrasts); no
   adjustment.
10. **A17 comparison:** A17 verdict from `a17_contrast.csv` (learned rule): "effect" = A17 primary with Holm < 0.05;
    "not estimable" = gate failed, no anchor, both at the boundary. A row is "changed" when its A17 category (effect /
    estimable / not estimable) differs from the A153 one (effect / not estimable / any other class).
11. **D159.1 flag** (descriptive): `clock_specific` of the healthy-atlas compartment with the same name (lung_hlca for
    lung targets, heart for heart targets), empty elsewhere.
12. **Pseudo-clock companion (registered step 3):** not available. No disease control arm has a genome-wide table;
    columns `pseudo_clock` say so. The pseudo-d* quantile in the specificity qualifier is therefore missing.
13. **Data:** every pseudobulk table of the 21 datasets used carries a rule stamp equal to the current builder hash
    (built 2026-10-28 17:03).
