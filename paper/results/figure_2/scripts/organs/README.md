# organs: the shared departure between organs, brain vs periphery, ventilated donors

GTEx v10 tissues as compartments of one donor (timed: time of death known), the brain-vs-periphery replications in sheep and
baboon, the ventilator natural experiment and the timed age test in prefrontal cortex. Clock states and hosts come from
`clocksync` (`Dataset.data`, `fit`, `agreement`, `timed`, `gate`, `stats`, `paired.projector / states`); GTEx through
`common/gtex.py`; sheep, baboon through the builders `paper/common/datasets/{sheep,baboon}`. No plotting code.
Outputs: `paper/results/figure_2/output/organs/` (`.../subsample/` for checks). GTEx time of death is controlled access: aggregates only.

```
run.sh [--subsample N]      all scripts in order (full run by default)
_shared.py                  cohorts, Dataset wrapper, timed per-tissue fit, departures, factor analysis, block correlation (local helpers)
```

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_organs_timed.py` | Aim 1 (pre-A1, `01_gtex_organs`), message 2 | per cohort (fast Hardy 1-2, ventilator Hardy 0): departures of every tissue from recorded time; rank of the shared departure (10-fold CV factor analysis vs donors permuted within tissue; all tissues and periphery); gains, shared fraction, ischaemic slope, untimed coupling / amplitude / offset; brain vs periphery r and r*; covariates (age, sex, ischaemic time, Hardy 2) pooled over tissues | `01_gtex_organs/scripts/gtex_organs.py` |
| `compute_brain_artefact.py` | A103: H103.1, H103.2, D1-D3 | is the brain/periphery split a collection or agonal artefact? r* after adjusting departures for the sample's ischaemic time and Hardy class, without PER1/PER2, exact and amplitude departures, Hardy 1 and 2 alone; donor slope on the brain-periphery ischaemic delay; immediate-early score vs brain-specific departure | `brain_artefact/scripts/a103_gtex.py` |
| `compute_brain_periphery_species.py` | A42 readouts 1-3; A103: H103.3 | sheep and baboon: r* of a brain-only and a periphery-only host (`Dataset.agreement`), study overlap, baboon timed departures (descriptive); brain couples less to the peripheral host (Mann-Whitney); coupling ranks vs matched GTEx tissues (Spearman); sheep without PER1/PER2 | `cross_species/scripts/{r1,r23,human_ext,common}.py`, `brain_artefact/scripts/a103_sheep.py` |
| `compute_brain_artefact_d4.py` | A103: D4 (descriptive) | the donor-level brain-minus-periphery slope on the ischaemic delay, and r*, in Hardy 0 (ventilated) donors; writes the feasibility table (brain tissues with >= 40 ventilated donors; needs >= 4) | `a103_gtex.py` (the `tv` block) |
| `compute_species_rstar_full.py` | A42 readout 1; A103: H103.3 (full size) | the sheep and baboon r* with the registered 200 / 100 bootstraps as six parallel, checkpointed tasks. **Documented, not run: hours** (a `--subsample` check runs in minutes) | `r1.py`, `a103_sheep.py` |
| `compute_chen_age.py` | A3 (`05_brain_age_timed`) | prefrontal cortex (BA11, BA47, Chen 2016, time of death known): log2 amplitude ratio old (> 60) / young (< 40), phase difference, per decade, SD of departures; sensitivity `--adjust ph rin pmi` | `05_brain_age_timed/scripts/chen_age.py` |

Part of the amendments of this figure is in `gtex_cues` (A46 ventilation dose, A47, A81/A81b cue reference,
A89-A93, A100-A102b, glucocorticoid classes D6). A27 part 1 (gene-level go / no-go) belongs to `gene_rhythms`; A27 D4 and
D5 (per-gene) are in S07.

## Outputs
`rank_cv, tissues, covariates, covariates_pooled, phase_loadings_q2, brain_periphery, phase_corr_<cohort>, phase_resid_corr_<cohort>` (01),
`a103_gtex_{rstar,delay,ieg,isch_tissue}`, `r1_{sheep,baboon}_readout, r1_study_overlap, r1_baboon_timed, r2_{couplings,readout},
r3_{readout,pairs}, a103_sheep_readout`, `a33_{tissues,summary,couplings}`, `a33_gc_{couplings,summary}`, `a103_d4_feasibility`, `r1_full/`, `r1_{sheep,baboon}_readout_full`, `a103_sheep_readout_full`, `walltime_hosts`, `chen_age[_adj-ph-rin-pmi]_{region,pooled}, chen_age_genes`.
Column names follow the legacy files. Per-donor departures (`departures.csv` in the legacy) are not written (controlled access).
Note: the files in `paper/results/figure_2/output/organs/` used for the checks below were produced with reduced resampling (`--n-boot`, see Checks);
`run.sh` without arguments overwrites them with the full run.

## Checks
Every script runs with `--subsample` (60-200 donors; `compute_organs_timed.py --cohorts ventilator` is slow because of the factor analysis).
Full point estimates against `legacy/desync/**/output/` (full data, bootstraps reduced: `--n-boot 5` for 01, 10 for A103, 2 for A33):

| readout | result here | legacy | note |
|---|---|---|---|
| 01 fast: donors, tissues | 244, 40 | 244, 40 | |
| 01 rank CV gain over q = 0, phase q = 1, 2, 3, 4 | 0.255, 0.310, 0.307, 0.297 | identical | amplitude 0.298, 0.334, 0.327, 0.315 (legacy 0.333 at q = 2); joint identical |
| 01 brain vs periphery, phase (205 donors): r, r* | 0.575, 0.637 | 0.575, 0.639 | reliabilities by random splits (0.916 / 0.888 vs 0.911 / 0.887) |
| 01 brain vs periphery, amplitude: r, r* | 0.680, 0.748 | 0.680, 0.749 | |
| 01 per tissue (40): |a|, phase gain, ischaemic slope | identical (max diff 0, 0.0002, 0) | | shared fraction max diff 0.002 |
| 01 per tissue: untimed log2 amplitude, offset | max diff 0.007 h-units, 0.0002 | | |
| 01 per tissue: coupling | max diff 0.056, median 0.570 vs 0.589 | | deviation (leave-one-out host, below) |
| A103 base r, r* (205 donors) | 0.5748, 0.6391 | 0.5748, 0.6391 | |
| A103 H103.1 adjusted r* | 0.6385 | 0.6400 | PASS in both (bound < 1 with 10 bootstraps) |
| A103 H103.2 no PER1/PER2 r* | 0.5844 | 0.5843 | |
| A103 exact / exact noPER12 / amplitude / Hardy 1 / Hardy 2 r* | 0.6905, 0.6996, 0.7464, 0.4120, 0.6682 | 0.6903, 0.7013, 0.7458, 0.4119, 0.6651 | |
| A103 donor slope (h), 14 genes | -0.0314 | -0.031 | |
| A42 sheep readout 2 (brain vs periphery coupling to the peripheral host) | 0.205 vs 0.416, p 0.0090 | identical | sensitivity 0.0022 identical |
| A42 baboon readout 2 | 0.420 vs 0.712, U 145, p 4e-5 | identical | |
| A42 readout 3 sheep primary | n 21, rho 0.331, p 0.067 | identical | MATCH + SECONDARY 0.399 / 0.0275 identical; non-brain 0.291 / 0.143; extended human fit 0.266 / 0.119 |
| A42 readout 3 baboon secondary | n 25, rho 0.469, p 0.011 | identical (0.011) | peripheral pairs 0.230 / 0.182 |
| A42 baboon timed departures (subsample run: all 12 animals) | r 0.455 (p 0.11), Fisher-Lee 0.455 (0.074), rotated state 0.665, R 0.876 / 0.990, circular SD 1.96 / 0.53 h | 0.455 (0.126), 0.455 (0.064), 0.665, 0.876 / 0.990, 1.96 / 0.53 | permutation p differ by Monte Carlo (500 here) |
| A33 loss without GC: D, p (tissue permutation, 3432 splits) | 0.0724, 0.1361 | 0.0724, 0.1361 | per-tissue losses identical; bootstrap not rerun |
| A33 loss with GC | 0.0741, p_perm 0.131 | 0.0769, 0.136 | approximate GC adjustment (deviation) |
| A103 D4 (Hardy 0): brain tissues with >= 40 ventilated donors | 0 of 4 needed (largest: 16 donors, frontal cortex) | not run (same rule) | infeasible, as legacy: `a103_d4_feasibility.csv` |
| A33 loss with GC, target-specific adjustment: D, p_perm (3432 splits; no bootstrap) | 0.0603, 0.1745 | 0.0769, 0.136 | between the approximate adjustment above (0.0741) and zero; see deviations |
| A3 Chen pooled (500 bootstraps): log2 amplitude old/young, phase, per decade, phase per decade | -0.532, +3.47 h, -0.053, +0.88 h | -0.525, +3.45 h, -0.053, +0.88 h | p 0.022 / 0.002 / 0.38 / 0.002 identical; SD of departures +0.057 vs +0.001 (p 0.83 vs 0.94) |
| A3 adjusted for pH, RIN, PMI | -0.514, +2.67 h | -0.499, +2.66 h | |
| walltime hosts (subsample, 60 donors per cohort) | R fast 0.62, ventilator 0.20 | R 0.63 (A46 context), 0.14 | full run not made |

Not rerun in full: the bootstraps (01, A103, A33), the ventilator cohort of 01, the r* of sheep and baboon (`Dataset.agreement`
with 50 splits x 200 bootstraps is serial: hours) and A103.3. Their subsample runs completed (sheep 200 animals, baboon 12 animals).

## Deviations
- **Untimed coupling in 01.** The legacy coupling came from the free host of all tissues (`factor.fit`); `Dataset.fit` reads an anchor on its
  leave-one-out host (as the A104 example). Amplitude and offset (read on the full host) are unchanged.
- **Likelihood check in the factor analysis** every 5 EM steps (the legacy: every step). Held-out gains agree to 3 decimals.
- **Bootstrap in 01.** `Dataset.resample` replicates (host refitted) instead of the legacy `factor.fit` refit; same design, other random stream.
- **A33 GC adjustment.** The legacy residualised the target and the host on [1, GC of target, GC of host tissues] after the host fit (internal
  `on_host`). Here the state of each tissue is residualised on its own GC score (per cohort) before the host fit (`adjust_gc`): D 0.0741 vs 0.0769.
  The unadjusted readout is exact. Bootstrap intervals reference 1000 replicates by default (as legacy).
- **A33 GC, target-specific (`compute_ventilator_gc.py`).** The legacy residualised the target's state and the host (the 2-column
  host posterior mean, fitted without any covariate) on [1, GC of the target, mean GC of the host tissues]. Here the states of all
  compartments are residualised on the same two covariates before the host is fitted (`Dataset.adjust`), so the host is built from adjusted
  anchors instead of being projected afterwards; D 0.0603 (p_perm 0.17) vs legacy 0.0769 (0.136). The legacy order (host first, then
  coupling given a supplied host) needs a public coupling-on-a-given-host call (API gaps). Bootstrap: `--n-boot` (100 by default, not run in full).
- **A42 sensitivity "no stratum centring"** (point estimate only in the legacy) is not ported: `Dataset.agreement` always centres on the modal stratum.
- **A3 projector.** `paired.projector` centres over genes; the legacy centred each region over samples only (differences of 0.01-0.02 in the ratios).
- **Brain tissues in the ventilator cohort.** Fewer than 4 brain tissues have >= 40 ventilated donors: no brain-periphery row (as legacy); A103 D4 not run.

## Local helpers (not in clocksync), `_shared.py` unless noted
- `timed_fit`: per-tissue Poisson-weighted complex least squares on cos / sin of the recorded time with ML own noise (the package has `timed.cosinor` for real series only) and `departures`.
- `fa`, `rank_cv`, `heldout`, `varimax`: factor analysis with missing entries and its cross-validation.
- `block_stat`, `block_boot`, `reliability`: brain-vs-periphery correlation of departures with split-half disattenuation.
- `sample_ischaemic_h`: ischaemic time of the sample behind each (donor, tissue) row (`gtex.load` keeps one value per donor).
- `adjust_gc` (A33), `timed` (complex least squares with age x e^{it}, `compute_chen_age.py`).

## Not ported
- 01: the `exact/` rerun as a separate output (exact departures are in A103 D1), per-donor `departures.csv`.
- A27 part 1 (go/no-go, `gonogo*`) is in figure_S08; the A27 diagnostics D1, D2 and D4 are in `gtex_cues/compute_gonogo_diagnostics.py`.
- A42: the `r1_baboon_boot.csv` / `r1_sheep_boot.csv` replicate files (`Dataset.agreement` returns the interval only).
- Run in full only as documented scripts: `compute_species_rstar_full.py` (hours), `compute_ventilator_gc.py --n-boot` beyond a handful.

## API gaps
- `Dataset.timed_loadings(hours)`: per-compartment amplitude and offset against known times with the Poisson-weighted fit (now local `timed_fit`).
- `Dataset.fit(..., covariates=...)` / `Dataset.adjust` now exist (A33 uses `adjust` per target); still missing: the coupling of a target on a *supplied* host after both are residualised on target-specific covariates (the legacy order), e.g. `Dataset.fit(host=...)`.
- `Dataset.agreement`: parallel bootstraps and the replicate table; an option without stratum centring.
- `Dataset.resample` replicates should expose the row index (to carry covariates; worked around by attaching columns to `ds.table`).

