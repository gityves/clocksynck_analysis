# prep_d2_gse27272: mother and fetus at term (HANDOFF D2), prep + gate + power only

Status 2026-11-01: PREP, GATE AND POWER DONE; no observed coupling, coherence or group difference computed, printed or saved. Not registered,
not committed. Script: `scripts/prep_d2_gse27272.py --stage prep|gate|power` (run with `python -I`; each stage < 1 min, no detached run needed).

## Data and design
GSE27272 (Illumina HumanRef-8 v3, GPL6883; Ceske Budejovice, 2008-09; maternal smoking study). 183 samples, 69 deliveries (`individual`):
maternal blood M 65, cord blood D 64, placenta PL 54. Deliveries with >= 2 compartments 67 (20 with 2, 47 with all three); M+D 60, M+PL 51,
D+PL 50; fetal target with maternal host 64. Smokers 20 deliveries (D 18, PL 16 with host), non-smokers 49 (D 42, PL 35).
Design factors in the metadata: maternal smoking (plus cotinine, age, BMI, parity, gestational age 37-41 w). **All 69 deliveries are
vaginal, so there is no labour / delivery-mode factor**, no labour duration, time of day or season. Values are deposited lumi-quantile log2
(each compartment normalised separately; only probes detected at p < 0.01 in all samples kept: 24,526 probes). The platform annotation was
not in the series matrix: `GPL6883_full.txt` (GEO text, sha256 in `fetal_gse27272/PROVENANCE.md`) was fetched; probes -> 18,630 symbols (mean
of log2 per symbol, `CIART` via the `C1orf51` alias).
Clock genes: 13 of 14 present; **PER3 is absent from the platform** (no probe). The GTEx waveform is restricted to the 13 (as A177). Several
genes are expressed near the floor in blood (PER1, NR1D1, NPAS2, CIART mean log2 1.5-2.3: `output/clock_gene_coverage.csv`).
Counts for clocksync: 1000 x 2^x (depth cancels in the centred log fraction; Poisson variance negligible), so the model-free gate is exact.

## Gate (`Dataset.gate`, 2000 permutations; `output/gates.csv`)
| scope | donors | p_pairing | p_acrophase | pass |
|---|---|---|---|---|
| all three compartments, all deliveries | 69 | 0.0005 | 0.0030 | yes |
| fetal (D, PL) x maternal host | 69 | 0.0165 | 0.0210 | yes |
| cord blood x maternal blood | 69 | 0.0220 | 0.0130 | yes |
| placenta x maternal blood | 69 | 0.133 | 0.190 | **no (uninformative)** |
| cord blood x placenta | 69 | 0.0005 | 0.0055 | yes |
| non-smoker arm, all compartments / fetal x host | 49 | 0.088 / 0.134 | 0.120 / 0.113 | **no / no** |
| smoker arm, all compartments / fetal x host | 20 | 0.001 / 0.321 | 0.003 / 0.246 | yes / **no** |
The shared state is clock-shaped overall and for cord blood x mother, so the fetal-follows-mother readout is gated for cord blood only;
placenta x mother fails (uninformative, not "decoupled"). **Each smoking arm must pass its own gate for a contrast: the non-smoker arm and
the smoker fetal x host pair fail, so a smoking contrast is not estimable** (the `require_host_gate_each_arm` rule would refuse it).
Caveat: a shared non-circadian structure between mother and cord blood (delivery stress, shared processing) is not excluded by a
clock-shaped gate alone; chip/batch layout is not in the metadata.

## Power (`output/power.csv`; 5000 label permutations; MDE = 2.49 x null SD; observed never evaluated)
Structural limit: the mother contributes ONE compartment, and the clocksync host (factor model) needs >= 2 anchor compartments, so
"fetal compartment coupling c / d* on the maternal host" is not available in the public API. Two substitutes:
- Pairwise estimand (A54/A115 `clocksync.paired`, fetal vs maternal blood per delivery), smoker vs non-smoker label permutations:
  cord blood (18 vs 42) coherence null SD 0.144, MDE 0.359 (offset: SD 3.9 h, MDE 9.6 h); placenta (16 vs 35) coherence SD 0.189, MDE 0.470
  (offset SD 4.8 h, MDE 11.8 h). One-group coherence null (fetal states shuffled across deliveries): SD 0.060 (cord), 0.065 (placenta), q95 0.224, 0.244.
  Coherence is bounded by 1, so smoking effects below about 0.36 (cord) are undetectable; d* is not defined for this estimand.
- d* on the free host (anchors = all three compartments, each target on the other two; the host then contains the other fetal tissue, so
  it is NOT the mother's time): smoking contrast, cord blood (19 vs 44) null SD 0.055, MDE 0.137; placenta (17 vs 37) SD 0.099, MDE 0.246.
  Both <= 0.3 (stop rule not triggered) but the arm gates fail, so the contrast stays not estimable.
No labour contrast exists in this series; GSE73685 (13 term women) is the only labour series noted in SCOUT.md.

## What the archive did with placenta (A105-A108c)
Only A106 touched placenta/mother: first-trimester decidua/placenta of the Arutyunyan 2023 single-cell atlas (16 pregnancies, 13 with
fetal VCT/SCT x maternal compartments). H106.1 met (maternal decidua shares a clock state, gate 0.012/0.030); H106.2 uninformative (fetal
VCT+SCT x maternal: p_pairing 0.51, p_acrophase 0.13; couplings on the anchored maternal host VCT 0.34, fibroblast 0.42, SCT 0.19,
Hofbauer 0.09, EVT 0, ambient maternal RNA a confound). Maternal blood was not used there (<= 4 pregnancies). A107/A107b/A108/A108c: fetal
organs (Suo, Cao) and mouse embryo, no mother. A105-A105c: COVID decedent organs. **No term bulk series was used**; GSE27272 does not
overlap with A106 (term, bulk, blood and placenta, 69 deliveries) and has a maternal host of one compartment, unlike A106's several.
(Source: archive `RESULTS_preliminary.md` A106 block; `SCOUT.md` D2.)

## Files
`output/counts.csv`, `gates.csv`, `power.csv`, `clock_gene_coverage.csv`, `gse27272_samples.csv` (metadata, no private data), `clock_table.csv`
(pseudo-counts of 13 clock genes per delivery x compartment).
