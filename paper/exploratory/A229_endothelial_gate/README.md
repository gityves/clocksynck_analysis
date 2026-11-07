# A229: do tumour endothelial adhesion genes and the tumour CD8 T share follow the tumour endothelial clock state? (registered; power stage only)

- `scripts/compute_a229.py --stage power`: eligibility (patients / cohorts per test and side) and the within-cohort permutation null of
  r = |sum y conj(z)| / sqrt(sum y^2 sum |z|^2) from y SHUFFLED across patients only (identity pairing never evaluated, no observed
  statistic computed or stored), pooled with patient-count weights; MDE = 2.49 x null SD; scored if >= 40 patients and MDE <= 0.25.
  Outputs `output/{eligibility,power}.csv` (`reason` in eligibility gives why a cohort drops out of a test).
- `--stage tests`: implemented, NOT run (observed pooled r, permutation p, 200 expression-matched pseudo-clock sets, Holm over primaries
  H229.1-2 and, separately, secondaries H229.3-5; pass = Holm p < 0.05 and r > pseudo q95). Not-scored tests are skipped unless `--force`.
  Descriptive: ICAM1 percentile among 200 expression-matched endothelial genes (drawn with a fixed seed from the 400 closest in mean
  log2(CPM+1); genes present in all cohorts of the test, test and waveform genes excluded), the same tests in adjacent-normal endothelium,
  cancer types with >= 12 patients, implied peak hour.
- z: endothelial complex clock state (`clocksync.paired.projector / states`, GTEx waveform, counts=True), projector fitted on the cohort's
  endothelial rows of that side only; z and y centred within cohort x chemistry over the patients of the test.
- y: ICAM1, VCAM1, SELE = log2(CPM + 1) from the genome-wide endothelial pseudobulk (`compute_a151.genome_pair`; a cohort lacking the gene
  drops out of that test); CD8 share = log((CD8T + .5) / (CD8neg_T + .5)) from A228 `subset_fractions.csv` (samples of the side, summed per
  patient, >= 50 T/NK cells); T/NK share = log(T_NK cells / cells of all compartment pseudobulks of the A227 table).
- Cohorts, tables, pseudo-clock construction: the 34 A227 cohorts (`compute_a227`, `compute_a226`, `a149_lib`). Patients are not
  de-duplicated across cohorts.
- Implied peak: theta = arg(sum y z) = -arg(sum y conj z), hours = (theta mod 2 pi) * 24 / 2 pi (`clocksync.stats.HOURS`, as `Dataset.host`
  `phase_h`). `tests.csv` carries both arg(sum y conj z) in hours and the implied peak.
- Run: `SCCC_DATA=/mnt/max8tb/data/CoPhaser_data/cancer_circadian RUN_CAP_GB=70 paper/common/run.sh
  paper/exploratory/A229_endothelial_gate/scripts/compute_a229.py --stage power|tests --jobs 3`.
