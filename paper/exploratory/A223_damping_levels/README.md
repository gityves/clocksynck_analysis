# A223: tumour-vs-paired-normal clock levels as a damped clock read at daytime surgery (registered)

`scripts/compute_a223.py --stage power|tests` (seconds). Inputs: A186 CPTAC-3 cache (8 types, 713 pairs; library = all
cached ENSG genes), A221 TCGA cache (13 projects), GTEx waveform `clocksync.gtex()`.

**Not scored (power rule).** Rotation-null SD of B 1.13-1.28, MDE 2.8-3.2 against the registered limit 0.5
(`output/power.csv`). With 14 genes, gene-specific tumour-normal changes unrelated to phase are large compared with the
spread of the damping predictor (|x_g| <= about 0.4 log2 at lam = 1), so a phase-organised pattern of the size full
damping would give cannot be separated from rotations. The tests stage is implemented and was not run.
