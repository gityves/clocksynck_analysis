# C8b: replication and proliferation test of published clock-survival claims (A192-A194)

Registered in `paper/docs/registered_experiments.md` (commit a6d5922), on the user's instruction (2026-11-02) to make the
A191 null convincing to the groups who reported clock-survival associations. Literature: `paper/docs/literature_positioning.md`
section 10.

| amendment | claim tested | cohorts | script |
|---|---|---|---|
| A192 | clock genes predict metastasis-free survival (Cadenas et al., Cell Cycle 2014), and whether proliferation carries it | Mainz GSE11121, Rotterdam GSE2034, Transbig GSE7390 (node-negative, untreated; fresh) | `scripts/score_a192.py` |
| A193 | high CYCLOPS magnitude predicts worse survival in luminal A (Li et al., PNAS 2024) | TCGA (reproduction), METABRIC and SCAN-B luminal A | `scripts/cyclops_magnitude.py` (readout), `scripts/score_a193.py` |
| A194 | TimeTeller clock dysfunction Theta is an independent prognostic factor (Vlachou et al., bioRxiv 2019) | SCAN-B, TCGA BRCA; METABRIC secondary | `scripts/timeteller_theta.py` (readout), `scripts/score_a194.py` |

Data: `paper/common/datasets/breast_prognosis_geo/build.py` (A192 cohorts; provenance in
`$SCCC_DATA/breast_prognosis_geo/PROVENANCE.md`). Li et al.'s code (read, not run): `$SCCC_DATA/li2024_repo/`.

Outputs (`output/`, not versioned): `a192_power.csv`, `a192_reproduction.csv`, `a192_tests.csv`, `a192_prolif_match.csv`,
`a192_secondary.csv`; A193 and A194 outputs listed in their scripts. Only summary statistics are written for outcomes; GTEx
enters only as the TimeTeller training set, and no per-donor GTEx value is written.

    python scripts/score_a192.py [--check] [--power-only]      # about 5 min
