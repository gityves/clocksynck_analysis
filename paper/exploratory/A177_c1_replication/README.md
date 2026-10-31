# A177: primary replication of A175 H175.1 (preparation, gate and power only)

Status: data preparation, clock gates and the label-only permutation null of D (MDE). The observed D, and any pairing
statistic by proliferation stratum or unstratified, has NOT been computed; the replication is not yet registered.

Cohorts: GSE14520 HCC (GPL3921 only, 213 tumour/non-tumour pairs) and GSE66229 ACRG gastric (GPL570, 98 paired patients of
300 tumours + 100 normals). Proliferation = A97 G2M score of the paired tumours, tertiles by rank within cohort.

- `scripts/prep_a177.py` (`python -I`): paired tumour/normal tables, G2M, tertiles -> `output/{key}_{tumour,normal}.parquet`, `{key}_pairs.csv`, `prep_log.csv`.
- `scripts/power_a177.py`: gates (`gates.csv`), tertile arms (`power_arms.csv`), permutation null of D and MDE (`power.csv`).
- Data: `$SCCC_DATA/gse66229/` (PROVENANCE.md; full series matrix, sha256 in it). GSE14520 from A171's cached table.

Microarray adaptation: `cohort_states(..., counts=False)` (unweighted projector, log2 x ln 2 states, Poisson variances 0), so
the A175 statistic is the plain |complex correlation|. GSE66229 values (0.57-4.47) are on the GSE62254 scale and are converted
x log2(10) as in A173 (scale-free for the statistic and the G2M ranks). CIART is absent on GPL3921 (13 clock genes).
The pair gate is the `Dataset.gate` logic written for two arrays; the reference gate is `paired.reference_gate` on the normals.
Run `python -I prep_a177.py`, then `python power_a177.py` (6 s).
