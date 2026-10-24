"""Breast-cancer mastectomies with recorded resection time: tumour and adjacent normal of 43 patients (GSE233242).

Source: `$SCCC_DATA/validation/GSE233242/GSE233242_Processed_file_tpm.csv.gz` (TPM and transcript length per gene, GEO
TPM only; Li et al. 2024, PNAS 121:e2311854121), joined to the recorded resection times (SI Appendix Table S1 and the GEO
sample key) in `resection_times.csv` next to this file (copy of the legacy tracked table
`desync/validation/GSE233242_resection_times.csv`: patient id, time of resection, GSM and TPM column of each library).
Legacy builder replaced: `desync/validation/scripts/v25_5_breast_timed.py:samples` and `pseudobulk` (V25.5, A25 post
hoc; also E4 of extended_clock/scripts/evaluate.py, A82 hospital-hours check).
Used by: V25.5 (internal time vs recorded time, hospital prior), E4.

Choices (as legacy): approximate counts = TPM x transcript length, rescaled per library to `n_reads` reads
(primary 30 M; 10 M and 60 M are sensitivities that only scale the Poisson covariance), first row per gene symbol,
clock-gene counts rounded; depth = n_reads. Compartments `normal` and `tumour` (both form the host), donor =
patient id (Table S1). One patient is spelled W006063 in GEO and W006053 in Table S1: excluded in the primary run
(`join_typo=False`), joined in the sensitivity (`join_typo=True`). Donors table: `resection_h` (h), `resection_time`,
`age` (y; SI Table S1, added by `add_age.py`, A130), `subtype` (GEO). The time of resection is a recorded hospital time (08-18 h window), not a circadian time.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import clocksync as cs
from common.paths import DATA

HERE = Path(__file__).resolve().parent
TPM = DATA / "validation" / "GSE233242" / "GSE233242_Processed_file_tpm.csv.gz"
TRUTH = HERE / "resection_times.csv"


def build(n_reads: float = 30e6, join_typo: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    t = pd.read_csv(TRUTH)
    if not join_typo:
        t = t[t.note.fillna("") == ""]
    key = pd.concat([pd.DataFrame(dict(patient=t.patient_id, compartment=c, SL=t[f"tpm_column_{c}"]))
                     for c in ("normal", "tumour")], ignore_index=True)
    x = pd.read_csv(TPM).drop_duplicates("gene_name").set_index("gene_name")
    rl = x[key.SL].mul(x.transcriptLen, axis=0)
    genes = list(cs.gtex().genes)
    cnt = (rl / rl.sum() * n_reads).loc[genes].T.round()
    cnt.index = range(len(key))
    pb = cs.from_bulk(cnt, key, donor="patient", compartment="compartment",
                      depth=pd.Series(float(n_reads), index=cnt.index))
    donors = t.set_index("patient_id")[["resection_h", "age", "resection_time", "subtype_geo"]].rename(
        columns={"subtype_geo": "subtype"})
    return pb, donors
