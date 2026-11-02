"""Outcome and covariate tables for A193 and A194 (SCAN-B, TCGA BRCA, METABRIC), keyed as the readout scripts key tumours.

Each loader returns a DataFrame indexed by the tumour key with `time` (months), `event` (0/1, overall survival), `g2m`,
`MKI67`, the cohort's covariates (numeric, complete-case filtering left to the caller) and `pam50`.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EXP = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(EXP / "A127_scanb" / "scripts"), str(EXP / "A166_metabric" / "scripts"),
                str(EXP / "C8_internal_time_outcomes" / "scripts")]
import a127_lib as A  # noqa: E402
import a166_lib as L  # noqa: E402
import score_a191 as S  # noqa: E402

TCGA_CLIN = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/tcga/pancan_subtypes/brca_tcga_pan_can_atlas_2018")
COVARIATES = {"SCAN-B": ["age", "log_size", "node_pos", "grade3", "grade1", "er", "her2"],
              "TCGA-BRCA": ["age", "t2", "t34", "node_pos", "er", "erbb2"],
              "METABRIC": ["age", "log_size", "nodes", "grade3", "grade1", "er", "erbb2"]}


def _z(x):
    return (x - x.mean()) / x.std()


def scanb():
    m, X = A.load()
    c = pd.DataFrame(index=m.index)
    c["time"] = pd.to_numeric(m["overall survival days"], errors="coerce") / 30.4375
    c["event"] = pd.to_numeric(m["overall survival event"], errors="coerce")
    c["age"] = pd.to_numeric(m["age at diagnosis"], errors="coerce")
    c["log_size"] = np.log1p(pd.to_numeric(m["tumor size"], errors="coerce"))
    c["node_pos"] = m["lymph node status"].map({"NodePositive": 1.0, "NodeNegative": 0.0})
    g = m["nhg"].astype(str)
    c["grade3"], c["grade1"] = (g == "G3").astype(float).where(g.isin(["G1", "G2", "G3"])), (g == "G1").astype(float).where(g.isin(["G1", "G2", "G3"]))
    c["grade12"] = g.isin(["G1", "G2"]).astype(float).where(g.isin(["G1", "G2", "G3"]))
    c["er"] = m.er.astype(float)
    c["her2"] = pd.to_numeric(m["her2 status"], errors="coerce")
    c["g2m"], c["MKI67"], c["pam50"] = m.prolif, _z(X["MKI67"]), m.pam50
    return c, X


def tcga():
    m = pd.read_parquet(L.OUT / "tcga_meta.parquet")
    X = pd.read_parquet(L.OUT / "tcga_expr.parquet").loc[m.index]
    p = pd.read_csv(TCGA_CLIN / "data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
    p = p.reindex(m.patient)
    c = pd.DataFrame(index=m.index)
    c["time"] = p.OS_MONTHS.to_numpy()
    c["event"] = p.OS_STATUS.astype(str).str.startswith("1").astype(float).where(p.OS_STATUS.notna()).to_numpy()
    c["age"] = p.AGE.to_numpy(float)
    t = p.PATH_T_STAGE.astype(str).str.upper().str[:2].to_numpy()
    ok = np.isin(t, ["T1", "T2", "T3", "T4"])
    c["t2"] = np.where(ok, (t == "T2").astype(float), np.nan)
    c["t34"] = np.where(ok, np.isin(t, ["T3", "T4"]).astype(float), np.nan)
    n = p.PATH_N_STAGE.astype(str).str.upper().to_numpy()
    c["node_pos"] = np.where(np.char.startswith(n.astype(str), "N0"), 0.0, np.where(np.char.startswith(n.astype(str), "N"), 1.0, np.nan))
    c.loc[np.char.startswith(n.astype(str), "NX"), "node_pos"] = np.nan
    c["m1"] = (p.PATH_M_STAGE.astype(str).str.upper().str.startswith("M1")).astype(float).to_numpy()
    c["er"] = m.er.astype(float)
    c["erbb2"] = _z(X["ERBB2"])
    c["g2m"], c["MKI67"], c["pam50"], c["patient"] = m.prolif, _z(X["MKI67"]), m.pam50, m.patient
    return c, X


def metabric():
    m, X = S.metabric()
    c = pd.DataFrame(index=m.index)
    c["time"], c["event"] = m.os_t, m.os_e
    c["rfs_t"], c["rfs_e"] = m.rfs_t, m.rfs_e
    c["age"], c["log_size"], c["nodes"] = m.AGE_AT_DIAGNOSIS, np.log1p(m.TUMOR_SIZE), m.nodes
    g = m.GRADE
    c["grade3"], c["grade1"] = (g == 3).astype(float).where(g.notna()), (g == 1).astype(float).where(g.notna())
    c["grade12"] = g.isin([1, 2]).astype(float).where(g.notna())
    c["er"], c["erbb2"] = m.er.astype(float), _z(X["ERBB2"])
    c["g2m"], c["MKI67"], c["pam50"], c["mb_cohort"] = m.prolif, _z(X["MKI67"]), m.pam50, m.cohort
    return c, X


LOADERS = {"SCAN-B": scanb, "TCGA-BRCA": tcga, "METABRIC": metabric}


def block_design(d, cov_by_cohort, extra=()):
    """Covariates as cohort-specific columns (zero outside their cohort), so each cohort keeps its own covariate effects in
    a model stratified by cohort; `extra` columns are shared."""
    cols = []
    for k, cov in cov_by_cohort.items():
        inside = (d.cohort == k).to_numpy()
        for v in cov:
            cols.append(np.where(inside, d[v].to_numpy(float), 0.0))
    for v in extra:
        cols.append(d[v].to_numpy(float))
    return np.column_stack(cols) if cols else None
