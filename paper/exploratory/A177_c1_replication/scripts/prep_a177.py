"""A177 data preparation (blind to every pairing statistic): paired tumour / normal tables of the two replication cohorts of
A175 H175.1. No pairing statistic, stratified or not, is computed here; only expression, the A97 G2M score of the paired
tumours and its tertile labels (rank within cohort).

  gse14520  HCC, GPL3921 only: tumour (A) and non-tumour (B) liver of one patient; cached A171 table (log2, probes averaged
            per symbol), pairs = patients with exactly one of each (A171 rule)
  gse66229  ACRG gastric, GPL570: 'tissue: Gastric tumor' vs 'normal gastric tissue' of one patient. Parsed with the A173
            conventions (A171 `expression_file`: probes averaged per symbol). Deposited values 0.57-4.47 are the same
            scale as GSE62254 (same ACRG tumour set), which A173 registered as log10 -> x log2(10).

Outputs (output/): {key}_tumour.parquet, {key}_normal.parquet (pairs x genes, rows aligned, log2), {key}_pairs.csv
(patient, prolif, tertile 0 low / 1 mid / 2 high), prep_log.csv.
Run: python -I prep_a177.py   (reads downloaded text files only; nothing downloaded is executed)"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "paper/exploratory/A171_luad_lihc/scripts"))
sys.path.insert(0, str(ROOT / "paper/exploratory/A168_endocrine/scripts"))
import a168_lib as L  # noqa: E402
import prep_a171 as Q  # noqa: E402

OUT = HERE.parent / "output"
OUT.mkdir(exist_ok=True)
SCCC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
PLATFORMS = SCCC / "scout_gastric_lung"  # GPL570.annot.gz
_probe_symbols = Q.probe_symbols


def probe_symbols(gpl):
    """A171 platform parser reading the platform table from the A173 folder; series files are read from SCCC / <series>."""
    Q.P.D = PLATFORMS
    try:
        return _probe_symbols(gpl)
    finally:
        Q.P.D = SCCC


Q.probe_symbols = probe_symbols
Q.P.D = SCCC
A171 = ROOT / "paper/exploratory/A171_luad_lihc/output"
LOG = []


def tertiles(x) -> np.ndarray:
    return pd.qcut(pd.Series(x).rank(method="first"), 3, labels=False).to_numpy()


def finish(key, patients, tum, nor):
    """Rows aligned by patient; G2M = A97 (mean z of HALLMARK_G2M_CHECKPOINT without clock genes) over the paired tumours."""
    g2m = L.sh._zmean(tum, L.sh.hallmark_g2m()).to_numpy()
    tum.to_parquet(OUT / f"{key}_tumour.parquet")
    nor.to_parquet(OUT / f"{key}_normal.parquet")
    pd.DataFrame({"patient": patients, "prolif": g2m, "tertile": tertiles(g2m)}).to_csv(OUT / f"{key}_pairs.csv", index=False)
    LOG.extend([(key, "n_pairs", len(patients)), (key, "n_genes", tum.shape[1]),
                (key, "clock_genes_present", ",".join(g for g in L.GENES if g in tum.columns)),
                (key, "g2m_genes_present", len([g for g in L.sh.hallmark_g2m() if g in tum.columns]))])
    print(key, len(patients), tum.shape, flush=True)


def gse14520():
    m = pd.read_parquet(A171 / "gse14520_meta.parquet")
    X = pd.read_parquet(A171 / "gse14520_expr.parquet")
    m["patient"] = m.patient.str.strip().str[:-1]  # Individual = patient id + A (tumour) / B (non-tumour)
    tis = m.tissue.str.strip()
    a = m[tis == "Non-Tumor"].drop_duplicates("patient", keep=False)
    b = m[tis == "Tumor"].drop_duplicates("patient", keep=False)
    both = sorted(set(a.patient) & set(b.patient))
    ga = a.reset_index().set_index("patient").loc[both, "gsm"]
    gb = b.reset_index().set_index("patient").loc[both, "gsm"]
    finish("gse14520", both, X.loc[gb.to_numpy()].reset_index(drop=True), X.loc[ga.to_numpy()].reset_index(drop=True))


def gse66229():
    f = SCCC / "gse66229" / "GSE66229_series_matrix.txt.gz"
    rows = {}
    with gzip.open(f, "rt") as fh:
        for l in fh:
            if l.startswith("!Sample_geo_accession"):
                gsm = [x.strip('"') for x in l.rstrip("\n").split("\t")[1:]]
            elif l.startswith("!Sample_characteristics_ch1"):
                v = [x.strip('"') for x in l.rstrip("\n").split("\t")[1:]]
                rows[v[0].split(":")[0]] = [x.split(":", 1)[1].strip() for x in v]
            elif l.startswith('"ID_REF"'):
                break
    m = pd.DataFrame(rows, index=gsm)
    assert set(m.tissue) == {"Gastric tumor", "normal gastric tissue"}, set(m.tissue)
    ex = Q.expression_file("gse66229", "GSE66229_series_matrix.txt.gz", "GPL570")
    mx = [r for r in Q.P.LOG if r[0] == "gse66229" and r[1] == "deposited_max"][0][2]
    assert mx < 10, mx  # same range as GSE62254 (log10, A173); a log2 scale would reach ~14
    ex = ex * np.log2(10)
    Q.P.LOG.append(("gse66229", "deposited_scale", "log10 -> log2 (A173 convention, same ACRG set as GSE62254)"))
    t = m[m.tissue == "Gastric tumor"].drop_duplicates("patient", keep=False).reset_index(names="gsm").set_index("patient")
    n = m[m.tissue == "normal gastric tissue"].drop_duplicates("patient", keep=False).reset_index(names="gsm").set_index("patient")
    both = sorted(set(t.index) & set(n.index), key=int)
    ex = ex[ex.isna().mean(axis=1) <= 0.2]
    X = ex.T.astype(np.float32)
    X = X.fillna(X.mean())
    finish("gse66229", both, X.loc[t.loc[both, "gsm"].to_numpy()].reset_index(drop=True), X.loc[n.loc[both, "gsm"].to_numpy()].reset_index(drop=True))
    LOG.append(("gse66229", "n_tumours_total", int((m.tissue == "Gastric tumor").sum())))
    LOG.append(("gse66229", "n_normals_total", int((m.tissue != "Gastric tumor").sum())))


if __name__ == "__main__":
    gse14520()
    gse66229()
    pd.DataFrame(LOG + Q.P.LOG, columns=["series", "item", "value"]).to_csv(OUT / "prep_log.csv", index=False)
