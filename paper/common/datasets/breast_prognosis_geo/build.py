"""Node-negative, systemically untreated breast cancer cohorts on Affymetrix U133A (GPL96), as in Cadenas et al. 2014:
Mainz (GSE11121, 200), Rotterdam (GSE2034, 286) and Transbig (GSE7390, 198).

Writes to `$SCCC_DATA/breast_prognosis_geo/`:
- `expr.parquet`: tumours x genes, log2 intensity (MAS5 values of GSE11121 and GSE2034 log2-transformed after flooring at 1;
  GSE7390 is already log2), one probe per gene (highest mean over the three cohorts, GPL96 annotation, single-symbol probes);
- `clinical.csv`: cohort, time to distant metastasis (months) and event, ER (recorded: GSE2034, GSE7390), grade and size
  (GSE11121, GSE7390), age (GSE7390). GSE2034's events are distant metastases (series summary: "106 ... developed a distant
  metastasis"; its clinical table is in the GEO text record, `GSE2034_clinical.txt`).

Untrusted downloads: parse only (run with `python -I`).

    python -I build.py
"""
import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd

D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
OUT = D / "breast_prognosis_geo"
GPL96 = D / "scout_endocrine" / "platforms" / "GPL96.annot.gz"
COHORTS = {"GSE11121": "Mainz", "GSE2034": "Rotterdam", "GSE7390": "Transbig"}


def series(acc):
    """(probe x sample matrix, sample characteristics as a DataFrame)."""
    with gzip.open(OUT / f"{acc}_series_matrix.txt.gz", "rt") as f:
        lines = f.read().splitlines()
    i0 = lines.index("!series_matrix_table_begin")
    i1 = lines.index("!series_matrix_table_end")
    X = pd.read_csv(io.StringIO("\n".join(lines[i0 + 1:i1])), sep="\t", index_col=0)
    gsm = next(l for l in lines if l.startswith("!Sample_geo_accession")).split("\t")[1:]
    gsm = [g.strip('"') for g in gsm]
    ch = {}
    for l in lines:
        if l.startswith("!Sample_characteristics_ch1"):
            vals = [v.strip('"') for v in l.split("\t")[1:]]
            key = vals[0].split(":")[0].strip()
            ch[key] = [v.split(":", 1)[1].strip() if ":" in v else np.nan for v in vals]
    return X, pd.DataFrame(ch, index=gsm)


def num(s):
    return pd.to_numeric(s, errors="coerce")


def main():
    mats, clin = [], []
    for acc, name in COHORTS.items():
        X, ch = series(acc)
        X = X.astype(float)
        if acc != "GSE7390":
            X = np.log2(X.clip(lower=1))
        c = pd.DataFrame(index=X.columns)
        c["cohort"] = name
        if acc == "GSE11121":
            c["time"], c["event"] = num(ch["t.dmfs"]), num(ch["e.dmfs"])
            c["grade"], c["size"] = num(ch["grade"]), num(ch["size_in_cm"])
        elif acc == "GSE7390":
            c["time"], c["event"] = num(ch["t.dmfs"]) / 30.4375, num(ch["e.dmfs"])
            c["grade"], c["size"], c["age"], c["er"] = num(ch["grade"]), num(ch["size"]), num(ch["age"]), num(ch["er"])
        else:
            t = pd.read_csv(OUT / "GSE2034_clinical.txt", sep="\t").set_index("GEO asscession number")
            t = t.reindex(c.index)
            c["time"], c["event"] = num(t["time to relapse or last follow-up (months)"]), num(t["relapse (1=True)"])
            c["er"] = t["ER Status"].map({"ER+": 1.0, "ER-": 0.0})
        mats.append(X)
        clin.append(c)
    X = pd.concat(mats, axis=1)
    a = pd.read_csv(GPL96, sep="\t", skiprows=27, usecols=["ID", "Gene symbol"]).dropna()
    a = a.assign(sym=a["Gene symbol"].str.split("///")).explode("sym")
    single = set(a.loc[~a["Gene symbol"].str.contains("///"), "sym"])
    # multi-symbol probes count for a symbol only when it has no single-symbol probe (PER1: MIR6883///PER1, NR1D1///THRA)
    a = a[~a["Gene symbol"].str.contains("///") | ~a.sym.isin(single)]
    a = a[a.ID.isin(X.index)]
    best = a.assign(m=X.mean(axis=1).reindex(a.ID).to_numpy()).sort_values("m").groupby("sym").ID.last()
    E = X.loc[best.values].set_axis(best.index).T.astype("float32")
    c = pd.concat(clin)
    assert E.index.equals(c.index)
    E.to_parquet(OUT / "expr.parquet")
    c.rename_axis("gsm").to_csv(OUT / "clinical.csv")
    print(E.shape, c.groupby("cohort").agg(n=("event", "size"), events=("event", "sum")).to_dict())


if __name__ == "__main__":
    main()
