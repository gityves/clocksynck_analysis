"""A171 data preparation (blind to the clock statistic): GEO series matrices of the LUAD and LIHC cohorts parsed with the A168
parser (`A168_endocrine/scripts/prep_a168.py`: log2 if linear, probes averaged per symbol, G2M / purity-proxy scores) into
../output/{key}_meta.parquet and {key}_expr.parquet. No coherence is computed here.
Run: python -I prep_a171.py   (reads downloaded text files only; nothing downloaded is executed)"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts"))
import a168_lib as L  # noqa: E402
import prep_a168 as P  # noqa: E402

OUT = HERE.parent / "output"
OUT.mkdir(exist_ok=True)
L.OUT = OUT  # prep_a168.finish writes to L.OUT
P.D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_luad_lihc")
META = HERE.parents[1] / "scouts_2026-10-31/luad_lihc/meta"
SYMBOL_COLUMNS = ("Gene symbol", "Gene Symbol", "GeneSymbol")


def probe_symbols(gpl):
    """As prep_a168.probe_symbols, for both GEO .annot.gz and full .txt platform tables (symbol column named per table)."""
    f = P.D / "platforms" / f"{gpl}.annot.gz"
    opener = (lambda: gzip.open(f, "rt", errors="replace")) if f.exists() else (lambda: open(P.D / "platforms" / f"{gpl}.txt", errors="replace"))
    with opener() as fh:
        lines = [l for l in fh if not l.startswith(("^", "!", "#"))]
    a = pd.read_csv(pd.io.common.StringIO("".join(lines)), sep="\t", dtype=str, low_memory=False)
    col = next(c for c in SYMBOL_COLUMNS if c in a.columns)
    a = a[["ID", col]].dropna().rename(columns={col: "Gene symbol"})
    a["sym"] = a["Gene symbol"].str.split("///")
    a = a.explode("sym"); a["sym"] = a.sym.str.strip().map(P.canon)
    a = a[a.sym != ""]
    P.LOG.append((gpl, "probes_with_symbol", a.ID.nunique()))
    P.LOG.append((gpl, "probes_with_several_symbols", int((a.groupby("ID").size() > 1).sum())))
    return a.set_index("ID")["sym"]


P.probe_symbols = probe_symbols


def series(gse, key, gpl, keep, cols, file=None):
    m = pd.read_csv(META / f"{gse}_samples.csv", dtype=str)
    if file:
        m = m[m.file == file]
    m = m[keep(m)].set_index("gsm")
    m = m[[c for c in cols if c in m]].rename(columns=cols)
    ex = expression_file(gse, file or f"{gse}_series_matrix.txt.gz", gpl)
    P.finish(gse, key, m, ex)


def expression_file(gse, file, gpl):
    """Linear (max > 100) -> log2 with the smallest positive value as floor; probes averaged per symbol (= prep_a168.expression)."""
    path = P.D / gse / file
    ex = pd.read_csv(path, sep="\t", comment="!", index_col=0, quotechar='"', na_values=["null", "NA", ""]).astype(np.float32)
    ex.index = ex.index.astype(str)
    mx, mn = float(np.nanmax(ex.values)), float(np.nanmin(ex.values))
    if mx > 100:
        ex = np.log2(ex.clip(lower=float(ex.values[ex.values > 0].min())))
    P.LOG += [(gse, "deposited_max", mx), (gse, "deposited_min", mn), (gse, "deposited_scale", "linear -> log2" if mx > 100 else "log2")]
    sym = probe_symbols(gpl); sym = sym[sym.index.isin(ex.index)]
    ex = ex.loc[sym.index].groupby(sym.to_numpy()).mean()
    P.LOG += [(gse, "probe_symbol_rows", len(sym)), (gse, "genes", len(ex))]
    return ex


def scores_by(m, X, by):
    """GSE68465: scores z-scored within site (registered)."""
    for _, i in m.groupby(by).indices.items():
        sub = P.scores(X.iloc[i], m.iloc[i].copy())
        m.loc[m.index[i], ["prolif", "purity"]] = sub[["prolif", "purity"]].to_numpy()
    return m


if __name__ == "__main__":
    tum = lambda col, val: (lambda m: m[col].str.strip() == val)
    series("GSE72094", "gse72094", "GPL15048", lambda m: m.gsm.notna(),
           {"patient_id": "patient", "smoking_status": "smoking", "Stage": "stage", "kras_status": "kras", "egfr_status": "egfr",
            "stk11_status": "stk11", "tp53_status": "tp53"})
    series("GSE68465", "gse68465", "GPL96", tum("disease_state", "Lung Adenocarcinoma"),
           {"source": "site", "smoking_history": "smoking", "disease_stage": "stage"})
    # scores within site (registered): recompute from the cached expression
    m, X = L.load("gse68465", OUT); m = scores_by(m, X, "site"); m.to_parquet(OUT / "gse68465_meta.parquet")
    series("GSE14520", "gse14520", "GPL3921", lambda m: m["Tissue Type"].str.strip().isin(["Tumor", "Non-Tumor"]),
           {"Individual": "patient", "Tissue Type": "tissue", "HBV viral status": "hbv", "Cirrhosis": "cirrhosis", "BCLC staging": "bclc"},
           file="GSE14520-GPL3921_series_matrix.txt.gz")
    series("GSE36376", "gse36376", "GPL10558", lambda m: m.gsm.notna(), {"tissue": "tissue", "ajcc stage": "stage"})
    series("GSE31210", "gse31210", "GPL570", tum("tissue", "primary lung tumor"),
           {"smoking status": "smoking", "gene alteration status": "alteration", "pathological stage": "stage"})
    series("GSE41271", "gse41271", "GPL6884", tum("histology", "Adenocarcinoma"), {"tobacco history": "smoking", "final patient stage": "stage"})
    series("GSE63898", "gse63898", "GPL13667", lambda m: m.gsm.notna(), {"tissue subtype": "tissue"})
    pd.DataFrame(P.LOG, columns=["series", "item", "value"]).to_csv(OUT / "prep_log.csv", index=False)
    cnt = []
    for k in ("gse72094", "gse68465", "gse14520", "gse36376", "gse31210", "gse41271", "gse63898"):
        m = pd.read_parquet(OUT / f"{k}_meta.parquet")
        cnt.append(m.assign(series=k).groupby(["series"] + (["tissue"] if "tissue" in m else [])).size().rename("n"))
    print(pd.concat(cnt).to_string())
