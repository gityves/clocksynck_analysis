"""A168 data preparation (blind to the clock statistic): parse the three GEO series matrices (GSE105777, GSE20181, GSE87411) with their
platform annotations into tidy tables cached under ../output/: {series}_meta.parquet (sample table: patient, arm, time point, label,
proliferation / purity / oestrogen-response scores) and {series}_expr.parquet (samples x gene symbols, log2, float32), plus
sample_counts.csv and prep_log.csv (log2 decisions, probe counts, imputed cells). No coherence is computed here.
Run: python -I prep_a168.py   (reads downloaded text files only; nothing downloaded is executed)"""
import gzip
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a168_lib as L  # noqa: E402

D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_endocrine")
PLATFORM = {"GSE105777": "GPL10558", "GSE20181": "GPL96", "GSE87411": "GPL6480"}
L.OUT.mkdir(exist_ok=True)
LOG = []


def header(path):
    """Series-matrix header lines up to the table: key -> list of value lists."""
    h = {}
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!"):
                k, *v = line.rstrip("\n").split("\t"); h.setdefault(k, []).append([x.strip('"') for x in v])
    return h


def probe_symbols(gpl):
    """Long table probe -> gene symbol of the platform annotation. A probe annotated to several symbols ('A /// B': overlapping genes,
    e.g. THRA /// NR1D1) counts towards each of them; probes without a symbol are dropped."""
    with gzip.open(D / "platforms" / f"{gpl}.annot.gz", "rt", errors="replace") as f:
        lines = [l for l in f if not l.startswith(("^", "!", "#"))]
    a = pd.read_csv(pd.io.common.StringIO("".join(lines)), sep="\t", usecols=["ID", "Gene symbol"], dtype=str).dropna()
    a["sym"] = a["Gene symbol"].str.split("///")
    a = a.explode("sym"); a["sym"] = a.sym.str.strip().map(canon)
    LOG.append((gpl, "probes_with_symbol", a.ID.nunique())); LOG.append((gpl, "probes_with_several_symbols", int((a.groupby("ID").size() > 1).sum())))
    return a.set_index("ID")["sym"]


def expression(gse):
    """(genes x samples) log2 expression: probes mapped to symbols, log2-transformed if linear (max > 100), averaged per symbol."""
    path = D / gse / f"{gse}_series_matrix.txt.gz"
    ex = pd.read_csv(path, sep="\t", comment="!", index_col=0, quotechar='"', na_values=["null", "NA", ""]).astype(np.float32)
    ex.index = ex.index.astype(str)
    mx, mn = float(np.nanmax(ex.values)), float(np.nanmin(ex.values))
    linear = mx > 100
    if linear:
        pos = ex.values[ex.values > 0]
        ex = ex.clip(lower=float(pos.min()))  # non-positive values (if any) set to the smallest positive one
        ex = np.log2(ex)
    LOG.append((gse, "deposited_max", mx)); LOG.append((gse, "deposited_min", mn)); LOG.append((gse, "deposited_scale", "linear -> log2" if linear else "log2 (log2 ratio for GSE87411)"))
    sym = probe_symbols(PLATFORM[gse]); sym = sym[sym.index.isin(ex.index)]
    ex = ex.loc[sym.index]  # one row per (probe, symbol)
    LOG.append((gse, "probe_symbol_rows", len(ex)))
    ex = ex.groupby(sym.to_numpy()).mean()  # A166 convention: probes averaged per symbol (log scale)
    LOG.append((gse, "genes", len(ex)))
    return ex


_canon = {}
def canon(s):  # canonical clock-gene symbol (aliases such as C1orf51 -> CIART, BMAL1 -> ARNTL); others unchanged
    if s not in _canon:
        c = L.canonical(s); _canon[s] = c if c in L.GENES else s
    return _canon[s]


def scores(X, m):
    """Proliferation (G2M), purity proxy (A97 11 genes) and oestrogen-response (early + late) scores: mean gene z over the series, z-scored."""
    def z(genes):
        s = L.sh._zmean(X, [g for g in genes if g in X.columns]); return ((s - s.mean()) / s.std()).to_numpy()
    m["prolif"] = z(L.sh.hallmark_g2m()); m["purity"] = z(L.sh.PURITY); m["estro"] = z(L.estrogen_genes())
    return m


def finish(gse, key, m, ex):
    ex = ex[m.index]  # samples in the table's order
    nan_clock = int(ex.loc[[g for g in L.GENES if g in ex.index]].isna().sum().sum())
    ex = ex[ex.isna().mean(axis=1) <= 0.2]  # genes missing in > 20 % of samples dropped
    nan_cells = int(ex.isna().sum().sum())
    X = ex.T.astype(np.float32); X = X.fillna(X.mean())  # remaining NaN cells: gene mean
    m = scores(X, m)
    m["pos"] = np.arange(len(m))
    X.to_parquet(L.OUT / f"{key}_expr.parquet"); m.to_parquet(L.OUT / f"{key}_meta.parquet")
    LOG.append((gse, "clock_genes_present", ",".join(g for g in L.GENES if g in X.columns)))
    LOG.append((gse, "nan_cells_clock_genes_before_filter", nan_clock)); LOG.append((gse, "nan_cells_imputed", nan_cells))
    LOG.append((gse, "n_samples", len(m))); LOG.append((gse, "n_genes_kept", X.shape[1]))
    print(key, len(m), X.shape, flush=True)


def poetic():
    gse = "GSE105777"; h = header(D / gse / f"{gse}_series_matrix.txt.gz")
    t = h["!Sample_title"][0]; gsm = h["!Sample_geo_accession"][0]
    rows = []
    for g, x in zip(gsm, t):
        a = re.fullmatch(r"AI\.Treated\.(\d+)([BS])", x); u = re.fullmatch(r"(\d+)([DS]) \[reanalysis\]", x)
        assert a or u, x
        rows.append((g, f"AI{a.group(1)}" if a else f"UNT{u.group(1)}", "AI" if a else "UNT", (a or u).group(2), x))
    m = pd.DataFrame(rows, columns=["gsm", "patient", "arm", "tp", "title"]).set_index("gsm")
    # consistency with the characteristics lines (group / timepoint / sampling time)
    flat = [r for r in h["!Sample_characteristics_ch1"]]
    grp = np.array([next((r[i] for r in flat if r[i].startswith("group:")), "") for i in range(len(t))])
    assert ((grp == "group: treated") == (m.arm.to_numpy() == "AI")).all() and ((grp == "group: untreated") == (m.arm.to_numpy() == "UNT")).all()
    finish(gse, "poetic", m, expression(gse))


def letrozole():
    gse = "GSE20181"; h = header(D / gse / f"{gse}_series_matrix.txt.gz")
    t = h["!Sample_title"][0]; gsm = h["!Sample_geo_accession"][0]
    rows = []
    for g, x in zip(gsm, t):
        r = re.match(r"(\d+)([ABC])[;:]", x); assert r, x
        resp = x.strip().split(";")[-1].strip()
        rows.append((g, r.group(1), "LET", r.group(2), resp, x))
    m = pd.DataFrame(rows, columns=["gsm", "patient", "arm", "tp", "response", "title"]).set_index("gsm")
    assert (m.groupby("patient").response.nunique() == 1).all()  # response label constant within patient
    finish(gse, "letrozole", m, expression(gse))


def z1031b():
    gse = "GSE87411"; h = header(D / gse / f"{gse}_series_matrix.txt.gz")
    t = h["!Sample_title"][0]; gsm = h["!Sample_geo_accession"][0]
    cls = {"BL": "base", "BL Redo": "base", "BL_Redo": "base", "re-analyzed": "reanalyzed", "w2": "w2", "W2": "w2", "w2 Redo": "w2",
           "w2_Redo": "w2", "W2_Redo": "w2", "M": "month", "M_Redo": "month"}
    rows = []
    for g, x in zip(gsm, t):
        p, tp = re.search(r"Patient ID_(\w+?)_(.*)$", x).groups()
        rows.append((g, p, "Z1031B", cls[tp], tp, x))
    m = pd.DataFrame(rows, columns=["gsm", "patient", "arm", "tp", "label", "title"]).set_index("gsm")
    finish(gse, "z1031b", m, expression(gse))


if __name__ == "__main__":
    poetic(); letrozole(); z1031b()
    pd.DataFrame(LOG, columns=["series", "item", "value"]).to_csv(L.OUT / "prep_log.csv", index=False)
    cnt = pd.concat([pd.read_parquet(L.OUT / f"{k}_meta.parquet").assign(series=k) for k in L.SERIES]).groupby(["series", "arm", "tp"]).size()
    cnt.rename("n_samples").to_csv(L.OUT / "sample_counts.csv"); print(cnt.to_string())
