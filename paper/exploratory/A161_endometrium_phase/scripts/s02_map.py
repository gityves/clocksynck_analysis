"""A161 step 2 (before any clock statistic): gene-level matrices and the 14-gene mapping (probe with the highest mean over all
samples of the series). Writes output/clock_gene_map_<series>.csv and caches gene x sample log2 matrices in $SCCC_DATA/.../derived.
python -I s02_map.py"""
import gzip, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
import a161_lib as L
from clocksync.genes import canonical
import clocksync as cs

GENES = list(cs.gtex().genes)
DER = L.EN / "derived"; DER.mkdir(exist_ok=True)


def annot(path, sym="Gene symbol"):
    with gzip.open(path, "rt", errors="replace") as f:
        lines = [ln for ln in f if not ln.startswith(("#", "!", "^"))]
    import io
    a = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str, usecols=["ID", sym]).fillna("")
    return a.set_index("ID")[sym]


def collapse(E: pd.DataFrame, probe_syms: pd.Series, series: str):
    """E probes x samples (log2). Clock genes: any probe whose symbol list contains the gene; highest mean. Other genes:
    single-symbol probes only, highest mean per symbol. Returns gene x sample frame (clock genes under their canonical name)."""
    syms = probe_syms.reindex(E.index).fillna("")
    lst = syms.map(lambda s: [canonical(x) for x in s.split("///") if x.strip()])
    mean = E.mean(axis=1)
    rows, pick = [], {}
    for g in GENES:
        c = [p for p, l in lst.items() if g in l]
        if not c:
            rows.append(dict(gene=g, probe="", mean=np.nan, n_probes=0, multi_symbol=False, symbols="")); continue
        p = mean[c].idxmax()
        pick[g] = p
        rows.append(dict(gene=g, probe=p, mean=mean[p], n_probes=len(c), multi_symbol=len(lst[p]) > 1, symbols="///".join(lst[p])))
    mp = pd.DataFrame(rows); mp.to_csv(L.OUT / f"clock_gene_map_{series}.csv", index=False)
    single = pd.Series({p: l[0] for p, l in lst.items() if len(l) == 1})
    best = mean[single.index].groupby(single).idxmax()
    G = E.loc[best.values].copy(); G.index = best.index
    for g, p in pick.items():
        G.loc[g] = E.loc[p].values
    return G, mp


def save(G, series):
    G.to_parquet(DER / f"{series}_gene_log2.parquet")


# GSE51981
E = L.matrix(L.EN / "GSE51981_series_matrix.txt.gz")
print("GSE51981 value range", np.nanmin(E.values), np.nanmax(E.values), E.shape, "NaN", int(E.isna().sum().sum()))
G, mp = collapse(E, annot(L.EN / "GPL570.annot.gz"), "GSE51981"); save(G, "GSE51981"); print(mp, G.shape)
# GSE6364
E = L.matrix(L.EN / "GSE6364_series_matrix.txt.gz")
print("GSE6364 value range (raw MAS5-type intensities, log2 taken with a floor of 1)", np.nanmin(E.values), np.nanmax(E.values), E.shape)
E = np.log2(E.clip(lower=1.0))
G, mp = collapse(E, annot(L.EN / "GPL570.annot.gz"), "GSE6364"); save(G, "GSE6364"); print(mp, G.shape)
# GSE141549 (supplementary normalised values; platform annotation = the ILMN_Gene column of the file)
for tag in ("HT12", "WG6"):
    d = pd.read_parquet(DER / f"GSE141549_{tag}.parquet")
    E = d.drop(columns="ILMN_Gene")
    print(tag, "range", np.nanmin(E.values), np.nanmax(E.values), E.shape)
    G, mp = collapse(E, d["ILMN_Gene"].astype(str), f"GSE141549_{tag}"); save(G, f"GSE141549_{tag}"); print(mp, G.shape)
