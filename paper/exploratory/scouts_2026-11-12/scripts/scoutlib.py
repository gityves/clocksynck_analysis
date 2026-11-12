"""Helpers for the lesion / multi-tissue scout: parse GEO series matrices (headers, expression tables), probe->symbol maps,
clock-gene presence and median percentile among all genes. No clock statistic beyond that is computed."""
import gzip, io, re, csv
import numpy as np, pandas as pd
CLOCK = ["ARNTL","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","DBP","RORA","TEF","CIART"]
ALIAS = {"BMAL1":"ARNTL","CHRONO":"CIART","C1orf51":"CIART","REV-ERBA":"NR1D1","REV-ERBB":"NR1D2"}
ROOT = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian"

def _open(p): return io.TextIOWrapper(gzip.open(p, "rb"), errors="replace")

def read_sheet(p):
    """Return DataFrame (one row per sample) with title, geo, and every 'key: value' characteristic as a column."""
    meta, chars = {}, []
    with _open(p) as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"): break
            if not l.startswith("!Sample_"): continue
            k, _, v = l.rstrip("\n").partition("\t"); vals = [x.strip('"') for x in v.split("\t")]
            if k == "!Sample_characteristics_ch1" or k == "!Sample_characteristics_ch2": chars.append(vals)
            else: meta.setdefault(k[8:], vals)
    d = pd.DataFrame({"geo": meta["geo_accession"], "title": meta["title"]})
    for k in ("platform_id", "source_name_ch1", "description"):
        if k in meta: d[k] = meta[k]
    for row in chars:
        keys = {r.split(":", 1)[0].strip() for r in row if ":" in r}
        if len(keys) == 1:
            k = keys.pop(); d[k] = [r.split(":", 1)[1].strip() if ":" in r else "" for r in row]
        else:  # mixed keys on one row (GEO quirk): split by key
            for r_i, r in enumerate(row):
                if ":" in r:
                    k, v = r.split(":", 1); d.loc[r_i, k.strip()] = v.strip()
    return d

def read_matrix(p):
    with _open(p) as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"): break
        df = pd.read_csv(f, sep="\t", index_col=0, quotechar='"', comment=None, low_memory=False)
    df = df[~df.index.astype(str).str.startswith("!")]
    return df

def read_annot(p, symcols=("Gene symbol", "Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_symbol")):
    """probe -> symbol from GEO .annot.gz (table between !platform_table_begin/end)."""
    with _open(p) as f:
        for l in f:
            if l.startswith("!platform_table_begin"): break
        a = pd.read_csv(f, sep="\t", index_col=0, low_memory=False, quoting=csv.QUOTE_NONE)
    a = a[~a.index.astype(str).str.startswith("!")]
    col = next((c for c in symcols if c in a.columns), None)
    if col is None: raise KeyError(f"no symbol column in {p}: {list(a.columns)[:15]}")
    return a[col].astype(str).str.split("///").str[0].str.strip()

def clock_report(expr, sym, higher_is_more=True):
    """expr: genes/probes x samples (numeric); sym: Series probe->symbol (or None if index already symbols).
    Returns DataFrame gene, present, n_probes, median_percentile (rank among all genes, per sample, median over samples)."""
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if sym is not None:
        s = sym.reindex(expr.index); expr = expr[s.notna() & (s != "") & (s != "nan")]; s = s.reindex(expr.index)
        g = expr.groupby(s.values).mean()
        nprobe = s.value_counts()
    else:
        g = expr.groupby(expr.index.astype(str)).mean(); nprobe = pd.Series(1, index=g.index)
    g = g[g.notna().sum(axis=1) > 0]
    pct = g.rank(axis=0, pct=True) * 100
    out = []
    for c in CLOCK:
        alt = [c] + [k for k, v in ALIAS.items() if v == c]
        hit = next((a for a in alt if a in g.index), None)
        out.append(dict(gene=c, present=hit is not None, n_probes=int(nprobe.get(hit, 0)) if hit else 0,
                        median_percentile=float(np.nanmedian(pct.loc[hit].values)) if hit else np.nan))
    return pd.DataFrame(out), g.shape[0]
