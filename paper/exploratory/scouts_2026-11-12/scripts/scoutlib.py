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

def _table(p, begin):
    """Read the first platform table: header line after `begin`, rows until '!platform_table_end'."""
    lines = []
    with _open(p) as f:
        for l in f:
            if l.startswith(begin): break
        for l in f:
            if l.startswith("!platform_table_end"): break
            lines.append(l)
    a = pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0, low_memory=False, quoting=csv.QUOTE_NONE, on_bad_lines="skip")
    return a

def read_annot(p):
    """probe -> Series of symbols (list) from a GEO .annot.gz ('Gene symbol', '///'-separated) or a GPL soft family file
    (Affymetrix gene_assignment: 'acc // SYMBOL // title // ...' entries separated by ' /// ')."""
    a = _table(p, "!platform_table_begin")
    if "Gene symbol" in a.columns:
        return a["Gene symbol"].dropna().astype(str).str.split("///").apply(lambda L: [x.strip() for x in L if x.strip()])
    if "gene_assignment" in a.columns:
        g = a["gene_assignment"].dropna().astype(str)
        return g.apply(lambda s: sorted({e.split("//")[1].strip() for e in s.split("///") if "//" in e and len(e.split("//")) > 1}))
    raise KeyError(f"no symbol column in {p}: {list(a.columns)[:12]}")

def _present_clock(symset):
    out = {}
    for c in CLOCK:
        alt = [c] + [k for k, v in ALIAS.items() if v == c]
        out[c] = next((a for a in alt if a in symset), None)
    return out

def presence_only(sym):
    """sym: probe -> list of symbols. Returns clock presence / probe counts (no expression values)."""
    cnt = {}
    for L in sym.values:
        for x in L: cnt[x] = cnt.get(x, 0) + 1
    hit = _present_clock(set(cnt))
    return pd.DataFrame([dict(gene=c, present=hit[c] is not None, n_probes=cnt.get(hit[c], 0) if hit[c] else 0,
                              median_percentile=np.nan) for c in CLOCK]), len(cnt)

def clock_report(expr, sym=None):
    """expr: probes/genes x samples; sym: Series probe -> list of symbols, or None when the index is already gene symbols.
    Collapses probes to genes by the mean, ranks each sample across all genes, reports median percentile over samples."""
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if sym is not None:
        rows = []
        s = sym.reindex(expr.index).dropna()
        ex = pd.DataFrame({"probe": np.repeat(s.index.values, s.str.len().values), "gene": np.concatenate(s.values) if len(s) else []})
        nprobe = ex.gene.value_counts()
        vals = expr.loc[ex.probe.values].copy(); vals.index = ex.gene.values
        g = vals.groupby(level=0).mean()
    else:
        g = expr.groupby(expr.index.astype(str)).mean(); nprobe = pd.Series(1, index=g.index)
    g = g[g.notna().sum(axis=1) > 0]
    pct = g.rank(axis=0, pct=True) * 100
    hit = _present_clock(set(g.index))
    out = [dict(gene=c, present=hit[c] is not None, n_probes=int(nprobe.get(hit[c], 0)) if hit[c] else 0,
                median_percentile=float(np.nanmedian(pct.loc[hit[c]].values)) if hit[c] else np.nan) for c in CLOCK]
    return pd.DataFrame(out), g.shape[0]

def read_matrix_rows(p, keep):
    """Stream a series matrix and return only rows whose probe id is in `keep` (set of str)."""
    rows = {}
    with _open(p) as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"): break
        header = next(f).rstrip("\n").split("\t")
        for l in f:
            if l.startswith("!series_matrix_table_end"): break
            k, _, v = l.partition("\t")
            k = k.strip('"')
            if k in keep: rows[k] = [float(x) if x not in ('', 'null', 'NA', '"null"') else np.nan for x in v.rstrip("\n").split("\t")]
    return pd.DataFrame(rows, index=[h.strip('"') for h in header[1:]]).T

def ensembl_map():
    """ENSG (versionless) -> gene_name from the local GENCODE v29 GTF (parsed as text)."""
    m = {}
    with _open(f"{ROOT}/gencode_v29/gencode.v29.basic.annotation.gtf.gz") as f:
        for l in f:
            if l.startswith("#"): continue
            x = l.split("\t", 8)
            if len(x) < 9 or x[2] != "gene": continue
            gi = re.search(r'gene_id "([^"]+)"', x[8]); gn = re.search(r'gene_name "([^"]+)"', x[8])
            if gi and gn: m[gi.group(1).split(".")[0]] = gn.group(1)
    return pd.Series(m)
