"""C6 blind cohorts: header-only summary (no expression value is read).
Reads from each GEO series matrix ONLY the '!' header lines (stopping at '!series_matrix_table_begin') and the first
column (ID_REF) of the table, never the value columns; probe-to-symbol maps come from GEO platform tables (GPL text pages /
.annot.gz). Output: ../output/header_summary.csv (+ header_overlap.csv).
Run: python -I header_c6.py   (parses downloaded text only; nothing downloaded is executed)"""
import csv
import gzip
import io
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "output"
sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts"))
import a168_lib as L  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402

D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
PLAT = {"GPL570": D / "scout_luad_lihc/platforms/GPL570.annot.gz", "GPL10558": D / "scout_luad_lihc/platforms/GPL10558.annot.gz",
        "GPL25318": D / "c6_platforms/GPL25318.txt", "GPL18109": D / "c6_platforms/GPL18109.txt"}
SYMCOL = ("Gene symbol", "Gene Symbol")
G14 = list(L.GENES)  # the pipeline's 14 clock genes (ARNTL CLOCK NPAS2 PER1-3 CRY1-2 NR1D1 NR1D2 DBP RORA TEF CIART)
EXTRA = ["BHLHE40", "BHLHE41"]  # in the brief's list, not in the pipeline's 14 (reported only)
G2M, PURITY = set(L.sh.hallmark_g2m()), list(L.sh.PURITY)


def tumour_rule(key):
    """(description of annotation field, tumour predicate on one sample's {field: value}, normal predicate)."""
    return {
        "gse39582": ("characteristics 'dataset' in {discovery, validation}; normal = 'Non Tumoral'",
                     lambda s: s.get("dataset") in ("discovery", "validation"), lambda s: s.get("dataset") == "Non Tumoral"),
        "gse17536": ("all samples (source 'Moffitt Cancer Center'; no tumour-type field; primary colon cancers by design)",
                     lambda s: True, lambda s: False),
        "gse116918": ("characteristics 'tissue' == 'Primary Prostate Cancer' (diagnostic biopsies)",
                      lambda s: s.get("tissue") == "Primary Prostate Cancer", lambda s: False),
        "gse65858": ("characteristics 'tumor_type' == 'Primary' (Secondary and Relapse dropped)",
                     lambda s: s.get("tumor_type") == "Primary", lambda s: False),
        "gse53625": ("source_name starts 'cancer tissue' (paired 'normal tissue' dropped)",
                     lambda s: s["_source"].startswith("cancer tissue"), lambda s: s["_source"].startswith("normal tissue")),
        "gse157009": ("characteristics 'diagnosis' == 'Squamous Cell Carcinoma of Lung'",
                      lambda s: s.get("diagnosis") == "Squamous Cell Carcinoma of Lung", lambda s: False),
        "gse157010": ("characteristics 'diagnosis' == 'Squamous Cell Carcinoma of Lung'",
                      lambda s: s.get("diagnosis") == "Squamous Cell Carcinoma of Lung", lambda s: False),
    }[key]


COHORTS = [("colorectal", "gse39582", "GSE39582", "c6_colorectal_gse39582"), ("colorectal", "gse17536", "GSE17536", "c6_colorectal_gse17536"),
           ("prostate", "gse116918", "GSE116918", "c6_prostate_gse116918"), ("hnsc", "gse65858", "GSE65858", "c6_hnsc_gse65858"),
           ("escc", "gse53625", "GSE53625", "c6_escc_gse53625"), ("lusc", "gse157009", "GSE157009", "c6_lusc_gse157009"),
           ("lusc", "gse157010", "GSE157010", "c6_lusc_gse157010")]


def read_header(path):
    """Header '!' lines and the ID_REF column. The value columns are never parsed: each table line is cut at its first tab."""
    head, ids = [], []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            head.append(line.rstrip("\n"))
        next(fh)  # the ID_REF / GSM header line
        for line in fh:
            if line.startswith("!series_matrix_table_end"):
                break
            ids.append(line.split("\t", 1)[0].strip('"'))
    return head, ids


def samples(head):
    """Per-sample dicts from the header: gsm, title, _source and the 'key: value' characteristics."""
    rows = {}
    for line in head:
        p = next(csv.reader([line], delimiter="\t", quotechar='"'), [""])
        if not p or not p[0]:
            continue
        rows.setdefault(p[0], []).append(p[1:])
    gsm = rows["!Sample_geo_accession"][0]
    out = [{"gsm": g, "title": t, "_source": s} for g, t, s in zip(gsm, rows["!Sample_title"][0], rows["!Sample_source_name_ch1"][0])]
    for vals in rows.get("!Sample_characteristics_ch1", []):
        for d, v in zip(out, vals):
            if ":" in v:
                k, x = v.split(":", 1); d.setdefault(k.strip(), x.strip())
    one = lambda k: sorted({v for r in rows.get(k, []) for v in r})
    return out, one("!Sample_platform_id"), one("!Sample_data_processing"), rows


def platform_symbols(gpl):
    """probe ID -> canonical symbols (list) from the platform table; None if the table has no gene symbol column."""
    f = PLAT[gpl]
    with (gzip.open(f, "rt", errors="replace") if f.suffix == ".gz" else open(f, errors="replace")) as fh:
        lines = [l for l in fh if not l.startswith(("^", "!", "#"))]
    a = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str, low_memory=False)
    col = next((c for c in SYMCOL if c in a.columns), None)
    if col is None:
        return None
    a = a[["ID", col]].dropna()
    a["sym"] = a[col].str.split("///")
    a = a.explode("sym"); a["sym"] = a.sym.str.strip()
    a = a[~a.sym.isin(["", "---", "NA"])]
    a["sym"] = a.sym.map(lambda s: L.canonical(s) if L.canonical(s) in L.GENES else s)
    return a


def main():
    OUT.mkdir(exist_ok=True)
    rows, gsms, heads = [], {}, {}
    for typ, key, gse, folder in COHORTS:
        f = next((D / folder).glob("*_series_matrix.txt.gz"))
        head, ids = read_header(f)
        s, gpl, proc, _ = samples(head)
        keep, norm = tumour_rule(key)[1:]
        n_t, n_n = sum(map(keep, s)), sum(map(norm, s))
        gsms[gse] = {d["gsm"] for d in s}
        a = platform_symbols(gpl[0])
        ids = set(ids)
        if a is None:
            cov = dict(n_probes_in_matrix=len(ids), n_symbols=0, clock_present="", n_clock=0, n_g2m=0, n_purity=0, pseudo_pool=0, symbol_map="none: GEO platform table has no gene symbols")
        else:
            a = a[a.ID.isin(ids)]
            syms = set(a.sym)
            pool = [g for g in syms if L.canonical(g) not in CLOCK_LIKE and g not in G2M and g not in PURITY and g not in L.GENES]
            present = [g for g in G14 if g in syms]
            cov = dict(n_probes_in_matrix=len(ids), n_symbols=len(syms), clock_present=" ".join(present), n_clock=len(present),
                       n_g2m=len(G2M & syms), n_purity=len(set(PURITY) & syms), pseudo_pool=len(pool), symbol_map=f"{gpl[0]} platform table")
            cov["clock_missing"] = " ".join(g for g in G14 if g not in syms)
            cov["purity_missing"] = " ".join(g for g in PURITY if g not in syms)
            cov["extra_present"] = " ".join(g for g in EXTRA if g in syms)
        rows.append(dict(type=typ, gse=gse, platform=";".join(gpl), n_samples=len(s), n_primary_tumours=n_t, n_normals=n_n,
                         tumour_field=tumour_rule(key)[0], value_type=" | ".join(proc)[:300], n_g2m_set=len(G2M), n_purity_set=len(PURITY),
                         n_clock_set=len(G14), **cov))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "header_summary.csv", index=False)
    # sample-identity overlap from header accessions (GSE14333 accessions: GEO brief page, metadata only)
    g14 = {l.strip() for l in open(OUT / "gse14333_gsm_ids.txt") if l.strip()}
    gsms["GSE14333"] = g14
    ov = [dict(a=a, b=b, n_a=len(gsms[a]), n_b=len(gsms[b]), shared_gsm=len(gsms[a] & gsms[b])) for a in gsms for b in gsms if a < b]
    pd.DataFrame(ov).to_csv(OUT / "header_overlap.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
    print(df.drop(columns=["tumour_field", "value_type", "clock_present"]).to_string())
    print(pd.DataFrame(ov).to_string())


if __name__ == "__main__":
    main()
