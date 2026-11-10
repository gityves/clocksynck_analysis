"""Kang et al. 2011 (Nature 478:483), human brain development, GEO GSE25219, Affymetrix Human Exon 1.0 ST (GPL5175), gene-level matrix: prep ONLY.
NO clock statistic (no correlation, phase, coherence, gate or fit): design tables and detection-level summaries.

Source (untrusted, parse only): `GSE25219-GPL5175_series_matrix.txt.gz` (1,340 samples x 17,565 transcript clusters, Partek-normalised: the series
header gives the processing; the value scale is reported at run time) and the platform table `GPL5175_full.txt` (`gene_assignment` column).

Collapse to symbols: the symbol of each transcript cluster is the second field of every `///`-separated `gene_assignment` entry; clusters whose entries name
more than one distinct symbol, or none, are dropped; clusters sharing a symbol are averaged (mean on the deposited scale). Counts are printed.

Ages: `age_pcw` = post-conception weeks (prenatal as given in PCW; postnatal (`PMonth`, `Y` in the series) 40 + weeks since birth: months x 4.348, years x 52.18);
`age_log` = log2(age_pcw); `stage` = the paper's period (1-15).

    python -I prep_kang.py
"""
from __future__ import annotations

import gzip
import io
import time
from pathlib import Path

import numpy as np
import pandas as pd

from clocksync.genes import CLOCK_GENES, canonical

DAT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/natexp/kang2011")
RAW = DAT / "raw"
OUT = Path(__file__).resolve().parents[1] / "output"
GROUPS = [("pre_le12", 0, 12.99), ("pre_13_22", 13, 22.99), ("pre_25_37", 23, 40), ("infant_lt1y", 40, 40 + 52.18), ("child_1_12y", 40 + 52.18, 40 + 52.18 * 13),
          ("ge13y", 40 + 52.18 * 13, 1e9)]


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def read_series():
    meta, rows, table, in_t = {}, [], [], False
    with gzip.open(RAW / "GSE25219-GPL5175_series_matrix.txt.gz", "rt", errors="replace") as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"):
                in_t = True
            elif l.startswith("!series_matrix_table_end"):
                in_t = False
            elif in_t:
                table.append(l)
            elif l.startswith("!Sample_"):
                k, *v = l.rstrip("\n").split("\t")
                rows.append((k[len("!Sample_"):], [x.strip('"') for x in v]))
    gsm = dict(rows)["geo_accession"]
    m = pd.DataFrame({"gsm": gsm, "title": dict(rows)["title"]})
    for k, v in rows:
        if k == "characteristics_ch1":
            key = v[0].partition(":")[0].strip().lower()
            m[key] = [x.partition(":")[2].strip() for x in v]
    ex = pd.read_csv(io.StringIO("".join(table)), sep="\t", index_col=0, quotechar='"')
    ex.index = ex.index.astype(str)
    return m, ex.astype("float32")


def cluster_symbols():
    with open(RAW / "GPL5175_full.txt", errors="replace") as fh:
        t = pd.read_csv(io.StringIO("".join(l for l in fh if not l.startswith(("^", "!", "#")))), sep="\t", dtype=str, usecols=["ID", "gene_assignment"])
    ga = t.set_index("ID").gene_assignment.fillna("---")
    def syms(x):
        if x == "---":
            return set()
        return {canonical(e.split(" // ")[1]) for e in x.split(" /// ") if len(e.split(" // ")) > 1 and e.split(" // ")[1].strip() not in ("", "---")}
    S = ga.map(syms)
    return S, ga


def prep():
    m, ex = read_series()
    say("series matrix", ex.shape, "samples in header", len(m), "header/table columns equal:", list(ex.columns) == list(m.gsm))
    m = m.rename(columns={"brain code": "brain", "sex": "sex", "postmortem interval": "pmi", "ph": "ph", "rna integrity number": "rin", "age": "age_string"})
    m["stage"] = m["stage"].astype(int)
    m["prenatal"] = m.age_string.str.contains("PCW")
    m["age_pcw"] = [float(a.split()[0]) if p else (40 + float(a.split()[0]) * {"PMonth": 4.348, "Y": 52.18}.get(a.split()[1], np.nan))
                    for a, p in zip(m.age_string, m.prenatal)]
    assert m.age_pcw.notna().all(), m.loc[m.age_pcw.isna(), "age_string"].unique()
    m["age_log"] = np.log2(m.age_pcw)
    m["age_value"] = m.age_string.str.split().str[0].astype(float)
    m["age_unit"] = m.age_string.str.split().str[1]
    for c in ("pmi", "ph", "rin"):
        m[c] = pd.to_numeric(m[c], errors="coerce")
    m["group"] = [next(g for g, a, b in GROUPS if a <= x <= b) for x in m.age_pcw]
    cols = ["gsm", "title", "brain", "region", "hemisphere", "sex", "age_string", "age_value", "age_unit", "age_pcw", "age_log", "prenatal", "stage", "group", "pmi", "ph", "rin"]
    m[cols].to_csv(DAT / "samples.csv", index=False)
    m[cols].to_csv(OUT / "kang_samples.csv", index=False)
    say("brains", m.brain.nunique(), "prenatal", m[m.prenatal].brain.nunique(), "postnatal", m[~m.prenatal].brain.nunique())
    v = ex.to_numpy()
    say("value scale (all clusters x samples): min", float(np.nanmin(v)), "median", float(np.nanmedian(v)), "max", float(np.nanmax(v)), "NaN frac", float(np.isnan(v).mean()))
    ex.to_parquet(DAT / "expr_cluster.parquet")
    S, ga = cluster_symbols()
    S = S.reindex(ex.index)
    n0 = int(S.map(len).eq(0).sum()); nmulti = int(S.map(len).gt(1).sum())
    keep = S.map(len).eq(1)
    sym = S[keep].map(lambda s: next(iter(s)))
    say("clusters", len(ex), "no symbol", n0, "several symbols (dropped)", nmulti, "single symbol", int(keep.sum()), "distinct symbols", sym.nunique(),
        "clusters per symbol >1:", int((sym.value_counts() > 1).sum()))
    g = ex.loc[sym.index].groupby(sym.to_numpy()).mean().astype("float32")
    g.to_parquet(DAT / "expr_gene.parquet")
    say("expr_gene", g.shape)
    design(m)
    detection(m, ex, g, S, sym)


def design(m):
    br = m.drop_duplicates("brain")
    t = br.groupby(["stage"]).agg(brains=("brain", "size"), age_min_pcw=("age_pcw", "min"), age_max_pcw=("age_pcw", "max"), prenatal=("prenatal", "first"))
    t.to_csv(OUT / "design_kang_brains_per_stage.csv")
    print("\nbrains per stage\n", t.to_string())
    g = br.groupby("group").size().reindex([x[0] for x in GROUPS]); g.to_csv(OUT / "design_kang_brains_per_group.csv")
    print("\nbrains per group\n", g.to_string())
    rpb = m.groupby("brain").agg(n_samples=("gsm", "size"), n_regions=("region", "nunique"), group=("group", "first"))
    d = rpb.groupby("group").n_regions.describe().reindex([x[0] for x in GROUPS]); d.to_csv(OUT / "design_kang_regions_per_brain.csv")
    print("\nregions per brain by group\n", d.round(1).to_string())
    print("\nregions-per-brain distribution (all)\n", rpb.n_regions.value_counts().sort_index().to_string())
    rs = m.groupby("stage").region.apply(lambda x: " ".join(sorted(set(x)))); rs.to_csv(OUT / "design_kang_region_labels_per_stage.csv")
    print("\nregion labels per stage\n", "\n".join(f"{k}: {v}" for k, v in rs.items()))
    bg = m.drop_duplicates(["brain", "region"]).groupby(["region", "group"]).brain.nunique().unstack(fill_value=0)[[x[0] for x in GROUPS]]
    bg.to_csv(OUT / "design_kang_brains_per_region_group.csv")
    print("\nbrains having each region, by group\n", bg.to_string())
    h = m.groupby(["brain", "region"]).hemisphere.nunique()
    both = h[h > 1].reset_index().groupby("brain").size().reindex(rpb.index, fill_value=0)
    hb = pd.DataFrame({"n_regions_both_hemispheres": both, "n_samples": rpb.n_samples})
    hb.to_csv(OUT / "design_kang_both_hemisphere_per_brain.csv")
    print("\nhemispheres:", m.hemisphere.value_counts().to_dict(), "; brains with a region sampled in both hemispheres:", int((both > 0).sum()),
          "; brains with samples in both hemispheres:", int((m.groupby("brain").hemisphere.nunique() > 1).sum()))
    print("sex per brain:", br.sex.value_counts().to_dict())


def detection(m, ex, g, S, sym):
    from clocksync import gtex
    print("\nclocksync CLOCK_GENES == gtex().genes:", list(CLOCK_GENES) == list(gtex().genes))
    flat = ex.to_numpy()
    q = np.quantile(flat, [.05, .25, .5, .75, .95]); say("background, all clusters x samples, quantiles 5/25/50/75/95:", q.round(2).tolist())
    cm = ex.median(axis=1)
    say("per-cluster median across samples, quantiles 5/25/50/75/95:", np.quantile(cm, [.05, .25, .5, .75, .95]).round(2).tolist())
    pre, post = m.gsm[m.prenatal], m.gsm[~m.prenatal]
    rows = []
    for gene in CLOCK_GENES:
        cl = [c for c in S.index if gene in S[c]]
        single = [c for c in cl if c in sym.index]
        r = {"gene": gene, "n_clusters_naming_gene": len(cl), "n_clusters_single_gene": len(single), "in_expr_gene": gene in g.index,
             "clusters_naming_gene": ";".join(cl)}
        if gene in g.index:
            r["median_log2_prenatal"] = float(g.loc[gene, pre].median()); r["median_log2_postnatal"] = float(g.loc[gene, post].median())
            r["pctl_of_cluster_medians_prenatal"] = float((cm.to_numpy() < r["median_log2_prenatal"]).mean() * 100)
            r["pctl_of_cluster_medians_postnatal"] = float((cm.to_numpy() < r["median_log2_postnatal"]).mean() * 100)
        if cl:  # cluster level (includes clusters shared with other genes, dropped from expr_gene): best-expressed cluster naming the gene
            c = max(cl, key=lambda c: ex.loc[c].median())
            r.update(cluster_used=c, cluster_genes="/".join(sorted(S[c])), cluster_median_log2_prenatal=float(ex.loc[c, pre].median()),
                     cluster_median_log2_postnatal=float(ex.loc[c, post].median()))
        rows.append(r)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "clock_detection_kang.csv", index=False); d.to_csv(DAT / "clock_detection.csv", index=False)
    print(d.drop(columns="clusters_naming_gene").round(2).to_string())


if __name__ == "__main__":
    prep()
