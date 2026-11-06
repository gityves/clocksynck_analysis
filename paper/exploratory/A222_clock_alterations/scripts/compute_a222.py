"""A222: are clock genes selected for alteration in TCGA tumours, or passengers? (registered in
paper/docs/registered_experiments.md, A222).

H222.1 mutation: R = log[(non-silent + 0.5) / (silent + 0.5)] summed over the 14 waveform clock genes (MC3 PASS-class calls,
hypermutators > 1,000 coding mutations excluded); null 1,000 random 14-gene sets, each clock gene replaced by a gene of the same
CDS-length decile x expression decile cell. H222.2 / H222.3 copy number: per gene and project, deep-deletion (-2) / amplification
(+2) frequency minus the median over its 25 flanking protein-coding genes on each side (same chromosome, GENCODE v29 order); mean
over clock genes and projects; null 1,000 random 14-gene sets of protein-coding genes with a full flank.

Stages: `power` writes the gene tables and the null distributions only (no clock-gene statistic); `tests` scores.
    python compute_a222.py --stage power|tests
Outputs: ../output/{genes.csv (background only in power), null_*.csv, power.csv, tests.csv, descriptive_*.csv}
"""
from __future__ import annotations

import argparse
import gzip
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "paper"))
import clocksync as cs  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parent / "output"
GTF = DATA / "gencode_v29" / "gencode.v29.basic.annotation.gtf.gz"
MC3 = DATA / "tcga" / "mc3" / "mc3.v0.2.8.PUBLIC.maf.gz"
GISTIC = DATA / "tcga" / "xena_gistic" / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
RC3 = DATA / "tcga" / "recount3"
MC3_FILTER_OK = {"PASS", "wga", "native_wga_mix"}  # as A215
NONSILENT = {"Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
             "In_Frame_Ins", "Splice_Site", "Translation_Start_Site"}
TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation"}
CODING = NONSILENT | {"Silent"}
HYPER, FLANK, NNULL, MIN_PROJ, N_EXPR = 1000, 25, 1000, 100, 30
ALIAS = {"CIART": "C1orf51"}  # MC3 / Xena GISTIC predate the rename
CLOCK = list(cs.gtex().genes)
SEED = zlib.crc32(b"A222") % (2 ** 31)
K = 2.49
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


# ---- gene table ------------------------------------------------------------------------------------------------------------
def gencode() -> pd.DataFrame:
    """Protein-coding genes: chromosome, start, longest CDS (bp) over basic transcripts."""
    genes, cds = {}, {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 8)
            if c[2] not in ("gene", "CDS"):
                continue
            a = dict(kv.strip().split(" ", 1) for kv in c[8].strip().rstrip(";").split(";") if kv.strip())
            name, gtype = a["gene_name"].strip('"'), a["gene_type"].strip('"')
            if gtype != "protein_coding" or c[0] in ("chrY", "chrM"):
                continue
            if c[2] == "gene":
                genes[name] = (c[0], int(c[3]), a["gene_id"].strip('"').split(".")[0])
            else:
                t = a["transcript_id"].strip('"')
                cds.setdefault(name, {}).setdefault(t, 0)
                cds[name][t] += int(c[4]) - int(c[3]) + 1
    g = pd.DataFrame([(n, *v, max(cds.get(n, {0: 0}).values())) for n, v in genes.items()],
                     columns=["gene", "chrom", "start", "ensg", "cds_bp"])
    return g[g.cds_bp > 0].reset_index(drop=True)


def project_of_patient() -> dict:
    m = {}
    for p in sorted(x.name for x in RC3.iterdir() if (x / f"tcga.tcga.{x.name}.MD.gz").exists()):
        md = pd.read_csv(RC3 / p / f"tcga.tcga.{p}.MD.gz", sep="\t", usecols=["tcga_barcode"])
        m.update({b[:12]: p for b in md.tcga_barcode})
    return m


def expression(g: pd.DataFrame) -> pd.Series:
    """Pan-TCGA median log2 CPM (protein-coding library) over N_EXPR primary tumours per project (first by barcode)."""
    ensg = set(g.ensg)
    med = []
    for p in sorted(x.name for x in RC3.iterdir() if (x / f"tcga.tcga.{x.name}.MD.gz").exists()):
        md = pd.read_csv(RC3 / p / f"tcga.tcga.{p}.MD.gz", sep="\t", usecols=["external_id", "tcga_barcode", "cgc_sample_sample_type"])
        md = md[md.cgc_sample_sample_type == "Primary Tumor"].sort_values("tcga_barcode").head(N_EXPR)
        x = pd.read_csv(RC3 / p / f"tcga.gene_sums.{p}.G029.gz", sep="\t", comment="#", usecols=["gene_id", *md.external_id])
        x.index = x.pop("gene_id").str.split(".").str[0]
        x = x[x.index.isin(ensg)]
        x = x[~x.index.duplicated()]
        med.append(np.log2(x / x.sum() * 1e6 + 1).median(1).rename(p))
        say("expression", p, len(md))
    e = pd.concat(med, axis=1).median(1)
    return g.ensg.map(e)


def mutations(g: pd.DataFrame, proj: dict) -> tuple[pd.DataFrame, dict]:
    sym = {ALIAS.get(n, n): n for n in g.gene}
    cols = ["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "FILTER"]
    rows = []
    for ch in pd.read_csv(MC3, sep="\t", usecols=cols, chunksize=500_000, dtype=str):
        ch = ch[ch.Variant_Classification.isin(CODING)]
        ch = ch[ch.FILTER.fillna("").map(lambda x: set(x.split(";")) <= MC3_FILTER_OK)]
        rows.append(ch.assign(s=ch.Tumor_Sample_Barcode.str[:15]).drop(columns=["FILTER", "Tumor_Sample_Barcode"]))
    m = pd.concat(rows)
    load = m.groupby("s").size()
    keep = load.index[load <= HYPER]
    info = dict(samples=int(len(load)), hypermutators=int((load > HYPER).sum()), kept=int(len(keep)))
    m = m[m.s.isin(keep) & m.Hugo_Symbol.isin(sym)]
    m["gene"] = m.Hugo_Symbol.map(sym)
    t = pd.DataFrame({"nonsilent": m[m.Variant_Classification.isin(NONSILENT)].groupby("gene").size(),
                      "silent": m[m.Variant_Classification == "Silent"].groupby("gene").size(),
                      "truncating": m[m.Variant_Classification.isin(TRUNC)].groupby("gene").size()}).reindex(g.gene).fillna(0)
    return t.reset_index(drop=True), info


def copy_number(g: pd.DataFrame, proj: dict) -> tuple[dict, dict]:
    """Per project: deep-deletion and amplification frequency per gene (rows aligned to g, NaN when absent)."""
    x = pd.read_csv(GISTIC, sep="\t", index_col=0)
    x = x[~x.index.duplicated()]
    x = x.reindex([ALIAS.get(n, n) for n in g.gene])
    pj = pd.Series([proj.get(c[:12]) for c in x.columns], index=x.columns)
    dele, amp, n = {}, {}, {}
    for p, cols in pj.dropna().groupby(pj.dropna()).groups.items():
        if len(cols) < MIN_PROJ:
            continue
        v = x[cols].to_numpy(float)
        ok = ~np.isnan(v)
        dele[p] = np.where(ok.any(1), (v == -2).sum(1) / np.maximum(ok.sum(1), 1), np.nan)
        amp[p] = np.where(ok.any(1), (v == 2).sum(1) / np.maximum(ok.sum(1), 1), np.nan)
        n[p] = len(cols)
    return dict(del_=pd.DataFrame(dele), amp=pd.DataFrame(amp)), n


def focal_excess(f: pd.DataFrame, g: pd.DataFrame) -> pd.DataFrame:
    """Frequency minus the median over FLANK genes on each side (same chromosome, start order); NaN without a full flank."""
    out = pd.DataFrame(np.nan, index=f.index, columns=f.columns)
    for ch, idx in g.sort_values("start").groupby("chrom").groups.items():
        idx = np.asarray(g.loc[idx].sort_values("start").index)
        v = f.loc[idx].to_numpy()
        for i in range(FLANK, len(idx) - FLANK):
            nb = np.r_[v[i - FLANK:i], v[i + 1:i + 1 + FLANK]]
            out.loc[idx[i]] = v[i] - np.nanmedian(nb, axis=0)
    return out


# ---- statistics and nulls ---------------------------------------------------------------------------------------------------
def R_stat(t, idx):
    return float(np.log((t.nonsilent.to_numpy()[idx].sum() + 0.5) / (t.silent.to_numpy()[idx].sum() + 0.5)))


def matched_null(g, clock_idx, rng, n=NNULL):
    """Random sets: each clock gene replaced by a non-clock gene of its CDS-length decile x expression decile cell."""
    cell = g.cds_dec.astype(str) + "_" + g.expr_dec.astype(str)
    pool = {c: np.flatnonzero((cell == c).to_numpy() & ~g.is_clock.to_numpy()) for c in cell.unique()}
    return np.array([[rng.choice(pool[cell.iloc[i]]) for i in clock_idx] for _ in range(n)])


def cn_stat(e, idx):
    return float(np.nanmean(e.to_numpy()[idx]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", required=True, choices=["power", "tests"])
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "gene_cache.pkl"
    if cache.exists():
        g, t, cn, ncn, info = pd.read_pickle(cache)
    else:
        g = gencode(); say("genes", len(g))
        proj = project_of_patient()
        g["expr"] = expression(g)
        t, info = mutations(g, proj); say("mutations", info)
        cn, ncn = copy_number(g, proj); say("copy number projects", len(ncn))
        pd.to_pickle((g, t, cn, ncn, info), cache)
    missing = [c for c in CLOCK if c not in set(g.gene)]
    assert not missing, missing
    g["is_clock"] = g.gene.isin(CLOCK)
    g = g[g.expr.notna()].copy()
    g["cds_dec"] = pd.qcut(g.cds_bp, 10, labels=False)
    g["expr_dec"] = pd.qcut(g.expr.rank(method="first"), 10, labels=False)
    t = t.loc[g.index]
    clock_idx = np.flatnonzero(g.is_clock.to_numpy())
    rng = np.random.default_rng(SEED)

    null_mut = matched_null(g.reset_index(drop=True), clock_idx, rng)
    tr = t.reset_index(drop=True)
    nR = np.array([R_stat(tr, s) for s in null_mut])
    ex = {k: focal_excess(v.loc[g.index], g) for k, v in cn.items()}
    full = np.flatnonzero(ex["del_"].notna().all(1).to_numpy() & ~g.is_clock.to_numpy())
    null_cn = np.array([rng.choice(full, len(CLOCK), replace=False) for _ in range(NNULL)])
    nD = np.array([cn_stat(ex["del_"], s) for s in null_cn])
    nA = np.array([cn_stat(ex["amp"], s) for s in null_cn])
    pd.DataFrame({"R": nR, "del_focal": nD, "amp_focal": nA}).to_csv(OUT / "null_draws.csv", index=False)
    pw = pd.DataFrame([
        dict(test="H222.1", stat="R log(nonsilent/silent)", null_mean=nR.mean(), null_sd=nR.std(ddof=1), mde=K * nR.std(ddof=1)),
        dict(test="H222.2", stat="deep-deletion focal excess", null_mean=nD.mean(), null_sd=nD.std(ddof=1), mde=K * nD.std(ddof=1)),
        dict(test="H222.3", stat="amplification focal excess", null_mean=nA.mean(), null_sd=nA.std(ddof=1), mde=K * nA.std(ddof=1))])
    elig = dict(**info, genes=len(g), cn_projects=len(ncn), cn_samples=int(sum(ncn.values())),
                clock_with_flank=int(ex["del_"].loc[g.index[clock_idx]].notna().all(1).sum()))
    pw = pw.assign(**{k: v for k, v in elig.items()})
    pw.to_csv(OUT / "power.csv", index=False)
    say("\n" + pw.round(5).to_string(index=False))
    if a.stage == "power":
        return

    # ---- tests -------------------------------------------------------------------------------------------------------------
    R = R_stat(tr, clock_idx)
    D = cn_stat(ex["del_"], clock_idx)
    A = cn_stat(ex["amp"], clock_idx)
    p1 = (1 + (nR >= R).sum()) / (1 + len(nR))
    p2 = (1 + (nD >= D).sum()) / (1 + len(nD))
    p3 = 2 * min((1 + (nA >= A).sum()), (1 + (nA <= A).sum())) / (1 + len(nA))
    res = pd.DataFrame([dict(test="H222.1", T=R, p=p1, pctile=100 * (nR < R).mean()),
                        dict(test="H222.2", T=D, p=p2, pctile=100 * (nD < D).mean()),
                        dict(test="H222.3", T=A, p=min(p3, 1.0), pctile=100 * (nA < A).mean())])
    res = res.merge(pw[["test", "null_mean", "null_sd", "mde"]], on="test")
    order = np.argsort(res.p.to_numpy())
    holm = np.maximum.accumulate([min(1, (3 - k) * res.p.iloc[i]) for k, i in enumerate(order)])
    res.loc[res.index[order], "holm"] = holm
    res.to_csv(OUT / "tests.csv", index=False)
    say("\n" + res.round(5).to_string(index=False))

    # descriptive per clock gene
    gc = g.iloc[clock_idx]
    tc = t.loc[gc.index]
    bg = t.loc[g.index[~g.is_clock.to_numpy()]]
    gbg = g[~g.is_clock]
    rate_bg = (bg.nonsilent + bg.silent) / (gbg.cds_bp / 1e6)
    desc = pd.DataFrame(dict(gene=gc.gene, chrom=gc.chrom, cds_bp=gc.cds_bp, expr=gc.expr.round(2), nonsilent=tc.nonsilent, silent=tc.silent,
                             truncating=tc.truncating, per_mb=((tc.nonsilent + tc.silent) / (gc.cds_bp / 1e6)).round(1),
                             per_mb_matched_median=[float(rate_bg[(gbg.cds_dec == r.cds_dec) & (gbg.expr_dec == r.expr_dec)].median()) for r in gc.itertuples()],
                             del_freq_mean=cn["del_"].loc[gc.index].mean(1).round(4), del_focal_mean=ex["del_"].loc[gc.index].mean(1).round(4),
                             amp_freq_mean=cn["amp"].loc[gc.index].mean(1).round(4), amp_focal_mean=ex["amp"].loc[gc.index].mean(1).round(4)))
    desc.to_csv(OUT / "descriptive_genes.csv", index=False)
    say("\n" + desc.to_string(index=False))


if __name__ == "__main__":
    main()
