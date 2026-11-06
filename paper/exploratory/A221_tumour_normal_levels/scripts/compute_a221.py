"""A221: are "clock genes downregulated in tumours" real or a design artefact? (registered in
paper/docs/registered_experiments.md, A221).
  --stage power : download-independent data build (cached), eligibility, 200 random decile-matched 14-gene sets per project,
                  null SD of L and Delta from the random sets only, MDE = 2.49 x null SD (eligibility.csv, power.csv,
                  random_sets.csv, null_draws.csv). No clock-gene tumour-vs-normal difference is computed or written.
  --stage tests : H221.1-3 with Holm, descriptive tables (NOT run before the power note is dated and committed).
    python -I compute_a221.py --stage power|tests [--jobs 4] [--smoke --out DIR]
Time of death and Hardy class are controlled (dbGaP): held in memory, only aggregates are written."""
import argparse
import gzip
import os
import sys
import zlib
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "paper"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm, spearmanr, ttest_ind  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import gene_columns  # noqa: E402
from common.datasets.tcga_paired import build as TB  # noqa: E402
from common.paths import DATA  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--jobs", type=int, default=4)
ap.add_argument("--nsets", type=int, default=200)
ap.add_argument("--nperm", type=int, default=2000)
ap.add_argument("--smoke", action="store_true", help="two small projects, 20 sets; use with --out")
ap.add_argument("--out", type=Path, default=HERE.parent / "output")
args = ap.parse_args()
args.jobs = min(args.jobs, 8)
OUT = args.out; OUT.mkdir(parents=True, exist_ok=True)
GT = DATA / "gtex_recount3"
CACHE = (OUT / "cache") if args.smoke else (GT / "cache_a221"); CACHE.mkdir(exist_ok=True)
K = norm.ppf(0.95) + norm.ppf(0.80)  # 2.49
MDE_MAX, MIN_PAIRS, MIN_GTEX = 0.5, 10, 50
SEED = zlib.crc32(b"A221")
GENES = list(cs.gtex().genes)
IEG = ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "EGR2", "EGR3", "IER2", "ATF3", "NR4A1"]
# TCGA project -> (recount3 GTEx SMTS file, SMTSD label)
MAP = {"BRCA": ("BREAST", "Breast - Mammary Tissue"), "LUAD": ("LUNG", "Lung"), "LUSC": ("LUNG", "Lung"),
       "KIRC": ("KIDNEY", "Kidney - Cortex"), "KIRP": ("KIDNEY", "Kidney - Cortex"), "KICH": ("KIDNEY", "Kidney - Cortex"),
       "LIHC": ("LIVER", "Liver"), "THCA": ("THYROID", "Thyroid"), "PRAD": ("PROSTATE", "Prostate"),
       "COAD": ("COLON", "Colon - Transverse"), "STAD": ("STOMACH", "Stomach"),
       "ESCA": ("ESOPHAGUS", "Esophagus - Mucosa"), "UCEC": ("UTERUS", "Uterus")}
PROJECTS = ["BRCA", "COAD", "ESCA", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC"]
if args.smoke:
    PROJECTS, args.nsets = ["KICH", "UCEC"], 20
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


# ---- data build (parse only) -------------------------------------------------------------------------------------------------------
def annotation() -> pd.DataFrame:
    """Protein-coding genes of GENCODE v29 (gene_id, symbol); duplicated symbols are summed downstream."""
    f = CACHE / "pc_genes.csv"
    if f.exists():
        return pd.read_csv(f)
    rows = []
    with gzip.open(TB.RECOUNT3 / "annotation" / "human.gene_sums.G029.gtf.gz", "rt") as h:
        for line in h:
            if line[0] == "#":
                continue
            c = line.split("\t", 8)
            if c[2] != "gene":
                continue
            a = dict(x.strip().split(" ", 1) for x in c[8].strip().strip(";").split(";") if x.strip())
            if a["gene_type"].strip('"') == "protein_coding":
                rows.append((a["gene_id"].strip('"'), a["gene_name"].strip('"')))
    d = pd.DataFrame(rows, columns=["gene_id", "symbol"]); d.to_csv(f, index=False)
    return d


def log_cpm(path, ids, ann) -> pd.DataFrame:
    """samples x protein-coding symbols of log2(CPM + 1), library = summed protein-coding coverage (recount3 gene_sums are
    base-pair coverage sums; the per-sample AUC scaling cancels in a CPM)."""
    want = set(ids)
    cov = pd.read_csv(path, sep="\t", comment="#", index_col=0, usecols=lambda c: c == "gene_id" or c in want)[list(ids)]
    cov = cov.reindex(ann.gene_id).dropna(how="all")
    cov.index = ann.set_index("gene_id").symbol.reindex(cov.index).to_numpy()
    cov = cov.groupby(level=0).sum().T
    lc = np.log2(cov.div(cov.sum(axis=1), axis=0) * 1e6 + 1)
    return lc.astype(np.float32)


def build_project(p: str, ann) -> Path:
    """Cache per project: log2 CPM of the genes with mean CPM >= 1 over tumours + normals (and the 14 clock genes) for
    primary tumours (one per patient), solid-tissue normals, the pairs, and the matched GTEx samples."""
    f = CACHE / f"{p}.npz"
    if f.exists():
        return f
    say(p, "build")
    s = TB.samples_of(p)
    t = s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")
    n = s[s.sample_type == "Solid Tissue Normal"].sort_values("barcode").drop_duplicates("patient")
    pts = sorted(set(t.patient) & set(n.patient))
    sm = pd.concat([t, n]).drop_duplicates("external_id")
    lc = log_cpm(TB.RECOUNT3 / p / f"tcga.gene_sums.{p}.G029.gz", sm.external_id.tolist(), ann)
    bar = sm.set_index("external_id")
    lt, ln = lc.loc[t.external_id], lc.loc[n.external_id]
    dir_, lab = MAP[p]
    md = pd.read_csv(GT / dir_ / f"gtex.gtex.{dir_}.MD.gz", sep="\t", low_memory=False)
    md = md[md.SMTSD == lab]
    lg = log_cpm(GT / dir_ / f"gtex.gene_sums.{dir_}.G029.gz", md.external_id.tolist(), ann).loc[md.external_id]
    genes = lt.columns.intersection(lg.columns)
    cpm = (2 ** pd.concat([lt, ln])[genes].astype(np.float64) - 1).mean()
    clock = [genes[i] for i in gene_columns(list(genes), GENES)]
    keep = genes[(cpm.to_numpy() >= 1) | genes.isin(clock)]
    ti = {x: i for i, x in enumerate(t.patient)}; ni = {x: i for i, x in enumerate(n.patient)}
    np.savez_compressed(
        f, genes=np.array(keep).astype(str), tum=lt[keep].to_numpy(), nor=ln[keep].to_numpy(), gtex=lg[keep].to_numpy(),
        pair_t=np.array([ti[x] for x in pts]), pair_n=np.array([ni[x] for x in pts]),
        tum_sample=bar.loc[t.external_id, "barcode"].to_numpy().astype(str),
        gtex_subject=md.SUBJID.to_numpy().astype(str), gtex_hardy=md.DTHHRDY.to_numpy(float))  # hardy: memory only, see load()
    return f


class Proj:
    """One project's expression matrices (in memory)."""

    def __init__(self, p):
        z = np.load(CACHE / f"{p}.npz", allow_pickle=False)
        self.p, self.genes = p, pd.Index(z["genes"])
        self.tum, self.nor, self.gtex = z["tum"], z["nor"], z["gtex"]
        self.pt, self.pn = z["pair_t"], z["pair_n"]
        self.tum_sample, self.subject, self.hardy = z["tum_sample"], z["gtex_subject"], z["gtex_hardy"]
        self.clock = np.array(gene_columns(list(self.genes), GENES))
        self.mean_expr = np.concatenate([self.tum, self.nor]).mean(0)  # mean log2 CPM over tumours and normals

    def lfc_paired(self):
        return (self.tum[self.pt] - self.nor[self.pn]).mean(0)

    def lfc_gtex(self, mask=None):
        g = self.gtex if mask is None else self.gtex[mask]
        return self.tum.mean(0) - g.mean(0)


def tod_by_subject() -> pd.Series:
    """Time of death (h, mod 24) per GTEx subject from the labelled h5ad obs. Controlled: in memory only."""
    import anndata as ad
    from common.gtex import LABELLED
    a = ad.read_h5ad(LABELLED, backed="r")
    try:
        o = a.obs[["subject", "tod"]].astype({"subject": str})
    finally:
        a.file.close()
    o = o.dropna().drop_duplicates("subject").set_index("subject").tod.astype(float) % 24
    return o


def is_day(tod): return (tod >= 8) & (tod < 18)


# ---- random sets ---------------------------------------------------------------------------------------------------------------------
def random_sets(P: Proj, nsets: int) -> np.ndarray:
    """nsets x 14 gene indices: each clock gene is replaced by a gene of its expression decile (deciles of the mean log2 CPM
    over tumours and normals, among expressed protein-coding genes outside the clock set); distinct within a set. The draw
    depends on the project and crc32('A221') only. Uses expression only, no tumour-vs-normal difference."""
    cand = np.setdiff1d(np.arange(len(P.genes)), P.clock)
    x = P.mean_expr
    edges = np.quantile(x[cand], np.linspace(0, 1, 11)[1:-1])
    dec = np.digitize(x, edges)
    pools = {d: cand[dec[cand] == d] for d in range(10)}
    need = pd.Series(dec[P.clock]).value_counts()
    rng = np.random.default_rng([SEED, zlib.crc32(P.p.encode())])
    return np.array([np.concatenate([rng.choice(pools[d], n, replace=False) for d, n in need.sort_index().items()])
                     for _ in range(nsets)])


# ---- power --------------------------------------------------------------------------------------------------------------------------
def eligibility(tod) -> pd.DataFrame:
    rows = []
    for p in PROJECTS:
        P = Proj(p); tt = tod.reindex(P.subject)
        rows.append(dict(project=p, pairs=len(P.pt), tumours=len(P.tum), normals=len(P.nor), gtex_tissue=MAP[p][1],
                         gtex_samples=len(P.gtex), gtex_with_tod=int(tt.notna().sum()),
                         gtex_daytime_08_18=int(is_day(tt).sum()),
                         eligible=len(P.pt) >= MIN_PAIRS and len(P.gtex) >= MIN_GTEX))
    return pd.DataFrame(rows)


def stage_power():
    ann = annotation()
    Parallel(n_jobs=args.jobs)(delayed(build_project)(p, ann) for p in PROJECTS)
    el = eligibility(tod_by_subject()); el.to_csv(OUT / "eligibility.csv", index=False); say("eligibility\n" + el.to_string())
    ok = el.project[el.eligible].tolist()
    m_pair, m_diff, sets_rows = [], [], []
    for p in ok:
        P = Proj(p); S = random_sets(P, args.nsets)
        lp, lg = P.lfc_paired(), P.lfc_gtex()
        m_pair.append(lp[S].mean(1)); m_diff.append((lg - lp)[S].mean(1))
        sets_rows += [(p, k, r, P.genes[g]) for k, s in enumerate(S) for r, g in enumerate(s)]
    pd.DataFrame(sets_rows, columns=["project", "set", "rank", "gene"]).to_csv(OUT / "random_sets.csv", index=False)
    m_pair, m_diff = np.array(m_pair), np.array(m_diff)  # projects x sets
    # null draw of set k: its mean over projects minus the mean of all sets (the reference the observed clock set is read against)
    Lk = (m_pair - m_pair.mean(1, keepdims=True)).mean(0)
    Dk = (m_diff - m_diff.mean(1, keepdims=True)).mean(0)
    pd.DataFrame({"set": range(len(Lk)), "L": Lk, "Delta": Dk}).to_csv(OUT / "null_draws.csv", index=False)
    rows = []
    for name, d, raw in [("L (H221.1)", Lk, m_pair.mean(0)), ("Delta (H221.2)", Dk, m_diff.mean(0))]:
        sd = d.std(ddof=1)
        rows.append(dict(statistic=name, n_projects=len(ok), n_sets=len(d), null_mean_random_sets=raw.mean(), null_mean_centred=d.mean(),
                         null_sd=sd, mde_log2=K * sd, scored=bool(K * sd <= MDE_MAX)))
    rows.append(dict(statistic="H221.3 (gene x project points)", n_projects=len(ok), n_sets=np.nan,
                     null_mean_random_sets=np.nan, null_mean_centred=np.nan, null_sd=np.nan, mde_log2=np.nan,
                     scored=True, n_points=len(GENES) * len(ok)))
    pw = pd.DataFrame(rows); pw.to_csv(OUT / "power.csv", index=False); say("power\n" + pw.to_string())


# ---- tests (registered; not run before the power note) ------------------------------------------------------------------------------
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); adj = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o])
    out = np.empty_like(p); out[o] = np.minimum(adj, 1); return out


def welch_down(a, b):
    """Share of genes with mean(a) < mean(b) and Welch p < 0.05 (columns = genes)."""
    t, pv = ttest_ind(a, b, axis=0, equal_var=False)
    return (t < 0) & (pv < 0.05)


def pooled_spearman(x, y): return spearmanr(x, y)[0]


def h2213(diff, offset, proj_id, nperm, rng):
    """Spearman(design difference, daytime offset) over gene x project points; offsets permuted over genes within project."""
    obs = pooled_spearman(diff, offset); idx = [np.flatnonzero(proj_id == u) for u in np.unique(proj_id)]
    null = np.empty(nperm)
    for i in range(nperm):
        o = offset.copy()
        for j in idx:
            o[j] = offset[rng.permutation(j)]
        null[i] = pooled_spearman(diff, o)
    return obs, (1 + (null >= obs).sum()) / (1 + nperm)


def stage_tests():
    el = pd.read_csv(OUT / "eligibility.csv"); ok = el.project[el.eligible].tolist()
    tod = tod_by_subject(); rng = np.random.default_rng(SEED)
    sets = pd.read_csv(OUT / "random_sets.csv")
    pur = pd.read_csv(DATA / "tcga_purity_absolute" / "mastercalls.txt", sep="\t")
    pur = pur.set_index("array").purity
    cm, cd, nm, nd, gene_rows, hardy_rows, ie_rows, share_rows, pur_rows = [], [], [], [], [], [], [], [], []
    h3_diff, h3_off, h3_proj = [], [], []
    for p in ok:
        P = Proj(p); gi = {g: i for i, g in enumerate(P.genes)}
        S = np.array([[gi[g] for g in d.sort_values("rank").gene] for _, d in sets[sets.project == p].groupby("set")])
        lp, lg = P.lfc_paired(), P.lfc_gtex(); diff = lg - lp
        cm.append(lp[P.clock].mean() - lp[S].mean()); cd.append(diff[P.clock].mean())
        nm.append(lp[S].mean(1)); nd.append(diff[S].mean(1))
        tp, np_ = P.tum[P.pt], P.nor[P.pn]
        dp, dg = welch_down(tp, np_), welch_down(P.tum, P.gtex)
        for g, i in zip(GENES, P.clock):
            gene_rows.append(dict(project=p, gene=g, lfc_paired=lp[i], lfc_gtex=lg[i], design_diff=diff[i],
                                  down_paired=bool(dp[i]), down_gtex=bool(dg[i])))
        share_rows.append(dict(project=p, clock_down_paired=dp[P.clock].mean(), clock_down_gtex=dg[P.clock].mean(),
                               random_down_paired=dp[S].mean(), random_down_gtex=dg[S].mean()))
        ie = [gi[g] for g in IEG if g in gi]
        ie_rows.append(dict(project=p, n_genes=len(ie), lfc_paired=lp[ie].mean(), lfc_gtex=lg[ie].mean()))
        tt = pd.Series(tod.reindex(P.subject).to_numpy()); has = tt.notna().to_numpy(); day = is_day(tt).to_numpy()
        off = P.gtex[day][:, P.clock].mean(0) - P.gtex[has][:, P.clock].mean(0)
        h3_diff += list(diff[P.clock]); h3_off += list(off); h3_proj += [p] * len(GENES)
        for h in range(5):  # GTEx by Hardy class (aggregates; classes with >= 10 samples)
            m = P.hardy == h
            if m.sum() >= 10:
                hardy_rows.append(dict(project=p, hardy=h, n=int(m.sum()), clock_lfc_gtex=P.lfc_gtex(m)[P.clock].mean(),
                                       random_lfc_gtex=P.lfc_gtex(m)[S].mean()))
        pr = pur.reindex([s[:15] for s in P.tum_sample]).to_numpy()
        lev = P.tum[:, P.clock].mean(1); k = ~np.isnan(pr)
        pur_rows.append(dict(project=p, n=int(k.sum()), spearman_clock_level_vs_purity=spearmanr(lev[k], pr[k])[0]))
    cm, cd = np.array(cm), np.array(cd); nm, nd = np.array(nm), np.array(nd)
    Lk = (nm - nm.mean(1, keepdims=True)).mean(0); Dk = (nd - nd.mean(1, keepdims=True)).mean(0)
    L, D = cm.mean(), cd.mean()  # Delta is read against the raw random-set values (two-sided around their mean)
    Draw = nd.mean(0)
    pL = (1 + (Lk <= L).sum()) / (1 + len(Lk)); pD = (1 + (np.abs(Draw - Draw.mean()) >= abs(D - Draw.mean())).sum()) / (1 + len(Draw))
    r3, p3 = h2213(np.array(h3_diff), np.array(h3_off), np.array(h3_proj), args.nperm, rng)
    res = pd.DataFrame([dict(hypothesis="H221.1", statistic="L", value=L, p=pL, sided="lower"),
                        dict(hypothesis="H221.2", statistic="Delta", value=D, p=pD, sided="two"),
                        dict(hypothesis="H221.3", statistic="spearman", value=r3, p=p3, sided="greater")])
    res["p_holm"] = holm(res.p); res.to_csv(OUT / "results.csv", index=False); say("results\n" + res.to_string())
    pd.DataFrame(gene_rows).to_csv(OUT / "per_gene_project.csv", index=False)
    pd.DataFrame(share_rows).to_csv(OUT / "share_down.csv", index=False)
    pd.DataFrame(ie_rows).to_csv(OUT / "immediate_early.csv", index=False)
    pd.DataFrame(hardy_rows).to_csv(OUT / "gtex_by_hardy.csv", index=False)
    pd.DataFrame(pur_rows).to_csv(OUT / "purity.csv", index=False)


if __name__ == "__main__":
    {"power": stage_power, "tests": stage_tests}[args.stage]()
