"""A213: is tumour clock coherence lower in heterogeneous / copy-number-burdened tumours at matched proliferation? (registered in
paper/docs/registered_experiments.md, A213; mirrors A210).
  --stage power : measures (ABSOLUTE subclonal genome fraction, FGA from GISTIC), cell eligibility, permutation nulls of T,
                  MDE = 2.49 x null SD (eligibility.csv, counts.csv, power.csv, null_draws.csv). No coherence of any real group is computed.
  --stage tests : observed T, one-sided permutation p, Holm over the scored H213.1-3, pseudo-clock companion (200 sets), S213.4
                  (FGA-matched rerun of A210 H210.2) and descriptives.
    python -I compute_a213.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 4]
Cells: stratum (TCGA project; METABRIC PAM50) x proliferation tertile (tertile within stratum over tumours with the measure); within a cell,
top vs bottom tertile of the measure; >= 10 per arm. T = mean over cells of coherence(top) - coherence(bottom); labels permuted within cell.
H213.1 TCGA subclonal fraction; H213.2 TCGA FGA; H213.3 METABRIC FGA. Coherence, expression, proliferation, pseudo sets as A210."""
import argparse
import csv
import gzip
import os
import sys
import zlib
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts"))
sys.path.insert(0, str(ROOT / "paper"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm, spearmanr  # noqa: E402

import a168_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.decompose import components  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=4)
args = ap.parse_args()
OUT = HERE.parent / "output"; OUT.mkdir(exist_ok=True)
A210 = HERE.parents[1] / "A210_myc_coherence" / "output"
CACHE = A210 / "tcga_cache"  # per-project clock counts, proliferation, MYC expression (A210; open-access TCGA)
K, MIN_ARM, MDE_MAX = norm.ppf(0.95) + norm.ppf(0.80), 10, 0.25
SEED = zlib.crc32(b"A213") % (2 ** 31)
NCHUNK = 20
PROJ_A97 = ["BLCA", "BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "OV", "PRAD", "THCA", "UCEC"]
PROJ_A172 = ["CESC", "ESCA", "GBM", "HNSC", "LGG", "LUSC", "PAAD", "PCPG", "READ", "SARC", "STAD", "TGCT"]
PROJECTS = PROJ_A97 + PROJ_A172
GENES, tertile, LN2 = list(L.GENES), L.A166.tertile, L.LN2
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)
XENA = DATA / "tcga" / "xena_gistic" / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
MB = DATA / "metabric" / "brca_metabric" / "data_cna.txt"
ABS = DATA / "tcga" / "absolute_pancanatlas" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"


# ---- measures (text parse only) ---------------------------------------------------------------------------------------------------
def fga(path, opener, skip):
    """Per-sample fraction of genes with GISTIC != 0, over genes non-missing in every sample of the matrix (`skip` annotation columns)."""
    with opener(path, "rt") as f:
        r = csv.reader(f, delimiter="\t"); head = next(r)
        M = np.array([pd.to_numeric(row[skip:], errors="coerce") for row in r], dtype=np.float32)
    ok = ~np.isnan(M).any(axis=1)
    say(path.name, "genes", M.shape[0], "complete in all samples", int(ok.sum()))
    return pd.Series((M[ok] != 0).mean(axis=0), index=head[skip:])


def cached(name, fn, idx):
    f = OUT / name
    if not f.exists():
        s = fn(); s.index.name = idx; s.rename("fga").to_csv(f)
    return pd.read_csv(f, index_col=0).fga


def measures():
    """(METABRIC patient -> FGA, TCGA 15-char sample -> FGA, ABSOLUTE table by 15-char barcode: subclonal fraction, purity, ploidy)."""
    fm = cached("fga_metabric.csv", lambda: fga(MB, open, 2), "patient")
    ft = cached("fga_tcga.csv", lambda: fga(XENA, gzip.open, 1), "sample")
    a = pd.read_csv(ABS, sep="\t"); a = a[a["call status"] == "called"].drop_duplicates("array").set_index("array")
    a = a[["Subclonal genome fraction", "purity", "ploidy"]].rename(columns={"Subclonal genome fraction": "subclonal"})
    return fm, ft, a


# ---- cohort tables ----------------------------------------------------------------------------------------------------------------
def primaries(B, p):
    s = B.samples_of(p)
    return s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")


def zmean(lc, genes):
    g = [x for x in genes if x in lc.columns]
    z = ((lc[g] - lc[g].mean()) / lc[g].std(ddof=0).replace(0, np.nan)).mean(axis=1)
    return ((z - z.mean()) / z.std()).to_numpy()


def load_lc(B, p):
    """log2 CPM of genes with mean CPM >= 1 for the primaries of p (A172 `load`, rows as the A210 cache)."""
    ct = B.counts_of(p, primaries(B, p)); depth = ct.sum(axis=1)
    lc = np.log2(ct.div(depth, axis=0) * 1e6 + 1)
    return lc.loc[:, (2 ** lc - 1).mean() >= 1].astype(np.float32)


def tcga_table(p, ft, ab, cna_myc):
    """A210 cache of project p (patient, 15-char sample, prolif, depth, 14 clock counts) + FGA, ABSOLUTE columns, MYC GISTIC state."""
    d = pd.read_parquet(CACHE / f"{p}.parquet").reset_index(drop=True)
    d["fga"] = ft.reindex(d["sample"]).to_numpy()
    for c in ("subclonal", "purity", "ploidy"):
        d[c] = ab[c].reindex(d["sample"]).to_numpy()
    d["state"] = cna_myc.reindex(d["sample"]).to_numpy()
    return d


def metabric(fm):
    m, X = L.A166.load("metabric"); m = m[m.pam50.notna()].copy(); X = X.loc[m.index]
    m["fga"] = fm.reindex(m.patient).to_numpy(); m = m[m.fga.notna()]; X = X.loc[m.index]
    m = m.reset_index(drop=True); m["stratum"] = m.pam50.astype(str)
    return m, X


# ---- cells and statistic ----------------------------------------------------------------------------------------------------------
def make_cells(m, measure):
    """Stratum x proliferation-tertile cells of the tumours of m with `measure` (positions of m kept); within a cell the measure's top
    (label 1) vs bottom (label 0) tertile, >= MIN_ARM per arm. Returns (eligibility rows, eligible cells)."""
    m = m[m[measure].notna()].copy(); m["t"] = m.groupby("stratum").prolif.transform(tertile); rows, cells = [], []
    for (s, t), g in m.groupby(["stratum", "t"]):
        pos = g.index.to_numpy(); e = tertile(g[measure].to_numpy(float)) if len(g) >= 3 else np.full(len(g), 1)
        top, bot = pos[e == 2], pos[e == 0]
        ok = min(len(top), len(bot)) >= MIN_ARM
        rows.append(dict(stratum=s, tertile=int(t), n=len(g), n_top=len(top), n_bottom=len(bot), eligible=ok))
        if ok:
            cells.append(dict(key=f"{s}|{int(t)}", stratum=s, tertile=int(t), idx=np.r_[top, bot],
                              lab=np.r_[np.ones(len(top), int), np.zeros(len(bot), int)]))
    return pd.DataFrame(rows), cells


def diffs(Y, wf, cells, perms=None):
    out = []
    for k, c in enumerate(cells):
        l = c["lab"] if perms is None else c["lab"][perms[k]]
        out.append(L.coh(Y, c["idx"][l == 1], wf) - L.coh(Y, c["idx"][l == 0], wf))
    return np.array(out)


def null_chunk(blocks, seed, n):
    rng = np.random.default_rng(seed); out = np.empty(n)
    for i in range(n):
        out[i] = np.mean(np.concatenate([diffs(Y, wf, cells, [rng.permutation(len(c["lab"])) for c in cells]) for Y, wf, cells in blocks]))
    return out


def swap_chunk(blocks, seed, n):
    """Null of S213.4: the labels of each matched pair swapped with probability 1/2. Cells: dict(a, z) = paired positions."""
    rng = np.random.default_rng(seed); out = np.empty(n)
    for i in range(n):
        d = []
        for Y, wf, cells in blocks:
            for c in cells:
                f = rng.random(len(c["a"])) < 0.5
                d.append(L.coh(Y, np.where(f, c["z"], c["a"]), wf) - L.coh(Y, np.where(f, c["a"], c["z"]), wf))
        out[i] = np.mean(d)
    return out


def draws(fn, blocks, tag, n):
    sizes = [n // NCHUNK + (j < n % NCHUNK) for j in range(NCHUNK)]
    return np.concatenate(Parallel(n_jobs=args.jobs)(delayed(fn)(blocks, np.random.SeedSequence([SEED, tag, j]), s) for j, s in enumerate(sizes)))


# ---- designs ----------------------------------------------------------------------------------------------------------------------
fm, ft, ab = measures()
cna_myc = pd.read_csv(A210 / "myc_cna_tcga.csv", index_col=0).myc_cna
mb, Xmb = metabric(fm)
tc = {p: tcga_table(p, ft, ab, cna_myc) for p in PROJECTS}
HYP = {"H213.1": ("TCGA", "subclonal"), "H213.2": ("TCGA", "fga"), "H213.3": ("METABRIC", "fga")}
el, cells_h = [], {}
for h, (coh_, meas) in HYP.items():
    if coh_ == "METABRIC":
        e, c = make_cells(mb, meas); e.insert(0, "hypothesis", h); el.append(e); cells_h[h] = {"METABRIC": c}
    else:
        cells_h[h] = {}
        for p in PROJECTS:
            e, c = make_cells(tc[p].assign(stratum=p), meas); e.insert(0, "hypothesis", h); el.append(e); cells_h[h][p] = c
            say(h, p, "tumours with the measure", int(tc[p][meas].notna().sum()), "eligible cells", len(c))
pd.concat(el, ignore_index=True).to_csv(OUT / "eligibility.csv", index=False)
cnt = [dict(hypothesis=h, stratum=k, n_with_measure=int((mb if k == "METABRIC" else tc[k])[HYP[h][1]].notna().sum()), n_cells=len(c),
            n_tumours=int(sum(len(x["lab"]) for x in c))) for h in HYP for k, c in cells_h[h].items()]
pd.DataFrame(cnt).to_csv(OUT / "counts.csv", index=False)
Ymb = LN2 * Xmb[GENES].to_numpy(float)
Yof = lambda counts, genes: np.log((counts[list(genes)].to_numpy(float) + 0.5) / counts.depth.to_numpy(float)[:, None])
Ytc = {p: Yof(tc[p], GENES) for p in PROJECTS}
blocks = {h: [(Ymb, L.WF, cells_h[h]["METABRIC"])] if HYP[h][0] == "METABRIC" else [(Ytc[p], L.WF, cells_h[h][p]) for p in PROJECTS if cells_h[h][p]]
          for h in HYP}
cells_all = {h: [c for cs_ in cells_h[h].values() for c in cs_] for h in HYP}
summ = lambda c: dict(n_cells=len(c), n_tumours=int(sum(len(x["lab"]) for x in c)), n_top=int(sum(x["lab"].sum() for x in c)))
say({h: summ(c) for h, c in cells_all.items()})

if args.stage == "power":
    nl = {h: draws(null_chunk, blocks[h], i, args.nperm) for i, h in enumerate(HYP) if cells_all[h]}
    pd.DataFrame(nl).to_csv(OUT / "null_draws.csv", index=False)
    pw = pd.DataFrame([dict(hypothesis=h, cohort=HYP[h][0], measure=HYP[h][1], **summ(cells_all[h]), n_perm=args.nperm,
                            n_projects=sum(bool(c) for c in cells_h[h].values()),
                            null_mean=nl[h].mean() if h in nl else np.nan, null_sd=nl[h].std(ddof=1) if h in nl else np.nan,
                            MDE_80pct=K * nl[h].std(ddof=1) if h in nl else np.nan) for h in HYP])
    pw["scored"] = pw.MDE_80pct <= MDE_MAX
    pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests (not run before the dated power note) ----------------------------------------------------------------------------------
pw = pd.read_csv(OUT / "power.csv").set_index("hypothesis"); nl = pd.read_csv(OUT / "null_draws.csv"); assert len(nl) == args.nperm
obs_cells = {h: np.concatenate([diffs(Y, wf, c) for Y, wf, c in blocks[h]]) for h in HYP}
obs = {h: float(obs_cells[h].mean()) for h in HYP}

# pseudo-clock companion (200 sets per project, one set family per project shared by H213.1 and H213.2; A210 construction)
B = data._builder("tcga_paired")
sets_mb = pd.read_csv(L.A166.OUT / "pseudo_sets_metabric.csv").to_numpy()[: args.npseudo].tolist()
pseudo = {"H213.3": np.array([np.mean(diffs(LN2 * Xmb[s].to_numpy(float), L.A166.wf_for(s), cells_h["H213.3"]["METABRIC"])) for s in sets_mb])}
pc = {"H213.1": [], "H213.2": []}
for j, p in enumerate(PROJECTS):
    if not (cells_h["H213.1"][p] or cells_h["H213.2"][p]):
        continue
    lc = load_lc(B, p); t = primaries(B, p)
    assert (t.patient.to_numpy() == tc[p].patient.to_numpy()).all()
    lc = lc.set_axis(t.patient.to_numpy())
    sets = L.pseudo_sets(lc, L.A166.pool_genes(lc), GENES, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 1000 + j])))
    dep = tc[p].depth.to_numpy(float)[:, None]
    for h in pc:
        if cells_h[h][p]:
            pc[h].append(np.array([diffs(np.log((2 ** lc[list(g)].to_numpy(float) - 1) / 1e6 + 0.5 / dep), L.wf_pseudo(L.WF, g), cells_h[h][p])
                                   for g in sets]))
    say(p, "pseudo done")
for h in pc:
    pseudo[h] = np.concatenate(pc[h], axis=1).mean(axis=1)
pd.DataFrame(pseudo).to_csv(OUT / "pseudo_T.csv", index=False)

rows = []
for h in HYP:
    sc = bool(pw.loc[h, "scored"])
    rows.append(dict(hypothesis=h, cohort=HYP[h][0], measure=HYP[h][1], scored=sc, T=obs[h], null_mean=pw.loc[h, "null_mean"], null_sd=pw.loc[h, "null_sd"],
                     MDE_80pct=pw.loc[h, "MDE_80pct"], p_one_lower=pl(nl[h].to_numpy(), obs[h]) if sc else np.nan,
                     pseudo_percentile=100 * np.mean(pseudo[h] <= obs[h]), pseudo_q05=np.quantile(pseudo[h], 0.05),
                     T_below_pseudo_q05=bool(obs[h] < np.quantile(pseudo[h], 0.05))))
res = pd.DataFrame(rows); s = res.scored.to_numpy(); ps = res.p_one_lower[s].to_numpy(); o = np.argsort(ps)
res["p_holm"] = np.nan
res.loc[s, "p_holm"] = np.minimum(1, np.maximum.accumulate(np.sort(ps) * (len(ps) - np.arange(len(ps)))))[np.argsort(o)]
res["verdict"] = np.where(~res.scored, "not estimable (MDE > 0.25)", np.where(res.p_holm < 0.05, "pass", "fail"))
res["specificity"] = np.where(res.verdict == "pass", np.where(res.T_below_pseudo_q05, "clock-specific", "not specific to the clock"), "not applicable")
res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True)
say("tests done")

# ---- S213.4: A210 H210.2 cells, amplified matched 1:1 to neutral on FGA -----------------------------------------------------------
def match(fa, fz, caliper=0.05):
    """1:1 nearest-neighbour matching without replacement: all (amplified, neutral) pairs within the caliper taken in increasing distance
    (ties by position), each tumour used once. Returns (positions in fa, positions in fz)."""
    D = np.abs(fa[:, None] - fz[None, :]); i, j = np.nonzero(D <= caliper)
    o = np.lexsort((j, i, D[i, j])); ua, uz, A, Z = set(), set(), [], []
    for k in o:
        if i[k] not in ua and j[k] not in uz:
            ua.add(i[k]); uz.add(j[k]); A.append(i[k]); Z.append(j[k])
    return np.array(A, int), np.array(Z, int)


mcells, mrows = [], []
for p in PROJECTS:
    m = tc[p][tc[p].state.notna()].reset_index(drop=True); m["stratum"] = p; m["t"] = m.groupby("stratum").prolif.transform(tertile)
    for t, g in m.groupby("t"):
        st = g.state.to_numpy()
        if min((st == 2).sum(), (st == 0).sum()) < MIN_ARM:
            continue
        # A210 eligible cell (positions in the project table tc[p]); matching needs FGA
        keep = np.flatnonzero(tc[p].state.notna().to_numpy()); ga = g[g.state == 2]; gz = g[g.state == 0]
        ga, gz = ga[ga.fga.notna()], gz[gz.fga.notna()]
        a, z = match(ga.fga.to_numpy(), gz.fga.to_numpy())
        mrows.append(dict(project=p, tertile=int(t), n_amp=int((st == 2).sum()), n_neutral=int((st == 0).sum()), n_amp_fga=len(ga), n_neutral_fga=len(gz),
                          n_pairs=len(a), eligible=len(a) >= MIN_ARM))
        # positions refer to the filtered-by-MYC table m, whose rows are the tc[p] rows with a MYC call (reset)
        if len(a) >= MIN_ARM:
            mcells.append(dict(project=p, a=keep[ga.index.to_numpy()[a]], z=keep[gz.index.to_numpy()[z]],
                               fa=ga.fga.to_numpy()[a], fz=gz.fga.to_numpy()[z]))
pd.DataFrame(mrows).to_csv(OUT / "s213_4_cells.csv", index=False)
blk4 = [(Ytc[p], L.WF, [c for c in mcells if c["project"] == p]) for p in PROJECTS if any(c["project"] == p for c in mcells)]
T4 = float(np.mean([L.coh(Y, c["a"], wf) - L.coh(Y, c["z"], wf) for Y, wf, cc in blk4 for c in cc]))
nl4 = draws(swap_chunk, blk4, 4, args.nperm)
a210 = pd.read_csv(A210 / "tests.csv").set_index("hypothesis").loc["H210.2"]
a210_pw = pd.read_csv(A210 / "power.csv").set_index("hypothesis").loc["H210.2"]
exc210 = float(a210["T"] - a210_pw["null_mean"])
r4 = dict(n_cells=len(mcells), n_pairs=int(sum(len(c["a"]) for c in mcells)), T_matched=T4, null_mean=nl4.mean(), null_sd=nl4.std(ddof=1),
          excess_matched=T4 - nl4.mean(), p_one_lower=pl(nl4, T4), A210_T=float(a210["T"]), A210_null_mean=float(a210_pw["null_mean"]), A210_excess=exc210,
          ratio_excess=(T4 - nl4.mean()) / exc210,
          mean_abs_fga_gap=float(np.mean(np.concatenate([np.abs(c["fa"] - c["fz"]) for c in mcells]))),
          mean_fga_amp=float(np.mean(np.concatenate([c["fa"] for c in mcells]))), mean_fga_neutral=float(np.mean(np.concatenate([c["fz"] for c in mcells]))))
r4["FGA_carries_A210"] = bool((r4["excess_matched"] > 0.5 * exc210) and (r4["p_one_lower"] >= 0.05))  # excess negative: "less than half" in magnitude
pd.DataFrame([r4]).to_csv(OUT / "s213_4.csv", index=False); print(pd.DataFrame([r4]).T.to_string(), flush=True)
say("S213.4 done")

# ---- descriptives -----------------------------------------------------------------------------------------------------------------
def parts(Y, wf, idx_top, idx_bot):
    """(dR, dlog_plane, dlog_off) of top vs bottom, each group centred (A176 Contrast.parts)."""
    c = []
    for i in (idx_top, idx_bot):
        y = Y[i]; y = y - y.mean(0); c.append(components(y, np.zeros(y.shape[1]), 1, wf))
    return [c[0]["coherence"] - c[1]["coherence"], c[0]["log_plane"] - c[1]["log_plane"], c[0]["log_off"] - c[1]["log_off"]]


dc = []
for h, (coh_, meas) in HYP.items():
    for k, cells in cells_h[h].items():
        tab = mb if k == "METABRIC" else tc[k]
        Y, wf = (Ymb, L.WF) if k == "METABRIC" else (Ytc[k], L.WF)
        for c in cells:
            top, bot = c["idx"][c["lab"] == 1], c["idx"][c["lab"] == 0]
            pr = parts(Y, wf, top, bot)
            r = dict(hypothesis=h, project=k, cell=c["key"], n_top=len(top), n_bottom=len(bot), T_cell=L.coh(Y, top, wf) - L.coh(Y, bot, wf),
                     dR_components=pr[0], dlog_plane=pr[1], dlog_off=pr[2], ds=pr[1] - pr[2],
                     prolif_top_minus_bottom=tab.prolif.to_numpy()[top].mean() - tab.prolif.to_numpy()[bot].mean(),
                     measure_top_mean=tab[meas].to_numpy()[top].mean(), measure_bottom_mean=tab[meas].to_numpy()[bot].mean())
            if k != "METABRIC":
                r["purity_top_minus_bottom"] = np.nanmean(tab.purity.to_numpy()[top]) - np.nanmean(tab.purity.to_numpy()[bot])
            dc.append(r)
dc = pd.DataFrame(dc); dc.to_csv(OUT / "descriptive_cells.csv", index=False)
num = ["T_cell", "dR_components", "dlog_plane", "dlog_off", "ds", "prolif_top_minus_bottom", "purity_top_minus_bottom"]
dsum = dc.groupby("hypothesis")[num].mean(); dsum.insert(0, "n_cells", dc.groupby("hypothesis").size())
dsum.to_csv(OUT / "descriptive.csv"); print(dsum.round(4).to_string(), flush=True)
pp = dc.groupby(["hypothesis", "project"]).agg(n_cells=("T_cell", "size"), T=("T_cell", "mean")).reset_index()
pp.to_csv(OUT / "descriptive_per_project.csv", index=False)

# Spearman: pooled over tumours and mean of within-project coefficients
allt = pd.concat([tc[p].assign(project=p) for p in PROJECTS], ignore_index=True)
sp = []
for a_, b_, name in (("subclonal", "fga", "subclonal~fga"), ("subclonal", "prolif", "subclonal~prolif"), ("fga", "prolif", "fga~prolif")):
    d = allt.dropna(subset=[a_, b_]); w = [spearmanr(g[a_], g[b_])[0] for _, g in d.groupby("project") if len(g) >= 20]
    sp.append(dict(cohort="TCGA", pair=name, n=len(d), rho_pooled=spearmanr(d[a_], d[b_])[0], rho_within_project_mean=np.mean(w), n_projects=len(w)))
d = mb.dropna(subset=["fga", "prolif"]); w = [spearmanr(g.fga, g.prolif)[0] for _, g in d.groupby("stratum") if len(g) >= 20]
sp.append(dict(cohort="METABRIC", pair="fga~prolif", n=len(d), rho_pooled=spearmanr(d.fga, d.prolif)[0], rho_within_project_mean=np.mean(w), n_projects=len(w)))
pd.DataFrame(sp).to_csv(OUT / "descriptive_spearman.csv", index=False); print(pd.DataFrame(sp).round(3).to_string(), flush=True)
say("done tests")
