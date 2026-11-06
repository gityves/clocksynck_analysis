"""A210: does MYC amplification lower tumour clock coherence at matched proliferation? (registered in
paper/docs/registered_experiments.md, A210).
  --stage power : MYC copy-number extraction (python -I parse of the downloaded files), cell eligibility, permutation nulls of T,
                  MDE = 2.49 x null SD (eligibility.csv, counts.csv, power.csv, null_draws.csv). No coherence of any real group is computed.
  --stage tests : observed T, one-sided permutation p, Holm over the scored hypotheses, pseudo-clock companion (200 sets),
                  descriptive rows (per-cell T, proliferation imbalance, gain vs neutral, MYC expression top vs bottom tertile).
    python -I compute_a210.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 4]
Cells: METABRIC = PAM50 x proliferation tertile (tertile within PAM50); TCGA = project x proliferation tertile (within project). Tertiles
are taken over all tumours of the stratum with a copy-number call (any GISTIC state), before gains / losses are set aside. A cell needs
>= 10 amplified (GISTIC 2) and >= 10 neutral (0) tumours; T = mean over cells of coherence(amplified) - coherence(neutral); the MYC labels
are permuted within cell. Coherence = A87 `decompose.components` of the 14 clock genes, group centred (a166 / a168 `coh`).
Expression, PAM50, G2M score and pseudo-clock sets as A166 (METABRIC; the pseudo sets are those of A166, pseudo_sets_metabric.csv) and as
A172 (TCGA; recount3 `tcga_paired` loader, counts-based Y, A168 pseudo construction, seeds from crc32('A210'))."""
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
from scipy.stats import norm  # noqa: E402

import a168_lib as L  # noqa: E402  (also exposes a166_lib as L.A166)
import clocksync as cs  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=4)
args = ap.parse_args()
OUT = HERE.parent / "output"; OUT.mkdir(exist_ok=True)
CACHE = OUT / "tcga_cache"; CACHE.mkdir(exist_ok=True)
K, MIN_ARM, MDE_MAX = norm.ppf(0.95) + norm.ppf(0.80), 10, 0.25
SEED = zlib.crc32(b"A210") % (2 ** 31)
NCHUNK = 20  # null draws are generated in fixed chunks, so the result does not depend on --jobs
PROJ_A97 = ["BLCA", "BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "OV", "PRAD", "THCA", "UCEC"]
PROJ_A172 = ["CESC", "ESCA", "GBM", "HNSC", "LGG", "LUSC", "PAAD", "PCPG", "READ", "SARC", "STAD", "TGCT"]
PROJECTS = PROJ_A97 + PROJ_A172
GENES, tertile, LN2 = list(L.GENES), L.A166.tertile, L.LN2
STATE = {-2: "loss", -1: "loss", 0: "neutral", 1: "gain", 2: "amplified"}
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)
XENA = DATA / "tcga" / "xena_gistic" / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
MB = DATA / "metabric" / "brca_metabric" / "data_cna.txt"


# ---- MYC copy-number rows (text parse only) ---------------------------------------------------------------------------------------
def myc_row(path, opener, skip):
    """The MYC row of a tab-separated CNA matrix (genes x samples; `skip` leading annotation columns) as a Series over samples."""
    with opener(path, "rt") as f:
        r = csv.reader(f, delimiter="\t"); head = next(r)
        for row in r:
            if row[0] == "MYC":
                return pd.Series(pd.to_numeric(row[skip:], errors="coerce").tolist(), index=head[skip:])
    raise KeyError("MYC")


def myc_cna():
    """(METABRIC patient -> GISTIC state, TCGA 15-character sample -> GISTIC state); derived tables cached in output/."""
    fm, ft = OUT / "myc_cna_metabric.csv", OUT / "myc_cna_tcga.csv"
    if not fm.exists():
        s = myc_row(MB, open, 2); s.index.name = "patient"; s.rename("myc_cna").to_csv(fm)  # columns: Hugo, Entrez, samples
    if not ft.exists():
        s = myc_row(XENA, gzip.open, 1); s.index.name = "sample"; s.rename("myc_cna").to_csv(ft)
    return (pd.read_csv(fm, index_col=0).myc_cna, pd.read_csv(ft, index_col=0).myc_cna)


# ---- cohort tables ----------------------------------------------------------------------------------------------------------------
def metabric(cna):
    m, X = L.A166.load("metabric"); m = m[m.pam50.notna()].copy(); X = X.loc[m.index]
    m["state"] = cna.reindex(m.patient).to_numpy(); m = m[m.state.notna()]; X = X.loc[m.index]
    m = m.reset_index(drop=True); m["stratum"] = m.pam50.astype(str); m["myc_expr"] = X["MYC"].to_numpy(float)
    return m, X


def primaries(B, p):
    s = B.samples_of(p)
    return s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")


def zmean(lc, genes):
    g = [x for x in genes if x in lc.columns]
    z = ((lc[g] - lc[g].mean()) / lc[g].std(ddof=0).replace(0, np.nan)).mean(axis=1)
    return ((z - z.mean()) / z.std()).to_numpy()


def load_tcga(B, p, need_lc=False):
    """A172 `load`: (14 clock counts + depth, log2 CPM of genes with mean CPM >= 1 (None if not need_lc), G2M proliferation z)."""
    ct = B.counts_of(p, primaries(B, p)); depth = ct.sum(axis=1)
    lc = np.log2(ct.div(depth, axis=0) * 1e6 + 1); lc = lc.loc[:, (2 ** lc - 1).mean() >= 1]
    return ct[GENES].assign(depth=depth.to_numpy()), (lc.astype(np.float32) if need_lc else None), zmean(lc, L.sh.hallmark_g2m()), \
        lc["MYC"].to_numpy(float)


def tcga_table(B, cna, p):
    """Cached per project: patient, 15-character barcode, GISTIC state, proliferation, MYC expression, depth, 14 clock counts."""
    f = CACHE / f"{p}.parquet"
    if not f.exists():
        t = primaries(B, p); clock, _, prolif, myc = load_tcga(B, p)
        d = clock.copy(); d.insert(0, "patient", t.patient.to_numpy()); d.insert(1, "sample", t.barcode.str[:15].to_numpy())
        d["prolif"], d["myc_expr"] = prolif, myc
        d.to_parquet(f)
    d = pd.read_parquet(f); d["state"] = cna.reindex(d["sample"]).to_numpy()
    return d[d.state.notna()].reset_index(drop=True)


# ---- cells and statistic ----------------------------------------------------------------------------------------------------------
def make_cells(m):
    """All stratum x proliferation-tertile cells of table m (stratum, prolif, state, row positions 0..n-1) with their counts; eligible
    cells carry the row positions (amplified first) and labels (1 = amplified)."""
    m = m.copy(); m["t"] = m.groupby("stratum").prolif.transform(tertile); rows, cells = [], []
    for (s, t), g in m.groupby(["stratum", "t"]):
        st, pos = g.state.to_numpy(), g.index.to_numpy()
        n = {k: int((st == v).sum()) for k, v in (("n_amplified", 2), ("n_gain", 1), ("n_neutral", 0))}; n["n_loss"] = int((st < 0).sum())
        ok = min(n["n_amplified"], n["n_neutral"]) >= MIN_ARM
        rows.append(dict(stratum=s, tertile=int(t), **n, eligible=ok))
        if ok:
            a, z = pos[st == 2], pos[st == 0]
            cells.append(dict(key=f"{s}|{int(t)}", stratum=s, tertile=int(t), idx=np.r_[a, z], lab=np.r_[np.ones(len(a), int), np.zeros(len(z), int)]))
    return pd.DataFrame(rows), cells


def diffs(Y, wf, cells, perms=None):
    """coherence(amplified) - coherence(neutral) per cell of one block (Y rows = block positions); perms: one label permutation per cell."""
    out = []
    for k, c in enumerate(cells):
        l = c["lab"] if perms is None else c["lab"][perms[k]]
        out.append(L.coh(Y, c["idx"][l == 1], wf) - L.coh(Y, c["idx"][l == 0], wf))
    return np.array(out)


def null_chunk(blocks, seed, n):
    """n permutation draws of T (mean over all cells of all blocks); blocks = [(Y, wf, cells)]."""
    rng = np.random.default_rng(seed); out = np.empty(n)
    for i in range(n):
        out[i] = np.mean(np.concatenate([diffs(Y, wf, cells, [rng.permutation(len(c["lab"])) for c in cells]) for Y, wf, cells in blocks]))
    return out


def null_draws(blocks, tag, n):
    sizes = [n // NCHUNK + (j < n % NCHUNK) for j in range(NCHUNK)]
    res = Parallel(n_jobs=args.jobs)(delayed(null_chunk)(blocks, np.random.SeedSequence([SEED, tag, j]), s) for j, s in enumerate(sizes))
    return np.concatenate(res)


# ---- designs ----------------------------------------------------------------------------------------------------------------------
cna_mb, cna_tc = myc_cna()
say("MYC rows: METABRIC", cna_mb.notna().sum(), "TCGA", cna_tc.notna().sum())
mb, Xmb = metabric(cna_mb)
el_mb, cells_mb = make_cells(mb); el_mb.insert(0, "cohort", "METABRIC")
B = data._builder("tcga_paired")
tc, el_tc, cells_tc = {}, [], {}
for j, p in enumerate(PROJECTS):
    tc[p] = tcga_table(B, cna_tc, p)
    e, c = make_cells(tc[p].assign(stratum=p)); e.insert(0, "cohort", "TCGA"); el_tc.append(e); cells_tc[p] = c
    say(p, "tumours with a call", len(tc[p]), "eligible cells", len(c))
el = pd.concat([el_mb] + el_tc, ignore_index=True); el.to_csv(OUT / "eligibility.csv", index=False)
cnt = pd.concat([mb.groupby("stratum").state.apply(lambda s: pd.Series(dict(n=len(s), n_amplified=(s == 2).sum(), n_gain=(s == 1).sum(),
                                                                          n_neutral=(s == 0).sum(), n_loss=(s < 0).sum()))).unstack().assign(cohort="METABRIC")]
                + [pd.DataFrame([dict(stratum=p, cohort="TCGA", n=len(t), n_amplified=(t.state == 2).sum(), n_gain=(t.state == 1).sum(),
                                      n_neutral=(t.state == 0).sum(), n_loss=(t.state < 0).sum())]).set_index("stratum") for p, t in tc.items()])
cnt.reset_index().to_csv(OUT / "counts.csv", index=False)
summ = lambda c: dict(n_cells=len(c), n_tumours=int(sum(len(x["lab"]) for x in c)), n_amplified=int(sum(x["lab"].sum() for x in c)))
Ymb = LN2 * Xmb[GENES].to_numpy(float)
Yof = lambda counts, genes: np.log((counts[list(genes)].to_numpy(float) + 0.5) / counts.depth.to_numpy(float)[:, None])
blocks = {"H210.1": [(Ymb, L.WF, cells_mb)],
          "H210.2": [(Yof(tc[p], GENES), L.WF, cells_tc[p]) for p in PROJECTS if cells_tc[p]]}
cells_all = {"H210.1": cells_mb, "H210.2": [c for p in PROJECTS for c in cells_tc[p]]}
say({h: summ(c) for h, c in cells_all.items()})

if args.stage == "power":
    nl = {h: null_draws(blocks[h], i, args.nperm) for i, h in enumerate(blocks) if cells_all[h]}
    pd.DataFrame(nl).to_csv(OUT / "null_draws.csv", index=False)
    pw = pd.DataFrame([dict(hypothesis=h, cohort="METABRIC" if h == "H210.1" else "TCGA", **summ(cells_all[h]), n_perm=args.nperm,
                            null_mean=nl[h].mean() if h in nl else np.nan, null_sd=nl[h].std(ddof=1) if h in nl else np.nan,
                            MDE_80pct=K * nl[h].std(ddof=1) if h in nl else np.nan) for h in blocks])
    pw["scored"] = pw.MDE_80pct <= MDE_MAX
    pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests (not run before the dated power note) ----------------------------------------------------------------------------------
pw = pd.read_csv(OUT / "power.csv").set_index("hypothesis"); nl = pd.read_csv(OUT / "null_draws.csv"); assert len(nl) == args.nperm
obs_cells, obs = {}, {}
for h in blocks:
    d = [diffs(Y, wf, c) for Y, wf, c in blocks[h]]; obs_cells[h] = np.concatenate(d); obs[h] = float(obs_cells[h].mean())

# pseudo-clock companion: per block, k-th set -> per-cell differences; T_k = mean over all cells
rng = np.random.default_rng(1660)
sets_mb = pd.read_csv(L.A166.OUT / "pseudo_sets_metabric.csv").to_numpy()[: args.npseudo].tolist()
pseudo = {"H210.1": np.array([np.mean(diffs(LN2 * Xmb[s].to_numpy(float), L.wf_for(s), cells_mb)) for s in sets_mb])}
pc = []
for j, p in enumerate(PROJECTS):
    if not cells_tc[p]:
        continue
    _, lc, _, _ = load_tcga(B, p, need_lc=True)
    t = primaries(B, p); keep = tc[p]["patient"].to_numpy(); lc = lc.set_axis(t.patient.to_numpy()).loc[keep]
    sets = L.pseudo_sets(lc, L.A166.pool_genes(lc), GENES, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 1000 + j])))
    dep = tc[p].depth.to_numpy(float)[:, None]
    pc.append(np.array([diffs(np.log((2 ** lc[list(g)].to_numpy(float) - 1) / 1e6 + 0.5 / dep), L.wf_pseudo(L.WF, g), cells_tc[p]) for g in sets]))
    say(p, "pseudo done")
pseudo["H210.2"] = np.concatenate(pc, axis=1).mean(axis=1)
pd.DataFrame(pseudo).to_csv(OUT / "pseudo_T.csv", index=False)

rows = []
for h in blocks:
    sc = bool(pw.loc[h, "scored"])
    r = dict(hypothesis=h, scored=sc, T=obs[h], null_sd=pw.loc[h, "null_sd"], MDE_80pct=pw.loc[h, "MDE_80pct"],
             p_one_lower=pl(nl[h].to_numpy(), obs[h]) if sc else np.nan, pseudo_percentile=100 * np.mean(pseudo[h] <= obs[h]),
             pseudo_q05=np.quantile(pseudo[h], 0.05), T_below_pseudo_q05=bool(obs[h] < np.quantile(pseudo[h], 0.05)))
    rows.append(r)
res = pd.DataFrame(rows); s = res.scored.to_numpy(); ps = res.p_one_lower[s].to_numpy(); o = np.argsort(ps)
res["p_holm"] = np.nan
res.loc[s, "p_holm"] = np.minimum(1, np.maximum.accumulate(np.sort(ps) * (len(ps) - np.arange(len(ps)))))[np.argsort(o)]
res["verdict"] = np.where(~res.scored, "not estimable (MDE > 0.25)", np.where(res.p_holm < 0.05, "pass", "fail"))
res["specificity"] = np.where(res.verdict == "pass", np.where(res.T_below_pseudo_q05, "clock-specific", "not specific to the clock"), "not applicable")
res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True)

# descriptive: per cell T, residual proliferation imbalance, gain vs neutral, MYC expression top vs bottom tertile
desc = []
for h, tabs in (("H210.1", [("METABRIC", mb, Ymb, L.WF, cells_mb)]),
                ("H210.2", [(p, tc[p].assign(stratum=p), Yof(tc[p], GENES), L.WF, cells_tc[p]) for p in PROJECTS])):
    for coh_name, m, Y, wf, cells in tabs:
        m["t"] = m.groupby("stratum").prolif.transform(tertile)
        for (sname, t), g in m.groupby(["stratum", "t"]):
            st, pos = g.state.to_numpy(), g.index.to_numpy(); r = dict(hypothesis=h, cell=f"{sname}|{int(t)}")
            ok = min((st == 2).sum(), (st == 0).sum()) >= MIN_ARM
            r.update(eligible=ok, n_amplified=int((st == 2).sum()), n_neutral=int((st == 0).sum()), n_gain=int((st == 1).sum()))
            if ok:
                r["T_cell"] = L.coh(Y, pos[st == 2], wf) - L.coh(Y, pos[st == 0], wf)
                r["prolif_amplified_minus_neutral"] = g.prolif.to_numpy()[st == 2].mean() - g.prolif.to_numpy()[st == 0].mean()
            if min((st == 1).sum(), (st == 0).sum()) >= MIN_ARM:
                r["T_gain_minus_neutral"] = L.coh(Y, pos[st == 1], wf) - L.coh(Y, pos[st == 0], wf)
            if len(g) >= 3 * MIN_ARM:
                e = tertile(g.myc_expr.to_numpy()); r["T_myc_expr_top_minus_bottom"] = L.coh(Y, pos[e == 2], wf) - L.coh(Y, pos[e == 0], wf)
            desc.append(r)
desc = pd.DataFrame(desc); desc.to_csv(OUT / "descriptive_cells.csv", index=False)
summ_d = desc.groupby("hypothesis")[["T_cell", "prolif_amplified_minus_neutral", "T_gain_minus_neutral", "T_myc_expr_top_minus_bottom"]].mean()
summ_d.to_csv(OUT / "descriptive.csv"); print(summ_d.round(4).to_string(), flush=True)
say("done tests")
