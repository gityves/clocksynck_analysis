"""A217 / A218: TCGA genomic-state losses and the proliferation loss at matched molecular subtype (registered in
paper/docs/registered_experiments.md, A217 and A218; mirrors A215 / A216).
  --stage power : PanCanAtlas SUBTYPE table (from the stored cBioPortal API responses), purity proxy, project eligibility, cell eligibility and counts,
                  arm subtype composition, permutation nulls (labels only) of H217.1-3 and of P_sub (H218.1), MDE = 2.49 x null SD
                  (subtypes_tcga.csv, purity_tcga.csv, projects.csv, eligibility.csv, counts.csv, arm_composition.csv, power.csv, null_draws.csv).
                  No coherence of any real group is computed. The A218 precondition (P_proj below its own pseudo 5th percentile) and the jackknife
                  need observed values and belong to the tests stage.
  --stage tests : observed T, one-sided permutation p, Holm over H217.1-3, pseudo-clock companion (200 sets per project and amendment), the original
                  design (project x proliferation cells) on the same tumours with its excess and the ratio of the excesses, per-cell / per-project /
                  per-subtype values, BRCA separately (descriptive); A218 P_proj, P_sub, precondition, H218.1, A = P_sub / P_proj with a
                  delete-one-project jackknife interval.
  --smoke --out DIR : both stages on random noise expression (real labels and cells) in DIR, to test the code path without touching real coherence.
    python -I compute_a217.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 4]
H217.1-3: cells = project x subtype x proliferation tertile (tertile within project x subtype over tumours with the measure); arms TP53 altered vs
intact, FGA top vs bottom tertile within cell, MYC GISTIC 2 vs 0; >= 10 per arm; T = mean over cells of coherence(arm 1) - coherence(arm 0); labels
permuted within cell. H218.1: P_sub = mean over project x subtype x purity-tertile cells of coherence(top) - coherence(bottom proliferation tertile within
cell); P_proj = the same over project x purity-tertile cells on the tumours of the eligible P_sub cells; proliferation arm labels permuted within cell.
Coherence, expression, proliferation, FGA, MYC, TP53 status and pseudo sets are those of A210 / A213 / A215."""
import argparse
import json
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

import a168_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=4)
ap.add_argument("--smoke", action="store_true"); ap.add_argument("--out", default=None)
args = ap.parse_args()
if args.smoke:
    assert args.out, "--smoke needs --out DIR"; args.nperm, args.npseudo = min(args.nperm, 40), min(args.npseudo, 3)
OUT = Path(args.out) if args.out else HERE.parent / "output"; OUT.mkdir(parents=True, exist_ok=True)
EXP = HERE.parents[1]
A210, A213, A215 = EXP / "A210_myc_coherence" / "output", EXP / "A213_heterogeneity_coherence" / "output", EXP / "A215_checkpoint_coherence" / "output"
CACHE = A210 / "tcga_cache"
API = DATA / "tcga" / "pancan_subtypes" / "api"
K, MIN_ARM, MDE_MAX, MIN_COVER = norm.ppf(0.95) + norm.ppf(0.80), 10, 0.25, 0.5
SEED = {"A217": zlib.crc32(b"A217") % (2 ** 31), "A218": zlib.crc32(b"A218") % (2 ** 31)}
NCHUNK = 20
PROJECTS = ["BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ",
            "SARC", "STAD", "TGCT", "THCA", "UCEC"]
GENES, LN2 = list(L.GENES), L.LN2


def tertile(x):
    """A166 tertile (rank method "first"); groups of fewer than 3 tumours are one middle group (they never reach the arm minimum)."""
    return L.A166.tertile(x) if len(x) >= 3 else np.ones(len(x), int)

PURITY = list(L.sh.PURITY)
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)
MEAS = {"H217.1": "tp53", "H217.2": "fga", "H217.3": "myc"}
LABEL = {"tp53": "TP53 altered vs intact", "fga": "FGA top vs bottom", "myc": "MYC GISTIC 2 vs 0"}


def study_of(p):
    return "coadread" if p in ("COAD", "READ") else p.lower()


# ---- subtype and purity tables (no coherence) -------------------------------------------------------------------------------------
def primaries(B, p):
    s = B.samples_of(p)
    return s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")


def load_lc(B, p):
    """log2 CPM of genes with mean CPM >= 1 for the primaries of p (A213 `load_lc`)."""
    ct = B.counts_of(p, primaries(B, p)); depth = ct.sum(axis=1)
    lc = np.log2(ct.div(depth, axis=0) * 1e6 + 1)
    return lc.loc[:, (2 ** lc - 1).mean() >= 1].astype(np.float32)


def zmean(lc, genes):
    g = [x for x in genes if x in lc.columns]
    z = (lc[g] - lc[g].mean()) / lc[g].std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1), len(g)


def subtype_table():
    """project, patient, subtype of the primaries of the 23 projects (empty = no label); parsed from the stored cBioPortal API JSON."""
    f = OUT / "subtypes_tcga.csv"
    if not f.exists():
        lab = {}
        for p in set(map(study_of, PROJECTS)):
            for r in json.load(open(API / f"{p}_tcga_pan_can_atlas_2018_SUBTYPE.json")):
                lab[r["patientId"]] = r["value"]
        rows = []
        for p in PROJECTS:
            for pt in pd.read_parquet(CACHE / f"{p}.parquet", columns=["patient"]).patient:
                rows.append(dict(project=p, patient=pt, subtype=lab.get(pt, "")))
        pd.DataFrame(rows).to_csv(f, index=False)
    return pd.read_csv(f, keep_default_na=False)


def purity_table():
    """A127 purity proxy per primary: mean z (within project, ddof 0) of log2 CPM of the 11 stromal / immune genes (genes with mean CPM >= 1)."""
    f = OUT / "purity_tcga.csv"
    if args.smoke and not f.exists():
        f = HERE.parent / "output" / "purity_tcga.csv"
    if not f.exists():
        B = data._builder("tcga_paired"); rows = []
        for p in PROJECTS:
            lc = load_lc(B, p); t = primaries(B, p); z, ng = zmean(lc, PURITY)
            rows.append(pd.DataFrame(dict(project=p, patient=t.patient.to_numpy(), purity=((z - z.mean()) / z.std()).to_numpy(), n_purity_genes=ng)))
            say(p, "purity genes", ng)
        pd.concat(rows).to_csv(f, index=False)
    return pd.read_csv(f)


SUB, PUR = subtype_table(), purity_table()
myc_s = pd.read_csv(A210 / "myc_cna_tcga.csv", index_col=0).myc_cna
fga_s = pd.read_csv(A213 / "fga_tcga.csv", index_col=0).fga
tp53_s = pd.read_csv(A215 / "status_tcga.csv", index_col=0).tp53


def project_table(p):
    """A210 cache of p (row order = primaries) with subtype, purity, MYC state, FGA, TP53 status."""
    d = pd.read_parquet(CACHE / f"{p}.parquet").reset_index(drop=True)
    s = SUB[SUB.project == p].set_index("patient").subtype.replace("", np.nan)
    d["sub"] = d.patient.map(s)
    d["purity"] = d.patient.map(PUR[PUR.project == p].set_index("patient").purity)
    d["myc"] = d["sample"].map(myc_s); d["fga"] = d["sample"].map(fga_s); d["tp53"] = d["sample"].map(tp53_s)
    return d


# ---- project eligibility ----------------------------------------------------------------------------------------------------------
tc = {p: project_table(p) for p in PROJECTS}
pj = []
for p in PROJECTS:
    d = tc[p]; lv = sorted(d["sub"].dropna().unique())
    pj.append(dict(project=p, n_primaries=len(d), n_subtype=int(d["sub"].notna().sum()), fraction=d["sub"].notna().mean(), n_levels=len(lv), levels="; ".join(lv)))
pj = pd.DataFrame(pj); pj["coverage_ok"] = pj.fraction >= MIN_COVER; pj["informative"] = pj.n_levels >= 2
pj["in_tests"] = pj.coverage_ok & pj.informative & (pj.project != "BRCA"); pj.to_csv(OUT / "projects.csv", index=False)
TP = list(pj.project[pj.in_tests])
say("projects entering:", TP)
Yof = lambda counts, genes: np.log((counts[list(genes)].to_numpy(float) + 0.5) / counts.depth.to_numpy(float)[:, None])
rng_noise = np.random.default_rng(0)
Ytc = {p: (rng_noise.normal(size=(len(tc[p]), len(GENES))) if args.smoke else Yof(tc[p], GENES)) for p in PROJECTS}


# ---- designs ----------------------------------------------------------------------------------------------------------------------
def diffs(Y, wf, cells, perms=None):
    out = []
    for k, c in enumerate(cells):
        l = c["lab"] if perms is None else c["lab"][perms[k]]
        out.append(L.coh(Y, c["idx"][l == 1], wf) - L.coh(Y, c["idx"][l == 0], wf))
    return np.array(out)


class CellD:
    """Labelled arms within cells (label 1 minus label 0); labels permuted within cell."""
    def __init__(s, cells): s.cells = cells
    def draw(s, rng): return [rng.permutation(len(c["lab"])) for c in s.cells]
    def values(s, Y, wf, st=None): return diffs(Y, wf, s.cells, st)


def arms(g, meas):
    """(arm 1, arm 0) positions of the tumours g (index = row positions of the project table)."""
    if meas == "tp53":
        return g.index[g.tp53 == 1].to_numpy(), g.index[g.tp53 == 0].to_numpy()
    if meas == "myc":
        return g.index[g.myc == 2].to_numpy(), g.index[g.myc == 0].to_numpy()
    e = tertile(g.fga.to_numpy(float)) if len(g) >= 3 else np.full(len(g), 1)
    return g.index[e == 2].to_numpy(), g.index[e == 0].to_numpy()


def cells_measure(d, p, meas, by_sub=True):
    """Project x subtype (by_sub) or project only x proliferation-tertile cells of the tumours of d with a subtype and the measure (tertile within
    the group over those tumours); >= MIN_ARM per arm. Returns (eligibility rows, eligible cells)."""
    m = d[d["sub"].notna() & d[meas].notna()].copy(); m["g"] = m["sub"] if by_sub else "all"
    m["t"] = m.groupby("g").prolif.transform(tertile); rows, cells = [], []
    for (s, t), g in m.groupby(["g", "t"]):
        a, b = arms(g, meas); ok = min(len(a), len(b)) >= MIN_ARM
        rows.append(dict(stratum=p, subtype=s, cell=int(t), n=len(g), n_arm1=len(a), n_arm0=len(b), eligible=ok))
        if ok:
            cells.append(dict(key=f"{p}|{s}|{int(t)}", stratum=p, subtype=s, idx=np.r_[a, b], lab=np.r_[np.ones(len(a), int), np.zeros(len(b), int)]))
    return pd.DataFrame(rows), cells


def top_bottom(g):
    e = tertile(g.prolif.to_numpy(float)) if len(g) >= 3 else np.full(len(g), 1)
    return g.index[e == 2].to_numpy(), g.index[e == 0].to_numpy()


def cells_prolif(d, p):
    """A218: purity tertile within project (over all primaries with a purity), P_sub cells = project x subtype x purity tertile, proliferation tertile
    within cell, top (label 1) vs bottom (label 0), >= MIN_ARM per arm; P_proj cells = project x purity tertile on the tumours of the eligible P_sub cells.
    Returns (eligibility rows, P_sub cells, P_proj cells, positions of the tumours of the eligible P_sub cells)."""
    d = d[d.purity.notna()].copy(); d["tp"] = tertile(d.purity.to_numpy(float)); m = d[d["sub"].notna()]
    rows, sub, keep = [], [], []
    for (s, t), g in m.groupby(["sub", "tp"]):
        a, b = top_bottom(g); ok = min(len(a), len(b)) >= MIN_ARM
        rows.append(dict(hypothesis="H218.1 P_sub", stratum=p, subtype=s, cell=int(t), n=len(g), n_arm1=len(a), n_arm0=len(b), eligible=ok))
        if ok:
            sub.append(dict(key=f"{p}|{s}|{int(t)}", stratum=p, subtype=s, idx=np.r_[a, b], lab=np.r_[np.ones(len(a), int), np.zeros(len(b), int)]))
            keep += list(g.index)
    proj = []
    for t, g in d.loc[sorted(keep)].groupby("tp"):
        a, b = top_bottom(g); ok = min(len(a), len(b)) >= MIN_ARM
        rows.append(dict(hypothesis="P_proj", stratum=p, subtype="", cell=int(t), n=len(g), n_arm1=len(a), n_arm0=len(b), eligible=ok))
        if ok:
            proj.append(dict(key=f"{p}|all|{int(t)}", stratum=p, subtype="all", idx=np.r_[a, b], lab=np.r_[np.ones(len(a), int), np.zeros(len(b), int)]))
    return pd.DataFrame(rows), sub, proj, np.array(sorted(keep), int)


def build(ps):
    """Cells of every hypothesis for the projects ps."""
    el, C = [], {h: {} for h in MEAS}; C["P_sub"], C["P_proj"], C["keep"] = {}, {}, {}
    for p in ps:
        for h, meas in MEAS.items():
            e, c = cells_measure(tc[p], p, meas); e.insert(0, "hypothesis", h); el.append(e); C[h][p] = c
        e, s, pr, kp = cells_prolif(tc[p], p); el.append(e); C["P_sub"][p], C["P_proj"][p], C["keep"][p] = s, pr, kp
    return pd.concat(el, ignore_index=True), C


# ---- nulls ------------------------------------------------------------------------------------------------------------------------
def null_chunk(blocks, seed, n):
    """blocks: (Y, wf, design); the statistic of a draw is the mean over all units of all blocks."""
    rng = np.random.default_rng(seed); out = np.empty(n)
    for i in range(n):
        out[i] = np.mean(np.concatenate([d.values(Y, wf, d.draw(rng)) for Y, wf, d in blocks]))
    return out


def draws(blocks, seed, tag, n):
    sizes = [n // NCHUNK + (j < n % NCHUNK) for j in range(NCHUNK)]
    return np.concatenate(Parallel(n_jobs=args.jobs)(delayed(null_chunk)(blocks, np.random.SeedSequence([seed, tag, j]), s) for j, s in enumerate(sizes)))


ncells = lambda cc: int(sum(len(c["lab"]) for c in cc))
elig, C = build(PROJECTS); elig.to_csv(OUT / "eligibility.csv", index=False)
HYPS = list(MEAS) + ["H218.1"]
CELLS = {**{h: C[h] for h in MEAS}, "H218.1": C["P_sub"]}
SEEDOF = lambda h: SEED["A218"] if h == "H218.1" else SEED["A217"]
blocks = {h: [(Ytc[p], L.WF, CellD(CELLS[h][p])) for p in TP if CELLS[h][p]] for h in HYPS}
nunits = lambda h: int(sum(len(d.cells) for _, _, d in blocks[h]))
ntum = lambda h: int(sum(ncells(d.cells) for _, _, d in blocks[h]))

cnt, comp = [], []
for p in TP + ["BRCA"]:
    d = tc[p]
    for h, meas in MEAS.items():
        t = d[d["sub"].notna() & d[meas].notna()]
        cnt.append(dict(hypothesis=h, project=p, n_with_subtype_and_measure=len(t), n_arm1_all=int(sum(len(arms(g, meas)[0]) for _, g in t.groupby("sub"))),
                        n_cells=len(C[h][p]), n_tumours_in_arms=ncells(C[h][p])))
        for c in C[h][p]:
            for lab, name in ((1, "arm1"), (0, "arm0")):
                comp.append(dict(hypothesis=h, project=p, subtype=c["subtype"], arm=name, n=int((c["lab"] == lab).sum())))
    cnt.append(dict(hypothesis="H218.1", project=p, n_with_subtype_and_measure=int((d["sub"].notna() & d.purity.notna()).sum()), n_cells=len(C["P_sub"][p]),
                    n_tumours_in_arms=ncells(C["P_sub"][p]), n_tumours_all_P_sub=len(C["keep"][p]), n_cells_P_proj=len(C["P_proj"][p]),
                    n_tumours_in_arms_P_proj=ncells(C["P_proj"][p])))
pd.DataFrame(cnt).to_csv(OUT / "counts.csv", index=False)
comp = pd.DataFrame(comp).groupby(["hypothesis", "project", "subtype", "arm"], as_index=False).n.sum(); comp.to_csv(OUT / "arm_composition.csv", index=False)
say({h: (len(blocks[h]), nunits(h), ntum(h)) for h in HYPS})

if args.stage == "power":
    nl = {}
    for i, h in enumerate(HYPS):
        if blocks[h]:
            nl[h] = draws(blocks[h], SEEDOF(h), i, args.nperm); say(h, "null done", nl[h].mean().round(5), nl[h].std(ddof=1).round(5))
    pd.DataFrame(nl).to_csv(OUT / "null_draws.csv", index=False)
    meas = {**{h: LABEL[m] for h, m in MEAS.items()}, "H218.1": "P_sub (proliferation top vs bottom within project x subtype x purity)"}
    pw = pd.DataFrame([dict(hypothesis=h, measure=meas[h], n_projects=len(blocks[h]), n_units=nunits(h), n_tumours=ntum(h), n_perm=args.nperm,
                            null_mean=nl[h].mean() if h in nl else np.nan, null_sd=nl[h].std(ddof=1) if h in nl else np.nan,
                            MDE_80pct=K * nl[h].std(ddof=1) if h in nl else np.nan) for h in HYPS])
    pw["scored"] = pw.MDE_80pct <= MDE_MAX
    pw["precondition"] = np.where(pw.hypothesis == "H218.1", "P_proj below its pseudo 5th percentile: evaluated in the tests stage", "")
    pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests (not run before the dated power note) ----------------------------------------------------------------------------------
pw = pd.read_csv(OUT / "power.csv").set_index("hypothesis"); nl = pd.read_csv(OUT / "null_draws.csv"); assert len(nl) == args.nperm
vals = lambda cells_of, ps: {p: CellD(cells_of[p]).values(Ytc[p], L.WF) for p in ps if cells_of[p]}
obs_c = {h: vals(CELLS[h], TP) for h in HYPS}                       # per project: per-cell T
obs = {h: float(np.concatenate(list(obs_c[h].values())).mean()) if obs_c[h] else np.nan for h in HYPS}
obs_proj_c = vals(C["P_proj"], TP)
P_proj = float(np.concatenate(list(obs_proj_c.values())).mean())

# original design (project x proliferation cells) on the same tumours: projects with an eligible matched cell, tumours with a subtype and the measure
orig_cells, orig_blocks, orig_null_mean, orig_obs = {}, {}, {}, {}
for i, (h, meas) in enumerate(MEAS.items()):
    pr = [p for p in TP if C[h][p]]; orig_cells[h] = {p: cells_measure(tc[p], p, meas, by_sub=False)[1] for p in pr}
    orig_blocks[h] = [(Ytc[p], L.WF, CellD(orig_cells[h][p])) for p in pr if orig_cells[h][p]]
    nlo = draws(orig_blocks[h], SEED["A217"], 100 + i, args.nperm); orig_null_mean[h] = float(nlo.mean())
    orig_obs[h] = float(np.mean(np.concatenate([d.values(Y, wf) for Y, wf, d in orig_blocks[h]])))

# pseudo-clock companion: 200 sets per project and amendment (A213 / A215 construction)
B = None if args.smoke else data._builder("tcga_paired")
pc = {h: [] for h in MEAS}; pc_sub, pc_proj = [], []
for j, p in enumerate(TP):
    if args.smoke:
        sets = [list(range(len(GENES)))] * args.npseudo
        mk = lambda g, rn=np.random.default_rng(j): (rn.normal(size=(len(tc[p]), len(GENES))), L.WF)
    else:
        lc = load_lc(B, p); t = primaries(B, p); assert (t.patient.to_numpy() == tc[p].patient.to_numpy()).all()
        lc = lc.set_axis(t.patient.to_numpy()); dep = tc[p].depth.to_numpy(float)[:, None]; pool = L.A166.pool_genes(lc)
        mk = lambda g: (np.log((2 ** lc[list(g)].to_numpy(float) - 1) / 1e6 + 0.5 / dep), L.wf_pseudo(L.WF, g))
    for amd, hs in (("A217", list(MEAS)), ("A218", ["P_sub", "P_proj"])):
        if args.smoke:
            sets = [None] * args.npseudo
        else:
            sets = L.pseudo_sets(lc, pool, GENES, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED[amd], 1000 + j])))
        for h in hs:
            cells = CELLS[h][p] if h in MEAS else C[h][p]
            if cells:
                arr = np.array([CellD(cells).values(*mk(g)) for g in sets])
                (pc[h] if h in MEAS else pc_sub if h == "P_sub" else pc_proj).append(arr)
    say(p, "pseudo done")
pseudo = {h: np.concatenate(pc[h], axis=1).mean(axis=1) for h in MEAS}
pseudo["H218.1"] = np.concatenate(pc_sub, axis=1).mean(axis=1); pseudo["P_proj"] = np.concatenate(pc_proj, axis=1).mean(axis=1)
pd.DataFrame(pseudo).to_csv(OUT / "pseudo_T.csv", index=False)
q05 = lambda h: float(np.quantile(pseudo[h], 0.05))

rows = []
for h in HYPS:
    sc = bool(pw.loc[h, "scored"])
    r = dict(hypothesis=h, measure=pw.loc[h, "measure"], scored_power=sc, stat=obs[h], null_mean=pw.loc[h, "null_mean"], null_sd=pw.loc[h, "null_sd"],
             MDE_80pct=pw.loc[h, "MDE_80pct"], excess=obs[h] - pw.loc[h, "null_mean"])
    if h == "H218.1":
        r["P_proj"] = P_proj; r["P_proj_pseudo_q05"] = q05("P_proj"); r["P_proj_below_pseudo_q05"] = bool(P_proj < q05("P_proj")); sc = sc and r["P_proj_below_pseudo_q05"]
    r["scored"] = sc; r["p_one_lower"] = pl(nl[h].to_numpy(), obs[h]) if sc else np.nan
    r.update(pseudo_percentile=100 * np.mean(pseudo[h] <= obs[h]), pseudo_q05=q05(h), below_pseudo_q05=bool(obs[h] < q05(h)))
    if h in MEAS:
        r.update(orig_T=orig_obs[h], orig_null_mean=orig_null_mean[h], orig_excess=orig_obs[h] - orig_null_mean[h],
                 excess_ratio=(obs[h] - pw.loc[h, "null_mean"]) / (orig_obs[h] - orig_null_mean[h]) if orig_obs[h] != orig_null_mean[h] else np.nan)
    rows.append(r)
res = pd.DataFrame(rows)
s = res.hypothesis.isin(list(MEAS)) & res.scored; ps = res.p_one_lower[s].to_numpy(); o = np.argsort(ps)
adj = np.empty(len(ps)); adj[o] = np.minimum(1, np.maximum.accumulate(np.sort(ps) * (len(ps) - np.arange(len(ps))))); res["p_holm"] = np.nan; res.loc[s, "p_holm"] = adj
res.loc[res.hypothesis == "H218.1", "p_holm"] = res.loc[res.hypothesis == "H218.1", "p_one_lower"]  # single hypothesis
res["verdict"] = np.where(~res.scored, "not scored", np.where(res.p_holm < 0.05, "pass", "fail"))
res["specificity"] = np.where(res.verdict == "pass", np.where(res.below_pseudo_q05, "clock-specific", "not specific to the clock"), "")

# A218 attenuation: A = P_sub / P_proj, delete-one-project jackknife (normal interval on the jackknife SE)
ps_c, pp_c = obs_c["H218.1"], obs_proj_c; pr = [p for p in TP if p in ps_c and p in pp_c]
tot = lambda c, ps_: np.concatenate([c[p] for p in ps_]).mean()
A = tot(ps_c, pr) / tot(pp_c, pr); Aj = np.array([tot(ps_c, [q for q in pr if q != p]) / tot(pp_c, [q for q in pr if q != p]) for p in pr])
se = np.sqrt((len(pr) - 1) / len(pr) * ((Aj - Aj.mean()) ** 2).sum())
att = pd.DataFrame([dict(P_sub=obs["H218.1"], P_proj=P_proj, A=A, A_jack_se=se, A_lo=A - 1.96 * se, A_hi=A + 1.96 * se, n_projects=len(pr),
                         A_jack_min=Aj.min(), A_jack_max=Aj.max())])
h18 = res[res.hypothesis == "H218.1"].iloc[0]
att["reading"] = ("not scored" if not h18.scored else
                  "holds within subtype" if (h18.verdict == "pass" and h18.specificity == "clock-specific" and A >= 0.5) else
                  "carried by subtype" if (h18.verdict == "fail" and A + 1.96 * se < 0.5) else "inconclusive")
att.to_csv(OUT / "attenuation.csv", index=False)
pd.DataFrame(dict(project_left_out=pr, A=Aj)).to_csv(OUT / "attenuation_jackknife.csv", index=False)
res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True); print(att.round(4).to_string(), flush=True)

# per-cell, per-project and per-subtype values (descriptive)
dc = []
for h in HYPS:
    for p in TP:
        for c, v in zip(CELLS[h][p], obs_c[h].get(p, [])):
            dc.append(dict(hypothesis=h, project=p, subtype=c["subtype"], cell=c["key"], n_arm1=int((c["lab"] == 1).sum()), n_arm0=int((c["lab"] == 0).sum()), T_cell=v))
for p in TP:
    for c, v in zip(C["P_proj"][p], obs_proj_c.get(p, [])):
        dc.append(dict(hypothesis="P_proj", project=p, subtype="all", cell=c["key"], n_arm1=int((c["lab"] == 1).sum()), n_arm0=int((c["lab"] == 0).sum()), T_cell=v))
dc = pd.DataFrame(dc); dc.to_csv(OUT / "descriptive_cells.csv", index=False)
dc.groupby(["hypothesis", "project"]).agg(n_cells=("T_cell", "size"), T=("T_cell", "mean")).reset_index().to_csv(OUT / "descriptive_per_project.csv", index=False)
dc.groupby(["hypothesis", "subtype"]).agg(n_cells=("T_cell", "size"), T=("T_cell", "mean")).reset_index().to_csv(OUT / "descriptive_per_subtype.csv", index=False)

# BRCA separately (descriptive; not blind): PAM50 x proliferation (matched) and proliferation (original) cells, same statistic, own null
br = []
for i, (h, meas) in enumerate(MEAS.items()):
    for design, by_sub in (("matched", True), ("original", False)):
        c = cells_measure(tc["BRCA"], "BRCA", meas, by_sub)[1]
        if c:
            b = [(Ytc["BRCA"], L.WF, CellD(c))]; nb = draws(b, SEED["A217"], 200 + 10 * i + by_sub, args.nperm); o = float(CellD(c).values(Ytc["BRCA"], L.WF).mean())
            br.append(dict(hypothesis=h, measure=LABEL[meas], design=design, n_cells=len(c), n_tumours=ncells(c), T=o, null_mean=nb.mean(), null_sd=nb.std(ddof=1),
                           excess=o - nb.mean(), p_one_lower=pl(nb, o)))
_, sb, pj_, _ = cells_prolif(tc["BRCA"], "BRCA")
for name, c in (("P_sub", sb), ("P_proj", pj_)):
    if c:
        b = [(Ytc["BRCA"], L.WF, CellD(c))]; nb = draws(b, SEED["A218"], 300 + (name == "P_sub"), args.nperm); o = float(CellD(c).values(Ytc["BRCA"], L.WF).mean())
        br.append(dict(hypothesis="A218", measure=name, design="matched" if name == "P_sub" else "original", n_cells=len(c), n_tumours=ncells(c), T=o, null_mean=nb.mean(),
                       null_sd=nb.std(ddof=1), excess=o - nb.mean(), p_one_lower=pl(nb, o)))
pd.DataFrame(br).to_csv(OUT / "brca_descriptive.csv", index=False); print(pd.DataFrame(br).round(4).to_string(), flush=True)
say("done tests")
