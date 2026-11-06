"""A219: clock coherence between molecular subtypes at matched proliferation (registered in paper/docs/registered_experiments.md, A219; mirrors A217).
  --stage power : eligibility (cells = project / cohort x proliferation tertile; subtypes with >= 15 tumours per cell; cells with >= 2 subtypes), counts,
                  permutation nulls (subtype labels permuted within cell) of S (H219.1) and D (H219.2, H219.3), MDE = 2.49 x null SD (H219.1 on S - null mean)
                  (projects.csv, eligibility.csv, counts.csv, power.csv, null_draws.csv). No coherence of any real group is computed.
  --stage tests : observed S / D, one-sided greater permutation p, Holm over H219.1-3, pseudo-clock companion (200 sets per TCGA project; METABRIC: stored A166
                  sets; SCAN-B: stored A127 mean- and variance-matched sets), descriptives (subtype coherences per cell / project, the same without proliferation
                  matching, residual proliferation imbalance, ER within LumB / Her2).
  --smoke --out DIR : both stages on random noise expression (real labels and cells) in DIR, to test the code path without touching real coherence.
    python -I compute_a219.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 4]
Cells: project (METABRIC, SCAN-B: cohort) x proliferation tertile (rank "first", within project over tumours with a subtype label). A subtype enters a cell with >= 15
tumours; a cell needs >= 2 entering subtypes (H219.1); for D, LumA and >= 1 of LumB / Her2 / Basal must enter (Normal-like left out), the entering others are pooled.
S = mean over cells of the tumour-weighted SD (population form) of the subtype coherences around the weighted mean; D = mean over cells of coh(LumA) - coh(others pooled).
Coherence, expression, proliferation and TCGA subtypes are those of A210 / A217; METABRIC and SCAN-B are the A166 / A127 cohort tables."""
import argparse
import os
import sys
import zlib
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
EXP = HERE.parents[1]
sys.path.insert(0, str(EXP / "A168_endocrine" / "scripts"))
sys.path.insert(0, str(ROOT / "paper"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a168_lib as L  # noqa: E402
import a127_lib as A127  # noqa: E402
import clocksync as cs  # noqa: E402
from common import data  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=4)
ap.add_argument("--smoke", action="store_true"); ap.add_argument("--out", default=None)
args = ap.parse_args()
if args.smoke:
    assert args.out, "--smoke needs --out DIR"; args.nperm, args.npseudo = min(args.nperm, 40), min(args.npseudo, 3)
OUT = Path(args.out) if args.out else HERE.parent / "output"; OUT.mkdir(parents=True, exist_ok=True)
A217 = EXP / "A217_subtype_matched" / "output"
CACHE = EXP / "A210_myc_coherence" / "output" / "tcga_cache"
K, MIN_SUB, MDE_MAX = norm.ppf(0.95) + norm.ppf(0.80), 15, 0.25
SEED = zlib.crc32(b"A219") % (2 ** 31)
NCHUNK = 20
GENES, LN2 = list(L.GENES), L.LN2
HYPS = ["H219.1", "H219.2", "H219.3"]
BREAST_D = ["LumA", "LumB", "Her2", "Basal"]
pg = lambda nl, o: cs.stats.perm_p(nl, o, "greater", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)
tertile = lambda x: L.A166.tertile(x) if len(x) >= 3 else np.ones(len(x), int)
rng_noise = np.random.default_rng(0)
noise = lambda n: rng_noise.normal(size=(n, len(GENES)))


# ---- cohorts: expression Y (n x 14), table with columns sub, prolif (+ er) --------------------------------------------------------
sub_tab = pd.read_csv(A217 / "subtypes_tcga.csv", keep_default_na=False)
pj = pd.read_csv(A217 / "projects.csv"); TP = list(pj.project[pj.in_tests])
TCGA_PROJECTS = TP + ["BRCA"]
CO = {}   # name -> dict(tab, Y, kind)
for p in TCGA_PROJECTS:
    d = pd.read_parquet(CACHE / f"{p}.parquet").reset_index(drop=True)
    s = sub_tab[sub_tab.project == p].set_index("patient").subtype.replace("", np.nan)
    d["sub"] = d.patient.map(s)
    if p == "BRCA":
        d["sub"] = d["sub"].str.replace("BRCA_", "", regex=False)
    Y = noise(len(d)) if args.smoke else np.log((d[GENES].to_numpy(float) + 0.5) / d.depth.to_numpy(float)[:, None])
    CO[p] = dict(tab=d[["sub", "prolif"]].copy(), Y=Y)
mb, Xmb = L.A166.load("metabric"); mb = mb.reset_index(drop=True); Xmb = Xmb.reset_index(drop=True)
mb["sub"] = mb.pam50.where(mb.pam50.notna()); mb["er"] = mb.er.astype(float)
CO["METABRIC"] = dict(tab=mb[["sub", "prolif", "er"]].copy(), Y=noise(len(mb)) if args.smoke else LN2 * Xmb[GENES].to_numpy(float), X=Xmb)
sb, Xsb = A127.load(); sb = sb.reset_index(drop=True); Xsb = Xsb.reset_index(drop=True)
sb["sub"] = sb.pam50
CO["SCANB"] = dict(tab=sb[["sub", "prolif", "er"]].copy(), Y=noise(len(sb)) if args.smoke else LN2 * Xsb[GENES].to_numpy(float), X=Xsb)
HYP_COH = {"H219.1": TP, "H219.2": ["METABRIC", "BRCA"], "H219.3": ["SCANB"]}
KIND = {"H219.1": "S", "H219.2": "D", "H219.3": "D"}


# ---- cells ------------------------------------------------------------------------------------------------------------------------
def make_cells(name, kind, matched=True):
    """Cells of cohort `name` for statistic kind 'S' or 'D'. Returns (eligibility rows, cells)."""
    t = CO[name]["tab"]; t = t[t["sub"].notna()].copy(); t["pos"] = t.index
    t["g"] = t.groupby(lambda i: 0).prolif.transform(tertile) if matched else 1
    rows, cells = [], []
    for g, c in t.groupby("g"):
        n = c["sub"].value_counts(); enter = [s for s in n.index if n[s] >= MIN_SUB and (kind == "S" or s in BREAST_D)]
        ok = len(enter) >= 2 and (kind == "S" or ("LumA" in enter and len(enter) >= 2))
        for s in sorted(n.index):
            rows.append(dict(kind=kind, cohort=name, cell=int(g), subtype=s, n=int(n[s]), enters=s in enter, cell_eligible=ok))
        if ok:
            k = c[c["sub"].isin(enter)]
            lab = k["sub"].to_numpy() if kind == "S" else (k["sub"] == "LumA").astype(int).to_numpy()
            cells.append(dict(key=f"{name}|{int(g)}", cohort=name, idx=k.pos.to_numpy(), lab=lab, sub=k["sub"].to_numpy()))
    return pd.DataFrame(rows), cells


def cell_value(Y, wf, c, lab, kind):
    if kind == "D":
        return L.coh(Y, c["idx"][lab == 1], wf) - L.coh(Y, c["idx"][lab == 0], wf)
    lv = np.unique(lab); coh = np.array([L.coh(Y, c["idx"][lab == s], wf) for s in lv]); w = np.array([(lab == s).sum() for s in lv], float)
    m = (w * coh).sum() / w.sum()
    return float(np.sqrt((w * (coh - m) ** 2).sum() / w.sum()))


class CellD:
    """Subtype labels within cells; labels permuted within cell. values -> one S (weighted SD) or D (LumA - others) per cell."""
    def __init__(s, cells, kind): s.cells, s.kind = cells, kind
    def draw(s, rng): return [rng.permutation(len(c["lab"])) for c in s.cells]
    def values(s, Y, wf, st=None):
        return np.array([cell_value(Y, wf, c, c["lab"] if st is None else c["lab"][st[k]], s.kind) for k, c in enumerate(s.cells)])


def null_chunk(blocks, seed, n):
    rng = np.random.default_rng(seed); out = np.empty(n)
    for i in range(n):
        out[i] = np.mean(np.concatenate([d.values(Y, wf, d.draw(rng)) for Y, wf, d in blocks]))
    return out


def draws(blocks, tag, n):
    sizes = [n // NCHUNK + (j < n % NCHUNK) for j in range(NCHUNK)]
    return np.concatenate(Parallel(n_jobs=args.jobs)(delayed(null_chunk)(blocks, np.random.SeedSequence([SEED, tag, j]), s) for j, s in enumerate(sizes)))


ncells = lambda cc: int(sum(len(c["lab"]) for c in cc))
elig, CELLS = [], {h: {} for h in HYPS}
for h in HYPS:
    for nm in HYP_COH[h]:
        e, c = make_cells(nm, KIND[h]); e.insert(0, "hypothesis", h); elig.append(e); CELLS[h][nm] = c
elig = pd.concat(elig, ignore_index=True); elig.to_csv(OUT / "eligibility.csv", index=False)
blocks = {h: [(CO[nm]["Y"], L.WF, CellD(CELLS[h][nm], KIND[h])) for nm in HYP_COH[h] if CELLS[h][nm]] for h in HYPS}
nunits = lambda h: int(sum(len(d.cells) for _, _, d in blocks[h]))
ntum = lambda h: int(sum(ncells(d.cells) for _, _, d in blocks[h]))
cnt = []
for h in HYPS:
    for nm in HYP_COH[h]:
        cc = CELLS[h][nm]; e = elig[(elig.hypothesis == h) & (elig.cohort == nm)]
        cnt.append(dict(hypothesis=h, cohort=nm, n_labelled=int(e.n.sum()), n_cells_possible=e.cell.nunique(), n_cells=len(cc), n_tumours_in_cells=ncells(cc),
                        subtypes_per_cell="; ".join(str(len(np.unique(c["sub"]))) for c in cc), levels="; ".join(sorted({s for c in cc for s in np.unique(c["sub"])}))))
pd.DataFrame(cnt).to_csv(OUT / "counts.csv", index=False)
say({h: (len(blocks[h]), nunits(h), ntum(h)) for h in HYPS})

if args.stage == "power":
    nl = {}
    for i, h in enumerate(HYPS):
        if blocks[h]:
            nl[h] = draws(blocks[h], i, args.nperm); say(h, "null done", nl[h].mean().round(5), nl[h].std(ddof=1).round(5))
    pd.DataFrame(nl).to_csv(OUT / "null_draws.csv", index=False)
    meas = {"H219.1": "S: SD of subtype coherences, TCGA non-BRCA", "H219.2": "D: LumA - others, METABRIC + TCGA BRCA", "H219.3": "D: LumA - others, SCAN-B"}
    pw = pd.DataFrame([dict(hypothesis=h, measure=meas[h], n_cohorts=len(blocks[h]), n_units=nunits(h), n_tumours=ntum(h), n_perm=args.nperm,
                            null_mean=nl[h].mean() if h in nl else np.nan, null_sd=nl[h].std(ddof=1) if h in nl else np.nan,
                            MDE_80pct=K * nl[h].std(ddof=1) if h in nl else np.nan) for h in HYPS])
    pw["scored"] = pw.MDE_80pct <= MDE_MAX
    pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests (not run before the dated power note) ----------------------------------------------------------------------------------
pw = pd.read_csv(OUT / "power.csv").set_index("hypothesis"); nl = pd.read_csv(OUT / "null_draws.csv"); assert len(nl) == args.nperm
obs_c = {h: {nm: CellD(CELLS[h][nm], KIND[h]).values(CO[nm]["Y"], L.WF) for nm in HYP_COH[h] if CELLS[h][nm]} for h in HYPS}
obs = {h: float(np.concatenate(list(obs_c[h].values())).mean()) if obs_c[h] else np.nan for h in HYPS}

# pseudo-clock companion: TCGA 200 sets per project (A213 / A217 construction, seeds crc32("A219")); METABRIC stored A166 sets; SCAN-B stored A127 sets
def primaries(B, p):
    s = B.samples_of(p)
    return s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")


def load_lc(B, p):
    """log2 CPM of genes with mean CPM >= 1 for the primaries of p (A213 `load_lc`)."""
    ct = B.counts_of(p, primaries(B, p)); depth = ct.sum(axis=1)
    lc = np.log2(ct.div(depth, axis=0) * 1e6 + 1)
    return lc.loc[:, (2 ** lc - 1).mean() >= 1].astype(np.float32)


sets_mb = pd.read_csv(L.A166.OUT / "pseudo_sets_metabric.csv").to_numpy()[: args.npseudo].tolist()
sets_sb = pd.read_csv(EXP / "A127_scanb" / "output" / "pseudo_sets_variance_matched.csv").to_numpy()[: args.npseudo].tolist()
pc = {h: [] for h in HYPS}
for j, nm in enumerate(TCGA_PROJECTS):
    hs = [h for h in HYPS if nm in HYP_COH[h] and CELLS[h][nm]]
    if not hs:
        continue
    if args.smoke:
        rn = np.random.default_rng(j); sets = [None] * args.npseudo
        mk = lambda g: (rn.normal(size=(len(CO[nm]["Y"]), len(GENES))), L.WF)
    else:
        B = data._builder("tcga_paired"); lc = load_lc(B, nm); t = primaries(B, nm)
        d = pd.read_parquet(CACHE / f"{nm}.parquet").reset_index(drop=True); assert (t.patient.to_numpy() == d.patient.to_numpy()).all()
        lc = lc.set_axis(t.patient.to_numpy()); dep = d.depth.to_numpy(float)[:, None]
        sets = L.pseudo_sets(lc, L.A166.pool_genes(lc), GENES, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 1000 + j])))
        mk = lambda g: (np.log((2 ** lc[list(g)].to_numpy(float) - 1) / 1e6 + 0.5 / dep), L.wf_pseudo(L.WF, g))
    for h in hs:
        pc[h].append(np.array([CellD(CELLS[h][nm], KIND[h]).values(*mk(g)) for g in sets]))
    say(nm, "pseudo done")
for h, nm, sets, key in (("H219.2", "METABRIC", sets_mb, "X"), ("H219.3", "SCANB", sets_sb, "X")):
    if CELLS[h][nm]:
        des = CellD(CELLS[h][nm], KIND[h]); X = CO[nm][key]
        if args.smoke:
            pc[h].append(np.array([des.values(rng_noise.normal(size=(len(CO[nm]["Y"]), len(GENES))), L.WF) for _ in range(args.npseudo)]))
        else:
            pc[h].append(np.array([des.values(LN2 * X[list(g)].to_numpy(float), L.wf_pseudo(L.WF, g)) for g in sets]))
pseudo = {h: np.concatenate(pc[h], axis=1).mean(axis=1) for h in HYPS if pc[h]}
pd.DataFrame(pseudo).to_csv(OUT / "pseudo_stat.csv", index=False)
q95 = {h: float(np.quantile(pseudo[h], 0.95)) for h in pseudo}

rows = []
for h in HYPS:
    sc = bool(pw.loc[h, "scored"]); r = dict(hypothesis=h, measure=pw.loc[h, "measure"], scored=sc, stat=obs[h], null_mean=pw.loc[h, "null_mean"], null_sd=pw.loc[h, "null_sd"],
                                             MDE_80pct=pw.loc[h, "MDE_80pct"], excess=obs[h] - pw.loc[h, "null_mean"])
    r["p_one_greater"] = pg(nl[h].to_numpy(), obs[h]) if sc else np.nan
    r.update(pseudo_percentile=100 * np.mean(pseudo[h] <= obs[h]), pseudo_q95=q95[h], pseudo_mean=float(pseudo[h].mean()), above_pseudo_q95=bool(obs[h] > q95[h]))
    rows.append(r)
res = pd.DataFrame(rows)
s = res.scored.to_numpy(); ps = res.p_one_greater[s].to_numpy(); o = np.argsort(ps)
adj = np.empty(len(ps)); adj[o] = np.minimum(1, np.maximum.accumulate(np.sort(ps) * (len(ps) - np.arange(len(ps))))); res["p_holm"] = np.nan; res.loc[s, "p_holm"] = adj
res["verdict"] = np.where(~res.scored, "not scored", np.where(res.p_holm < 0.05, "pass", "fail"))
res["specificity"] = np.where(res.verdict == "pass", np.where(res.above_pseudo_q95, "clock-specific", "not specific to the clock"), "")
res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True)

# ---- descriptives -----------------------------------------------------------------------------------------------------------------
dc, imb = [], []
for h in HYPS:
    for nm in HYP_COH[h]:
        Y = CO[nm]["Y"]; t = CO[nm]["tab"]
        for c in CELLS[h][nm]:
            lab = c["lab"]; g = c["sub"] if KIND[h] == "S" else np.where(lab == 1, "LumA", "others")
            for s_ in np.unique(g):
                dc.append(dict(hypothesis=h, cohort=nm, cell=c["key"], group=s_, n=int((g == s_).sum()), coherence=L.coh(Y, c["idx"][g == s_], L.WF),
                               mean_prolif=float(t.prolif.to_numpy()[c["idx"][g == s_]].mean())))
            imb.append(dict(hypothesis=h, cohort=nm, cell=c["key"], stat_cell=cell_value(Y, L.WF, c, lab, KIND[h]),
                            sd_prolif_between_groups=float(np.std([t.prolif.to_numpy()[c["idx"][g == s_]].mean() for s_ in np.unique(g)]))))
dc = pd.DataFrame(dc); dc.to_csv(OUT / "descriptive_cells.csv", index=False)
dc.groupby(["hypothesis", "cohort", "group"]).apply(lambda x: pd.Series(dict(n_cells=len(x), n=x.n.sum(), coherence_w=np.average(x.coherence, weights=x.n))),
                                                    include_groups=False).reset_index().to_csv(OUT / "descriptive_per_group.csv", index=False)
imb = pd.DataFrame(imb); imb.to_csv(OUT / "descriptive_cells_stat_and_prolif_imbalance.csv", index=False)
imb.groupby(["hypothesis", "cohort"]).agg(n_cells=("stat_cell", "size"), stat=("stat_cell", "mean"), prolif_sd_between_groups=("sd_prolif_between_groups", "mean")) \
    .reset_index().to_csv(OUT / "descriptive_per_cohort.csv", index=False)

# the same without proliferation matching: one cell per project / cohort, own null
un = []
for i, h in enumerate(HYPS):
    bl, cc_all = [], []
    for nm in HYP_COH[h]:
        e, c = make_cells(nm, KIND[h], matched=False)
        if c:
            bl.append((CO[nm]["Y"], L.WF, CellD(c, KIND[h]))); cc_all.append(c)
    if bl:
        nu = draws(bl, 100 + i, args.nperm); o = float(np.mean(np.concatenate([d.values(Y, wf) for Y, wf, d in bl])))
        un.append(dict(hypothesis=h, n_cells=sum(len(d.cells) for _, _, d in bl), n_tumours=sum(ncells(d.cells) for _, _, d in bl), stat=o, null_mean=nu.mean(), null_sd=nu.std(ddof=1),
                       excess=o - nu.mean(), p_one_greater=pg(nu, o)))
pd.DataFrame(un).to_csv(OUT / "descriptive_unmatched.csv", index=False); print(pd.DataFrame(un).round(4).to_string(), flush=True)

# ER within LumB and Her2 (breast cohorts with ER: METABRIC, SCAN-B; arms >= MIN_SUB; descriptive)
er = []
for nm in ("METABRIC", "SCANB"):
    t = CO[nm]["tab"]
    for s_ in ("LumB", "Her2"):
        g = t[(t["sub"] == s_) & t.er.isin([0.0, 1.0])]; a, b = g.index[g.er == 1].to_numpy(), g.index[g.er == 0].to_numpy()
        er.append(dict(cohort=nm, subtype=s_, n_er_pos=len(a), n_er_neg=len(b), coh_er_pos=L.coh(CO[nm]["Y"], a, L.WF) if len(a) >= MIN_SUB else np.nan,
                       coh_er_neg=L.coh(CO[nm]["Y"], b, L.WF) if len(b) >= MIN_SUB else np.nan))
pd.DataFrame(er).to_csv(OUT / "descriptive_er_within_subtype.csv", index=False); print(pd.DataFrame(er).round(4).to_string(), flush=True)
say("done tests")
