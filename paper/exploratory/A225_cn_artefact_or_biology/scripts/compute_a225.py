"""A225: copy number at the clock loci, readout artefact or clock biology? (registered in paper/docs/registered_experiments.md, A225)
  --stage sim   : part A (gate). Semi-synthetic METABRIC and TCGA BRCA on the real GISTIC profiles; scenarios N (identical clocks) and D (non-LumA lam x 0.7);
                  estimators raw (A219), adjusted (A224, betas re-estimated), masked (GISTIC != 0 set missing, pairwise-complete coherence); 200 replicates.
                  sim_calibration.csv, sim_replicates.csv, sim_summary.csv (gate, V225.1).
  --stage power : MDEs (permutation nulls, labels permuted; no observed statistic is computed) of H225.1 (gate-passing estimator), H225.2 (ARNTL), H225.3
                  (other clock genes): power.csv, eligibility.csv, null_draws.csv.
  --stage tests : observed H225.1 / H225.2 / H225.3 with permutation p, pseudo companions, Holm. NOT RUN before the dated power note.
  --smoke --out DIR : every stage on noise expression (real labels, cells and copy numbers), few replicates / permutations, in DIR.
    python -I compute_a225.py --stage sim|power|tests [--jobs 8] [--nrep 200] [--nperm 2000] [--npseudo 200]"""
import argparse
import os
import sys
import tempfile
import zlib
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
EXP = HERE.parents[1]
A224S = EXP / "A224_cis_copy_number" / "scripts"
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(EXP / "A168_endocrine" / "scripts")); sys.path.insert(0, str(ROOT / "paper"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a168_lib as L  # noqa: E402
import a225_lib as A  # noqa: E402
import clocksync as cs  # noqa: E402
from common import data  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["sim", "power", "tests"])
ap.add_argument("--nrep", type=int, default=200); ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=8); ap.add_argument("--smoke", action="store_true"); ap.add_argument("--out", default=None)
args = ap.parse_args()
assert args.jobs <= 8, "budget: at most 8 workers"
if args.smoke:
    assert args.out, "--smoke needs --out DIR"; args.nrep, args.nperm, args.npseudo = min(args.nrep, 4), min(args.nperm, 60), min(args.npseudo, 3)
OUT = Path(args.out) if args.out else HERE.parent / "output"; OUT.mkdir(parents=True, exist_ok=True)
SEED = zlib.crc32(b"A225") % (2 ** 31)
SEED219 = zlib.crc32(b"A219") % (2 ** 31)
K, MDE_H1, MDE_H2, MIN_ARM, MIN_PAIR, DAMP = norm.ppf(0.95) + norm.ppf(0.80), 0.25, 0.15, 10, 5, 0.7
OBS_D_A219 = 0.276                                       # V225.1 denominator: A219 H219.2 observed LumA - others (registered)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)

# ---- reuse A224 (cells, copy numbers, dosage fits, adjusted expression): its code up to the stage branch, run in a scratch output dir -----------------
src = (A224S / "compute_a224.py").read_text(); src = src[: src.index('HYPS = ["H224.1"')]
sys.argv = ["compute_a224.py", "--stage", "power", "--out", tempfile.mkdtemp(prefix="a225_a224_"), "--jobs", "1"] + (["--smoke"] if args.smoke else [])
ns = {"__file__": str(A224S / "compute_a224.py"), "__name__": "a224"}
exec(compile(src, str(A224S / "compute_a224.py"), "exec"), ns)
GENES, WF, LN2, CW = ns["GENES"], L.WF, ns["LN2"], ns["CW"]
tc, Y0, CNt, CNm, mb, Ymb0, Ymb, Yt, c_mb, c_br, cn_of, dos = (ns[k] for k in ("tc", "Y0", "CNt", "CNm", "mb", "Ymb0", "Ymb", "Yt", "c_mb", "c_br", "cn_of", "dos"))
PROJECTS, tertile, TCGA219 = ns["PROJECTS"], ns["tertile"], ns["TCGA219"]
CNM = lambda CN, ids: np.column_stack([cn_of(CN, g, ids) for g in GENES])          # n x 14, NaN where no call
COH = {"METABRIC": dict(Y=Ymb0, Yadj=Ymb, CN=CNM(CNm, mb.patient.to_numpy()), cells=c_mb, sub=mb["sub"], dcoh="METABRIC"),
       "BRCA": dict(Y=Y0["BRCA"], Yadj=Yt["BRCA"], CN=CNM(CNt, tc["BRCA"]["sample"].to_numpy()), cells=c_br, sub=tc["BRCA"]["sub"], dcoh="TCGA_BRCA")}
say("A224 objects loaded; cells:", {k: len(v["cells"]) for k, v in COH.items()})

def cells_for(name, est):
    """Cells of cohort `name` usable by the estimator: all for raw / adjusted; for masked those whose two arms have >= MIN_PAIR complete rows for every gene pair."""
    c = COH[name]
    if est != "masked":
        return c["cells"]
    M = (np.nan_to_num(c["CN"]) == 0).astype(float)
    return [cl for cl in c["cells"] if A.min_pair_n(M, cl["idx"], cl["lab"]) >= MIN_PAIR]


# ---- cohort parameters of the simulation ---------------------------------------------------------------------------------------------------------
def make_par(name):
    c = COH[name]; Y, CN = c["Y"], c["CN"]; ph = A.leave_out_phases(Y, CW); d = dos[dos.cohort == c["dcoh"]].set_index("gene")
    fits = [A.fit_cn(Y[:, j], CN[:, j], ph[:, j]) for j in range(len(GENES))]
    for j, g in enumerate(GENES):                          # our OLS reproduces the A224 betas
        if g in d.index and d.loc[g, "note"] == "":
            assert abs(fits[j][0] - d.loc[g, "beta_ln"]) < 1e-6 * max(1, abs(fits[j][0])) + 1e-9, (name, g, fits[j][0], d.loc[g, "beta_ln"])
    col = lambda k: np.array([d.loc[g, k] * LN2 if (g in d.index and d.loc[g, k] == d.loc[g, k]) else 0.0 for g in GENES])
    sub = c["sub"]
    return dict(n=len(Y), amp=np.asarray(WF.amp, float), psi=np.asarray(WF.psi, float), cw=CW, wf=WF, cn=CN, cn0=np.nan_to_num(CN), beta=col("beta_inter"), ic=col("int_cos"),
                sd=np.array([f[1] for f in fits]), is_other=(sub.notna() & (sub != "LumA")).to_numpy(), nfit=np.array([f[2] for f in fits]), **{"is": col("int_sin")})


# ---- pseudo genes (own CN, own beta, own residual SD; no clock term) -------------------------------------------------------------------------------
def pseudo_source(name):
    """(sets, gy, gcn) of the cohort's first args.npseudo A219 pseudo sets (stored A166 sets for METABRIC; A219-seed sets on the BRCA counts for BRCA):
    gy(g) = ln expression of pseudo gene g (n), gcn(g) = its GISTIC state (n, NaN where no call)."""
    c = COH[name]
    if args.smoke:
        rg = np.random.default_rng(5); return [GENES] * args.npseudo, (lambda g: rg.normal(size=len(c["Y"]))), (lambda g: c["CN"][:, GENES.index("PER1")])
    if name == "METABRIC":
        sets = pd.read_csv(L.A166.OUT / "pseudo_sets_metabric.csv").to_numpy()[: args.npseudo].tolist(); Xmb = L.A166.load("metabric")[1].reset_index(drop=True)
        CNp = ns["cn_metabric"](sorted({g for s in sets for g in s}))
        return sets, (lambda g: L.LN2 * Xmb[g].to_numpy(float)), (lambda g: cn_of(CNp, g, mb.patient.to_numpy()))
    lc = load_lc("BRCA"); assert (lc.index.to_numpy() == tc["BRCA"].patient.to_numpy()).all()
    sets = L.pseudo_sets(lc, L.A166.pool_genes(lc), GENES, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED219, 1000 + TCGA219.index("BRCA")])))
    CNp = ns["cn_tcga"](sorted({g for s in sets for g in s})); dep = tc["BRCA"].depth.to_numpy(float)
    return sets, (lambda g: np.log((2 ** lc[g].to_numpy(float) - 1) / 1e6 + 0.5 / dep)), (lambda g: cn_of(CNp, g, tc["BRCA"]["sample"].to_numpy()))


def pseudo_ctx(name):
    """Simulation inputs of the pseudo sets: list of dicts (cn, cn0, beta, sd) with own CN, own beta (cohort fit, phase of the 14 real clock genes) and residual SD."""
    c = COH[name]; ph14 = np.angle((c["Y"] - c["Y"].mean(0)) @ CW); sets, gy, gcn = pseudo_source(name); cache = {}
    def fit(g):
        if g not in cache:
            cn = gcn(g); b, sd, _ = A.fit_cn(gy(g), cn, ph14); cache[g] = (cn, b, sd)
        return cache[g]
    out = []
    for s in sets:
        f = [fit(g) for g in s]; cn = np.column_stack([x[0] for x in f])
        out.append(dict(cn=cn, cn0=np.nan_to_num(cn), beta=np.array([x[1] for x in f]), sd=np.array([x[2] for x in f])))
    return out


def primaries(B, p):
    s = B.samples_of(p); return s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")


def load_lc(p):
    """log2 CPM (genes with mean CPM >= 1) of the primaries of project p, rows = patients (A213 `load_lc`)."""
    B = data._builder("tcga_paired"); pr = primaries(B, p); ct = B.counts_of(p, pr); lc = np.log2(ct.div(ct.sum(axis=1), axis=0) * 1e6 + 1)
    return lc.loc[:, (2 ** lc - 1).mean() >= 1].astype(np.float32).set_axis(pr.patient.to_numpy())


# =============================================================================================================================================
if args.stage == "sim":
    say("sim: parameters"); target = pd.read_csv(EXP / "A219_subtype_coherence" / "output" / "descriptive_per_group.csv")
    target = {n: float(target[(target.hypothesis == "H219.2") & (target.cohort == n) & (target.group == "LumA")].coherence_w.iloc[0]) for n in COH}
    ctx, cal, tru = {}, [], {}
    for k, name in enumerate(COH):
        par = make_par(name); lam, f, extra = A.calibrate(par, COH[name]["cells"], WF, target[name], SEED + k, nrep=5 if args.smoke else 20)
        assert np.isfinite(lam), (name, "calibration bracket", f, extra)
        cb = {e: cells_for(name, e) for e in A.EST}; ntrue = 20 if args.smoke else 200
        tru[name] = {e: (A.true_D(par, cb[e], lam, DAMP, WF, ntrue, SEED + 100 + k) if cb[e] else np.array([])) for e in ("raw", "masked")}   # adjusted uses the raw cells
        tru[name]["adjusted"] = tru[name]["raw"]
        cal.append(dict(cohort=name, target_luma_raw_coherence=target[name], lam=lam, achieved=f, true_D_under_D_mean_cells=tru[name]["raw"].mean(), true_D_under_D_masked_cells=tru[name]["masked"].mean() if len(tru[name]["masked"]) else np.nan, n=par["n"], n_cells=len(COH[name]["cells"]),
                        resid_sd_mean=par["sd"].mean(), **{f"sd_{g}": s for g, s in zip(GENES, par["sd"])}))
        say(name, "lam", round(lam, 4), "achieved", round(f, 4), "target", round(target[name], 4), "true D under D (all cells)", round(tru[name]["raw"].mean(), 4), "(masked cells)", len(tru[name]["masked"]))
        ctx[name] = dict(par=par, cells_by={e: cells_for(name, e) for e in A.EST}, pseudo=pseudo_ctx(name)); say(name, "pseudo ctx", len(ctx[name]["pseudo"]))
    pd.DataFrame(cal).to_csv(OUT / "sim_calibration.csv", index=False)
    lams = {c["cohort"]: c["lam"] for c in cal}
    per = 2; chunks = [list(range(i, min(i + per, args.nrep))) for i in range(0, args.nrep, per)]; rows = []
    gen = Parallel(n_jobs=args.jobs, return_as="generator")(delayed(A.replicate_chunk)(rs, ctx, lams, DAMP, args.nperm, SEED) for rs in chunks)
    for i, r in enumerate(gen):
        rows += r
        if i % 5 == 0 or i == len(chunks) - 1:
            say(f"replicates done {min((i + 1) * per, args.nrep)}/{args.nrep}")
            pd.DataFrame(rows).to_csv(OUT / "sim_replicates.csv", index=False)
    rep = pd.DataFrame(rows); rep.to_csv(OUT / "sim_replicates.csv", index=False)
    nce = lambda e, n: len(ctx[n]["cells_by"][e])
    tcoh = lambda e, n: tru[n][e].mean() if len(tru[n][e]) else np.nan
    tpool = {e: np.concatenate([tru[n][e] for n in COH]).mean() for e in A.EST}
    out = []
    for (scen, e, lev), g in rep.groupby(["scenario", "estimator", "level"], sort=False):
        td = 0.0 if scen == "N" else (tpool[e] if lev == "pooled" else tcoh(e, lev))
        out.append(dict(scenario=scen, estimator=e, level=lev, n_rep=len(g), n_cells=sum(nce(e, n) for n in (COH if lev == "pooled" else [lev])), mean_D=g.D.mean(), sd_D=g.D.std(ddof=1), reject_rate_p05=(g.p < 0.05).mean(), true_D=td,
                        relative_recovery=g.D.mean() / td if scen == "D" else np.nan, pseudo_mean_D=g.pseudo_D.mean(), pseudo_sd_D=g.pseudo_D.std(ddof=1)))
    sm = pd.DataFrame(out)
    gate = []
    for e in A.EST:
        n_ = sm[(sm.scenario == "N") & (sm.estimator == e) & (sm.level == "pooled")].iloc[0]; d_ = sm[(sm.scenario == "D") & (sm.estimator == e) & (sm.level == "pooled")].iloc[0]
        ok = dict(N_bias=abs(n_.mean_D) <= 0.03, N_reject=n_.reject_rate_p05 <= 0.07, D_recovery=abs(d_.relative_recovery - 1) <= 0.25)
        gate.append(dict(estimator=e, **ok, gate_pass=all(ok.values())))
    yn = lambda v: "yes" if v else "no"
    gate = pd.DataFrame(gate).set_index("estimator").map(yn); sm = sm.join(gate, on="estimator"); sm.loc[sm.level != "pooled", gate.columns] = ""
    sm["V225_1_raw_D_N_over_0.276"] = np.where((sm.scenario == "N") & (sm.estimator == "raw") & (sm.level == "pooled"), sm.mean_D / OBS_D_A219, np.nan)
    sm.to_csv(OUT / "sim_summary.csv", index=False); print(sm.round(4).to_string(), flush=True); say("done sim"); sys.exit(0)


# =============================================================================================================================================
# real-data blocks
def h1_blocks(est):
    """H225.1 blocks [(Y, M, cells, wf)] of the estimator 'masked' | 'adjusted' | 'raw'."""
    out = []
    for name, c in COH.items():
        Y = c["Yadj"] if est == "adjusted" else c["Y"]
        M = (np.nan_to_num(c["CN"]) == 0).astype(float) if est == "masked" else np.ones_like(Y)
        out.append((Y, M, cells_for(name, est), WF))
    return out


ft = pd.read_csv(EXP / "A213_heterogeneity_coherence" / "output" / "fga_tcga.csv", index_col=0).fga
OTH = {g: [x for x in GENES if x != g] for g in GENES}


def trans_cells(g):
    """Cells of gene g in every project: GISTIC(g) <= -1 (label 1) vs == 0 (label 0), within-project FGA tertile x proliferation tertile (over tumours in an arm with FGA and
    proliferation); arms >= MIN_ARM; every gene pair of the 13 others complete in >= MIN_PAIR tumours of the cell (after masking). Returns eligibility rows, {project: (Y, M, cells)}."""
    rows, blocks = [], {}
    oi = [GENES.index(x) for x in OTH[g]]
    for p in PROJECTS:
        d = tc[p]; cn = cn_of(CNt, g, d["sample"].to_numpy()); fga = ft.reindex(d["sample"]).to_numpy(float) if not args.smoke else np.random.default_rng(1).random(len(d))
        arm = np.where(cn <= -1, 1.0, np.where(cn == 0, 0.0, np.nan)); ok = np.isfinite(arm) & np.isfinite(fga) & d.prolif.notna().to_numpy()
        t = pd.DataFrame(dict(pos=np.where(ok)[0], arm=arm[ok].astype(int), f=tertile(fga[ok]), q=tertile(d.prolif.to_numpy(float)[ok])))
        CN13 = CNM(CNt, d["sample"].to_numpy())[:, oi]; M = (np.nan_to_num(CN13) == 0).astype(float); cells = []
        for (f_, q_), gg in t.groupby(["f", "q"]):
            n1, n0 = int((gg.arm == 1).sum()), int((gg.arm == 0).sum()); mp = A.min_pair_n(M, gg.pos.to_numpy(), gg.arm.to_numpy()) if min(n1, n0) > 0 else 0
            el = min(n1, n0) >= MIN_ARM and mp >= MIN_PAIR
            rows.append(dict(gene=g, project=p, fga_tertile=int(f_), prolif_tertile=int(q_), n_loss=n1, n_neutral=n0, min_pair_complete=mp, eligible=el))
            if el:
                cells.append(dict(key=f"{p}|{f_}|{q_}", idx=gg.pos.to_numpy(), lab=gg.arm.to_numpy()))
        if cells:
            blocks[p] = (Y0[p][:, oi], M, cells)
    return pd.DataFrame(rows), blocks


def trans_blocks(g, blocks):
    wf = WF.subset(OTH[g]); return [(Y, M, c, wf) for Y, M, c in blocks.values()]


def null_trans(g, nperm, seed):
    el, bl = trans_cells(g); b = trans_blocks(g, bl)
    return g, (A.block_null(b, nperm, seed) if b else np.full(nperm, np.nan))


# =============================================================================================================================================
if args.stage == "power":
    sm = pd.read_csv(OUT / "sim_summary.csv"); gp = sm[sm.level == "pooled"].drop_duplicates("estimator").set_index("estimator").gate_pass
    est = "masked" if gp["masked"] == "yes" else ("adjusted" if gp["adjusted"] == "yes" else None); say("gate:", gp.to_dict(), "-> H225.1 estimator:", est)
    rows, nulls = [], {}
    nc = lambda b: (sum(len(x[2]) for x in b), sum(len(c["lab"]) for x in b for c in x[2]))
    for e in A.EST:                                          # H225.1: gate-passing estimator is primary, the others are reference rows
        b = h1_blocks(e); nl = A.block_null(b, args.nperm, np.random.SeedSequence([SEED, 1, A.EST.index(e)])); nulls[f"H225.1_{e}"] = nl
        sd = np.nanstd(nl, ddof=1); ncell, ntum = nc(b)
        rows.append(dict(hypothesis="H225.1", role="primary" if e == est else "reference", estimator=e, gene="", n_cells=ncell, n_tumours=ntum, n_perm=args.nperm, n_nan_draws=int(np.isnan(nl).sum()),
                         null_mean=np.nanmean(nl), null_sd=sd, MDE_80pct=K * sd, limit=MDE_H1, scored=(e == est) and K * sd <= MDE_H1))
        say("H225.1", e, round(np.nanmean(nl), 5), round(sd, 5))
    res = Parallel(n_jobs=args.jobs)(delayed(null_trans)(g, args.nperm, np.random.SeedSequence([SEED, 2, i])) for i, g in enumerate(GENES))
    elig = []
    for name, c in COH.items():
        Mc = (np.nan_to_num(c["CN"]) == 0).astype(float)
        for cl in c["cells"]:
            mp = A.min_pair_n(Mc, cl["idx"], cl["lab"])
            elig.append(dict(hypothesis="H225.1", cohort=name, gene="", project=name, cell=cl["key"], n_loss=int((cl["lab"] == 1).sum()), n_neutral=int((cl["lab"] == 0).sum()),
                             min_pair_complete=mp, eligible=mp >= MIN_PAIR, note="n_loss = LumA arm, n_neutral = others arm; eligible = usable by the masked estimator (raw / adjusted use all cells)"))
    for g, nl in res:
        el, bl = trans_cells(g); elig += [dict(hypothesis="H225.2" if g == "ARNTL" else "H225.3", **r) for r in el.to_dict("records")]; b = trans_blocks(g, bl); nulls[f"trans_{g}"] = nl
        sd = np.nanstd(nl, ddof=1) if b else np.nan; ncell, ntum = nc(b) if b else (0, 0)
        h = "H225.2" if g == "ARNTL" else "H225.3"
        rows.append(dict(hypothesis=h, role="primary" if g == "ARNTL" else "secondary", estimator="masked-other-13", gene=g, n_cells=ncell, n_tumours=ntum, n_perm=args.nperm, n_nan_draws=int(np.isnan(nl).sum()),
                         null_mean=np.nanmean(nl) if b else np.nan, null_sd=sd, MDE_80pct=K * sd, limit=MDE_H2,
                         scored=(K * sd <= MDE_H2) if g == "ARNTL" else (ncell >= 5)))
    pw = pd.DataFrame(rows); pw["mde_le_0.15"] = pw.MDE_80pct <= MDE_H2; pw.to_csv(OUT / "power.csv", index=False)
    pd.DataFrame(nulls).to_csv(OUT / "null_draws.csv", index=False); pd.DataFrame(elig).to_csv(OUT / "eligibility.csv", index=False)
    print(pw.round(4).to_string(), flush=True); say("done power"); sys.exit(0)


# =============================================================================================================================================
# ---- tests (NOT run before the dated power note) -------------------------------------------------------------------------------------------------
pw = pd.read_csv(OUT / "power.csv"); nl = pd.read_csv(OUT / "null_draws.csv"); assert len(nl) == args.nperm
est = pw[(pw.hypothesis == "H225.1") & (pw.role == "primary")].estimator
est = est.iloc[0] if len(est) else None
if est == "adjusted":
    print("adjusted passed the gate and masked did not: H225.1 equals A224 H224.1 (A224 output/tests.csv); not recomputed here"); est = None
pg = lambda null, o: cs.stats.perm_p(null, o, "greater", tie=0.0)
pl = lambda null, o: cs.stats.perm_p(null, o, "less", tie=0.0)
rows, pstat = [], {}

# H225.1 (masked): observed D, p one-sided greater, 200 pseudo-clock sets (A219 sets, masked by their own CN) --------------------------------------------
if est == "masked" and bool(pw[(pw.hypothesis == "H225.1") & (pw.role == "primary")].scored.iloc[0]):
    obs = float(A.block_obs(h1_blocks("masked")).mean()); psv = []
    for name, c in COH.items():
        sets, gy, gcn = pseudo_source(name); v = []
        for s in sets:
            Yp = np.column_stack([gy(g) for g in s]); Mp = (np.nan_to_num(np.column_stack([gcn(g) for g in s])) == 0).astype(float)
            v.append(A.cell_stats(Yp, Mp, cells_for(name, "masked"), L.wf_pseudo(WF, s), 0, None)[0])
        psv.append(np.array(v)); say("H225.1 pseudo", name)
    pstat["H225.1"] = np.nanmean(np.concatenate(psv, axis=1), axis=1); n1 = nl["H225.1_masked"].to_numpy()
    rows.append(dict(hypothesis="H225.1", gene="", stat=obs, null_mean=np.nanmean(n1), null_sd=np.nanstd(n1, ddof=1), p=pg(n1, obs), alt="greater"))

# H225.2 / H225.3: observed D_g and pseudo companions on the same arms -----------------------------------------------------------------------------
SCORED = set(pw[pw.hypothesis.isin(["H225.2", "H225.3"]) & pw.scored].gene)          # H225.2 needs MDE <= 0.15, H225.3 >= 5 eligible cells
TC = {g: trans_cells(g) for g in GENES if g in SCORED}
LC, PSETS = {}, {}
for p in PROJECTS:
    if args.smoke or not any(p in TC[g][1] for g in TC):
        continue
    LC[p] = load_lc(p); assert (LC[p].index.to_numpy() == tc[p].patient.to_numpy()).all()
for gi, g in enumerate(GENES):
    for pi, p in enumerate(PROJECTS):
        if g in TC and p in TC[g][1]:
            PSETS[(g, p)] = ([OTH[g]] * args.npseudo if args.smoke else
                             L.pseudo_sets(LC[p], L.A166.pool_genes(LC[p]), OTH[g], args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 3, gi, pi]))))
allg = sorted({x for v in PSETS.values() for s in v for x in s})
CNp = CNt if args.smoke else ns["cn_tcga"](allg); say("pseudo sets generated for", len(PSETS), "gene x project units;", len(allg), "pseudo genes")
for g in GENES:
    if g not in TC:
        continue
    bl = TC[g][1]; b = trans_blocks(g, bl); k = f"trans_{g}"
    if not b:
        continue
    cellsum = np.zeros(args.npseudo); ncell = 0   # pseudo sets with a non-estimable cell (pair N < 3) become NaN and are left out of the quantiles
    obs = float(A.block_obs(b).mean())
    for p, (Y, M, cells) in bl.items():
        d = tc[p]; dep = d.depth.to_numpy(float); samples = d["sample"].to_numpy()
        for si in range(args.npseudo):
            s = PSETS[(g, p)][si]
            Yp = np.random.default_rng([SEED, si]).normal(size=(len(d), 13)) if args.smoke else np.log((2 ** LC[p][list(s)].to_numpy(float) - 1) / 1e6 + 0.5 / dep[:, None])
            Mp = np.ones_like(Yp) if args.smoke else (np.nan_to_num(np.column_stack([cn_of(CNp, x, samples) for x in s])) == 0).astype(float)
            cellsum[si] += A.cell_stats(Yp, Mp, cells, L.wf_pseudo(WF.subset(OTH[g]), s), 0, None)[0].sum()
        ncell += len(cells)
    pstat[k] = cellsum / ncell; n1 = nl[k].to_numpy()
    rows.append(dict(hypothesis="H225.2" if g == "ARNTL" else "H225.3", gene=g, stat=obs, null_mean=np.nanmean(n1), null_sd=np.nanstd(n1, ddof=1), p=pl(n1, obs), alt="less"))
    say(g, "done")
pd.DataFrame(pstat).to_csv(OUT / "pseudo_stat.csv", index=False)
res = pd.DataFrame(rows)
for i, r in res.iterrows():
    k = "H225.1" if r.hypothesis == "H225.1" else f"trans_{r.gene}"; ps = pstat[k]
    ps = ps[np.isfinite(ps)]; res.loc[i, "n_pseudo_finite"] = len(ps)
    res.loc[i, "pseudo_percentile"] = 100 * np.mean(ps <= r.stat); res.loc[i, "pseudo_q"] = np.quantile(ps, 0.95 if r.hypothesis == "H225.1" else 0.05)
    res.loc[i, "specific"] = bool(r.stat > res.loc[i, "pseudo_q"]) if r.hypothesis == "H225.1" else bool(r.stat < res.loc[i, "pseudo_q"])
holm = lambda ps: np.minimum(1, np.maximum.accumulate(np.sort(ps) * (len(ps) - np.arange(len(ps)))))[np.argsort(np.argsort(ps))]
res["p_holm"] = np.nan
m12 = res.hypothesis.isin(["H225.1", "H225.2"]); m3 = res.hypothesis == "H225.3"
for m in (m12, m3):
    if m.any():
        res.loc[m, "p_holm"] = holm(res.p[m].to_numpy())
res["verdict"] = np.where(res.p_holm < 0.05, np.where(res.specific, "pass", "significant, not clock-specific"), "fail")
res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True); say("done tests")
