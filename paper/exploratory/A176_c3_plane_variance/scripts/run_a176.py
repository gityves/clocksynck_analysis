"""A176: is the coherence loss a swamped clock (clock-plane variance down relative to off-plane) or a disordered one?
Seven registered coherence contrasts (A127 H127.1/H127.2, A166 H166.1, A171 H171.1/2, A173 H173.1/3), low-coherence arm minus
reference arm, mean over the registered strata; statistics from clocksync.decompose.components on the same arm-centred covariance
as the registered coherence: dR, dlog_plane, dlog_off, ds = dlog_plane - dlog_off. One label permutation (registered scheme) gives all
four statistics.
  --stage calibration : reproduces dR with the registered coherence code (asserted to 1e-6 against the amendment tables), the
                        permutation nulls of all statistics, the pseudo-clock sets (registered 200 sets per cohort, observed labels),
                        regression dR ~ ds and the MDE (2.49 x null SD). The observed dlog_plane / dlog_off / ds of the clock genes are
                        never computed in this stage.
  --stage tests       : observed statistics, one-sided p of ds < 0, Holm over seven, pseudo percentile, H176.2 calls, descriptives.
    python run_a176.py --stage calibration|tests [--nperm 2000] [--jobs 4]"""
import argparse
import os
import sys
import zlib
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
sys.path.insert(0, str(EXP / "A168_endocrine/scripts"))
sys.path.insert(0, str(EXP / "A127_scanb/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a127_lib as A127  # noqa: E402
import a168_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.decompose import components  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["calibration", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--jobs", type=int, default=4)
args = argparse.Namespace(stage=None, nperm=2000, jobs=4)  # set by main(); importers (A178, A179) may set args.jobs
OUT = HERE.parent / "output"; OUT.mkdir(exist_ok=True)
K, SEED = norm.ppf(0.95) + norm.ppf(0.80), zlib.crc32(b"A176") % (2 ** 31)
A166, tertile = L.A166, L.A166.tertile
STATS = ["dR", "dlp", "dlo", "ds"]
NAMES = ["H127.1", "H127.2", "H166.1", "H171.1", "H171.2", "H173.1", "H173.3"]
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
EXPR = lambda d: EXP / d / "output"


class Contrast:
    """Groups (index arrays); labels per group; case arm = label `case`, reference arm = label `ref`. draw(rng) returns permuted
    labels (registered scheme), obs = observed labels."""
    def __init__(self, name, G, obs, draw, case, ref):
        self.name, self.G, self.obs, self.draw, self.case, self.ref = name, G, obs, draw, case, ref

    def parts(self, Y, wf, labs=None):
        """Per-group (dR, dlp, dlo) via components (coherence, log_plane, log_off)."""
        out = []
        for i, l in zip(self.G, self.obs if labs is None else labs):
            c = []
            for v in (self.case, self.ref):
                y = Y[i[l == v]]; y = y - y.mean(0)
                c.append(components(y, np.zeros(y.shape[1]), 1, wf))
            out.append([c[0]["coherence"] - c[1]["coherence"], c[0]["log_plane"] - c[1]["log_plane"], c[0]["log_off"] - c[1]["log_off"]])
        return out

    def stats(self, Y, wf, labs=None):
        d = np.mean(self.parts(Y, wf, labs), 0)
        return np.array([d[0], d[1], d[2], d[1] - d[2]])

    def dR_registered(self, Y, wf):
        """Registered coherence statistic (a127_lib.coh)."""
        return float(np.mean([A127.coh(Y, i[l == self.case], wf) - A127.coh(Y, i[l == self.ref], wf) for i, l in zip(self.G, self.obs)]))


def lab_contrast(name, G, lab, case, ref):
    lab = np.asarray(lab, float)
    G = list(G)
    return Contrast(name, G, [lab[i] for i in G], lambda rng: [lab[i][rng.permutation(len(i))] for i in G], case, ref)


def prolif_contrast(name, strata, prolif):
    """Proliferation tertile labels (registered: tertile of proliferation within stratum); null permutes proliferation within stratum."""
    p = np.asarray(prolif, float)
    strata = list(strata)

    def draw(rng):
        return [tertile(rng.permutation(p[i])) for i in strata]
    return Contrast(name, strata, [tertile(p[i]) for i in strata], draw, 2, 0)


# ---------------------------------------------------------------------------------------------------------------------------------
def zscore_scores(X, m, by=None):  # A171 definition
    m = m.copy()
    groups = m.groupby(by).indices.values() if by else [np.arange(len(m))]
    for i in groups:
        for col, genes in (("prolif", L.sh.hallmark_g2m()), ("purity", L.sh.PURITY)):
            s = L.sh._zmean(X.iloc[i], [g for g in genes if g in X.columns])
            m.iloc[i, m.columns.get_loc(col)] = ((s - s.mean()) / s.std()).to_numpy()
    return m


def purity_strata(purity, min_n=30):  # A171 / A173 Prolif
    pt = tertile(purity)
    return [np.flatnonzero(pt == k) for k in range(3) if (pt == k).sum() >= min_n]


def build_all():
    """dict name -> dict(contrast, X, m, cg, Y, wf, genes csv of the registered pseudo sets, site-centring labels)."""
    R = {}
    # SCAN-B ------------------------------------------------------------------------------------------------------------------
    m, X = A127.load(); LN2 = L.LN2
    m["pur_t"] = m.groupby("pam50").purity.transform(lambda s: A127.tertile(s))
    m["pur_t_cohort"] = A127.tertile(m.purity)
    S1 = []
    for _, g in m.groupby(["pam50", "pur_t"]):
        if len(g) >= 30:
            S1.append(np.array([m.index.get_loc(i) for i in g.index]))
    prolif = m.prolif.to_numpy()
    lab1 = [A127.tertile(prolif[i]) for i in S1]
    ptert = [np.flatnonzero(m.pur_t_cohort.to_numpy() == t) for t in range(3)]
    er = m.er.to_numpy()
    cg = list(A127.GENES)
    d1 = Contrast("H127.1", S1, lab1, lambda rng: [l[rng.permutation(len(l))] for l in lab1], 2, 0)
    d2 = lab_contrast("H127.2", ptert, er, 0, 1)
    for d in (d1, d2):
        R[d.name] = dict(c=d, X=X, m=m, cg=cg, sets=EXPR("A127_scanb") / "pseudo_sets_variance_matched.csv", center=None)
    # METABRIC ----------------------------------------------------------------------------------------------------------------
    mb, Xb = A166.load("metabric"); mb["pos"] = np.arange(len(mb))
    h = A166.build(mb)[0]
    R["H166.1"] = dict(c=lab_contrast("H166.1", h.G, h.er, 0, 1), X=Xb, m=mb, cg=list(A166.GENES),
                       sets=EXPR("A166_metabric") / "pseudo_sets_metabric.csv", center=None)
    # lung ----------------------------------------------------------------------------------------------------------------------
    for n, key in (("H171.1", "gse72094"), ("H171.2", "gse68465")):
        mm, XX = L.load(key, EXPR("A171_luad_lihc"))
        by = "site" if key == "gse68465" else None
        mm = zscore_scores(XX, mm, by).assign(pos=np.arange(len(mm)))
        R[n] = dict(c=prolif_contrast(n, purity_strata(mm.purity), mm.prolif), X=XX, m=mm, cg=L.clock_genes(XX),
                    sets=EXPR("A171_luad_lihc") / f"pseudo_sets_{key}.csv", center=mm[by].to_numpy() if by else None)
    for n, key in (("H173.1", "gse62254"), ("H173.3", "gse15459")):
        mm, XX = L.load(key, EXPR("A173_gastric_lung_histology"))
        R[n] = dict(c=prolif_contrast(n, purity_strata(mm.purity), mm.prolif), X=XX, m=mm, cg=L.clock_genes(XX),
                    sets=EXPR("A173_gastric_lung_histology") / f"pseudo_sets_{key}.csv", center=None)
    return R


def Yof(r, genes):
    Y = L.Yof(r["X"], genes)
    if r["center"] is not None:  # GSE68465: Y centred per gene within site (registered)
        Y = Y - pd.DataFrame(Y).groupby(r["center"]).transform("mean").to_numpy()
    return Y


def wf_for(genes, base):
    return L.wf_pseudo(base, genes)


def registered_dR():
    t = {}
    a = pd.read_csv(EXPR("A127_scanb") / "a127_tests.csv").set_index("test")["T"]
    t["H127.1"] = a[[k for k in a.index if k.startswith("H127.1")][0]]; t["H127.2"] = a[[k for k in a.index if k.startswith("H127.2")][0]]
    t["H166.1"] = pd.read_csv(EXPR("A166_metabric") / "a166_tests.csv").set_index("hypothesis")["T"]["H166.1"]
    b = pd.read_csv(EXPR("A171_luad_lihc") / "tests.csv").set_index("hypothesis")["T"]
    t["H171.1"], t["H171.2"] = b["H171.1"], b["H171.2"]
    c = pd.read_csv(EXPR("A173_gastric_lung_histology") / "tests.csv").set_index("hypothesis")["stat"]
    t["H173.1"], t["H173.3"] = c["H173.1"], c["H173.3"]
    return t


def null_chunk(c, Y, wf, seed, n, tag=""):
    rng = np.random.default_rng(seed)
    return np.array([c.stats(Y, wf, c.draw(rng)) for _ in range(n)])


def null_of(c, Y, wf, j, n_total):
    """Permutation null (n_total x 4); chunks of 500 with SeedSequence([SEED, j]).spawn."""
    nch = int(np.ceil(n_total / 500)); seeds = np.random.SeedSequence([SEED, j]).spawn(nch)
    sizes = [min(500, n_total - 500 * k) for k in range(nch)]
    return np.vstack(Parallel(n_jobs=args.jobs)(delayed(null_chunk)(c, Y, wf, s, n) for s, n in zip(seeds, sizes)))


def pseudo_stats(c, r, sets):
    return np.array([c.stats(Yof(r, g), wf_for(g, r["wf"])) for g in sets])


def read_sets(path, n_genes):
    s = pd.read_csv(path)
    return [[g for g in row if isinstance(g, str)] for row in s.values.tolist()]


def main():
    global args
    args = ap.parse_args()
    say("stage", args.stage)
    R = build_all()
    for n, r in R.items():
        r["wf"] = L.wf_of(r["cg"]); r["Y"] = Yof(r, r["cg"])
        r["sets_list"] = read_sets(r["sets"], len(r["cg"]))
        assert all(len(s) == len(r["cg"]) for s in r["sets_list"]), n
        say(n, "groups", [len(i) for i in r["c"].G], "genes", len(r["cg"]), "pseudo sets", len(r["sets_list"]))

    # ---- reproduction of the registered coherence statistic --------------------------------------------------------------------------
    reg = registered_dR()
    rep = []
    for n in NAMES:
        r = R[n]; v = r["c"].dR_registered(r["Y"], r["wf"])
        assert abs(v - reg[n]) < 1e-6, (n, v, reg[n])
        rep.append(dict(contrast=n, dR_reproduced=v, dR_registered_table=reg[n], abs_diff=abs(v - reg[n])))
    rep = pd.DataFrame(rep); say("dR reproduced to 1e-6 for all seven"); print(rep.to_string(), flush=True)

    if args.stage == "calibration":
        rep.to_csv(OUT / "reproduction.csv", index=False)
        nulls, prow, cal = {}, [], []
        for j, n in enumerate(NAMES):
            r = R[n]; c = r["c"]
            nl = null_of(c, r["Y"], r["wf"], j, args.nperm); nulls[n] = nl
            say(n, "null done")
            ps_ = pseudo_stats(c, r, r["sets_list"])
            for k, v in enumerate(ps_):
                prow.append(dict(contrast=n, set=k, **dict(zip(STATS, v))))
            sl, ic = np.polyfit(ps_[:, 3], ps_[:, 0], 1)
            slo = (ps_[:, 3] @ ps_[:, 0]) / (ps_[:, 3] @ ps_[:, 3])
            sd = nl[:, 3].std(ddof=1)
            cal.append(dict(contrast=n, dR_registered=reg[n], n_pseudo=len(ps_), slope_dR_on_ds=sl, intercept=ic,
                            r_dR_ds=np.corrcoef(ps_[:, 3], ps_[:, 0])[0, 1], pseudo_ds_sd=ps_[:, 3].std(ddof=1), pseudo_dR_sd=ps_[:, 0].std(ddof=1),
                            ds_swamping_needed=(reg[n] - ic) / sl, ds_swamping_needed_through_origin=reg[n] / slo,
                            null_sd_ds=sd, MDE_ds=K * sd, null_sd_dR=nl[:, 0].std(ddof=1), null_sd_dlp=nl[:, 1].std(ddof=1),
                            null_sd_dlo=nl[:, 2].std(ddof=1), n_perm=args.nperm))
        pd.DataFrame({f"{n}|{s}": nulls[n][:, k] for n in NAMES for k, s in enumerate(STATS)}).to_csv(OUT / "null_draws.csv", index=False)
        pd.DataFrame(prow).to_csv(OUT / "pseudo_stats.csv", index=False)
        cal = pd.DataFrame(cal); cal.to_csv(OUT / "calibration.csv", index=False); print(cal.round(4).to_string(), flush=True)
        say("done calibration"); sys.exit(0)

    # ---- tests ---------------------------------------------------------------------------------------------------------------------
    nl = pd.read_csv(OUT / "null_draws.csv"); assert len(nl) == args.nperm
    ps_all = pd.read_csv(OUT / "pseudo_stats.csv"); cal = pd.read_csv(OUT / "calibration.csv").set_index("contrast")
    obs = {n: R[n]["c"].stats(R[n]["Y"], R[n]["wf"]) for n in NAMES}
    for n in NAMES:
        assert abs(obs[n][0] - reg[n]) < 1e-6, n
    p = np.array([pl(nl[f"{n}|ds"].to_numpy(), obs[n][3]) for n in NAMES])
    ph = np.minimum(1, np.maximum.accumulate(np.sort(p) * (len(p) - np.arange(len(p)))))[np.argsort(np.argsort(p))]
    rows = []
    for j, n in enumerate(NAMES):
        o = obs[n]; ps_ = ps_all[ps_all.contrast == n]; cr = cal.loc[n]
        h1 = bool(ph[j] < 0.05)
        need = abs(cr.ds_swamping_needed)
        h2 = "swamped (H176.1 passes)" if h1 else ("disordered" if cr.MDE_ds <= need else "unresolved")
        q05 = np.quantile(ps_.ds, 0.05)
        spec = ("clock-gene specific" if o[3] < q05 else "general feature of expression (not below pseudo 5th percentile)") if o[3] < 0 else "not applicable (ds >= 0)"
        rows.append(dict(contrast=n, dR=o[0], dlog_plane=o[1], dlog_off=o[2], ds=o[3], null_mean_ds=nl[f"{n}|ds"].mean(), null_sd_ds=cr.null_sd_ds,
                         p_ds_less=p[j], p_holm=ph[j], H176_1_pass=h1, MDE_ds=cr.MDE_ds, ds_swamping_needed=cr.ds_swamping_needed,
                         abs_needed=need, H176_2=h2, pseudo_percentile_ds=100 * np.mean(ps_.ds.to_numpy() <= o[3]), pseudo_q05_ds=q05,
                         pseudo_median_ds=np.median(ps_.ds), specificity=spec,
                         p_dlp_two_sided=cs.stats.perm_p(nl[f"{n}|dlp"].to_numpy(), o[1], "two-sided", tie=0.0),
                         p_dlo_two_sided=cs.stats.perm_p(nl[f"{n}|dlo"].to_numpy(), o[2], "two-sided", tie=0.0),
                         pseudo_percentile_dR=100 * np.mean(ps_.dR.to_numpy() <= o[0]),
                         pseudo_percentile_dlp=100 * np.mean(ps_.dlp.to_numpy() <= o[1]), pseudo_percentile_dlo=100 * np.mean(ps_.dlo.to_numpy() <= o[2]),
                         n_perm=args.nperm))
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "tests.csv", index=False)
    n_pass, n_dis = int(res.H176_1_pass.sum()), int((res.H176_2 == "disordered").sum())
    summ = dict(n_H176_1_pass=n_pass, H176_1_supported=n_pass >= 5, n_disordered=n_dis, disordered_ge5=n_dis >= 5,
                n_unresolved=int((res.H176_2 == "unresolved").sum()))
    pd.Series(summ).to_csv(OUT / "verdict.csv", header=False)
    print(res.round(4).T.to_string(), flush=True); print(summ, flush=True)

    # ---- descriptives: without PER1 and NR1D1 ----------------------------------------------------------------------------------------
    dd = []
    for j, n in enumerate(NAMES):
        r = R[n]; c = r["c"]; cg = r["cg"]
        keep = [g for g in cg if g not in ("PER1", "NR1D1")]
        Yk = r["Y"][:, [cg.index(g) for g in keep]]; wk = L.wf_of(keep)
        o = c.stats(Yk, wk); x = null_of(c, Yk, wk, 100 + j, args.nperm)
        dd.append(dict(contrast=n, n_genes=len(keep), dR=o[0], dlog_plane=o[1], dlog_off=o[2], ds=o[3], null_sd_ds=x[:, 3].std(ddof=1),
                       p_ds_less=pl(x[:, 3], o[3]), full_ds=obs[n][3]))
        say(n, "without PER1/NR1D1 done")
    pd.DataFrame(dd).to_csv(OUT / "descriptive_without_per1_nr1d1.csv", index=False)
    print(pd.DataFrame(dd).round(4).to_string(), flush=True)

    # ---- descriptives: proliferation-responsive pseudo pool (C4 pool) ----------------------------------------------------------------
    cc = []
    for j, n in enumerate(NAMES):
        r = R[n]; c = r["c"]; X = r["X"]
        pool = A166.pool_genes(X)
        cor = X[pool].corrwith(pd.Series(r["m"].prolif.to_numpy(), index=X.index)).abs()
        top = cor[cor >= cor.quantile(0.9)].index.tolist()
        sets = L.pseudo_sets(X, top, r["cg"], 200, np.random.default_rng(np.random.SeedSequence([SEED, 200 + j])))
        pst = pseudo_stats(c, r, sets); o = obs[n]
        pd.DataFrame(sets).to_csv(OUT / f"pseudo_sets_prolif_pool_{n}.csv", index=False)
        pd.DataFrame(pst, columns=STATS).assign(contrast=n).to_csv(OUT / f"pseudo_stats_prolif_pool_{n}.csv", index=False)
        cc.append(dict(contrast=n, n_pool=len(top), n_sets=len(sets), pool_ds_median=np.median(pst[:, 3]), pool_ds_q05=np.quantile(pst[:, 3], .05),
                       pool_dR_median=np.median(pst[:, 0]), ds_percentile=100 * np.mean(pst[:, 3] <= o[3]), dR_percentile=100 * np.mean(pst[:, 0] <= o[0]),
                       dlp_percentile=100 * np.mean(pst[:, 1] <= o[1]), dlo_percentile=100 * np.mean(pst[:, 2] <= o[2]),
                       ds_below_pool_q05=bool(o[3] < np.quantile(pst[:, 3], .05))))
        say(n, "prolif pool done")
    pd.DataFrame(cc).to_csv(OUT / "descriptive_prolif_pool.csv", index=False)
    print(pd.DataFrame(cc).round(4).to_string(), flush=True)
    say("done tests")


if __name__ == "__main__":
    main()
