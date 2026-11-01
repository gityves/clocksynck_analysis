"""A182 (registered in paper/docs/registered_experiments.md): is the proliferation clock-coherence loss specific to breast, lung
adenocarcinoma and gastric cancer? H182.1 GSE39582, H182.2 GSE17536 (colon), H182.3 GSE116918 (prostate) = Family 1;
H182.4 GSE157009, H182.5 GSE157010 (lung squamous) = Family 2 (Holm within family). GSE65858 (head and neck): descriptive, two-sided.
Statistic = H171 (purity-proxy tertile strata >= 30 tumours, T = mean of coh(high G2M tertile) - coh(low); proliferation permuted
within stratum, one-sided for a loss). Code reuses a168_lib / A171 / A173 conventions. Data from prep_c6.py.
  --stage power : permutation nulls, null SD, MDE = 2.49 x SD (power.csv, null_draws.csv); no observed statistic is computed.
  --stage tests : T, p, Holm, pseudo-clock companion (200 sets), equivalence T - 1.645 x null SD > -0.19, MDE > 0.30 -> not estimable
                  (tests.csv); GSE65858 and descriptive rows (descriptive.csv, strata.csv).
    python run_c6.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 4]"""
import argparse
import os
import sys
import zlib
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a168_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=4)
args = ap.parse_args()
OUT = HERE.parent / "output"
K, SEED, MIN_STRATUM = norm.ppf(0.95) + norm.ppf(0.80), zlib.crc32(b"A182") % (2 ** 31), 30
EQ_MARGIN, EQ_Z, MDE_MAX = -0.19, norm.ppf(0.95), 0.30
FAMILY = {"H182.1": (1, "gse39582"), "H182.2": (1, "gse17536"), "H182.3": (1, "gse116918"),
          "H182.4": (2, "gse157009"), "H182.5": (2, "gse157010")}
DESC_SERIES = "gse65858"
TESTS = {n: s for n, (_, s) in FAMILY.items()}
MAIN = list(TESTS.values())
ALL = MAIN + [DESC_SERIES]
tertile = L.A166.tertile
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


def tertile_effect(Y, wf, strata, p):
    return float(np.mean([L.coh(Y, i[tertile(p[i]) == 2], wf) - L.coh(Y, i[tertile(p[i]) == 0], wf) for i in strata]))


def permute_within(p, strata, rng):
    q = p.copy()
    for i in strata:
        q[i] = rng.permutation(q[i])
    return q


class Prolif:
    """H171 statistic on purity-proxy tertile strata (>= 30 tumours); null: proliferation permuted within stratum."""
    def __init__(self, m):
        pt = tertile(m.purity)
        self.strata = [np.flatnonzero(pt == k) for k in range(3) if (pt == k).sum() >= MIN_STRATUM]
        self.prolif = m.prolif.to_numpy(float)

    def T(self, Y, wf, p=None):
        return tertile_effect(Y, wf, self.strata, self.prolif if p is None else p)

    def null(self, Y, wf, n, rng):
        return np.array([self.T(Y, wf, permute_within(self.prolif, self.strata, rng)) for _ in range(n)])


def series(key):
    m, X = L.load(key, OUT)
    cg = L.clock_genes(X)
    return dict(m=m, X=X, cg=cg, Y=L.Yof(X, cg), wf=L.wf_of(cg), d=Prolif(m))


DATA = {k: series(k) for k in ALL}
for k, s in DATA.items():
    say(k, "tumours", len(s["m"]), "genes", len(s["cg"]), "strata", [len(i) for i in s["d"].strata])


def null_of(j, key, genes=None, d=None):
    s = DATA[key]
    g = s["cg"] if genes is None else genes
    Y = s["Y"] if genes is None else L.Yof(s["X"], g)
    return (d or s["d"]).null(Y, L.wf_of(g), args.nperm, np.random.default_rng(np.random.SeedSequence([SEED, j])))


if args.stage == "power":
    res = dict(zip(ALL, Parallel(n_jobs=args.jobs)(delayed(null_of)(j + 1, k) for j, k in enumerate(ALL))))
    pd.DataFrame(res).to_csv(OUT / "null_draws.csv", index=False)
    names = {**{s: n for n, s in TESTS.items()}, DESC_SERIES: "descriptive HNSC"}
    pw = pd.DataFrame([dict(analysis=names[k], series=k.upper(), n=len(DATA[k]["m"]), n_genes=len(DATA[k]["cg"]),
                            strata=";".join(str(len(i)) for i in DATA[k]["d"].strata), null_mean=res[k].mean(),
                            null_sd=res[k].std(ddof=1), MDE_80pct=K * res[k].std(ddof=1),
                            not_estimable=K * res[k].std(ddof=1) > MDE_MAX) for k in ALL])
    pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests --------------------------------------------------------------------------------------------------------------------------
nl = pd.read_csv(OUT / "null_draws.csv")
assert len(nl) == args.nperm
obs = {k: DATA[k]["d"].T(DATA[k]["Y"], DATA[k]["wf"]) for k in ALL}


def pseudo(j, key):
    s = DATA[key]
    sets = L.pseudo_sets(s["X"], L.A166.pool_genes(s["X"]), s["cg"], args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 10 + j])))
    return sets, np.array([s["d"].T(L.Yof(s["X"], g), L.wf_pseudo(s["wf"], g)) for g in sets])


ps = dict(zip(ALL, Parallel(n_jobs=args.jobs)(delayed(pseudo)(j, k) for j, k in enumerate(ALL))))
for k in ALL:
    pd.DataFrame(ps[k][0]).to_csv(OUT / f"pseudo_sets_{k}.csv", index=False)
    pd.DataFrame({"T": ps[k][1]}).to_csv(OUT / f"pseudo_stats_{k}.csv", index=False)

p = {n: pl(nl[s].to_numpy(), obs[s]) for n, s in TESTS.items()}
holm = {}
for fam in (1, 2):
    nm = [n for n in TESTS if FAMILY[n][0] == fam]
    pv = np.array([p[n] for n in nm])
    ph = np.minimum(1, np.maximum.accumulate(np.sort(pv) * (len(pv) - np.arange(len(pv)))))[np.argsort(np.argsort(pv))]
    holm.update(dict(zip(nm, ph)))


def verdict(passed, specific, mde, equiv):
    if passed:
        return "passes, clock-specific (type joins the claim / histology wording dropped)" if specific else "passes, not clock-specific"
    if mde > MDE_MAX:
        return "not estimable (MDE > 0.30)"
    return "no loss of the replicated size (equivalent)" if equiv else "not resolved at this n"


rows = []
for n, s in TESTS.items():
    T, pT, sd = obs[s], ps[s][1], nl[s].std(ddof=1)
    mde = K * sd
    passed, specific, equiv = holm[n] < 0.05, T < np.quantile(pT, 0.05), (T - EQ_Z * sd) > EQ_MARGIN
    rows.append(dict(hypothesis=n, family=FAMILY[n][0], series=s.upper(), T=T, p_one_lower=p[n], p_holm=holm[n], null_mean=nl[s].mean(), null_sd=sd,
                     MDE_80pct=mde, n_tumours=len(DATA[s]["m"]), n_genes=len(DATA[s]["cg"]), pseudo_percentile=100 * np.mean(pT <= T),
                     pseudo_q05=np.quantile(pT, 0.05), clock_specific=bool(passed and specific), equiv_bound=T - EQ_Z * sd,
                     equivalent=bool(equiv), not_estimable=bool(mde > MDE_MAX), verdict=verdict(passed, specific, mde, equiv),
                     n_perm=args.nperm, n_pseudo=args.npseudo))
res = pd.DataFrame(rows); res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True)

# ---- descriptive --------------------------------------------------------------------------------------------------------------------
desc = []
s = DESC_SERIES; x = nl[s].to_numpy(); sd = x.std(ddof=1)
desc.append(dict(analysis="GSE65858 HNSC T (two-sided)", series=s.upper(), n=len(DATA[s]["m"]), stat=obs[s], p_two_sided=cs.stats.perm_p(x, obs[s], "two-sided", tie=0.0),
                 p_one_lower=pl(x, obs[s]), null_sd=sd, MDE_80pct=K * sd, pseudo_percentile=100 * np.mean(ps[s][1] <= obs[s]),
                 pseudo_q05=np.quantile(ps[s][1], 0.05), equiv_bound=obs[s] - EQ_Z * sd))
jobs = []
for key in ALL:
    d = DATA[key]; keep = [g for g in d["cg"] if g not in ("PER1", "NR1D1")]
    jobs.append((f"{key.upper()} without PER1 and NR1D1", key, keep, d["d"]))
m39 = DATA["gse39582"]["m"]
for ds in ("discovery", "validation"):
    i = np.flatnonzero((m39.dataset == ds).to_numpy())
    jobs.append((f"GSE39582 {ds} only (strata recomputed within the dataset)", "gse39582", None, Prolif(m39.iloc[i]), i))


def one(j, name, key, genes, d, idx=None):
    s = DATA[key]; g = s["cg"] if genes is None else genes
    Y = L.Yof(s["X"], g); wf = L.wf_of(g)
    if idx is not None:
        Y = Y[idx]
    x = d.null(Y, wf, args.nperm, np.random.default_rng(np.random.SeedSequence([SEED, 100 + j])))
    o = d.T(Y, wf)
    return dict(analysis=name, series=key.upper(), n=len(Y), n_genes=len(g), stat=o, p_one_lower=pl(x, o), null_sd=x.std(ddof=1), MDE_80pct=K * x.std(ddof=1))


for r in Parallel(n_jobs=args.jobs)(delayed(one)(j, *jb) for j, jb in enumerate(jobs)):
    desc.append(r)
dd = pd.DataFrame(desc); dd.to_csv(OUT / "descriptive.csv", index=False); print(dd.round(4).to_string(), flush=True)

st = []
for k in ALL:
    s = DATA[k]; m = s["m"]
    for q, i in enumerate(s["d"].strata):
        t = tertile(s["d"].prolif[i])
        c = {lv: L.coh(s["Y"], i[t == z], s["wf"]) for z, lv in enumerate(("low", "mid", "high"))}
        st.append(dict(series=k.upper(), stratum=q, n=len(i), **{f"coh_prolif_{lv}": v for lv, v in c.items()}, T_stratum=c["high"] - c["low"]))
pd.DataFrame(st).to_csv(OUT / "strata.csv", index=False); print(pd.DataFrame(st).round(4).to_string(), flush=True)
say("done tests")
