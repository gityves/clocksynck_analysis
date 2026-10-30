"""A169: does tumour clock coherence fall progressively with the duration of letrozole? H169.1 (GSE111563, > 120 d vs pre), H169.2
(GSE59515, 3 months vs pre); descriptive progression on patients with all three time points; A168 definitions and library.
  --stage power : permutation-null SDs and MDE = 2.49 x SD (power.csv, null_draws.csv); no observed statistic is computed.
  --stage tests : observed T, one-sided p, Holm over the two, pseudo-clock companion, progression, sensitivities.
    python -I run_a169.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 2] [--quick]"""
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
from scipy import stats  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a168_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=2); ap.add_argument("--quick", action="store_true", help="50 permutations, 10 pseudo sets, output/quick/")
args = ap.parse_args()
if args.quick:
    args.nperm, args.npseudo = 50, 10
DATA = HERE.parent / "output"
OUT = DATA / "quick" if args.quick else DATA
OUT.mkdir(parents=True, exist_ok=True)
K, CHUNK, SEED = norm.ppf(0.95) + norm.ppf(0.80), 250, zlib.crc32(b"A169") % (2 ** 31)
METABRIC, SCANB = -0.165, -0.315
SERIES = {"gse111563": "GSE111563", "gse59515": "GSE59515", "gse55374": "GSE55374"}
TEST_SERIES = {"H169.1": "gse111563", "H169.2": "gse59515"}
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


class Trend:
    """Patients with all three time points (rows pos[:, 0..2]); T = least-squares slope of coherence on the time rank 0, 1, 2;
    null: the three labels permuted within each patient."""
    def __init__(self, name, pos):
        self.name, self.pos = name, np.asarray(pos)

    @property
    def n(self):
        return len(self.pos)

    def draw(self, rng):
        return np.argsort(rng.random(self.pos.shape), axis=1)

    def T(self, Y, wf, cohf=L.coh, lab=None):
        P = self.pos if lab is None else np.take_along_axis(self.pos, lab, axis=1)
        c = np.array([cohf(Y, P[:, k], wf) for k in range(3)])
        return float(np.polyfit([0, 1, 2], c, 1)[0])

    def null(self, Y, wf, n, rng, cohf=L.coh):
        return np.array([self.T(Y, wf, cohf, self.draw(rng)) for _ in range(n)])

    def counts(self):
        return dict(n_patients=self.n)


def triplets(m):
    w = m.pivot_table(index="patient", columns="tp", values="pos", aggfunc="first").dropna()
    return w[[0, 1, 2]].to_numpy(int)


M, XX, YY, WFS, DES = {}, {}, {}, {}, {}
for s in SERIES:
    M[s], XX[s] = L.load(s, DATA); cg = L.clock_genes(XX[s])
    YY[s], WFS[s] = L.Yof(XX[s], cg), L.wf_of(cg)
for n, s in TEST_SERIES.items():
    x, y, _ = L.pairs_from_meta(M[s], 0, 2); DES[n] = L.Pairs(n, "swap", x, y)
for n, d in DES.items():
    say(n, SERIES[TEST_SERIES[n]], "genes", YY[TEST_SERIES[n]].shape[1], d.counts())


def chunked_null(k, d, Y, wf, n, cohf=L.coh):
    sizes = [CHUNK] * (n // CHUNK) + ([n % CHUNK] if n % CHUNK else [])
    rngs = [np.random.default_rng(np.random.SeedSequence([SEED, k, c])) for c in range(len(sizes))]
    return np.concatenate([d.null(Y, wf, m, r, cohf) for m, r in zip(sizes, rngs)])


def job(k, d, s, drop=(), resid=False):
    cg = [g for g in L.clock_genes(XX[s]) if g not in drop]
    Y, wf = L.Yof(XX[s], cg), L.wf_of(cg)
    cohf = L.residual_coh(M[s][["prolif", "purity"]].to_numpy()) if resid else L.coh
    return chunked_null(k, d, Y, wf, args.nperm, cohf), Y, wf, cohf, len(cg)


if args.stage == "power":
    res = Parallel(n_jobs=args.jobs)(delayed(job)(i + 1, DES[n], TEST_SERIES[n]) for i, n in enumerate(DES))
    nulls = {n: r[0] for n, r in zip(DES, res)}
    pd.DataFrame(nulls).to_csv(OUT / "null_draws.csv", index=False)
    pw = pd.DataFrame([dict(hypothesis=n, series=SERIES[TEST_SERIES[n]], n_patients=DES[n].n, n_genes=YY[TEST_SERIES[n]].shape[1],
                            null_mean=nulls[n].mean(), null_sd=nulls[n].std(ddof=1), MDE_80pct=K * nulls[n].std(ddof=1),
                            metabric_effect=METABRIC, scanb_effect=SCANB, n_perm=args.nperm) for n in DES])
    pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests --------------------------------------------------------------------------------------------------------------------------
nulls = pd.read_csv(OUT / "null_draws.csv")
assert len(nulls) == args.nperm
names = list(DES)
obs = {n: DES[n].T(YY[TEST_SERIES[n]], WFS[TEST_SERIES[n]]) for n in names}
p = np.array([pl(nulls[n].to_numpy(), obs[n]) for n in names])
ph = np.minimum(1, np.maximum.accumulate(np.sort(p) * (len(p) - np.arange(len(p)))))[np.argsort(np.argsort(p))]
est, rows = L.estrogen_genes(), []
for j, n in enumerate(names):
    s = TEST_SERIES[n]; X, cg = XX[s], L.clock_genes(XX[s])
    sets = L.pseudo_sets(X, L.pool_genes(X, est), cg, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 10 + j])))
    pd.DataFrame(sets).to_csv(OUT / f"pseudo_sets_{s}.csv", index=False)
    pT = np.array([DES[n].T(L.Yof(X, g), L.wf_pseudo(WFS[s], g)) for g in sets])
    T, q05, sd = obs[n], np.quantile(pT, 0.05), nulls[n].std(ddof=1)
    spec = ("clock-specific" if T < q05 else "coherence loss not specific to the clock") if ph[j] < 0.05 else "not applicable (test failed)"
    rows.append(dict(hypothesis=n, series=SERIES[s], T=T, p_one_lower=p[j], p_holm=ph[j], null_mean=nulls[n].mean(), null_sd=sd,
                     MDE_80pct=K * sd, n_patients=DES[n].n, n_genes=len(cg), pseudo_percentile=100 * np.mean(pT <= T), pseudo_q05=q05,
                     verdict="pass" if ph[j] < 0.05 else "fail", specificity=spec, n_perm=args.nperm, n_pseudo=args.npseudo))
res = pd.DataFrame(rows); res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True)

# ---- descriptive: coherence per time point (all samples and the complete triplets), progression slope, treatment-effect scores -------
cg_rows, prog, ch = [], [], []
for i, s in enumerate(SERIES):
    m, Y, wf = M[s], YY[s], WFS[s]
    for t, g in m.groupby("tp"):
        cg_rows.append(dict(series=SERIES[s], set="all samples", tp=t, n=len(g), coherence=L.coh(Y, g.pos.to_numpy(), wf)))
    tri = Trend("progression", triplets(m))
    for t in range(3):
        cg_rows.append(dict(series=SERIES[s], set="complete triplets", tp=t, n=tri.n, coherence=L.coh(Y, tri.pos[:, t], wf)))
    T = tri.T(Y, wf); nl = chunked_null(50 + i, tri, Y, wf, args.nperm)
    prog.append(dict(series=SERIES[s], n_patients=tri.n, slope_per_step=T, p_one_lower=pl(nl, T), null_sd=nl.std(ddof=1), MDE_80pct=K * nl.std(ddof=1)))
    for t in (1, 2):
        a, b, _ = L.pairs_from_meta(m, 0, t)
        r = dict(series=SERIES[s], contrast=f"pre -> tp {t}", n_pairs=len(a))
        for c in ("prolif", "purity", "estro"):
            d = m[c].to_numpy()[b] - m[c].to_numpy()[a]
            r.update({f"d_{c}": d.mean(), f"sd_{c}": d.std(ddof=1), f"p_{c}": stats.ttest_1samp(d, 0).pvalue})
        ch.append(r)
for f, t in (("coherence_tp.csv", cg_rows), ("progression.csv", prog), ("treatment_effect_scores.csv", ch)):
    pd.DataFrame(t).to_csv(OUT / f, index=False); print(pd.DataFrame(t).round(4).to_string(), flush=True)

# ---- sensitivities --------------------------------------------------------------------------------------------------------------------
specs = []
for n in names:
    s = TEST_SERIES[n]
    specs += [(f"{n} without PER1", DES[n], s, ["PER1"], False), (f"{n} without PER1 and NR1D1", DES[n], s, ["PER1", "NR1D1"], False),
              (f"{n} Y residualised on proliferation and purity proxy within time point", DES[n], s, [], True)]
m = M["gse59515"]
for r in ("Responder", "Non-responder"):
    x, y, _ = L.pairs_from_meta(m[m.response == r], 0, 2); specs.append((f"H169.2 {r}s", L.Pairs("H169.2", "swap", x, y), "gse59515", [], False))
x, y, _ = L.pairs_from_meta(M["gse55374"], 0, 2); specs.append(("GSE55374 3 months vs pre", L.Pairs("GSE55374", "swap", x, y), "gse55374", [], False))


def sens(j, lab, d, s, drop, resid):
    nl, Y, wf, cohf, ng = job(100 + j, d, s, drop, resid)
    T = d.T(Y, wf, cohf)
    return dict(analysis=lab, series=SERIES[s], T=T, p_one_lower=pl(nl, T), null_sd=nl.std(ddof=1), MDE_80pct=K * nl.std(ddof=1), n_genes=ng,
                **d.counts())


sn = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(sens)(j, *sp) for j, sp in enumerate(specs)))
sn.to_csv(OUT / "sensitivity.csv", index=False); print(sn.round(4).to_string(), flush=True)
say("done tests")
