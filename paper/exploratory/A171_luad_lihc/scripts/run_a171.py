"""A171: does tumour clock coherence fall with proliferation in lung adenocarcinoma and hepatocellular carcinoma? (registered in
paper/docs/registered_experiments.md, A171). H171.1 GSE72094, H171.2 GSE68465 (LUAD); H171.3 GSE14520, H171.4 GSE36376 (LIHC).
Statistic = A166 H166.3: within purity-proxy tertiles (>= 30 tumours), mean of coh(high proliferation tertile) - coh(low);
proliferation permuted within stratum.
  --stage power : permutation nulls, MDE = 2.49 x SD (power.csv, null_draws.csv), and the TCGA LUAD / LIHC values of the same
                  statistic (seen data, tcga.csv). No statistic of the four test series is computed.
  --stage tests : observed T, one-sided p, Holm over four, pseudo-clock companion; descriptive rows.
    python run_a171.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 4]"""
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
K, SEED, MIN_STRATUM, MIN_ARM = norm.ppf(0.95) + norm.ppf(0.80), zlib.crc32(b"A171") % (2 ** 31), 30, 10
TESTS = {"H171.1": "gse72094", "H171.2": "gse68465", "H171.3": "gse14520", "H171.4": "gse36376"}
TUMOUR = {"gse14520": "Tumor", "gse36376": "liver tumor", "gse63898": "primary hepatocellular carcinoma"}
DESCRIPTIVE = ("gse31210", "gse41271", "gse63898")
tertile = L.A166.tertile
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


def zscore_scores(X, m, by=None):
    """G2M and purity-proxy scores (mean gene z, z-scored) over the tumours, within `by` groups if given (GSE68465: site)."""
    m = m.copy()
    groups = m.groupby(by).indices.values() if by else [np.arange(len(m))]
    for i in groups:
        for col, genes in (("prolif", L.sh.hallmark_g2m()), ("purity", L.sh.PURITY)):
            s = L.sh._zmean(X.iloc[i], [g for g in genes if g in X.columns])
            m.iloc[i, m.columns.get_loc(col)] = ((s - s.mean()) / s.std()).to_numpy()
    return m


def tumours(key):
    """(meta of the tumours with scores over the tumours, expression, Y of the clock genes, waveform, clock genes)."""
    m, X = L.load(key, OUT)
    if key in TUMOUR:
        keep = (m.tissue.str.strip() == TUMOUR[key]).to_numpy()
        m, X = m[keep], X[keep]
    by = "site" if key == "gse68465" else None
    m = zscore_scores(X, m, by).assign(pos=np.arange(len(m)))
    cg = L.clock_genes(X)
    Y = L.Yof(X, cg)
    if by:  # registered: Y centred per gene within site
        Y = Y - pd.DataFrame(Y).groupby(m[by].to_numpy()).transform("mean").to_numpy()
    return m, X, Y, L.wf_of(cg), cg


class Prolif:
    """H166.3 statistic on purity-proxy tertile strata (>= 30 tumours); null: proliferation permuted within stratum."""
    def __init__(self, purity, prolif):
        pt = tertile(purity)
        self.strata = [np.flatnonzero(pt == k) for k in range(3) if (pt == k).sum() >= MIN_STRATUM]
        self.prolif = np.asarray(prolif, float)

    def T(self, Y, wf, prolif=None):
        p = self.prolif if prolif is None else prolif
        d = []
        for i in self.strata:
            t = tertile(p[i])
            d.append(L.coh(Y, i[t == 2], wf) - L.coh(Y, i[t == 0], wf))
        return float(np.mean(d))

    def null(self, Y, wf, n, rng):
        out = np.empty(n)
        for k in range(n):
            p = self.prolif.copy()
            for i in self.strata:
                p[i] = rng.permutation(p[i])
            out[k] = self.T(Y, wf, p)
        return out


class Groups:
    """coh(group) - coh(rest) within strata (each arm >= 10), mean over strata; null: labels permuted within stratum."""
    def __init__(self, label, strata_codes):
        self.lab = np.asarray(label)
        self.strata = [np.flatnonzero(strata_codes == k) for k in np.unique(strata_codes)]
        self.strata = [i for i in self.strata if min(self.lab[i].sum(), (~self.lab[i]).sum()) >= MIN_ARM]

    def T(self, Y, wf, lab=None):
        a = self.lab if lab is None else lab
        return float(np.mean([L.coh(Y, i[a[i]], wf) - L.coh(Y, i[~a[i]], wf) for i in self.strata])) if self.strata else np.nan

    def null(self, Y, wf, n, rng):
        out = np.empty(n)
        for k in range(n):
            a = self.lab.copy()
            for i in self.strata:
                a[i] = rng.permutation(a[i])
            out[k] = self.T(Y, wf, a)
        return out


def nulls(design, Y, wf, k):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, k]))
    return design.null(Y, wf, args.nperm, rng)


def summary(name, series, design, Y, wf, nl, n, ng, observe=True):
    r = dict(analysis=name, series=series, n=n, n_genes=ng, null_mean=nl.mean(), null_sd=nl.std(ddof=1), MDE_80pct=K * nl.std(ddof=1))
    if observe:
        t = design.T(Y, wf)
        r.update(T=t, T_minus_null_mean=t - nl.mean(), p_one_lower=pl(nl, t))
    return r


def tcga(project, k):
    """The H171 statistic in TCGA (seen data, descriptive): Y = log fractions of the 14 clock genes, A97 scores."""
    sys.path.insert(0, str(HERE.parents[2] / "results/figure_5/scripts/cancer"))
    import _shared as S
    t, _ = S.load_tumours(project)
    Y = np.log((t[S.GENES].to_numpy(float) + 0.5) / t.depth.to_numpy(float)[:, None])
    d = Prolif(t.purity, t.prolif)
    nl = nulls(d, Y, L.WF, k)
    return summary(f"TCGA {project}", project, d, Y, L.WF, nl, len(t), len(S.GENES))


DATA = {s: tumours(s) for s in set(TESTS.values()) | set(DESCRIPTIVE)}
DES = {n: Prolif(DATA[s][0].purity, DATA[s][0].prolif) for n, s in TESTS.items()}
for n, s in TESTS.items():
    say(n, s, "tumours", len(DATA[s][0]), "genes", len(DATA[s][4]), "strata", [len(i) for i in DES[n].strata])

if args.stage == "power":
    res = Parallel(n_jobs=args.jobs)(delayed(nulls)(DES[n], DATA[s][2], DATA[s][3], i + 1) for i, (n, s) in enumerate(TESTS.items()))
    nl = dict(zip(TESTS, res))
    pd.DataFrame(nl).to_csv(OUT / "null_draws.csv", index=False)
    pw = pd.DataFrame([summary(n, s, DES[n], *DATA[s][2:4], nl[n], len(DATA[s][0]), len(DATA[s][4]), observe=False) for n, s in TESTS.items()])
    tc = pd.DataFrame(Parallel(n_jobs=2)(delayed(tcga)(p, 20 + j) for j, p in enumerate(("LUAD", "LIHC"))))
    pw.to_csv(OUT / "power.csv", index=False); tc.to_csv(OUT / "tcga.csv", index=False)
    print(pw.round(4).to_string()); print(tc.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests ----------------------------------------------------------------------------------------------------------------------
nl = pd.read_csv(OUT / "null_draws.csv")
assert len(nl) == args.nperm
names = list(TESTS)
obs = {n: DES[n].T(*DATA[TESTS[n]][2:4]) for n in names}
p = np.array([pl(nl[n].to_numpy(), obs[n]) for n in names])
ph = np.minimum(1, np.maximum.accumulate(np.sort(p) * (len(p) - np.arange(len(p)))))[np.argsort(np.argsort(p))]


def pseudo(j, n):
    s = TESTS[n]; m, X, Y, wf, cg = DATA[s]
    sets = L.pseudo_sets(X, L.A166.pool_genes(X), cg, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 10 + j])))
    Yc = lambda g: L.Yof(X, g) - (pd.DataFrame(L.Yof(X, g)).groupby(m.site.to_numpy()).transform("mean").to_numpy() if s == "gse68465" else 0)
    return sets, np.array([DES[n].T(Yc(g), L.wf_pseudo(wf, g)) for g in sets])


ps = Parallel(n_jobs=args.jobs)(delayed(pseudo)(j, n) for j, n in enumerate(names))
rows = []
for j, n in enumerate(names):
    sets, pT = ps[j]; s = TESTS[n]
    pd.DataFrame(sets).to_csv(OUT / f"pseudo_sets_{s}.csv", index=False)
    T, q05, sd = obs[n], np.quantile(pT, 0.05), nl[n].std(ddof=1)
    spec = ("clock-specific" if T < q05 else "coherence loss not specific to the clock") if ph[j] < 0.05 else "not applicable (test failed)"
    rows.append(dict(hypothesis=n, series=s.upper(), T=T, p_one_lower=p[j], p_holm=ph[j], null_mean=nl[n].mean(), null_sd=sd,
                     MDE_80pct=K * sd, n_tumours=len(DATA[s][0]), n_genes=len(DATA[s][4]), pseudo_percentile=100 * np.mean(pT <= T),
                     pseudo_q05=q05, verdict="pass" if ph[j] < 0.05 else "fail", specificity=spec, n_perm=args.nperm, n_pseudo=args.npseudo))
res = pd.DataFrame(rows); res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True)

# per-stratum coherences (descriptive)
st = []
for n, s in list(TESTS.items()) + [(f"desc {d}", d) for d in DESCRIPTIVE]:
    m, X, Y, wf, cg = DATA[s]
    for k, i in enumerate(Prolif(m.purity, m.prolif).strata):
        t = tertile(m.prolif.to_numpy()[i])
        st.append(dict(analysis=n, series=s, stratum=k, n=len(i), **{f"coh_prolif_{lv}": L.coh(Y, i[t == q], wf) for q, lv in enumerate(("low", "mid", "high"))}))
pd.DataFrame(st).to_csv(OUT / "strata.csv", index=False)

# ---- descriptive --------------------------------------------------------------------------------------------------------------
jobs = []
for d in DESCRIPTIVE:
    m, X, Y, wf, cg = DATA[d]
    jobs.append((f"{d.upper()} proliferation (H171 statistic)", d, Prolif(m.purity, m.prolif), Y, wf, len(m), len(cg)))
for n, s in TESTS.items():
    m, X, Y, wf, cg = DATA[s]
    keep = [g for g in cg if g not in ("PER1", "NR1D1")]
    Yk = Y[:, [cg.index(g) for g in keep]]
    jobs.append((f"{n} without PER1 and NR1D1", s, DES[n], Yk, L.wf_of(keep), len(m), len(keep)))
m, X, Y, wf, cg = DATA["gse72094"]
pt, ft = tertile(m.purity), tertile(m.prolif)
for col, pos in (("kras", "Mut"), ("egfr", "Mut"), ("stk11", "Mut"), ("tp53", "Mut"), ("smoking", "Ever")):
    ok = m[col].isin(["Mut", "WT"] if pos == "Mut" else ["Ever", "Never"]).to_numpy()
    lab = (m[col] == pos).to_numpy()
    for sname, codes in (("purity tertile", pt), ("proliferation tertile", ft)):
        idx = np.flatnonzero(ok)
        jobs.append((f"GSE72094 {col} {pos} vs rest, within {sname}", "gse72094", Groups(lab[idx], codes[idx]), Y[idx], wf, int(ok.sum()), len(cg)))
# LIHC tumour vs non-tumour (all samples of the series; scores not used)
for key, normal, tum in (("gse14520", "Non-Tumor", "Tumor"), ("gse36376", "adjacent non-tumor liver", "liver tumor")):
    mm, XX = L.load(key, OUT); cg = L.clock_genes(XX); YY = L.Yof(XX, cg); mm = mm.assign(pos=np.arange(len(mm)))
    tis = mm.tissue.str.strip()
    if key == "gse14520":
        mm["patient"] = mm.patient.str.strip().str[:-1]  # 'Individual' = patient id + A (tumour) / B (non-tumour)
        a = mm[tis == normal].drop_duplicates("patient", keep=False).set_index("patient")
        b = mm[tis == tum].drop_duplicates("patient", keep=False).set_index("patient")
        both = a.index.intersection(b.index).sort_values()
        jobs.append(("GSE14520 tumour vs paired non-tumour (swap within patient)", key,
                     L.Pairs("pairs", "swap", a.loc[both, "pos"].to_numpy(), b.loc[both, "pos"].to_numpy()), YY, L.wf_of(cg), len(both), len(cg)))
    else:
        sel = tis.isin([normal, tum]).to_numpy()
        jobs.append(("GSE36376 tumour vs adjacent (unpaired)", key, Groups((tis == tum).to_numpy()[sel], np.zeros(sel.sum(), int)),
                     YY[sel], L.wf_of(cg), int(sel.sum()), len(cg)))


def desc(j, name, s, d, Y, wf, n, ng):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 100 + j]))
    nl_ = d.null(Y, wf, args.nperm, rng)
    return summary(name, s, d, Y, wf, nl_, n, ng)


dd = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(desc)(j, *jb) for j, jb in enumerate(jobs)))
dd.to_csv(OUT / "descriptive.csv", index=False); print(dd.round(4).to_string(), flush=True)
say("done tests")
