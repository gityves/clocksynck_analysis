"""A173: proliferation clock-coherence loss in gastric cancer (H173.1-3) and adenocarcinoma vs squamous in lung (H173.4)
(registered in paper/docs/registered_experiments.md, A173). Data from prep_a173.py.
  --stage power : permutation nulls, MDE = 2.49 x SD (power.csv, null_draws.csv); no observed statistic is computed.
  --stage tests : observed statistics, one-sided p, Holm (H173.1-3), pseudo-clock companion; descriptive rows.
    python run_a173.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 4]"""
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
K, SEED, MIN_STRATUM = norm.ppf(0.95) + norm.ppf(0.80), zlib.crc32(b"A173") % (2 ** 31), 30
GASTRIC = {"H173.1": "gse62254", "H173.2": "gse84437", "H173.3": "gse15459"}
LUNG = {"gse37745": ("adeno", "squamous"), "gse41271": ("Adenocarcinoma", "Squamous"), "gse30219": ("ADC", "SQC")}
tertile = L.A166.tertile
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


def tertile_effect(Y, wf, strata, p):
    """Mean over strata of coh(high proliferation tertile) - coh(low) (H166.3 / H171 statistic)."""
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

    def stats(self, Y, wf, p=None):
        return {"T": tertile_effect(Y, wf, self.strata, self.prolif if p is None else p)}

    def null(self, Y, wf, n, rng):
        return pd.DataFrame([self.stats(Y, wf, permute_within(self.prolif, self.strata, rng)) for _ in range(n)])


class Histology:
    """H173.4 per series: T_h on purity-proxy halves within histology h (median split, each half >= 30); D = T_ADC - T_SQC.
    Null: proliferation permuted within histology x half (one permutation gives D, T_ADC and T_SQC)."""
    def __init__(self, m, adc, sqc):
        h = m.histology.str.strip().to_numpy()
        pur = m.purity.to_numpy(float)
        self.strata = {}
        for lab, name in ((adc, "ADC"), (sqc, "SQC")):
            i = np.flatnonzero(h == lab)
            lo = pur[i] <= np.median(pur[i])
            self.strata[name] = [j for j in (i[lo], i[~lo]) if len(j) >= MIN_STRATUM]
        self.prolif = m.prolif.to_numpy(float)

    def stats(self, Y, wf, p=None):
        p = self.prolif if p is None else p
        a, s = (tertile_effect(Y, wf, self.strata[k], p) for k in ("ADC", "SQC"))
        return {"D": a - s, "T_ADC": a, "T_SQC": s}

    def null(self, Y, wf, n, rng):
        allst = self.strata["ADC"] + self.strata["SQC"]
        return pd.DataFrame([self.stats(Y, wf, permute_within(self.prolif, allst, rng)) for _ in range(n)])


def series(key):
    m, X = L.load(key, OUT)
    cg = L.clock_genes(X)
    d = Histology(m, *LUNG[key]) if key in LUNG else Prolif(m)
    return dict(m=m, X=X, cg=cg, Y=L.Yof(X, cg), wf=L.wf_of(cg), d=d)


KEYS = list(GASTRIC.values()) + list(LUNG) + ["gse26253"]
DATA = {k: series(k) for k in KEYS}
for k, s in DATA.items():
    st = s["d"].strata
    say(k, "tumours", len(s["m"]), "genes", len(s["cg"]),
        "strata", {h: [len(i) for i in v] for h, v in st.items()} if isinstance(st, dict) else [len(i) for i in st])


def null_of(j, key, genes=None):
    s = DATA[key]
    g = s["cg"] if genes is None else genes
    Y = s["Y"] if genes is None else L.Yof(s["X"], g)
    return s["d"].null(Y, L.wf_of(g), args.nperm, np.random.default_rng(np.random.SeedSequence([SEED, j])))


MAIN = list(GASTRIC.values()) + list(LUNG)
if args.stage == "power":
    res = dict(zip(MAIN, Parallel(n_jobs=args.jobs)(delayed(null_of)(j + 1, k) for j, k in enumerate(MAIN))))
    nl = pd.concat({k: v for k, v in res.items()}, axis=1)
    nl.columns = [f"{a}|{b}" for a, b in nl.columns]
    nl.to_csv(OUT / "null_draws.csv", index=False)
    rows = [dict(analysis=n, series=k, n=len(DATA[k]["m"]), null_sd=nl[f"{k}|T"].std(ddof=1), MDE_80pct=K * nl[f"{k}|T"].std(ddof=1))
            for n, k in GASTRIC.items()]
    Dm = nl[[f"{k}|D" for k in LUNG]].mean(axis=1)
    rows.append(dict(analysis="H173.4 mean D", series=";".join(LUNG), n=sum(len(DATA[k]["m"]) for k in LUNG),
                     null_sd=Dm.std(ddof=1), MDE_80pct=K * Dm.std(ddof=1)))
    for k in LUNG:
        for c in ("D", "T_ADC", "T_SQC"):
            rows.append(dict(analysis=f"{c} (descriptive)", series=k, n=len(DATA[k]["m"]), null_sd=nl[f"{k}|{c}"].std(ddof=1),
                             MDE_80pct=K * nl[f"{k}|{c}"].std(ddof=1)))
    pw = pd.DataFrame(rows); pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests --------------------------------------------------------------------------------------------------------------------------
nl = pd.read_csv(OUT / "null_draws.csv")
assert len(nl) == args.nperm
obs = {k: DATA[k]["d"].stats(DATA[k]["Y"], DATA[k]["wf"]) for k in MAIN}


def pseudo(j, key):
    s = DATA[key]
    sets = L.pseudo_sets(s["X"], L.A166.pool_genes(s["X"]), s["cg"], args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 10 + j])))
    return sets, pd.DataFrame([s["d"].stats(L.Yof(s["X"], g), L.wf_pseudo(s["wf"], g)) for g in sets])


ps = dict(zip(MAIN, Parallel(n_jobs=args.jobs)(delayed(pseudo)(j, k) for j, k in enumerate(MAIN))))
for k in MAIN:
    pd.DataFrame(ps[k][0]).to_csv(OUT / f"pseudo_sets_{k}.csv", index=False)
    ps[k][1].to_csv(OUT / f"pseudo_stats_{k}.csv", index=False)

p = np.array([pl(nl[f"{k}|T"].to_numpy(), obs[k]["T"]) for k in GASTRIC.values()])
ph = np.minimum(1, np.maximum.accumulate(np.sort(p) * (len(p) - np.arange(len(p)))))[np.argsort(np.argsort(p))]
rows = []
for j, (n, k) in enumerate(GASTRIC.items()):
    T, pT, x = obs[k]["T"], ps[k][1]["T"].to_numpy(), nl[f"{k}|T"]
    rows.append(dict(hypothesis=n, series=k.upper(), stat=T, p_one_lower=p[j], p_holm=ph[j], null_sd=x.std(ddof=1), MDE_80pct=K * x.std(ddof=1),
                     n_tumours=len(DATA[k]["m"]), n_genes=len(DATA[k]["cg"]), pseudo_percentile=100 * np.mean(pT <= T),
                     pseudo_q05=np.quantile(pT, 0.05), verdict="pass" if ph[j] < 0.05 else "fail"))
D = float(np.mean([obs[k]["D"] for k in LUNG])); Dn = nl[[f"{k}|D" for k in LUNG]].mean(axis=1).to_numpy()
Dp = np.mean([ps[k][1]["D"].to_numpy() for k in LUNG], axis=0)
pD = pl(Dn, D)
rows.append(dict(hypothesis="H173.4", series=";".join(k.upper() for k in LUNG), stat=D, p_one_lower=pD, p_holm=pD, null_sd=Dn.std(ddof=1),
                 MDE_80pct=K * Dn.std(ddof=1), n_tumours=sum(len(DATA[k]["m"]) for k in LUNG), n_genes=14,
                 pseudo_percentile=100 * np.mean(Dp <= D), pseudo_q05=np.quantile(Dp, 0.05), verdict="pass" if pD < 0.05 else "fail"))
res = pd.DataFrame(rows)
res["specificity"] = np.where(res.verdict == "pass", np.where(res.stat < res.pseudo_q05, "clock-specific", "not specific"), "")
res["n_perm"], res["n_pseudo"] = args.nperm, args.npseudo
res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True)

# ---- descriptive --------------------------------------------------------------------------------------------------------------------
desc = []
for k in LUNG:
    for c in ("D", "T_ADC", "T_SQC"):
        x = nl[f"{k}|{c}"].to_numpy()
        desc.append(dict(analysis=f"{k.upper()} {c}", stat=obs[k][c], p_one_lower=pl(x, obs[k][c]), null_sd=x.std(ddof=1), MDE_80pct=K * x.std(ddof=1),
                         pseudo_percentile=100 * np.mean(ps[k][1][c].to_numpy() <= obs[k][c])))
for c in ("T_ADC", "T_SQC"):
    x = nl[[f"{k}|{c}" for k in LUNG]].mean(axis=1).to_numpy(); o = float(np.mean([obs[k][c] for k in LUNG]))
    desc.append(dict(analysis=f"lung mean {c}", stat=o, p_one_lower=pl(x, o), null_sd=x.std(ddof=1), MDE_80pct=K * x.std(ddof=1),
                     pseudo_percentile=100 * np.mean(np.mean([ps[k][1][c].to_numpy() for k in LUNG], axis=0) <= o)))
x = nl[[f"{k}|T" for k in GASTRIC.values()]].mean(axis=1).to_numpy(); o = float(np.mean([obs[k]["T"] for k in GASTRIC.values()]))
desc.append(dict(analysis="gastric mean of H173.1-3", stat=o, p_one_lower=pl(x, o), null_sd=x.std(ddof=1), MDE_80pct=K * x.std(ddof=1),
                 pseudo_percentile=100 * np.mean(np.mean([ps[k][1]["T"].to_numpy() for k in GASTRIC.values()], axis=0) <= o)))
extra = [("GSE26253 (FFPE DASL) T", "gse26253", None, "T")]
extra += [(f"{k.upper()} without PER1 and NR1D1 {c}", k, [g for g in DATA[k]["cg"] if g not in ("PER1", "NR1D1")], c)
          for k, c in [(k, "T") for k in GASTRIC.values()] + [(k, "D") for k in LUNG]]
xn = Parallel(n_jobs=args.jobs)(delayed(null_of)(100 + j, k, g) for j, (_, k, g, _) in enumerate(extra))
for (name, k, g, c), x in zip(extra, xn):
    s = DATA[k]; gg = s["cg"] if g is None else g
    o = s["d"].stats(L.Yof(s["X"], gg), L.wf_of(gg))[c]
    desc.append(dict(analysis=name, stat=o, p_one_lower=pl(x[c].to_numpy(), o), null_sd=x[c].std(ddof=1), MDE_80pct=K * x[c].std(ddof=1)))
dd = pd.DataFrame(desc); dd.to_csv(OUT / "descriptive.csv", index=False); print(dd.round(4).to_string(), flush=True)
say("done tests")
