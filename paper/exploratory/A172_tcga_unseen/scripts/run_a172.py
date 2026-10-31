"""A172: the proliferation clock-coherence loss in the TCGA projects not among the 11 of A97 (registered in
paper/docs/registered_experiments.md, A172). Data: recount3 (fetch_a172.sh), loaded as `common/datasets/tcga_paired`.
  --stage power : eligibility (metadata only), per-project permutation nulls, MDE (eligibility.csv, power.csv, null_draws.csv).
  --stage tests : per-project T, H172.1 (mean over projects), H172.2 (Holm), pseudo-clock companion, descriptive rows.
    python run_a172.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 4]"""
import argparse
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

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=4)
args = ap.parse_args()
OUT = HERE.parent / "output"
OUT.mkdir(exist_ok=True)
K, SEED, MIN_TUMOURS, MIN_STRATUM = norm.ppf(0.95) + norm.ppf(0.80), zlib.crc32(b"A172") % (2 ** 31), 150, 30
PROJECTS = ["ACC", "CESC", "CHOL", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "LAML", "LGG", "LUSC", "MESO", "PAAD", "PCPG", "READ",
            "SARC", "SKCM", "STAD", "TGCT", "THYM", "UCS", "UVM"]
CARCINOMAS = ["CESC", "ESCA", "HNSC", "LUSC", "PAAD", "STAD"]
PRIMARY = {"LAML": "Primary Blood Derived Cancer - Peripheral Blood"}
B = data._builder("tcga_paired")
tertile = L.A166.tertile
GENES = list(L.GENES)
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


def primaries(p):
    s = B.samples_of(p)
    return s[s.sample_type == PRIMARY.get(p, "Primary Tumor")].sort_values("barcode").drop_duplicates("patient")


def zmean(lc, genes):
    g = [x for x in genes if x in lc.columns]
    z = ((lc[g] - lc[g].mean()) / lc[g].std(ddof=0).replace(0, np.nan)).mean(axis=1)
    return ((z - z.mean()) / z.std()).to_numpy()


def load(p):
    """(counts of the 14 clock genes + depth, log2 CPM of genes with mean CPM >= 1, proliferation, purity)."""
    t = primaries(p)
    ct = B.counts_of(p, t)
    depth = ct.sum(axis=1)
    lc = np.log2(ct.div(depth, axis=0) * 1e6 + 1)
    lc = lc.loc[:, (2 ** lc - 1).mean() >= 1]
    clock = ct[GENES].assign(depth=depth.to_numpy())
    return clock, lc.astype(np.float32), zmean(lc, L.sh.hallmark_g2m()), zmean(lc, L.sh.PURITY)


def Yof(counts, genes):
    return np.log((counts[list(genes)].to_numpy(float) + 0.5) / counts.depth.to_numpy(float)[:, None])


class Prolif:
    """H166.3 / H171 statistic: purity-proxy tertile strata (>= 30 tumours), coh(high proliferation tertile) - coh(low)."""
    def __init__(self, purity, prolif):
        pt = tertile(purity)
        self.strata = [np.flatnonzero(pt == k) for k in range(3) if (pt == k).sum() >= MIN_STRATUM]
        self.prolif = np.asarray(prolif, float)

    def T(self, Y, wf, prolif=None):
        p = self.prolif if prolif is None else prolif
        return float(np.mean([L.coh(Y, i[tertile(p[i]) == 2], wf) - L.coh(Y, i[tertile(p[i]) == 0], wf) for i in self.strata]))

    def null(self, Y, wf, n, rng):
        out = np.empty(n)
        for k in range(n):
            q = self.prolif.copy()
            for i in self.strata:
                q[i] = rng.permutation(q[i])
            out[k] = self.T(Y, wf, q)
        return out


def project_power(j, p, clock, prolif, purity):
    d = Prolif(purity, prolif)
    nl = d.null(Yof(clock, GENES), L.WF, args.nperm, np.random.default_rng(np.random.SeedSequence([SEED, j])))
    say(p, "null done")
    return p, nl, len(clock), [len(i) for i in d.strata]


def project_tests(j, p, clock, lc, prolif, purity):
    d = Prolif(purity, prolif)
    T = d.T(Yof(clock, GENES), L.WF)
    sets = L.pseudo_sets(lc, L.A166.pool_genes(lc), GENES, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 1000 + j])))
    ctp = lambda g: np.log((2 ** lc[list(g)].to_numpy(float) - 1) / 1e6 + 0.5 / clock.depth.to_numpy(float)[:, None])
    pT = np.array([d.T(ctp(g), L.wf_pseudo(L.WF, g)) for g in sets])
    keep = [g for g in GENES if g not in ("PER1", "NR1D1")]
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 2000 + j]))
    Yk = Yof(clock, keep); nk = d.null(Yk, L.wf_of(keep), args.nperm, rng); Tk = d.T(Yk, L.wf_of(keep))
    say(p, "tests done")
    return dict(project=p, T=T, pseudo=pT, T_noPER1NR1D1=Tk, p_noPER1NR1D1=pl(nk, Tk), MDE_noPER1NR1D1=K * nk.std(ddof=1))


if args.stage == "power":
    el = []
    for p in PROJECTS:
        n = len(primaries(p))
        el.append(dict(project=p, n_primary_patients=n, eligible=n >= MIN_TUMOURS))
    el = pd.DataFrame(el); el.to_csv(OUT / "eligibility.csv", index=False); print(el.to_string(), flush=True)
    ok = el.project[el.eligible].tolist()
    res = Parallel(n_jobs=args.jobs)(delayed(project_power)(j, p, c, g, u) for j, p in enumerate(ok)
                                        for c, _, g, u in [load(p)])  # loaded in the parent: the builder module does not pickle
    nl = pd.DataFrame({p: x for p, x, _, _ in res}); nl.to_csv(OUT / "null_draws.csv", index=False)
    rows = [dict(analysis=p, n=n, strata=str(s), null_mean=x.mean(), null_sd=x.std(ddof=1), MDE_80pct=K * x.std(ddof=1)) for p, x, n, s in res]
    for name, cols in (("H172.1 (mean over eligible)", ok), ("carcinomas only", [c for c in CARCINOMAS if c in ok]),
                       ("without LAML", [c for c in ok if c != "LAML"])):
        m = nl[cols].mean(axis=1)
        rows.append(dict(analysis=name, n=len(cols), strata="", null_mean=m.mean(), null_sd=m.std(ddof=1), MDE_80pct=K * m.std(ddof=1)))
    pw = pd.DataFrame(rows); pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests ------------------------------------------------------------------------------------------------------------------------
nl = pd.read_csv(OUT / "null_draws.csv")
assert len(nl) == args.nperm
ok = list(nl.columns)
res = {r["project"]: r for r in Parallel(n_jobs=args.jobs)(delayed(project_tests)(j, p, *load(p)) for j, p in enumerate(ok))}
pd.DataFrame({p: res[p]["pseudo"] for p in ok}).to_csv(OUT / "pseudo_T.csv", index=False)
per = pd.DataFrame([dict(project=p, T=res[p]["T"], null_mean=nl[p].mean(), null_sd=nl[p].std(ddof=1), MDE_80pct=K * nl[p].std(ddof=1),
                         p_one_lower=pl(nl[p].to_numpy(), res[p]["T"]), pseudo_percentile=100 * np.mean(res[p]["pseudo"] <= res[p]["T"]),
                         pseudo_q05=np.quantile(res[p]["pseudo"], 0.05), T_noPER1NR1D1=res[p]["T_noPER1NR1D1"],
                         p_noPER1NR1D1=res[p]["p_noPER1NR1D1"]) for p in ok])
o = np.argsort(per.p_one_lower.to_numpy())
per["p_holm"] = np.minimum(1, np.maximum.accumulate(np.sort(per.p_one_lower) * (len(per) - np.arange(len(per)))))[np.argsort(o)]
per["verdict"] = np.where(per.p_holm < 0.05, "pass", "fail")
per["specificity"] = np.where(per.p_holm < 0.05, np.where(per["T"] < per.pseudo_q05, "clock-specific", "not specific"), "")
per.to_csv(OUT / "per_project.csv", index=False); print(per.round(4).to_string(), flush=True)
rows = []
for name, cols in (("H172.1 (mean over eligible)", ok), ("carcinomas only", [c for c in CARCINOMAS if c in ok]),
                   ("without LAML", [c for c in ok if c != "LAML"])):
    T = float(np.mean([res[c]["T"] for c in cols])); m = nl[cols].mean(axis=1).to_numpy()
    ps = np.mean([res[c]["pseudo"] for c in cols], axis=0)
    rows.append(dict(analysis=name, projects=";".join(cols), T=T, null_mean=m.mean(), null_sd=m.std(ddof=1), MDE_80pct=K * m.std(ddof=1),
                     p_one_lower=pl(m, T), pseudo_percentile=100 * np.mean(ps <= T), pseudo_q05=np.quantile(ps, 0.05)))
pooled = pd.DataFrame(rows); pooled.to_csv(OUT / "pooled.csv", index=False); print(pooled.round(4).to_string(), flush=True)
say("done tests")
