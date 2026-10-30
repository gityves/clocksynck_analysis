"""A170: is the ER− loss of breast-tumour clock coherence carried by the malignant cells? (registered in
paper/docs/registered_experiments.md, A170). Wu 2021 + Pal 2021 per-patient pseudobulks from prep_a170.py.
  --stage power : permutation null of H170.1, MDE = 2.49 x SD (power.csv); no observed statistic is computed.
  --stage tests : observed T, one-sided p, pseudo-clock companion, descriptive rows (per cohort, whole tumour, Wu CAF and
                  endothelial, sensitivities).
    python run_a170.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 2]"""
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
from clocksync.decompose import logfrac  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=2)
args = ap.parse_args()
OUT = HERE.parent / "output"
K, SEED, MIN_CELLS = norm.ppf(0.95) + norm.ppf(0.80), zlib.crc32(b"A170") % (2 ** 31), 50
EXCLUDE = {"ER_0319", "ER_0114_T3", "TN_0114_T2"}  # ER-low; patient 0114 has one tumour in each arm
METABRIC, SCANB = -0.165, -0.315
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)

units = pd.read_csv(OUT / "units.csv").set_index("unit")
counts = []
for c in ("Wu", "Pal"):
    X = pd.read_parquet(OUT / f"counts_{c}.parquet")
    X.columns = [canonical(g) for g in X.columns]
    counts.append(X.T.groupby(level=0).sum().T)
common = sorted(set(counts[0].columns) & set(counts[1].columns))
C = pd.concat([x[common] for x in counts])
C["depth"] = units.depth.reindex(C.index)  # depth over all genes of the cohort's matrix
CLOCK = L.clock_genes(C)
WF = L.wf_of(CLOCK)


def select(cls, cohorts=("Wu", "Pal"), drop=()):
    u = units[(units.cls == cls) & units.cohort.isin(cohorts) & (units.n_cells >= MIN_CELLS)
              & ~units.patient.isin(EXCLUDE | set(drop))]
    return u


def design(u):
    """(units table, Y = log fractions of the clock genes centred within cohort, ER− index, ER+ index, cohort codes)."""
    y, _, _ = logfrac(C.loc[u.index], CLOCK, u.cohort.to_numpy())
    neg, pos = np.flatnonzero(u.er.to_numpy() == "ER-"), np.flatnonzero(u.er.to_numpy() == "ER+")
    return y, neg, pos, u.cohort.to_numpy()


def T(Y, neg, pos, wf):
    return L.coh(Y, neg, wf) - L.coh(Y, pos, wf)


def null(Y, er_neg, coh_codes, wf, n, k):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, k]))
    lab, out = er_neg.copy(), np.empty(n)
    for i in range(n):
        for c in np.unique(coh_codes):  # ER labels permuted within cohort
            m = coh_codes == c
            lab[m] = rng.permutation(lab[m])
        out[i] = T(Y, np.flatnonzero(lab), np.flatnonzero(~lab), wf)
    return out


def run(name, u, k, observe=True):
    Y, neg, pos, co = design(u)
    is_neg = np.zeros(len(u), bool); is_neg[neg] = True
    nl = null(Y, is_neg, co, WF, args.nperm, k)
    r = dict(analysis=name, n_er_neg=len(neg), n_er_pos=len(pos), n_genes=len(CLOCK), null_mean=nl.mean(),
             null_sd=nl.std(ddof=1), MDE_80pct=K * nl.std(ddof=1), n_perm=args.nperm)
    if observe:
        t = T(Y, neg, pos, WF)
        r.update(T=t, T_minus_null_mean=t - nl.mean(), p_one_lower=pl(nl, t), coh_er_neg=L.coh(Y, neg, WF),
                 coh_er_pos=L.coh(Y, pos, WF))
    return r, nl


prim = select("malignant")
say("H170.1 units:", prim.groupby(["cohort", "er"]).size().to_dict(), "clock genes", len(CLOCK), CLOCK)

if args.stage == "power":
    r, nl = run("H170.1", prim, 1, observe=False)
    pd.DataFrame({"H170.1": nl}).to_csv(OUT / "null_draws.csv", index=False)
    pw = pd.DataFrame([dict(r, metabric_effect=METABRIC, scanb_effect=SCANB)])
    pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- H170.1 and its pseudo-clock companion ------------------------------------------------------------------------------
nl = pd.read_csv(OUT / "null_draws.csv")["H170.1"].to_numpy()
assert len(nl) == args.nperm
Y, neg, pos, co = design(prim)
t = T(Y, neg, pos, WF)
p = pl(nl, t)
frac = C.loc[prim.index, common].div(C.loc[prim.index, "depth"], axis=0)
Xl = np.log(frac + 0.5 / C.loc[prim.index, "depth"].to_numpy()[:, None])  # uncentred log fractions (for matching)
pool = [g for g in L.pool_genes(Xl, L.estrogen_genes()) if (frac[g] >= 1e-5).all()]
sets = L.pseudo_sets(Xl, pool, CLOCK, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 10])))
pd.DataFrame(sets).to_csv(OUT / "pseudo_sets.csv", index=False)


def pseudo_T(g):
    y, _, _ = logfrac(C.loc[prim.index], list(g), co)
    return T(y, neg, pos, L.wf_pseudo(WF, g))


pT = np.array(Parallel(n_jobs=args.jobs)(delayed(pseudo_T)(g) for g in sets))
q05 = np.quantile(pT, 0.05)
spec = ("clock-specific" if t < q05 else "coherence loss not specific to the clock") if p < 0.05 else "not applicable (test failed)"
res = pd.DataFrame([dict(hypothesis="H170.1", T=t, p_one_lower=p, null_mean=nl.mean(), null_sd=nl.std(ddof=1),
                         MDE_80pct=K * nl.std(ddof=1), coh_er_neg=L.coh(Y, neg, WF), coh_er_pos=L.coh(Y, pos, WF),
                         n_er_neg=len(neg), n_er_pos=len(pos), n_genes=len(CLOCK), pool_size=len(pool),
                         pseudo_percentile=100 * np.mean(pT <= t), pseudo_q05=q05, verdict="pass" if p < 0.05 else "fail",
                         specificity=spec, n_perm=args.nperm, n_pseudo=args.npseudo)])
res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).T.to_string(), flush=True)

# ---- descriptive --------------------------------------------------------------------------------------------------------
same = set(prim.patient)
specs = [("malignant, Wu only", select("malignant", ("Wu",))), ("malignant, Pal only", select("malignant", ("Pal",))),
         ("whole tumour, same patients", select("whole").query("patient in @same")),
         ("whole tumour, all patients", select("whole")),
         ("Wu CAF", select("caf", ("Wu",))), ("Wu endothelial", select("endothelial", ("Wu",))),
         ("malignant without CID3963", select("malignant", drop={"CID3963"})),
         ("malignant without Pal BRCA1 TNBC", select("malignant", drop={p for p in units.patient if p.startswith("TN_B1_")}))]
rows = Parallel(n_jobs=args.jobs)(delayed(run)(n, u, 100 + j) for j, (n, u) in enumerate(specs))
desc = pd.DataFrame([r for r, _ in rows])
desc.to_csv(OUT / "descriptive.csv", index=False); print(desc.round(4).to_string(), flush=True)
say("done tests")
