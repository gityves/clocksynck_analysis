"""A168: does aromatase-inhibitor ER blockade reduce breast-tumour clock coherence? POETIC difference-in-differences (H168.1), POETIC treated
change (H168.2), GSE20181 replication (H168.3).
  --stage power : permutation-null SDs and MDE = 2.49 x SD (power.csv, null_draws.csv). No observed statistic is computed or printed.
  --stage tests : observed T, one-sided permutation p, Holm over the three, pseudo-clock companion (200 variance-matched sets per series),
                  registered sensitivities and descriptive tables (needs only the cached tables of prep_a168.py).
Seeds: crc32("A168"); null of test k: SeedSequence([SEED, k, chunk]) (chunks of 250 permutations: independent of --jobs); pseudo sets
SeedSequence([SEED, 10 + series]); sensitivity j: SeedSequence([SEED, 100 + j]).
    python -I run_a168.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 2] [--quick] [--out DIR]"""
import argparse
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
sys.path.insert(0, str(Path(__file__).parent))
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
ap.add_argument("--out", default=None)
args = ap.parse_args()
if args.quick:
    args.nperm, args.npseudo = 50, 10
OUT = Path(args.out or (L.OUT / "quick" if args.quick else L.OUT)); OUT.mkdir(parents=True, exist_ok=True)
K, CHUNK, SEED = norm.ppf(0.95) + norm.ppf(0.80), 250, L.SEED  # K = 2.49
METABRIC, SCANB = -0.165, -0.315  # ER effect beyond the null mean (registration)
SERIES_OF = {"H168.1": "poetic", "H168.2": "poetic", "H168.3": "letrozole"}
NAMES = {"H168.1": "POETIC difference-in-differences: [coh(AI,S)-coh(AI,B)] - [coh(untreated,S)-coh(untreated,D)]",
         "H168.2": "POETIC treated change: coh(AI,S)-coh(AI,B)", "H168.3": "GSE20181 (letrozole): coh(day 10-14)-coh(pre)"}
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
holm = lambda p: np.minimum(1, np.maximum.accumulate(np.sort(p) * (len(p) - np.arange(len(p)))))[np.argsort(np.argsort(p))]
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)

# ---- data and registered designs ---------------------------------------------------------------------------------------------------
M, XX, YY, WFS, DES = {}, {}, {}, {}, {}
for s in ("poetic", "letrozole"):
    M[s], XX[s] = L.load(s); cg = L.clock_genes(XX[s])
    YY[s], WFS[s] = L.Yof(XX[s], cg), L.wf_of(cg)
    DES.update(L.build(s, M[s]))
for k, d in DES.items():
    say(k, "series", SERIES_OF[k], "genes", YY[SERIES_OF[k]].shape[1], d.counts())


def chunked_null(k, d, Y, wf, n, cohf=L.coh, stream=None):
    """Permutation null of design d with n draws in chunks of CHUNK (each with its own SeedSequence child; results do not depend on jobs)."""
    sizes = [CHUNK] * (n // CHUNK) + ([n % CHUNK] if n % CHUNK else [])
    rngs = [np.random.default_rng(np.random.SeedSequence([SEED, k if stream is None else stream, c])) for c in range(len(sizes))]
    return np.concatenate([d.null(Y, wf, m, r, cohf) for m, r in zip(sizes, rngs)])


def job_null(k, name):
    s = SERIES_OF[name]
    return name, chunked_null(k, DES[name], YY[s], WFS[s], args.nperm)

if args.stage == "power":
    # ---- (1) power: null draws only; no observed statistic is computed or printed ------------------------------------------------------
    res = Parallel(n_jobs=args.jobs)(delayed(job_null)(i + 1, n) for i, n in enumerate(DES))
    nulls = dict(res)
    pd.DataFrame(nulls).to_csv(OUT / "null_draws.csv", index=False)
    pw = pd.DataFrame([dict(hypothesis=n, description=NAMES[n], series=L.SERIES[SERIES_OF[n]], n_patients=DES[n].n, null_mean=nulls[n].mean(),
                            null_sd=nulls[n].std(ddof=1), MDE_80pct=K * nulls[n].std(ddof=1), metabric_effect=METABRIC, scanb_effect=SCANB,
                            n_perm=args.nperm, **({"n_ai": int((DES[n].arm == 1).sum()), "n_untreated": int((DES[n].arm == 0).sum())} if DES[n].kind == "did" else {}))
                        for n in DES])
    pw["mde_le_metabric"] = pw.MDE_80pct <= abs(METABRIC); pw["mde_le_scanb"] = pw.MDE_80pct <= abs(SCANB)
    pw.to_csv(OUT / "power.csv", index=False); print(pw.drop(columns="description").round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ======================================= stage tests ==============================================================================
# ---- (2) observed T, p, Holm -----------------------------------------------------------------------------------------------------
nd = OUT / "null_draws.csv"
nulls = pd.read_csv(nd) if nd.exists() else None
if nulls is None or len(nulls) != args.nperm:
    say("recomputing nulls (same seeds)"); nulls = pd.DataFrame(dict(Parallel(n_jobs=args.jobs)(delayed(job_null)(i + 1, n) for i, n in enumerate(DES))))
names = list(DES)
obs = {n: DES[n].T(YY[SERIES_OF[n]], WFS[SERIES_OF[n]]) for n in names}
p = np.array([pl(nulls[n].to_numpy(), obs[n]) for n in names]); ph = holm(p)

# ---- (3) specificity: pseudo-clock sets ------------------------------------------------------------------------------------------------
est = L.estrogen_genes(); pseudo, pT = {}, {}
for j, s in enumerate(("poetic", "letrozole")):
    X, cg = XX[s], L.clock_genes(XX[s])
    used = X.loc[M[s].tp.isin({"poetic": ["B", "S", "D"], "letrozole": ["A", "B"]}[s]).to_numpy()]
    sets = L.pseudo_sets(used, L.pool_genes(used, est), cg, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED, 10 + j])))
    pseudo[s] = sets; pd.DataFrame(sets).to_csv(OUT / f"pseudo_sets_{s}.csv", index=False)
    for n in (x for x in names if SERIES_OF[x] == s):
        pT[n] = np.array([DES[n].T(L.Yof(X, g), L.wf_pseudo(WFS[s], g)) for g in sets])
pd.DataFrame(pT).to_csv(OUT / "pseudo_T_draws.csv", index=False)
rows = []
for n, pp, hh in zip(names, p, ph):
    T, q05 = obs[n], np.quantile(pT[n], 0.05); sd = nulls[n].std(ddof=1)
    spec = (("clock-specific" if T < q05 else "coherence loss not specific to the clock") if hh < 0.05 else
            ("not applicable (test failed); T below pseudo q05" if T < q05 else "not applicable (test failed)"))
    rows.append(dict(hypothesis=n, description=NAMES[n], series=L.SERIES[SERIES_OF[n]], T=T, p_one_lower=pp, p_holm=hh, null_mean=nulls[n].mean(), null_sd=sd,
                     MDE_80pct=K * sd, metabric_effect=METABRIC, scanb_effect=SCANB, n_patients=DES[n].n,
                     pseudo_percentile=100 * np.mean(pT[n] <= T), pseudo_q05=q05, pseudo_median=np.median(pT[n]), T_below_pseudo_q05=bool(T < q05),
                     verdict="pass" if hh < 0.05 else "fail", specificity=spec, n_perm=args.nperm, n_pseudo=args.npseudo))
res = pd.DataFrame(rows); res.to_csv(OUT / "tests.csv", index=False); print(res.drop(columns="description").round(4).to_string(), flush=True)

# ---- descriptive tables: coherence per arm x time point; treatment-effect scores ----------------------------------------------------------
Mz, Xz = L.load("z1031b"); Yz, wz = L.Yof(Xz, L.clock_genes(Xz)), L.wf_of(L.clock_genes(Xz))
ALL = {"poetic": (M["poetic"], YY["poetic"], WFS["poetic"]), "letrozole": (M["letrozole"], YY["letrozole"], WFS["letrozole"]), "z1031b": (Mz, Yz, wz)}
cg_rows = []
for s, (m, Y, wf) in ALL.items():
    for (a, t), g in m.groupby(["arm", "tp"]):
        cg_rows.append(dict(series=L.SERIES[s], arm=a, tp=t, n_samples=len(g), coherence=L.coh(Y, g.pos.to_numpy(), wf)))
pd.DataFrame(cg_rows).to_csv(OUT / "coherence_groups.csv", index=False)
pd.DataFrame(cg_rows).pipe(lambda d: print(d.round(4).to_string()))


def change_row(s, label, m, pre, post, arms=None):
    a, b, _ = L.pairs_from_meta(m, pre, post, arms)
    r = dict(series=L.SERIES[s], contrast=label, n_pairs=len(a))
    for c in ("prolif", "purity", "estro"):
        d = m[c].to_numpy()[b] - m[c].to_numpy()[a]
        r.update({f"d_{c}": d.mean(), f"sd_{c}": d.std(ddof=1), f"p_{c}": stats.ttest_1samp(d, 0).pvalue})
    return r


mz2 = Mz.assign(tp=Mz.tp.replace({"reanalyzed": "base"}))
ch = [change_row("poetic", "AI B->S", M["poetic"], "B", "S", ["AI"]), change_row("poetic", "untreated D->S", M["poetic"], "D", "S", ["UNT"]),
      change_row("letrozole", "A->B (pre -> day 10-14)", M["letrozole"], "A", "B"), change_row("letrozole", "A->C (pre -> day 90)", M["letrozole"], "A", "C"),
      change_row("z1031b", "BL -> 2-4 wk", Mz, "base", "w2"), change_row("z1031b", "BL or re-analyzed -> 2-4 wk", mz2, "base", "w2")]
ch = pd.DataFrame(ch); ch.to_csv(OUT / "treatment_effect_scores.csv", index=False); print(ch.round(4).to_string(), flush=True)

# ---- (4) sensitivities -------------------------------------------------------------------------------------------------------------------
def spec_list():
    """(label, series, design, clock genes dropped, residualise) of every sensitivity; Y, waveform and cohf are built in sens_job."""
    S = []
    for lab, drop in (("without PER1", ["PER1"]), ("without PER1 and NR1D1", ["PER1", "NR1D1"])):
        S += [(lab, SERIES_OF[n], DES[n], drop, False) for n in names]
    S += [("Y residualised on proliferation and purity proxy within arm x time point", SERIES_OF[n], DES[n], [], True) for n in names]
    pm = M["poetic"]; ai = DES["H168.2"]; dpro = pm.prolif.to_numpy()[ai.post] - pm.prolif.to_numpy()[ai.pre]
    third = L.A166.tertile(dpro)
    for k, nm in enumerate(("lowest", "middle", "highest")):
        S.append((f"POETIC treated, {nm} third of change in proliferation score", "poetic", L.Pairs("H168.2", "swap", ai.pre[third == k], ai.post[third == k]), [], False))
    ml = M["letrozole"]
    for resp in ("responder", "nonresponder"):
        x, y, _ = L.pairs_from_meta(ml[ml.response == resp], "A", "B"); S.append((f"GSE20181 {resp}s", "letrozole", L.Pairs("H168.3", "swap", x, y), [], False))
    x, y, _ = L.pairs_from_meta(ml, "A", "C"); S.append(("GSE20181 pre vs day 90", "letrozole", L.Pairs("H168.3", "swap", x, y), [], False))
    return S


def zspecs():
    out = []
    for lab, m in (("GSE87411 baseline (BL) vs 2-4 wk, without re-analyzed baselines", Mz), ("GSE87411 baseline incl. re-analyzed vs 2-4 wk", mz2)):
        x, y, _ = L.pairs_from_meta(m, "base", "w2"); out.append((lab, "z1031b", L.Pairs("H168.3", "swap", x, y), [], False))
    return out


def sens_job(j, lab, s, d, drop, resid):
    m, Y0, wf0 = ALL[s]; X = XX[s] if s != "z1031b" else Xz
    cg = [g for g in L.clock_genes(X) if g not in drop]
    Y, wf = L.Yof(X, cg), L.wf_of(cg)
    cohf = L.residual_coh(m[["prolif", "purity"]].to_numpy()) if resid else L.coh
    T = d.T(Y, wf, cohf); nl = chunked_null(0, d, Y, wf, args.nperm, cohf, stream=100 + j)
    return dict(analysis=lab, statistic=d.name, series=L.SERIES[s], T=T, p_one_lower=pl(nl, T), null_sd=nl.std(ddof=1), MDE_80pct=K * nl.std(ddof=1),
                n_genes=len(cg), n_perm=args.nperm, **d.counts())


specs = spec_list() + zspecs()
sens = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(sens_job)(j, *sp) for j, sp in enumerate(specs)))
sens.to_csv(OUT / "sensitivity.csv", index=False); print(sens.round(4).to_string(), flush=True)
say("done tests")
