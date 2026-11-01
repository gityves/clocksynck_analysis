"""A178 (C7): robustness of the coherence losses to single genes (V178.1) and to the pseudo-pool size (V178.2).
Definitions are those of A176 (run_a176: Contrast, build_all, registered permutation schemes); pseudo sets by the A168 construction.
    python run_a178.py [--nperm 2000] [--nsets 1000] [--jobs 4]"""
import argparse
import sys
import time
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A176_c3_plane_variance/scripts"))
ap = argparse.ArgumentParser()
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--nsets", type=int, default=1000); ap.add_argument("--jobs", type=int, default=4)
args = ap.parse_args()
import run_a176 as D  # noqa: E402,F401
D.args.jobs = args.jobs
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import beta  # noqa: E402

L, A166, EXPR, say, pl = D.L, D.A166, D.EXPR, D.say, D.pl
OUT = HERE.parent / "output"; OUT.mkdir(exist_ok=True)
SEED = zlib.crc32(b"A178") % (2 ** 31)
LOGO = ["H127.1", "H127.2", "H166.1", "H166.2"]
PSEUDO = ["H171.1", "H173.1", "H127.1", "H127.2", "H166.1", "H171.2", "H173.3"]  # first two = V178.2, rest descriptive
COHORT = {"H127.1": "SCAN-B", "H127.2": "SCAN-B", "H166.1": "METABRIC", "H166.2": "METABRIC", "H171.1": "GSE72094", "H171.2": "GSE68465",
          "H173.1": "GSE62254", "H173.3": "GSE15459"}
t0 = time.time()

R = D.build_all()
mb = R["H166.1"]["m"]  # carries pos
h2 = A166.build(mb)[1]
R["H166.2"] = dict(c=D.lab_contrast("H166.2", h2.G, h2.er, 0, 1), X=R["H166.1"]["X"], m=mb, cg=R["H166.1"]["cg"], sets=None, center=None)
for n, r in R.items():
    r["wf"] = L.wf_of(r["cg"]); r["Y"] = D.Yof(r, r["cg"])
    say(n, "groups", [len(i) for i in r["c"].G], "genes", len(r["cg"]))

# ---- first: reproduce every full-gene statistic against the amendment tables (1e-6) -------------------------------------------------
reg = D.registered_dR()
reg["H166.2"] = pd.read_csv(EXPR("A166_metabric") / "a166_tests.csv").set_index("hypothesis")["T"]["H166.2"]
for n in LOGO + ["H171.1", "H171.2", "H173.1", "H173.3"]:
    v = R[n]["c"].dR_registered(R[n]["Y"], R[n]["wf"])
    assert abs(v - reg[n]) < 1e-6, (n, v, reg[n])
say("full-gene statistics reproduced to 1e-6 for all", len(reg) + 0, "contrasts:", {k: round(v, 6) for k, v in reg.items()})


# ---- V178.1 leave-one-gene-out -----------------------------------------------------------------------------------------------------
def logo_task(n, gi, ss):
    r = R[n]; c = r["c"]; cg = r["cg"]
    keep = [g for g in cg if g != cg[gi]]
    Yk = r["Y"][:, [cg.index(g) for g in keep]]; wk = L.wf_of(keep)
    o = c.stats(Yk, wk)[0]
    rng = np.random.default_rng(ss)
    nl = np.array([c.stats(Yk, wk, c.draw(rng))[0] for _ in range(args.nperm)])
    return dict(contrast=n, dropped_gene=cg[gi], statistic=o, p=pl(nl, o), null_sd=nl.std(ddof=1), n_perm=args.nperm)


tasks = [(n, gi, np.random.SeedSequence([SEED, 1, k, gi])) for k, n in enumerate(LOGO) for gi in range(len(R[n]["cg"]))]
logo = pd.DataFrame(Parallel(n_jobs=args.jobs, verbose=0)(delayed(logo_task)(*t) for t in tasks))
full = {n: R[n]["c"].dR_registered(R[n]["Y"], R[n]["wf"]) for n in LOGO}
logo["full_statistic"] = logo.contrast.map(full)
logo.to_csv(OUT / "logo.csv", index=False)
say("LOGO done", f"{time.time() - t0:.0f}s")

# ---- V178.2 pseudo pool of 1000 new sets (A168 construction, new seed) -----------------------------------------------------------------
sets, ps = {}, {}
for j, cohort in enumerate(dict.fromkeys(COHORT[n] for n in PSEUDO)):
    first = next(n for n in PSEUDO if COHORT[n] == cohort); r = R[first]; X = r["X"]
    sets[cohort] = L.pseudo_sets(X, A166.pool_genes(X), r["cg"], args.nsets, np.random.default_rng(np.random.SeedSequence([SEED, 2, j])))
    pd.DataFrame(sets[cohort]).to_csv(OUT / f"pseudo_sets_{cohort}.csv", index=False)
for n in PSEUDO:
    r = R[n]
    ps[n] = np.array(Parallel(n_jobs=args.jobs)(delayed(lambda s: r["c"].stats(D.Yof(r, s), D.wf_for(s, r["wf"]))[0])(s) for s in sets[COHORT[n]]))
    say(n, "pseudo stats done", f"{time.time() - t0:.0f}s")
rows = []
for n in PSEUDO:
    o = full[n] if n in full else R[n]["c"].dR_registered(R[n]["Y"], R[n]["wf"])
    k = int((ps[n] <= o).sum()); N = len(ps[n])
    lo = beta.ppf(0.025, k, N - k + 1) if k > 0 else 0.0; hi = beta.ppf(0.975, k + 1, N - k) if k < N else 1.0
    rows.append(dict(contrast=n, cohort=COHORT[n], role="V178.2" if n in PSEUDO[:2] else "descriptive", statistic=o, n_sets=N, n_below=k,
                     percentile=100 * k / N, ci95_lo=100 * lo, ci95_hi=100 * hi, pseudo_q05=np.quantile(ps[n], 0.05)))
pc = pd.DataFrame(rows); pc.to_csv(OUT / "pseudo1000.csv", index=False)

# ---- verdict ---------------------------------------------------------------------------------------------------------------------------
v = []
for n in LOGO:
    d = logo[logo.contrast == n]
    d2 = d.assign(shift=d.statistic - d.full_statistic)
    big = d2.loc[d2["shift"].abs().idxmax()]
    neg, sig = int((d.statistic < 0).sum()), int((d.p < 0.05).sum())
    v.append(dict(test="V178.1", contrast=n, n_genes=len(d), n_negative=neg, n_p_lt_005=sig,
                  robust=bool(neg == len(d) and sig >= 12), max_p=d.p.max(), most_influential_gene=big.dropped_gene,
                  statistic_without_it=big.statistic, full_statistic=full[n], note="" if (neg == len(d) and sig >= 12) else
                  "failing genes (p >= 0.05 or statistic >= 0): " + ",".join(d[(d.p >= 0.05) | (d.statistic >= 0)].dropped_gene)))
for _, p in pc[pc.role == "V178.2"].iterrows():
    v.append(dict(test="V178.2", contrast=p.contrast, percentile=p.percentile, ci95_lo=p.ci95_lo, ci95_hi=p.ci95_hi,
                  clock_specific=bool(p.percentile < 5), at_the_line=bool(p.ci95_lo <= 5 <= p.ci95_hi)))
v = pd.DataFrame(v)
v1 = v[v.test == "V178.1"].robust.all(); v2 = v[v.test == "V178.2"].clock_specific.all()
v = pd.concat([v, pd.DataFrame([dict(test="overall", contrast="V178.1", robust=bool(v1)), dict(test="overall", contrast="V178.2", clock_specific=bool(v2)),
                                dict(test="overall", contrast="V178.1 and V178.2", robust=bool(v1 and v2))])], ignore_index=True)
v.to_csv(OUT / "verdict.csv", index=False)
print(logo.groupby("contrast").agg(min_stat=("statistic", "min"), max_stat=("statistic", "max"), max_p=("p", "max"), n_p05=("p", lambda s: (s < .05).sum())).to_string())
print(pc.round(3).to_string()); print(v.to_string())
say("done A178", f"{time.time() - t0:.0f}s")
