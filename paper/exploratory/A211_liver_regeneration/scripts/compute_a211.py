"""A211: does synchronous physiological proliferation (mouse liver regeneration after partial hepatectomy, GSE180012) swamp the liver clock?

Registered in paper/docs/registered_experiments.md (A211). Data: $SCCC_DATA/gse180012_phx (parse_gse180012.py).
  --stage gate  : manipulation gate (a) G2M PHx > sham over 20-32 h (label swaps within time point), (b) reference_gate on the sham samples
                  -> output/samples.csv, output/gate.csv. No coherence of either arm is computed.
  --stage power : null of H211.1 (arm labels swapped within time point), MDE = 2.49 x null SD -> output/power.csv, output/null_draws.csv.
                  No observed coherence or T is computed or printed.
  --stage tests : observed T of H211.1, one-sided permutation p, pseudo-clock companion (200 sets), descriptive splits 4-16 h / 20-32 h,
                  clock-gene cosinor per arm, day/night series -> output/tests.csv, pseudo_T_draws.csv, cosinor.csv, dn_descriptive.csv.
Seeds: crc32("A211"); stream k of SeedSequence([SEED, k]).
    python -I compute_a211.py --stage gate|power|tests [--nperm 2000] [--npseudo 200]
"""
import argparse
import zlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "A168_endocrine/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a168_lib as L  # noqa: E402  (coh, WF, GENES, hallmark, pseudo_sets, pool_genes, wf_pseudo)
import clocksync as cs  # noqa: E402
from clocksync.paired import reference_gate  # noqa: E402
from common.paths import DATA  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["gate", "power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
args = ap.parse_args()
D = DATA / "gse180012_phx"
OUT = Path(__file__).resolve().parents[1] / "output"; OUT.mkdir(exist_ok=True)
SEED, K, MDE_MAX = zlib.crc32(b"A211") % 2 ** 31, norm.ppf(0.95) + norm.ppf(0.80), 0.3
SURGERY_H = 8.5  # surgery 08:00-09:00 (GEO extract protocol); time_h = hours after surgery
WINDOWS = {"all_4-32h": (4, 32), "injury_4-16h": (4, 16), "replicative_20-32h": (20, 32)}
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)
rng_of = lambda k: np.random.default_rng(np.random.SeedSequence([SEED, k]))

S = pd.read_csv(D / "samples.csv")
S = S[S.series == "main"].set_index("title")
X = pd.read_parquet(D / "expression_main.parquet").loc[S.index]          # 48 samples x genes, log2(RPKM + 1)
S["pos"] = np.arange(len(S))
S.reset_index().to_csv(OUT / "samples.csv", index=False)
phx = (S.arm == "PHx").to_numpy(); tpt = S.time_h.to_numpy()
CG = L.GENES
Y = L.Yof(X, CG)                                                           # ln2 x log2 expression (A168 convention)


def win(w):
    return np.where((tpt >= WINDOWS[w][0]) & (tpt <= WINDOWS[w][1]))[0]


def swap(idx, rng):
    """Arm labels (bool PHx) with the labels permuted within each time point of the samples idx."""
    lab = phx.copy()
    for t in np.unique(tpt[idx]):
        j = idx[tpt[idx] == t]
        lab[j] = rng.permutation(phx[j])
    return lab


def T_coh(Yc, wf, idx, lab):
    return L.coh(Yc, idx[lab[idx]], wf) - L.coh(Yc, idx[~lab[idx]], wf)


def null_T(w, n, rng):
    idx = win(w)
    return np.array([T_coh(Y, L.WF, idx, swap(idx, rng)) for _ in range(n)])


if args.stage == "gate":
    g2m = [g for g in L.hallmark("HALLMARK_G2M_CHECKPOINT") if g in X.columns]
    Z = X[g2m]; prolif = ((Z - Z.mean()) / Z.std(ddof=0).replace(0, np.nan)).mean(axis=1).to_numpy()   # A97 _zmean over the 48 samples
    idx = win("replicative_20-32h")
    Tg = prolif[idx][phx[idx]].mean() - prolif[idx][~phx[idx]].mean()
    rng = rng_of(1)
    nl = []
    for _ in range(args.nperm):
        lab = swap(idx, rng); nl.append(prolif[idx][lab[idx]].mean() - prolif[idx][~lab[idx]].mean())
    pa = cs.stats.perm_p(np.array(nl), Tg, "greater", tie=0.0)
    byt = {f"g2m_diff_{t}h": float(prolif[(tpt == t) & phx].mean() - prolif[(tpt == t) & ~phx].mean()) for t in np.unique(tpt)}
    sham = np.where(~phx)[0]
    g = reference_gate(X.iloc[sham][CG].to_numpy(float), L.WF, n_perm=200, rng=211)
    gate = pd.DataFrame([
        dict(gate="a_proliferation", n_genes=len(g2m), stat=Tg, p=pa, n_swaps=args.nperm, passes=bool(pa < 0.05), **byt),
        dict(gate="b_sham_clock_visible", n=len(sham), **{k: v for k, v in g.items() if np.ndim(v) == 0}, passes=bool(g["p"] < 0.05))])
    gate["both_pass"] = bool(gate.passes.all())
    gate.to_csv(OUT / "gate.csv", index=False)
    print(gate.round(4).T.to_string())

elif args.stage == "power":
    rows, draws = [], {}
    for k, w in enumerate(WINDOWS):
        say(w); nl = null_T(w, args.nperm, rng_of(100 + k)); draws[w] = nl
        mde = K * nl.std(ddof=1)
        rows.append(dict(window=w, n_per_arm=len(win(w)) // 2, null_mean=nl.mean(), null_sd=nl.std(ddof=1), mde=mde,
                         registered=(w == "all_4-32h"), estimable=bool(mde <= MDE_MAX), n_swaps=args.nperm))
    P = pd.DataFrame(rows); P.to_csv(OUT / "power.csv", index=False); pd.DataFrame(draws).to_csv(OUT / "null_draws.csv", index=False)
    print(P.round(4).to_string())

else:  # tests
    assert (OUT / "power.csv").exists() and (OUT / "gate.csv").exists(), "run --stage gate and power first"
    pw, gt = pd.read_csv(OUT / "power.csv").set_index("window"), pd.read_csv(OUT / "gate.csv")
    nd = pd.read_csv(OUT / "null_draws.csv")
    pool = L.pool_genes(X)
    sets = L.pseudo_sets(X, pool, CG, args.npseudo, rng_of(10))
    pd.DataFrame(sets).to_csv(OUT / "pseudo_sets.csv", index=False)
    rows, pT = [], {}
    for w in WINDOWS:
        idx = win(w)
        T = T_coh(Y, L.WF, idx, phx)
        cP, cS = L.coh(Y, idx[phx[idx]], L.WF), L.coh(Y, idx[~phx[idx]], L.WF)
        p = cs.stats.perm_p(nd[w].to_numpy(), T, "less", tie=0.0)
        pT[w] = np.array([T_coh(L.Yof(X, g), L.wf_pseudo(L.WF, g), idx, phx) for g in sets])
        q05 = np.quantile(pT[w], 0.05)
        rows.append(dict(window=w, registered=(w == "all_4-32h"), n_per_arm=len(idx) // 2, coh_PHx=cP, coh_sham=cS, T=T, p_less=p,
                         pseudo_percentile=100 * np.mean(pT[w] <= T), pseudo_q05=q05, pseudo_median=np.median(pT[w]), T_below_pseudo_q05=bool(T < q05),
                         mde=pw.loc[w, "mde"], estimable=bool(pw.loc[w, "estimable"]), gates_pass=bool(gt.both_pass.iloc[0]), n_swaps=len(nd), n_pseudo=len(sets)))
    R = pd.DataFrame(rows)
    r = R.iloc[0]
    R["verdict"] = ("uninformative (gate failed)" if not r.gates_pass else "not estimable (MDE > 0.3)" if not r.estimable else
                    ("pass, clock-specific" if r.T_below_pseudo_q05 else "pass, not clock-specific") if r.p_less < 0.05 else "fail (powered)")
    R.to_csv(OUT / "tests.csv", index=False); pd.DataFrame(pT).to_csv(OUT / "pseudo_T_draws.csv", index=False)
    print(R.round(4).T.to_string())
    # clock-gene cosinor per arm over the 4-32 h series (24 h period; acrophase as clock hour assuming surgery at 08:30)
    cos = []
    for arm, m in (("sham", ~phx), ("PHx", phx)):
        th = tpt[m] * 2 * np.pi / 24
        A = np.column_stack([np.ones(m.sum()), np.cos(th), np.sin(th)])
        B, *_ = np.linalg.lstsq(A, X[CG].to_numpy(float)[m], rcond=None)
        res = X[CG].to_numpy(float)[m] - A @ B
        for j, gname in enumerate(CG):
            ph = np.arctan2(B[2, j], B[1, j]) * 24 / (2 * np.pi) % 24
            cos.append(dict(arm=arm, gene=gname, mean=B[0, j], amp_log2=float(np.hypot(B[1, j], B[2, j])), acrophase_h_after_surgery=float(ph),
                            acrophase_clock_h=float((ph + SURGERY_H) % 24), r2=float(1 - res[:, j].var() / X[CG].to_numpy(float)[m][:, j].var())))
    pd.DataFrame(cos).to_csv(OUT / "cosinor.csv", index=False)
    # day / night surgery series (file 2, log2(TPM + 1) of the clock genes): group means and PHx - sham
    dn = pd.read_csv(D / "expression_dn_clock.csv", index_col="sample")
    grp = pd.Series(dn.index.str.replace(r"-\d+$", "", regex=True), index=dn.index)
    mean = dn.groupby(grp).mean().T
    for t in ("4", "24"):
        for dd in ("D", "N"):
            mean[f"{dd}{t}_PHx_minus_sham"] = mean[f"{dd}{t}P"] - mean[f"{dd}{t}S"]
    mean.to_csv(OUT / "dn_descriptive.csv")
    say("done")
