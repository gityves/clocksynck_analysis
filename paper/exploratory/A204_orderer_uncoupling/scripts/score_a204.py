"""A204 scoring: simulator gate, H204.1-H204.4 (Holm over four), descriptives. Reads <out>/runs/ from compute_a204.py.

    $(conda info --base)/envs/sccc/bin/python paper/exploratory/A204_orderer_uncoupling/scripts/score_a204.py

Outputs (<out>/): gate.csv, tests.csv, cofe_per_rep.csv, cofe_cells.csv, auc.csv, tt_per_rep.csv, tt_cells.csv.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import compute_a204 as A  # noqa: E402

RAD = 2 * np.pi / 24


def resultant(h) -> float:
    return float(np.abs(np.mean(np.exp(1j * np.asarray(h, float) * RAD))))


def rayleigh_p(h) -> float:
    """Rayleigh test of uniformity (Zar's approximation)."""
    n = len(h)
    Rn = n * resultant(h)
    return float(min(1.0, np.exp(np.sqrt(1 + 4 * n + 4 * (n * n - Rn * Rn)) - (1 + 2 * n))))


def largest_gap(h) -> float:
    s = np.sort(np.asarray(h, float) % 24)
    return float(np.max(np.diff(np.r_[s, s[0] + 24])))


def circ_corr(a_h, b_h) -> float:
    """Fisher-Lee circular correlation."""
    a, b = np.asarray(a_h) * RAD, np.asarray(b_h) * RAD
    a0, b0 = stats.circmean(a), stats.circmean(b)
    sa, sb = np.sin(a - a0), np.sin(b - b0)
    return float(np.sum(sa * sb) / np.sqrt(np.sum(sa ** 2) * np.sum(sb ** 2)))


def circ_lin(x, th_h) -> float:
    """Circular-linear correlation (Mardia) of a linear x with a circular time (h)."""
    th = np.asarray(th_h) * RAD
    rc, rs = np.corrcoef(x, np.cos(th))[0, 1], np.corrcoef(x, np.sin(th))[0, 1]
    rcs = np.corrcoef(np.cos(th), np.sin(th))[0, 1]
    return float(np.sqrt((rc ** 2 + rs ** 2 - 2 * rc * rs * rcs) / (1 - rcs ** 2)))


def binom_ge(k, n, p=0.5) -> float:
    return float(stats.binom.sf(k - 1, n, p))


def auc(pos, neg) -> float:
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    return float(np.mean((pos[:, None] > neg[None, :]) + 0.5 * (pos[:, None] == neg[None, :])))


def auc_ci(pos, neg, B=2000, seed=204):
    r = np.random.default_rng(seed)
    v = [auc(r.choice(pos, len(pos)), r.choice(neg, len(neg))) for _ in range(B)]
    return np.quantile(v, [0.025, 0.975])


def cofe_rows(runs: Path, reps: int) -> pd.DataFrame:
    rows = []
    for w in A.COFE_WINDOWS:
        for tr in A.COFE_TRUTHS:
            for rep in range(reps):
                r = pd.read_pickle(runs / f"cofe_{w}_{tr}_{rep:03d}.pkl")
                m = r["matched"]
                diff = (r["phase_T"].loc[m].to_numpy() - r["phase_N"].loc[m].to_numpy()) % 24
                rows.append(dict(window=w, truth=tr, rep=rep, n_matched=len(m),
                                 R_true=resultant(r["t"]), gap_true_h=largest_gap(r["t"]),
                                 R_labels=resultant(r["phase_T"]), gap_labels_h=largest_gap(r["phase_T"]),
                                 abs_cc_labels_true=abs(circ_corr(r["phase_T"].loc[r["t"].index], r["t"])),
                                 R_diff=resultant(diff), rayleigh_p_diff=rayleigh_p(diff),
                                 interior_min=r["interior_min"], cs_corr=r["cs_corr"], cs_offset_h=r["cs_offset_h"],
                                 cs_p=r["cs_p"]))
    return pd.DataFrame(rows)


def tt_rows(runs: Path, reps: int) -> pd.DataFrame:
    rows = []
    for w in A.WINDOWS:
        for c in A.TT_CS:
            for rep in range(reps):
                r = pd.read_pickle(runs / f"tt_{w}_{c:.2f}_{rep:03d}.pkl")
                e = dict(window=w, c_T=c, rep=rep, failed=r["pred"] is None)
                if r["pred"] is not None:
                    t, p = r["t"].to_numpy(), r["pred"].time_pred.to_numpy(float)
                    ok = np.isfinite(p)
                    lr = stats.linregress(t[ok], p[ok])
                    err = (p[ok] - (t[ok] + A.OFFSET_T_H) + 12) % 24 - 12
                    e.update(n=int(ok.sum()), slope=lr.slope, slope_p=lr.pvalue, r=lr.rvalue, circ_lin_r=circ_lin(t[ok], p[ok]),
                             median_abs_err_h=float(np.median(np.abs(err))),
                             share_7_20=float(np.mean((p[ok] % 24 >= 7) & (p[ok] % 24 < 20))),
                             median_theta=float(np.nanmedian(r["pred"].theta)))
                rows.append(e)
    return pd.DataFrame(rows)


def main(out: Path = A.OUT, reps: int = 30):
    runs = out / "runs"
    cf, tt = cofe_rows(runs, reps), tt_rows(runs, reps)
    cf.to_csv(out / "cofe_per_rep.csv", index=False)
    tt.to_csv(out / "tt_per_rep.csv", index=False)

    g = cf.groupby(["window", "truth"])
    cells = g.agg(R_true=("R_true", "median"), gap_true_h=("gap_true_h", "median"), R_labels=("R_labels", "median"),
                  gap_labels_h=("gap_labels_h", "median"), n_gap_lt4=("gap_labels_h", lambda x: int((x < 4).sum())),
                  abs_cc_labels_true=("abs_cc_labels_true", "median"), R_diff=("R_diff", "median"),
                  n_uniform_diff=("rayleigh_p_diff", lambda x: int((x > 0.05).sum())),
                  interior_min_share=("interior_min", "mean"), cs_corr=("cs_corr", "median"),
                  cs_offset_h=("cs_offset_h", lambda x: A.HOURS * float(np.angle(np.mean(np.exp(1j * np.asarray(x) / A.HOURS))))),
                  cs_detect=("cs_p", lambda x: int((x < 0.05).sum())), n_matched=("n_matched", "median"), reps=("rep", "size"))
    cells.reset_index().to_csv(out / "cofe_cells.csv", index=False)

    tc = tt[~tt.failed].groupby(["window", "c_T"]).agg(
        power=("slope_p", lambda x: float((x < 0.05).mean())), n_detect=("slope_p", lambda x: int((x < 0.05).sum())),
        median_slope=("slope", "median"), median_r=("r", "median"), median_circ_lin_r=("circ_lin_r", "median"),
        median_abs_err_h=("median_abs_err_h", "median"), share_7_20=("share_7_20", "median"),
        median_theta=("median_theta", "median"), reps=("rep", "size"))
    tc["failed"] = tt.groupby(["window", "c_T"]).failed.sum()
    tc.reset_index().to_csv(out / "tt_cells.csv", index=False)

    # simulator gate (arm F)
    c_cc = float(cf.query("window == 'F' and truth == 'c0.9'").abs_cc_labels_true.median())
    t_err = float(tt.query("window == 'F' and c_T == 0.9").median_abs_err_h.median())
    gate = pd.DataFrame([dict(check="cofe_F_c0.9_median_abs_circ_corr", value=c_cc, threshold=">= 0.7", passes=c_cc >= 0.7),
                         dict(check="timeteller_F_c0.9_median_abs_err_h", value=t_err, threshold="<= 3", passes=t_err <= 3)])
    gate.to_csv(out / "gate.csv", index=False)
    print(gate.to_string(), flush=True)

    # hypotheses
    W9, WU = cf.query("window == 'W' and truth == 'c0.9'"), cf.query("window == 'W' and truth == 'uncoupled'")
    n = len(W9)
    k1 = int((W9.gap_labels_h < 4).sum())
    k2 = int((W9.rayleigh_p_diff > 0.05).sum())
    k3, f3 = int((W9.cs_p < 0.05).sum()), int((WU.cs_p < 0.05).sum())
    T9 = tt.query("window == 'W' and c_T == 0.9 and not failed")
    k4 = int((T9.slope_p < 0.05).sum())
    tests = pd.DataFrame([
        dict(h="H204.1", readout="W c0.9: COFE largest label gap < 4 h", k=k1, n=n, rule=">= 20", p=binom_ge(k1, n), passes_rule=k1 >= 20),
        dict(h="H204.2", readout="W c0.9: Rayleigh p > 0.05 on COFE tumour-normal differences", k=k2, n=n, rule=">= 20",
             p=binom_ge(k2, n), passes_rule=k2 >= 20),
        dict(h="H204.3", readout=f"W c0.9: clocksync pairing p < 0.05 (companion W uncoupled: {f3} of {len(WU)})", k=k3, n=n,
             rule=">= 24 and companion <= 3", p=binom_ge(k3, n), passes_rule=(k3 >= 24) and (f3 <= 3)),
        dict(h="H204.4", readout="W c0.9 n 108: TimeTeller slope p < 0.05", k=k4, n=len(T9), rule=">= 24",
             p=binom_ge(k4, len(T9)), passes_rule=k4 >= 24),
    ])
    m = len(tests)
    order = np.argsort(tests.p.to_numpy())
    adj, run = np.empty(m), 0.0
    for i, j in enumerate(order):
        run = max(run, min(1.0, (m - i) * tests.p.iloc[j]))
        adj[j] = run
    tests["p_holm"] = adj
    tests["verdict"] = np.where(tests.passes_rule & (tests.p_holm < 0.05), "PASS", "FAIL")
    if not gate.passes.all():
        tests["verdict"] = "not scored (simulator gate failed)"
    tests.to_csv(out / "tests.csv", index=False)
    print(tests.to_string(), flush=True)

    # discrimination on the same replicates (W, c0.9 vs uncoupled; also c0.5 vs uncoupled)
    au = []
    for w in A.COFE_WINDOWS:
        for pos in ("c0.9", "c0.5"):
            P, N = cf.query("window == @w and truth == @pos"), cf.query("window == @w and truth == 'uncoupled'")
            for meth, col in (("cofe_R_diff", "R_diff"), ("clocksync_corr", "cs_corr")):
                lo, hi = auc_ci(P[col].to_numpy(), N[col].to_numpy())
                au.append(dict(window=w, positive=pos, readout=meth, auc=auc(P[col], N[col]), lo95=lo, hi95=hi))
    pd.DataFrame(au).to_csv(out / "auc.csv", index=False)
    print(pd.DataFrame(au).round(3).to_string(), flush=True)
    print(cells.round(3).to_string(), flush=True)
    print(tc.round(3).to_string(), flush=True)
    mde = tc.reset_index().query("power >= 0.8").groupby("window").c_T.min()
    print("TimeTeller MDE (smallest c_T with power >= 0.8):", mde.to_dict(), flush=True)


if __name__ == "__main__":
    main()
