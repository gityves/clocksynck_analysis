"""A145 scoring: H145.1-3 (Holm), D145.1 (all readouts), D145.2 (A129 known answer), failures, from compute_a145.py outputs.

Rules (registration + the choices in ../README.md):
- coupling error per replicate = mean over the 5 targets of |c_hat - c_true|; clocksync's estimates are used whether or not the
  gate passes (forced), the gate is reported (H145.3 counts a claim only if it passes).
- best two-step method of an arm = lowest mean error (H145.1, H145.2) or lowest false-coupling rate (H145.3) over the
  replicates in which its orderer fit succeeded (the timed oracle is a reference, not a pipeline).
- paired comparison over the replicates where that method succeeded; H145.1/2 one-sided exact sign test (ties dropped),
  H145.3 one-sided exact McNemar test over targets x replicates (discordant pairs). Sensitivity: orderer failure counted as a
  clocksync win (column p_failure_as_win). Holm over the three.
- offsets and amplitudes are scored on targets with true c > 0 (the phase of an uncoupled compartment is undefined).
    python score_a145.py [--out DIR]
Outputs: h145.csv, d145_1.csv, d145_2.csv, failures.csv.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import lib as A  # noqa: E402

TWO_STEP = list(A.METHODS)
ARMS_ORDER = list(A.ARMS)


def circ_corr(a_h, b_h) -> float:
    """Jammalamadaka-SenGupta circular correlation of two angle series (hours)."""
    a, b = np.asarray(a_h, float) / A.HOURS, np.asarray(b_h, float) / A.HOURS
    sa = np.sin(a - np.angle(np.exp(1j * a).sum()))
    sb = np.sin(b - np.angle(np.exp(1j * b).sum()))
    den = np.sqrt((sa ** 2).sum() * (sb ** 2).sum())
    return float((sa * sb).sum() / den) if den > 0 else np.nan


def wrap(h):
    return (np.asarray(h, float) + 12.0) % 24.0 - 12.0


def sign_p(wins: int, n: int) -> float:
    """One-sided exact binomial P(X >= wins), X ~ Bin(n, 1/2); 1 if n = 0."""
    return float(stats.binom.sf(wins - 1, n, 0.5)) if n > 0 else 1.0


def holm(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p)
    adj, run = np.empty(len(p)), 0.0
    for j, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - j) * p[i]))
        adj[i] = run
    return adj


def readouts(d: pd.DataFrame) -> dict:
    """Readouts of one (arm, method) over its non-failed rows."""
    ok = d[~d.failed.astype(bool)]
    n_rep, n_ok = d.rep.nunique(), ok.rep.nunique()
    r = dict(n_reps=n_rep, n_failed=n_rep - n_ok)
    if ok.empty:
        return r
    r["spearman_c"] = stats.spearmanr(ok.c_hat, ok.c_true).statistic
    r["mae_c"] = (ok.c_hat - ok.c_true).abs().mean()
    pos = ok[ok.c_true > 0]
    r["circ_corr_offset"] = circ_corr(pos.offset_hat_h, pos.offset_true_h)
    r["median_abs_offset_err_h"] = float(np.median(np.abs(wrap(pos.offset_hat_h - pos.offset_true_h))))
    r["corr_log_amp"] = float(np.corrcoef(pos.log_amp_hat, np.log(pos.amp_true))[0, 1])
    r["mae_log2_amp"] = float((np.abs(pos.log_amp_hat - np.log(pos.amp_true)) / np.log(2)).mean())
    r["false_coupling_share"] = float((ok.c_hat > 0.2).mean())
    if ok.method.iloc[0] == "clocksync":
        g = ok.groupby("rep").gate_pass.first().astype(bool)
        r["gate_pass_rate"] = float(g.mean())
        r["false_claim_share"] = float(((ok.c_hat > 0.2) & ok.gate_pass.astype(bool)).mean())
        gp = ok[ok.gate_pass.astype(bool)]
        r["mae_c_gate_passing"] = float((gp.c_hat - gp.c_true).abs().mean()) if len(gp) else np.nan
    return r


def per_rep_mae(d: pd.DataFrame) -> pd.Series:
    return (d.c_hat - d.c_true).abs().groupby(d.rep).mean()


def compare_mae(L: pd.DataFrame, arm: str, hid: str) -> dict:
    a = L[L.arm == arm]
    reps = sorted(a.rep.unique())
    cs_mae = per_rep_mae(a[a.method == "clocksync"]).reindex(reps)
    best, best_mae, tab = None, np.inf, {}
    for m in TWO_STEP:
        d = a[(a.method == m) & ~a.failed.astype(bool)]
        tab[m] = d.c_hat.sub(d.c_true).abs().mean() if len(d) else np.nan
        if len(d) and tab[m] < best_mae:
            best, best_mae = m, tab[m]
    d = a[a.method == best]
    ok_reps = [r for r in reps if not d[d.rep == r].failed.astype(bool).any()]
    ts = per_rep_mae(d[d.rep.isin(ok_reps)]).reindex(ok_reps)
    diff = cs_mae.reindex(ok_reps) - ts
    wins, losses = int((diff < 0).sum()), int((diff > 0).sum())
    n_fail = len(reps) - len(ok_reps)
    return dict(hypothesis=hid, arm=arm, best_two_step=best, n_reps=len(reps), n_paired=len(ok_reps), n_failed_best=n_fail,
                mae_clocksync=float(cs_mae.reindex(ok_reps).mean()), mae_best_two_step=float(ts.mean()),
                mae_all_two_step=";".join(f"{m}={tab[m]:.3f}" for m in TWO_STEP),
                clocksync_wins=wins, two_step_wins=losses, ties=len(ok_reps) - wins - losses,
                p_one_sided=sign_p(wins, wins + losses), p_failure_as_win=sign_p(wins + n_fail, wins + losses + n_fail))


def compare_false(L: pd.DataFrame, arm: str, hid: str) -> dict:
    a = L[L.arm == arm]
    cs_ = a[a.method == "clocksync"].set_index(["rep", "target"])
    cs_claim = ((cs_.c_hat > 0.2) & cs_.gate_pass.astype(bool))
    best, best_rate, rates = None, np.inf, {}
    for m in TWO_STEP:
        d = a[(a.method == m) & ~a.failed.astype(bool)]
        rates[m] = float((d.c_hat > 0.2).mean()) if len(d) else np.nan
        if len(d) and rates[m] < best_rate:
            best, best_rate = m, rates[m]
    d = a[(a.method == best) & ~a.failed.astype(bool)].set_index(["rep", "target"])
    ts_false = d.c_hat > 0.2
    cc = cs_claim.reindex(ts_false.index)
    b = int((ts_false & ~cc).sum())  # two-step false, clocksync not: clocksync better
    c = int((~ts_false & cc).sum())
    n_fail = a[a.method == best].rep.nunique() - d.index.get_level_values(0).nunique()
    return dict(hypothesis=hid, arm=arm, best_two_step=best, n_reps=a.rep.nunique(), n_paired=d.index.get_level_values(0).nunique(),
                n_failed_best=n_fail, rate_clocksync_claim=float(cc.mean()), rate_best_two_step=float(ts_false.mean()),
                rate_all_two_step=";".join(f"{m}={rates[m]:.3f}" for m in TWO_STEP),
                discordant_clocksync_better=b, discordant_two_step_better=c, p_one_sided=sign_p(b, b + c),
                p_failure_as_win=np.nan)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(A.OUT))
    args = ap.parse_args(argv)
    out = Path(args.out)
    L = pd.read_csv(out / "results_long.csv")
    L["failed"] = L.failed.astype(bool)
    # D145.1
    rows = []
    for arm in ARMS_ORDER:
        for m in ["clocksync"] + TWO_STEP + ["oracle"]:
            d = L[(L.arm == arm) & (L.method == m)]
            if len(d):
                rows.append(dict(arm=arm, method=m, **readouts(d)))
    d1 = pd.DataFrame(rows)
    d1.to_csv(out / "d145_1.csv", index=False)
    # H145
    h = pd.DataFrame([compare_mae(L, "W", "H145.1"), compare_mae(L, "S_large", "H145.2"), compare_false(L, "N", "H145.3")])
    h["p_holm"] = holm(h.p_one_sided.to_numpy())
    h["verdict"] = np.where(h.p_holm < 0.05, "pass", "fail")
    h.to_csv(out / "h145.csv", index=False)
    # failures
    fr = []
    rt = pd.read_csv(out / "runtime.csv")
    for (arm, m), g in rt[rt.method.isin(TWO_STEP)].groupby(["arm", "method"]):
        errs = g[g.failed.astype(bool)].error.fillna("").astype(str).str.slice(0, 80)
        fr.append(dict(arm=arm, method=m, n_reps=len(g), n_failed=int(g.failed.astype(bool).sum()),
                       errors="; ".join(sorted(set(errs))), median_runtime_s=float(g.runtime_s.median())))
    fl = pd.DataFrame(fr)
    fl["arm"] = pd.Categorical(fl.arm, ARMS_ORDER)
    fl.sort_values(["arm", "method"]).to_csv(out / "failures.csv", index=False)
    # D145.2
    if (out / "real_long.csv").exists():
        R = pd.read_csv(out / "real_long.csv")
        C = pd.read_csv(out / "real_clocksync.csv").rename(columns={"coupling": "c_hat", "offset_h": "offset_hat_h"}).assign(method="clocksync")
        R = pd.concat([C[["genotype", "method", "c_hat", "offset_hat_h", "passes"]],
                       R.rename(columns={"offset_h": "offset_hat_h"})[["genotype", "method", "failed", "c_hat", "offset_hat_h", "log_amp_rel"]]],
                      ignore_index=True)
        rows = []
        for m, g in R.groupby("method", sort=False):
            g = g.set_index("genotype")
            if "AlbCre" in g.index and "BLKO" in g.index and not g["failed"].fillna(False).astype(bool).any():
                c0, c1 = g.c_hat["AlbCre"], g.c_hat["BLKO"]
                rows.append(dict(method=m, c_albcre=c0, c_blko=c1, drop=c1 - c0, drop_relative=(c1 - c0) / c0 if c0 > 0 else np.nan,
                                 sees_drop=bool(c1 < c0), offset_albcre_h=g.offset_hat_h["AlbCre"], offset_blko_h=g.offset_hat_h["BLKO"],
                                 failed=False))
            else:
                rows.append(dict(method=m, failed=True))
        rows.append(dict(method="clocksync_A111_contrast_reference", c_albcre=0.967103, c_blko=0.676373, drop=-0.290730,
                         drop_relative=-0.290730 / 0.967103, sees_drop=True, failed=False,
                         note="a111_h111_3.csv (key Liver, AL arm): pooled-host contrast, d* -0.287 [-0.411, -0.118]"))
        pd.DataFrame(rows).to_csv(out / "d145_2.csv", index=False)
    pd.set_option("display.width", 250)
    print(d1.round(3).to_string()); print(h.round(4).T.to_string()); print(fl.to_string())
    if (out / "d145_2.csv").exists():
        print(pd.read_csv(out / "d145_2.csv").round(3).to_string())


if __name__ == "__main__":
    main()
