"""A315 (registered): within-tissue clock-plane amplitude of brain-dead organ donors vs ventilated post-mortem donors.

A109 `Tissues(200, seed 0)` (12 GTEx body tissues), groups as A132 (Hardy 0, COHORT organ donor vs post-mortem). Per
tissue and group the debiased clock-plane variance W = mean |z - mean z|^2 - mean Poisson variance; L_t = log(W_OD /
W_PM); Q = median over tissues. Matched (200) and responsive (200, score stage) sets give Q_k. Aggregates only.

    python score_a315.py --stage power   # matched-set Q spread and the clock's relabelling null (no clock contrast)
    python score_a315.py --stage score   # H315.1 and the descriptive readouts
Outputs (`../output/`): power.csv | tissues.csv, sets.csv, tests.csv, descriptive.csv.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A109_coupling_nulls" / "scripts"))
sys.path.insert(0, str(EXP / "A132_brain_death" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

import _helpers as H  # noqa: E402
import compute_brain_death as B  # noqa: E402
import compute_h109_responsive as R  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
MIN_N, N_BOOT, N_PERM = 20, 1000, 1000
TARGET = np.log(0.7)


def W(T, kind, t, donors):
    """Debiased clock-plane variance of tissue t over `donors` (repeats allowed), for every set of `kind` (K,)."""
    r = T.R[t][donors]
    r = r[r >= 0]
    if len(r) < MIN_N:
        return None
    z, v = T.Z[kind][t][:, r], T.V[kind][t][:, r]
    return np.mean(np.abs(z - z.mean(1, keepdims=True)) ** 2, 1) - v.mean(1)


def L_table(T, kind, a, b):
    """log(W_a / W_b) per tissue (rows) and set (columns); NaN where a group is too small or W <= 0."""
    out = {}
    for t in H.BODY:
        wa, wb = W(T, kind, t, a), W(T, kind, t, b)
        if wa is None or wb is None:
            continue
        with np.errstate(invalid="ignore", divide="ignore"):
            out[t] = np.where((wa > 0) & (wb > 0), np.log(wa / wb), np.nan)
    return pd.DataFrame(out).T


def Q(Lt):
    return np.nanmedian(Lt.to_numpy(), axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    if a.stage == "score":
        R.add_sets(T, 200, 0)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    h0 = T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    print(f"states built ({time.time() - t0:.0f}s): organ donors {len(od)}, post-mortem {len(pm)}", flush=True)
    rng = np.random.default_rng(315)
    Qm = Q(L_table(T, "pseudo", od, pm))
    if a.stage == "power":
        perm = []
        for _ in range(N_PERM):
            p = rng.permutation(h0)
            perm.append(Q(L_table(T, "clock", p[:len(od)], p[len(od):]))[0])
        sd = max(float(np.nanstd(Qm, ddof=1)), float(np.nanstd(perm, ddof=1)))
        res = dict(n_od=len(od), n_pm=len(pm), sd_matched_Q=float(np.nanstd(Qm, ddof=1)),
                   sd_relabel_Q=float(np.nanstd(perm, ddof=1)), sd_used=sd, target=TARGET,
                   power=float(norm.cdf(-TARGET / sd - 1.645)),
                   tissues=int(L_table(T, "pseudo", od, pm).notna().any(axis=1).sum()))
        res["gate"] = bool(res["power"] >= 0.8)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).to_string(), flush=True)
        return
    Lc = L_table(T, "clock", od, pm)
    q = float(Q(Lc)[0])
    Qr = Q(L_table(T, "resp", od, pm))
    boot, tb = [], []
    for _ in range(N_BOOT):
        bo, bp = rng.choice(od, len(od)), rng.choice(pm, len(pm))
        Lb = L_table(T, "clock", bo, bp)
        boot.append(Q(Lb)[0])
        tb.append(Lb[0].reindex(Lc.index))
    boot = np.array(boot)
    tb = pd.concat(tb, axis=1)
    hi = float(np.nanquantile(boot, 0.95))
    q05 = float(np.nanquantile(Qm, 0.05))
    tests = [dict(test="H315.1", Q=q, ratio=float(np.exp(q)), boot_lo90=float(np.nanquantile(boot, 0.05)), boot_hi90=hi,
                  matched_median=float(np.nanmedian(Qm)), matched_q05=q05, percentile_vs_matched=float(np.nanmean(Qm < q)),
                  Q_star=float(q - np.nanmedian(Qm)), ratio_star=float(np.exp(q - np.nanmedian(Qm))),
                  responsive_median=float(np.nanmedian(Qr)), percentile_vs_responsive=float(np.nanmean(Qr < q)),
                  tissues=int(Lc[0].notna().sum()), passes=bool(q < 0 and q < q05 and hi < 0))]
    tab = pd.DataFrame(dict(L=Lc[0], lo90=tb.quantile(0.05, axis=1), hi90=tb.quantile(0.95, axis=1),
                            matched_median=L_table(T, "pseudo", od, pm).median(axis=1)))
    tab["n_od"] = [int((T.R[t][od] >= 0).sum()) for t in tab.index]
    tab["n_pm"] = [int((T.R[t][pm] >= 0).sum()) for t in tab.index]
    tab.to_csv(OUT / "tissues.csv")
    pd.DataFrame(dict(matched_Q=Qm)).to_csv(OUT / "sets.csv", index=False)
    pd.DataFrame(tests).to_csv(OUT / "tests.csv", index=False)
    # descriptive: fast vs post-mortem; ventilation-duration tertiles within organ donors (A132 cuts 61 / 101 h)
    fast = T.group_donors("fast")
    dur = pd.to_numeric(ph.DTHVNTD, errors="coerce").to_numpy()
    top, bot = od[dur[od] > 101], od[dur[od] <= 61]
    desc = []
    for name, x, y in (("fast_vs_postmortem", fast, pm), ("od_top_vs_bottom_tertile", top, bot)):
        q_c, q_m = float(Q(L_table(T, "clock", x, y))[0]), Q(L_table(T, "pseudo", x, y))
        desc.append(dict(contrast=name, n_a=len(x), n_b=len(y), Q=q_c, matched_median=float(np.nanmedian(q_m)),
                         Q_star=float(q_c - np.nanmedian(q_m)), percentile_vs_matched=float(np.nanmean(q_m < q_c))))
    pd.DataFrame(desc).to_csv(OUT / "descriptive.csv", index=False)
    print(tab.round(3).to_string())
    print(pd.DataFrame(tests).round(4).T.to_string())
    print(pd.DataFrame(desc).round(3).to_string(index=False))
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
