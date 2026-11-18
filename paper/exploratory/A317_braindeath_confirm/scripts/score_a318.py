"""A318 (registered): is the A317 brain-death arrangement change explained by ischaemic time or RNA quality?

H318.1: per tissue, Hardy 0 clock states (and every matched / responsive set) residualised on standardised ischaemic
time and RIN (complex least squares, pooled groups); D_rho as A317 H317.1 on the adjusted states. H318.2 (descriptive):
within organ donors, top vs bottom tertile of donor-level median ischaemic time. Aggregates only.

    python score_a318.py --stage power | score
Outputs (`../output/`): a318_power.csv | a318_tests.csv, a318_sets.csv.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a317 as A  # noqa: E402
from score_a317 import B, H, RS, gtex  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"


def adjust(T, h0, attrs):
    """Residualise Hardy 0 rows of every kind on (ischaemic, RIN); returns R with rows lacking attributes removed."""
    R = {}
    for t in H.BODY:
        r = T.R[t].copy()
        a = attrs[t]
        isch, rin = a.ischemic_time.to_numpy(float), a.rin.to_numpy(float)
        rows = r[h0][r[h0] >= 0]
        ok = rows[np.isfinite(isch[rows]) & np.isfinite(rin[rows])]
        bad = np.setdiff1d(rows, ok)
        drop = np.isin(r, bad) & (r >= 0)
        r[drop] = -1
        R[t] = r
        X = np.column_stack([np.ones(len(ok)), (isch[ok] - isch[ok].mean()) / isch[ok].std(), (rin[ok] - rin[ok].mean()) / rin[ok].std()])
        P = X @ np.linalg.pinv(X)
        for kind in ("clock", "pseudo", "resp"):
            Z = T.Z[kind][t]
            sub = Z[:, ok]
            T.Z[kind][t][:, ok] = sub - (P @ sub.T).T + sub.mean(1, keepdims=True)
    return R


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    RS.add_sets(T, 200, 0)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    attrs = A.sample_attrs(T)
    # aggregates of the confounds per group (medians over samples)
    agg = {}
    for g, d in (("organ_donor", od), ("postmortem", pm)):
        v = {c: np.concatenate([attrs[t][c].to_numpy(float)[T.R[t][d][T.R[t][d] >= 0]] for t in H.BODY]) for c in ("ischemic_time", "rin")}
        agg.update({f"{g}_ischemic_median": float(np.nanmedian(v["ischemic_time"])), f"{g}_rin_median": float(np.nanmedian(v["rin"]))})
    h2 = None
    if a.stage == "score":
        # H318.2 inputs (unadjusted states)
        don_isch = pd.Series({d: np.nanmedian([attrs[t].ischemic_time.to_numpy(float)[T.R[t][d]] for t in H.BODY if T.R[t][d] >= 0])
                              for d in od})
        lo_c, hi_c = don_isch.quantile([1 / 3, 2 / 3])
        top, bot = np.array(don_isch[don_isch > hi_c].index), np.array(don_isch[don_isch <= lo_c].index)
        Zc0, Mc0 = A.raw(T, "clock", 0, T.R)
        m2 = B.pair_mask(Mc0, fast, top, bot)
        rng = np.random.default_rng(318)
        h2 = A.rho(Zc0, Mc0, fast, top, m2) - A.rho(Zc0, Mc0, fast, bot, m2)
        n2 = []
        for _ in range(A.N_PERM):
            p = rng.permutation(np.concatenate([top, bot]))
            n2.append(A.rho(Zc0, Mc0, fast, p[:len(top)], m2) - A.rho(Zc0, Mc0, fast, p[len(top):], m2))
    # H318.1: adjusted states
    R = adjust(T, h0, attrs)
    od_ok = od[np.array([any(R[t][d] >= 0 for t in H.BODY) for d in od])]
    pm_ok = pm[np.array([any(R[t][d] >= 0 for t in H.BODY) for d in pm])]
    rng = np.random.default_rng(A.SEED)
    draws = [rng.choice(od_ok, len(pm_ok), replace=False) for _ in range(A.N_DRAWS)]
    Zc, Mc = A.raw(T, "clock", 0, R)
    mask = B.pair_mask(Mc, fast, pm_ok, draws[0]) & B.pair_mask(Mc, od_ok)
    h0_ok = np.concatenate([od_ok, pm_ok])
    rn = []
    for _ in range(A.N_PERM):
        p = rng.permutation(h0_ok)
        rn.append(A.rho(Zc, Mc, fast, p[:len(pm_ok)], mask) - A.rho(Zc, Mc, fast, p[len(pm_ok):2 * len(pm_ok)], mask))
    rn = np.array(rn)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "power":
        res = dict(pairs=int(mask.sum()), n_od=len(od_ok), n_pm=len(pm_ok), sd_D_rho=float(np.nanstd(rn, ddof=1)),
                   mde_D_rho=float(2.49 * np.nanstd(rn, ddof=1)), **agg)
        pd.Series(res).to_csv(OUT / "a318_power.csv", header=["value"])
        print(pd.Series(res).to_string(), flush=True)
        return
    D = A.d_rho(Zc, Mc, fast, pm_ok, draws, mask)
    sets = []
    for kind, lab in (("pseudo", "matched"), ("resp", "responsive")):
        for k in range(200):
            Z, M = A.raw(T, kind, k, R)
            sets.append(dict(family=lab, set=k, D_rho=A.d_rho(Z, M, fast, pm_ok, draws, mask)))
    sets = pd.DataFrame(sets)
    sets.to_csv(OUT / "a318_sets.csv", index=False)
    q = sets.groupby("family").D_rho.quantile(0.95)
    p_rho = float((1 + np.sum(rn >= D)) / (1 + len(rn)))
    ex = A.excess_C(T, pm_ok, draws, R)
    tests = [dict(test="H318.1", pairs=int(mask.sum()), D_rho=D, p_labels=p_rho, matched_q95=float(q["matched"]),
                  responsive_q95=float(q["responsive"]),
                  pct_matched=float(np.mean(sets[sets.family == "matched"].D_rho < D)),
                  pct_responsive=float(np.mean(sets[sets.family == "responsive"].D_rho < D)),
                  passes=bool(D > 0 and p_rho < 0.05 and D > q["matched"] and D > q["responsive"]),
                  D_C_adjusted=ex["pm"]["excess"] - ex["od"]["excess"]),
             dict(test="H318.2", n_top=len(top), n_bottom=len(bot), pairs=int(m2.sum()), rho_top_minus_bottom=float(h2),
                  p_two_sided=float((1 + np.sum(np.abs(n2) >= abs(h2))) / (1 + len(n2)))),
             dict(test="aggregates", **agg)]
    pd.DataFrame(tests).to_csv(OUT / "a318_tests.csv", index=False)
    print(pd.DataFrame(tests).round(4).T.to_string())
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
