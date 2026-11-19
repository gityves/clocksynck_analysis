"""A327 (registered, check B3): does the brain-death arrangement change depend on whether a tissue can be reached by
the procurement flush? A317 data, groups, draws and estimator; pair classes P (both possibly reached), N (both not
reached), X (mixed); Delta = D_rho(N) - D_rho(P u X) with a donor-label relabelling null. Aggregates only.

    python score_a327.py --stage power | score
Outputs (`../output/`): a327_power.csv | a327_tests.csv
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
from score_a317 import B, H, gtex  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

REACHED = ["Heart - Atrial Appendage", "Lung", "Adipose - Visceral (Omentum)", "Esophagus - Mucosa"]
NOT_REACHED = ["Artery - Tibial", "Nerve - Tibial", "Skin - Sun Exposed (Lower leg)", "Muscle - Skeletal",
               "Adipose - Subcutaneous"]
N_PERM, SEED = 2000, 3270


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    rng0 = np.random.default_rng(A.SEED)
    od_draws = [rng0.choice(od, len(pm), replace=False) for _ in range(A.N_DRAWS)]
    Z, M = A.raw(T, "clock", 0, T.R)
    mask = B.pair_mask(M, fast, pm, *od_draws[:1]) & B.pair_mask(M, od)
    cls = np.array(["" for _ in H.PAIRS], dtype=object)
    for i, (x, y) in enumerate(H.PAIRS):
        tx, ty = H.BODY[x], H.BODY[y]
        if tx in REACHED and ty in REACHED:
            cls[i] = "P"
        elif tx in NOT_REACHED and ty in NOT_REACHED:
            cls[i] = "N"
        elif {tx, ty} <= set(REACHED + NOT_REACHED):
            cls[i] = "X"
    mN, mPX = mask & (cls == "N"), mask & np.isin(cls, ["P", "X"])
    masks = {c: mask & (cls == c) for c in ("P", "N", "X")} | {"PX": mPX}
    rng = np.random.default_rng(SEED)
    null_delta, null_N = [], []
    for _ in range(N_PERM):
        p = rng.permutation(h0)
        gp, go = p[:len(pm)], p[len(pm):2 * len(pm)]
        dN = A.rho(Z, M, fast, gp, mN) - A.rho(Z, M, fast, go, mN)
        dPX = A.rho(Z, M, fast, gp, mPX) - A.rho(Z, M, fast, go, mPX)
        null_N.append(dN)
        null_delta.append(dN - dPX)
    null_N, null_delta = np.array(null_N), np.array(null_delta)
    base = dict(pairs_P=int(masks["P"].sum()), pairs_N=int(masks["N"].sum()), pairs_X=int(masks["X"].sum()),
                sd_delta=float(np.std(null_delta, ddof=1)), mde_delta=float(2.49 * np.std(null_delta, ddof=1)),
                sd_D_N=float(np.std(null_N, ddof=1)), mde_D_N=float(2.49 * np.std(null_N, ddof=1)))
    if a.stage == "power":
        pd.Series(base).to_csv(A.OUT / "a327_power.csv", header=["value"])
        print(pd.Series(base).to_string(), f"\ndone {time.time() - t0:.0f}s")
        return
    rows = []
    for c, m in masks.items():
        rows.append(dict(cls=c, pairs=int(m.sum()), rho_pm=A.rho(Z, M, fast, pm, m),
                         rho_od=float(np.median([A.rho(Z, M, fast, d, m) for d in od_draws])),
                         D_rho=A.d_rho(Z, M, fast, pm, od_draws, m)))
    out = pd.DataFrame(rows).set_index("cls")
    delta = out.loc["N", "D_rho"] - out.loc["PX", "D_rho"]
    out["delta_N_minus_PX"] = delta
    out["p_delta_two_sided"] = float((1 + np.sum(np.abs(null_delta) >= abs(delta))) / (1 + N_PERM))
    out["p_D_N"] = float((1 + np.sum(null_N >= out.loc["N", "D_rho"])) / (1 + N_PERM))
    for k, v in base.items():
        out[k] = v
    out.to_csv(A.OUT / "a327_tests.csv")
    print(out.round(3).T.to_string(), f"\ndone {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
