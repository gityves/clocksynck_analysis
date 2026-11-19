"""Living reference for brain death (prompt cell_route_2026-11-19, Task 3): precondition, aggregate run before registering.

Pair-level offsets (phase of E[z_a conj z_b] over shared donors, A317 estimator, in hours; positive = first tissue's
clock-plane angle larger) of the GTEx pairs that living surgical cohorts also carry, per group: fast deaths, Hardy 4,
ventilated post-mortem, organ donors (median over the A317 draws of 79, and all 298). Donor bootstrap (500, seed 3250)
for the organ-donor minus post-mortem difference. Aggregates only (GTEx controlled access).

    python living_ref_precondition.py
Output: ../output/living_ref_precondition.csv
"""
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

SUB, VIS, MUS = "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Muscle - Skeletal"
SKL, SKS = "Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"
PAIRS = [(SUB, VIS), (MUS, SUB), (SUB, SKL), (SUB, SKS)]
N_BOOT, SEED = 500, 3250


def offset_h(Z, M, idx, a, b):
    d = idx[(M[idx, a] > 0) & (M[idx, b] > 0)]
    Zc = (Z[idx] - (Z[idx] * M[idx]).sum(0) / np.maximum(M[idx].sum(0), 1)) * M[idx]
    za, zb = Zc[np.isin(idx, d), a], Zc[np.isin(idx, d), b]
    return float(np.angle(np.sum(za * zb.conj())) * H.HOURS), len(d)


def wrap(x):
    return (x + 12) % 24 - 12


def main():
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h4, h0 = T.group_donors("fast"), T.group_donors("H4"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    rng = np.random.default_rng(A.SEED)
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(A.N_DRAWS)]
    Z, M = A.raw(T, "clock", 0, T.R)
    j = {t: i for i, t in enumerate(H.BODY)}
    groups = {"fast": fast, "hardy4": h4, "postmortem_vent": pm, "organ_donors_all": od}
    brng = np.random.default_rng(SEED)
    rows = []
    for a, b in PAIRS:
        ia, ib = j[a], j[b]
        row = dict(pair=f"{a} | {b}")
        for g, idx in groups.items():
            row[f"{g}_h"], row[f"{g}_n"] = offset_h(Z, M, idx, ia, ib)
        dd = [offset_h(Z, M, d, ia, ib)[0] for d in od_draws]
        row["organ_donors_draws_h"] = float(np.angle(np.mean(np.exp(1j * np.array(dd) / H.HOURS))) * H.HOURS)
        row["od_minus_pm_h"] = wrap(row["organ_donors_draws_h"] - row["postmortem_vent_h"])
        row["od_minus_fast_h"] = wrap(row["organ_donors_all_h"] - row["fast_h"])
        row["pm_minus_fast_h"] = wrap(row["postmortem_vent_h"] - row["fast_h"])
        boot = []
        for _ in range(N_BOOT):
            bo, bp = brng.choice(od, len(od)), brng.choice(pm, len(pm))
            boot.append(wrap(offset_h(Z, M, bo, ia, ib)[0] - offset_h(Z, M, bp, ia, ib)[0]))
        boot = np.array(boot)
        c = np.angle(np.mean(np.exp(1j * boot / H.HOURS))) * H.HOURS
        dev = wrap(boot - c)
        row["od_all_minus_pm_h"] = wrap(row["organ_donors_all_h"] - row["postmortem_vent_h"])
        row["boot_lo"], row["boot_hi"] = float(c + np.quantile(dev, 0.025)), float(c + np.quantile(dev, 0.975))
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(A.OUT / "living_ref_precondition.csv", index=False)
    print(out.round(2).T.to_string())
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
