"""A317 post hoc (2026-11-25, after the slide 4.2c scatter showed the per-tissue agreement of ventilated donors resting
on the lung): leave-one-tissue-out of H317.1. For each of the 12 body tissues, the pairs involving it are removed
(states are per-tissue projections, no host, so this removes the tissue entirely) and D_rho = rho(post-mortem, fast) -
median rho(organ-donor draws, fast) is recomputed with the A317 draws, its donor-label permutation p (1000) and the
matched / responsive 95th percentiles (200 sets each). The full-mask row must reproduce A317 (D_rho 0.621).
Descriptive, not registered. Aggregates only.

    python posthoc_loto.py
Output: ../output/posthoc_loto.csv
"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a317 as A  # noqa: E402
from score_a317 import B, H, RS, gtex, rho, d_rho, raw  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

N_PERM = 1000


def main():
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    RS.add_sets(T, 200, 0)
    coh = gtex.phenotypes(B.PH_COLS).reindex(T.donors).COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    rng = np.random.default_rng(A.SEED)  # same draws as A317
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(A.N_DRAWS)]
    Zc, Mc = raw(T, "clock", 0, T.R)
    full = B.pair_mask(Mc, fast, pm, *od_draws[:1]) & B.pair_mask(Mc, od)
    sets = {(kind, k): raw(T, kind, k, T.R) for kind in ("pseudo", "resp") for k in range(200)}
    print(f"states built ({time.time() - t0:.0f}s)", flush=True)
    rows = []
    for drop in [None, *H.BODY]:
        j = None if drop is None else H.BODY.index(drop)
        mask = full & np.array([j not in p for p in H.PAIRS])
        D = d_rho(Zc, Mc, fast, pm, od_draws, mask)
        null = []
        for _ in range(N_PERM):
            p = rng.permutation(h0)
            null.append(rho(Zc, Mc, fast, p[:len(pm)], mask) - rho(Zc, Mc, fast, p[len(pm):2 * len(pm)], mask))
        comp = {kind: np.array([d_rho(*sets[kind, k], fast, pm, od_draws, mask) for k in range(200)])
                for kind in ("pseudo", "resp")}
        rows.append(dict(dropped=drop or "none", pairs=int(mask.sum()), D_rho=D,
                         rho_pm=rho(Zc, Mc, fast, pm, mask),
                         rho_od=float(np.median([rho(Zc, Mc, fast, d, mask) for d in od_draws])),
                         p_labels=float((1 + np.sum(np.array(null) >= D)) / (1 + N_PERM)),
                         matched_q95=float(np.nanquantile(comp["pseudo"], 0.95)),
                         responsive_q95=float(np.nanquantile(comp["resp"], 0.95))))
        print(pd.Series(rows[-1]).to_dict(), f"{time.time() - t0:.0f}s", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(A.OUT / "posthoc_loto.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
