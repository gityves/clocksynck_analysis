"""A318 post hoc note (2026-11-19): same-pairs decomposition of A317's common-ischaemic-range sensitivity.

A317's `ischaemic_common_range` row gave D_rho 0.239 on 36 pairs against 0.621 on 65: a change of pair set and of donor
set at once. On the same 36 pairs this reports, for post-mortem ventilated donors and organ donors (median of 20 draws),
the Fisher-Lee rho with fast deaths and D_rho:
  U36  unrestricted samples (A317 main draws; also the sensitivity's draws),
  R36  restricted samples, the sensitivity's draws (reproduces 0.239),
  S36  unrestricted samples subsampled per pair and group to R36's shared-donor counts (20 draws, seed 3180),
and the same on the 29 dropped pairs (U29). Descriptive, post hoc. Aggregates only (GTEx controlled access).

    python same_pairs_a318.py
Outputs (`../output/`): same_pairs.csv (one row per readout and draw), same_pairs_summary.csv, dropped_pairs.csv.
"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a317 as A  # noqa: E402  (sets the A109 / A132 paths)
from score_a317 import B, H, RS, ST, gtex  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = A.OUT
N_SUB, SEED_SUB = 20, 3180


def replay_draws(T, ph, fast, pm, od, h0):
    """Re-create A317's random stream: main organ-donor draws and the ischaemic sensitivity's draws."""
    rng = np.random.default_rng(A.SEED)
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(A.N_DRAWS)]
    for _ in range(A.N_PERM):  # perm_null
        rng.permutation(h0)
    steroid = pd.to_numeric(ph.MHSTRDLT, errors="coerce").to_numpy() == 1
    pm_s, od_s = pm[~steroid[pm]], od[~steroid[od]]
    for _ in range(A.N_DRAWS):  # steroids_excluded sensitivity
        rng.choice(od_s, min(len(pm_s), len(od_s)), replace=False)
    dr = [rng.choice(od, min(len(pm), len(od)), replace=False) for _ in range(A.N_DRAWS)]
    return od_draws, dr


def pair_phases_sub(Z, M, idx, pairs, n_take=None, rng=None):
    """Pair phases over each pair's shared donors among idx (optionally subsampled to n_take[p]); tissues centred on them."""
    out = np.empty(len(pairs))
    for i, (a, b) in enumerate(pairs):
        d = idx[(M[idx, a] > 0) & (M[idx, b] > 0)]
        if n_take is not None:
            d = rng.choice(d, int(n_take[i]), replace=False)
        za, zb = Z[d, a] - Z[d, a].mean(), Z[d, b] - Z[d, b].mean()
        out[i] = np.angle(np.sum(za * zb.conj()))
    return out


def shared(M, idx, pairs):
    return np.array([int(((M[idx, a] > 0) & (M[idx, b] > 0)).sum()) for a, b in pairs])


def readout(name, draw, th_fast, th_pm, th_od):
    r_pm = float(ST.fisher_lee(th_fast, th_pm))
    r_od = float(np.median([ST.fisher_lee(th_fast, t) for t in th_od]))
    return dict(readout=name, draw=draw, rho_pm=r_pm, rho_od=r_od, D_rho=r_pm - r_od)


def main():
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    RS.add_sets(T, 200, 0)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    od_draws, dr = replay_draws(T, ph, fast, pm, od, h0)

    Zu, Mu = A.raw(T, "clock", 0, T.R)
    mask65 = B.pair_mask(Mu, fast, pm, *od_draws[:1]) & B.pair_mask(Mu, od)
    attrs = A.sample_attrs(T)
    isch_pm = np.concatenate([attrs[t].ischemic_time.to_numpy(float)[T.R[t][pm][T.R[t][pm] >= 0]] for t in H.BODY])
    lo, hi = np.nanquantile(isch_pm, [0.05, 0.95])
    rows = {t: (attrs[t].ischemic_time.to_numpy(float) >= lo) & (attrs[t].ischemic_time.to_numpy(float) <= hi)
            for t in H.BODY}
    Rr = A.restrict(T, rows)
    Zr, Mr = A.raw(T, "clock", 0, Rr)
    m36 = mask65 & B.pair_mask(Mr, fast, pm, dr[0])
    m29 = mask65 & ~m36
    print(f"states ({time.time() - t0:.0f}s): pairs 65={mask65.sum()} 36={m36.sum()} 29={m29.sum()}; "
          f"ischaemic range {lo:.0f}-{hi:.0f} min", flush=True)

    res = []
    # A317 estimator (group-centred) on fixed pair sets
    def std(name, Z, M, draws, m):
        th = lambda g: B.phase(Z, M, g, m)  # noqa: E731
        res.append(readout(name, -1, th(fast), th(pm), [th(d) for d in draws]))
    std("R36", Zr, Mr, dr, m36)
    std("U36", Zu, Mu, od_draws, m36)
    std("U36_sens_draws", Zu, Mu, dr, m36)
    std("U65", Zu, Mu, od_draws, mask65)
    std("U29", Zu, Mu, od_draws, m29)
    std("R29_available", Zr, Mr, dr, m29)
    # per-pair estimator (shared-donor-centred): check against U36, then noise-matched subsamples
    pairs36 = [H.PAIRS[i] for i in np.flatnonzero(m36)]
    th = lambda g, Z, M, n=None, rng=None: pair_phases_sub(Z, M, g, pairs36, n, rng)  # noqa: E731
    res.append(readout("R36_pairwise", -1, th(fast, Zr, Mr), th(pm, Zr, Mr), [th(d, Zr, Mr) for d in dr]))
    res.append(readout("U36_pairwise", -1, th(fast, Zu, Mu), th(pm, Zu, Mu), [th(d, Zu, Mu) for d in dr]))
    n_fast, n_pm = shared(Mr, fast, pairs36), shared(Mr, pm, pairs36)
    n_od = [shared(Mr, d, pairs36) for d in dr]
    rng = np.random.default_rng(SEED_SUB)
    for k in range(N_SUB):
        res.append(readout("S36", k, th(fast, Zu, Mu, n_fast, rng), th(pm, Zu, Mu, n_pm, rng),
                           [th(d, Zu, Mu, n, rng) for d, n in zip(dr, n_od)]))
    res = pd.DataFrame(res)
    res.to_csv(OUT / "same_pairs.csv", index=False)

    s = res[res.readout == "S36"]
    r36 = res.set_index("readout").loc["R36"]
    summ = res[res.readout != "S36"].drop(columns="draw").copy()
    summ = pd.concat([summ, pd.DataFrame([dict(readout=f"S36_{q}", **s[["rho_pm", "rho_od", "D_rho"]].quantile(v).to_dict())
                                          for q, v in (("q05", 0.05), ("median", 0.5), ("q95", 0.95))])])
    summ["R36_pct_in_S36"] = float(np.mean(s.D_rho < r36.D_rho))
    summ.to_csv(OUT / "same_pairs_summary.csv", index=False)

    names = [f"{H.BODY[a]} | {H.BODY[b]}" for a, b in H.PAIRS]
    drop = []
    for i in np.flatnonzero(mask65):
        a, b = H.PAIRS[i]
        drop.append(dict(pair=names[i], kept=bool(m36[i]),
                         **{f"n_{g}_{w}": int(((M[idx, a] > 0) & (M[idx, b] > 0)).sum())
                            for w, M in (("all", Mu), ("restricted", Mr))
                            for g, idx in (("fast", fast), ("pm", pm), ("od_draw0", dr[0]))}))
    pd.DataFrame(drop).to_csv(OUT / "dropped_pairs.csv", index=False)
    print(summ.round(3).to_string(index=False))
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
