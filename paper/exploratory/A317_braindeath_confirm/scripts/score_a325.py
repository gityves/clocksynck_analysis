"""A325 (registered, check B2): does post-mortem ischaemia shape the reference organ arrangement?

Within fast deaths (H325.1) and all post-mortem deaths (H325.2), the top and bottom tertiles of donor ischaemic time
(median over the donor's samples), cut within hour-of-death strata (08-20 h, 20-08 h; x death group for H325.2), are
compared by the Fisher-Lee rho of their pair phases (A317 estimator), against random disjoint splits of the same sizes
within the same strata. Power stage: the organ-donor minus post-mortem tissue shifts injected into the top group of
random splits (f = 1, 0.5). Aggregates only (GTEx controlled access).

    python score_a325.py --stage power | score
Outputs (`../output/`): a325_power.csv | a325_tests.csv, a325_sets.csv, a325_hour.csv
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as SS  # noqa: E402

import score_a317 as A  # noqa: E402
from score_a317 import B, H, RS, gtex  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = A.OUT
N_NULL, N_INJ, SEED = 2000, 200, 3251
SETS = list(range(0, 200, 4))


def tertiles(isch, strata):
    """Top / bottom tertile donors, cut within each stratum (strata: dict label -> donor array)."""
    top, bot = [], []
    for d in strata.values():
        x = isch.reindex(d).dropna()
        if len(x) < 3:
            continue
        lo, hi = x.quantile([1 / 3, 2 / 3])
        top.append(x[x > hi].index.to_numpy())
        bot.append(x[x <= lo].index.to_numpy())
    return np.concatenate(top), np.concatenate(bot)


def random_split(strata, top, bot, rng):
    """Disjoint random groups with the observed per-stratum sizes of top and bottom."""
    gt, gb = [], []
    for d in strata.values():
        nt, nb = np.isin(d, top).sum(), np.isin(d, bot).sum()
        p = rng.permutation(d)
        gt.append(p[:nt])
        gb.append(p[nt:nt + nb])
    return np.concatenate(gt), np.concatenate(gb)


def tissue_shifts(Z, M, od, pm, mask):
    """Least-squares tissue shifts s (rad, sum 0) with s_a - s_b = wrapped (od - pm) pair phase difference."""
    d = np.angle(np.exp(1j * (B.phase(Z, M, od, mask) - B.phase(Z, M, pm, mask))))
    X = np.zeros((int(mask.sum()), len(H.BODY)))
    for r, i in enumerate(np.flatnonzero(mask)):
        a, b = H.PAIRS[i]
        X[r, a], X[r, b] = 1, -1
    return np.linalg.lstsq(X, d, rcond=None)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    if a.stage == "score":
        RS.add_sets(T, 200, 0)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h3, h4, h0 = (T.group_donors(g) for g in ("fast", "H3", "H4", "vent"))
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    attrs = A.sample_attrs(T)
    every = np.concatenate([fast, h3, h4, pm])
    isch = pd.Series({d: np.nanmedian([attrs[t].ischemic_time.to_numpy(float)[T.R[t][d]] for t in H.BODY
                                       if T.R[t][d] >= 0] or [np.nan]) for d in every})
    day = (T.tod >= 8) & (T.tod < 20)
    designs = {
        "H325.1": {f"{s}": fast[day[fast] == (s == "day")] for s in ("day", "night")},
        "H325.2": {f"{g}_{s}": d[day[d] == (s == "day")] for g, d in (("fast", fast), ("H3", h3), ("H4", h4), ("pm", pm))
                   for s in ("day", "night")},
        "H325.1_unstratified": {"all": fast},
    }
    Z, M = A.raw(T, "clock", 0, T.R)
    rng = np.random.default_rng(SEED)

    # hour-of-death association (aggregates)
    hr = []
    for g, d in (("fast", fast), ("all_postmortem", every)):
        x, w = isch.reindex(d).to_numpy(), 2 * np.pi * T.tod[d] / 24
        ok = np.isfinite(x)
        hr.append(dict(group=g, n=int(ok.sum()), spearman_cos=SS.spearmanr(x[ok], np.cos(w[ok]))[0],
                       spearman_sin=SS.spearmanr(x[ok], np.sin(w[ok]))[0],
                       isch_median_day=float(np.nanmedian(x[day[d]])), isch_median_night=float(np.nanmedian(x[~day[d]])),
                       n_day=int(day[d].sum()), n_night=int((~day[d]).sum())))
    pd.DataFrame(hr).to_csv(OUT / "a325_hour.csv", index=False)

    rows, setrows = [], []
    od_draw0 = np.random.default_rng(A.SEED).choice(od, len(pm), replace=False)
    m65 = B.pair_mask(M, fast, pm, od_draw0) & B.pair_mask(M, od)
    s = tissue_shifts(Z, M, od, pm, m65)
    for name, strata in designs.items():
        top, bot = tertiles(isch, strata)
        mask = B.pair_mask(M, top, bot)
        null = []
        splits = [random_split(strata, top, bot, rng) for _ in range(N_NULL)]
        null = np.array([A.rho(Z, M, gt, gb, mask) for gt, gb in splits])
        q05 = float(np.quantile(null, 0.05))
        row = dict(test=name, n_top=len(top), n_bottom=len(bot), pairs=int(mask.sum()),
                   isch_top_median=float(isch.reindex(top).median()), isch_bottom_median=float(isch.reindex(bot).median()),
                   null_median=float(np.median(null)), null_q05=q05)
        if a.stage == "power":
            for f in (1.0, 0.5):
                hits = []
                for gt, gb in splits[:N_INJ]:
                    Zi = Z.copy()
                    Zi[gt] = Zi[gt] * np.exp(1j * f * s)[None, :]
                    hits.append(A.rho(Zi, M, gt, gb, mask) < q05)
                row[f"power_f{f:g}"] = float(np.mean(hits))
            row["shift_sd_h"] = float(np.std(s) * H.HOURS)
        else:
            r = A.rho(Z, M, top, bot, mask)
            row.update(rho_tb=r, p=float((1 + np.sum(null <= r)) / (1 + len(null))))
            if name == "H325.1":
                for kind, lab in (("pseudo", "matched"), ("resp", "responsive")):
                    for k in SETS:
                        Zk, Mk = A.raw(T, kind, k, T.R)
                        rk = A.rho(Zk, Mk, top, bot, mask)
                        nk = np.array([A.rho(Zk, Mk, gt, gb, mask) for gt, gb in splits])
                        setrows.append(dict(family=lab, set=k, rho_tb=rk, p=float((1 + np.sum(nk <= rk)) / (1 + len(nk)))))
                sr = pd.DataFrame(setrows)
                for lab in ("matched", "responsive"):
                    row[f"{lab}_p_q05"] = float(sr[sr.family == lab].p.quantile(0.05))
                    row[f"{lab}_share_p_below"] = float(np.mean(sr[sr.family == lab].p <= row["p"]))
        rows.append(row)
        print(pd.Series(row).to_string(), f"\n({time.time() - t0:.0f}s)", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / ("a325_power.csv" if a.stage == "power" else "a325_tests.csv"), index=False)
    if setrows:
        pd.DataFrame(setrows).to_csv(OUT / "a325_sets.csv", index=False)
    print(pd.DataFrame(hr).round(3).to_string(index=False))
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
