"""A349 part 1 (registered, power only): detection limit of a change in the organ order at the ventilated design.

As A336 part 3 (`a336_detectability.order_sim`): the arm under test has one organ's clock states rotated by x h (each
organ in turn; 0.25 h grid to 12 h) or every organ by N(0, x^2); detection = Fisher-Lee rho(full fast-death map, arm map)
below the 5th percentile of the unrotated draws; limit = smallest x with detection >= 0.8. 200 draws, pairs with >= 20
donors in the fast map and the arm of each draw.
  ventilated arm : the 79 ventilated post-mortem donors resampled with replacement.
  companion arm  : 79 fast deaths drawn without replacement; reference = the remaining fast deaths (the registered
                   "full fast map" cannot contain the test donors).
Gene sets: 14 genes and core loop. Aggregates only (no donor-level output).

    python a349_part1.py [--check] [--n-draw 200]
Outputs (`../output/`): part1_sim.csv, part1_limits.csv
"""
import argparse
import time

import numpy as np
import pandas as pd

import a349_common as C
from a349_common import B, H, S
import compute_coupling_nulls as N
from clocksync import stats as ST
from clocksync.config import limit_threads

W = 2 * np.pi / 24
GRID = np.arange(0.25, 12.01, 0.25)


def theta(Z, M, idx, mask):
    return H.pair_phase(N._H_of(Z, M, idx))[mask]


def sim(Z, M, fast, pm, arm, n_draw, seed):
    rng = np.random.default_rng(seed)
    n = len(pm)
    draws = []
    for _ in range(n_draw):
        if arm == "ventilated":
            ref, g = fast, rng.choice(pm, n)
        else:
            p = rng.permutation(fast)
            g, ref = p[:n], p[n:]
        mask = B.pair_mask(M, ref, g)
        draws.append((ref, g, mask, H.pair_phase(N._H_of(Z, M, ref))[mask]))
    rho = lambda Zx, d: float(ST.fisher_lee(d[3], theta(Zx, M, d[1], d[2])))  # noqa: E731
    base = np.array([rho(Z, d) for d in draws])
    thr = float(np.quantile(base, 0.05))
    rows = []
    for x in GRID:
        rot = np.exp(1j * x * W)
        for o, t in enumerate(H.BODY):
            det = []
            for d in draws:
                Zs = Z.copy()
                Zs[d[1], o] *= rot
                det.append(rho(Zs, d) < thr)
            rows.append(dict(design="one organ shifted", organ=t, x_h=x, detection=float(np.mean(det))))
        det = []
        for d in draws:
            Zs = Z.copy()
            Zs[d[1], :] *= np.exp(1j * rng.normal(0, x * W, len(H.BODY)))[None, :]
            det.append(rho(Zs, d) < thr)
        rows.append(dict(design="organ phases dispersed (SD)", organ="all", x_h=x, detection=float(np.mean(det))))
    return pd.DataFrame(rows), dict(base_median=float(np.median(base)), threshold=thr,
                                    pairs_median=float(np.median([d[2].sum() for d in draws])), n_draw=n_draw)


def limits(s, design):
    out = {}
    for o, g in s[s.design == design].groupby("organ"):
        ok = g[g.detection >= 0.8]
        out[o] = float(ok.x_h.min()) if len(ok) else np.nan
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="5 draws, coarse output to scratch; not a run")
    ap.add_argument("--n-draw", type=int, default=200)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    out = C.OUT if not a.check else C.OUT / "check"
    out.mkdir(parents=True, exist_ok=True)
    T, fast, pm = C.tissues()
    print(f"states built ({time.time() - t0:.0f}s): fast {len(fast)}, ventilated post-mortem {len(pm)}", flush=True)
    sims, lims = [], []
    for gs, k in (("14 genes", "clock"), ("core", "core")):
        Z, M = S.raw(T, k, 0, T.R)
        for arm in ("ventilated", "fast companion"):
            s, att = sim(Z, M, fast, pm, "ventilated" if arm == "ventilated" else "fast", 5 if a.check else a.n_draw,
                         [349, 1, C.SEEDS[gs], arm == "ventilated"])
            s.insert(0, "arm", arm), s.insert(0, "gene_set", gs)
            sims.append(s)
            one, disp = limits(s, "one organ shifted"), limits(s, "organ phases dispersed (SD)")["all"]
            v = np.array(list(one.values()), float)
            lims.append(dict(gene_set=gs, arm=arm, n_arm=len(pm), **att, one_organ_median_h=float(np.nanmedian(v)) if np.isfinite(v).any() else np.nan,
                             one_organ_min_h=float(np.nanmin(v)) if np.isfinite(v).any() else np.nan,
                             one_organ_max_h=float(np.nanmax(v)) if np.isfinite(v).any() else np.nan,
                             n_organs_never=int((~np.isfinite(v)).sum()), dispersion_sd_h=disp,
                             **{f"limit_{t}": one[t] for t in H.BODY}))
            print(f"{gs} / {arm} done ({time.time() - t0:.0f}s)", flush=True)
            pd.concat(sims).to_csv(out / "part1_sim.csv", index=False)
            pd.DataFrame(lims).to_csv(out / "part1_limits.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(lims).drop(columns=[c for c in pd.DataFrame(lims).columns if c.startswith("limit_")]).round(3).to_string(index=False))
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
