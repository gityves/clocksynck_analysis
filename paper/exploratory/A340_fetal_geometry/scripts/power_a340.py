"""A340 power stage (registered in `paper/docs/registered_experiments.md`), before any fetal improperness iota is computed.

iota = |sum_i (z_i^2 - nu_i)| / sum_i (|z_i|^2 - t_i): 0 for a clock-plane cloud with no preferred axis, 1 for points
on a line; z_i the spliced clock-plane state (log fractions residualised as A328 part b, one projector), nu_i and t_i
its counting pseudo-variance and variance (exact log-count noise, `lognoise`).
Noise matching (after the score code check, 2026-11-22, before any fetal iota): genes permuted independently gave iota
0.56-0.93, because each gene's biological variance projected through the unequal projector weights is itself
elongated. Every reference draw is therefore given independent per-gene log-normal noise up to the fetal unit's
per-gene residual variance of the spliced log fractions (`match_noise`; the fetal clock signal is part of that
variance, so the added noise is an upper bound and the separation conservative); a no-clock reference (genes permuted
across the drawn samples) shows where pure noise falls.

Reference distributions per unit (del Valle, prenatal liver; 10x cortex descriptive): n distinct GTEx samples of the
matching tissue, fast deaths (Hardy 1-2) at 08-18 h (time of death in memory only), clock at 1, 0.75 and 0.5 of the
adult size (`simulate.scale_clock`), thinned to the unit's spliced clock-gene depths (`power_a328b` design and
targets): iota_daily. The same draws with every clock-plane state replaced by its projection on the draw's main axis,
total clock-plane variance kept (`simulate.thin_clock`), then thinned: iota_line. Decidable if q95(iota_daily) <
q05(iota_line) at 0.75.

    python power_a340.py [--n-sub 1000] [--jobs 4] [--check]
Outputs (`../output/power/`): pools.csv, subsamples.csv, quantiles.csv, decisions.csv (`_check` with --check).
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "A328_fetal_running_loop" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import power_a328b as P  # noqa: E402  (sets the import paths of the analysis repo)
import lognoise  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from clocksync.simulate import scale_clock, thin_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402

OUT = HERE.parents[1] / "output" / "power"
SEED = 3400
SCALES = (1.0, 0.75, 0.5)
UNITS = ("delvalle", "liver", "cortex")
S = {}


def iota(s, depth, Z, wf, groups) -> float:
    """Noise-corrected improperness of the spliced clock-plane states (s: n x genes of wf, raw counts)."""
    L = projector(pd.DataFrame(s, columns=list(wf.genes)), wf)
    y = dev.residual(np.log((s + 0.5) / depth[:, None]), Z)
    b = (y - y.mean(0)) @ L.T
    z = b[:, 0] + 1j * b[:, 1]
    dof = 1 - (1 + Z.shape[1]) / len(s)
    w = dof * lognoise.var_log(lognoise.expected(s, depth, groups))  # n x G counting variance of log(s + 0.5)
    V = np.einsum("ag,ng,bg->nab", L, w, L)              # n x 2 x 2
    t = V[:, 0, 0] + V[:, 1, 1]
    nu = V[:, 0, 0] - V[:, 1, 1] + 2j * V[:, 0, 1]
    return float(np.abs(np.sum(z ** 2) - nu.sum()) / (np.sum(np.abs(z) ** 2) - t.sum()))


def gene_var(s, depth, Z) -> np.ndarray:
    """Per-gene variance of the residualised log fractions (all 14 genes)."""
    y = dev.residual(np.log((s + 0.5) / depth[:, None]), Z)
    return y.var(0, ddof=1 + Z.shape[1])


def match_noise(s, depth, Z, target, rng):
    """Counts with independent log-normal per-gene noise added until each gene's residual variance reaches the fetal
    unit's (`target`); genes already above it are left as they are."""
    extra = np.clip(target - gene_var(s, depth, Z), 0, None)
    eps = rng.normal(0, np.sqrt(extra), s.shape) - extra / 2
    return rng.poisson(s * np.exp(eps))


def line(x, depth, rng):
    """Counts whose centred clock-plane states are projected on their main axis (total variance kept)."""
    df = pd.DataFrame(x, columns=P.GENES).assign(depth=depth)
    beta, _ = states(df, projector(df, P.WF), P.WF, depth)
    bc = beta - beta.mean(0)
    ev, evec = np.linalg.eigh(bc.T @ bc)
    u = evec[:, -1]
    bl = np.outer(bc @ u, u) * np.sqrt(ev.sum() / ev[-1])
    return thin_clock(x.astype(float), P.WF.X, bl - bc, rng).astype(np.int64)


def one(task) -> list:
    unit, i = task
    rng = np.random.default_rng([SEED, i, sum(map(ord, unit))])
    obs, ts, pool = S["obs"][unit], S["ts"][unit], S["pools"][S["tissue"][unit]]
    n = min(len(obs), len(pool["depth"]))
    if n < len(obs):
        keep = np.sort(rng.choice(len(obs), n, replace=False))
        obs, ts = obs.iloc[keep].reset_index(drop=True), ts[keep]
    pick = rng.choice(len(pool["depth"]), n, replace=False)
    depth = pool["depth"][pick]
    Z = P.covariates(obs)
    ds = obs.depth_s.to_numpy(float)
    out = []
    for a in SCALES:
        x = pool["x"][a][pick]
        null = np.column_stack([rng.permutation(c) for c in x.T])  # no between-gene covariance
        for kind, xx in (("daily", x), ("line", line(x, depth, rng)), ("noclock", null)):
            s, _ = P.thin(xx, ts, depth, rng)
            s = match_noise(s, ds, Z, S["fvar"][unit], rng)
            for w, wf in P.WAVEFORMS.items():
                idx = [P.GENES.index(g) for g in wf.genes]
                out.append(dict(unit=unit, scale=a, kind=kind, waveform=w, sub=i, n=n,
                                iota=iota(s[:, idx], ds, Z, wf, obs.cohort.to_numpy())))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0, sfx = time.time(), "_check" if a.check else ""
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 220)
    design, det = P.prep()
    rng = np.random.default_rng([SEED, 0])
    pools, prow = {}, []
    for key in ("adrenal", "liver"):
        pb, d = gtex.load([P.TISSUE[key]], hardy=gtex.FAST, genes="clock", depth="library_size")
        pb = pb[pb.donor.isin(set(d.index[in_window(d.tod.to_numpy())]))].reset_index(drop=True)
        clock = pb[["donor", "compartment", "depth", *P.GENES]].copy()
        pools[key] = dict(depth=clock.depth.to_numpy(float),
                          x={s: scale_clock(clock, P.WF, s, depth=clock.depth, rng=rng)[P.GENES].to_numpy(np.int64)
                             for s in SCALES})
        prow.append(dict(tissue=P.TISSUE[key], fast_deaths_08_18=len(pb)))
    pd.DataFrame(prow).to_csv(OUT / f"pools{sfx}.csv", index=False)
    print(pd.DataFrame(prow).to_string(index=False), flush=True)
    S.update(pools=pools, obs={u: design[design.compartment == u].reset_index(drop=True) for u in UNITS},
             tissue={u: design[design.compartment == u].tissue.iloc[0] for u in UNITS})
    S["ts"] = {u: P.targets(S["obs"][u], det[det.compartment == u], "spliced") for u in UNITS}
    import score_a328b as B
    fv = {}
    for u in B.load_units():
        if u["name"] in UNITS:
            s = u["s"][P.GENES].to_numpy(np.int64)
            fv[u["name"]] = gene_var(s, u["obs"].depth_s.to_numpy(float), B.covariates(u["obs"]))
    S["fvar"] = fv
    print("fetal per-gene residual variance (median over genes):", {k: round(float(np.median(v)), 3) for k, v in fv.items()},
          flush=True)
    n = 5 if a.check else a.n_sub
    tasks, res = [(u, i) for i in range(n) for u in UNITS], []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(one, tasks)):
            res.extend(r)
            if (k + 1) % 300 == 0 or k + 1 == len(tasks):
                print(f"{k + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res)
    res.to_csv(OUT / f"subsamples{sfx}.csv", index=False)
    q = (res.groupby(["unit", "waveform", "scale", "kind"]).iota.quantile([.05, .5, .95]).unstack()
         .rename(columns={.05: "q05", .5: "median", .95: "q95"}).reset_index())
    q.to_csv(OUT / f"quantiles{sfx}.csv", index=False)
    dec = []
    for (u, w, s), r in q.groupby(["unit", "waveform", "scale"]):
        v = r.set_index("kind")
        dec.append(dict(unit=u, waveform=w, scale=s, daily_q95=v.loc["daily", "q95"], line_q05=v.loc["line", "q05"],
                        separated=bool(v.loc["daily", "q95"] < v.loc["line", "q05"])))
    dec = pd.DataFrame(dec)
    dec.to_csv(OUT / f"decisions{sfx}.csv", index=False)
    print(q.round(3).to_string(index=False), dec.round(3).to_string(index=False),
          f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
