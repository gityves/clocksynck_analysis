"""A339 power stage (registered in `paper/docs/registered_experiments.md`), before any damping ratio rho is computed on
fetal or postnatal data.

rho = |sum_i z_u,i conj(z_s,i)| / (sum_i |z_u,i|^2 - N_u) on the matched estimator (spliced thinned gene by gene to
the unspliced level within cohort; 5 replicates averaged here, 20 in the scoring), N_u the Poisson variance of z_u
(sum over genes of L_kg^2 x Var[log(Poisson(lam) + 0.5)], lam the expected count, `lognoise`; times the residual
degrees-of-freedom factor). 90 % bootstrap interval within
cohort. H339.1: lower bound > 0.43 (3 h at tau = 1 h); H339.2: lower bound > 0.69 (6 h).

Stand-ins and design: those of `power_a328b.py` (fetal designs and cohort clock-gene medians; GTEx samples of the
matching tissue, deaths 08-18 h; Greenwell AL + NR real counts thinned to each unit). Settings: `a338_24h` (unspliced
state turned by phi and amplified by 1 / cos phi, phi from the 1.27 h lead: true rho 0.945), `a338_6h` (true rho
cos(arctan(2 pi tau / 6 h))), `a338_3h` (0.34), `steady` (1.0), `greenwell` (real 24 h clock).
Decidable (per unit and criterion): pass rate >= 0.8 under a338_24h and Greenwell, <= 0.1 under a338_3h (H339.1) or
a338_6h (H339.2).

    python power_a339.py [--n-sub 200] [--jobs 4] [--check]
Outputs (`../output/power_a339/`): subsamples.csv, passrates.csv, decisions.csv (`_check` with --check).
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import power_a328b as P  # noqa: E402  (sets the import paths of the analysis repo)
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector  # noqa: E402
from common import development as dev  # noqa: E402
from score_a337 import boot  # noqa: E402
import lognoise  # noqa: E402

OUT = HERE.parents[1] / "output" / "power_a339"
SEED, N_BOOT, REPS = 3390, 500, 5
PHI_6H = np.arctan(2 * np.pi * P.TAU / 6.0)
P.SETTINGS["a338_6h"] = (1.0, PHI_6H * P.HOURS, 1 / np.cos(PHI_6H))
SETTINGS = ("a338_24h", "a338_6h", "a338_3h", "steady")
TRUE = {"a338_24h": np.cos(P.PHI_REF), "a338_6h": np.cos(PHI_6H), "a338_3h": np.cos(P.PHI_3H), "steady": 1.0}
THRESH = {"H339_1": 0.43, "H339_2": 0.69}
S = {}


def rho_stats(s, u, obs, Z, groups, wf, rng, depth_groups=None) -> dict:
    """Matched-estimator rho with its bootstrap interval."""
    idx = [P.GENES.index(g) for g in wf.genes]
    s, u = s[:, idx], u[:, idx]
    L = projector(pd.DataFrame(s, columns=list(wf.genes)), wf)
    ds, du = obs.depth_s.to_numpy(float), obs.depth_u.to_numpy(float)
    dof = 1 - (1 + Z.shape[1]) / len(s)

    def z(x, depth):
        y = dev.residual(np.log((x + 0.5) / np.broadcast_to(depth, x.shape)), Z)
        b = (y - y.mean(0)) @ L.T
        return b[:, 0] + 1j * b[:, 1]

    zu = z(u, du[:, None])
    zs = np.zeros(len(s), complex)
    for _ in range(REPS):
        sm, dsm = P.depth_match(s, u, ds, du, groups if depth_groups is None else depth_groups, rng)
        zs += z(sm, dsm)
    zs /= REPS
    nu = dof * lognoise.var_log(lognoise.expected(u, du, groups if depth_groups is None else depth_groups)) @ (L ** 2).sum(0)  # counting variance of each z_u
    cross, var_u = zu * np.conj(zs), np.abs(zu) ** 2

    def rho(i):
        den = var_u[i].sum(-1) - nu[i].sum(-1)
        return np.abs(cross[i].sum(-1)) / np.where(den > 0, den, np.nan)

    B = boot(groups, rng, N_BOOT)
    rb = rho(B)
    return dict(rho=float(rho(np.arange(len(s)))), rho_lo=float(np.nanquantile(rb, .05)),
                rho_hi=float(np.nanquantile(rb, .95)), noise_share=float(nu.sum() / var_u.sum()))


def one(task) -> list:
    comp, i = task
    rng = np.random.default_rng([SEED, i, sum(map(ord, comp))])
    obs, ts, tu, pool = S["obs"][comp], S["ts"][comp], S["tu"][comp], S["pools"][S["tissue"][comp]]
    n = min(len(obs), len(pool["depth"]))
    if n < len(obs):
        keep = np.sort(rng.choice(len(obs), n, replace=False))
        obs, ts, tu = obs.iloc[keep].reset_index(drop=True), ts[keep], tu[keep]
    pick = rng.choice(len(pool["depth"]), n, replace=False)
    depth = pool["depth"][pick]
    Z, groups = P.covariates(obs), obs.cohort.to_numpy()
    out = []
    for name in SETTINGS:
        s, _ = P.thin(pool["xs"][name][pick], ts, depth, rng)
        u, _ = P.thin(pool["xu"][name][pick], tu, depth, rng)
        for w, wf in P.WAVEFORMS.items():
            out.append(dict(unit=comp, setting=name, waveform=w, sub=i, n=n, true_rho=TRUE[name],
                            **rho_stats(s, u, obs, Z, groups, wf, rng)))
    g = S["greenwell"]
    m = min(len(obs), len(g["obs"]))
    gp, fp = rng.choice(len(g["obs"]), m, replace=False), rng.choice(len(obs), m, replace=False)
    fo = obs.iloc[fp].reset_index(drop=True)
    s, _ = P.thin(g["s"][gp], ts[fp], g["ds"][gp], rng)
    u, _ = P.thin(g["u"][gp], tu[fp], g["du"][gp], rng)
    go = fo.assign(regimen=g["obs"].group.to_numpy()[gp])
    Zg = dev.design(go, age=None, factors=("regimen", "cohort"))
    for w, wf in P.WAVEFORMS.items():
        out.append(dict(unit=comp, setting="greenwell", waveform=w, sub=i, n=m, true_rho=np.nan,
                        **rho_stats(s, u, go, Zg, (go.regimen + "|" + go.cohort).to_numpy(), wf, rng,
                                    fo.cohort.to_numpy())))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0, sfx = time.time(), "_check" if a.check else ""
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 220)
    design, det = P.prep()
    S["pools"], _ = P.pools(np.random.default_rng([SEED, 0]))
    gw = P.score_a337.load("greenwell")[0]
    S["greenwell"] = dict(obs=gw["obs"], s=gw["spliced"].reindex(columns=P.GENES, fill_value=0).to_numpy(np.int64),
                          u=gw["unspliced"].reindex(columns=P.GENES, fill_value=0).to_numpy(np.int64),
                          ds=gw["spliced"].sum(axis=1).to_numpy(float), du=gw["unspliced"].sum(axis=1).to_numpy(float))
    comps = list(design.compartment.unique())
    S.update(obs={c: design[design.compartment == c].reset_index(drop=True) for c in comps},
             tissue={c: design[design.compartment == c].tissue.iloc[0] for c in comps})
    S["ts"] = {c: P.targets(S["obs"][c], det[det.compartment == c], "spliced") for c in comps}
    S["tu"] = {c: P.targets(S["obs"][c], det[det.compartment == c], "unspliced") for c in comps}
    print(f"prep {time.time() - t0:.0f}s; true rho", {k: round(v, 3) for k, v in TRUE.items()}, flush=True)
    n = 3 if a.check else a.n_sub
    tasks, res = [(c, i) for i in range(n) for c in comps], []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, r in enumerate(pool.imap_unordered(one, tasks)):
            res.extend(r)
            if (k + 1) % 50 == 0 or k + 1 == len(tasks):
                print(f"{k + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res).sort_values(["unit", "setting", "waveform", "sub"])
    res.to_csv(OUT / f"subsamples{sfx}.csv", index=False)
    rows = []
    for (c, s, w), r in res.groupby(["unit", "setting", "waveform"], sort=False):
        rows.append(dict(unit=c, setting=s, waveform=w, n=int(r.n.iloc[0]), true_rho=r.true_rho.iloc[0],
                         median_rho=r.rho.median(), median_lo=r.rho_lo.median(), noise_share=r.noise_share.median(),
                         **{h: (r.rho_lo > t).mean() for h, t in THRESH.items()}))
    rates = pd.DataFrame(rows)
    rates.to_csv(OUT / f"passrates{sfx}.csv", index=False)
    dec = []
    for (c, w), r in rates.groupby(["unit", "waveform"], sort=False):
        v = r.set_index("setting")
        for h, alt in (("H339_1", "a338_3h"), ("H339_2", "a338_6h")):
            sens, gwr, false = v.loc["a338_24h", h], v.loc["greenwell", h], v.loc[alt, h]
            dec.append(dict(unit=c, waveform=w, hypothesis=h, sensitivity=sens, greenwell=gwr, false_rate=false,
                            decidable=bool(sens >= 0.8 and gwr >= 0.8 and false <= 0.1)))
    dec = pd.DataFrame(dec)
    dec.to_csv(OUT / f"decisions{sfx}.csv", index=False)
    print(rates.round(3).to_string(index=False), dec.round(3).to_string(index=False),
          f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
