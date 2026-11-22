"""A339 scoring (registered in `paper/docs/registered_experiments.md`, A339 and its power note): the damping ratio rho
of the fetal clock-plane states. One run.

rho = |sum_i z_u,i conj(z_s,i)| / (sum_i |z_u,i|^2 - N_u) on the matched estimator (20 thinning replicates averaged),
N_u the counting variance of z_u (exact Var[log(Poisson(lam) + 0.5)] at the expected count, `lognoise`, times the
residual degrees-of-freedom factor). Units and residualisation of A328 part b (`score_a328b.load_units`). 90 %
interval and one-sided p (share of bootstrap rho <= threshold) from 2,000 resamples within cohort.
H339.1: rho > 0.43 (periods <= 3 h excluded), H339.2: rho > 0.69 (<= 6 h); scored on the decidable readouts of the
power note (H339.1: del Valle and liver, 14 genes and core; H339.2: del Valle and liver, 14 genes), Holm over the two
units per hypothesis and gene set; the rest descriptive. `--check`: unspliced := thinned spliced (true rho about 1).

    python score_a339.py [--jobs 4] [--check]
Outputs (`../output/score_a339[_check]/`): rho.csv, verdicts.csv.
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

import score_a328b as B  # noqa: E402  (sets the import paths of the analysis repo)
import lognoise  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector  # noqa: E402
from common import development as dev  # noqa: E402
from score_a337 import boot  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, N_BOOT, REPS = 3391, 2000, 20
THRESH = {"H339_1": 0.43, "H339_2": 0.69}
SCORED = {("H339_1", "all"), ("H339_1", "core"), ("H339_2", "all")}
S = {}


def task(k):
    unit, wname, which = S["tasks"][k]
    wf, obs = B.WAVEFORMS[wname], unit["obs"]
    genes = list(wf.genes)
    rng = np.random.default_rng([SEED, k])
    s, u, ds, du = B.counts(unit, genes, which)
    Z, groups = B.covariates(obs), obs.cohort.to_numpy()
    L = projector(pd.DataFrame(s, columns=genes), wf)
    to_z = lambda y: (y @ L.T)[:, 0] + 1j * (y @ L.T)[:, 1]  # noqa: E731
    zu = to_z(B.logfrac(u, du[:, None], Z))
    zs = np.zeros(len(s), complex)
    for _ in range(REPS):
        sm, dsm = B.depth_match(s, u, ds, du, groups, rng)
        zs += to_z(B.logfrac(sm, dsm, Z))
    zs /= REPS
    dof = 1 - (1 + Z.shape[1]) / len(s)
    nu = dof * lognoise.var_log(lognoise.expected(u, du, groups)) @ (L ** 2).sum(0)
    cross, var_u = zu * np.conj(zs), np.abs(zu) ** 2

    def rho(i):
        den = var_u[i].sum(-1) - nu[i].sum(-1)
        return np.abs(cross[i].sum(-1)) / np.where(den > 0, den, np.nan)

    rb = rho(boot(groups, rng, N_BOOT))
    row = dict(unit=unit["name"], role=unit["role"], waveform=wname, unspliced="intron" if which == "u" else "direct",
               n=len(obs), rho=float(rho(np.arange(len(s)))), rho_lo=float(np.nanquantile(rb, .05)),
               rho_hi=float(np.nanquantile(rb, .95)), noise_share=float(nu.sum() / var_u.sum()),
               boot_nan=float(np.isnan(rb).mean()))
    for h, t in THRESH.items():
        row[f"p_{h}"] = float((1 + np.sum(~(rb > t))) / (1 + N_BOOT))  # nan (non-positive denominator) counts as fail
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    global N_BOOT
    limit_threads(1)
    t0 = time.time()
    out = OUT / ("score_a339_check" if a.check else "score_a339")
    out.mkdir(parents=True, exist_ok=True)
    units = B.load_units()
    if a.check:
        N_BOOT = 300
        rng = np.random.default_rng(0)
        for u in units:
            th = lambda x: pd.DataFrame(rng.binomial(x.to_numpy(np.int64), 0.15), columns=x.columns)  # noqa: E731
            u["u"] = th(u["s"])
            u["obs"]["depth_u"] = u["u"].sum(axis=1).to_numpy(float)
            if u["d"] is not None:
                u["d"] = th(u["s"])
                u["obs"]["depth_d"] = u["d"].sum(axis=1).to_numpy(float)
    print(f"load {time.time() - t0:.0f}s", flush=True)
    S["tasks"] = [(u, w, which) for u in units for w in B.WAVEFORMS
                  for which in (("u", "d") if u["d"] is not None else ("u",))]
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        rows = list(pool.imap_unordered(task, range(len(S["tasks"]))))
    t = pd.DataFrame(rows).sort_values(["role", "unit", "waveform", "unspliced"])
    t.to_csv(out / "rho.csv", index=False)
    v = []
    for (h, w) in sorted(SCORED):
        x = t[t.unit.isin(B.SCORED) & (t.waveform == w) & (t.unspliced == "intron")].copy()
        x["hypothesis"], x["p"] = h, x[f"p_{h}"]
        x["holm"] = dev.holm(x.p.to_numpy())
        x["pass"] = x.holm < .05
        v.append(x[["hypothesis", "unit", "waveform", "n", "rho", "rho_lo", "rho_hi", "p", "holm", "pass"]])
    v = pd.concat(v)
    v.to_csv(out / "verdicts.csv", index=False)
    pd.set_option("display.width", 220)
    print(t.round(3).to_string(index=False), v.round(4).to_string(index=False), f"runtime {time.time() - t0:.0f}s",
          sep="\n", flush=True)


if __name__ == "__main__":
    main()
