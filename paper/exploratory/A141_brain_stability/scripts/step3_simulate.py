"""A141 step 3 (decision): is the anchor dependence of the brain contrasts what noise alone produces?

Per brain contrast (12, the WP0 item 8 sets), `--nsim` simulations (registered 200) at the fitted amplitude, phase offset
and own noise psi of every compartment (fitted on the contrast's donors on the A17 host), the real donors, groups,
rows, measurement covariances and strata, host phases uniform on 24 h, no case effect. Each simulated data set is
refitted with the A17 anchors and with the all-LOO host; the statistic is the number of target compartments whose
case-minus-control d changes sign between the two hosts. Noise-driven: observed count <= the 95th percentile of the
simulated counts. The same simulations give the uniform-sampling reference of step 1 (Rayleigh R, spread and
eccentricity of the host fitted on the control donors and on all donors, A17 and all-compartment hosts).
Also written: the fitted parameters (`params_<contrast>.csv`), the observed d on both hosts, and the observed host
statistics. Choices not in the registration: see README.

    run.sh step3_simulate.py [--subsample N] [--nsim 200] [--jobs 6] [--chunk 10] [--only c1 c2 ...]
Output: output/step3/ (subsample: output/subsample/step3/)
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _a141 as H  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

SETS = ("control", "all")


def host_rows(d, C) -> dict:
    """Host statistics of one data set: A17 and all-compartment host, control donors and all donors."""
    lab = C.lab
    sets = {"control": lab.index[lab == "control"].to_numpy(), "all": lab.index.to_numpy()}
    out = {}
    for hn, anc in (("a17", C.anchors), ("allcomp", C.ds.compartments)):
        for sn, don in sets.items():
            for k, v in H.host_stats(d, anc, don).items():
                out[f"{hn}_{sn}_{k}"] = v
    return out


def task(name, chunk, sims, subsample, seed, outdir):
    limit_threads(1)
    f = outdir / f"sims_{name}_{chunk:03d}.csv"
    if f.exists():
        return str(f)
    C = H.Contrast(name, subsample, seed)
    par = H.fitted_params(C.ds)
    S = H.Simulator(C.ds, par)
    rows = []
    for s in sims:
        rng = np.random.default_rng([seed, H.CONTRASTS.index(name), s])
        d = S.draw(rng)
        da = H.d_values(H.with_data(C.ds, d), C.targets, C.blocks)
        dl = H.d_values(H.with_data(C.ds, d, C.ds.compartments), C.targets, C.blocks)
        rows.append(dict(contrast=name, sim=s, flips=H.n_flips(da, dl), n_targets=len(C.targets),
                         n_finite=int((np.isfinite(da) & np.isfinite(dl)).sum()), **host_rows(d, C)))
    pd.DataFrame(rows).to_csv(f, index=False)
    return str(f)


def observed_one(name, subsample, seed, outdir):
    """Observed d on both hosts (recomputed and compared to the WP0 table), observed host statistics, fitted params."""
    limit_threads(1)
    C = H.Contrast(name, subsample, seed)
    par = H.fitted_params(C.ds)
    par.to_csv(outdir / f"params_{name}.csv")
    da = H.d_values(C.ds, C.targets, C.blocks)
    dl = H.d_values(C.loo(), C.targets, C.blocks)
    o = C.obs.reindex(C.targets)
    t = pd.DataFrame(dict(contrast=name, target=C.targets, d_anchored=da.to_numpy(), d_loo=dl.to_numpy(),
                          d_anchored_wp0=o.contrast_anchored_recomputed.to_numpy(), d_loo_wp0=o.contrast_all_loo.to_numpy(),
                          se_anchored_wp0=o.se_anchored_recomputed.to_numpy(), se_loo_wp0=o.se_all_loo.to_numpy(),
                          n_control=o.n_control.to_numpy(), n_case=o.n_case.to_numpy()))
    t.to_csv(outdir / f"observed_d_{name}.csv", index=False)
    return dict(contrast=name, n_targets=len(C.targets), n_control_donors=C.n_groups[0], n_case_donors=C.n_groups[1],
                flips_recomputed=H.n_flips(da, dl), flips_wp0=int((o.sign_change == True).sum()),  # noqa: E712
                max_abs_d_diff_wp0=float(max(np.nanmax(np.abs(da - o.contrast_anchored_recomputed)),
                                             np.nanmax(np.abs(dl - o.contrast_all_loo)))),
                **host_rows(C.ds.data, C))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--nsim", type=int, default=200)
    ap.add_argument("--chunk", type=int, default=10)
    ap.add_argument("--only", nargs="+", default=H.CONTRASTS)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    outdir = H.OUT / ("subsample" if run["quick"] else "") / "step3"
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print(f"observed ({len(args.only)} contrasts)", flush=True)
    obs = Parallel(n_jobs=args.jobs)(delayed(observed_one)(n, args.subsample, args.seed, outdir) for n in args.only)
    pd.DataFrame(obs).to_csv(outdir / "observed.csv", index=False)
    print(f"observed done {time.time() - t0:.0f}s", flush=True)
    jobs = [(n, c, list(range(c * args.chunk, min((c + 1) * args.chunk, args.nsim))))
            for n in args.only for c in range(-(-args.nsim // args.chunk))]
    done = 0
    gen = Parallel(n_jobs=args.jobs, return_as="generator_unordered")(
        delayed(task)(n, c, s, args.subsample, args.seed, outdir) for n, c, s in jobs)
    for f in gen:
        done += 1
        print(f"[{time.strftime('%H:%M:%S')}] {done}/{len(jobs)} chunks {Path(f).name} ({time.time() - t0:.0f}s)", flush=True)
    sims = pd.concat([pd.read_csv(outdir / f"sims_{n}_{c:03d}.csv") for n, c, _ in jobs], ignore_index=True)
    sims.to_csv(outdir / "sims.csv", index=False)
    o = pd.DataFrame(obs).set_index("contrast")
    rows = []
    for n in args.only:
        x = sims[sims.contrast == n].flips.to_numpy()
        p95 = float(np.quantile(x, 0.95, method="inverted_cdf"))
        ob = int(o.loc[n, "flips_wp0"])
        rows.append(dict(contrast=n, n_targets=int(o.loc[n, "n_targets"]), observed_flips=ob,
                         observed_flips_recomputed=int(o.loc[n, "flips_recomputed"]), n_sim=len(x), sim_mean=x.mean(),
                         sim_median=float(np.median(x)), sim_p95=p95, sim_max=int(x.max()), p_upper=float((x >= ob).mean()),
                         noise_driven=bool(ob <= p95)))
    s = pd.DataFrame(rows)
    s.to_csv(outdir / "step3_summary.csv", index=False)
    print(s.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
