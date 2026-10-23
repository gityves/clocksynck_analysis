"""Sampling-window checks of the method (A18 part 1): does the shared clock variance between donors scale with
1 - R_W^2, R_W = |mean exp(i omega t)| over the collection times?

Pre-registered as A18 part 1 (2026-10-12; docs/PLAN_internal_desync.md (archive)). Statistic for two donor sets A and B of one
table: per compartment pair (k, l), the complex cross-covariance gamma_kl of the stratum-centred clock coefficients
(b = beta_1 + i beta_2, the clock-plane states of `Dataset`) over the donors the pair shares in the set;
ratio = sum n_kl Re(gamma_kl^A u_kl*) / sum n_kl Re(gamma_kl^B u_kl*), u_kl the unit phase of gamma_kl on all donors
(a signed projection: no positive bias from noise). Donor bootstrap within sets. Pass: the 95 % interval of
observed / predicted contains 1.
  1a GTEx bulk, fast deaths, recorded time of death: died 08-18 h (A) vs all others (B); predicted ratio
     (1 - R_A^2) / (1 - R_B^2) from the recorded times. (Also 08-12, 12-16, 20-24 h windows: post hoc additions of the run.)
  1b HLCA healthy lung (A4 table): hospital (surgical resection, biopsy, brush) vs any time (donor lung, explant,
     autopsy); predicted 0.455 for 08-18 h vs uniform. Descriptive: the sampling class is nested in study and platform.
Controlled access: the times of death stay in memory; the output has the window sums and ratios only.

    python compute_window_check.py [--subsample N] [--seed S] [--n-boot N]
Outputs (paper/results/figure_1/output/windows/): window_check.csv
Legacy: legacy/desync/common/window_check.py -> common_output/a18_window_check.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "windows"
HOSPITAL = {"surgical_resection", "biopsy", "brush", "biopsy;brush"}
ANYTIME = {"donor_lung", "donor_lung_cryopreserved", "lung_explant", "autopsy"}


def states(ds: cs.Dataset):
    """Complex clock coefficients per row, donors, compartments, strata (the states `Dataset` reads)."""
    d = ds.data
    return d.beta[:, 0] + 1j * d.beta[:, 1], d.resection, d.compartment, d.stratum


def _centred(b, don, comp, strat, keep):
    """Complex coefficients per (donor, compartment), centred per compartment x stratum over `keep` donors."""
    d = pd.DataFrame(dict(re=b.real, im=b.imag, don=don, comp=comp, strat=strat))
    d = d[d.don.isin(keep)]
    g = d.groupby(["comp", "strat"])
    d["re"], d["im"] = d.re - g.re.transform("mean"), d.im - g.im.transform("mean")
    d = d[g.re.transform("size") >= 2]
    w = lambda v: d.pivot_table(index="don", columns="comp", values=v, aggfunc="first")
    return w("re") + 1j * w("im")


def gammas(Wd: pd.DataFrame, min_shared: int = 5) -> dict:
    out, cols = {}, list(Wd.columns)
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            x = Wd[[cols[i], cols[j]]].dropna()
            if len(x) >= min_shared:
                out[(cols[i], cols[j])] = (len(x), np.mean(x.iloc[:, 0].to_numpy() * np.conj(x.iloc[:, 1].to_numpy())))
    return out


def ratio(b, don, comp, strat, set_a, set_b, u) -> float:
    ga, gb = gammas(_centred(b, don, comp, strat, set_a)), gammas(_centred(b, don, comp, strat, set_b))
    common = [p for p in u if p in ga and p in gb]
    na = sum(ga[p][0] * np.real(ga[p][1] * np.conj(u[p])) for p in common)
    nb = sum(gb[p][0] * np.real(gb[p][1] * np.conj(u[p])) for p in common)
    return na / nb if nb > 0 else np.nan


def run(label, st, set_a, set_b, predicted, n_boot, rng) -> dict:
    b, don, comp, strat = st
    allg = gammas(_centred(b, don, comp, strat, np.r_[set_a, set_b]))
    u = {p: g / abs(g) for p, (n, g) in allg.items() if abs(g) > 0}
    r = ratio(b, don, comp, strat, set_a, set_b, u)
    reps = []
    for _ in range(n_boot):
        ra, rb = rng.choice(set_a, len(set_a)), rng.choice(set_b, len(set_b))
        idx, dd = [], []
        for j, x in enumerate(np.r_[ra, rb]):  # duplicated donors: relabel copies so that each draw is its own row
            rows = np.flatnonzero(don == x)
            idx.append(rows), dd.append(np.full(len(rows), f"{x}#{j}"))
        idx, dd = np.concatenate(idx), np.concatenate(dd)
        na = np.array([f"{x}#{j}" for j, x in enumerate(ra)])
        nb = np.array([f"{x}#{j + len(ra)}" for j, x in enumerate(rb)])
        reps.append(ratio(b[idx], dd, comp[idx], strat[idx], na, nb, u))
    reps = np.array(reps, float)
    reps = reps[np.isfinite(reps)]
    lo, hi = np.quantile(reps, [0.025, 0.975]) if len(reps) > 1 else (np.nan, np.nan)
    return dict(check=label, donors_a=len(set_a), donors_b=len(set_b), pairs=len(u), ratio=r, lo=lo, hi=hi,
                predicted=predicted, obs_over_pred=r / predicted, obs_over_pred_lo=lo / predicted,
                obs_over_pred_hi=hi / predicted, interval_contains_prediction=bool(lo <= predicted <= hi),
                n_boot=len(reps))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None, help="bootstrap draws (default 500 as pre-registered; 20 with --subsample)")
    args = ap.parse_args(argv)
    run_ = cli.setup(args)
    out = out_dir(FOLDER, run_["quick"])
    nb = args.n_boot if args.n_boot is not None else (20 if run_["quick"] else 500)
    rng = np.random.default_rng(args.seed)
    rows = []
    pb, donors = sh.gtex_fast(30, args.subsample, args.seed)
    ds = cs.Dataset(pb, anchors=None, donors=donors, name="gtex_fast")
    st = states(ds)
    tod = donors.tod.dropna()
    for name, lo, hi in (("08-18 h", 8, 18), ("08-12 h", 8, 12), ("12-16 h", 12, 16), ("20-24 h", 20, 24)):
        a = tod.index[(tod >= lo) & (tod < hi)].to_numpy()
        rest = tod.index[~tod.index.isin(a)].to_numpy()
        if len(a) < 10 or len(rest) < 10:
            continue
        pred = sh.window_factor(tod[a]) / sh.window_factor(tod[rest])
        rows.append(run(f"GTEx fast deaths: died {name} vs others", st, a, rest, pred, nb, rng))
        print(rows[-1], flush=True)
    pbl, dl = data.load("lung_hlca", subsample=args.subsample, seed=args.seed)
    dl = dl[dl.healthy.astype(bool)]
    pbl = pbl[pbl.donor.isin(dl.index)].reset_index(drop=True)
    dsl = cs.Dataset(pbl, anchors=None, donors=dl, name="lung_hlca")
    meth = dl.tissue_sampling_method_major.astype(str)
    h, a = meth.index[meth.isin(HOSPITAL)].to_numpy(), meth.index[meth.isin(ANYTIME)].to_numpy()
    if len(h) >= 5 and len(a) >= 5:
        rows.append(run("HLCA healthy lung: hospital (resection, biopsy, brush) vs any time (donor, explant, autopsy)",
                        states(dsl), h, a, 0.455, nb, rng))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(out / "window_check.csv", index=False)


if __name__ == "__main__":
    main()
