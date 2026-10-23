"""A17 acceptance simulations of the anchored host (V12): identity, coverage, boundary share, power.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17 step 5, written before any anchored coupling on real data). Cohorts from
`clocksync.simulate.simulate` with the global GTEx waveform, hospital-hours sampling (08-18 h plus 1 h jitter), depth 4e6. Anchors:
three compartments coupled to host time (loadings 1.0, 0.9, 1.1).
1. identity: a competing module of three compartments with its own clock state, independent of host time, same size as the anchors;
   100 cohorts of 42 donors. Share of cohorts in which the free host (all compartments, reference in the anchor or in the competing
   module) is closer to the competing state than to host time; share with anchored host similarity to host time >= 0.9; anchored
   stability (share of 100 donor-bootstrap hosts with similarity > 0.7 to the full host) over the first 20 cohorts. Criteria 0.9, 0.9.
2. coverage: 30 donors, target with true coupling c in {0.3, 0.6, 0.9}, anchors precise (psi 0.002) or noisy (0.05); 200 cohorts
   each, `Dataset.fit(targets=["k"])` with a 200-draw donor bootstrap. Criteria: coverage >= 0.9, |median bias| <= 0.05.
3. boundary: 20 donors, c in {0.95, 1.0}: share of estimates at the boundary (ML psi at its floor). The test-inversion lower bound
   (`anchored.lower_bound`) is not in the public API: not ported.
4. power: 30 + 30 donors, target coupling 0.8 in both groups, damping s in {1, 0.6, 0} (`clocksync.power.mde`, 100 replicates, 200
   permutations). Criteria: power rises with damping; rejection rate at s = 1 <= 0.07.
The per-fit bias prediction (V12 part 5, `anchored.design_bias`) is not ported (API gap).

    python compute_anchored_simulations.py [--subsample N] [--only identity coverage boundary power] [--jobs 4]
`--subsample N`: check mode (N cohorts per part, 10 bootstraps, 20 permutations).
Outputs: paper/results/extended_data/output/calibration/anchored_sim_{identity,coverage,boundary,power,summary}.csv
Legacy: legacy/scripts/v12_anchored.py (outputs/v12_anchored/)
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import power  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "calibration"
WF = cs.gtex()
HOURS = 24 / (2 * np.pi)
ANCH = ["a0", "a1", "a2"]
LOAD = [1.0, 0.9, 1.1]


def f_of(t):
    return np.column_stack([np.cos(np.asarray(t) / HOURS), np.sin(np.asarray(t) / HOURS)])


def similarity(fa, fb) -> float:
    """local: not in clocksync. |complex correlation| of two host states over the same donors (invariant to rotation and scale)."""
    a, b = np.asarray(fa, float), np.asarray(fb, float)
    za, zb = (a[:, 0] - a[:, 0].mean()) + 1j * (a[:, 1] - a[:, 1].mean()), (b[:, 0] - b[:, 0].mean()) + 1j * (b[:, 1] - b[:, 1].mean())
    den = np.sqrt(np.sum(np.abs(za) ** 2) * np.sum(np.abs(zb) ** 2))
    return float(np.abs(np.sum(za * np.conj(zb))) / den) if den > 0 else np.nan


def cohort(n, seed, c=None, apsi=0.002, module=False, zk=0.8, prefix="r"):
    """Anchors (+ target 'k' with true coupling c to host time) (+ competing module b0-b2). Returns table, t."""
    df, t = cs.simulate(WF, {f"a{i}": cs.Compartment(z, psi=apsi, tau=0.0) for i, z in enumerate(LOAD)}, n=n, seed=seed, tau=0.0)
    parts = [df]
    if c is not None:
        Vt = np.trace(np.cov(f_of(t).T))
        psi = zk ** 2 * Vt * (1 - c) / (2 * c)
        parts.append(cs.simulate(WF, {"k": cs.Compartment(zk, psi=psi, tau=0.0)}, n=n, seed=seed, tau=0.0)[0])
    if module:
        parts.append(cs.simulate(WF, {f"b{i}": cs.Compartment(z, psi=apsi, tau=0.0) for i, z in enumerate(LOAD)},
                                 n=n, seed=seed + 10 ** 6, tau=0.0)[0])
    out = pd.concat(parts, ignore_index=True)
    out["resection"] = prefix + out.resection
    return out, t


def host_of(ds: cs.Dataset, seed=0) -> pd.DataFrame:
    return ds.fit(targets=[ds.anchors[0]], n_boot=0, gate=False).host[["x", "y"]]


def dataset(df, anchors):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cs.Dataset(df, anchors=anchors)


# ---- 1. identity -------------------------------------------------------------------------------------------------------

def identity_one(j, n_stab):
    cs.config.limit_threads(1)
    df, t = cohort(42, 100 + j, module=True)
    tb = cs.simulate(WF, {"x": cs.Compartment()}, n=42, seed=100 + j + 10 ** 6)[1]
    comps = sorted(df.compartment.unique())
    out = dict(cohort=j)
    for name, anc in (("free_ref_anchor", ["a0"] + [c for c in comps if c != "a0"]),
                      ("free_ref_competing", ["b0"] + [c for c in comps if c != "b0"]), ("anchored", ANCH)):
        h = host_of(dataset(df, anc))
        idx = [int(r[2:]) for r in h.index]
        out[f"{name}_host_time"] = similarity(h.to_numpy(), f_of(t)[idx])
        out[f"{name}_competing"] = similarity(h.to_numpy(), f_of(tb)[idx])
    if n_stab:
        ds = dataset(df, ANCH)
        base = host_of(ds)
        rng = np.random.default_rng(j)
        sims = []
        for _ in range(n_stab):
            hb = host_of(ds.resample(rng))
            orig = np.array([r.split("#")[0] for r in hb.index])
            _, first = np.unique(orig, return_index=True)
            sims.append(similarity(hb.to_numpy()[first], base.loc[orig[first]].to_numpy()))
        out["anchored_stability"] = float(np.mean(np.asarray(sims) > 0.7))
    return out


def identity(jobs, n, n_stab):
    r = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(identity_one)(j, n_stab if j < min(20, n) else 0) for j in range(n)))
    rows = [dict(check=f"identity: free host (reference in {m} module) closer to the competing module (share of cohorts)",
                 value=float((r[f"free_ref_{m}_competing"] > r[f"free_ref_{m}_host_time"]).mean()), criterion="descriptive")
            for m in ("anchor", "competing")]
    rows += [dict(check="identity: anchored host similarity to host time >= 0.9 (share of cohorts)",
                  value=float((r.anchored_host_time >= 0.9).mean()), criterion=">= 0.9"),
             dict(check="identity: anchored stability (median over cohorts)", value=float(r.anchored_stability.median()), criterion=">= 0.9")]
    return r, rows


# ---- 2. coverage ------------------------------------------------------------------------------------------------------

def coverage_one(c, apsi, rep, n_boot):
    cs.config.limit_threads(1)
    df, _ = cohort(30, 1000 + rep, c=c, apsi=apsi)
    f = dataset(df, ANCH).fit(targets=["k"], n_boot=n_boot, seed=rep, gate=False).table.loc["k"]
    return dict(c=c, apsi=apsi, rep=rep, est=f.coupling, lo=f.coupling_lo, hi=f.coupling_hi)


def coverage(jobs, n, n_boot):
    r = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(coverage_one)(c, a, rep, n_boot) for c in (0.3, 0.6, 0.9) for a in (0.002, 0.05)
                                           for rep in range(n)))
    rows = []
    for (c, a), g in r.groupby(["c", "apsi"]):
        rows += [dict(check=f"coverage c={c} anchor psi={a}", value=float(((g.lo <= c) & (g.hi >= c)).mean()), criterion=">= 0.9"),
                 dict(check=f"bias c={c} anchor psi={a}", value=float(g.est.median() - c), criterion="|.| <= 0.05")]
    return r, rows


# ---- 3. boundary ------------------------------------------------------------------------------------------------------

def boundary_one(c, rep):
    cs.config.limit_threads(1)
    df, _ = cohort(20, 5000 + rep, c=c)
    f = dataset(df, ANCH).fit(targets=["k"], n_boot=0, gate=False).table.loc["k"]
    return dict(c=c, rep=rep, boundary=bool(f.boundary), est=f.coupling)


def boundary(jobs, n):
    r = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(boundary_one)(c, rep) for c in (0.95, 1.0) for rep in range(n)))
    return r, [dict(check=f"boundary c={c}: share of estimates at the boundary", value=float(g.boundary.mean()), criterion="descriptive")
               for c, g in r.groupby("c")]


# ---- 4. power ---------------------------------------------------------------------------------------------------------

def power_part(jobs, n_rep, n_perm):
    a, _ = cohort(30, 77, c=0.8, prefix="ctl")
    b, _ = cohort(30, 78, c=0.8, prefix="cas")
    ds = dataset(pd.concat([a, b], ignore_index=True), ANCH)
    lab = pd.Series({r: ("control" if r.startswith("ctl") else "case") for r in ds.table.resection.unique()})
    m = power.mde(ds.study, lab, ds.anchors, "k", scales=(1.0, 0.6, 0.0), n_rep=n_rep, n_perm=n_perm, seed=79, jobs=jobs)
    curve = m["curve"]
    pw = curve.set_index("s").power
    return curve, [dict(check="power: rejection at s = 1", value=float(pw[1.0]), criterion="<= 0.07"),
                   dict(check="power: monotone in damping (s 1 <= 0.6 <= 0)", value=float(pw[1.0] <= pw[0.6] <= pw[0.0]), criterion="1"),
                   dict(check="power: MDE (coupling drop at 80 % power)", value=m["mde"], criterion="descriptive")]


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="*", default=["identity", "coverage", "boundary", "power"])
    args = ap.parse_args()
    opt = cli.setup(args)
    quick = opt["quick"]
    N = args.subsample
    n_id, n_stab = (N, 10) if quick else (100, 100)
    n_cov, n_boot = (N, 10) if quick else (200, 200)
    n_bound = N if quick else 200
    n_rep, n_perm = (N, 20) if quick else (100, 200)
    out = out_dir(FOLDER, quick)
    summary = []
    for p in args.only:
        r, rows = {"identity": lambda: identity(args.jobs, n_id, n_stab), "coverage": lambda: coverage(args.jobs, n_cov, n_boot),
                   "boundary": lambda: boundary(args.jobs, n_bound), "power": lambda: power_part(args.jobs, n_rep, n_perm)}[p]()
        r.to_csv(out / f"anchored_sim_{p}.csv", index=False)
        rows = pd.DataFrame(rows).assign(part=p)
        print(rows.to_string(index=False), flush=True)
        summary.append(rows)
    pd.concat(summary).to_csv(out / "anchored_sim_summary.csv", index=False)


if __name__ == "__main__":
    main()
