"""A19 part 2, validations 1 and 2: null calibration and planted effects of the age model at real designs (acceptance simulations).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A19 part 2, "Validation before the real data", and the A19 note of 2026-10-13: Wald statistics,
primary p from a null simulated at each compartment's real design, estimates centred by that null's median). At the real designs of
`psychad_mssm` and `lung_hlca` (healthy donors, ages, strata, learned anchors, as `resilience/compute_age.py`), for the compartments
with the highest, median and lowest anchored coupling, `Dataset.covariate` is run on:
  1. null (w = 0, gamma = 0): the donors' ages are permuted within stratum (the real clock data stay, so the design, host and noise are real and
     there is no association with the covariate); rejection rate of every estimand at alpha 0.05 (`_p` from the package's null, `--n-null`;
     `--null permutation` (default since clocksync 613e224, PI Q10 2026-10-26) or `simulated` (the archive A19 note; fails here));
     acceptance <= 0.07.
  2. planted Re(w/z) = -0.1 per decade and gamma = +0.15 per decade: on the permuted ages, the compartment's clock plane is changed by
     binomial thinning of its counts (`clocksync.simulate.thin_clock`): host-locked part + w x z f_t with w/z = -0.1, own noise raised to
     psi exp(gamma (x - x_min)) by adding N(0, psi (exp(gamma (x - x_min)) - 1)) in the plane; z, psi from the known-host regression on the
     fit's host (local). Acceptance: median estimate (centred, `_c`) within 25 % of the truth; power = rejection rate.
Estimands: locked_amp, drift_h, gamma, coupling_slope (`Dataset.covariate`).

    python compute_age_acceptance.py [--subsample N] [--n-sim 200] [--n-null 200] [--null permutation|simulated] [--tag T] [--jobs 4]
`--subsample N`: N donors per design (`min_age` N / 2), 3 simulations, 20 null draws.
Outputs: paper/results/extended_data/output/calibration/age_acceptance_{runs,summary}.csv
Legacy: legacy/desync/20_age/scripts/acceptance.py (outputs: desync/20_age/output/acceptance_summary.csv).
Deviation: the legacy simulated every data set from the fitted state-level model (beta = mu + M(z) f + M(w) x f + e, host posterior noise); here the real
count tables are modified (permuted ages; thinning), a count-level simulation at the same design. Legacy numbers are not reproduced draw by draw.
"""

from __future__ import annotations

import sys
import warnings
import zlib
from pathlib import Path

PAPER = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("resilience")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import thin_clock  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import compute_age as CA  # noqa: E402  (resilience)

import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location("s12_shared", group_dir("species") / "_shared.py")
S12 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S12)  # local helper of figure_S12 (known_host); loaded under its own name, `_shared` is figure_5's

FOLDER = "calibration"
DESIGNS = ("psychad_mssm", "lung_hlca")
ESTIMANDS = CA.ESTIMANDS
PLANT = dict(locked_amp=-0.1, gamma=0.15)
GENES = list(cs.gtex().genes)


def design(name, subsample, seed):
    tg = CA.sh.targets("healthy")
    pb, donors = CA.healthy_fit(next(tg[tg.name == name].itertuples()), subsample, seed)
    age = donors["age"].astype(float).dropna()
    age = age[age.index.isin(set(pb.donor))]
    rule = "learned" if CA.sh.anchor_set(name, "learned") else "fixed"
    anc = [a for a in CA.sh.anchor_set(name, rule) if a in set(pb.compartment)]
    return pb, donors, age, anc


def picks(ds, age, min_age):
    """Compartments with the highest, median and lowest anchored coupling among those with >= `min_age` aged donors."""
    n = ds.table[ds.table.resection.isin(age.index)].groupby("compartment").resection.nunique()
    c = ds.fit(n_boot=0, gate=False).table.coupling
    c = c[c.index.isin(n.index[n >= min_age])].dropna().sort_values()
    return [c.index[-1], c.index[len(c) // 2], c.index[0]], c


def permuted_ages(ds, age, rng):
    """Ages shuffled among the donors within stratum (donors without age keep none)."""
    st = ds.table.groupby("resection").stratum.first().reindex(age.index)
    out = age.copy()
    for s in st.unique():
        i = age.index[(st == s).to_numpy()]
        out.loc[i] = rng.permutation(age.loc[i].to_numpy())
    return out


def plant(ds, host, k, x, rng, w_ratio, gamma):
    """Count table with the age effects planted in compartment k (x: centred age in decades, indexed by donor)."""
    d = ds.data
    i = np.flatnonzero((d.compartment == k) & np.isin(d.resection, x.index) & np.isin(d.resection, host.index))
    f = host.loc[d.resection[i], ["x", "y"]].to_numpy(float)
    kh = S12.known_host(d.beta[i], d.V[i], f, np.asarray(d.group)[i])
    xi = x.reindex(d.resection[i]).to_numpy(float)
    z = kh["z"]
    locked = w_ratio * xi[:, None] * np.column_stack([z.real * f[:, 0] - z.imag * f[:, 1], z.imag * f[:, 0] + z.real * f[:, 1]])
    sd = np.sqrt(kh["psi"] * (np.exp(gamma * (xi - xi.min())) - 1.0))
    delta = locked + rng.standard_normal((len(i), 2)) * sd[:, None]
    t = ds.table.copy()
    t.iloc[i, t.columns.get_indexer(GENES)] = thin_clock(t.iloc[i][GENES].to_numpy(float), ds.waveform.X, delta, rng)
    return t


def one(name, k, scenario, sim, n_null, subsample, seed, min_age, null="permutation"):
    limit_threads(1)
    warnings.simplefilter("ignore")
    pb, donors, age, anc = design(name, subsample, seed)
    ds = cs.Dataset(pb, anchors=anc, donors=donors, name=name)
    rng = np.random.default_rng([seed, sim, scenario == "planted", zlib.crc32(str(k).encode()) % 10**6])  # crc32: stable across processes (str hash is salted)
    a = permuted_ages(ds, age, rng) / 10.0
    x = a - a.median()
    if scenario == "planted":
        host = ds.fit(n_boot=0, gate=False).host
        ds = cs.Dataset(plant(ds, host, k, x, rng, PLANT["locked_amp"], PLANT["gamma"]), anchors=anc, donors=donors, name=name)
    r = ds.covariate(x, targets=[k], center=False, n_null=n_null, seed=sim, null=null)
    if "error" in r:
        return None
    r = r.apply(pd.to_numeric, errors="coerce").iloc[0]
    row = dict(design=name, compartment=k, scenario=scenario, sim=sim, coupling=r.get("coupling"))
    for e in ESTIMANDS:
        row[e], row[f"{e}_c"], row[f"{e}_p"] = r.get(e), r.get(f"{e}_c"), r.get(f"{e}_p")
    return row


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-sim", type=int, default=200)
    ap.add_argument("--n-null", type=int, default=None)
    ap.add_argument("--null", choices=("permutation", "simulated"), default="permutation")
    ap.add_argument("--tag", default="", help="suffix of the output files (e.g. _simulated)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_sim = 3 if run["quick"] else args.n_sim
    n_null = args.n_null or (20 if run["quick"] else 200)
    min_age = args.subsample // 2 if run["quick"] else CA.MIN_AGE_DONORS
    jobs = []
    for name in DESIGNS:
        pb, donors, age, anc = design(name, args.subsample, args.seed)
        ks, c = picks(cs.Dataset(pb, anchors=anc, donors=donors), age, min_age)
        print(name, "donors with age", len(age), "picked", {k: round(float(c[k]), 3) for k in ks}, flush=True)
        jobs += [(name, k, sc, s, n_null, args.subsample, args.seed, min_age, args.null) for k in ks for sc in ("null", "planted") for s in range(n_sim)]
    rows = [r for r in Parallel(n_jobs=args.jobs)(delayed(one)(*j) for j in jobs) if r]
    df = pd.DataFrame(rows)
    df.to_csv(out / f"age_acceptance_runs{args.tag}.csv", index=False)
    summ = []
    for (n, k, sc), g in df.groupby(["design", "compartment", "scenario"]):
        r = dict(design=n, compartment=k, scenario=sc, n=len(g), coupling=g.coupling.median())
        for e in ESTIMANDS:
            r |= {f"{e}_median": g[e].median(), f"{e}_c_median": g[f"{e}_c"].median(), f"{e}_reject": float((g[f"{e}_p"] < 0.05).mean())}
        summ.append(r)
    s = pd.DataFrame(summ)
    s.to_csv(out / f"age_acceptance_summary{args.tag}.csv", index=False)
    pd.set_option("display.width", 250)
    print(s.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
