"""A234: does internal time drift from external time with age, and does the brain drift differently from the body?

Registered in paper/docs/registered_experiments.md (A234). GTEx fast deaths (Hardy 1-2) with time of death. Two untimed
hosts per donor (body: the 12 A81 / A92 body tissues; brain: 11 regions); Δ = host phase − time of death (h) or brain −
body host phase; circular-linear regression of Δ on exact age (h per decade, > 0 = advance), strata Hardy x sex.
Two stages:
  power  D234.0 direction check, wall-time alignment R, age-permutation null of β (no observed β is computed or printed):
         null SD and MDE_80 (2.80 x null SD) per test -> output/{power,direction,counts}.csv.
  score  (after the power note is dated) observed β, permutation p, donor bootstrap (hosts refitted), sensitivities
         (Hardy 2 only, bracket age, brain host without the cerebellar regions), 50 pseudo-clock sets for H234.2
         -> output/score_{tests,sensitivity,pseudo}.csv. `--placebo`: ages shuffled over donors (code test).

    python compute_a234.py --stage power|score [--subsample N] [--jobs 8]
Controlled access: exact age and time of death stay in memory; only aggregates are written.

Choices not fixed by the registration
- Hosts: `cs.Dataset(pb, anchors=<all tissues of the host>).host()` phase_h (default depth "sum", GTEx waveform); a donor
  enters a host with >= 5 body tissues / >= 3 brain regions; H234.2 uses donors with both.
- Estimator: β maximises Σ_s n_s R_s(β), R_s the resultant length of (Δ − β·age) in stratum s (μ_s profiled out); grid
  over [-3, 3] h per decade (step 0.05), then bounded refinement within ±0.05 of the grid maximum (β is bounded at
  ±3.05: a host without time-of-day signal, e.g. a pseudo-clock host, gives β at the bound). Age in decades, centred.
- Strata with < 3 donors are merged into the other sex of the same Hardy class (permutation within the merged stratum).
- D234.0: Fisher-Lee(host phase, time of death) > 0 required; a negative host is negated (phase -> -phase) before Δ.
- Bootstrap: donors resampled with replacement (copies renamed), both hosts refitted on the replicate, β recomputed.
- Pseudo-clock: `compute_a208.pseudo_sets` (expression-matched, seed from "A234"), both hosts refitted on the pseudo genes;
  the D234.0 orientation of each pseudo host is set the same way (Fisher-Lee sign with time of death).
"""

from __future__ import annotations

import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "exploratory" / "A208_sex_coupling" / "scripts"))
sys.path.insert(0, str(ROOT / "results" / "figure_3" / "scripts" / "gtex_cues"))
import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync import stats as CS  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402

OUT = HERE.parent / "output"
W = 2 * np.pi / 24
BODY = list(S.BODY)
BRAIN = ["Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate (basal ganglia)",
         "Brain - Cerebellar Hemisphere", "Brain - Cerebellum", "Brain - Cortex", "Brain - Frontal Cortex (BA9)",
         "Brain - Hippocampus", "Brain - Hypothalamus", "Brain - Nucleus accumbens (basal ganglia)",
         "Brain - Putamen (basal ganglia)"]
NO_CEREBELLUM = [t for t in BRAIN if "Cerebell" not in t]
MIN_T = {"body": 5, "brain": 3}
GRID = np.arange(-3, 3.0001, 0.05)
N_PERM, N_BOOT, N_PSEUDO = 10_000, 200, 50


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def seed_of(*parts) -> int:
    return zlib.crc32("|".join(map(str, parts)).encode()) % (2 ** 31)


# ---- data and hosts ---------------------------------------------------------------------------------------------------

def load(genes="clock"):
    """{host: pseudobulk filtered to donors with enough tissues}, donors (hardy, tod, male, bracket age, exact age)."""
    pb, donors = gtex.load(tissues=BODY + BRAIN, hardy=gtex.FAST, genes=genes)
    donors["age_exact"] = gtex.phenotypes(["AGE"]).AGE.reindex(donors.index).astype(float)
    out = {}
    for h, tis in (("body", BODY), ("brain", BRAIN)):
        p = pb[pb.compartment.isin(tis)]
        n = p.groupby("donor").compartment.nunique()
        out[h] = p[p.donor.isin(n.index[n >= MIN_T[h]])].reset_index(drop=True)
    return out, donors


def host_phase(pb, anchors) -> pd.Series:
    """Host phase (h) per donor, anchored on every tissue of the host present in pb."""
    ds = cs.Dataset(pb, anchors=[a for a in anchors if a in set(pb.compartment)])
    return ds.host().phase_h


def orient(ph: pd.Series, tod: pd.Series) -> tuple[pd.Series, float]:
    """D234.0: phases negated if their Fisher-Lee correlation with time of death is negative; returns (phases, FL)."""
    fl = CS.fisher_lee(W * ph.to_numpy(), W * tod.reindex(ph.index).to_numpy())
    return (ph if fl >= 0 else -ph), float(fl)


def deltas(pbs, donors, brain_anchors=BRAIN, hosts=None) -> tuple[dict, dict]:
    """{test: Δ series (h)} for H234.1 (body − tod), H234.2 (brain − body), H234.3 (brain − tod); D234.0 rows."""
    hosts = hosts or {"body": host_phase(pbs["body"], BODY),
                      "brain": host_phase(pbs["brain"][pbs["brain"].compartment.isin(brain_anchors)], brain_anchors)}
    ph, fl = {}, {}
    for h, v in hosts.items():
        ph[h], fl[h] = orient(v, donors.tod)
    both = ph["body"].index.intersection(ph["brain"].index)
    d = {"H234.1": ph["body"] - donors.tod.reindex(ph["body"].index),
         "H234.2": ph["brain"].reindex(both) - ph["body"].reindex(both),
         "H234.3": ph["brain"] - donors.tod.reindex(ph["brain"].index)}
    return {k: pd.Series(CS.wrap(v.to_numpy()), index=v.index) for k, v in d.items()}, fl


# ---- estimator --------------------------------------------------------------------------------------------------------

def strata(donors: pd.DataFrame, idx) -> np.ndarray:
    """Hardy x sex codes; strata with < 3 donors merged into the other sex of the same Hardy class."""
    d = donors.loc[idx]
    s = np.array(d.hardy.astype(int) * 2 + d.male.astype(int), dtype=int)  # writable copy
    for v, n in zip(*np.unique(s, return_counts=True)):
        if n < 3:
            s[s == v] = v ^ 1
    return s


def _objective(beta, z, age, codes, ncode):
    """Σ_s |Σ_{i in s} exp(i(Δ_i − β a_i))| for a vector of β (rows) and one or many age vectors."""
    e = z * np.exp(-1j * W * np.multiply.outer(beta, age))          # (..., n)
    tot = np.zeros(e.shape[:-1] + (ncode,), complex)
    for k in range(ncode):
        tot[..., k] = e[..., codes == k].sum(-1)
    return np.abs(tot).sum(-1)


def slope(delta_h, age, codes) -> float:
    """β (h per decade) of the circular-linear regression with stratum intercepts (age in decades, centred)."""
    _, codes = np.unique(codes, return_inverse=True)
    ncode, z = codes.max() + 1, np.exp(1j * W * np.asarray(delta_h, float))
    a = (np.asarray(age, float) - np.mean(age)) / 10
    g = GRID[int(np.argmax(_objective(GRID, z, a, codes, ncode)))]
    r = minimize_scalar(lambda b: -_objective(np.array([b]), z, a, codes, ncode)[0], bounds=(g - .05, g + .05),
                        method="bounded", options={"xatol": 1e-4})
    return float(r.x)


def perm_null(delta_h, age, codes, n, seed) -> np.ndarray:
    """β under age permuted within strata."""
    rng = np.random.default_rng(seed)
    idx = [np.flatnonzero(codes == c) for c in np.unique(codes)]
    out = np.empty(n)
    for k in range(n):
        a = np.array(age, float)
        for i in idx:
            a[i] = a[rng.permutation(i)]
        out[k] = slope(delta_h, a, codes)
    return out


def _null_chunk(delta_h, age, codes, n, seed):
    limit_threads(1)
    return perm_null(delta_h, age, codes, n, seed)


def null_of(delta_h, age, codes, n, jobs, seed) -> np.ndarray:
    k = max(jobs, 1)
    parts = Parallel(n_jobs=jobs)(delayed(_null_chunk)(delta_h, age, codes, len(c), seed + j)
                                  for j, c in enumerate(np.array_split(np.arange(n), k)))
    return np.concatenate(parts)


def frame(d: pd.Series, donors, age_col="age_exact", keep=None):
    idx = d.index if keep is None else d.index[keep(donors.loc[d.index]).to_numpy()]
    return d.loc[idx].to_numpy(), donors.loc[idx, age_col].to_numpy(float), strata(donors, idx)


# ---- bootstrap and pseudo-clock ---------------------------------------------------------------------------------------

def _boot_chunk(pbs, donors, seeds):
    limit_threads(1)
    ids = np.array(sorted(set(pbs["body"].donor) | set(pbs["brain"].donor)))
    out = []
    for sd in seeds:
        draw = ids[np.random.default_rng(sd).choice(len(ids), len(ids))]
        m = pd.DataFrame({"orig": draw, "donor_new": [f"{d}@{j}" for j, d in enumerate(draw)]})
        rep = {h: p.merge(m, left_on="donor", right_on="orig").drop(columns=["donor", "orig"])
               .rename(columns={"donor_new": "donor"}) for h, p in pbs.items()}
        dn = donors.loc[m.orig].set_axis(m.donor_new)
        try:
            dl, _ = deltas(rep, dn)
            out.append([slope(*frame(dl[t], dn)) for t in ("H234.1", "H234.2", "H234.3")])
        except Exception as e:  # noqa: BLE001  (a failed host refit is recorded, not fatal)
            say(f"boot seed {sd} failed: {e}")
            out.append([np.nan] * 3)
    return out


def boots(pbs, donors, n, jobs, seed) -> np.ndarray:
    seeds = [seed_of("A234boot", seed, i) for i in range(n)]
    parts = Parallel(n_jobs=jobs)(delayed(_boot_chunk)(pbs, donors, list(c)) for c in np.array_split(seeds, max(jobs, 1)))
    return np.array([r for p in parts for r in p], float)


def _pseudo_one(w, donors_index):
    limit_threads(1)
    pbs, donors = load(genes=list(w.genes))
    pbs = {h: p[p.donor.isin(donors_index)] for h, p in pbs.items()}
    hosts = {}
    for h, tis in (("body", BODY), ("brain", BRAIN)):
        ds = cs.Dataset(pbs[h], anchors=[a for a in tis if a in set(pbs[h].compartment)], waveform=w)
        hosts[h] = ds.host().phase_h
    dl, fl = deltas(pbs, donors, hosts=hosts)
    return [slope(*frame(dl[t], donors)) for t in ("H234.1", "H234.2", "H234.3")] + [fl["body"], fl["brain"]]


def pseudo(donors_index, n, jobs, seed) -> pd.DataFrame:
    import compute_a208 as B
    ws = B.pseudo_sets(gtex.FAST, n, seed_of("A234pseudo", seed))
    rows = Parallel(n_jobs=jobs)(delayed(_pseudo_one)(w, donors_index) for w in ws)
    return pd.DataFrame(rows, columns=["beta_H234.1", "beta_H234.2", "beta_H234.3", "fl_body", "fl_brain"])


# ---- stages -----------------------------------------------------------------------------------------------------------

def subsample(pbs, donors, n, seed):
    keep = np.random.default_rng(seed).permutation(sorted(donors.index))[:n]
    return {h: p[p.donor.isin(keep)].reset_index(drop=True) for h, p in pbs.items()}, donors.loc[keep]


def power(pbs, donors, n_perm, jobs, seed):
    dl, fl = deltas(pbs, donors)
    pd.DataFrame([dict(host=h, fisher_lee_tod=fl[h], negated=fl[h] < 0, n=len(dl["H234.1" if h == "body" else "H234.3"]),
                       wall_R=float(abs(np.exp(1j * W * dl["H234.1" if h == "body" else "H234.3"]).mean())))
                  for h in fl]).to_csv(OUT / "direction.csv", index=False)
    rows, cnt = [], []
    for t, d in dl.items():
        x, a, s = frame(d, donors)
        null = null_of(x, a, s, n_perm, jobs, seed_of("A234null", t, seed))
        rows.append(dict(test=t, n=len(x), age_min=a.min(), age_max=a.max(), null_sd=null.std(),
                         null_q025=np.quantile(null, .025), null_q975=np.quantile(null, .975), mde_80=2.80 * null.std()))
        dd = donors.loc[d.index]
        cnt.append(dict(test=t, n=len(x), hardy1=int((dd.hardy == 1).sum()), hardy2=int((dd.hardy == 2).sum()),
                        male=int(dd.male.sum()), under40=int((dd.age_exact < 40).sum())))
        say(t, "null SD", round(null.std(), 3))
    pd.DataFrame(rows).to_csv(OUT / "power.csv", index=False)
    pd.DataFrame(cnt).to_csv(OUT / "counts.csv", index=False)


def score(pbs, donors, n_perm, n_boot, n_pseudo, jobs, seed):
    dl, _ = deltas(pbs, donors)
    bt = boots(pbs, donors, n_boot, jobs, seed)
    rows = []
    for j, (t, d) in enumerate(dl.items()):
        x, a, s = frame(d, donors)
        b = slope(x, a, s)
        null = null_of(x, a, s, n_perm, jobs, seed_of("A234null", t, seed))
        se = max(np.nanstd(bt[:, j]), null.std())
        rows.append(dict(test=t, n=len(x), beta_h_per_decade=b, p_two_sided=(1 + np.sum(np.abs(null) >= abs(b))) / (1 + n_perm),
                         null_sd=null.std(), boot_se=np.nanstd(bt[:, j]), n_boot_ok=int(np.isfinite(bt[:, j]).sum()),
                         lo90=b - 1.645 * se, hi90=b + 1.645 * se, mde_80=2.80 * null.std()))
        say(t, "scored")
    tab = pd.DataFrame(rows)
    conf = tab.test.isin(["H234.1", "H234.3"])
    p = tab.loc[conf, "p_two_sided"].to_numpy()
    o = np.argsort(p)
    holm = np.maximum.accumulate(np.minimum(1, p[o] * (len(p) - np.arange(len(p)))))
    tab.loc[conf, "p_holm"] = holm[np.argsort(o)]

    sens = []
    for t, label, kw in (("H234.1", "hardy2_only", dict(keep=lambda q: q.hardy == 2)),
                         ("H234.3", "hardy2_only", dict(keep=lambda q: q.hardy == 2)),
                         ("H234.1", "bracket_age", dict(age_col="age"))):
        x, a, s = frame(dl[t], donors, **kw)
        null = perm_null(x, a, s, 2000, seed_of("A234sens", t, label))
        b = slope(x, a, s)
        sens.append(dict(test=t, variant=label, n=len(x), beta_h_per_decade=b, null_sd=null.std(),
                         p_two_sided=(1 + np.sum(np.abs(null) >= abs(b))) / 2001))
    dl_nc, _ = deltas(pbs, donors, brain_anchors=NO_CEREBELLUM)
    for t in ("H234.2", "H234.3"):
        x, a, s = frame(dl_nc[t], donors)
        null = perm_null(x, a, s, 2000, seed_of("A234sens", t, "nocb"))
        b = slope(x, a, s)
        sens.append(dict(test=t, variant="brain_host_without_cerebellum", n=len(x), beta_h_per_decade=b,
                         null_sd=null.std(), p_two_sided=(1 + np.sum(np.abs(null) >= abs(b))) / 2001))
    pd.DataFrame(sens).to_csv(OUT / "score_sensitivity.csv", index=False)

    ps = pseudo(donors.index, n_pseudo, jobs, seed)
    ps.to_csv(OUT / "score_pseudo.csv", index=False)
    for t in dl:
        q = float(np.quantile(np.abs(ps[f"beta_{t}"]), .95))
        tab.loc[tab.test == t, "pseudo_abs_q95"] = q
        tab.loc[tab.test == t, "clock_specific"] = bool(abs(tab.loc[tab.test == t, "beta_h_per_decade"].iloc[0]) > q)
    tab.to_csv(OUT / "score_tests.csv", index=False)
    say("done")


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--placebo", action="store_true", help="ages shuffled over donors (code test; output/placebo/)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    global OUT
    OUT = OUT / "placebo" if args.placebo else OUT / "quick" if run["quick"] else OUT
    OUT.mkdir(parents=True, exist_ok=True)
    pbs, donors = load()
    if run["quick"]:
        pbs, donors = subsample(pbs, donors, args.subsample, args.seed)
    if args.placebo:
        perm = np.random.default_rng(seed_of("A234placebo", args.seed)).permutation(len(donors))
        donors[["age_exact", "age"]] = donors[["age_exact", "age"]].to_numpy()[perm]
    n_perm = 200 if run["quick"] else N_PERM
    say("donors body", pbs["body"].donor.nunique(), "brain", pbs["brain"].donor.nunique())
    if args.stage == "power":
        power(pbs, donors, n_perm, args.jobs, args.seed)
    else:
        score(pbs, donors, n_perm, 20 if run["quick"] else N_BOOT, 4 if run["quick"] else N_PSEUDO, args.jobs, args.seed)


if __name__ == "__main__":
    main()
