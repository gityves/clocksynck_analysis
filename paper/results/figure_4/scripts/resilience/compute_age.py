"""A19 part 2: does a cell type follow the host clock less with age? Locked amplitude, phase drift and own noise per
decade across the healthy fits of A17, pooled across organs (H19.1, H19.2).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A19 part 2 and the A19 note of 2026-10-13). Eligible fits (rule, at
run time): healthy A17 fits whose gate passes, with anchors (learned; fixed as fallback, flagged), >= 30 donors with an
age spanning >= 30 years. Per compartment on the anchored host (leave-one-out for anchors), x = (age - median) / 10:
locked amplitude Re(w/z), phase drift Im(w/z) (h per decade), own noise gamma, coupling slope, each minus its median
under a null simulated at the real design (`_c`) with that null's two-sided p (`_p`; the primary p of the A19 note);
donor bootstrap SE. Pooling: one value per fit and family (mean over the fit's eligible compartments, mean SE; the two
SEA-AD fits are one study), DerSimonian-Laird across fits.
  H19.1: pooled gamma over all compartments > 0 (two-sided p from the pooled z).
  H19.2: pooled locked amplitude, lymphocyte + myeloid minus fibroblast_stromal + endothelial + mural, < 0
         (within-fit difference, then pooled).
Holm over the two. Secondary: drift, gamma and locked amplitude per family; neurons minus glia (brain, retina); both
primaries without brain.

    python compute_age.py [--subsample N] [--seed S] [--jobs 4]
Outputs (paper/results/figure_4/output/resilience/): a19_age_registry.csv, a19_age_compartments.csv, a19_age_pooled.csv
Legacy: legacy/desync/20_age/scripts/{age_fit,ages,pool}.py. Donor age: the builders' `age` column (exact years or bin
midpoints, as legacy `ages.py`). Aggregates only (no controlled-access data).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "resilience"
ESTIMANDS = ("locked_amp", "drift_h", "gamma", "coupling_slope")
MIN_AGE_DONORS, MIN_SPAN = 30, 30.0
SHARED = {"seaad_dlpfc": "seaad", "seaad_mtg": "seaad"}


def healthy_fit(row, subsample, seed):
    """(pseudobulk, donor table) of a healthy target: donors flagged `healthy` (all donors when the table has no flag)."""
    pb, donors = data.load(row.dataset)
    if "healthy" in donors:
        donors = donors[donors.healthy.astype(bool)]
    pb = pb[pb.donor.isin(donors.index)]
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb = pb[pb.donor.isin(keep)]
        donors = donors.loc[donors.index.isin(keep)]
    return pb.reset_index(drop=True), donors


def boot_se(pb: pd.DataFrame, anchors, x: pd.Series, targets, n_boot: int, seed: int) -> pd.DataFrame:
    """local: not in clocksync (`Dataset.covariate` does not run on `Dataset.resample` replicates). Donor-bootstrap SD of
    the covariate estimands: donors drawn with replacement, repeats renamed `<id>.<j>`, the Dataset (projection, host)
    rebuilt in every replicate (legacy: host fixed, so these SEs are the larger, host-refitted ones)."""
    rng = np.random.default_rng(seed)
    donors = np.array(sorted(pb.donor.unique()))
    by = {d: g for d, g in pb.groupby("donor")}
    reps = []
    for _ in range(n_boot):
        pick = rng.choice(donors, len(donors))
        pbb = pd.concat([by[d].assign(donor=f"{d}.{j}") for j, d in enumerate(pick)], ignore_index=True)
        xb = pd.Series({f"{d}.{j}": x.get(d, np.nan) for j, d in enumerate(pick)}).dropna()
        try:
            r = cs.Dataset(pbb, anchors=anchors).covariate(xb, targets=targets, center=False, n_null=0)
            reps.append(np.asarray(r[list(ESTIMANDS)].reindex(targets), float))
        except (ValueError, KeyError, np.linalg.LinAlgError):
            continue
    a = np.stack(reps) if reps else np.empty((0, len(targets), 4))
    return pd.DataFrame(np.nanstd(a, axis=0, ddof=1), index=targets, columns=[f"{e}_se" for e in ESTIMANDS])


def one_fit(row, run, subsample, seed, min_age):
    limit_threads(1)
    warnings.simplefilter("ignore", RuntimeWarning)
    pb, donors = healthy_fit(row, subsample, seed)
    age = donors["age"].astype(float).dropna() if "age" in donors else pd.Series(dtype=float)
    age = age[age.index.isin(set(pb.donor))]
    reg = dict(name=row.name, tissue=row.tissue, study=row.study, donors=pb.donor.nunique(), donors_age=len(age),
               age_min=age.min(), age_median=age.median(), age_max=age.max())
    rule = "learned" if sh.anchor_set(row.name, "learned") else "fixed" if sh.anchor_set(row.name, "fixed") else None
    reg.update(rule=rule, fallback=rule == "fixed")
    if rule is None or len(age) < min_age or (age.max() - age.min()) < MIN_SPAN:
        return dict(reg, gate=np.nan, eligible=False), None
    anc = [a for a in sh.anchor_set(row.name, rule) if a in set(pb.compartment)]
    ds = cs.Dataset(pb, anchors=anc, donors=donors, name=row.name)
    g = ds.gate(n_perm=run["n_perm"], seed=seed)
    reg.update(gate=bool(g.passes), p_pairing=g.p_pairing, p_acrophase=g.p_acrophase)
    reg["eligible"] = bool(g.passes)
    if not g.passes:
        return reg, None
    x = age / 10.0
    cov = ds.covariate(x, n_null=run["n_null"], seed=seed)
    if "error" in cov:
        cov = cov[cov["error"].isna()]
    cov = cov.apply(pd.to_numeric, errors="coerce")
    cov = cov[cov.n >= min_age]
    se = boot_se(pb, anc, x - x.median(), list(cov.index), run["n_boot"], seed + 1)
    out = cov.join(se).reset_index().assign(target=row.name, rule=rule, tissue=row.tissue)
    out["anchor"] = out.compartment.isin(anc)
    out["family"] = out.compartment.map(sh.family_map())
    return reg, out


def per_fit(df, e, sel):
    """One value per study: mean of `e`_c over the selected compartments of each fit, SE the mean SE."""
    u = df[np.isfinite(df[f"{e}_c"]) & np.isfinite(df[f"{e}_se"])]
    u = u[sel(u)]
    g = u.groupby("target").agg(y=(f"{e}_c", "mean"), se=(f"{e}_se", "mean"), k=(f"{e}_c", "size"),
                                tissue=("tissue", "first")).reset_index()
    g["study"] = g.target.map(SHARED).fillna(g.target)
    return g.groupby("study").agg(y=("y", "mean"), se=("se", "mean"), k=("k", "sum"),
                                  tissue=("tissue", "first")).reset_index()


def diff_fit(df, e, a, b):
    ga, gb = per_fit(df, e, a), per_fit(df, e, b)
    m = ga.merge(gb, on="study", suffixes=("_a", "_b"))
    return pd.DataFrame(dict(study=m.study, y=m.y_a - m.y_b, se=np.sqrt(m.se_a ** 2 + m.se_b ** 2),
                             k=m.k_a + m.k_b, tissue=m.tissue_a))


def pooled(g, label, **kw) -> dict:
    r = cs.stats.meta(g.y, g.se)
    z = r["random"] / r["se_random"] if r["k"] and r["se_random"] > 0 else np.nan
    return dict(test=label, studies=r["k"], compartments=int(g.k.sum()) if len(g) else 0, estimate=r["random"],
                se=r["se_random"], tau2=r["tau2"], p=2 * norm.sf(abs(z)) if np.isfinite(z) else np.nan, **kw)


def pool_all(df: pd.DataFrame) -> pd.DataFrame:
    fam = lambda *f: (lambda u: u.family.isin(f).to_numpy())
    allc = lambda u: np.ones(len(u), bool)
    imm, strv = fam(*sh.IMMUNE), fam(*sh.STROMAL_VASCULAR)
    rows = [pooled(per_fit(df, "gamma", allc), "H19.1 gamma, all compartments", primary=True),
            pooled(diff_fit(df, "locked_amp", imm, strv), "H19.2 locked amplitude, immune - stromal-vascular",
                   primary=True)]
    p = np.array([r["p"] for r in rows], float)
    h = np.full(2, np.nan)
    if np.isfinite(p).all():
        h = cs.stats.holm(p)
    for r, v in zip(rows, h):
        r["p_holm"] = v
    nb = lambda g: g[g.tissue != "brain"]
    rows += [pooled(nb(per_fit(df, "gamma", allc)), "H19.1 without brain", primary=False),
             pooled(nb(diff_fit(df, "locked_amp", imm, strv)), "H19.2 without brain", primary=False)]
    for f in sorted(df.family.dropna().unique()):
        for e in ("drift_h", "gamma", "locked_amp"):
            rows.append(pooled(per_fit(df, e, fam(f)), f"{e}, {f}", primary=False))
    cns = df[df.tissue.isin(["brain", "retina"])]
    for e in ("locked_amp", "gamma", "drift_h"):
        rows.append(pooled(diff_fit(cns, e, fam("neuron"), fam("glia")), f"{e}, neuron - glia", primary=False))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="+", default=None, help="healthy fits (default: all 34)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    min_age = min(MIN_AGE_DONORS, args.subsample // 2) if args.subsample else MIN_AGE_DONORS
    tg = sh.targets("healthy")
    if args.only:
        tg = tg[tg.name.isin(args.only)]
    res = Parallel(n_jobs=args.jobs)(delayed(one_fit)(r, run, args.subsample, args.seed + i, min_age)
                                     for i, r in enumerate(tg.itertuples()))
    reg = pd.DataFrame([r for r, _ in res])
    comp = pd.concat([c for _, c in res if c is not None], ignore_index=True)
    reg.to_csv(out / "a19_age_registry.csv", index=False)
    comp.to_csv(out / "a19_age_compartments.csv", index=False)
    pool = pool_all(comp)
    pool.to_csv(out / "a19_age_pooled.csv", index=False)
    print(reg[reg.eligible][["name", "donors_age", "rule"]].to_string(index=False))
    print(pool.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
