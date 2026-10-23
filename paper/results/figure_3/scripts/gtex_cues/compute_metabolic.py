"""A47: body mass and type 2 diabetes in GTEx. Do the metabolic organs follow the shared clock less in heavier or
diabetic donors?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A47; run note there).
  H47.1 (primary): module {Liver, subcutaneous and visceral adipose, Pancreas, Skeletal muscle} minus companion
        {Tibial nerve, skin x2}: Re(w/z) per SD of BMI < 0, one-sided p < 0.05. Fitted separately in Hardy 0 and Hardy
        1-2 (own hosts, `Dataset.covariate` per tissue, bias-corrected `locked_amp_c`), pooled by inverse null
        variance. The null here is a permutation of BMI over donors (see README: deviations).
  Secondary: T2D case-control coupling loss (control minus case; own host spread d and common spread d*) per tissue and
  module with the companion (`Dataset.contrast`), pooled over the two cohorts; liver disease (descriptive, counts).
  Sensitivities: BMI without T2D donors, without long-term steroids. Aggregates only.

    python compute_metabolic.py [--subsample N] [--jobs 2] [--n-perm 100]
Outputs: paper/results/figure_3/output/gtex_cues/a47_{design_counts,bmi_tissues,bmi_cohort_modules,bmi_pooled,
t2d_tissues,t2d_cohort_modules,t2d_pooled,primary}.csv
Legacy: legacy/desync/validation/gtex_covariates/scripts/a47_metabolic.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import _covariates as K  # noqa: E402
import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

MET = ["Liver", "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Pancreas", "Muscle - Skeletal"]
COMP = ["Nerve - Tibial", "Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"]
TISS = MET + COMP
MIN_DONORS, MIN_CASES = 20, 5
COH = {"vent": (0,), "fast": (1, 2)}


def donors_ok(ph, ds, variant):
    q = ph.reindex(ds.donor_ids)
    ok = pd.Series(True, index=q.index)
    if variant == "excl_t2d":
        ok &= q.t2d != 1
    if variant == "excl_steroid":
        ok &= q.steroid != 1
    return q.index[ok.to_numpy()]


def tissue_donors(ds, t, ok):
    return np.intersect1d(sorted(set(ds.table.resection[ds.table.compartment == t])), ok)


def bmi_fits(ds, ph, variant, n_null, seed) -> dict:
    """Per tissue: standardised BMI over the tissue's donors, `Dataset.covariate` (no further centring)."""
    ok = donors_ok(ph, ds, variant)
    out = {}
    for i, t in enumerate(TISS):
        b = ph.bmi.reindex(tissue_donors(ds, t, ok)).dropna()
        if len(b) < MIN_DONORS:
            continue
        r = ds.covariate((b - b.mean()) / b.std(), targets=[t], center=False, n_null=n_null, seed=seed + i)
        if t in r.index and "locked_amp" in r.columns:
            out[t] = r.loc[t]
    return out


def bmi_perm(ds, ph, seed) -> dict:
    limit_threads(1)
    rng = np.random.default_rng(seed)
    b = ph.bmi.reindex(donors_ok(ph, ds, "primary")).dropna()
    x = pd.Series(rng.permutation(((b - b.mean()) / b.std()).to_numpy()), index=b.index)
    r = ds.covariate(x, targets=TISS, center=False, n_null=0)
    return {c: K.module({t: r.loc[t, c] for t in r.index if r.loc[t].get("n", 0) >= MIN_DONORS}, MET, COMP)
            for c in ("locked_amp", "drift_h")}


def t2d_labels(ph, ds):
    t = ph.t2d.reindex(ds.donor_ids).dropna()
    return t.map({0.0: "control", 1.0: "case"})


def t2d_perm(ds, lab, seed) -> dict:
    limit_threads(1)
    rng = np.random.default_rng(seed)
    pl = pd.Series(rng.permutation(lab.to_numpy()), index=lab.index)
    return t2d_stats(ds.contrast(pl, targets=TISS, n_perm=0, n_boot=0).table)


def t2d_stats(tab) -> dict:
    """Module statistics of the coupling loss (control minus case) from a contrast table."""
    out = {}
    for q, col in (("d", "d"), ("dstar", "d_star")):
        v = {t: -tab.loc[t, col] for t in tab.index if np.isfinite(tab.loc[t].get(col, np.nan))}
        out[q, "module mean"] = float(np.mean([v[t] for t in MET if t in v])) if any(t in v for t in MET) else np.nan
        out[q, "module minus companion"] = K.module(v, MET, COMP)
    return out


def pool(stats_by_cohort: dict, nulls: dict, upper: bool) -> dict:
    """Inverse-variance pooling of cohort statistics; null = weighted sum of independently shuffled cohort nulls."""
    rng = np.random.default_rng(471)
    cs_ = [c for c in stats_by_cohort if np.isfinite(stats_by_cohort[c]) and len(nulls[c]) > 1]
    w = np.array([1 / np.var(nulls[c]) for c in cs_])
    w = w / w.sum()
    D = float(sum(wi * stats_by_cohort[c] for wi, c in zip(w, cs_)))
    nl = sum(wi * rng.permutation(np.asarray(nulls[c])) for wi, c in zip(w, cs_))
    return dict(D=D, null_sd=float(np.std(nl)), p=K.p_upper(D, nl) if upper else K.p_lower(D, nl),
                weights=";".join(f"{c}:{wi:.3f}" for c, wi in zip(cs_, w)))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(K.FOLDER, run["quick"])
    n_perm = args.n_perm or (8 if run["quick"] else 100)
    ph = K.protected()
    ds = {}
    for c, h in COH.items():
        pb, _ = K.cohort(h, subsample=args.subsample, seed=args.seed)
        ds[c] = cs.Dataset(pb, anchors=list(S.TISSUES), name=f"gtex_{c}")

    # ---- design counts ---------------------------------------------------------------------------------------------
    cnt = []
    for c, d in ds.items():
        ok = donors_ok(ph, d, "primary")
        for t in TISS:
            dn = tissue_donors(d, t, ok)
            cnt.append(dict(cohort=c, tissue=t, n_bmi=int(ph.bmi.reindex(dn).notna().sum()),
                            n_t2d=int(ph.t2d.reindex(dn).notna().sum()), n_t2d_case=int((ph.t2d.reindex(dn) == 1).sum()),
                            n_liver_dis_case=int((ph.liver_dis.reindex(dn) == 1).sum())))
    pd.DataFrame(cnt).to_csv(out / "a47_design_counts.csv", index=False)

    # ---- BMI: tissues, cohort modules, pooled ----------------------------------------------------------------------
    rows, mods = [], []
    for vi, v in enumerate(("primary", "excl_t2d", "excl_steroid")):
        for ci, (c, d) in enumerate(ds.items()):
            fits = bmi_fits(d, ph, v, run["n_null"], 4700 + 100 * vi + 10 * ci)
            for col in ("locked_amp", "drift_h"):
                vals = {t: r[f"{col}_c"] for t, r in fits.items()}
                mods.append(dict(variant=v, estimand=col, cohort=c, D=K.module(vals, MET, COMP),
                                 n_met=sum(t in vals for t in MET), n_comp=sum(t in vals for t in COMP)))
                for t, r in fits.items():
                    rows.append(dict(variant=v, estimand=col, cohort=c, tissue=t, n=int(r["n"]), est=r[col],
                                     est_c=r[f"{col}_c"], p_two=r[f"{col}_p"], coupling=r["coupling"]))
    T, M = pd.DataFrame(rows), pd.DataFrame(mods)
    # permutation null of the module (BMI shuffled over donors), per cohort, primary design
    perm = {}
    for c, d in ds.items():
        r = Parallel(n_jobs=args.jobs)(delayed(bmi_perm)(d, ph, 4710 + 1000 * (c == "fast") + i) for i in range(n_perm))
        perm[c] = {col: np.array([x[col] for x in r], float) for col in ("locked_amp", "drift_h")}
        print("BMI permutations", c, {k: round(float(np.nanstd(v)), 4) for k, v in perm[c].items()}, flush=True)
    pooled = []
    for v in ("primary", "excl_t2d", "excl_steroid"):
        for col in ("locked_amp", "drift_h"):
            s = {c: M[(M.variant == v) & (M.estimand == col) & (M.cohort == c)].D.iloc[0] for c in ds}
            for name, cs_ in (("pooled", list(ds)), ("vent only", ["vent"]), ("fast only", ["fast"])):
                sub = {c: s[c] for c in cs_}
                r = pool(sub, {c: perm[c][col][np.isfinite(perm[c][col])] for c in cs_}, upper=False)
                pooled.append(dict(variant=v, estimand=col, set=name, D=r["D"], null_sd=r["null_sd"],
                                   p_lower=r["p"] if v == "primary" else np.nan, weights=r["weights"]))
    P = pd.DataFrame(pooled)

    # ---- T2D: coupling loss with its companion --------------------------------------------------------------------
    t2d_t, t2d_m, tn = [], [], {}
    obs_stats, null_stats = {}, {}
    for c, d in ds.items():
        lab = t2d_labels(ph, d)
        con = d.contrast(lab, targets=TISS, companion=COMP, n_perm=n_perm, n_boot=0, blocks="none")
        tab = con.table
        for t in tab.index:
            if not np.isfinite(tab.loc[t].get("d", np.nan)):
                continue
            t2d_t.append(dict(cohort=c, tissue=t, n_control=int(tab.loc[t, "donors_control"]),
                              n_case=int(tab.loc[t, "donors_case"]), loss=-tab.loc[t, "d"], loss_star=-tab.loc[t, "d_star"],
                              p_loss=tab.loc[t, "p_less"], p_loss_star=np.nan, V0=tab.loc[t, "V0"],
                              c_control=tab.loc[t, "c_control"], c_case=tab.loc[t, "c_case"]))
        obs_stats[c] = t2d_stats(tab)
        r = Parallel(n_jobs=args.jobs)(delayed(t2d_perm)(d, lab, 4720 + 1000 * (c == "fast") + i) for i in range(n_perm))
        null_stats[c] = {k: np.array([x[k] for x in r], float) for k in obs_stats[c]}
        for k, v in obs_stats[c].items():
            nl = null_stats[c][k][np.isfinite(null_stats[c][k])]
            t2d_m.append(dict(cohort=c, quantity=k[0], module=k[1], D=v, null_sd=float(nl.std()), p_upper=K.p_upper(v, nl)))
        if con.companion is not None:
            con.companion.assign(cohort=c).to_csv(out / f"a47_t2d_companion_{c}.csv", index=False)
    t2d_p = []
    for k in obs_stats["vent"]:
        r = pool({c: obs_stats[c][k] for c in ds}, {c: null_stats[c][k][np.isfinite(null_stats[c][k])] for c in ds}, True)
        t2d_p.append(dict(quantity=k[0], module=k[1], D=r["D"], null_sd=r["null_sd"], p_upper=r["p"],
                          lo=r["D"] - 1.96 * r["null_sd"], hi=r["D"] + 1.96 * r["null_sd"], weights=r["weights"]))

    h = P[(P.variant == "primary") & (P.estimand == "locked_amp") & (P.set == "pooled")].iloc[0]
    prim = pd.DataFrame([dict(hypothesis="H47.1", statistic=h.D, lo=h.D - 1.96 * h.null_sd, hi=h.D + 1.96 * h.null_sd,
                              p=h.p_lower, null_sd=h.null_sd, passed=bool(h.D < 0 and h.p_lower < .05))])
    T.to_csv(out / "a47_bmi_tissues.csv", index=False)
    M.to_csv(out / "a47_bmi_cohort_modules.csv", index=False)
    P.to_csv(out / "a47_bmi_pooled.csv", index=False)
    pd.DataFrame(t2d_t).to_csv(out / "a47_t2d_tissues.csv", index=False)
    pd.DataFrame(t2d_m).to_csv(out / "a47_t2d_cohort_modules.csv", index=False)
    pd.DataFrame(t2d_p).to_csv(out / "a47_t2d_pooled.csv", index=False)
    prim.to_csv(out / "a47_primary.csv", index=False)
    pd.set_option("display.width", 220)
    print(prim.round(4).to_string(index=False), "\n", P.round(4).to_string(index=False), "\n",
          pd.DataFrame(t2d_p).round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
