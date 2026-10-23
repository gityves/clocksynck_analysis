"""A46: ventilation dose-response in GTEx ventilated donors (Hardy 0). Does the body's clock drift from clock time and
do the metabolic organs come loose from the rest as ventilation lasts?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A46; run note there). Hypotheses:
  H46.1 (internal): module {Pancreas, Skeletal muscle, Liver} minus companion {Tibial nerve, skin x2}: Re(w/z), the
        relative change of the host-locked amplitude per day of ventilation (x = days capped at 14, centred), < 0,
        one-sided p < 0.05. Estimand per tissue: `Dataset.covariate` on the leave-one-out host of the other 15 tissues,
        bias-corrected by the median of its simulated nulls (`locked_amp_c`).
  H46.2 (external): per donor cos(host phase - time of death - common offset); slope on x < 0, one-sided permutation
        p (hours permuted within cause-of-death group).
  Holm over the two. Variants: primary, cap7, log2, excl_steroid_sepsis. Secondary: drift (h/day), own noise gamma,
  host spread and alignment by tertile of days. Aggregates only (protected phenotypes held in memory).

    python compute_ventilation_dose.py [--subsample N] [--jobs 2] [--n-perm 100] [--n-boot 50]
Outputs: paper/results/figure_3/output/gtex_cues/a46_{tissues,modules,h462,primary,tertiles,permutation,
fast_context}.csv
Legacy: legacy/desync/validation/gtex_covariates/scripts/a46_ventilation.py (+ common.py)
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

TARGET = ["Pancreas", "Muscle - Skeletal", "Liver"]
COMP = ["Nerve - Tibial", "Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"]
TISS = TARGET + COMP
VARIANTS = ["primary", "cap7", "log2", "excl_steroid_sepsis"]
MIN_DONORS = 30
COLS = ("locked_amp", "drift_h", "gamma")


def make_x(hours: pd.Series, variant: str) -> pd.Series:
    """Ventilation days (capped at 14, 7 for cap7), centred at the median; log2(1 + days) standardised for log2."""
    days = np.minimum(hours / 24.0, 7.0 if variant == "cap7" else 14.0)
    if variant == "log2":
        lg = np.log2(1 + days)
        return (lg - lg.mean()) / lg.std()
    return days - days.median()


def eligible(ph: pd.DataFrame, variant: str, donors) -> pd.Index:
    q = ph.reindex(donors)
    ok = q.hours.notna() & (q.vent != 0)
    if variant == "excl_steroid_sepsis":
        ok &= (q.steroid != 1) & (q.sepsis != 1)
    return q.index[ok.to_numpy()]


def tissue_fits(ds: cs.Dataset, ph, variant, hours, n_null, seed) -> dict:
    """Per tissue: `Dataset.covariate` with that tissue's own design (x centred over the donors that have it)."""
    ok = eligible(ph, variant, ds.donor_ids)
    hrs = (ph.hours if hours is None else hours).reindex(ok)
    out = {}
    for i, t in enumerate(TISS):
        dn = np.intersect1d(sorted(set(ds.table.resection[ds.table.compartment == t])), ok)
        if len(dn) < MIN_DONORS:
            continue
        x = make_x(hrs.reindex(dn), variant)
        r = ds.covariate(x, targets=[t], center=False, n_null=n_null, seed=seed + i)
        if t in r.index and "locked_amp" in r.columns:
            out[t] = r.loc[t]
    return out


def global_fit(ds, ph, hours, variant="primary") -> dict:
    """One call for all tissues with a single design (used inside permutations and the bootstrap)."""
    ok = eligible(ph, variant, ds.donor_ids)
    x = make_x(hours.reindex(ok), variant)
    r = ds.covariate(x, targets=TISS, center=False, n_null=0)
    return {t: r.loc[t, "locked_amp"] for t in r.index if "locked_amp" in r.columns and r.loc[t].get("n", 0) >= MIN_DONORS}


def perm_D(ds, ph, hours, ok, cause, seed) -> float:
    limit_threads(1)
    rng = np.random.default_rng(seed)
    h = pd.Series(hours.reindex(ok).to_numpy()[K.perm_within(cause, rng)], index=ok)
    v = global_fit(ds, ph, h)
    return K.module(v, TARGET, COMP)


def boot_D(pb, ph, seed) -> float:
    """Donor bootstrap (host refitted): donors drawn with replacement; the first copy keeps its id and further copies are
    `<id>#k`, so that the covariate (indexed by id) applies to every copy (local: `Dataset.resample` replicates
    cannot be passed to `covariate`)."""
    limit_threads(1)
    rng = np.random.default_rng(seed)
    donors = np.array(sorted(pb.donor.unique()))
    pick = rng.choice(donors, len(donors), replace=True)
    seen, lab = {}, []
    for p in pick:
        seen[p] = seen.get(p, 0) + 1
        lab.append(p if seen[p] == 1 else f"{p}#{seen[p]}")
    parts = []
    g = {d: x for d, x in pb.groupby("donor")}
    for p, l in zip(pick, lab):
        parts.append(g[p].assign(donor=l))
    bd = cs.Dataset(pd.concat(parts, ignore_index=True), anchors=list(S.TISSUES), name="boot")
    q = ph.reindex(donors)
    ok = (q.hours.notna() & (q.vent != 0)).to_numpy()
    x = make_x(pd.Series(q.hours.to_numpy()[ok], index=donors[ok]), "primary")
    try:
        r = bd.covariate(x, targets=TISS, center=False, n_null=0)
    except (ValueError, np.linalg.LinAlgError):
        return np.nan
    return K.module({t: r.loc[t, "locked_amp"] for t in r.index if r.loc[t].get("n", 0) >= MIN_DONORS}, TARGET, COMP)


def cos_stat(phase, tod_h, delta=None):
    ang = phase - S.W * tod_h
    if delta is None:
        delta = np.angle(np.exp(1j * ang).mean())
    return np.cos(ang - delta), delta


def slope(a, x) -> float:
    xc = x - x.mean()
    return float((xc * a).sum() / (xc ** 2).sum())


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(K.FOLDER, run["quick"])
    n_perm = args.n_perm or (10 if run["quick"] else 100)
    n_boot = args.n_boot or (6 if run["quick"] else 50)
    n_h2 = 200 if run["quick"] else 1000
    ph = K.protected()

    pb, don = K.cohort((0,), subsample=args.subsample, seed=args.seed)
    ds = cs.Dataset(pb, anchors=list(S.TISSUES), name="gtex_vent")
    ok0 = eligible(ph, "primary", ds.donor_ids)
    print({v: len(eligible(ph, v, ds.donor_ids)) for v in VARIANTS}, "donors analysed", flush=True)

    # ---- H46.1: per tissue and module, observed (bias-corrected) -------------------------------------------------
    rows, mods = [], []
    for vi, v in enumerate(VARIANTS):
        fits = tissue_fits(ds, ph, v, None, run["n_null"], 4600 + 10 * vi)
        for col in COLS:
            vals = {t: r[f"{col}_c"] for t, r in fits.items()}
            mods.append(dict(variant=v, estimand=col, D=K.module(vals, TARGET, COMP),
                             n_target=sum(t in vals for t in TARGET), n_comp=sum(t in vals for t in COMP)))
            for t, r in fits.items():
                rows.append(dict(variant=v, estimand=col, tissue=t, n=int(r["n"]), est=r[col], est_c=r[f"{col}_c"],
                                 p_two=r[f"{col}_p"], coupling=r["coupling"]))
        print(v, pd.DataFrame(mods)[lambda d: d.variant == v].D.round(4).tolist(), flush=True)
    T, M = pd.DataFrame(rows), pd.DataFrame(mods)

    # ---- within-cause permutation and donor bootstrap of the primary module ---------------------------------------
    cause = ph.cause.reindex(ok0).to_numpy()
    D_unc = K.module(global_fit(ds, ph, ph.hours), TARGET, COMP)
    Dp = np.array(Parallel(n_jobs=args.jobs)(delayed(perm_D)(ds, ph, ph.hours, ok0, cause, 4620 + i)
                                             for i in range(n_perm)), float)
    perm = dict(D_obs_uncorrected=D_unc, p_lower=K.p_lower(D_unc, Dp), perm_mean=float(np.nanmean(Dp)),
                perm_sd=float(np.nanstd(Dp)), n=int(np.isfinite(Dp).sum()))
    Db = np.array(Parallel(n_jobs=args.jobs)(delayed(boot_D)(pb, ph, 4640 + i) for i in range(n_boot)), float)
    Db = Db[np.isfinite(Db)]
    bias = M[(M.variant == "primary") & (M.estimand == "locked_amp")].D.iloc[0] - D_unc
    lo, hi = np.quantile(Db + bias, [0.025, 0.975]) if len(Db) > 1 else (np.nan, np.nan)
    M.loc[(M.variant == "primary") & (M.estimand == "locked_amp"), ["boot_lo", "boot_hi", "n_boot", "perm_p_lower"]] = \
        [lo, hi, len(Db), perm["p_lower"]]

    # ---- H46.2: external alignment of the host with the time of death --------------------------------------------
    host_v = S.host_states(pb, S.TISSUES, name="vent")
    tod_v = don.tod.reindex(host_v.index)
    pbf, donf = K.cohort((1, 2), tissues=S.TISSUES, subsample=args.subsample, seed=args.seed, require_tod=True)
    host_f = S.host_states(pbf, S.TISSUES, name="fast")
    fast = pd.DataFrame({"ph": np.arctan2(host_f.y, host_f.x), "tod": donf.tod.reindex(host_f.index)}).dropna()
    cf, delta_fast = cos_stat(fast.ph.to_numpy(), fast.tod.to_numpy())
    ctx = dict(fast_R=float(abs(np.exp(1j * (fast.ph - S.W * fast.tod)).mean())), fast_mean_cos=float(cf.mean()),
               fast_offset_h=float(delta_fast / S.W), n_fast=len(fast))
    rng = np.random.default_rng(463)
    h2 = []
    for v in VARIANTS:
        ok = eligible(ph, v, ds.donor_ids).intersection(tod_v.dropna().index)
        x = make_x(ph.hours.reindex(ok), v).to_numpy()
        ang = np.arctan2(host_v.y.reindex(ok), host_v.x.reindex(ok)).to_numpy()
        td = tod_v.reindex(ok).to_numpy()
        av, delta = cos_stat(ang, td)
        g = ph.cause.reindex(ok).to_numpy()
        s_obs = slope(av, x)
        sp = np.array([slope(av, x[K.perm_within(g, rng)]) for _ in range(n_h2)])
        bs = []
        for _ in range(n_h2):
            i = rng.integers(0, len(x), len(x))
            bs.append(slope(cos_stat(ang[i], td[i])[0], x[i]))
        a2, _ = cos_stat(ang, td, delta_fast)
        s2 = slope(a2, x)
        sp2 = np.array([slope(a2, x[K.perm_within(g, rng)]) for _ in range(n_h2)])
        h2.append(dict(variant=v, n=len(ok), slope=s_obs, boot_lo=np.quantile(bs, .025), boot_hi=np.quantile(bs, .975),
                       p_lower=K.p_lower(s_obs, sp), perm_sd=float(sp.std()), mean_cos=float(av.mean()),
                       offset_h=float(delta / S.W), R=float(abs(np.exp(1j * (ang - S.W * td)).mean())),
                       slope_fastoffset=s2, p_lower_fastoffset=K.p_lower(s2, sp2), mean_cos_fastoffset=float(a2.mean())))
        if v == "primary":
            q = K.tertiles(pd.Series(x, index=ok))
            hv = host_v.reindex(ok)
            tert = []
            for k in range(3):
                m = (q == k).to_numpy()
                tert.append(dict(tertile=f"t{k + 1}", n=int(m.sum()),
                                 days_median=float(np.median(ph.hours.reindex(ok)[m] / 24)),
                                 host_spread=float(np.trace(np.cov(hv[["x", "y"]].to_numpy()[m].T))),
                                 mean_cos=float(av[m].mean())))
    H2 = pd.DataFrame(h2)

    # ---- primaries, Holm over H46.1 and H46.2 ----------------------------------------------------------------------
    m1 = M[(M.variant == "primary") & (M.estimand == "locked_amp")].iloc[0]
    p46 = [perm["p_lower"], float(H2.p_lower.iloc[0])]
    adj = cs.stats.holm(p46)
    prim = pd.DataFrame([
        dict(hypothesis="H46.1", statistic=m1.D, lo=m1.boot_lo, hi=m1.boot_hi, p=p46[0], p_holm=adj[0],
             null_sd=perm["perm_sd"], passed=bool(m1.D < 0 and p46[0] < .05), passed_holm=bool(m1.D < 0 and adj[0] < .05)),
        dict(hypothesis="H46.2", statistic=H2.slope.iloc[0], lo=H2.boot_lo.iloc[0], hi=H2.boot_hi.iloc[0], p=p46[1],
             p_holm=adj[1], null_sd=H2.perm_sd.iloc[0], passed=bool(H2.slope.iloc[0] < 0 and p46[1] < .05),
             passed_holm=bool(H2.slope.iloc[0] < 0 and adj[1] < .05))])
    T.to_csv(out / "a46_tissues.csv", index=False)
    M.to_csv(out / "a46_modules.csv", index=False)
    H2.to_csv(out / "a46_h462.csv", index=False)
    prim.to_csv(out / "a46_primary.csv", index=False)
    pd.DataFrame(tert).to_csv(out / "a46_tertiles.csv", index=False)
    pd.Series(perm).to_csv(out / "a46_permutation.csv", header=["value"])
    pd.Series(ctx).to_csv(out / "a46_fast_context.csv", header=["value"])
    pd.set_option("display.width", 220)
    print(prim.round(4).to_string(index=False), "\n", H2.round(4).to_string(index=False), "\n", pd.Series(perm).round(4),
          "\n", pd.Series(ctx).round(4), flush=True)


if __name__ == "__main__":
    main()
