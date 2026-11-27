"""A403b (POST HOC, DESCRIPTIVE; registered in paper/docs/registered_experiments.md, "A403b"): tumour clock loudness and
inflammation across cancer types.

Same pooled table, P401 cancers, hosts, loudness, readouts and adjustment as A403 (`a403_a405_score.py`).
  D1  per cancer slopes (adjusted) with 2,000 within-cohort patient bootstraps, permutation p, companion range; I^2
  D2  immune-hot slope with the malignant G2M index as an extra covariate; within-cohort r(immune-hot, proliferation)
  D3  core-loop loudness; slope; the 364 three-gene drops of the 14 genes; paired 14 - core on shared bootstraps
  D4  stromal inflammatory score (A363 gene set) on loudness, stroma sum and fibroblasts only
                                                                        -> output/a403b/{d1,d2,d3,d3_drops,d4}.csv

    python a403b.py [--jobs 20] [--genome-jobs 4]                          (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import a403_a405_score as S
from a402_a405 import BASE, CORE, META, MDE_K, U, EXP, cs, log, trim, centred, cancer_hosts, loudness_frame, load_slim, readouts
from clocksync.config import limit_threads

sys.path.insert(0, str(EXP / "A361_A363_volume_mechanism/scripts"))
import vm_common as V  # noqa: E402

OUT = U.OUT / "a403b"
N_BOOT, N_PERM = 2000, 10_000
BASECOV = ["n_anchor", "log_depth"]


def units(df, col, covs, loud_col="loud"):
    """Per cohort (>= 4 patients): stratum-centred readout and loudness residualised on `covs`, z-scored."""
    out = []
    for s, g in df.dropna(subset=[col, loud_col] + covs).groupby("cohort"):
        if len(g) < 4:
            continue
        y = centred(g[col].to_numpy(float), g.stratum)
        l_ = centred(g[loud_col].to_numpy(float), g.stratum)
        X = np.column_stack([centred(g[c].to_numpy(float), g.stratum) for c in covs]) if covs else np.empty((len(g), 0))
        if len(g) > X.shape[1] + 2:
            y, l_ = S.residual(y, X), S.residual(l_, X)
        if y.std() < 1e-12 or l_.std() < 1e-12:
            continue
        out.append(dict(cohort=s, n=len(g), y=y / y.std(), l=l_ / l_.std()))
    return out


def slope(df, col, covs, loud_col="loud"):
    un = units(df, col, covs, loud_col)
    return (S.W.pooled_slope(un) if un else np.nan), un


def boot_index(df, rng):
    """Patients resampled with replacement within cohort."""
    idx = [rng.choice(g.index.to_numpy(), len(g)) for _, g in df.groupby("cohort")]
    return np.concatenate(idx)


def describe(df, col, covs, comp, seed, loud_col="loud"):
    obs, un = slope(df, col, covs, loud_col)
    if not un:
        return None
    rng = np.random.default_rng(seed)
    perms = [np.argsort(rng.random((N_PERM, u["n"])), axis=1) for u in un]
    nul = S.W.pooled_slope(un, perms)
    bs = np.array([slope(df.loc[boot_index(df, rng)].reset_index(drop=True), col, covs, loud_col)[0] for _ in range(N_BOOT)])
    r = dict(readout=col, patients=int(sum(u["n"] for u in un)), cohorts=len(un), slope=obs, lo90=float(np.nanquantile(bs, .05)),
             hi90=float(np.nanquantile(bs, .95)), boot_sd=float(np.nanstd(bs)), null_sd=float(nul.std()),
             p_two_sided=float((np.sum(np.abs(nul) >= abs(obs) - 1e-12) + 1) / (len(nul) + 1)))
    for fam, sets in (comp or {}).items():
        sl = []
        for lj in sets:
            d2 = df.drop(columns=loud_col).merge(pd.Series(lj, name=loud_col, dtype=float), left_on="donor", right_index=True)
            sl.append(slope(d2, col, covs, loud_col)[0] if len(d2) else np.nan)
        sl = np.array(sl)[np.isfinite(sl)]
        if len(sl):
            r.update({f"{fam}_n": len(sl), f"{fam}_median": float(np.median(sl)), f"{fam}_q05": float(np.quantile(sl, .05)),
                      f"{fam}_q95": float(np.quantile(sl, .95))})
    return r


def i2(rows):
    s = np.array([r["slope"] for r in rows])
    w = 1 / np.array([r["boot_sd"] for r in rows]) ** 2
    m = (w * s).sum() / w.sum()
    Q = float((w * (s - m) ** 2).sum())
    return dict(Q=Q, k=len(s), I2=max(0.0, (Q - (len(s) - 1)) / Q) if Q > 0 else 0.0, fixed_effect=float(m))


def per_cancer(df, col, covs, comp, tag, loud_col="loud"):
    rows = []
    r = describe(df, col, covs, comp, [4032, sum(map(ord, tag + col))], loud_col)
    if r:
        rows.append(dict(r, cancer="pooled"))
    pc = []
    for c, d in df.groupby("cancer"):
        r = describe(d.reset_index(drop=True), col, covs, comp, [4032, sum(map(ord, tag + col + c))], loud_col)
        if r:
            pc.append(dict(r, cancer=c))
    if len(pc) >= 2:
        h = i2(pc)
        rows[0].update(h)
    return rows + pc


def host_loud(tab, classes, w, stratum):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        h = cs.Dataset(tab, anchors=classes, waveform=w, mixed_strata=True).host()
        return pd.Series(centred(np.log(h.radius.to_numpy()), stratum.reindex(h.index)), index=h.index)
    except Exception:  # noqa: BLE001
        return pd.Series(dtype=float)


def louds_for(hosts, wfs, args):
    """{label: Series donor -> loudness} for each waveform, pooled over cancers."""
    tasks = [(lab, c, h["tab"][META + list(w.genes)], h["classes"], w, h["stratum"]) for lab, w in wfs for c, h in hosts.items()]
    res = Parallel(n_jobs=args.jobs)(delayed(host_loud)(*t[2:]) for t in tasks)
    out = {}
    for t, r in zip(tasks, res):
        out.setdefault(t[0], []).append(r)
    return {k: pd.concat(v) for k, v in out.items()}


def inflammation(hosts):
    """Per donor: A363 inflammatory score in the summed stroma and in fibroblasts (mean within-cohort z of log2 CPM)."""
    rows = []
    for c, h in hosts.items():
        Gc = h["Gc"]
        genes = V.inflammatory_genes([g for g in Gc.columns if g not in META])
        for lab, sub in (("infl_stroma", Gc), ("infl_fibro", Gc[Gc.compartment == "tumour:fibroblast"])):
            if sub.empty:  # no fibroblast class in this cancer's host (e.g. KIRC)
                continue
            t = sub.groupby("donor")[genes + ["depth"]].sum()
            lc = np.log2(t[genes].div(t.depth, axis=0) * 1e6 + 1)
            coh = np.array([d.split("|")[0] for d in lc.index])
            g = lc.groupby(coh)
            z = (lc - g.transform("mean")) / g.transform(lambda x: x.std(ddof=0)).replace(0, np.nan)
            z = z.loc[:, z.notna().mean() > 0.5]
            rows.append(pd.DataFrame({"donor": z.index, "score": z.mean(axis=1).to_numpy(), "kind": lab, "n_genes": len(genes)}))
    R = pd.concat(rows).pivot_table(index="donor", columns="kind", values="score").reset_index()
    return R


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=20)
    ap.add_argument("--genome-jobs", type=int, default=4)
    args = ap.parse_args()
    args.smoke = False
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    frames, donors, f388 = U.pool_sources(args, OUT)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    p401 = pd.read_csv(U.OUT / "power" / "p401.csv")
    cancers = [c for c in p401[(p401.readout.astype(str) == "14") & p401.p401_pass].cancer if c in set(donors.cancer) and c != "unspecified"]
    donors = donors[donors.cancer.isin(cancers)].reset_index(drop=True)
    hosts = cancer_hosts(frames, donors, cancers)
    sources, f388 = list(frames), set(f388)
    frames.clear()
    trim()
    loud = loudness_frame(hosts)
    assert len(hosts) == 12 and len(loud) == 574, "pool differs from A403"
    slim, errs = load_slim(sources, donors, f388, args)
    assert len(slim) == 51 and errs.empty
    ro = readouts(slim, loud)
    df = loud.merge(ro.drop(columns="cohort"), on="donor").merge(S.covariates(hosts), on="donor", how="left")
    df["cohort"] = df.donor.str.split("|").str[0]
    df = df.merge(inflammation(hosts), on="donor", how="left")
    log(f"table: {len(df)} patients; inflammatory score in {df.infl_stroma.notna().sum()}")
    # check: A403 reproduced
    a403 = pd.read_csv(U.OUT / "score_a403_a405" / "a403_tests.csv").query("adjusted and readout == 'immune_hot'").slope.iloc[0]
    rep = slope(df, "immune_hot", BASECOV)[0]
    assert abs(rep - a403) < 1e-9, (rep, a403)
    log("A403 pooled slope reproduced")
    comp = S.companion_louds(hosts, args, 200)
    get_reusable_executor().shutdown(wait=True)
    # D1
    D1 = [r for col in S.ROUT for r in per_cancer(df, col, BASECOV, comp, "d1")]
    pd.DataFrame(D1).to_csv(OUT / "d1.csv", index=False)
    log("D1 done")
    # D2
    D2 = per_cancer(df, "immune_hot", BASECOV + ["prolif"], comp, "d2")
    rr = []
    for s, g in df.dropna(subset=["immune_hot", "prolif"]).groupby("cohort"):
        if len(g) >= 4:
            a, b = centred(g.immune_hot.to_numpy(float), g.stratum), centred(g.prolif.to_numpy(float), g.stratum)
            if a.std() > 0 and b.std() > 0:
                rr.append((len(g), float(np.corrcoef(a, b)[0, 1])))
    n_, r_ = np.array(rr).T
    D2 = pd.DataFrame(D2).assign(r_immune_prolif_within_cohort=float((n_ * r_).sum() / n_.sum()))
    D2.to_csv(OUT / "d2.csv", index=False)
    log("D2 done")
    # D4
    D4 = [dict(r, score=k) for k in ("infl_stroma", "infl_fibro") for r in per_cancer(df, k, BASECOV, comp, "d4")]
    pd.DataFrame(D4).to_csv(OUT / "d4.csv", index=False)
    log("D4 done")
    del comp
    for h in hosts.values():
        h.pop("Gc", None)
    trim()
    # D3: core loop and the 364 three-gene drops
    drops, trio = U.C.drops(BASE)
    L = louds_for(hosts, drops + [("14", BASE)], args)
    assert np.allclose(L["14"].reindex(df.donor).to_numpy(), df.loud.to_numpy(), equal_nan=True)
    vals = {lab: slope(df.drop(columns="loud").merge(L[lab].rename("loud"), left_on="donor", right_index=True), "immune_hot", BASECOV)[0]
            for lab, _ in drops}
    Dd = pd.DataFrame({"drop": list(vals), "slope": list(vals.values())})
    Dd["trio"] = Dd["drop"] == trio
    Dd.to_csv(OUT / "d3_drops.csv", index=False)
    core = vals[trio]
    null = Dd[~Dd.trio].slope.to_numpy()
    df["loud_core"] = L[trio].reindex(df.donor).to_numpy()
    rng = np.random.default_rng([4033])
    pb = []
    for _ in range(N_BOOT):
        d = df.loc[boot_index(df, rng)].reset_index(drop=True)
        a = slope(d, "immune_hot", BASECOV)[0]
        b = slope(d, "immune_hot", BASECOV, "loud_core")[0]
        pb.append((a, b, a - b))
    pb = np.array(pb)
    D3 = pd.DataFrame([dict(slope_14=rep, slope_core=core, core_percentile_towards_zero=float(np.mean(np.abs(null) < abs(core)) * 100),
                            n_drops=len(null), drop_median=float(np.median(null)), drop_q10=float(np.quantile(null, .10)),
                            drop_q90=float(np.quantile(null, .90)), paired_diff=float(rep - core),
                            paired_lo90=float(np.quantile(pb[:, 2], .05)), paired_hi90=float(np.quantile(pb[:, 2], .95)),
                            core_lo90=float(np.quantile(pb[:, 1], .05)), core_hi90=float(np.quantile(pb[:, 1], .95)))])
    D3.to_csv(OUT / "d3.csv", index=False)
    pd.set_option("display.width", 300, "display.max_columns", 40)
    for nm, T in (("D1", pd.DataFrame(D1)), ("D2", D2), ("D3", D3), ("D4", pd.DataFrame(D4))):
        log(f"{nm}\n" + T.round(3).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
