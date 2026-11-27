"""A403, A404, A405 scoring run (registered in paper/docs/registered_experiments.md, "A402-A405", power notes and
"A403-A405 note 1"). A402 closed at the power stage and is not scored.

Same pooled tumour table, P401 cancers, tumour-stroma hosts and slim tables as the power stage (`a402_a405.py`).

  A403  pooled within-cohort slope of each readout on clock loudness (adjusted: stratum-centred, residualised on the
        number of stromal anchors and the log median stromal depth; unadjusted beside), 10,000 within-cohort shuffles,
        companions: the observed readout on the loudness of 200 matched and 200 narrow responsive pseudo-hosts per cancer;
        H403.1 immune-hot (two-sided), H403.2 proliferation (one-sided, slope < 0), Holm; descriptive: shares, per cancer.
                                                                       -> a403_tests.csv, a403_per_cancer.csv, a403_companions.csv
  A404  immune-compartment coupling, amplitude and offset (fibroblast reference) on (a) the A387 normal-stroma table,
        (b) the pooled tumour-stroma table, (b17) the 17-source companion table; companion q95 from the power stage;
        H404.1 / H404.2 (myeloid, 14 genes), Holm; per cancer descriptive.   -> a404_fits.csv, a404_tests.csv, a404_per_cancer.csv
  A405  BIOKEY: pre-treatment loudness, immune-hot score and cancer-cell coupling contribution by expansion (medians,
        bootstrap intervals, label-shuffle MDE); pre vs on change (joint host), depth-adjusted, by expansion.
                                                                       -> a405_patients.csv, a405_pre.csv, a405_change.csv

    python a403_a405_score.py [--jobs 20] [--genome-jobs 4] [--smoke]          (long: paper/common/run.sh)

--smoke: a code check on a small pool with readouts shuffled within cohort (A403), the immune rows' donor labels permuted
within stratum (A404) and the expansion labels and on-treatment patient labels shuffled (A405): no target value is printed.
"""
from __future__ import annotations

import argparse
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import a402_a405 as W
from a402_a405 import (ANCH, BASE, CORE, IMMUNE, MALT, META, MDE_K, U, M, cs, log, trim, centred, pseudo_host,
                       cancer_hosts, loudness_frame, load_slim, readouts, a387_immune_table)
from clocksync.config import limit_threads

OUT = U.OUT / "score_a403_a405"
N_PERM = 10_000
ROUT = ("immune_hot", "prolif", "tnk_share", "myeloid_share")
HYP = {"immune_hot": "H403.1", "prolif": "H403.2", "tnk_share": "descriptive", "myeloid_share": "descriptive"}


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o])
    out = np.empty_like(p)
    out[o] = np.minimum(adj, 1)
    return out


# ---- A403 --------------------------------------------------------------------------------------------------------------

def covariates(hosts):
    """Per donor: number of stromal anchor classes and log median stromal depth (the host's stroma rows)."""
    rows = []
    for h in hosts.values():
        t = h["tab"]
        g = t.groupby("donor")
        rows.append(pd.DataFrame({"n_anchor": g.compartment.nunique(), "log_depth": np.log(g.depth.median())}))
    return pd.concat(rows).rename_axis("donor").reset_index()


def residual(v, X):
    """OLS residual of v on [1, X] (X centred columns; constant columns dropped)."""
    X = X[:, X.std(0) > 0] if X.size else X
    A = np.column_stack([np.ones(len(v)), X])
    return v - A @ np.linalg.lstsq(A, v, rcond=None)[0]


def units(df, col, adjust, loud_col="loud"):
    """Per cohort (>= 4 patients): stratum-centred, optionally covariate-residualised, z-scored readout and loudness."""
    out = []
    for s, g in df.dropna(subset=[col, loud_col]).groupby("cohort"):
        if len(g) < 4:
            continue
        y = centred(g[col].to_numpy(float), g.stratum)
        l = centred(g[loud_col].to_numpy(float), g.stratum)
        if adjust:
            X = np.column_stack([centred(g[c].to_numpy(float), g.stratum) for c in ("n_anchor", "log_depth")])
            if len(g) > X.shape[1] + 2:
                y, l = residual(y, X), residual(l, X)
        if y.std() < 1e-12 or l.std() < 1e-12:
            continue
        out.append(dict(cohort=s, n=len(g), y=y / y.std(), l=l / l.std(), donor=g.donor.to_numpy()))
    return out


def slope_test(un, seed, n_perm=N_PERM):
    obs = W.pooled_slope(un)
    rng = np.random.default_rng(seed)
    perms = [np.argsort(rng.random((n_perm, u["n"])), axis=1) for u in un]
    nul = W.pooled_slope(un, perms)
    return obs, nul


def companion_louds(hosts, args, n_set):
    """{family: [{donor: loudness} per set]} on each cancer's pooled stroma; matched seed [401, 0] (= power), responsive [401, 1]."""
    from clocksync.paired import pseudo_waveforms
    out = {"matched": [dict() for _ in range(n_set)], "responsive": [dict() for _ in range(n_set)]}
    for c, h in hosts.items():
        Gc = h["Gc"]
        genes = [g for g in Gc.columns if g not in META]
        frac = Gc[genes].div(Gc.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(Gc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        for code, fam in enumerate(("matched", "responsive")):
            ws = pseudo_waveforms(lm, n_set, np.random.default_rng([U.SEED, code]), BASE, fam)
            used = sorted({g for w in ws for g in w.genes} | set(BASE.genes))
            res = Parallel(n_jobs=args.jobs)(delayed(safe_host)(Gc[META + used], h["classes"], w, h["stratum"]) for w in ws)
            for j, r in enumerate(res):
                out[fam][j].update(r.to_dict())
        log(f"A403 companion loudness {c}")
    return out


def safe_host(Gu, classes, w, stratum):
    try:
        return pseudo_host(Gu, classes, w, stratum)
    except Exception:  # noqa: BLE001  (a pseudo-host that cannot be fitted: no companion loudness for that set)
        return pd.Series(dtype=float)


def stage_a403(df, comp, out):
    rows, crow = [], []
    for adjust in (True, False):
        for col in ROUT:
            un = units(df, col, adjust)
            obs, nul = slope_test(un, [403, int(adjust), sum(map(ord, col))])
            p2 = float((np.sum(np.abs(nul) >= abs(obs) - 1e-12) + 1) / (len(nul) + 1))
            p_lo = float((np.sum(nul <= obs + 1e-12) + 1) / (len(nul) + 1))
            r = dict(adjusted=adjust, readout=col, hypothesis=HYP[col], patients=int(sum(u["n"] for u in un)), cohorts=len(un),
                     slope=obs, null_sd=float(nul.std()), mde=MDE_K * float(nul.std()), p_two_sided=p2, p_one_sided_negative=p_lo)
            for fam, sets in comp.items():
                sl = []
                for j, lj in enumerate(sets):
                    d2 = df.drop(columns="loud").merge(pd.Series(lj, name="loud", dtype=float), left_on="donor", right_index=True)
                    u2 = units(d2, col, adjust) if len(d2) else []
                    sl.append(W.pooled_slope(u2) if u2 else np.nan)
                    if adjust:
                        crow.append(dict(readout=col, family=fam, set=j, slope=sl[-1]))
                sl = np.array(sl)[np.isfinite(sl)]
                r.update({f"{fam}_n": len(sl), f"{fam}_median": float(np.median(sl)), f"{fam}_q05": float(np.quantile(sl, .05)),
                          f"{fam}_q95": float(np.quantile(sl, .95))})
            rows.append(r)
    T = pd.DataFrame(rows)
    def beyond(r, f):
        return (r[f"{f}_q95"] < r["slope"] or r[f"{f}_q05"] > r["slope"]) if r["readout"] != "prolif" else r[f"{f}_q05"] > r["slope"]
    T["beyond_companions"] = [all(beyond(r, f) for f in comp) for _, r in T.iterrows()]
    prim = T.adjusted & T.readout.isin(["immune_hot", "prolif"])
    pv = np.where(T.readout == "prolif", T.p_one_sided_negative, T.p_two_sided)
    T["p_registered"] = pv
    T["p_holm"] = np.nan
    T.loc[prim, "p_holm"] = holm(pv[prim.to_numpy()])
    T["passes"] = prim & (T.p_holm < 0.05) & T.beyond_companions
    T.to_csv(out / "a403_tests.csv", index=False)
    pd.DataFrame(crow).to_csv(out / "a403_companions.csv", index=False)
    # per cancer, descriptive (adjusted)
    pc = []
    for c, d in df.groupby("cancer"):
        for col in ROUT:
            un = units(d, col, True)
            if not un:
                continue
            obs, nul = slope_test(un, [403, 9, sum(map(ord, c + col))])
            pc.append(dict(cancer=c, readout=col, patients=int(sum(u["n"] for u in un)), cohorts=len(un), slope=obs,
                           null_sd=float(nul.std()), p_two_sided=float((np.sum(np.abs(nul) >= abs(obs) - 1e-12) + 1) / (len(nul) + 1))))
    P = pd.DataFrame(pc)
    P.to_csv(out / "a403_per_cancer.csv", index=False)
    log("A403\n" + T.round(4).to_string(index=False))
    return T


# ---- A404 --------------------------------------------------------------------------------------------------------------

def fit_task(tab, anchors, rd, ref, n_boot, n_perm, tag):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    targets = [t for t in IMMUNE if t in set(tab.compartment)]
    ds = cs.Dataset(tab, anchors=anchors, waveform=U.READOUTS[rd], mixed_strata=True)
    mode = "loo" if len(anchors) >= 3 else "full"
    kw = dict(seed=U.SEED, n_perm=n_perm, gate=False, level=0.90, attenuation="df", anchor_readout=mode,
              targets=targets + ([a for a in anchors] if ref else []), reference=ref)
    try:
        t = ds.fit(n_boot=n_boot, **kw).table
        t["interval_fallback"] = ""
    except ValueError as e:
        t = ds.fit(n_boot=0, **kw).table
        t["interval_fallback"] = f"no intervals ({str(e)[:80]})"
    t = t.loc[[k for k in targets if k in t.index]].reset_index().rename(columns={"index": "compartment"})
    if "compartment" not in t:
        t = t.rename(columns={t.columns[0]: "compartment"})
    return t.assign(table=tag, readout=rd, kind="reference" if ref else "coupling")


def stage_a404(ta, anch_a, tb, tb17, hosts, args, out, n_boot, n_perm):
    tasks = []
    for tag, tab, anc in (("A387 normal-stroma", ta, anch_a), ("tumour-stroma pooled", tb, ANCH), ("tumour-stroma 17 sources", tb17, ANCH)):
        anc = [a for a in anc if tab[tab.compartment == a].donor.nunique() >= 5]
        ref = next(a for a in anc if a.endswith("fibroblast"))
        for rd in ("14", "core"):
            if tag.endswith("17 sources"):
                tasks.append((tab, anc, rd, None, 0, 0, tag))
                continue
            tasks.append((tab, anc, rd, None, n_boot, n_perm, tag))
            tasks.append((tab, anc, rd, ref, n_boot, 0, tag))
    for c, h in hosts.items():
        ids = set(h["tab"].donor)
        tab = tb[tb.donor.isin(ids)]
        if not tab.compartment.isin(IMMUNE).any():
            continue
        for rd in ("14", "core"):
            tasks.append((tab, h["classes"], rd, None, min(n_boot, 200), min(n_perm, 200), f"cancer {c}"))
    res = Parallel(n_jobs=min(args.jobs, len(tasks)))(delayed(fit_task)(*t) for t in tasks)
    F = pd.concat(res, ignore_index=True)
    F.to_csv(out / "a404_fits.csv", index=False)
    pw = pd.read_csv(W.OUT / "a404_power.csv")
    pw["key"] = np.where(pw.host.str.startswith("A387"), "A387", "A401")
    q = pw.set_index(["key", "readout", "target"])
    rows = []
    for tag, key in (("A387 normal-stroma", "A387"), ("tumour-stroma pooled", "A401")):
        for rd in ("14", "core"):
            for t in IMMUNE:
                f = F[(F.table == tag) & (F.readout == rd) & (F.compartment == t) & (F.kind == "coupling")]
                g = F[(F.table == tag) & (F.readout == rd) & (F.compartment == t) & (F.kind == "reference")]
                if f.empty:
                    continue
                f, qq = f.iloc[0], q.loc[(key, rd, t)]
                r = dict(host=tag, readout=rd, target=t, n=int(f.n), coupling=float(f.coupling), c_lo90=f.coupling_lo, c_hi90=f.coupling_hi,
                         p_coupling=float(f.p_coupling), host_reliability=f.get("host_reliability", np.nan),
                         q95_matched=float(qq.q95_matched), q95_responsive=float(qq.q95_responsive),
                         median_matched=float(qq.median_matched), median_responsive=float(qq.median_responsive),
                         interval_fallback=f.interval_fallback)
                if not g.empty:
                    g = g.iloc[0]
                    r.update(amp_rel=float(np.exp(g.log_amp_rel)), amp_lo=float(np.exp(g.log_amp_lo)), amp_hi=float(np.exp(g.log_amp_hi)),
                             offset_h=float(g.offset_h), offset_lo=float(g.offset_lo), offset_hi=float(g.offset_hi))
                if key == "A401":
                    f17 = F[(F.table == "tumour-stroma 17 sources") & (F.readout == rd) & (F.compartment == t)]
                    r["coupling_17_sources"] = float(f17.coupling.iloc[0]) if len(f17) else np.nan
                    r["n_17_sources"] = int(f17.n.iloc[0]) if len(f17) else 0
                comp_c = r.get("coupling_17_sources", r["coupling"])
                r["above_both_q95"] = bool(comp_c > r["q95_matched"] and comp_c > r["q95_responsive"])
                r["full_above_both_q95"] = bool(r["coupling"] > r["q95_matched"] and r["coupling"] > r["q95_responsive"])
                r["hypothesis"] = {("A387", "14", "tumour:myeloid"): "H404.1", ("A401", "14", "tumour:myeloid"): "H404.2"}.get((key, rd, t), "descriptive")
                rows.append(r)
    T = pd.DataFrame(rows)
    h = T.hypothesis.str.startswith("H404")
    T["p_holm"] = np.nan
    T.loc[h, "p_holm"] = holm(T.loc[h, "p_coupling"].to_numpy())
    T["passes"] = h & (T.p_holm < 0.05) & T.above_both_q95
    T.to_csv(out / "a404_tests.csv", index=False)
    F[F.table.str.startswith("cancer")].to_csv(out / "a404_per_cancer.csv", index=False)
    log("A404\n" + T.round(3).to_string(index=False))
    return T


# ---- A405 --------------------------------------------------------------------------------------------------------------

def bassez_on(genes):
    """BIOKEY on-treatment malignant and stromal rows (>= 50 cells), donor `<patient>|On`."""
    G = pd.read_parquet(U.BASSEZ)
    G = G[G.timepoint == "On"].drop(columns=[c for c in ("timepoint", "expansion", "cohort", "BC_type") if c in G])
    G = U.tidy(G, only=genes)
    G, _ = U.finish(G, "bassez2021", {d: "BRCA" for d in set(G.donor)})
    G["donor"] = G.donor + "|On"
    return G


def contribution(stab, host):
    """Per donor Re(m conj(b h)) / mean |m|^2 (malignant state m, host h, b = sum m conj(h) / sum |h|^2)."""
    from clocksync.deviation import states
    s = states(stab[stab.compartment == MALT], BASE)
    m = pd.Series(s.to_numpy(), index=s.index.get_level_values("donor"))
    hz = pd.Series((host.x + 1j * host.y).to_numpy(), index=host.index.astype(str))
    ids = m.index.intersection(hz.index)
    m, hz = m[ids].to_numpy(), hz[ids].to_numpy()
    b = (m * np.conj(hz)).sum() / (np.abs(hz) ** 2).sum()
    return pd.Series(np.real(m * np.conj(b * hz)) / np.mean(np.abs(m) ** 2), index=ids)


def boot_median(y, rng, n=2000):
    b = np.array([np.median(rng.choice(y, len(y))) for _ in range(n)])
    return float(np.median(y)), float(np.quantile(b, .05)), float(np.quantile(b, .95))


def stage_a405(loud, ro, hosts, slim, args, out, smoke):
    meta = pd.read_parquet(U.BASSEZ, columns=["donor", "timepoint", "compartment", "depth", "n_cells", "expansion", "cohort"])
    meta = meta[meta.n_cells >= U.MIN_CELLS]
    lab = meta.drop_duplicates("donor").set_index("donor")[["expansion", "cohort"]]
    rng0 = np.random.default_rng([405, 0])
    if smoke:
        lab["expansion"] = lab.expansion.to_numpy()[rng0.permutation(len(lab))]
    br = hosts["BRCA"]
    mal_pre = pd.concat([G[(G.compartment == MALT) & G.donor.isin(set(br["tab"].donor))][META + list(BASE.genes)]
                         for G in slim.values() if set(BASE.genes) <= set(G.columns)], ignore_index=True)
    c_pre = contribution(pd.concat([br["tab"], mal_pre], ignore_index=True), br["host"])
    P = loud[loud.cancer == "BRCA"].set_index("donor")[["loud"]].join(ro.set_index("donor")[["immune_hot"]]).join(c_pre.rename("contribution"))
    P = P[P.index.str.startswith("bassez2021|")].join(lab)
    P.rename_axis("donor").reset_index().drop(columns=[]).to_csv(out / "a405_patients.csv", index=False)  # no per-donor clinical data beyond the public labels
    rows = []
    rng = np.random.default_rng([405, 1])
    for q in ("loud", "immune_hot", "contribution"):
        d = P.dropna(subset=[q])
        d = d[d.expansion.isin(["E", "NE"])]
        r = dict(quantity=q, n_E=int((d.expansion == "E").sum()), n_NE=int((d.expansion == "NE").sum()))
        for g in ("E", "NE"):
            r[f"median_{g}"], r[f"lo90_{g}"], r[f"hi90_{g}"] = boot_median(d[d.expansion == g][q].to_numpy(), rng)
        y, e = d[q].to_numpy(), (d.expansion == "E").to_numpy()
        r["mean_diff_E_minus_NE"] = float(y[e].mean() - y[~e].mean())
        sh = np.array([y[(p := rng.permutation(e))].mean() - y[~p].mean() for _ in range(2000)])
        r["null_sd"], r["mde"] = float(sh.std()), MDE_K * float(sh.std())
        r["shuffle_share_abs_ge"] = float(np.mean(np.abs(sh) >= abs(r["mean_diff_E_minus_NE"])))
        rows.append(r)
    Pre = pd.DataFrame(rows)
    Pre.to_csv(out / "a405_pre.csv", index=False)
    # (2) pre vs on: joint host on the BRCA stroma with the on-treatment rows
    on = bassez_on(set(BASE.genes))
    on = on[on.donor.str[:-3].isin(set(br["tab"].donor))]
    if smoke:  # on-treatment patient labels shuffled
        ids = sorted(set(on.donor))
        mp = dict(zip(ids, np.array(ids)[rng0.permutation(len(ids))]))
        on["donor"] = on.donor.map(mp)
    st_on = on[on.compartment.isin(br["classes"])][META + list(BASE.genes)]
    tab = pd.concat([br["tab"], st_on], ignore_index=True)
    host = cs.Dataset(tab, anchors=br["classes"], waveform=BASE, mixed_strata=True).host()
    strat = tab.drop_duplicates("donor").set_index("donor").stratum.reindex(host.index)
    L_ = pd.Series(centred(np.log(host.radius.to_numpy()), strat), index=host.index)
    cj = contribution(pd.concat([tab, mal_pre, on[on.compartment == MALT][META + list(BASE.genes)]], ignore_index=True), host)
    dep = pd.concat([br["tab"], st_on]).groupby("donor").depth.sum()
    pts = [d for d in P.index if d + "|On" in L_.index and d in L_.index]
    C = pd.DataFrame({"d_loud": [L_[p + "|On"] - L_[p] for p in pts],
                      "d_contribution": [cj.get(p + "|On", np.nan) - cj.get(p, np.nan) for p in pts],
                      "ldr": [np.log(dep[p + "|On"] / dep[p]) for p in pts]}, index=pd.Index(pts, name="donor")).join(lab)
    rows = []
    for q in ("d_loud", "d_contribution"):
        d = C.dropna(subset=[q])

        def adj(d):
            X = np.column_stack([np.ones(len(d)), d.ldr.to_numpy()] + [(d.cohort == k).to_numpy(float) for k in sorted(set(d.cohort))[1:]])
            beta = np.linalg.lstsq(X, d[q].to_numpy(), rcond=None)[0]
            a = d[q].to_numpy() - beta[1] * d.ldr.to_numpy()  # change at equal stromal depth
            e = (d.expansion == "E").to_numpy()
            ne = (d.expansion == "NE").to_numpy()
            return dict(all=a.mean(), E=a[e].mean() if e.any() else np.nan, NE=a[ne].mean() if ne.any() else np.nan,
                        E_minus_NE=(a[e].mean() - a[ne].mean()) if e.any() and ne.any() else np.nan, depth_slope=beta[1])
        est = adj(d)
        bs = pd.DataFrame([adj(d.iloc[rng.integers(0, len(d), len(d))]) for _ in range(2000)])
        for k in ("all", "E", "NE", "E_minus_NE"):
            rows.append(dict(quantity=q, group=k, n=len(d) if k in ("all", "E_minus_NE") else int((d.expansion == k).sum()),
                             estimate=est[k], lo90=float(bs[k].quantile(.05)), hi90=float(bs[k].quantile(.95)), depth_slope=est["depth_slope"]))
        e = d.expansion.isin(["E", "NE"])
        de = d[e]
        a = de[q].to_numpy()
        lab_e = (de.expansion == "E").to_numpy()
        sh = np.array([a[(p := rng.permutation(lab_e))].mean() - a[~p].mean() for _ in range(2000)])
        rows[-1].update(null_sd=float(sh.std()), mde=MDE_K * float(sh.std()))
    Ch = pd.DataFrame(rows)
    Ch.to_csv(out / "a405_change.csv", index=False)
    log("A405 pre\n" + Pre.round(3).to_string(index=False))
    log("A405 change\n" + Ch.round(3).to_string(index=False))


# ---- main --------------------------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=20)
    ap.add_argument("--genome-jobs", type=int, default=4)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    out = OUT.parent / ("score_a403_a405_smoke" if args.smoke else OUT.name)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n_set, n_boot, n_perm = (6, 20, 50) if args.smoke else (200, 500, 1000)
    frames, donors, f388 = U.pool_sources(args, out)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    p401 = pd.read_csv(U.OUT / "power" / "p401.csv")
    cancers = [c for c in p401[(p401.readout.astype(str) == "14") & p401.p401_pass].cancer if c in set(donors.cancer) and c != "unspecified"]
    log(f"P401 cancers: {cancers}")
    donors = donors[donors.cancer.isin(cancers)].reset_index(drop=True)
    hosts = cancer_hosts(frames, donors, cancers)
    sources, f388 = list(frames), set(f388)
    frames.clear()
    trim()
    loud = loudness_frame(hosts)
    log(f"hosts: {len(hosts)} cancers, {len(loud)} patients with a host")
    if not args.smoke:
        assert len(hosts) == 12 and len(loud) == 574, "pool differs from the power stage"
    slim, errs = load_slim(sources, donors, f388, args)
    log(f"slim tables: {len(slim)} sources, errors: {len(errs)}")
    if not args.smoke:
        assert len(slim) == 51 and errs.empty, "slim tables differ from the power stage"
    # A403
    ro = readouts(slim, loud)
    df = loud.merge(ro.drop(columns="cohort"), on="donor").merge(covariates(hosts), on="donor", how="left")
    df["cohort"] = df.donor.str.split("|").str[0]
    if args.smoke:
        rng = np.random.default_rng(1)
        for col in ROUT:
            df[col] = df.groupby("cohort")[col].transform(lambda x: x.to_numpy()[rng.permutation(len(x))])
    comp = companion_louds(hosts, args, n_set)
    get_reusable_executor().shutdown(wait=True)
    stage_a403(df, comp, out)
    del comp
    # A405 (before the genome-wide stroma tables are dropped)
    stage_a405(loud, ro, hosts, slim, args, out, args.smoke)
    # A404
    ta, anch_a = a387_immune_table(args.smoke)
    stroma = pd.concat([h["tab"] for h in hosts.values()], ignore_index=True)
    imm = pd.concat([G[G.compartment.isin(IMMUNE) & G.donor.isin(set(stroma.donor))][META + list(BASE.genes)] for G in slim.values()
                     if all(g in G.columns for g in BASE.genes)], ignore_index=True)
    census = {p.name[:-8] for p in U.CEN.glob("*.parquet")}
    src17 = {s for s in slim if s in U.GEO_CANCER or s == "bassez2021" or s in census}
    for h in hosts.values():
        h.pop("Gc", None)
    trim()
    if args.smoke:  # immune rows' donor labels permuted within stratum
        rng = np.random.default_rng(2)
        for T_ in (ta, imm):
            for _, g in T_[T_.compartment.isin(IMMUNE)].groupby(["stratum", "compartment"]):
                T_.loc[g.index, "donor"] = g.donor.to_numpy()[rng.permutation(len(g))]
    tb = pd.concat([stroma, imm], ignore_index=True)
    tb17 = tb[tb.donor.str.split("|").str[0].isin(src17)].reset_index(drop=True)
    log(f"A404 tables: A387 {ta.donor.nunique()} patients, pooled {tb.donor.nunique()}, 17 sources {tb17.donor.nunique()}")
    stage_a404(ta, anch_a, tb, tb17, hosts, args, out, n_boot, n_perm)
    log(f"score done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
