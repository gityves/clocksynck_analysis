"""A102 and A102b: covariate sensitivity of the genome-wide GTEx tests.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A102, A102b). Extended covariates = D2 covariates + RIN + collection centre +
the top 10 genome-wide PCs of the tissue (log CPM residualised on the D2 covariates) that are associated neither with the host
nor with the time of death (both F-test p > 0.01). Steps (`--step`, default all):
- `leakage`   A102.1: R^2 of the donor-level host (fast deaths) on Hardy class, age bracket, sex, mean ischaemic time, mean RIN and
              centre, minus its permutation null mean (200); material if the median excess > 0.05. Also `interactions` (A102.5):
              log-amplitude change per SD of ischaemic time / RIN (Re(interaction / main)), core genes and A91 clock-following outputs.
- `a91`       A102.2: A91 (compute_internal_external) with the extended covariates; robust if H91 keeps p < 0.05 and the clock-following
              excess changes by < 25 %.
- `age`       A102.3: A101 with extended covariates and noise-debiased amplitudes (|coef|^2 minus its sampling variance); D = half log2 of
              the ratio of summed debiased amplitudes, old / young, outputs and core; pooled, Hardy-matched, fast; residual noise ratio.
- `counts`    A102.4: A101b counts (q < 0.2, amplitude > 1) with old and young matched on Hardy class x RIN tertile (10 draws).
- `targets`   A102.6: A100 wall-time retention (ventilated / fast) with debiased ventilated amplitude and extended covariates.
- `chiral_a91`, `chiral_a101`, `chiral_a100`  A102b: the A91, A101 and A100 tests on CHIRAL's preprocessing (RIN > 4, log2 CPM with prior 0.25,
              genes with mean count > 10, residuals of ischaemic-time bins, age, sex and Hardy class; entries > 5 log2 from the gene mean
              excluded per gene). A102b A101b is not ported.
Needs compute_internal_external.py, compute_drug_targets.py and compute_chiral_age.py outputs. Hosts: `clocksync`. Aggregates only.

    python compute_covariate_sensitivity.py [--step STEP ...] [--subsample N] [--jobs 2]
Outputs: a102_leakage.csv, a102_interactions.csv, a102_a91_tissues.csv, a102_a91_genes.csv.gz, a102_fast_summary.csv,
a102_age_{tissues,summary}.csv, a102_counts_{tissues,summary}.csv, a102_targets_{pairs.csv.gz,summary.csv},
a102b_a91_{tissues,summary}.csv, a102b_a91_genes.csv.gz, a102b_a101_{tissues,summary}.csv, a102b_a100_summary.csv.
Legacy: legacy/desync/validation/paired_damping/scripts/a102_covariates.py, a102b_chiral_prep.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _genes as G  # noqa: E402
import _shared as S  # noqa: E402
import compute_chiral_age as CA  # noqa: E402
import compute_internal_external as IE  # noqa: E402
from _genes import FOLDER, W  # noqa: E402

SEED = 102
CORE = CA.CORE


# ---- extended covariates -------------------------------------------------------------------------------------------

def ext_covariates(ob: pd.DataFrame) -> np.ndarray:
    """D2 covariates + RIN (z) + centre dummies + the PC columns (`pc_*`) attached to `ob`."""
    Z = S.d2_covariates(ob)
    rin = ob.rin.astype(float).fillna(ob.rin.astype(float).median()).to_numpy()
    cen = pd.get_dummies(ob.center.astype(str), drop_first=True).to_numpy(float)
    pcs = ob[[c for c in ob.columns if c.startswith("pc_")]].to_numpy(float)
    return np.c_[Z, (rin - rin.mean()) / (rin.std() + 1e-9), cen, pcs]


def assoc_p(x, M):
    D1 = np.c_[np.ones(len(x)), M]
    r1 = x - D1 @ np.linalg.lstsq(D1, x, rcond=None)[0]
    r0 = x - x.mean()
    k, dfr = M.shape[1], len(x) - D1.shape[1]
    return stats.f.sf(((r0 @ r0 - r1 @ r1) / k) / (r1 @ r1 / dfr), k, dfr)


def add_pcs(Y, ob, f=None, n_pc=10):
    """Attach the top PCs of the residualised expression not associated with f or the time of death (if known)."""
    Z = S.d2_covariates(ob)
    Z = np.c_[np.ones(len(Y)), Z[:, Z.std(0) > 0]]
    R = Y - Z @ np.linalg.lstsq(Z, Y, rcond=None)[0]
    R = (R - R.mean(0)) / (R.std(0) + 1e-9)
    U, s, _ = np.linalg.svd(R, full_matrices=False)
    t = ob.tod.to_numpy(float) % 24
    E = np.c_[np.cos(W * t), np.sin(W * t)]
    ob = ob.drop(columns=[c for c in ob.columns if c.startswith("pc_")])
    kept = 0
    for j in range(min(n_pc, U.shape[1])):
        x = U[:, j] * s[j]
        ok = assoc_p(x, E) > 0.01 if np.isfinite(t).all() else True
        if f is not None:
            ok = ok and assoc_p(x, f) > 0.01
        if ok:
            ob[f"pc_{j}"] = x
            kept += 1
    return ob, kept


def ext_cov(Y, ob, f):
    obp, _ = add_pcs(Y, ob, f)
    return ext_covariates(obp)


def expressed(tis, hardy, sub, need_tod=True):
    X, g, ob = G.load_rows(tis, hardy, sub, need_tod=need_tod)
    ok = np.median(X, 0) >= 5
    return np.log(X[:, ok] + 0.5), g[ok], ob


def ols(Y, D):
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    s2 = ((Y - D @ B) ** 2).sum(0) / max(len(Y) - D.shape[1], 1)
    return B, s2, np.linalg.pinv(D.T @ D)


def debiased_sq(B, s2, V, cols):
    return (B[cols] ** 2).sum(0) - s2 * np.trace(V[np.ix_(cols, cols)])


# ---- leakage, A91 with extended covariates, interactions ----------------------------------------------------------

def interactions(tis, fv, outputs, sub):
    Y, genes, ob = expressed(tis, (1, 2), sub)
    subj = ob.subject.to_numpy()
    m = np.isin(subj, fv.index)
    Y, ob = Y[m], ob[m]
    f = fv[["x", "y"]].loc[subj[m]].to_numpy()
    f = f - f.mean(0)
    Zc = S.d2_covariates(ob)
    Zc = Zc[:, Zc.std(0) > 0]
    out = dict(tissue=tis)
    for name, z in (("isch", ob.ischemic_min.astype(float)), ("rin", ob.rin.astype(float))):
        z = z.fillna(z.median()).to_numpy()
        z = (z - z.mean()) / (z.std() + 1e-9)
        D = np.c_[np.ones(len(Y)), Zc, z, f, f * z[:, None]]
        B, *_ = np.linalg.lstsq(D, Y, rcond=None)
        r = np.real((B[-2] + 1j * B[-1]) / (B[-4] + 1j * B[-3]))
        for cls, sel in (("core", np.isin(genes, CORE)), ("outputs", np.isin(genes, outputs))):
            out[f"{name}_{cls}_dlogamp_per_sd"] = float(np.median(r[sel])) if sel.any() else np.nan
    return out


def step_fast(out, o, sub, jobs, quick, steps):
    ts = G.ALL26
    hosts = G.hosts_for(ts, (1, 2), donors=None if sub is None else sub)
    fast = None if sub is None else sub
    calls = pd.read_csv(out / "a91_genes.csv.gz", usecols=["tissue", "gene", "cls"])
    rng = np.random.default_rng(SEED)
    full = CA.read_obs()
    if "leakage" in steps:
        leak = []
        for tis, fv in zip(ts, (hosts[t] for t in ts)):
            Y, genes, ob = expressed(tis, (1, 2), fast)
            m = np.isin(ob.subject.to_numpy(), fv.index)
            f = np.full((len(ob), 2), np.nan)
            f[m] = fv[["x", "y"]].loc[ob.subject.to_numpy()[m]].to_numpy()
            _, kept = add_pcs(Y[m], ob[m], f[m])
            don = full[full.DTHHRDY.isin([1, 2]) & full.subject.isin(fv.index)]
            don = don.assign(ischemic_time=don.ischemic_time.astype(float), rin=don.rin.astype(float))
            g = don.groupby("subject")
            dc = pd.DataFrame(dict(hardy=g.DTHHRDY.first().astype(str), age=g.age_bracket.first().astype(str),
                                   sex=g.sex.first().astype(str), isch=g.ischemic_time.mean(),
                                   rin=g.rin.mean()))
            dc["center"] = ob.drop_duplicates("subject").set_index("subject").center.reindex(dc.index).fillna("na")
            dc = dc.reindex(fv.index).dropna(subset=["hardy"])
            M = pd.get_dummies(dc[["hardy", "age", "sex", "center"]], drop_first=True).astype(float)
            M[["isch", "rin"]] = dc[["isch", "rin"]].fillna(dc[["isch", "rin"]].median())
            M = M.to_numpy(float)
            F = fv[["x", "y"]].loc[dc.index].to_numpy()

            def r2(Mx):
                D = np.c_[np.ones(len(Mx)), Mx]
                res = F - D @ np.linalg.lstsq(D, F, rcond=None)[0]
                return float(1 - (res ** 2).sum() / ((F - F.mean(0)) ** 2).sum())

            obs = r2(M)
            null = np.mean([r2(M[rng.permutation(len(M))]) for _ in range(20 if quick else 200)])
            leak.append(dict(tissue=tis, donors=len(M), pcs_kept=kept, r2=obs, r2_null=null, excess=obs - null))
        leak = pd.DataFrame(leak)
        leak.to_csv(out / "a102_leakage.csv", index=False)
        inter = pd.DataFrame(Parallel(n_jobs=jobs)(
            delayed(interactions)(t, hosts[t], calls[(calls.tissue == t) & calls.cls.isin(["clock", "both"])].gene, fast)
            for t in ts))
        inter.to_csv(out / "a102_interactions.csv", index=False)
        v = dict(leak_median_excess_r2=float(leak.excess.median()), leak_max_excess_r2=float(leak.excess.max()),
                 leak_material=bool(leak.excess.median() > 0.05))
        for cls in ("core", "outputs"):
            for z in ("isch", "rin"):
                v[f"interaction_{z}_{cls}_median"] = float(inter[f"{z}_{cls}_dlogamp_per_sd"].median())
        pd.Series(v).to_csv(out / "a102_leakage_summary.csv", header=False)
        print(pd.Series(v).to_string(), flush=True)
    if "a91" in steps:
        res = Parallel(n_jobs=jobs)(delayed(IE.a91_tissue)(t, hosts[t], fast, SEED + i, ext_cov) for i, t in enumerate(ts))
        summ = pd.DataFrame([r[0] for r in res])
        per = pd.concat([r[1] for r in res], ignore_index=True)
        summ.to_csv(out / "a102_a91_tissues.csv", index=False)
        per.to_csv(out / "a102_a91_genes.csv.gz", index=False)
        g26 = pd.read_csv(out / "selection_genes_26.csv.gz")
        h, _, _ = IE.h91(per, summ, g26, 100 if quick else 2000, rng)
        old = pd.read_csv(out / "a91_tissues.csv")
        ex = lambda s, k: float((s[f"n_{k}"] - s[f"null_{k}"]).clip(lower=0).sum())
        v = dict(H91_ext_diff=h["median_diff"], H91_ext_p=h["p_upper"], H91_ext_tissues=h["tissues"],
                 clock_excess_a91=ex(old, "clock"), clock_excess_ext=ex(summ, "clock"),
                 time_excess_a91=ex(old, "time"), time_excess_ext=ex(summ, "time"))
        v["clock_excess_change"] = v["clock_excess_ext"] / v["clock_excess_a91"] - 1
        v["A102_2_robust"] = bool(v["H91_ext_p"] < 0.05 and abs(v["clock_excess_change"]) < 0.25)
        pd.Series(v).to_csv(out / "a102_fast_summary.csv", header=False)
        print(pd.Series(v).to_string(), flush=True)


# ---- A102.3: age with debiased amplitudes ---------------------------------------------------------------------------

def tissue_D_deb(Y, Zc, f, out_cols, core_cols, rows_old, rows_young, h_old, h_young, version, n_sub, rng):
    res = []
    for _ in range(n_sub):
        io, iy = CA.draw(rows_old, rows_young, h_old, h_young, version, rng)
        if min(len(io), len(iy)) < 20:
            return (np.nan,) * 3
        vals = []
        for idx in (io, iy):
            Z = Zc[idx][:, Zc[idx].std(0) > 0]
            D = np.c_[np.ones(len(idx)), Z, f[idx] - f[idx].mean(0)]
            B, s2, V = ols(Y[idx], D)
            k = D.shape[1]
            vals.append((debiased_sq(B, s2, V, [k - 2, k - 1]), s2))
        (do, so), (dy, sy) = vals
        res.append((0.5 * np.log2(do[out_cols].sum() / dy[out_cols].sum()),
                    0.5 * np.log2(do[core_cols].sum() / dy[core_cols].sum()),
                    np.log2(np.median(so[out_cols]) / np.median(sy[out_cols]))))
    return tuple(np.nanmedian(np.array(res, float), axis=0))


def one_age(tis, fv, out_genes, age, seed, sub):
    Y, genes, ob = expressed(tis, (0, 1, 2, 3, 4), sub, need_tod=False)
    ob = ob.assign(age=age.reindex(ob.subject).to_numpy())
    rng = np.random.default_rng(seed)
    subj = ob.subject.to_numpy()
    m = np.isin(subj, fv.index) & ob.age.notna().to_numpy() & ~pd.Series(subj).duplicated().to_numpy()
    Y, ob, subj = Y[m], ob[m], subj[m]
    f = fv[["x", "y"]].loc[subj].to_numpy()
    obp, kept = add_pcs(Y, ob.assign(tod=np.nan), f)
    Zc = ext_covariates(obp)
    old = (ob.age >= 50).to_numpy()
    hardy = ob.hardy.to_numpy(int)
    out_cols, core_cols = np.isin(genes, out_genes), np.isin(genes, CORE)
    o_, y_ = np.flatnonzero(old), np.flatnonzero(~old)
    row = dict(tissue=tis, pcs_kept=kept)
    for v in CA.VERSIONS:
        d, dc, noise = tissue_D_deb(Y, Zc, f, out_cols, core_cols, o_, y_, hardy[o_], hardy[y_], v, 20, rng)
        row.update({f"{v}_D": d, f"{v}_Dcore": dc, f"{v}_log2_noise_ratio": noise})
    return row


def step_age(out, o, sub, jobs, quick):
    age = pd.to_numeric(gtex.phenotypes(["AGE"]).AGE, errors="coerce")
    calls = pd.read_csv(out / "a91_genes.csv.gz", usecols=["tissue", "gene", "cls"])
    calls = calls[calls.cls.isin(["clock", "both"])]
    ts = [t for t in G.ALL26 if (calls.tissue == t).sum() >= 50]
    hosts = G.hosts_for(ts, (0, 1, 2, 3, 4), donors=sub)
    per = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(one_age)(t, hosts[t], calls[calls.tissue == t].gene.unique(), age,
                                                              SEED + i, sub) for i, t in enumerate(ts)))
    per.to_csv(out / "a102_age_tissues.csv", index=False)
    v = {c: float(per[c].median()) for c in per.columns if c not in ("tissue", "pcs_kept")}
    v["tissues"] = len(per)
    pd.Series(v).to_csv(out / "a102_age_summary.csv", header=False)
    print(pd.Series(v).to_string(), flush=True)


# ---- A102.4: counts matched on Hardy x RIN ---------------------------------------------------------------------------

def one_counts(tis, seed, meta, sub):
    cpm, genes, ob = G.load_rows(tis, (0, 1, 2, 3, 4), sub, need_tod=False)
    ob = ob.join(meta[["CT"]])
    ok = (ob.CT.notna() & (ob.rin > 4) & ob.ischemic_min.notna()).to_numpy()
    cpm, ob = cpm[ok], ob[ok]
    if ob.age_bracket.isin(CA.OLD).sum() < 25 or ob.age_bracket.isin(CA.YOUNG).sum() < 25:
        return None
    X = cpm * ob.library_size.to_numpy(float)[:, None] / 1e6
    keep = X.mean(0) > 10
    E = CA.residualise(np.log2((X[:, keep] + 0.25) / (X.sum(1, keepdims=True) + 0.5) * 1e6), ob).astype(np.float32)
    rng = np.random.default_rng(seed)
    phi = ob.CT.to_numpy(float) * np.pi / 12
    grp = np.where(ob.age_bracket.isin(CA.OLD), "old", np.where(ob.age_bracket.isin(CA.YOUNG), "young", "mid"))
    hardy = ob.hardy.to_numpy(int)
    rin_t = pd.qcut(ob.rin.astype(float).rank(method="first"), 3, labels=False).to_numpy()
    key = {"pooled": hardy, "matched": hardy, "matched_rin": hardy * 10 + rin_t}
    ro, ry = np.flatnonzero(grp == "old"), np.flatnonzero(grp == "young")
    row = dict(tissue=tis)
    for v, k in key.items():
        acc = []
        for _ in range(10):
            io, iy = CA.draw_b(ro, ry, k[ro], k[ry], "pooled" if v == "pooled" else "matched", rng)
            if min(len(io), len(iy)) < 25:
                break
            co, _ = CA.counts(E, phi, io)
            cy, _ = CA.counts(E, phi, iy)
            so = np.median(((E[io] - E[io].mean(0)) ** 2).mean(0))
            sy = np.median(((E[iy] - E[iy].mean(0)) ** 2).mean(0))
            acc.append((co[1], cy[1], np.log2(so / sy), len(io)))
        if acc:
            r = np.median(np.array(acc, float), axis=0)
            row.update({f"{v}_old": r[0], f"{v}_young": r[1], f"{v}_log2_var_ratio": r[2], f"{v}_n": r[3]})
    row["rin_old_mean"], row["rin_young_mean"] = float(ob.rin[grp == "old"].mean()), float(ob.rin[grp == "young"].mean())
    return row


def step_counts(out, o, sub, jobs):
    meta = pd.read_csv(CA.META, sep=" ")
    meta["sample"] = meta.fullID.str.replace(".", "-", regex=False)
    meta = meta.set_index("sample")[["CT"]]
    tissues = sorted(t for t in o.tissue.unique() if not t.startswith(("Cells", "Whole Blood")))
    per = pd.DataFrame([r for r in Parallel(n_jobs=jobs)(delayed(one_counts)(t, SEED + i, meta, sub)
                                                         for i, t in enumerate(tissues)) if r is not None])
    per.to_csv(out / "a102_counts_tissues.csv", index=False)
    v = {}
    for ver in ("pooled", "matched", "matched_rin"):
        x = per.dropna(subset=[f"{ver}_old"])
        x = x[(x[f"{ver}_old"] + x[f"{ver}_young"]) >= 20]
        v[f"{ver}_tissues"] = len(x)
        v[f"{ver}_median_tissue_logratio"] = float(np.median(np.log((x[f"{ver}_old"] + 1) / (x[f"{ver}_young"] + 1))))
        v[f"{ver}_total_logratio"] = float(np.log(x[f"{ver}_old"].sum() / x[f"{ver}_young"].sum()))
        v[f"{ver}_median_log2_var_ratio"] = float(per[f"{ver}_log2_var_ratio"].median())
    v["rin_old_minus_young_median"] = float((per.rin_old_mean - per.rin_young_mean).median())
    pd.Series(v).to_csv(out / "a102_counts_summary.csv", header=False)
    print(pd.Series(v).to_string(), flush=True)


# ---- A102.6: retention with debiased amplitudes -----------------------------------------------------------------------

def one_ret(tis, genes, sub):
    amp = {}
    for g, hardy in (("vent", (0,)), ("fast", (1, 2))):
        Y, names, ob = expressed(tis, hardy, sub)
        t = ob.tod.to_numpy(float) % 24
        obp, _ = add_pcs(Y, ob)
        Zc = ext_covariates(obp)
        Zc = Zc[:, Zc.std(0) > 0]
        D = np.c_[np.ones(len(Y)), Zc, np.cos(W * t), np.sin(W * t)]
        idx = pd.Index(names).get_indexer(genes)
        B, s2, V = ols(Y[:, idx], D)
        k = D.shape[1]
        amp[g] = (np.sqrt(np.maximum(debiased_sq(B, s2, V, [k - 2, k - 1]), 0)), names[idx[idx >= 0]], idx)
    keep = (amp["vent"][2] >= 0) & (amp["fast"][2] >= 0)
    return pd.DataFrame(dict(tissue=tis, gene=np.asarray(genes)[keep],
                             ret_vent_deb=amp["vent"][0][keep] / amp["fast"][0][keep]))


def step_targets(out, sub, jobs):
    pairs = pd.read_csv(out / "a100_rhythmic_pairs.csv.gz")
    res = pd.concat(Parallel(n_jobs=jobs)(delayed(one_ret)(t, x.gene.unique(), sub) for t, x in pairs.groupby("tissue")),
                    ignore_index=True)
    r = pairs.merge(res, on=["tissue", "gene"], how="left")
    v = {}
    for who, x in (("targets", r[r.target]), ("others", r[~r.target])):
        v[f"{who}_ret_vent_a100"] = float(x.ret_vent.median())
        v[f"{who}_ret_vent_debiased_ext"] = float(x.ret_vent_deb.median())
        v[f"{who}_frac_zero_debiased"] = float((x.ret_vent_deb == 0).mean())
        for c in ("locked", "intermediate", "cue-like"):
            v[f"{who}_{c}_ret_vent_debiased_ext"] = float(x[x.lock_class == c].ret_vent_deb.median())
    r.to_csv(out / "a102_targets_pairs.csv.gz", index=False)
    pd.Series(v).to_csv(out / "a102_targets_summary.csv", header=False)
    print(pd.Series(v).to_string(), flush=True)


# ---- A102b: CHIRAL preprocessing ---------------------------------------------------------------------------------------

def prep(tis, sub):
    """CHIRAL residuals (n x genes), obs and gene names of the Hardy 0-4 samples of one tissue."""
    cpm, genes, ob = G.load_rows(tis, (0, 1, 2, 3, 4), sub, need_tod=False)
    ok = (ob.rin > 4).to_numpy() & ob.ischemic_min.notna().to_numpy()
    cpm, ob = cpm[ok], ob[ok]
    X = cpm * ob.library_size.to_numpy(float)[:, None] / 1e6
    sel = X.mean(0) > 10
    L = np.log2((X[:, sel] + 0.25) / (X.sum(1, keepdims=True) + 0.5) * 1e6)
    return CA.residualise(L, ob).astype(np.float32), ob, genes[sel]


def mask5(E):
    E = E.astype(np.float64).copy()
    E[np.abs(E - np.nanmean(E, 0)) > 5] = np.nan
    return E


def ols_nan(Y, D, cols):
    """Per gene OLS on the finite rows: coefficients, RSS, residual df, trace of (D'D)^-1 over `cols`."""
    g = Y.shape[1]
    Bm = np.full((D.shape[1], g), np.nan)
    rss, df, tr = np.full(g, np.nan), np.full(g, np.nan), np.full(g, np.nan)
    fin = np.isfinite(Y)
    full = fin.all(0)
    if full.any():
        Bf, *_ = np.linalg.lstsq(D, Y[:, full], rcond=None)
        Bm[:, full] = Bf
        rss[full] = ((Y[:, full] - D @ Bf) ** 2).sum(0)
        df[full] = len(Y) - D.shape[1]
        tr[full] = np.trace(np.linalg.pinv(D.T @ D)[np.ix_(cols, cols)])
    # genes sharing one pattern of excluded entries are solved together (the patterns are few: single outliers)
    bad = np.flatnonzero(~full)
    if len(bad):
        keys = np.packbits(fin[:, bad], axis=0)
        _, inv = np.unique(keys, axis=1, return_inverse=True)
        inv = np.asarray(inv).ravel()
        for u in np.unique(inv):
            cols_u = bad[inv == u]
            m = fin[:, cols_u[0]]
            Dj = D[m]
            if m.sum() <= D.shape[1] + 2:
                continue
            bj, *_ = np.linalg.lstsq(Dj, Y[np.ix_(m, cols_u)], rcond=None)
            Bm[:, cols_u] = bj
            rss[cols_u] = ((Y[np.ix_(m, cols_u)] - Dj @ bj) ** 2).sum(0)
            df[cols_u] = m.sum() - D.shape[1]
            tr[cols_u] = np.trace(np.linalg.pinv(Dj.T @ Dj)[np.ix_(cols, cols)])
    return Bm, rss, df, tr


def partial_p(Y, D, idx):
    _, rss1, df1, _ = ols_nan(Y, D, [0])
    _, rss0, _, _ = ols_nan(Y, np.delete(D, idx, axis=1), [0])
    return stats.f.sf(((rss0 - rss1) / len(idx)) / (rss1 / df1), len(idx), df1)


def chiral_calls(Y, f, E):
    D = np.c_[np.ones(len(Y)), f, E]
    p_int, p_ext = partial_p(Y, D, [1, 2]), partial_p(Y, D, [3, 4])
    cls = np.full(Y.shape[1], "none", dtype=object)
    cls[(p_int < 0.01) & (p_ext > 0.05)] = "clock"
    cls[(p_ext < 0.01) & (p_int > 0.05)] = "time"
    cls[(p_int < 0.01) & (p_ext < 0.01)] = "both"
    return cls


def chiral_a91_tissue(tis, fv, seed, sub):
    R, ob, genes = prep(tis, sub)
    rng = np.random.default_rng(seed)
    subj = ob.subject.to_numpy()
    m = ob.hardy.isin([1, 2]).to_numpy() & ob.tod.notna().to_numpy() & np.isin(subj, fv.index)
    Y, ob = mask5(R[m]), ob[m]
    f = fv[["x", "y"]].loc[subj[m]].to_numpy()
    f = f - f.mean(0)
    t = ob.tod.to_numpy(float) % 24
    E = np.c_[np.cos(W * t), np.sin(W * t)]
    nonclock = ~np.isin(genes, list(G.CLOCK))
    count = lambda c: {k: int(((c == k) & nonclock).sum()) for k in ("clock", "time", "both")}
    cls = chiral_calls(Y, f, E)
    nf = pd.DataFrame([count(chiral_calls(Y, f[rng.permutation(len(f))], E)) for _ in range(20)]).median()
    nt = pd.DataFrame([count(chiral_calls(Y, f, E[rng.permutation(len(E))])) for _ in range(20)]).median()
    row = dict(tissue=tis, donors=len(Y), **{f"n_{k}": v for k, v in count(cls).items()},
               null_clock=float(nf["clock"]), null_time=float(nt["time"]), null_both=float(max(nf["both"], nt["both"])),
               outlier_frac=float(np.isnan(Y).mean()))
    return row, pd.DataFrame(dict(tissue=tis, gene=genes, cls=cls))


def step_chiral_a91(out, sub, jobs, quick, rng):
    hosts = G.hosts_for(G.ALL26, (1, 2), donors=sub)
    res = Parallel(n_jobs=jobs)(delayed(chiral_a91_tissue)(t, hosts[t], 1022 + i, sub) for i, t in enumerate(G.ALL26))
    summ = pd.DataFrame([r[0] for r in res])
    per = pd.concat([r[1] for r in res], ignore_index=True)
    summ.to_csv(out / "a102b_a91_tissues.csv", index=False)
    per.to_csv(out / "a102b_a91_genes.csv.gz", index=False)
    g26 = pd.read_csv(out / "selection_genes_26.csv.gz")
    h, _, _ = IE.h91(per, summ, g26, 100 if quick else 2000, np.random.default_rng(1022))
    old = pd.read_csv(out / "a91_tissues.csv")
    ex = lambda s, k: float((s[f"n_{k}"] - s[f"null_{k}"]).clip(lower=0).sum())
    v = dict(H91_diff=h["median_diff"], H91_p=h["p_upper"], H91_tissues=h["tissues"],
             clock_excess_a91=ex(old, "clock"), clock_excess_chiral=ex(summ, "clock"),
             time_excess_a91=ex(old, "time"), time_excess_chiral=ex(summ, "time"),
             median_outlier_frac=float(summ.outlier_frac.median()))
    v["clock_excess_change"] = v["clock_excess_chiral"] / v["clock_excess_a91"] - 1
    v["robust"] = bool(v["H91_p"] < 0.05 and abs(v["clock_excess_change"]) < 0.25)
    pd.Series(v).to_csv(out / "a102b_a91_summary.csv", header=False)
    print(pd.Series(v).to_string(), flush=True)


def chiral_a101_tissue(tis, fv, out_genes, age, seed, sub):
    R, ob, genes = prep(tis, sub)
    rng = np.random.default_rng(seed)
    ob = ob.assign(age=age.reindex(ob.subject).to_numpy())
    subj = ob.subject.to_numpy()
    m = np.isin(subj, fv.index) & ob.age.notna().to_numpy() & ~pd.Series(subj).duplicated().to_numpy()
    R, ob, subj = R[m], ob[m], subj[m]
    f = fv[["x", "y"]].loc[subj].to_numpy()
    old = (ob.age >= 50).to_numpy()
    hardy = ob.hardy.to_numpy(int)
    oc, cc = np.isin(genes, out_genes), np.isin(genes, CORE)
    o_, y_ = np.flatnonzero(old), np.flatnonzero(~old)
    row = dict(tissue=tis)
    for v in CA.VERSIONS:
        acc = []
        for _ in range(20):
            io, iy = CA.draw(o_, y_, hardy[o_], hardy[y_], v, rng)
            if min(len(io), len(iy)) < 20:
                break
            amp, deb = [], []
            for idx in (io, iy):
                Y = mask5(R[idx][:, oc | cc])
                D = np.c_[np.ones(len(idx)), f[idx] - f[idx].mean(0)]
                Bm, rss, df, tr = ols_nan(Y, D, [1, 2])
                a2 = Bm[1] ** 2 + Bm[2] ** 2
                amp.append(np.sqrt(a2))
                deb.append(a2 - rss / df * tr)
            so_, sc_ = oc[oc | cc], cc[oc | cc]
            r = np.log2(amp[0] / amp[1])
            acc.append((np.nanmedian(r[so_]), np.nanmedian(r[sc_]),
                        0.5 * np.log2(np.nansum(deb[0][so_]) / np.nansum(deb[1][so_])),
                        0.5 * np.log2(np.nansum(deb[0][sc_]) / np.nansum(deb[1][sc_]))))
        if acc:
            r = np.nanmedian(np.array(acc, float), axis=0)
            row.update({f"{v}_D": r[0], f"{v}_Dcore": r[1], f"{v}_Ddeb": r[2], f"{v}_Dcore_deb": r[3]})
    return row


def step_chiral_a101(out, sub, jobs):
    age = pd.to_numeric(gtex.phenotypes(["AGE"]).AGE, errors="coerce")
    calls = pd.read_csv(out / "a91_genes.csv.gz", usecols=["tissue", "gene", "cls"])
    calls = calls[calls.cls.isin(["clock", "both"])]
    ts = [t for t in G.ALL26 if (calls.tissue == t).sum() >= 50]
    hosts = G.hosts_for(ts, (0, 1, 2, 3, 4), donors=sub)
    per = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(chiral_a101_tissue)(t, hosts[t], calls[calls.tissue == t].gene.unique(),
                                                                         age, 1022 + i, sub) for i, t in enumerate(ts)))
    per.to_csv(out / "a102b_a101_tissues.csv", index=False)
    v = {c: float(per[c].median()) for c in per.columns if c != "tissue"}
    v["tissues"] = len(per)
    pd.Series(v).to_csv(out / "a102b_a101_summary.csv", header=False)
    print(pd.Series(v).to_string(), flush=True)


def chiral_a100_tissue(tis, pair_genes, sub):
    R, ob, genes = prep(tis, sub)
    idx = pd.Index(genes).get_indexer(pair_genes)
    pair_genes, idx = pair_genes[idx >= 0], idx[idx >= 0]
    amp = {}
    for g, hardy in (("vent", [0]), ("fast", [1, 2])):
        m = ob.hardy.isin(hardy).to_numpy() & ob.tod.notna().to_numpy()
        t = ob.tod.to_numpy(float)[m] % 24
        Y = mask5(R[m][:, idx])
        D = np.c_[np.ones(m.sum()), np.cos(W * t), np.sin(W * t)]
        Bm, rss, df, tr = ols_nan(Y, D, [1, 2])
        amp[g] = np.sqrt(np.maximum(Bm[1] ** 2 + Bm[2] ** 2 - rss / df * tr, 0))
    return pd.DataFrame(dict(tissue=tis, gene=pair_genes, ret_vent_chiral=amp["vent"] / amp["fast"]))


def step_chiral_a100(out, sub, jobs):
    pairs = pd.read_csv(out / "a100_rhythmic_pairs.csv.gz")
    res = pd.concat(Parallel(n_jobs=jobs)(delayed(chiral_a100_tissue)(t, x.gene.unique(), sub)
                                          for t, x in pairs.groupby("tissue")), ignore_index=True)
    r = pairs.merge(res, on=["tissue", "gene"], how="left")
    v = {}
    for who, x in (("targets", r[r.target]), ("others", r[~r.target])):
        v[f"{who}_ret_vent_chiral"] = float(x.ret_vent_chiral.median())
        v[f"{who}_frac_zero"] = float((x.ret_vent_chiral == 0).mean())
        v[f"{who}_genes_measured"] = int(x.ret_vent_chiral.notna().sum())
        for c in ("locked", "intermediate", "cue-like"):
            v[f"{who}_{c}_ret_vent_chiral"] = float(x[x.lock_class == c].ret_vent_chiral.median())
    pd.Series(v).to_csv(out / "a102b_a100_summary.csv", header=False)
    print(pd.Series(v).to_string(), flush=True)


STEPS = ("leakage", "a91", "age", "counts", "targets", "chiral_a91", "chiral_a101", "chiral_a100")


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--step", nargs="+", default=list(STEPS), choices=STEPS)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    o = CA.read_obs()
    sub = G.pick_donors(o, args.subsample, args.seed)
    st = set(args.step)
    if st & {"leakage", "a91"}:
        step_fast(out, o, sub, args.jobs, run["quick"], st)
    if "age" in st:
        step_age(out, o, sub, args.jobs, run["quick"])
    if "counts" in st:
        step_counts(out, o, sub, args.jobs)
    if "targets" in st:
        step_targets(out, sub, args.jobs)
    if "chiral_a91" in st:
        step_chiral_a91(out, sub, args.jobs, run["quick"], None)
    if "chiral_a101" in st:
        step_chiral_a101(out, sub, args.jobs)
    if "chiral_a100" in st:
        step_chiral_a100(out, sub, args.jobs)


if __name__ == "__main__":
    main()
