"""A101 and A101b: does GTEx show CHIRAL's age dampening of rhythmic outputs, and is it the death-type composition of age groups?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A101, A101b).
A101 (all donors, Hardy 0-4; exact age from the protected phenotypes, in memory): internal phase = host from the clock genes
of the other `v2_gtex` tissues, fitted on all donors (`clocksync` `Dataset.fit(...).host`). Per tissue (A91 tissues with >= 50
clock-following or both calls), amplitude |coefficient| of log CPM on the host (2 columns, D2 covariates) of each output gene
(A91 clock-following or both) and of the 14 core genes, old (>= 50) vs young (< 50) at equal n, three versions: (a) pooled over
death classes (CHIRAL-like), (b) equal n within each Hardy class (composition-matched), (c) fast deaths only (20 subsamples).
D = median over tissues of the median log2(old / young) over output genes; D_core the same for core genes; D_rel = D - D_core.
Donor bootstrap (500, within age x Hardy cells, one draw shared by all tissues, 4 subsamples per replicate).
H101.1: upper 95 % bound of D_pooled < 0. H101.2: upper bound of D_pooled - D_matched < 0.
A101b: CHIRAL's own age measure on CHIRAL's published per-donor phases (naef-lab/CHIRAL, commit 1f78a9a, sample_metadata.txt
under $SCCC_DATA/chiral_repo/raw, read as text only; the phases are per-donor controlled data and stay in memory): per tissue
log2 CPM (prior 0.25, genes with mean count > 10, RIN > 4), residuals of ischaemic-time bins, age, sex and Hardy class,
harmonic regression on the phase per group (|deviation| > 5 set to the gene mean), amplitude = log2 peak-to-trough, LRT p,
BH q; counts of rhythmic genes (q < 0.2, amplitude > 0.5 / 1 / 2). Young < 50 y, old 60-79 y. R1: CHIRAL's age subsample (ageSS);
pooled / Hardy-matched / fast / ventilated versions at equal n (10 draws); H101b: log(N_old / N_young) matched minus pooled,
donor bootstrap (200, within age x Hardy cells), PASS if the interval excludes 0 and the difference is > 0.
Needs compute_internal_external.py (A91 classes). Aggregates only.

    python compute_chiral_age.py [--subsample N] [--jobs 2] [--boot B]
Outputs: a101_tissues.csv, a101_summary.csv, a101b_tissues.csv, a101b_summary.csv.
Legacy: legacy/desync/validation/paired_damping/scripts/a101_chiral_age.py, a101b_chiral_measure.py.
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

from common import cli, gtex  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

import _genes as G  # noqa: E402
import _shared as S  # noqa: E402
from _genes import FOLDER  # noqa: E402

SEED, N_SUB, N_BOOT, N_SUB_BOOT = 101, 20, 500, 4
VERSIONS = ("pooled", "matched", "fast")
SEED_B, N_DRAW, N_BOOT_B = 1012, 10, 200
META = DATA / "chiral_repo" / "raw" / "sample_metadata.txt"
THR = (0.5, 1.0, 2.0)
YOUNG, OLD = ("20-29", "30-39", "40-49"), ("60-69", "70-79")
CORE = list(G.CORE_GENES)  # BMAL1 is not a symbol of the GTEx table (ARNTL is)


def read_obs() -> pd.DataFrame:
    """All sample annotations used by both analyses (local: not in clocksync)."""
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs[["DTHHRDY", "tissue", "tod", "subject", "age_bracket", "ischemic_time", "rin", "sex"]].copy()
    finally:
        a.file.close()
    o["subject"], o["tissue"] = o.subject.astype(str), o.tissue.astype(str)
    return o


# ---- A101 ----------------------------------------------------------------------------------------------------------

def amps(Y, Zc, f):
    Zc = Zc[:, Zc.std(0) > 0]
    D = np.c_[np.ones(len(Y)), Zc, f - f.mean(0)]
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    return np.hypot(B[-2], B[-1])


def draw(pool_old, pool_young, hardy_old, hardy_young, version, rng):
    """Indices (into the tissue's rows) of an equal-n old and young sample under `version`."""
    if version == "fast":
        pool_old, pool_young = pool_old[np.isin(hardy_old, [1, 2])], pool_young[np.isin(hardy_young, [1, 2])]
    if version in ("pooled", "fast"):
        n = min(len(pool_old), len(pool_young))
        return rng.choice(pool_old, n, replace=False), rng.choice(pool_young, n, replace=False)
    o, y = [], []
    for h in np.intersect1d(hardy_old, hardy_young):
        po, py = pool_old[hardy_old == h], pool_young[hardy_young == h]
        n = min(len(po), len(py))
        o.append(rng.choice(po, n, replace=False))
        y.append(rng.choice(py, n, replace=False))
    return np.concatenate(o), np.concatenate(y)


def tissue_D(Y, Zc, f, out_cols, core_cols, rows_old, rows_young, h_old, h_young, version, n_sub, rng):
    res = []
    for _ in range(n_sub):
        io, iy = draw(rows_old, rows_young, h_old, h_young, version, rng)
        if min(len(io), len(iy)) < 20:
            return np.nan, np.nan
        r = np.log2(amps(Y[io], Zc[io], f[io]) / amps(Y[iy], Zc[iy], f[iy]))
        res.append((np.median(r[out_cols]), np.median(r[core_cols])))
    return tuple(np.median(np.array(res), axis=0))


def a101_tissue(tis, fv, out_genes, boots, age, seed, sub):
    # as the legacy run: the 14 core genes are appended even when A91 already lists them as outputs (they then count twice)
    genes = list(out_genes) + list(CORE)
    X, gsel, ob = G.load_rows(tis, (0, 1, 2, 3, 4), sub, need_tod=False, genes=genes)
    Y = np.log(X + 0.5)
    rng = np.random.default_rng(seed)
    subj = ob.subject.to_numpy()
    ob = ob.assign(age=age.reindex(subj).to_numpy())
    m = np.isin(subj, fv.index) & ob.age.notna().to_numpy() & ~pd.Series(subj).duplicated().to_numpy()
    Y, ob, subj = Y[m], ob[m], subj[m]
    f = fv[["x", "y"]].loc[subj].to_numpy()
    Zc = S.d2_covariates(ob)
    old = (ob.age >= 50).to_numpy()
    hardy = ob.hardy.to_numpy(int)
    out_cols = np.isin(gsel, out_genes)
    core_cols = np.isin(gsel, CORE)
    row_of = pd.Series(np.arange(len(subj)), index=subj)

    def run(rows, n_sub):
        rows = np.asarray(rows)
        o_, y_ = rows[old[rows]], rows[~old[rows]]
        return {v: tissue_D(Y, Zc, f, out_cols, core_cols, o_, y_, hardy[o_], hardy[y_], v, n_sub, rng)
                for v in VERSIONS}

    obs = run(np.arange(len(subj)), N_SUB)
    bs = [run(row_of.reindex(d).dropna().astype(int).to_numpy(), N_SUB_BOOT) for d in boots]
    comp = {f"old_frac_hardy{h}": float(np.mean(hardy[old] == h)) for h in range(5)}
    comp.update({f"young_frac_hardy{h}": float(np.mean(hardy[~old] == h)) for h in range(5)})
    return tis, obs, bs, dict(n_old=int(old.sum()), n_young=int((~old).sum()), outputs=int(out_cols.sum()), **comp)


def a101(out, o, sub, jobs, n_boot):
    age = gtex.phenotypes(["AGE"]).AGE
    age = pd.to_numeric(age, errors="coerce")
    calls = pd.read_csv(out / "a91_genes.csv.gz", usecols=["tissue", "gene", "cls"])
    calls = calls[calls.cls.isin(["clock", "both"])]
    ts = [t for t in G.ALL26 if (calls.tissue == t).sum() >= 50]
    hosts = G.hosts_for(ts, (0, 1, 2, 3, 4), donors=sub)
    don = o.drop_duplicates("subject").assign(age=lambda d: age.reindex(d.subject).to_numpy()).dropna(subset=["age"])
    don = don[don.DTHHRDY.isin([0, 1, 2, 3, 4])]
    if sub is not None:
        don = don[don.subject.isin(sub)]
    cell = (don.age >= 50).astype(str) + "|" + don.DTHHRDY.astype(int).astype(str)
    rng = np.random.default_rng(SEED)
    boots = [np.concatenate([rng.choice(g.subject.to_numpy(), len(g)) for _, g in don.groupby(cell.to_numpy())])
             for _ in range(n_boot)]
    res = Parallel(n_jobs=jobs)(delayed(a101_tissue)(t, hosts[t], calls[calls.tissue == t].gene.unique(), boots, age,
                                                     SEED + i, sub) for i, t in enumerate(ts))
    per = pd.DataFrame([dict(tissue=t, **{f"{k}_{q}": ob_[k][j] for k in VERSIONS for j, q in enumerate(("D", "Dcore"))},
                             **comp) for t, ob_, _, comp in res])
    per.to_csv(out / "a101_tissues.csv", index=False)

    def agg(table):
        o_ = {}
        for k in VERSIONS:
            d, c = np.array([t[k][0] for t in table]), np.array([t[k][1] for t in table])
            o_[k] = (np.nanmedian(d), np.nanmedian(c), np.nanmedian(d - c))
        return o_

    obs = agg([r[1] for r in res])
    boot = [agg([r[2][b] for r in res]) for b in range(n_boot)]
    v = {}
    for k in VERSIONS:
        for j, q in enumerate(("D", "Dcore", "Drel")):
            x = np.array([b[k][j] for b in boot])
            lo_hi = np.nanquantile(x, [0.025, 0.975]) if n_boot else (np.nan, np.nan)
            v[f"{q}_{k}"], v[f"{q}_{k}_lo"], v[f"{q}_{k}_hi"] = obs[k][j], *lo_hi
    diff = np.array([b["pooled"][0] - b["matched"][0] for b in boot])
    v["pooled_minus_matched"] = obs["pooled"][0] - obs["matched"][0]
    v["pooled_minus_matched_lo"], v["pooled_minus_matched_hi"] = np.nanquantile(diff, [0.025, 0.975]) if n_boot else (np.nan,) * 2
    v["H101_1"] = "PASS" if v["D_pooled_hi"] < 0 else "FAIL"
    v["H101_2"] = "PASS" if v["pooled_minus_matched_hi"] < 0 else "FAIL"
    v["tissues"], v["n_boot"] = len(per), n_boot
    pd.Series(v).to_csv(out / "a101_summary.csv", header=False)
    print(pd.Series(v).to_string(), flush=True)


# ---- A101b ---------------------------------------------------------------------------------------------------------

def residualise(L, ob):
    """Residuals of ischaemic-time bins (0-200-500-800-2100 min), age bracket midpoint, sex and Hardy class (CHIRAL)."""
    isch = pd.cut(ob.ischemic_min.clip(lower=1), [0, 200, 500, 800, 2100]).astype(str)
    age = ob.age_bracket.astype(str).str[:2].astype(float) + 5
    D = pd.get_dummies(pd.DataFrame(dict(isch=isch, sex=ob.sex.astype(str), hardy=ob.hardy.astype(str))),
                       drop_first=True).astype(float)
    D = np.c_[np.ones(len(L)), age.to_numpy(), D.to_numpy()]
    B, *_ = np.linalg.lstsq(D, L, rcond=None)
    return L - D @ B


def bh(p):
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty_like(p)
    q[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return q


def harmonic(E, phi):
    """Per gene: log2 peak-to-trough amplitude and BH q of the harmonic regression on phi (radians); |dev| > 5 set to the mean."""
    E = E - E.mean(0)
    E = np.where(np.abs(E) > 5, 0.0, E)
    n = len(E)
    X = np.c_[np.ones(n), np.cos(phi), np.sin(phi)]
    B, *_ = np.linalg.lstsq(X, E, rcond=None)
    rss1 = ((E - X @ B) ** 2).sum(0)
    rss0 = ((E - E.mean(0)) ** 2).sum(0)
    p = stats.chi2.sf(n * np.log(rss0 / rss1), 2)
    return 2 * np.hypot(B[1], B[2]), bh(p)


def counts(E, phi, idx):
    amp, q = harmonic(E[idx], phi[idx])
    return np.array([np.sum((q < 0.2) & (amp > t)) for t in THR]), amp


def draw_b(rows_o, rows_y, h_o, h_y, version, rng):
    if version == "matched":
        o, y = [], []
        for h in np.intersect1d(h_o, h_y):
            po, py = rows_o[h_o == h], rows_y[h_y == h]
            n = min(len(po), len(py))
            o.append(rng.choice(po, n, replace=False))
            y.append(rng.choice(py, n, replace=False))
        return np.concatenate(o), np.concatenate(y)
    if version == "fast":
        rows_o, rows_y = rows_o[np.isin(h_o, [1, 2])], rows_y[np.isin(h_y, [1, 2])]
    if version == "vent":
        rows_o, rows_y = rows_o[h_o == 0], rows_y[h_y == 0]
    n = min(len(rows_o), len(rows_y))
    return rng.choice(rows_o, n, replace=False), rng.choice(rows_y, n, replace=False)


def a101b_tissue(tis, boots, seed, meta, sub):
    cpm, genes, ob = G.load_rows(tis, (0, 1, 2, 3, 4), sub, need_tod=False)
    ob = ob.join(meta[["CT", "ageSS"]])
    ok = (ob.CT.notna() & (ob.rin > 4) & ob.ischemic_min.notna()).to_numpy()
    cpm, ob = cpm[ok], ob[ok]
    if ob.age_bracket.isin(OLD).sum() < 25 or ob.age_bracket.isin(YOUNG).sum() < 25:
        return None
    X = cpm * ob.library_size.to_numpy(float)[:, None] / 1e6
    keep = X.mean(0) > 10
    lib = X.sum(1, keepdims=True)
    E = residualise(np.log2((X[:, keep] + 0.25) / (lib + 0.5) * 1e6), ob).astype(np.float32)
    genes = genes[keep]
    del X, cpm
    rng = np.random.default_rng(seed)
    phi = ob.CT.to_numpy(float) * np.pi / 12
    grp = np.where(ob.age_bracket.isin(OLD), "old", np.where(ob.age_bracket.isin(YOUNG), "young", "mid"))
    hardy = ob.hardy.to_numpy(int)
    core = np.isin(genes, CORE)
    subj = ob.subject.to_numpy()
    out = dict(tissue=tis, genes=len(genes), n_old=int((grp == "old").sum()), n_young=int((grp == "young").sum()))
    for g in ("young", "old"):
        idx = np.flatnonzero((grp == g) & (ob.ageSS.to_numpy() == 1))
        if len(idx) < 25:
            return None
        c, amp = counts(E, phi, idx)
        out.update({f"R1_{g}_{t}": int(x) for t, x in zip(THR, c)})
        out[f"R1_{g}_n"], out[f"R1_{g}_core_amp"] = len(idx), float(np.median(amp[core]))
    rows_o, rows_y = np.flatnonzero(grp == "old"), np.flatnonzero(grp == "young")

    def version_counts(ro, ry, version, k):
        acc = []
        for _ in range(k):
            io, iy = draw_b(ro, ry, hardy[ro], hardy[ry], version, rng)
            if min(len(io), len(iy)) < 25:
                return None
            co, ao = counts(E, phi, io)
            cy, ay = counts(E, phi, iy)
            acc.append(np.r_[co, cy, np.median(ao[core]), np.median(ay[core]), len(io)])
        return np.median(np.array(acc), axis=0)

    for v in ("pooled", "matched", "fast", "vent"):
        r = version_counts(rows_o, rows_y, v, N_DRAW)
        if r is None:
            continue
        for j, t in enumerate(THR):
            out[f"{v}_old_{t}"], out[f"{v}_young_{t}"] = r[j], r[len(THR) + j]
        out[f"{v}_core_amp_old"], out[f"{v}_core_amp_young"], out[f"{v}_n"] = r[-3], r[-2], r[-1]
    row_of = pd.Series(np.arange(len(subj)), index=subj)
    bs = []
    for donors in boots:
        rows = row_of.reindex(donors).dropna().astype(int).to_numpy()
        ro, ry = rows[grp[rows] == "old"], rows[grp[rows] == "young"]
        bs.append({v: version_counts(ro, ry, v, 1) for v in ("pooled", "matched")})
    return out, bs


def a101b(out, o, sub, jobs, n_boot):
    meta = pd.read_csv(META, sep=" ")
    meta["sample"] = meta.fullID.str.replace(".", "-", regex=False)
    meta = meta.set_index("sample")[["CT", "ageSS"]]
    meta["ageSS"] = meta.ageSS.fillna(0)
    oo = o.join(meta)
    ok = oo.CT.notna() & oo.DTHHRDY.isin([0, 1, 2, 3, 4]) & (oo.rin > 4) & oo.ischemic_time.notna()
    ok &= ~oo.tissue.str.startswith(("Cells", "Whole Blood"))
    oo = oo[ok]
    if sub is not None:
        oo = oo[oo.subject.isin(sub)]
    don = oo.drop_duplicates("subject")
    don = don[don.age_bracket.astype(str).isin(YOUNG + OLD)]
    cell = don.age_bracket.astype(str).isin(OLD).astype(str) + "|" + don.DTHHRDY.astype(int).astype(str)
    rng = np.random.default_rng(SEED_B)
    boots = [np.concatenate([rng.choice(g.subject.to_numpy(), len(g)) for _, g in don.groupby(cell.to_numpy())])
             for _ in range(n_boot)]
    tissues = sorted(oo.tissue.unique())
    res = [(t, r) for t, r in zip(tissues, Parallel(n_jobs=jobs)(
        delayed(a101b_tissue)(t, boots, SEED_B + i, meta, sub) for i, t in enumerate(tissues))) if r is not None]
    per = pd.DataFrame([r[0] for _, r in res])
    per.to_csv(out / "a101b_tissues.csv", index=False)
    v = dict(tissues=len(per))
    lr = lambda o_, y_: float(np.log(o_.sum() / y_.sum()))
    for t in THR:
        v[f"R1_old_{t}"], v[f"R1_young_{t}"] = int(per[f"R1_old_{t}"].sum()), int(per[f"R1_young_{t}"].sum())
        v[f"R1_logratio_{t}"] = lr(per[f"R1_old_{t}"], per[f"R1_young_{t}"])
        for ver in ("pooled", "matched", "fast", "vent"):
            if f"{ver}_old_{t}" in per:
                x = per.dropna(subset=[f"{ver}_old_{t}"])
                v[f"{ver}_logratio_{t}"] = lr(x[f"{ver}_old_{t}"], x[f"{ver}_young_{t}"])
                v[f"{ver}_tissues"] = len(x)
    v["R1_core_amp_log2ratio"] = float(np.median(np.log2(per.R1_old_core_amp / per.R1_young_core_amp)))
    v["pooled_core_amp_log2ratio"] = float(np.median(np.log2(per.pooled_core_amp_old / per.pooled_core_amp_young)))
    j = THR.index(1.0)
    diffs = []
    for b in range(n_boot):
        tot = {}
        for ver in ("pooled", "matched"):
            rr = [r[1][b][ver] for _, r in res if r[1][b][ver] is not None]
            o_ = sum(x[j] for x in rr)
            y_ = sum(x[len(THR) + j] for x in rr)
            tot[ver] = np.log(o_ / y_) if o_ > 0 and y_ > 0 else np.nan
        diffs.append(tot["matched"] - tot["pooled"])
    v["H101b_matched_minus_pooled"] = v["matched_logratio_1.0"] - v["pooled_logratio_1.0"]
    if n_boot:
        v["H101b_lo"], v["H101b_hi"] = np.nanquantile(np.array(diffs), [0.025, 0.975])
        v["H101b"] = "PASS" if (v["H101b_matched_minus_pooled"] > 0 and v["H101b_lo"] > 0) else "FAIL"
    v["n_boot"] = n_boot
    pd.Series(v).to_csv(out / "a101b_summary.csv", header=False)
    print(pd.Series(v).to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--boot", type=int, default=None, help="bootstrap replicates (default 500 / 200; 10 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    o = read_obs()
    sub = G.pick_donors(o, args.subsample, args.seed)
    nb = lambda full: (10 if run["quick"] else full) if args.boot is None else args.boot
    a101(out, o, sub, args.jobs, nb(N_BOOT))
    a101b(out, o, sub, args.jobs, nb(N_BOOT_B))


if __name__ == "__main__":
    main()
