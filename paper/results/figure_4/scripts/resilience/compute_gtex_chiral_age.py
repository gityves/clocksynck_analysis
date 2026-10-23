"""A101 (core genes only): does GTEx show CHIRAL's age dampening, and is it death-type composition of the age groups?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A100-A101). All donors with Hardy 0-4. Internal phase: the host of the
clock genes of the other v2_gtex tissues (target tissue left out; `Dataset.fit().host`, all donors). Per target tissue,
amplitude of each of the 14 core clock genes = |loading| of log CPM on the host (2 columns) with the D2 covariates (sex,
age bracket, Hardy class, ischaemic time and its square); old (>= 50 y, exact age) vs young (< 50) at equal n per
tissue, in three versions: (a) pooled over death classes (CHIRAL-like), (b) composition-matched (equal n old and young
within each Hardy class), (c) fast deaths only (Hardy 1-2). D_core = median over tissues of the median over core genes of
log2(amplitude old / young), median of 20 subsamples; donor bootstrap (500, within age x Hardy cells, one draw shared
by all tissues, 4 subsamples per replicate).
  Here: D_core per version with its bootstrap interval, and pooled minus matched (H101.2 for the core genes).
  NOT here: H101.1 / H101.2 on the output genes (D, D_rel), which need the A91 clock-following gene classes and a
  covariate-adjusted per-gene regression on a fixed host (`Dataset.gene_rhythms` has neither).
Controlled access: exact age, Hardy class, ischaemic time held in memory; aggregates only.

    python compute_gtex_chiral_age.py [--subsample N] [--seed S] [--jobs 4]      (--subsample: random donors)
Outputs (paper/results/figure_4/output/resilience/): a101_core_tissues.csv, a101_core_summary.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a101_chiral_age.py (core-gene part).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "resilience"
HOSTS = ["Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)",
         "Esophagus - Mucosa", "Nerve - Tibial", "Skin - Not Sun Exposed (Suprapubic)",
         "Skin - Sun Exposed (Lower leg)", "Adipose - Subcutaneous", "Lung", "Pituitary", "Brain - Cerebellum",
         "Brain - Cortex", "Whole Blood"]
EXTRA = ["Liver", "Adrenal Gland", "Pancreas"]
GENES = list(cs.gtex().genes)
VERSIONS = ("pooled", "matched", "fast")
N_SUB, N_BOOT, N_SUB_BOOT, MIN_N = 20, 500, 4, 20


def amps(Y, Zc, f):
    Zc = Zc[:, Zc.std(0) > 0]
    D = np.c_[np.ones(len(Y)), Zc, f - f.mean(0)]
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    return np.hypot(B[-2], B[-1])


def draw(pool_old, pool_young, h_old, h_young, version, rng):
    if version == "fast":
        pool_old, pool_young = pool_old[np.isin(h_old, [1, 2])], pool_young[np.isin(h_young, [1, 2])]
    if version in ("pooled", "fast"):
        n = min(len(pool_old), len(pool_young))
        return rng.choice(pool_old, n, replace=False), rng.choice(pool_young, n, replace=False)
    o, y = [], []
    for h in np.intersect1d(h_old, h_young):
        po, py = pool_old[h_old == h], pool_young[h_young == h]
        n = min(len(po), len(py))
        o.append(rng.choice(po, n, replace=False))
        y.append(rng.choice(py, n, replace=False))
    return np.concatenate(o), np.concatenate(y)


def tissue_D(Y, Zc, f, rows_old, rows_young, h_old, h_young, version, n_sub, rng):
    """Median over subsamples of the median over core genes of log2(old / young) amplitude."""
    res = []
    for _ in range(n_sub):
        io, iy = draw(rows_old, rows_young, h_old, h_young, version, rng)
        if min(len(io), len(iy)) < MIN_N:
            return np.nan
        res.append(np.median(np.log2(amps(Y[io], Zc[io], f[io]) / amps(Y[iy], Zc[iy], f[iy]))))
    return float(np.median(res))


def host_of(pb: pd.DataFrame, target: str) -> pd.DataFrame:
    """Host state of the clock genes of the other v2_gtex tissues (all donors), donor -> (x, y)."""
    other = [t for t in HOSTS if t != target and t in set(pb.compartment)]
    ds = cs.Dataset(pb[pb.compartment.isin(other)], anchors=other, name=f"host_wo_{target}")
    return ds.fit(n_boot=0, gate=False).host[["x", "y"]]


def one(tis, pb, don, boots, seed, n_sub, n_sub_boot):
    limit_threads(1)
    rng = np.random.default_rng(seed)
    fv = host_of(pb, tis)
    x = pb[pb.compartment == tis].drop_duplicates("donor").set_index("donor")
    x = x[x.index.isin(fv.index) & x.index.isin(don.index)]
    d = don.loc[x.index]
    ok = d.age_exact.notna().to_numpy()
    x, d = x[ok], d[ok]
    f = fv.loc[x.index].to_numpy()
    Y = np.log(x[GENES].to_numpy(float) / x.depth.to_numpy(float)[:, None] * 1e6 + 0.5)
    isch = (d.ischaemic_h - d.ischaemic_h.mean()) / (d.ischaemic_h.std() or 1)
    hardy = d.hardy.to_numpy(int)
    Zc = np.column_stack([1 - d.male, d.age, (isch.fillna(0)), isch.fillna(0) ** 2] +
                         [(hardy == h).astype(float) for h in (1, 2, 3, 4)])
    old = (d.age_exact >= 50).to_numpy()
    row_of = pd.Series(np.arange(len(x)), index=x.index)

    def run(rows, ns):
        rows = np.asarray(rows)
        o_, y_ = rows[old[rows]], rows[~old[rows]]
        return {v: tissue_D(Y, Zc, f, o_, y_, hardy[o_], hardy[y_], v, ns, rng) for v in VERSIONS}

    obs = run(np.arange(len(x)), n_sub)
    bs = [run(row_of.reindex(b).dropna().astype(int).to_numpy(), n_sub_boot) for b in boots]
    comp = {f"{g}_frac_hardy{h}": float(np.mean(hardy[m] == h)) for g, m in (("old", old), ("young", ~old))
            for h in range(5)}
    return tis, obs, bs, dict(n_old=int(old.sum()), n_young=int((~old).sum()), **comp)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = 20 if run["quick"] else N_BOOT
    n_sub = 5 if run["quick"] else N_SUB
    pb, don = gtex.load(tissues=HOSTS + EXTRA, hardy=(0, 1, 2, 3, 4), depth="library_size", require_tod=False)
    don = don.join(gtex.phenotypes(["AGE"]).AGE.astype(float).rename("age_exact"))
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[:args.subsample]
        pb = pb[pb.donor.isin(keep)]
        don = don.loc[don.index.isin(keep)]
    tissues = [t for t in HOSTS + EXTRA if (pb.compartment == t).sum() >= (60 if not run["quick"] else 30)]
    dd = don.dropna(subset=["age_exact"])
    cell = (dd.age_exact >= 50).astype(str) + "|" + dd.hardy.astype(int).astype(str)
    rng = np.random.default_rng(101 + args.seed)
    boots = [np.concatenate([rng.choice(g.index.to_numpy(), len(g)) for _, g in dd.groupby(cell.to_numpy())])
             for _ in range(n_boot)]
    res = Parallel(n_jobs=args.jobs)(delayed(one)(t, pb, don, boots, 101 + i, n_sub, N_SUB_BOOT)
                                     for i, t in enumerate(tissues))
    per = pd.DataFrame([dict(tissue=t, **{f"Dcore_{v}": o[v] for v in VERSIONS}, **c) for t, o, _, c in res])
    per.to_csv(out / "a101_core_tissues.csv", index=False)
    agg = lambda tab: {v: float(np.nanmedian([t[v] for t in tab])) for v in VERSIONS}
    obs = agg([r[1] for r in res])
    boot = [agg([r[2][b] for r in res]) for b in range(n_boot)]
    s = {}
    for v in VERSIONS:
        x = np.array([b[v] for b in boot])
        s[f"Dcore_{v}"], s[f"Dcore_{v}_lo"], s[f"Dcore_{v}_hi"] = obs[v], *np.nanquantile(x, [0.025, 0.975])
    diff = np.array([b["pooled"] - b["matched"] for b in boot])
    s["pooled_minus_matched"] = obs["pooled"] - obs["matched"]
    s["pooled_minus_matched_lo"], s["pooled_minus_matched_hi"] = np.nanquantile(diff, [0.025, 0.975])
    s["tissues"], s["n_boot"] = len(per), n_boot
    pd.Series(s).to_csv(out / "a101_core_summary.csv", header=False)
    print(per.round(3).to_string(index=False), "\n", pd.Series(s).round(3).to_string())


if __name__ == "__main__":
    main()
