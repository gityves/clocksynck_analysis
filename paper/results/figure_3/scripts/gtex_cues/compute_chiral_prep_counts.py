"""A102b, rerun 3 (A101b on CHIRAL's preprocessing, with the NaN rule): rhythmic-gene counts on CHIRAL's phases.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A102b item 3: "A101b: counts above log2 peak-to-trough 1 on CHIRAL's phases,
pooled, Hardy-matched and Hardy x RIN-matched, with the NA rule"; compared with A101b and A102.4; descriptive). Per tissue
(donors with a CHIRAL phase CT, Hardy 0-4, RIN > 4, ischaemic time known; not Whole Blood or cultured cells; >= 25 old and young):
log2 CPM (prior 0.25), genes with mean count > 10, residuals of ischaemic-time bins, age, sex and Hardy class
(`compute_chiral_age.residualise`), entries more than 5 log2 units from the gene's mean excluded per gene (`mask5`, NaN rule),
harmonic regression on CT per group by least squares on the finite rows, LRT p, BH q; rhythmic = q < 0.2 and log2 peak-to-trough
> 1. Groups: young (< 50) and old (60-79) at equal n, 10 draws (median): R1 (CHIRAL's age subsample ageSS), pooled over death
classes, Hardy-matched, Hardy x RIN-tertile-matched. Statistic: log(total old / total young) over tissues and median tissue
log ratio (tissues with >= 20 calls). Aggregates only (the phases are per-donor controlled data and stay in memory).

    python compute_chiral_prep_counts.py [--subsample N] [--jobs 4]
Outputs (paper/results/figure_3/output/gtex_cues/): a102b_a101b_tissues.csv, a102b_a101b_summary.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a102b_chiral_prep.py (`step_a101b`). Local: the NaN-rule least squares
(`compute_covariate_sensitivity.ols_nan`), not in clocksync.
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

from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _genes as G  # noqa: E402
import compute_chiral_age as CA  # noqa: E402
import compute_covariate_sensitivity as CS  # noqa: E402
from _genes import FOLDER  # noqa: E402

SEED = 102


def counts_nan(E, phi, idx) -> int:
    """Rhythmic genes (BH q < 0.2, log2 peak-to-trough > 1) of the rows `idx` with entries > 5 log2 from the mean excluded."""
    Y = CS.mask5(E[idx])
    D = np.c_[np.ones(len(idx)), np.cos(phi[idx]), np.sin(phi[idx])]
    Bm, rss1, df, _ = CS.ols_nan(Y, D, [1, 2])
    _, rss0, _, _ = CS.ols_nan(Y, D[:, :1], [0])
    p = stats.chi2.sf((df + 3) * np.log(rss0 / rss1), 2)
    ok = np.isfinite(p)
    q = np.full(len(p), 1.0)
    pp = p[ok]
    order = np.argsort(pp)
    ranked = pp[order] * len(pp) / np.arange(1, len(pp) + 1)
    qq = np.empty_like(pp)
    qq[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    q[ok] = qq
    return int(np.sum((q < 0.2) & (2 * np.hypot(Bm[1], Bm[2]) > 1.0)))


def one(tis, seed, meta, sub):
    limit_threads(1)
    cpm, genes, ob = G.load_rows(tis, (0, 1, 2, 3, 4), sub, need_tod=False)
    ob = ob.join(meta[["CT", "ageSS"]])
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
    row = dict(tissue=tis)
    for g in ("young", "old"):
        idx = np.flatnonzero((grp == g) & (ob.ageSS.fillna(0).to_numpy() == 1))
        row[f"R1_{g}"] = counts_nan(E, phi, idx) if len(idx) >= 25 else np.nan
    ro, ry = np.flatnonzero(grp == "old"), np.flatnonzero(grp == "young")
    for v, key in (("pooled", hardy), ("matched", hardy), ("matched_rin", hardy * 10 + rin_t)):
        acc = []
        for _ in range(10):
            io, iy = CA.draw_b(ro, ry, key[ro], key[ry], "pooled" if v == "pooled" else "matched", rng)
            if min(len(io), len(iy)) < 25:
                break
            acc.append((counts_nan(E, phi, io), counts_nan(E, phi, iy)))
        if acc:
            row[f"{v}_old"], row[f"{v}_young"] = np.median(np.array(acc, float), axis=0)
    return row


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    o = G.obs_table()
    sub = G.pick_donors(o, args.subsample, args.seed)
    meta = pd.read_csv(CA.META, sep=" ")
    meta["sample"] = meta.fullID.str.replace(".", "-", regex=False)
    meta = meta.set_index("sample")[["CT", "ageSS"]]
    tissues = sorted(t for t in o.tissue.unique() if not t.startswith(("Cells", "Whole Blood")))
    per = pd.DataFrame([r for r in Parallel(n_jobs=args.jobs)(delayed(one)(t, SEED + i, meta, sub)
                                                              for i, t in enumerate(tissues)) if r is not None])
    per.to_csv(out / "a102b_a101b_tissues.csv", index=False)
    v = {}
    for ver in ("R1", "pooled", "matched", "matched_rin"):
        x = per.dropna(subset=[f"{ver}_old"])
        v[f"{ver}_total_logratio"] = float(np.log(x[f"{ver}_old"].sum() / x[f"{ver}_young"].sum()))
        x = x[(x[f"{ver}_old"] + x[f"{ver}_young"]) >= 20]
        v[f"{ver}_median_tissue_logratio"] = float(np.median(np.log((x[f"{ver}_old"] + 1) / (x[f"{ver}_young"] + 1))))
        v[f"{ver}_tissues"] = len(x)
    pd.Series(v).to_csv(out / "a102b_a101b_summary.csv", header=False)
    print(per.round(1).to_string(index=False), "\n", pd.Series(v).to_string())


if __name__ == "__main__":
    main()
