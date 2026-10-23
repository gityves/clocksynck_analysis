"""WP0 item 5 (integrity check, no new hypothesis): robustness of the CHIRAL reanalysis (archive A101b, A102b).

(1) The death-type-matched rhythmic-gene count comparison (old 60-79 vs young < 50; R1, pooled, Hardy-matched, Hardy x
    RIN-matched) of `compute_chiral_prep_counts.py`, with CHIRAL's normalisation (edgeR TMM on the genes with mean count
    > 10, log2 CPM, prior 0.25; `_tmm.py`) beside library-size CPM. Same tissues, seeds (SEED + tissue index) and draws
    (10) for both, so the CPM rows reproduce a102b_a101b_tissues.csv and the two normalisations are paired.
(2) CHIRAL's sex comparison (female vs male rhythmic-gene counts; R1 = CHIRAL's sex subsample mfSS) under the same
    matching: random equal n (pooled), equal n within Hardy class, Hardy x RIN tertile. Statistic log(N_female / N_male).
Everything else as `compute_chiral_prep_counts.counts_nan` (residuals of ischaemic-time bins, age, sex and Hardy; entries > 5
log2 from the gene mean excluded; harmonic regression on CHIRAL's CT; LRT; BH q < 0.2 and log2 peak-to-trough > 1).
Statistics: per-tissue log ratio log((a + 1) / (b + 1)) Wilcoxon signed-rank p; total log ratio with a tissue bootstrap
(2000, seed 5) two-sided p; matched minus pooled (same tissues) with the tissue bootstrap. Aggregates only: CT, age, sex,
Hardy are per-donor controlled data held in memory.

    python compute_chiral_tmm.py [--subsample N] [--jobs 4]
Outputs (paper/exploratory/WP0_integrity/output/chiral_tmm/, .../subsample/ for checks): count_comparison_tmm.csv,
count_comparison_tmm_summary.csv, sex_comparison_matched.csv, sex_comparison_matched_summary.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "paper"))
sys.path.insert(0, str(REPO / "paper" / "results" / "figure_3" / "scripts" / "gtex_cues"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

import _genes as G  # noqa: E402
import compute_chiral_age as CA  # noqa: E402
import compute_chiral_prep_counts as PC  # noqa: E402
import _tmm  # noqa: E402

SEED, N_DRAW, N_BOOT, MIN_N = PC.SEED, 10, 2000, 25
OUT = REPO / "paper" / "exploratory" / "WP0_integrity" / "output" / "chiral_tmm"
NORMS = ("CPM", "TMM")


def expression(X, keep, ob, norm):
    """Residualised log2 CPM (prior 0.25) of the kept genes under `norm`."""
    if norm == "CPM":  # as compute_chiral_prep_counts (library = all genes of the sample)
        L = np.log2((X[:, keep] + 0.25) / (X.sum(1, keepdims=True) + 0.5) * 1e6)
    else:  # CHIRAL To_CPM.R: filter genes, DGEList, calcNormFactors, cpm(log = TRUE, prior.count = 0.25)
        Xk = X[:, keep]
        L = _tmm.log2_cpm(Xk, _tmm.tmm_factors(Xk), 0.25)
    return CA.residualise(L, ob).astype(np.float32)


def one(tis, seed, meta, sub):
    limit_threads(1)
    cpm, genes, ob = G.load_rows(tis, (0, 1, 2, 3, 4), sub, need_tod=False)
    ob = ob.join(meta[["CT", "ageSS", "mfSS"]])
    ok = (ob.CT.notna() & (ob.rin > 4) & ob.ischemic_min.notna()).to_numpy()
    cpm, ob = cpm[ok], ob[ok]
    sex = ob.sex.to_numpy(str)
    n_age = [ob.age_bracket.isin(CA.OLD).sum(), ob.age_bracket.isin(CA.YOUNG).sum()]
    n_sex = [(sex == "female").sum(), (sex == "male").sum()]
    do_age, do_sex = min(n_age) >= MIN_N, min(n_sex) >= MIN_N
    if not (do_age or do_sex):
        return None
    X = np.round(cpm * ob.library_size.to_numpy(float)[:, None] / 1e6)
    keep = X.mean(0) > 10
    E = {k: expression(X, keep, ob, k) for k in NORMS}
    del X, cpm
    phi = ob.CT.to_numpy(float) * np.pi / 12
    hardy = ob.hardy.to_numpy(int)
    rin_t = pd.qcut(ob.rin.astype(float).rank(method="first"), 3, labels=False).to_numpy()
    keys = dict(pooled=hardy, matched=hardy, matched_rin=hardy * 10 + rin_t)
    cnt = lambda k, idx: PC.counts_nan(E[k], phi, idx)

    def comparison(g1, g2, flag, rng):
        """Rows per normalisation: R1 and the three matched versions, counts (g1, g2) at equal n, median of draws."""
        rows = {k: {} for k in NORMS}
        for lab, g in ((0, g1), (1, g2)):
            idx = np.flatnonzero(g & (ob[flag].fillna(0).to_numpy() == 1))
            for k in NORMS:
                rows[k][("R1", lab)] = cnt(k, idx) if len(idx) >= MIN_N else np.nan
        r1, r2 = np.flatnonzero(g1), np.flatnonzero(g2)
        for v in ("pooled", "matched", "matched_rin"):
            acc = {k: [] for k in NORMS}
            for _ in range(N_DRAW):
                i1, i2 = CA.draw_b(r1, r2, keys[v][r1], keys[v][r2], "pooled" if v == "pooled" else "matched", rng)
                if min(len(i1), len(i2)) < MIN_N:
                    break
                for k in NORMS:
                    acc[k].append((cnt(k, i1), cnt(k, i2), len(i1)))
            if acc["CPM"]:
                for k in NORMS:
                    a = np.array(acc[k], float)
                    rows[k][(v, 0)], rows[k][(v, 1)], rows[k][(v, "n")] = np.median(a[:, 0]), np.median(a[:, 1]), a[0, 2]
        return rows

    out = {}
    rng = np.random.default_rng(seed)  # age draws first, exactly as compute_chiral_prep_counts
    if do_age:
        grp_old = ob.age_bracket.isin(CA.OLD).to_numpy()
        grp_young = ob.age_bracket.isin(CA.YOUNG).to_numpy()
        out["age"] = comparison(grp_old, grp_young, "ageSS", rng)
    if do_sex:
        out["sex"] = comparison(sex == "female", sex == "male", "mfSS", rng)
    return tis, out


def table(res, kind, names):
    """One row per (tissue, normalisation); columns R1_<a/b>, <version>_<a/b> (+ _n) like a102b_a101b_tissues.csv."""
    a, b = names  # age: (old, young), sex: (female, male)
    rows = []
    for tis, out in res:
        if kind not in out:
            continue
        for k in NORMS:
            r = out[kind][k]
            row = dict(tissue=tis, normalisation=k)
            for ver in ("R1", "pooled", "matched", "matched_rin"):
                if (ver, 0) in r:
                    row[f"{ver}_{a}"], row[f"{ver}_{b}"] = r[(ver, 0)], r[(ver, 1)]
                    if kind == "sex" and ver != "R1":
                        row[f"{ver}_n"] = r[(ver, "n")]
            rows.append(row)
    d = pd.DataFrame(rows)
    if kind == "age":  # column order of the CPM table: R1_young, R1_old, pooled_old, pooled_young, ...
        order = ["tissue", f"R1_{b}", f"R1_{a}"] + [f"{v}_{x}" for v in ("pooled", "matched", "matched_rin") for x in (a, b)]
        d = d.reindex(columns=order + ["normalisation"])
    return d


def boot_p(x, rng):
    """Tissue-bootstrap interval and two-sided p of a statistic stored as a vector of draws."""
    lo, hi = np.quantile(x, [0.025, 0.975])
    return lo, hi, float(min(1.0, 2 * min(np.mean(x <= 0), np.mean(x >= 0))))


def summarise(d, a, b, tag):
    """Per normalisation and version: total log(a / b), tissue-bootstrap CI and p, median tissue log ratio with Wilcoxon p;
    matched minus pooled and matched_rin minus pooled (same tissues) with the tissue bootstrap."""
    rng = np.random.default_rng(5)
    rows = []
    for k in NORMS:
        x = d[d.normalisation == k].set_index("tissue")
        idx_boot = {}
        for ver in ("R1", "pooled", "matched", "matched_rin"):
            if f"{ver}_{a}" not in x:
                continue
            y = x.dropna(subset=[f"{ver}_{a}"])
            A, B = y[f"{ver}_{a}"].to_numpy(float), y[f"{ver}_{b}"].to_numpy(float)
            yy = y[(A + B) >= 20]
            lr = np.log((yy[f"{ver}_{a}"] + 1) / (yy[f"{ver}_{b}"] + 1)).to_numpy(float)
            bt = rng.integers(0, len(y), (N_BOOT, len(y))) if len(y) else None
            tot = np.log(A.sum() / B.sum()) if len(y) else np.nan
            bs = np.log(A[bt].sum(1) / B[bt].sum(1)) if len(y) else np.array([np.nan])
            lo, hi, p = boot_p(bs[np.isfinite(bs)], rng)
            wp = stats.wilcoxon(lr).pvalue if len(lr) >= 6 and np.any(lr != 0) else np.nan
            rows.append(dict(comparison=tag, normalisation=k, version=ver, statistic="total_logratio", tissues=len(y),
                             estimate=tot, lo=lo, hi=hi, p=p))
            rows.append(dict(comparison=tag, normalisation=k, version=ver, statistic="median_tissue_logratio",
                             tissues=len(lr), estimate=float(np.median(lr)) if len(lr) else np.nan, lo=np.nan, hi=np.nan,
                             p=wp))
        for ver in ("matched", "matched_rin"):
            if f"{ver}_{a}" not in x:
                continue
            y = x.dropna(subset=[f"{ver}_{a}", f"pooled_{a}"])
            n = len(y)
            if n < 3:
                continue
            P = np.c_[y[f"pooled_{a}"], y[f"pooled_{b}"]].astype(float)
            M = np.c_[y[f"{ver}_{a}"], y[f"{ver}_{b}"]].astype(float)
            stat = lambda ix: np.log(M[ix, 0].sum() / M[ix, 1].sum()) - np.log(P[ix, 0].sum() / P[ix, 1].sum())
            bs = np.array([stat(rng.integers(0, n, n)) for _ in range(N_BOOT)])
            lo, hi, p = boot_p(bs[np.isfinite(bs)], rng)
            rows.append(dict(comparison=tag, normalisation=k, version=f"{ver}_minus_pooled",
                             statistic="total_logratio_difference", tissues=n, estimate=stat(np.arange(n)), lo=lo, hi=hi,
                             p=p))
    return pd.DataFrame(rows)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = OUT / "subsample" if run["quick"] else OUT
    out.mkdir(parents=True, exist_ok=True)
    o = G.obs_table()
    sub = G.pick_donors(o, args.subsample, args.seed)
    meta = pd.read_csv(CA.META, sep=" ")
    meta["sample"] = meta.fullID.str.replace(".", "-", regex=False)
    meta = meta.set_index("sample")[["CT", "ageSS", "mfSS"]]
    tissues = sorted(t for t in o.tissue.unique() if not t.startswith(("Cells", "Whole Blood")))
    res = [r for r in Parallel(n_jobs=args.jobs)(delayed(one)(t, SEED + i, meta, sub) for i, t in enumerate(tissues))
           if r is not None]
    age = table(res, "age", ("old", "young"))
    sex = table(res, "sex", ("female", "male"))
    age.to_csv(out / "count_comparison_tmm.csv", index=False)
    sex.to_csv(out / "sex_comparison_matched.csv", index=False)
    s_age = summarise(age, "old", "young", "age_old_vs_young") if len(age) else pd.DataFrame()
    s_sex = summarise(sex, "female", "male", "sex_female_vs_male") if len(sex) else pd.DataFrame()
    s_age.to_csv(out / "count_comparison_tmm_summary.csv", index=False)
    s_sex.to_csv(out / "sex_comparison_matched_summary.csv", index=False)
    pd.set_option("display.width", 200)
    print(s_age.round(3).to_string(index=False), "\n", s_sex.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
