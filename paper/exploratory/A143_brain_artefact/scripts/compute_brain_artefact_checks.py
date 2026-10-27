"""A143: is the brain-periphery split a processing artefact, and does the hypothalamus follow the body?

Registered in paper/docs/registered_experiments.md (A143). Cohort, departures and r* as A103 base
(`results/figure_2/scripts/organs/compute_brain_artefact.py`, `_shared.block_stat` / `block_boot`): GTEx v10 fast deaths,
tissues with >= 40 donors, linear phase departures of the 14-gene plane, standardised per tissue. Sample attributes (v8,
controlled access) joined on the sample behind each (donor, tissue) row; samples without attributes dropped from H143.1-2.
- H143.1: departures adjusted per tissue for RIN, autolysis (factor), ischaemic time, Hardy 2 vs 1 and centre (where they
  vary); r* with 200 donor bootstraps; base r* on the same samples for comparison.
- H143.2: per brain-region pair, correlation of departures over donors whose two samples share an isolation batch
  (SMNABTCH) minus over donors whose samples do not; weighted mean over pairs with >= 10 donors per class (weight: harmonic
  mean of the counts); labels permuted across donors within pair (2,000). Peripheral pairs as reference.
- H143.3: alignment ratio of each brain region, [r(k, P) / sqrt(rel_P)] / [r(k, B_-k) / sqrt(rel_B-k)]; hypothalamus minus
  the median of the other regions, 1,000 donor bootstraps, one-sided p. Descriptive: all regions, and Hardy 3-4 deaths.
Holm over H143.1-3. Aggregates only: no per-donor value is written or printed.

    python compute_brain_artefact_checks.py [--subsample N] [--jobs J]
Outputs: paper/exploratory/A143_brain_artefact/output/{coverage,h143_1,h143_2,h143_2_pairs,h143_3_regions,tests}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
ORGANS = PAPER / "results/figure_2/scripts/organs"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(ORGANS))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import GTEX_PHENOTYPES, exploratory_dir  # noqa: E402

import _shared as S  # noqa: E402
import compute_brain_artefact as A103  # noqa: E402

AMENDMENT = "A143_brain_artefact"
ATTRIBUTES = GTEX_PHENOTYPES / "phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt"
HYPOTHALAMUS = "Brain - Hypothalamus"
N_SPLIT, SEED, MIN_CLASS = A103.N_SPLIT, 143, 10


def sample_attributes(pb: pd.DataFrame) -> pd.DataFrame:
    """RIN, autolysis, centre and isolation batch of the sample behind each row of `pb` (first sample of a duplicate, as
    `gtex.load` and `_shared.sample_ischaemic_h`)."""
    import anndata as ad

    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs[["subject", "tissue"]].copy()
        o["SAMPID"] = a.obs_names.astype(str)
    finally:
        a.file.close()
    o["donor"], o["compartment"] = o.subject.astype(str), o.tissue.astype(str)
    o = o.drop_duplicates(["donor", "compartment"])
    att = pd.read_csv(ATTRIBUTES, sep="\t", comment="#", low_memory=False,
                      usecols=["SAMPID", "SMRIN", "SMATSSCR", "SMCENTER", "SMNABTCH"])
    o = o.merge(att, on="SAMPID", how="left").set_index(["donor", "compartment"])
    return o.reindex(pd.MultiIndex.from_frame(pb[["donor", "compartment"]])).reset_index(drop=True)


def adjust_quality(values, t: pd.DataFrame, att: pd.DataFrame, autolysis: bool = False):
    """Per-tissue OLS residuals (plus intercept) on RIN, ischaemic time, Hardy 2 vs 1 and centre (note 1), plus the autolysis
    factor with `autolysis` (only tissues with >= 70 % of samples scored); NaN where any covariate is missing."""
    out = np.full(len(values), np.nan)
    v = np.asarray(values, float)
    for _, idx in t.groupby("compartment").indices.items():
        X = pd.DataFrame({"rin": att.SMRIN.to_numpy()[idx], "isch": t.isch.to_numpy()[idx],
                          "h2": (t.hardy.to_numpy()[idx] == 2).astype(float)})
        cat = pd.DataFrame({"centre": att.SMCENTER.to_numpy()[idx]})
        if autolysis:
            if att.SMATSSCR.iloc[idx].notna().mean() < 0.7:
                continue
            cat["auto"] = att.SMATSSCR.to_numpy()[idx]
        ok = X.notna().all(1).to_numpy() & cat.notna().all(1).to_numpy() & np.isfinite(v[idx])
        if ok.sum() < 10:
            continue
        D = pd.get_dummies(cat[ok].astype(str), drop_first=True, dtype=float)
        M = pd.concat([X[ok].reset_index(drop=True), D.reset_index(drop=True)], axis=1)
        M = M.loc[:, M.std() > 0]
        M.insert(0, "one", 1.0)
        B, *_ = np.linalg.lstsq(M.to_numpy(), v[idx][ok], rcond=None)
        out[idx[ok]] = v[idx][ok] - M.to_numpy() @ B + B[0]
    return out


def rstar_row(Z, label, rng, n_boot):
    pt = S.block_stat(Z, rng, N_SPLIT)
    lo, hi = S.block_boot(Z, rng, n_boot)
    return dict(analysis=label, donors=pt["donors"], r=pt["r"], rel_periph=pt["rel_periph"], rel_brain=pt["rel_brain"],
                r_star=pt["r_star"], r_star_lo=lo, r_star_hi=hi, n_boot=n_boot)


def batch_pairs(Z: pd.DataFrame, Bt: pd.DataFrame, tissues: list) -> list:
    """Per tissue pair: (x, y, same-batch labels) over donors with both departures and both batches."""
    pairs = []
    for i, a in enumerate(tissues):
        for b in tissues[i + 1:]:
            ok = Z[a].notna() & Z[b].notna() & Bt[a].notna() & Bt[b].notna()
            same = (Bt.loc[ok, a] == Bt.loc[ok, b]).to_numpy()
            if same.sum() >= MIN_CLASS and (~same).sum() >= MIN_CLASS:
                pairs.append((a, b, Z.loc[ok, a].to_numpy(float), Z.loc[ok, b].to_numpy(float), same))
    return pairs


def _corr(x, y):
    return np.corrcoef(x, y)[0, 1]


def batch_stat(pairs, labels=None):
    num = den = 0.0
    rows = []
    for j, (a, b, x, y, same) in enumerate(pairs):
        s = same if labels is None else labels[j]
        ns, nd = s.sum(), (~s).sum()
        d = _corr(x[s], y[s]) - _corr(x[~s], y[~s])
        w = 2 / (1 / ns + 1 / nd)
        num, den = num + w * d, den + w
        rows.append(dict(tissue_a=a, tissue_b=b, n_same=int(ns), n_diff=int(nd), diff=d))
    return num / den, rows


def batch_test(pairs, n_perm, rng):
    obs, rows = batch_stat(pairs)
    null = np.array([batch_stat(pairs, [rng.permutation(p[4]) for p in pairs])[0] for _ in range(n_perm)])
    return dict(n_pairs=len(pairs), stat=obs, p_one_sided=(1 + np.sum(null >= obs)) / (1 + n_perm), null_mean=null.mean(),
                null_sd=null.std(ddof=1), null_q95=np.quantile(null, .95)), rows


def alignment(Z: pd.DataFrame, brain: list, periph: list, rng, n_split: int) -> dict:
    """Alignment ratio of every brain region: its disattenuated correlation with the peripheral score over that with the
    score of the other brain regions (the region's own reliability cancels)."""
    P = S._score(Z[periph])
    out = {}
    for k in brain:
        others = [c for c in brain if c != k]
        Bk = S._score(Z[others])
        df = pd.concat([Z[k], P, Bk], axis=1, keys=["k", "P", "B"]).dropna()
        if len(df) < 20:
            continue
        rp = S.reliability(Z.loc[df.index, periph], rng, n_split)
        rb = S.reliability(Z.loc[df.index, others], rng, n_split)
        if rp <= 0 or rb <= 0:
            continue
        out[k] = (df.k.corr(df.P) / np.sqrt(rp)) / (df.k.corr(df.B) / np.sqrt(rb))
    return out


def hypothalamus_stat(al: dict) -> float:
    others = [v for k, v in al.items() if k != HYPOTHALAMUS and np.isfinite(v)]
    return al.get(HYPOTHALAMUS, np.nan) - np.median(others)


def _boot_align(Z, brain, periph, seed):
    limit_threads(1)
    rng = np.random.default_rng(seed)
    Zb = Z.iloc[rng.integers(0, len(Z), len(Z))].reset_index(drop=True)
    return hypothalamus_stat(alignment(Zb, brain, periph, rng, 20))


def departures_wide(pb, donors):
    d = donors.reindex(pb.donor)
    t = pd.DataFrame({"donor": pb.donor.to_numpy(), "compartment": pb.compartment.to_numpy(), "tod": d.tod.to_numpy(),
                      "hardy": d.hardy.to_numpy(), "isch": S.sample_ischaemic_h(pb)})
    ph, _ = A103.phase_departures(pb, t, list(S.CLOCK))
    return t, ph


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    n_boot, n_perm, n_boot3 = (10, 100, 20) if quick else (200, 2000, 1000)
    out = exploratory_dir(AMENDMENT, quick)
    rng = np.random.default_rng(SEED)

    pb, donors = S.cohort("fast", subsample=args.subsample, seed=args.seed)
    t, ph = departures_wide(pb, donors)
    att = sample_attributes(pb)
    zw = lambda v: S.zwide(v, t.donor.to_numpy(), t.compartment.to_numpy())  # noqa: E731
    isb = t.compartment.str.startswith("Brain").to_numpy()
    has = att.SMNABTCH.notna().to_numpy()
    cov = pd.DataFrame([dict(block=b, samples=int(m.sum()), with_attributes=int((m & has).sum()),
                             share=float((m & has).sum() / m.sum())) for b, m in (("brain", isb), ("periphery", ~isb))])
    cov.to_csv(out / "coverage.csv", index=False)
    print(cov.to_string(index=False), flush=True)
    tests = []

    # ---- H143.1 ---------------------------------------------------------------------------------------------------
    rows = [rstar_row(zw(ph), "A103 base, all samples", rng, n_boot),
            rstar_row(zw(np.where(has, ph, np.nan)), "base, samples with attributes", rng, n_boot),
            rstar_row(zw(adjust_quality(ph, t, att)), "H143.1: quality adjusted", rng, n_boot)]
    za = zw(adjust_quality(ph, t, att, autolysis=True))
    n_brain_auto = int(za.columns[za.notna().any()].str.startswith("Brain").sum())
    if n_brain_auto >= 2:
        rows.append(rstar_row(za, f"descriptive: + autolysis ({n_brain_auto} brain regions)", rng, n_boot))
    else:
        print(f"autolysis sensitivity not computed: {n_brain_auto} brain regions with >= 70 % scored", flush=True)
    h1 = pd.DataFrame(rows)
    h1.to_csv(out / "h143_1.csv", index=False)
    print(h1.round(3).to_string(index=False), flush=True)
    estimable = cov.set_index("block").share["brain"] >= 0.5
    r1 = h1.iloc[2]
    tests.append(dict(hypothesis="H143.1", statistic=r1.r_star, lo=r1.r_star_lo, hi=r1.r_star_hi, p_one_sided=np.nan,
                      pass_rule="r* < 1 and upper < 1",
                      verdict=("not estimable" if not estimable else
                               "pass" if (r1.r_star < 1 and r1.r_star_hi < 1) else "fail")))

    # ---- H143.2 ---------------------------------------------------------------------------------------------------
    Z = zw(ph)
    Bt = S.wide(att.SMNABTCH.to_numpy(), t.donor.to_numpy(), t.compartment.to_numpy()).reindex(index=Z.index, columns=Z.columns)
    brain = [c for c in Z.columns if c.startswith("Brain")]
    periph = [c for c in Z.columns if c not in brain]
    h2, pair_rows = [], []
    for name, tis in (("brain", brain), ("periphery", periph)):
        res, pr = batch_test(batch_pairs(Z, Bt, tis), n_perm, rng)
        h2.append(dict(block=name, **res))
        pair_rows += [dict(block=name, **r) for r in pr]
    h2 = pd.DataFrame(h2)
    h2.to_csv(out / "h143_2.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(out / "h143_2_pairs.csv", index=False)
    print(h2.round(4).to_string(index=False), flush=True)
    b2 = h2.set_index("block").loc["brain"]
    tests.append(dict(hypothesis="H143.2", statistic=b2.stat, lo=np.nan, hi=np.nan, p_one_sided=b2.p_one_sided,
                      pass_rule="p < 0.05 (a pass = batch-borne departure)",
                      verdict="not estimable" if (not estimable or b2.n_pairs == 0) else None))

    # ---- H143.3 ---------------------------------------------------------------------------------------------------
    al = alignment(Z, brain, periph, rng, N_SPLIT)
    obs = hypothalamus_stat(al)
    bs = np.array(Parallel(n_jobs=args.jobs)(delayed(_boot_align)(Z, brain, periph, SEED * 10_000 + b) for b in range(n_boot3)))
    bs = bs[np.isfinite(bs)]
    reg = pd.DataFrame({"region": list(al), "alignment_fast": list(al.values())})
    try:  # descriptive: Hardy 3-4 deaths, same tissues where >= 20 donors
        pbs, dons = gtex.load(hardy=(3, 4), genes="clock", depth="sum")
        n = pbs.groupby("compartment").donor.nunique()
        pbs = pbs[pbs.compartment.isin([c for c in Z.columns if n.get(c, 0) >= 20])].reset_index(drop=True)
        ts, phs = departures_wide(pbs, dons)
        Zs = S.zwide(phs, ts.donor.to_numpy(), ts.compartment.to_numpy())
        bs_, ps_ = [c for c in Zs.columns if c.startswith("Brain")], [c for c in Zs.columns if not c.startswith("Brain")]
        reg["alignment_hardy34"] = reg.region.map(alignment(Zs, bs_, ps_, rng, N_SPLIT))
    except Exception as e:  # descriptive only: reported, never silent
        print(f"Hardy 3-4 descriptive not computed: {type(e).__name__}: {e}", flush=True)
    reg.to_csv(out / "h143_3_regions.csv", index=False)
    print(reg.round(3).to_string(index=False), flush=True)
    p3 = (1 + np.sum(bs <= 0)) / (1 + len(bs))
    tests.append(dict(hypothesis="H143.3", statistic=obs, lo=np.quantile(bs, .025), hi=np.quantile(bs, .975), p_one_sided=p3,
                      pass_rule="one-sided bootstrap p < 0.05", verdict=None))

    # ---- Holm over the p values (H143.1 is a CI rule: p set to 0 if pass else 1 for the ordering) -----------------------
    T = pd.DataFrame(tests)
    pp = T.p_one_sided.fillna(T.verdict.map({"pass": 0.0, "fail": 1.0})).fillna(1.0)
    order, run_max, adj = np.argsort(pp.to_numpy()), 0.0, np.empty(len(T))
    for i, j in enumerate(order):
        run_max = max(run_max, min(1.0, (len(T) - i) * pp.iloc[j]))
        adj[j] = run_max
    T["p_holm"] = adj
    for j in range(len(T)):
        if pd.isna(T.verdict.iloc[j]):
            T.loc[T.index[j], "verdict"] = "pass" if T.p_holm.iloc[j] < 0.05 else "fail"
    T.to_csv(out / "tests.csv", index=False)
    print(T.round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
