"""A143 H143.1 r* and H143.3 region alignment re-read with bounded statistics. DESCRIPTIVE / POST HOC (2026-10-29,
roadmap 5.3): registered verdicts are not changed.

r* (H143.1) is the donor correlation of brain and peripheral block scores over sqrt(rel_brain rel_periph); it is not a
`stats.ccorr`. The ccorr guard is transplanted onto the reliability (see `A142 .../rescore_a142_bounded.py`): NaN when a
split-half reliability is <= 2 SE.
H143.3 registered statistic: alignment ratio [r(k,P)/sqrt(rel_P)] / [r(k,B_-k)/sqrt(rel_B)], unbounded. Bounded
alternatives (same data, same splits):
  share  a_P / (a_P + a_B), a = r / sqrt(rel); in [0, 1] when both a > 0 (a monotone transform of the ratio: same ranking
         while a_B > 0, but no explosion at a_B -> 0); undefined when a_P or a_B <= 0 (reported NaN)
  r_P    the raw (uncorrected, |r| <= 1) correlation of region k with the peripheral score
  ratio_unc  r(k,P) / r(k,B_-k), no reliability
Hypothalamus minus median of the other 12 regions, 1,000 donor bootstraps, one-sided share <= 0, for each. Fast deaths
(registered) and Hardy 3-4 (descriptive). Aggregates only.

    python rescore_a143_bounded.py [--jobs J] [--n-boot 1000]
Outputs: output/{h143_1,h143_3_regions,h143_3_tests,breakdown}_bounded.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
for p in (PAPER, PAPER / "results/figure_2/scripts/organs", HERE.parents[0] / "scripts",
          PAPER / "exploratory/A142_organ_blocks/scripts", PAPER / "exploratory/A142_organ_blocks/bounded_2026-10-29"):
    sys.path.insert(0, str(p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402

import _shared as S  # noqa: E402
import compute_brain_artefact_checks as C  # noqa: E402
import rescore_a142_bounded as R142  # noqa: E402  (reliability with SE)

OUT = HERE / "output"
K = 2.0
HYP = C.HYPOTHALAMUS


def rstar(Z, rng, n_split=C.N_SPLIT):
    """Brain vs periphery r* as `_shared.block_stat` (same rng use), with SE-based flag."""
    isb = Z.columns.str.startswith("Brain")
    both = pd.concat([S._score(Z.loc[:, ~isb]), S._score(Z.loc[:, isb])], axis=1).dropna()
    r = both.corr().iloc[0, 1]
    (rp, sp), (rb, sb) = (R142.reliability(Z.loc[both.index, ~isb], rng, n_split), R142.reliability(Z.loc[both.index, isb], rng, n_split))
    ok = rp > 0 and rb > 0
    rs = r / np.sqrt(rp * rb) if ok else np.nan
    flag = (not ok) or rp <= K * sp or rb <= K * sb
    return dict(donors=len(both), r=r, rel_periph=rp, rel_brain=rb, r_star=rs, flag=bool(flag),
                r_star_bounded=np.nan if flag else rs, r_star_clip=min(rs, 1.0) if np.isfinite(rs) else np.nan)


def align_table(Z, brain, periph, rng, n_split):
    """Per region: old ratio and the bounded alternatives (NaN where a reliability is not > 0 or flagged)."""
    P = S._score(Z[periph])
    rows = {}
    for k in brain:
        others = [c for c in brain if c != k]
        Bk = S._score(Z[others])
        df = pd.concat([Z[k], P, Bk], axis=1, keys=["k", "P", "B"]).dropna()
        if len(df) < 20:
            continue
        (rp, sp) = R142.reliability(Z.loc[df.index, periph], rng, n_split)
        (rb, sb) = R142.reliability(Z.loc[df.index, others], rng, n_split)
        if rp <= 0 or rb <= 0:
            continue
        rP, rB = df.k.corr(df.P), df.k.corr(df.B)
        aP, aB = rP / np.sqrt(rp), rB / np.sqrt(rb)
        flag = rp <= K * sp or rb <= K * sb
        rows[k] = dict(ratio_old=aP / aB, flag=bool(flag), a_P=aP, a_B=aB,
                       share=aP / (aP + aB) if (aP > 0 and aB > 0) else np.nan, r_P=rP, ratio_unc=rP / rB if rB != 0 else np.nan,
                       ratio_old_bounded=np.nan if flag else aP / aB)
    return pd.DataFrame(rows).T.astype(float)


def hyp_stats(T: pd.DataFrame) -> dict:
    out = {}
    for col in ("ratio_old", "share", "r_P", "ratio_unc"):
        o = T.drop(index=HYP, errors="ignore")[col].dropna()
        out[col] = T.at[HYP, col] - o.median() if HYP in T.index and np.isfinite(T.at[HYP, col]) else np.nan
    return out


def _boot(Z, brain, periph, seed):
    limit_threads(1)
    rng = np.random.default_rng(seed)
    Zb = Z.iloc[rng.integers(0, len(Z), len(Z))].reset_index(drop=True)
    T = align_table(Zb, brain, periph, rng, 20)
    return dict(**hyp_stats(T), share_regions_flagged=float(T.flag.mean()) if len(T) else np.nan)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=1000)
    args = ap.parse_args(argv)
    cli.setup(args)
    rng = np.random.default_rng(C.SEED)
    pb, donors = C.S.cohort("fast")
    t, ph = C.departures_wide(pb, donors)
    att = C.sample_attributes(pb)
    zw = lambda v: S.zwide(v, t.donor.to_numpy(), t.compartment.to_numpy())  # noqa: E731
    has = att.SMNABTCH.notna().to_numpy()
    h1 = pd.DataFrame([dict(analysis="A103 base, all samples", **rstar(zw(ph), rng)),
                       dict(analysis="base, samples with attributes", **rstar(zw(np.where(has, ph, np.nan)), rng)),
                       dict(analysis="H143.1: quality adjusted", **rstar(zw(C.adjust_quality(ph, t, att)), rng))])
    old = pd.read_csv(HERE.parents[0] / "output/h143_1.csv").iloc[:3]
    h1["registered_r_star"] = old.r_star.to_numpy()
    h1["registered_hi"] = old.r_star_hi.to_numpy()
    h1.to_csv(OUT / "h143_1_bounded.csv", index=False)
    print(h1.round(4).T.to_string(), flush=True)

    Z = zw(ph)
    brain = [c for c in Z.columns if c.startswith("Brain")]
    periph = [c for c in Z.columns if c not in brain]
    T = align_table(Z, brain, periph, rng, C.N_SPLIT)
    reg = T.copy()
    for c in ("ratio_old", "share", "r_P", "ratio_unc"):
        reg[f"rank_{c}"] = reg[c].rank(ascending=False)
    reg.index.name = "region"
    # descriptive Hardy 3-4
    pbs, dons = gtex.load(hardy=(3, 4), genes="clock", depth="sum")
    n = pbs.groupby("compartment").donor.nunique()
    pbs = pbs[pbs.compartment.isin([c for c in Z.columns if n.get(c, 0) >= 20])].reset_index(drop=True)
    ts, phs = C.departures_wide(pbs, dons)
    Zs = S.zwide(phs, ts.donor.to_numpy(), ts.compartment.to_numpy())
    bs_, ps_ = [c for c in Zs.columns if c.startswith("Brain")], [c for c in Zs.columns if not c.startswith("Brain")]
    T34 = align_table(Zs, bs_, ps_, rng, C.N_SPLIT)
    for c in ("ratio_old", "share", "r_P", "ratio_unc", "flag"):
        reg[f"hardy34_{c}"] = T34[c]
        if c != "flag":
            reg[f"hardy34_rank_{c}"] = T34[c].rank(ascending=False)
    old3 = pd.read_csv(HERE.parents[0] / "output/h143_3_regions.csv").set_index("region")
    reg["registered_alignment_fast"] = old3.alignment_fast
    reg["registered_alignment_hardy34"] = old3.alignment_hardy34
    reg.to_csv(OUT / "h143_3_regions_bounded.csv")
    print(reg.round(3).to_string(), flush=True)

    obs = hyp_stats(T)
    bs = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(_boot)(Z, brain, periph, C.SEED * 10_000 + b) for b in range(args.n_boot)))
    rows = []
    for c in ("ratio_old", "share", "r_P", "ratio_unc"):
        x = bs[c].dropna().to_numpy()
        rows.append(dict(statistic=c, hyp_minus_median=obs[c], boot_lo=np.quantile(x, .025), boot_hi=np.quantile(x, .975),
                         boot_sd=x.std(ddof=1), n_boot=len(x), p_one_sided=(1 + np.sum(x <= 0)) / (1 + len(x)),
                         hyp_rank=int(T[c].rank(ascending=False)[HYP]), hyp_value=T.at[HYP, c],
                         top_region=T[c].idxmax(), share_regions_flagged_boot=float(bs.share_regions_flagged.mean())))
    Tt = pd.DataFrame(rows)
    Tt.to_csv(OUT / "h143_3_tests_bounded.csv", index=False)
    pd.DataFrame([dict(h143_1_flagged=int(h1.flag.sum()), regions=len(T), regions_flagged=int(T.flag.sum()),
                       regions_ratio_gt1=int((T.ratio_old > 1).sum()), regions_share_nan=int(T.share.isna().sum()),
                       hardy34_regions_flagged=int(T34.flag.sum()), hardy34_ratio_gt1=int((T34.ratio_old > 1).sum()))]).to_csv(
        OUT / "breakdown_bounded.csv", index=False)
    print(Tt.round(4).T.to_string(), flush=True)


if __name__ == "__main__":
    main()
