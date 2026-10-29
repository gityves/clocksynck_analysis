"""A142 block statistic r* re-read with a bounded statistic. DESCRIPTIVE / POST HOC (2026-10-29, roadmap 5.3): the registered
verdicts are not changed.

A142's r* is NOT a `clocksync.stats.ccorr`: it is the donor correlation r of two block scores divided by
sqrt(rel_A rel_B), the split-half (Spearman-Brown) reliabilities (`_shared.reliability`). The analogue of the ccorr guard is
applied to the reliability: with rr the split-half correlation over n donors, SE(rr) = (1 - rr^2) / sqrt(n - 2) and
SE(rel) = 2 SE(rr) / (1 + rr)^2; a block is FLAGGED (r* = NaN) when rel <= 2 SE(rel) (the correction is indistinguishable
from 0, as in `ccorr`). Bounded alternatives reported next to r*: the uncorrected donor correlation r (|r| <= 1), r*
truncated at 1. Same cohort, plane (base), seeds and reference blocks as `compute_organ_blocks.py` (the rng streams are
reproduced; r* must equal the registered value).

    python rescore_a142_bounded.py [--n-ref 1000] [--jobs J]
Outputs: output/{targets,reference,reference_summary,tests}_bounded.csv  (aggregates only)
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SC = HERE.parents[0] / "scripts"
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_2/scripts/organs"))
sys.path.insert(0, str(SC))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

import _shared as S  # noqa: E402
import compute_organ_blocks as M  # noqa: E402
import compute_brain_artefact as A103  # noqa: E402

OUT = HERE / "output"
K = 2.0


def reliability(X: pd.DataFrame, rng, n_split: int):
    """As `_shared.reliability` (same rng use), also the mean SE of the Spearman-Brown value (delta method)."""
    r, se = [], []
    for _ in range(n_split):
        c = rng.permutation(X.columns)
        h = pd.concat([S._score(X[c[: len(c) // 2]]), S._score(X[c[len(c) // 2:]])], axis=1).dropna()
        rr = h.corr().iloc[0, 1]
        r.append(2 * rr / (1 + rr))
        se.append(2 * (1 - rr ** 2) / np.sqrt(max(len(h) - 2, 1)) / (1 + rr) ** 2)
    return float(np.mean(r)), float(np.mean(se))


def block_stat(Z, a, b, rng):
    both = pd.concat([S._score(Z[a]), S._score(Z[b])], axis=1).dropna()
    r = both.corr().iloc[0, 1]
    (ra, sa), (rb, sb) = reliability(Z.loc[both.index, a], rng, M.N_SPLIT), reliability(Z.loc[both.index, b], rng, M.N_SPLIT)
    ok = ra > 0 and rb > 0
    rs = r / np.sqrt(ra * rb) if ok else np.nan
    flag = (not ok) or ra <= K * sa or rb <= K * sb
    return dict(donors=len(both), r=r, rel_a=ra, rel_b=rb, r_star=rs, flag=bool(flag),
                r_star_bounded=np.nan if flag else rs, r_star_clip=min(rs, 1.0) if np.isfinite(rs) else np.nan)


def _ref_one(Z, periph, k, seed):
    limit_threads(1)
    rng = np.random.default_rng(seed)
    a = list(rng.choice(periph, k, replace=False))
    return block_stat(Z, a, [t for t in periph if t not in a], rng)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-ref", type=int, default=1000)
    args = ap.parse_args(argv)
    cli.setup(args)
    rng = np.random.default_rng(M.SEED)
    pb, donors = S.cohort("fast")
    d = donors.reindex(pb.donor)
    t = pd.DataFrame({"donor": pb.donor.to_numpy(), "compartment": pb.compartment.to_numpy(), "tod": d.tod.to_numpy(),
                      "hardy": d.hardy.to_numpy()})
    ph14, _ = A103.phase_departures(pb, t, list(S.CLOCK))
    Z = S.zwide(ph14, t.donor.to_numpy(), t.compartment.to_numpy())
    tissues = list(Z.columns)
    brain = [c for c in tissues if c.startswith("Brain")]
    periph = [c for c in tissues if c not in brain]
    targets = {"brain": (brain, periph), "digestive": (list(M.DIGESTIVE), [c for c in periph if c not in M.DIGESTIVE])}
    for name, blk in M.DESCRIPTIVE.items():
        targets[name] = (list(blk), [c for c in periph if c not in blk])
    rows, refs, summ = [], {}, []
    for name, (a, b) in targets.items():
        st = block_stat(Z, a, b, rng)  # same order of rng use as the registered script (boot in between is skipped: values differ in the stream only)
        k = len(a)
        if k not in refs:
            refs[k] = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(_ref_one)(Z, periph, k, M.SEED * 100_000 + 1000 * k + i)
                                                              for i in range(args.n_ref)))
            R = refs[k]
            fin = R.r_star.dropna()
            summ.append(dict(k=k, n_ref=len(R), share_flagged=float(R.flag.mean()), share_rstar_gt1=float((fin > 1).mean()),
                             share_rstar_nan=float(R.r_star.isna().mean()), rstar_median=float(fin.median()),
                             rstar_q05=float(fin.quantile(.05)), bounded_n=int(R.r_star_bounded.notna().sum()),
                             bounded_q05=float(R.r_star_bounded.quantile(.05)), bounded_median=float(R.r_star_bounded.median()),
                             clip_q05=float(R.r_star_clip.quantile(.05)), r_q05=float(R.r.quantile(.05)), r_median=float(R.r.median())))
        R = refs[k]
        p = lambda col, v: float((1 + np.sum(R[col].dropna() <= v)) / (1 + R[col].notna().sum())) if np.isfinite(v) else np.nan  # noqa: E731
        rows.append(dict(block=name, k=k, **st, ref_q05_old=float(R.r_star.dropna().quantile(.05)),
                         ref_q05_bounded=float(R.r_star_bounded.quantile(.05)), p_old=p("r_star", st["r_star"]),
                         p_bounded=p("r_star_bounded", st["r_star_bounded"]), p_clip=p("r_star_clip", st["r_star_clip"]),
                         ref_r_q05=float(R.r.quantile(.05)), ref_r_median=float(R.r.median()), p_uncorrected_r=p("r", st["r"])))
        print({kk: (round(v, 4) if isinstance(v, float) else v) for kk, v in rows[-1].items()}, flush=True)
    res = pd.DataFrame(rows)
    old = pd.read_csv(PAPER / "exploratory/A142_organ_blocks/output/targets.csv")
    old = old[old.plane == "base"].set_index("block")
    res["registered_r_star"] = res.block.map(old.r_star)
    res["registered_p"] = res.block.map(old.p_one_sided)
    res["registered_ref_q05"] = res.block.map(old.ref_q05)
    res.to_csv(OUT / "targets_bounded.csv", index=False)
    pd.DataFrame(summ).to_csv(OUT / "reference_summary_bounded.csv", index=False)
    for k, R in refs.items():
        R.to_csv(OUT / f"reference_k{k}_bounded.csv", index=False)  # per random peripheral block: aggregates of tissues only
    print(res.round(4).T.to_string(), flush=True)
    print(pd.DataFrame(summ).round(4).T.to_string(), flush=True)


if __name__ == "__main__":
    main()
