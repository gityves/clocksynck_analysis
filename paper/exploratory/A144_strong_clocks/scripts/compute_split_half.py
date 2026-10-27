"""A144 part 1: is the A104 amplitude-coupling association noise shared within the same donors? (GTEx split halves)

Registered in paper/docs/registered_experiments.md (A144). A104 design (`results/figure_2/scripts/backbone/
compute_oscillator_gtex.py`): GTEx fast deaths, the 27 peripheral tissues with >= 40 donors, all tissues anchors, leave-one-out
host. 100 random donor halves (seed 144): per half and tissue c, A = own amplitude, s = z2 V (locked), u = 2 psi (unlocked).
H144.1: cross-half rho = mean of Spearman(A_H1, c_H2) and Spearman(A_H2, c_H1), averaged over splits, > 0; null: tissue labels
of A permuted against c (the same permutation in every split and half), 10,000; one-sided p. D144.1: same-half rho, the
difference, SD over tissues of log s and log u and their Spearman with c (full sample). Aggregates only (GTEx controlled).

    python compute_split_half.py [--subsample N] [--jobs J]
Outputs: paper/exploratory/A144_strong_clocks/output/{h144_1,d144_1_components,splits}.csv
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

_spec = importlib.util.spec_from_file_location("a104", PAPER / "results/figure_2/scripts/backbone/compute_oscillator_gtex.py")
A104 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(A104)

AMENDMENT = "A144_strong_clocks"
SEED, N_SPLITS, N_PERM = 144, 100, 10_000


def fit(pb: pd.DataFrame, comps: list) -> pd.DataFrame:
    t = A104.per_tissue(cs.Dataset(pb, anchors=comps, name="gtex_periphery")).reindex(comps)
    return t.assign(s=t.z2 * t.V, u=2 * t.psi)


def _split(pb, comps, donors, seed):
    limit_threads(1)
    rng = np.random.default_rng(seed)
    h = rng.permutation(donors)
    halves = (h[: len(h) // 2], h[len(h) // 2:])
    try:
        return [fit(pb[pb.donor.isin(d)], comps) for d in halves]
    except (ValueError, np.linalg.LinAlgError):
        return None


def rank_pairs(t1: pd.DataFrame, t2: pd.DataFrame):
    """Ranks of (A_H1, c_H2), (A_H2, c_H1), (A_H1, c_H1), (A_H2, c_H2) over the tissues estimated in both halves."""
    ok = t1[["own_amplitude", "coupling"]].notna().all(1) & t2[["own_amplitude", "coupling"]].notna().all(1)
    r = {k: rankdata(v[ok.to_numpy()]) for k, v in (("a1", t1.own_amplitude.to_numpy()), ("c1", t1.coupling.to_numpy()),
                                                      ("a2", t2.own_amplitude.to_numpy()), ("c2", t2.coupling.to_numpy()))}
    return ok.to_numpy(), r


def spear(x, y):
    return np.corrcoef(x, y)[0, 1]  # inputs are ranks


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    n_splits, n_perm = (10, 200) if quick else (N_SPLITS, N_PERM)
    out = exploratory_dir(AMENDMENT, quick)
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS)
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[:args.subsample]
        pb = pb[pb.donor.isin(keep)]
    comps = sorted(c for c in pb.compartment.unique() if not c.startswith("Brain"))
    pb = pb[pb.compartment.isin(comps)]
    donors = np.array(sorted(pb.donor.unique()))
    print(f"{len(donors)} donors, {len(comps)} peripheral tissues", flush=True)

    full = fit(pb, comps)
    comp = dict(tissues=int(full[["own_amplitude", "coupling"]].notna().all(1).sum()),
                rho_full=spearmanr(full.own_amplitude, full.coupling, nan_policy="omit").statistic,
                sd_log_s=np.log(full.s).std(ddof=1), sd_log_u=np.log(full.u).std(ddof=1),
                rho_log_s_c=spearmanr(np.log(full.s), full.coupling, nan_policy="omit").statistic,
                rho_log_u_c=spearmanr(np.log(full.u), full.coupling, nan_policy="omit").statistic,
                rho_log_s_log_u=spearmanr(np.log(full.s), np.log(full.u), nan_policy="omit").statistic)
    pd.DataFrame([comp]).to_csv(out / "d144_1_components.csv", index=False)
    print(pd.DataFrame([comp]).round(3).to_string(index=False), flush=True)

    res = [r for r in Parallel(n_jobs=args.jobs)(delayed(_split)(pb, comps, donors, SEED * 1000 + i) for i in range(n_splits))
           if r is not None]
    rows, ranked = [], []
    for i, (t1, t2) in enumerate(res):
        ok, r = rank_pairs(t1, t2)
        ranked.append(r)
        rows.append(dict(split=i, tissues=int(ok.sum()), cross=(spear(r["a1"], r["c2"]) + spear(r["a2"], r["c1"])) / 2,
                         same=(spear(r["a1"], r["c1"]) + spear(r["a2"], r["c2"])) / 2))
    sp = pd.DataFrame(rows)
    sp.to_csv(out / "splits.csv", index=False)
    obs = sp.cross.mean()
    # null: one permutation of the tissue positions of A, applied in every split (tissue sets equal across splits when no
    # tissue is lost; otherwise the permutation is taken over each split's own tissues with the same generator state)
    rng = np.random.default_rng(SEED)
    null = np.empty(n_perm)
    for b in range(n_perm):
        seed_b = rng.integers(2**32)
        vals = []
        for r in ranked:
            p = np.random.default_rng(seed_b).permutation(len(r["a1"]))
            vals.append((spear(r["a1"][p], r["c2"]) + spear(r["a2"][p], r["c1"])) / 2)
        null[b] = np.mean(vals)
    p1 = (1 + np.sum(null >= obs)) / (1 + n_perm)
    h1 = dict(n_splits=len(sp), tissues_median=float(sp.tissues.median()), cross_rho=obs, cross_rho_q05=sp.cross.quantile(.05),
              cross_rho_q95=sp.cross.quantile(.95), same_rho=sp.same.mean(), same_minus_cross=(sp.same - sp.cross).mean(),
              same_minus_cross_q05=(sp.same - sp.cross).quantile(.05), same_minus_cross_q95=(sp.same - sp.cross).quantile(.95),
              null_mean=null.mean(), null_sd=null.std(ddof=1), null_q95=np.quantile(null, .95), p_one_sided=p1,
              verdict="pass" if p1 < 0.05 else "fail")
    pd.DataFrame([h1]).to_csv(out / "h144_1.csv", index=False)
    print(pd.DataFrame([h1]).round(4).T.to_string(), flush=True)


if __name__ == "__main__":
    main()
