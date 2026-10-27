"""A142: is the brain-periphery split specific to the brain, and do digestive organs carry a departure of their own?

Registered in paper/docs/registered_experiments.md (A142). Statistic of A103 base (`results/figure_2/scripts/organs/
compute_brain_artefact.py`, `_shared.block_stat`), generalised from "Brain*" to any tissue set: per-tissue standardised
linear phase departures from the timed fit (GTEx v10 fast deaths, tissues with >= 40 donors), block score = donor mean over
the block's tissues (>= 2), r* = r / sqrt(rel_A rel_B) with split-half reliabilities (50 splits). Each target block is read
against a reference of 1,000 random size-matched blocks of peripheral tissues (each against the rest of the periphery).
H142.1 brain (k = 13), H142.2 digestive (k = 7); one-sided p = (1 + #{reference <= target}) / (1 + n_ref); Holm over the two.
Sensitivities (no PER1/PER2; ischaemic time + Hardy adjusted) and descriptive blocks as registered. Aggregates only.

    python compute_organ_blocks.py [--subsample N] [--n-ref 1000] [--n-boot 200] [--jobs J]
`--subsample N`: check mode, N donors, 50 reference blocks, 10 bootstraps.
Outputs: paper/exploratory/A142_organ_blocks/output/{targets,reference_summary,tests}.csv
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
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

import _shared as S  # noqa: E402
import compute_brain_artefact as A103  # noqa: E402

AMENDMENT = "A142_organ_blocks"
N_SPLIT = A103.N_SPLIT
SEED = 142
DIGESTIVE = ("Stomach", "Colon - Transverse", "Colon - Sigmoid", "Esophagus - Gastroesophageal Junction",
             "Esophagus - Muscularis", "Esophagus - Mucosa", "Liver")
DESCRIPTIVE = {
    "food_exposed": ("Stomach", "Colon - Transverse", "Esophagus - Mucosa", "Liver"),
    "cardiovascular": ("Artery - Aorta", "Artery - Coronary", "Artery - Tibial", "Heart - Atrial Appendage",
                       "Heart - Left Ventricle"),
    "endocrine_adipose": ("Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Adrenal Gland", "Thyroid", "Pituitary"),
}


def block_stat(Z: pd.DataFrame, a: list, b: list, rng) -> dict:
    """`_shared.block_stat` with block A = `a` and block B = `b` (tissue lists) instead of brain vs the rest."""
    both = pd.concat([S._score(Z[a]), S._score(Z[b])], axis=1).dropna()
    r = both.corr().iloc[0, 1]
    ra, rb = S.reliability(Z.loc[both.index, a], rng, N_SPLIT), S.reliability(Z.loc[both.index, b], rng, N_SPLIT)
    return dict(donors=len(both), r=r, rel_a=ra, rel_b=rb, r_star=r / np.sqrt(ra * rb) if ra > 0 and rb > 0 else np.nan)


def boot(Z, a, b, rng, n_boot):
    el = Z.index[(Z[a].notna().sum(1) >= 2) & (Z[b].notna().sum(1) >= 2)]
    Ze = Z.loc[el]
    bs = [block_stat(Ze.iloc[rng.integers(0, len(Ze), len(Ze))].reset_index(drop=True), a, b, rng)["r_star"]
          for _ in range(n_boot)]
    bs = np.where(np.isfinite(bs), bs, np.inf)
    return tuple(np.quantile(bs, [0.025, 0.975])) if n_boot else (np.nan, np.nan)


def _ref_one(Z, periph, k, seed):
    limit_threads(1)
    rng = np.random.default_rng(seed)
    a = list(rng.choice(periph, k, replace=False))
    return block_stat(Z, a, [t for t in periph if t not in a], rng)["r_star"]


def reference(Z, periph, k, n_ref, jobs):
    return np.array(Parallel(n_jobs=jobs)(delayed(_ref_one)(Z, periph, k, SEED * 100_000 + 1000 * k + i)
                                          for i in range(n_ref)), float)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-ref", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    n_ref = args.n_ref or (50 if quick else 1000)
    n_boot = args.n_boot if args.n_boot is not None else (10 if quick else 200)
    out = exploratory_dir(AMENDMENT, quick)
    rng = np.random.default_rng(SEED)

    pb, donors = S.cohort("fast", subsample=args.subsample, seed=args.seed)
    d = donors.reindex(pb.donor)
    t = pd.DataFrame({"donor": pb.donor.to_numpy(), "compartment": pb.compartment.to_numpy(), "tod": d.tod.to_numpy(),
                      "hardy": d.hardy.to_numpy(), "isch": S.sample_ischaemic_h(pb)})
    zw = lambda v: S.zwide(v, t.donor.to_numpy(), t.compartment.to_numpy())  # noqa: E731
    genes14 = list(S.CLOCK)
    genes12 = [g for g in genes14 if g not in A103.ACUTE]
    ph14, _ = A103.phase_departures(pb, t, genes14)
    ph12, _ = A103.phase_departures(pb, t, genes12)
    planes = {"base": zw(ph14), "noPER12": zw(ph12), "adjusted": zw(A103.adjust(ph14, t))}

    tissues = list(planes["base"].columns)
    brain = [c for c in tissues if c.startswith("Brain")]
    periph = [c for c in tissues if c not in brain]
    print(f"{len(planes['base'])} donors, {len(brain)} brain, {len(periph)} peripheral tissues", flush=True)
    targets = {"brain": (brain, periph), "digestive": (list(DIGESTIVE), [c for c in periph if c not in DIGESTIVE])}
    for name, blk in DESCRIPTIVE.items():
        targets[name] = (list(blk), [c for c in periph if c not in blk])
    missing = [c for a, b in targets.values() for c in a + b if c not in tissues]
    assert not missing, f"tissues not in the cohort: {missing}"

    rows, refs, tests = [], [], []
    for plane, Z in planes.items():
        ref_cache = {}
        for name, (a, b) in targets.items():
            k = len(a) if name != "brain" else len(a)
            st = block_stat(Z, a, b, rng)
            lo, hi = boot(Z, a, b, rng, n_boot)
            if k not in ref_cache:
                ref_cache[k] = reference(Z, periph, k, n_ref, args.jobs)
                rr = ref_cache[k][np.isfinite(ref_cache[k])]
                refs.append(dict(plane=plane, k=k, n_ref=len(rr), mean=rr.mean(), sd=rr.std(ddof=1), q05=np.quantile(rr, .05),
                                 q50=np.median(rr), q95=np.quantile(rr, .95), share_below_1=float((rr < 1).mean())))
            rr = ref_cache[k][np.isfinite(ref_cache[k])]
            p = (1 + np.sum(rr <= st["r_star"])) / (1 + len(rr))
            rows.append(dict(plane=plane, block=name, k=k, n_other=len(b), **st, r_star_lo=lo, r_star_hi=hi,
                             ref_q05=np.quantile(rr, .05), ref_median=np.median(rr), p_one_sided=p))
            print({kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in rows[-1].items()}, flush=True)
    res = pd.DataFrame(rows)
    prim = res[(res.plane == "base") & res.block.isin(["brain", "digestive"])].copy()
    prim["hypothesis"] = prim.block.map({"brain": "H142.1", "digestive": "H142.2"})
    order = prim.p_one_sided.sort_values().index
    m, adj, run_max = len(prim), {}, 0.0
    for i, ix in enumerate(order):
        run_max = max(run_max, min(1.0, (m - i) * prim.p_one_sided[ix]))
        adj[ix] = run_max
    prim["p_holm"] = pd.Series(adj)
    prim["verdict"] = np.where(prim.p_holm < 0.05, "pass", "fail")
    res.to_csv(out / "targets.csv", index=False)
    pd.DataFrame(refs).to_csv(out / "reference_summary.csv", index=False)
    prim.to_csv(out / "tests.csv", index=False)
    print(prim[["hypothesis", "r_star", "r_star_lo", "r_star_hi", "ref_q05", "ref_median", "p_one_sided", "p_holm", "verdict"]]
          .round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
