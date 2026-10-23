"""A56 H56.1 and A73: is the psoriatic lesion delay larger where the lesion is more active (IL-17) or more
proliferative?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A56 H56.1, A73 H73.1). Gated A54 series (discovery and held-out, from
`compute_replication.py`), weights 1 / se^2 of the A54 variance ratio. Per pair the lesional-minus-uninvolved
difference of a z-scored mean-expression score; pairs split at the within-series median. Estimand: pooled offset of the
high half minus pooled offset of the low half (h, wrapped; < 0: the high half is more delayed). Null: high/low labels
permuted within series (n_perm); one-sided p (lower). H56.1: IL-17 score (IFN-gamma secondary), high - low < 0 and
p < 0.05. H73.1: proliferation score (differentiation loss, basal minus differentiated, secondary).

    python compute_dose.py [--subsample N]      (needs compute_replication.py run with the same --subsample)
Outputs: paper/results/figure_5/output/psoriasis/a56_dose.csv, a56_dose_series.csv, a73_proliferation.csv, a73_series.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a56_dose.py, a73_proliferation.py (the descriptive treatment
part of a73 is in compute_response.py)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402

PD = H.PD


def split_stat(states, w, highs):
    hi = PD.pooled_offset([PD.ccorr(*(x[m] for x in b)) for b, m in zip(states, highs)], w)
    lo = PD.pooled_offset([PD.ccorr(*(x[~m] for x in b)) for b, m in zip(states, highs)], w)
    return H.wrap1(hi - lo), hi, lo


def run(items, weights, scorer, label, rng, n_perm):
    """items: {series: (case, ref, counts)}; scorer(case, ref, counts) -> score difference per pair."""
    states, highs, w, rows = [], [], [], []
    for s, (case, ref, counts) in items.items():
        d = scorer(case, ref, counts)
        states.append(H.states(case, ref, counts))
        highs.append(d > np.median(d))
        w.append(weights[s])
        rows.append(dict(series=s, n=len(d), median_diff=float(np.median(d))))
    w = np.asarray(w)
    stat, hi, lo = split_stat(states, w, highs)
    null = np.array([split_stat(states, w, [rng.permutation(m) for m in highs])[0] for _ in range(n_perm)])
    p = float((1 + np.sum(null <= stat)) / (1 + n_perm))
    return (dict(score=label, n_series=len(rows), n_pairs=int(sum(r["n"] for r in rows)), offset_high_h=hi,
                 offset_low_h=lo, high_minus_low_h=stat, p_lower=p), pd.DataFrame(rows))


def score(genes):
    return lambda case, ref, counts: H.score_diff(case, ref, counts, genes)


def diffloss(case, ref, counts):
    """Basal minus differentiated score difference (A73 secondary)."""
    return H.score_diff(case, ref, counts, H.BASAL) - H.score_diff(case, ref, counts, H.DIFF)


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run_ = cli.setup(args)
    out = H.out_dir(H.FOLDER, run_["quick"])
    n_perm = H.sizes(args, run_, 2000)["n_perm"]
    items, w = H.psoriasis_items(out, args.subsample, args.seed)
    print(len(items), "gated series:", ", ".join(items), flush=True)
    rng = np.random.default_rng(56)
    rows = []
    for genes, label in ((H.IL17, "IL17_primary"), (H.IFNG, "IFNG_secondary")):
        v, d = run(items, w, score(genes), label, rng, n_perm)
        if label == "IL17_primary":
            v["H56_1"] = "PASS" if (v["high_minus_low_h"] < 0 and v["p_lower"] < 0.05) else "FAIL"
            d.to_csv(out / "a56_dose_series.csv", index=False)
        rows.append(v)
        print(v, flush=True)
    pd.DataFrame(rows).to_csv(out / "a56_dose.csv", index=False)
    rng = np.random.default_rng(73)
    rows = []
    v, d = run(items, w, score(H.PROLIF), "proliferation_primary", rng, n_perm)
    v["H73_1"] = "PASS" if (v["high_minus_low_h"] < 0 and v["p_lower"] < 0.05) else "FAIL"
    d.to_csv(out / "a73_series.csv", index=False)
    rows.append(v)
    v2, _ = run(items, w, diffloss, "differentiation_loss_secondary", rng, n_perm)
    rows.append(v2)
    pd.DataFrame(rows).to_csv(out / "a73_proliferation.csv", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
