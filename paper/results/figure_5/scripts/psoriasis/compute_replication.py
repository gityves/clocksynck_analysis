"""A54: psoriatic lesions against uninvolved skin of the same visit, held-out replication (H54.1-H54.3), plus the
discovery series.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A54). Estimands, per series: r = log tr Cov(beta_lesional) - log tr
Cov(beta_uninvolved) (Poisson-corrected, one projector for both) and the phase offset arg(ccorr(lesional,
uninvolved)) (h, < 0: lesional delayed). Gate on the uninvolved samples; pooling over gated series (r: fixed-effect
inverse-variance, DerSimonian-Laird reported; offset: weighted resultant with weights 1/se^2 of r). Null: per-pair
swaps of the two labels (10,000). Interval: pair bootstrap SE of r (1,000). H54.1 r > 0, swap p < 0.05 and r above
the pseudo-clock q95; H54.2 offset < 0, one-sided swap p < 0.05; H54.3 leave-one-gene-out offsets all < 0 and >= 10
of 14 per-gene shifts negative. Composition (descriptive): offset in pairs below/above the median keratinocyte-minus-
dermal score difference.

    python compute_replication.py [--subsample N] [--jobs 4]
Outputs: paper/results/figure_5/output/psoriasis/a54_{discovery,heldout}_{per_series,pooled,genes}.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a54_skin.py (+ paired.py)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402

SEED = 54


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(H.FOLDER, run["quick"])
    n = args.subsample
    for which, load in (("heldout", H.heldout), ("discovery", H.discovery)):
        samples, means = load(n, args.seed)
        _, d, v = H.analyse(samples, means, f"a54_{which}", out, SEED, "H54", min_pairs=5 if run["quick"] else 10,
                            **H.sizes(args, run))
        print(which, "\n", d.round(3).to_string(index=False), "\n", pd.Series(v).to_string(), flush=True)


if __name__ == "__main__":
    main()
