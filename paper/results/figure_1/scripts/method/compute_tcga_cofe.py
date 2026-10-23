"""COFE on TCGA as published (Ananthasubramaniam et al., PLoS Biol 2025; S4 Text): the comparison method on tumour bulk, figure 1.

Post hoc (legacy `fig1_method`, no amendment). 11 adenocarcinomas (`TCGA_TYPES`), one histological type each, one primary tumour per
patient, from the same recount3 data; protein-coding genes with >= 64 reads per sample and mean log2 normalised expression >= 7
(`cofe_filter`); COFE with 3-repeated 5-fold cross-validation and 5 restarts (`run_cofe`, seed 0). Per cancer: the sample phases
(orderings are arbitrary in origin and direction), the cross-validation curve, the clock genes' log2 expression, sample and gene numbers,
best sparsity, whether the cross-validation minimum is interior (COFE's validity signal) and the run time. The orderings are not scored
against any time (TCGA has none): they are the inputs of the eigen-core clock comparison and of the COFE-versus-`clocksync` figure.
Method code: `common/benchmark/` (paper-only, ported unchanged from `legacy/src/sccc/benchmark/`); no clocksync call.

    python compute_tcga_cofe.py [--subsample N] [--cancers BLCA ...]     (N: random tumours of the first cancer only)
Outputs: paper/results/figure_1/output/method/tcga_cofe/<cancer>_{phase,cv,clock_expr}.csv, tcga_cofe_meta.csv
Legacy: legacy/paper_code/fig1_method/run_tcga_cofe.py (outputs/paper/fig1_method/tcga/)
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402

from common import cli  # noqa: E402
from common.benchmark.cofe import run_cofe  # noqa: E402
from common.benchmark.data import TCGA_TYPES, cofe_filter, tcga_counts, tcga_tumours  # noqa: E402
from common.benchmark.tcga import COFE_CLOCK  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "method"


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--cancers", nargs="*", default=None)
    ap.add_argument("--ncores", type=int, default=4)
    args = ap.parse_args()
    opt = cli.setup(args)
    out = out_dir(FOLDER, opt["quick"]) / "tcga_cofe"
    out.mkdir(parents=True, exist_ok=True)
    cancers = args.cancers or (["BLCA"] if opt["quick"] else list(TCGA_TYPES))
    rows = []
    for p in cancers:
        s = tcga_tumours(p)
        if opt["quick"]:
            s = s.sample(min(args.subsample, len(s)), random_state=args.seed)
        X = cofe_filter(tcga_counts(p, s))
        t0 = time.time()
        r = run_cofe(X, seed=0, ncores=args.ncores)
        r["phase"].rename("phase").to_csv(out / f"{p}_phase.csv")
        r["cv"].to_csv(out / f"{p}_cv.csv", index=False)
        X[[g for g in COFE_CLOCK if g in X.columns]].to_csv(out / f"{p}_clock_expr.csv")
        rows.append(dict(cancer=p, n_samples=X.shape[0], n_genes=X.shape[1], best_s=r["best_s"], interior_min=bool(r["interior_min"]),
                         runtime_s=round(time.time() - t0)))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(out.parent / "tcga_cofe_meta.csv", index=False)


if __name__ == "__main__":
    main()
