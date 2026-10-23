"""A27 part 1 / A31b / A50: timed per-gene truths of GTEx tissues (the answer key of every gene-rhythm test).

Post hoc (A27 part 1). Per tissue and gene, cosinor of log(CPM + 0.5) on the recorded time of death (v1: fast deaths;
v2: + RIN > 6 and covariates ischaemic time, sex, age bracket, Hardy; v3: CHIRAL internal phase (DIP) instead of the time of
death, all death types), genes with median CPM >= 5; 2-df F p and Benjamini-Hochberg q. The tissues are those matched to the
A27 and A31b atlases and the A50 cortex truths (BA9, Cortex). Aggregates per gene and tissue only (GTEx controlled access).

    python compute_gtex_truth.py [--subsample N] [--versions v1 v2 v3]
Outputs: paper/results/figure_2/output/gene_rhythms/gtex_cosinor{,_v2,_v3}.csv.gz, gtex_cosinor_summary.csv
Legacy: legacy/desync/validation/gonogo/scripts/gtex_truth{,_v2,_v3}.py, gonogo_new/scripts/gtex_truth_new.py,
gene_map/scripts/gtex_truth_gm.py (`--subsample N` uses N random samples per tissue).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402

import _shared as S  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "gene_rhythms"
SUFFIX = {"v1": "", "v2": "_v2", "v3": "_v3"}


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--versions", nargs="+", default=["v1", "v2", "v3"])
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    summ = []
    for v in args.versions:
        tabs = []
        for tis in S.TRUTH_TISSUES + S.BRAIN_TISSUES:
            t = S.gtex_truth(tis, v, n=args.subsample, seed=args.seed)
            tabs.append(t)
            summ.append(dict(version=v, tissue=tis, n=int(t.n.iloc[0]), genes=len(t), p01=int((t.p < 0.01).sum()),
                             q05=int((t.q < 0.05).sum())))
            print(v, tis, summ[-1], flush=True)
        pd.concat(tabs).to_csv(out / f"gtex_cosinor{SUFFIX[v]}.csv.gz", index=False)
    pd.DataFrame(summ).to_csv(out / "gtex_cosinor_summary.csv", index=False)


if __name__ == "__main__":
    main()
