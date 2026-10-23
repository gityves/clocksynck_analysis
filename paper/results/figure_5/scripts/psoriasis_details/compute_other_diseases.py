"""A53 (H53.2), A55 and A59: variance ratio and phase offset of lesions against uninvolved tissue in psoriasis, atopic
dermatitis, other skin diseases and inflamed gut.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A53, A55, A59). Estimators as A54 (see psoriasis/
compute_replication.py): per series the reference gate on the uninvolved samples, r = log tr Cov(beta_lesion) -
log tr Cov(beta_uninvolved), the lesion-vs-uninvolved offset (h, < 0: lesion delayed), pooled over gated series, per-pair
label swaps (10,000), pseudo-clock control.
- A53 (H53.2): damping, the hypothesis was r < 0 (psoriasis primary, atopic dermatitis secondary): lower-tail swap p and
  pseudo-clock q05. The data gave r > 0; A54 re-pre-registered r > 0 and the offset.
- A55: inflamed vs non-inflamed gut mucosa (six series; CD, UC; same bowel segment), composition genes EPCAM/VIL1/CDH17/
  KRT20 (epithelial) vs PTPRC/CD3E/CD68 (immune).
- A59: other skin diseases (all pooled: H59.1 offset < 0; without mycosis fungoides; each condition; linear-scale series
  descriptive).

    python compute_other_diseases.py [--subsample N] [--only a53,a55,a59]
Outputs: paper/results/figure_5/output/psoriasis_details/a53_skin_{per_series,pooled}.csv,
  a55_ibd[_CD|_UC|_samesegment]_{per_series,pooled,genes}.csv, a59_{all,no_mf,linear_descriptive,<condition>}_*.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a53_skin.py, a55_ibd.py, a59_other.py (+ a54_skin.analyse)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402


def a53(args, run, out):
    samples, means = D.discovery_all(args.subsample, args.seed)
    sz = H.sizes(args, run)
    d, pooled = D.damping(samples, means, 532, sz["n_perm"], sz["n_boot"],
                          n_pseudo=50 if run["quick"] else 200, min_pairs=5)
    pooled.insert(0, "analysis", ["H53_2_primary" if x == "psoriasis" else "H53_2_secondary_AD" for x in pooled.disease])
    d.to_csv(out / "a53_skin_per_series.csv", index=False)
    pooled.to_csv(out / "a53_skin_pooled.csv", index=False)
    print(d.round(3).to_string(index=False), "\n", pooled.round(3).T.to_string(), flush=True)


def analyse(samples, means, tag, out, seed, hyp, run, args, **kw):
    try:
        _, d, v = H.analyse(samples, means, tag, out, seed, hyp, min_pairs=5 if run["quick"] else 10,
                            **H.sizes(args, run), **kw)
    except ValueError as e:  # too few gated series for a pooled estimate
        print(tag, "skipped:", e, flush=True)
        return None
    print(tag, "\n", d.round(3).to_string(index=False), "\n", pd.Series(v).to_string(), flush=True)
    return v


def a55(args, run, out):
    samples, means = D.ibd(args.subsample, args.seed)
    kw = dict(pos=D_EPI, neg=D_IMM)
    analyse(samples, means, "a55_ibd", out, 55, "H55", run, args, **kw)
    for dis in ("CD", "UC"):
        analyse(samples[samples.disease == dis], means, f"a55_ibd_{dis}", out, 55, "H55", run, args, **kw)
    ss = samples[samples.same_segment.astype(bool)]
    analyse(ss, means, "a55_ibd_samesegment", out, 55, "H55", run, args, **kw)


def a59(args, run, out):
    samples, means = D.other_skin(args.subsample, args.seed)
    std = samples[samples.source == "standard"]
    analyse(std, means, "a59_all", out, 59, "H59", run, args)
    mf = std[~std.disease.astype(str).str.contains("mycosis")]
    analyse(mf, means, "a59_no_mf", out, 59, "H59", run, args)
    lin = samples[samples.source == "linear"]
    if len(lin):
        analyse(lin, means, "a59_linear_descriptive", out, 59, "H59", run, args)
    for c in sorted(std.disease.astype(str).unique()):
        analyse(std[std.disease == c], means, f"a59_{c}", out, 59, "H59", run, args)


D_EPI = ["EPCAM", "VIL1", "CDH17", "KRT20"]
D_IMM = ["PTPRC", "CD3E", "CD68"]


def main(argv=None):
    ap = H.parser(__doc__)
    ap.add_argument("--only", default="a53,a55,a59")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    for name, f in (("a53", a53), ("a55", a55), ("a59", a59)):
        if name in args.only.split(","):
            f(args, run, out)


if __name__ == "__main__":
    main()
