"""A27 part 2, D8 part B: the glucocorticoid score in Li 2024 breast cancer with recorded resection times (post hoc, exploratory).

Post hoc, exploratory (docs/PLAN_internal_desync.md (archive), A27 part 2, D8; no hypothesis ids). Li et al. 2024 (GSE233242): tumour and adjacent
normal bulk RNA-seq of 42 patients with a recorded hospital resection time (08-18 h). GC score per library = mean over FKBP5,
TSC22D3, ZBTB16, KLF9, KLF15, DDIT4 of log(TPM + 1) z-scored across the patients of the same compartment; per patient the mean of
normal and tumour. Estimand: Spearman correlation of the GC score with the recorded time, the internal time (V25.5 per-patient
estimate with the hospital-hours prior, `est_day`) and the internal-time error (circular difference est_day - recorded time), and the
normal-tumour correlation of the score. Reads the per-patient internal times of windows (`v25_5_patients.csv`, run
`windows/compute_breast_timed.py` first; the file holds recorded clinical hours, not controlled data).
D8 part A (21 single-cell cancer cohorts: GC against the A25 per-patient internal time) needs the per-patient A25 times and the
cohort panel pseudobulks of `tme_time.py` / `panel_pb.py`, which no folder of the refactor produces: not ported.

    python compute_gc_li_breast.py [--subsample N]
Outputs (paper/results/figure_3/output/gtex_cues/): diag8_li_gc.csv
Legacy: legacy/desync/validation/glucocorticoid/scripts/diag8_gc_cancer.py (part B)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from common import cli  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

import _genes as G  # noqa: E402
from _shared import GC, W  # noqa: E402

TPM = DATA / "validation" / "GSE233242" / "GSE233242_Processed_file_tpm.csv.gz"
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS  # noqa: E402
TRUTH = DATASETS / "breast_timed_gse233242" / "resection_times.csv"


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(G.FOLDER, run["quick"])
    t = pd.read_csv(TRUTH)
    t = t[t.note.fillna("") == ""]
    x = pd.read_csv(TPM).drop_duplicates("gene_name").set_index("gene_name")
    genes = [g for g in GC if g in x.index]
    sc = {}
    for comp in ("normal", "tumour"):
        L = np.log(x.loc[genes, t[f"tpm_column_{comp}"]].T.to_numpy() + 1)
        sc[comp] = pd.Series(((L - L.mean(0)) / L.std(0)).mean(1), index=t.patient_id.to_numpy())
    li_file = out_dir("windows", False) / "v25_5_patients.csv"
    li = pd.read_csv(li_file).set_index("patient")
    b = pd.DataFrame(dict(gc_normal=sc["normal"], gc_tumour=sc["tumour"])).join(li, how="inner")
    b["gc"] = b[["gc_normal", "gc_tumour"]].mean(axis=1)
    b["err"] = np.angle(np.exp(1j * W * (b.est_day - b.true_h))) / W
    rows = [dict(quantity="GC normal vs tumour", n=len(b), rho=stats.spearmanr(b.gc_normal, b.gc_tumour)[0],
                 p=stats.spearmanr(b.gc_normal, b.gc_tumour)[1])]
    for c, lab in (("true_h", "recorded resection time"), ("est_day", "internal time (hospital prior)"),
                   ("err", "internal-time error")):
        r, p = stats.spearmanr(b.gc, b[c])
        rows.append(dict(quantity=f"GC vs {lab}", n=len(b), rho=r, p=p))
    res = pd.DataFrame(rows)
    res.to_csv(out / "diag8_li_gc.csv", index=False)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
