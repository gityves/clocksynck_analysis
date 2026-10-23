"""A28 readout E4: the extended waveforms on a living hospital cohort with recorded times (Li 2024 breast, V25.5 estimator).

Pre-registered as A28 (2026-10-15, before the evaluation; docs/PLAN_internal_desync.md (archive)): "E4: V25.5 (Li 2024 breast, A25 estimator, hospital
prior). Spearman and median error. Genes absent from the TPM are dropped." Criterion (no harm): the V25.5 Spearman drops by <= 0.10 with
W18 relative to W14. Waveforms W14, W18, Wcons, Wgen are fitted on the GTEx selection half exactly as `compute_extended_clock.py` (aggregates only),
the original `gtex()` waveform is reported for reference. Cohort and estimator as `windows/compute_breast_timed.py` (two
compartments, normal and tumour, form the host; the frame shift is read from their cross-covariance, local: not in clocksync). Estimands: Spearman of
the tangent coordinate with the recorded hour (and its permutation p), median absolute error after the common offset under the hospital
prior (and its permutation p; 2,000 shuffles).

    python compute_extended_e4.py [--subsample N] [--seed S] [--n-perm P]
Outputs (paper/results/extended_data/output/gate/): extended_clock_e4.csv
Legacy: legacy/desync/validation/extended_clock/scripts/evaluate.py (e4) -> output/e4_li.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS, group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("windows")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import compute_extended_clock as X  # noqa: E402
import compute_breast_timed as BT  # noqa: E402

FOLDER = "gate"


def breast_counts(genes, n_reads=30e6):
    """local: not in clocksync. The V25.5 pseudobulk of `genes` (counts approximated from TPM x length, rescaled to `n_reads`;
    genes absent from the table are dropped), as the builder `breast_timed_gse233242` does for the 14 clock genes."""
    sys.path.insert(0, str(DATASETS / "breast_timed_gse233242"))
    import importlib.util
    spec = importlib.util.spec_from_file_location("bt_build", DATASETS / "breast_timed_gse233242" / "build.py")
    b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b)
    t = pd.read_csv(b.TRUTH)
    t = t[t.note.fillna("") == ""]
    key = pd.concat([pd.DataFrame(dict(patient=t.patient_id, compartment=c, SL=t[f"tpm_column_{c}"])) for c in BT.COMPS],
                    ignore_index=True)
    x = pd.read_csv(b.TPM).drop_duplicates("gene_name").set_index("gene_name")
    have = [g for g in genes if g in x.index]
    rl = x[key.SL].mul(x.transcriptLen, axis=0)
    cnt = (rl / rl.sum() * n_reads).loc[have].T.round()
    cnt.index = range(len(key))
    pb = cs.from_bulk(cnt, key, donor="patient", compartment="compartment", depth=pd.Series(float(n_reads), index=cnt.index),
                      genes=have)
    donors = t.set_index("patient_id")[["resection_h"]]
    return pb, donors, have


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    limit_threads(1)
    out = out_dir(FOLDER, run["quick"])
    n_perm = args.n_perm or (200 if run["quick"] else 2000)
    df = X.load(args)
    wfs = X.waveforms(df)
    rows = []
    for arm, wf in wfs.items():
        pb, donors, have = breast_counts(list(wf.genes))
        if args.subsample:
            keep = np.random.default_rng(args.seed).permutation(sorted(donors.index))[: args.subsample]
            pb, donors = pb[pb.donor.isin(keep)].reset_index(drop=True), donors.loc[keep]
        row = BT.run(pb, donors, waveform=wf.subset(have), n_boot=0, n_perm=n_perm, seed=args.seed)[0]
        rows.append(dict(arm=arm, genes=len(have), missing=";".join(sorted(set(wf.genes) - set(have))), spearman=row["spearman_t"],
                         p_spearman=row["p_spearman"], mae_day_h=row["mae_day"], p_mae=row["p_mae_day"]))
    e4 = pd.DataFrame(rows)
    base = e4.set_index("arm").spearman
    e4["spearman_change_vs_W14"] = e4.spearman - base["W14"]
    e4["no_harm"] = e4.spearman_change_vs_W14 >= -0.10
    e4.to_csv(out / "extended_clock_e4.csv", index=False)
    print(e4.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
