"""A87: why is the tumour clock less observable? Amplitude x spread, off-plane variance, gene coherence, alignment.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A87; descriptive, no single primary). Components of the case against
the paired reference (`clocksync.decompose`), unit bootstrap (1,000) 95 % intervals; a component contributes if its
interval excludes no change in the direction of loss (reading rule):
  d_log_plane  log Poisson-corrected clock-plane variance, case minus reference (amplitude x spread);
  d_log_off    log off-plane variance of the 14 centred clock genes;
  d_coherence  gene coherence R (leave-gene-out loadings against the waveform phases);
  alignment    TCGA: |complex correlation| of the tumour state with the same patient's normal;
               single cell: |ccorr| of malignant and of normal epithelium with the patient's stromal host (d_align).
(a) TCGA, the ten paired cancers (per cancer, then the mean over cancers).
(b) Single-cell cohorts: the patients with malignant and adjacent-normal epithelium and a stromal host (fibroblast,
    endothelial, mural pseudobulks, both sites; mean state per patient), centred within study x stratum, pooled.

    python compute_decomposition.py [--subsample N] [--seed S]
Outputs (paper/results/figure_5/output/cancer/): a87_tcga.csv, a87_tcga_mean.csv, a87_single_cell.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a87_decomposition.py. Names: alignment = legacy align_corr.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer"
SEED = 87


def single_cell_samples(pb: pd.DataFrame):
    """(malignant and normal epithelium rows of the paired patients, stromal host rows)."""
    b = sh.sc_builder()
    pts = b.paired_epithelium(pb)
    cls = pb.compartment.str.split(":").str[1]
    ok = pb.donor.isin(pts)
    cond = pb.compartment.map({sh.MAL: "malignant", sh.NORM: "normal"})
    samples = pb[ok & cond.notna()].assign(condition=cond[ok & cond.notna()])
    host = pb[ok & cls.isin(sh.STROMA)]
    return samples, host


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = 30 if run["quick"] else 1000
    seed = SEED + args.seed

    pb, _ = sh.load_sc(args.subsample, args.seed)
    samples, host = single_cell_samples(pb)
    s = cs.decompose(samples, unit="donor", condition="condition", levels=("normal", "malignant"), strata="stratum",
                     host=host, n_boot=n_boot, seed=seed)
    s = s.assign(patients=int(s.units.iloc[0]), studies=samples.study.nunique()).drop(columns="units")
    s.to_csv(out / "a87_single_cell.csv")
    print(s.round(3).to_string(), flush=True)

    tc = sh.load_tcga(args.subsample, args.seed)
    rows = []
    for c in sh.CANCERS:
        r = cs.decompose(tc[tc.cohort == c], unit="unit", condition="condition", levels=sh.LEVELS, n_boot=n_boot,
                         seed=seed)
        r = r.rename(index={"alignment": "align_corr"}).assign(cancer=c, n=int(r.units.iloc[0])).drop(columns="units")
        rows.append(r)
        print(c, r.round(3).to_string(), flush=True)
    d = pd.concat(rows)
    d.to_csv(out / "a87_tcga.csv")
    d.groupby(level=0).estimate.mean().to_csv(out / "a87_tcga_mean.csv", header=False)


if __name__ == "__main__":
    main()
