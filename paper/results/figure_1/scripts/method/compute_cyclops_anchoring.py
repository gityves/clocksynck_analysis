"""Figure 1: can whole-day samples anchor CYCLOPS on hospital-hours samples? (Li, Hammarlund, Wu et al., PNAS 2024, CYCLOPS 2.0)

Post hoc (legacy `fig1_method`, notebook 4; no amendment). Per GTEx tissue (Muscle, Thyroid, Lung): W = 80 donors who died 08-18 h and
R = 80 other donors (whole day), drawn with one generator (seed 7); five CYCLOPS 2.0 runs: `window` (W alone), `pooled` (W + R, no
batch effect), `covariate` (W + R, dataset as covariate: the design of Li et al.), `covariate_batch` (as `covariate` with a per-gene
log-normal batch factor, sd 0.3, on W), `covariate_batch_timed` (as `covariate_batch` with 15 of the W donors time-stamped during
training; they are left out of the scores). The paper's own case: GTEx breast + TCGA normal breast + 200 TCGA ductal tumours,
covariates dataset and tumour status. Scores on the W donors after the best rotation and reflection (oracle orientation): median error,
skill over the best constant prediction, concentration of the predicted times (true about 0.75), and the clock-gene arrangement
(`evaluate.arrangement_error`) against the timed truth; for the breast case the median error of the GTEx donors and the concentration and
centre of each group's anchored phases. Only aggregates are written (GTEx time of death is controlled access).
Method code: `common/benchmark/` (CYCLOPS in Julia under `$SCCC_DATA/tools`); no clocksync call.

    python compute_cyclops_anchoring.py [--subsample N] [--tissues Muscle - Skeletal ...]
`--subsample N`: check mode (N donors for W and for R, 3 time-stamped, muscle only, no breast case).
Outputs: paper/results/figure_1/output/method/cyclops_anchoring/{scores,genes,breast_summary,breast_genes}.csv
Legacy: legacy/paper_code/fig1_method/{run_cyclops_anchoring.py,4_cyclops_anchoring.ipynb} (outputs/paper/fig1_method/anchoring/)
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import cli  # noqa: E402
from common.benchmark import evaluate as ev  # noqa: E402
from common.benchmark.cyclops import run_cyclops, seeds  # noqa: E402
from common.benchmark.data import cofe_filter, in_window, log2_normalised, tcga_counts, tcga_normals, tcga_tumours  # noqa: E402
from common.benchmark.tcga import COFE_CLOCK  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402

FOLDER = "method"
TISSUES = ("Muscle - Skeletal", "Thyroid", "Lung")
N, N_TIMED, BATCH_SD = 80, 15, 0.3
CONDITIONS = ("window", "pooled", "covariate", "covariate_batch", "covariate_batch_timed")
SEED = 7


def designs(tissues, n, n_timed):
    """{tissue: dict(counts, tod, W, R, timed, shifted, dataset)}: the donor sets and the batch-shifted counts."""
    rng = np.random.default_rng(SEED)
    out = {}
    data = S.load(tissues)
    for tissue in tissues:
        counts, tod = data[tissue]
        W = rng.choice(tod.index[in_window(tod)], n, replace=False)
        R = rng.choice(tod.index.difference(W), n, replace=False)
        timed = rng.choice(W, n_timed, replace=False)
        factor = pd.Series(rng.lognormal(0, BATCH_SD, counts.shape[1]), index=counts.columns)
        shifted = counts.copy()
        shifted.loc[W] = np.round(counts.loc[W] * factor)
        dataset = pd.DataFrame({"dataset": np.where(np.isin(counts.index, W), "W", "R")}, index=counts.index)
        out[tissue] = dict(counts=counts, tod=tod, W=W, R=R, timed=timed, shifted=shifted, dataset=dataset)
    return out


def breast_data():
    """GTEx breast, TCGA normal breast and 200 TCGA tumours on their common genes, with the covariates."""
    g_counts, g_tod = S.load(["Breast - Mammary Tissue"])["Breast - Mammary Tissue"]
    nor, tum = tcga_normals("BRCA"), tcga_tumours("BRCA").sample(200, random_state=0)
    t_counts = tcga_counts("BRCA", pd.concat([nor, tum]))
    common = g_counts.columns.intersection(t_counts.columns)
    both = pd.concat([g_counts[common], t_counts[common]])
    cov = pd.DataFrame({"dataset": ["GTEx"] * len(g_counts) + ["TCGA"] * len(t_counts),
                        "tumour": ["no"] * (len(g_counts) + len(nor)) + ["yes"] * len(tum)}, index=both.index)
    return g_counts, g_tod, both, cov


def cyclops(counts, tmp, ncores, **kw):
    """CYCLOPS run on non-log expression; phases in hours (raw: arbitrary origin and direction)."""
    linear = 2 ** log2_normalised(counts) - 1
    r = run_cyclops(linear.T, seeds(), Path(tmp), seed=0, threads=ncores, **kw)
    return (r.phase_raw * 24 / (2 * np.pi) % 24).rename("phase")


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--tissues", nargs="*", default=None)
    ap.add_argument("--ncores", type=int, default=4)
    args = ap.parse_args()
    opt = cli.setup(args)
    quick = opt["quick"]
    out = out_dir(FOLDER, quick) / "cyclops_anchoring"
    out.mkdir(parents=True, exist_ok=True)
    n, n_timed = (args.subsample, 3) if quick else (N, N_TIMED)
    tissues = tuple(args.tissues) if args.tissues else TISSUES[:1] if quick else TISSUES
    scores, genes = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for tissue, d in designs(tissues, n, n_timed).items():
            counts, tod, W, R, timed = d["counts"], d["tod"], d["W"], d["R"], d["timed"]
            WR, dataset = np.concatenate([W, R]), d["dataset"]
            runs = {"window": cyclops(counts.loc[W], f"{tmp}/w", args.ncores),
                    "pooled": cyclops(counts.loc[WR], f"{tmp}/p", args.ncores),
                    "covariate": cyclops(counts.loc[WR], f"{tmp}/c", args.ncores, covariates=dataset.loc[WR]),
                    "covariate_batch": cyclops(d["shifted"].loc[WR], f"{tmp}/cb", args.ncores, covariates=dataset.loc[WR]),
                    "covariate_batch_timed": cyclops(d["shifted"].loc[WR], f"{tmp}/cbt", args.ncores, covariates=dataset.loc[WR],
                                                     train_times=tod.loc[timed])}
            truth = ev.cosinor(cofe_filter(counts), tod, COFE_CLOCK)
            score_on = np.setdiff1d(W, timed)
            logW = cofe_filter(counts.loc[W])
            for cond, ph in runs.items():
                p, t = ph.loc[score_on], tod.loc[score_on]
                al, mae, _ = ev.align(p, t)
                arr = ev.arrangement_error(ev.cosinor(logW.loc[score_on], al, COFE_CLOCK), truth)
                genes.append(arr.assign(tissue=tissue, condition=cond))
                row = dict(tissue=tissue, condition=cond, mae_h=mae, skill=ev.skill(mae, t), R_true=ev.spread(t),
                           R_pred=ev.spread(p), gene_err_h=float(np.median(np.abs(arr.phase_err))))
                if cond != "window":
                    pr = ph.loc[R]
                    _, mae_r, _ = ev.align(pr, tod.loc[R])
                    row.update(skill_R=ev.skill(mae_r, tod.loc[R]), R_pred_R=ev.spread(pr))
                scores.append(row)
            print(tissue, flush=True)
        pd.DataFrame(scores).to_csv(out / "scores.csv", index=False)
        pd.concat(genes).rename_axis("gene").reset_index().to_csv(out / "genes.csv", index=False)
        print(pd.DataFrame(scores).round(2).to_string(index=False))
        if quick:
            return
        g_counts, g_tod, both, cov = breast_data()
        ph = cyclops(both, f"{tmp}/breast", args.ncores, covariates=cov)
        aligned, mae_g, dr = ev.align(ph.loc[g_counts.index], g_tod)
        off = ev.wrap(aligned - dr * ph.loc[g_counts.index].to_numpy())[0]
        breast = pd.DataFrame({"phase_h": (dr * ph + off) % 24, "group": np.select(
            [cov.dataset == "GTEx", cov.tumour == "no"], ["GTEx breast", "TCGA normal"], "TCGA tumour")}, index=both.index)
        truth_b = ev.cosinor(cofe_filter(g_counts), g_tod, COFE_CLOCK)
        bg_ = []
        for grp, idx in breast.groupby("group").groups.items():
            X = cofe_filter(both.loc[idx])
            est = ev.cosinor(X, breast.phase_h.loc[idx], COFE_CLOCK)
            bg_.append(est.join(ev.arrangement_error(est, truth_b)).assign(group=grp))
        pd.concat(bg_).rename_axis("gene").reset_index().to_csv(out / "breast_genes.csv", index=False)
        summ = breast.groupby("group").phase_h.agg(n="size", R=ev.spread,
                                                   centre_h=lambda h: np.angle(np.exp(1j * ev.W * h).mean()) / ev.W % 24)
        summ["gtex_mae_h"], summ["gtex_skill"] = mae_g, ev.skill(mae_g, g_tod)
        summ.round(3).to_csv(out / "breast_summary.csv")
        print(summ.round(2))


if __name__ == "__main__":
    main()
