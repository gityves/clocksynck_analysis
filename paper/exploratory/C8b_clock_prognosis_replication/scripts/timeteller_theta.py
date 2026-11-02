"""A194 readout (registered in `paper/docs/registered_experiments.md`): per-tumour TimeTeller Theta, no outcome data touched.

TimeTeller (`common.benchmark.timeteller`) is trained on GTEx Breast - Mammary Tissue fast deaths (Hardy 1-2) with recorded
time of death (3 h bins), on the published human panel (ARNTL, NPAS2, PER1-3, NR1D1/2, CIART, TEF, DBP), log2 normalised counts
(`common.benchmark.data.log2_normalised`, computed on all genes, then the panel kept). It is applied to SCAN-B (all non-replicate
tumours, log2(FPKM + 0.1)), TCGA BRCA (one primary per patient, recount3 log2 CPM) and METABRIC (log2 intensity, cross-platform,
secondary). TimeTeller z-scores each sample's gene vector (intergene normalisation), so the cohort scales only enter through the
shape of the panel profile. Two variants: "raw" (registered: tumour log2 values as they are) and "aligned" (deviation, see
below: each gene shifted so that the cohort mean equals the training mean; cohorts differ from GTEx by constant per-gene
offsets of platform, gene length and depth scaling, up to 10 log2 units, and the raw tumours fall outside TimeTeller's training
manifold: Theta saturates at 1 and the likelihood is flat). Check: 5-fold out-of-fold time error on the training tissue (aggregate only).
GTEx controlled access: per-donor times only in a temporary directory deleted at the end; the repo gets aggregates.
Also written: Spearman of Theta with the G2M `prolif` score (aggregate).

    python timeteller_theta.py [--cores 4]
"""
import argparse
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[4]
EXP = ROOT / "paper" / "exploratory"
sys.path[:0] = [str(ROOT), str(ROOT / "paper"), str(EXP / "A127_scanb" / "scripts"), str(EXP / "A166_metabric" / "scripts")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import a127_lib as A  # noqa: E402
import a166_lib as L  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import log2_normalised  # noqa: E402
from common.benchmark.timeteller import BIN_H, MIN_PER_BIN, cross_fit, run_timeteller  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
PANEL = ["ARNTL", "NPAS2", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "CIART", "TEF", "DBP"]
TISSUE = "Breast - Mammary Tissue"
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:5.1f} min]", *a, flush=True)


def scanb():
    """All non-replicate SCAN-B tumours (one per patient); outcome columns are not read."""
    head = pd.read_csv(A.OUT / "scanb_samples_parsed.csv", nrows=0).columns[0]
    meta = pd.read_csv(A.OUT / "scanb_samples_parsed.csv", index_col=0, usecols=[head, "is_repl", "col_order", "patient"])
    E = pd.read_parquet(A.OUT / "scanb_genes_kept.parquet")
    m = meta[~meta.is_repl.astype(bool)].sort_values("col_order").drop_duplicates("patient")
    X = E.loc[m.index]
    g2m = [g for g in A.sh.hallmark_g2m() if g in X.columns]
    return X[PANEL], pd.Series(A.sh._zmean(X, g2m).to_numpy(), index=X.index)


def tcga():
    m = pd.read_parquet(L.OUT / "tcga_meta.parquet", columns=["prolif"])
    return pd.read_parquet(L.OUT / "tcga_expr.parquet", columns=PANEL).loc[m.index], m.prolif


def metabric():
    m = pd.read_parquet(L.OUT / "metabric_meta.parquet", columns=["prolif"])
    return pd.read_parquet(L.OUT / "metabric_expr.parquet", columns=PANEL).loc[m.index], m.prolif


def circ_err(pred, true):
    d = np.abs(pred - true) % 24
    return np.minimum(d, 24 - d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cores", type=int, default=4)
    args = ap.parse_args()
    tmp = DATA / "tmp" / "timeteller_a194"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    try:
        pb, donors = gtex.load(tissues=[TISSUE], hardy=(1, 2), genes="all", depth="sum")
        meta = {"donor", "compartment", "depth"}
        c = pb.set_index("donor")[[x for x in pb.columns if x not in meta]]
        if c.columns.duplicated().any():
            c = c.T.groupby(level=0).sum().T
        tod = donors.tod.reindex(c.index).dropna()
        c = c.loc[tod.index]
        X = log2_normalised(c)[PANEL]
        log("train", X.shape)
        # out-of-fold check (aggregate only)
        oof = cross_fit(X, tod, tmp / "oof", cores=args.cores)
        err = circ_err(oof.time_pred.to_numpy(), tod.loc[oof.index].to_numpy())
        tb = (np.floor(tod % 24 / BIN_H)).astype(int)
        used = tb.map(tb.value_counts()) >= MIN_PER_BIN
        pd.Series(dict(n_donors=len(X), bin_h=BIN_H, n_bins_used=int(tb[used].nunique()), min_per_bin=MIN_PER_BIN,
                       oof_median_abs_err_h=round(float(np.median(err)), 2), oof_mean_abs_err_h=round(float(err.mean()), 2),
                       oof_theta_median=round(float(oof.theta.median()), 3), genes=",".join(PANEL))
                  ).to_csv(OUT / "timeteller_training_check.csv", header=["value"], index_label="metric")
        log("oof median abs error h", np.median(err))
        # tumours
        sets = {"SCAN-B": scanb(), "TCGA-BRCA": tcga(), "METABRIC": metabric()}
        test, prolif = {}, {}
        for k, (x, p) in sets.items():
            for v, y in (("raw", x), ("aligned", x - x.mean() + X.mean())):
                test[k, v] = y.astype(float).set_axis([f"{k}|{v}|{i}" for i in y.index.astype(str)])
            prolif[k] = p
        log("test", {k: len(v[0]) for k, v in sets.items()})
        allt = pd.concat(test.values())
        assert allt.notna().all().all()
        res = run_timeteller(X, tod, allt, tmp / "test", cores=args.cores)
        rows = []
        for (k, v), y in test.items():
            r = res.loc[y.index]
            rows.append(pd.DataFrame({"cohort": k, "variant": v, "sample": [i.split("|", 2)[2] for i in y.index],
                                      "theta": r.theta.to_numpy(), "time_pred": r.time_pred.to_numpy(),
                                      "max_loglik": r.max_loglik.to_numpy(), "prolif": prolif[k].to_numpy()}))
        D = pd.concat(rows, ignore_index=True)
        D.drop(columns="prolif").to_csv(OUT / "timeteller_theta.csv", index=False)
        summ = []
        for (k, v), g in D.groupby(["cohort", "variant"], sort=False):
            q = g.theta.quantile([0, .05, .25, .5, .75, .95, 1])
            ok = g.theta.notna() & g.max_loglik.notna()
            ok_t = g.theta.notna()
            summ.append(dict(cohort=k, variant=v, n=len(g), n_theta=int(ok_t.sum()), **{f"q{int(p * 100)}": round(v, 4) for p, v in q.items()},
                             spearman_theta_prolif=round(spearmanr(g.theta[ok_t], g.prolif[ok_t])[0], 3),
                             n_loglik=int(ok.sum()), n_theta_lt1=int((g.theta < 1).sum()),
                             spearman_loglik_prolif=round(spearmanr(g.max_loglik[ok], g.prolif[ok])[0], 3) if ok.sum() > 3 else np.nan))
        pd.DataFrame(summ).to_csv(OUT / "timeteller_theta_summary.csv", index=False)
        log(pd.DataFrame(summ).to_string())
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
