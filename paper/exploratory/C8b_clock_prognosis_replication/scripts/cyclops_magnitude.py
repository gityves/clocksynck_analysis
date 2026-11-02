"""A193 (C8b) readout: per-tumour CYCLOPS 2.0 magnitude in luminal A breast cancer (registered in
`paper/docs/registered_experiments.md`, A193). Readout only: no survival or outcome column is read here.

Cohorts (CYCLOPS human seed list, template hyperparameters, 80 models, three fits per cohort = sample-order seeds 0-2):
  tcga      TCGA BRCA recount3 CPM (non-log); PanCanAtlas BRCA_LumA primary tumours (one per patient, lowest barcode) fitted
            together with adjacent normal breast (sample type 11, one per patient); discrete covariate = tumour status.
  metabric  PAM50 LumA tumours, raw illumina intensities (2**log2), duplicate probes averaged, NaN -> gene mean as in A166;
            all genes (the cached A166 matrix lacks 8 seed genes); covariate = METABRIC cohort.
  scanb     PAM50 LumA tumours, FPKM = 2**x - 0.1 clipped at 0 (`a127_lib.load`, 17.9k genes, all seeds present);
            covariate = platform.
Output (../output/): cyclops_magnitude.csv (cohort, sample, is_tumour, magnitude_s{0,1,2}, phase_s{0,1,2}, magnitude = mean of
within-fit z-scored magnitudes, z over the luminal A tumours of that fit), cyclops_runs.csv, work/<cohort>_s<seed>/ (CYCLOPS
files; a finished fit is reused on rerun). `--check` fits a stratified subsample with train_n_models=2 (outputs suffixed _check).

    python cyclops_magnitude.py [--check] [--threads 8] [--cohorts tcga metabric scanb]
"""
import argparse
import json
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[4]
EXP = ROOT / "paper" / "exploratory"
sys.path[:0] = [str(ROOT), str(ROOT / "paper"), str(EXP / "A166_metabric" / "scripts"), str(EXP / "A127_scanb" / "scripts")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import a127_lib as S  # noqa: E402
from common.benchmark import cyclops as C  # noqa: E402
from common.benchmark.data import RECOUNT3  # noqa: E402,F401  (path check)
from common.datasets.tcga_paired import build as TB  # noqa: E402

OUT = HERE.parents[1] / "output"
SCCC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
SEEDS = (0, 1, 2)
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:6.1f} min]", *a, flush=True)


def load_tcga():
    """(expr genes x samples CPM, meta: is_tumour, covariate). Samples = barcodes."""
    st = pd.read_csv(SCCC / "tcga/pancan_subtypes/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", sep="\t",
                     comment="#", usecols=["PATIENT_ID", "SUBTYPE"]).drop_duplicates("PATIENT_ID").set_index("PATIENT_ID").SUBTYPE
    s = TB.samples_of("BRCA")
    t = s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")
    t = t[t.patient.map(st) == "BRCA_LumA"]
    n = s[s.sample_type == "Solid Tissue Normal"].sort_values("barcode").drop_duplicates("patient")
    meta_s = pd.concat([t, n], ignore_index=True)
    ct = TB.counts_of("BRCA", meta_s)  # samples x gene symbols, scaled counts, index barcode
    ct = ct.loc[:, ct.sum() > 0]
    cpm = ct.div(ct.sum(axis=1), axis=0) * 1e6
    meta = pd.DataFrame({"sample": meta_s.barcode.to_numpy(), "is_tumour": [True] * len(t) + [False] * len(n)}).set_index("sample")
    meta["cov"] = np.where(meta.is_tumour, "tumour", "normal")
    return cpm.T, meta


def load_metabric():
    import a166_lib as L
    meta, _ = L.load("metabric")
    P = SCCC / "metabric/brca_metabric"
    ex = pd.read_csv(P / "data_mrna_illumina_microarray.txt", sep="\t", index_col=0).drop(columns="Entrez_Gene_Id")
    X = ex.groupby(level=0).mean().T.astype(np.float32)
    X = X.fillna(X.mean())
    m = meta[meta.pam50 == "LumA"]
    meta = pd.DataFrame({"is_tumour": True, "cov": m.cohort.astype(int).astype(str)}, index=m.index).rename_axis("sample")
    return (2.0 ** X.loc[meta.index].astype(float)).T, meta


def load_scanb():
    m, X = S.load()
    m = m[m.pam50 == "LumA"]
    meta = pd.DataFrame({"is_tumour": True, "cov": m.platform.astype(str)}, index=m.index).rename_axis("sample")
    return (2.0 ** X.loc[meta.index].astype(float) - 0.1).clip(lower=0).T, meta


LOADERS = {"tcga": load_tcga, "metabric": load_metabric, "scanb": load_scanb}


def fit_one(name, expr, meta, seed, threads, check):
    d = OUT / "work" / f"{name}_s{seed}{'_check' if check else ''}"
    fo = d / "fit_output.csv"
    if not fo.exists():
        log(f"{name} seed {seed}: fit start, {expr.shape[1]} samples x {expr.shape[0]} genes")
        hyper = {"train_n_models": 2} if check else {}
        C.run_cyclops(expr, C.seeds(), d, seed=seed, threads=threads, covariates=meta[["cov"]].rename(columns={"cov": "covariate"}),
                      **hyper)
    info = json.loads(re.sub(r"Float\d+\[", "[", (d / "run_info.json").read_text()))  # wrapper writes "Float32[...]"
    fit = pd.read_csv(fo)
    back = {f"S{i}": c for i, c in enumerate(expr.columns)}  # wrapper's id mapping
    fit["sample"] = fit.ID.map(back)
    assert fit["sample"].notna().all() and len(fit) == expr.shape[1]
    fit = fit.set_index("sample")
    log(f"{name} seed {seed}: done in {info['runtime_s'] / 60:.1f} min, {info['n_seed_genes_used']} seed genes, "
        f"{info['n_eigengenes']} eigengenes, aligned_by {info['aligned_by']}")
    return fit, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--cohorts", nargs="+", default=list(LOADERS))
    a = ap.parse_args()
    sfx = "_check" if a.check else ""
    OUT.mkdir(exist_ok=True)
    res, runs = [], []
    for name in a.cohorts:
        expr, meta = LOADERS[name]()
        log(f"{name}: {int(meta.is_tumour.sum())} luminal A tumours, {int((~meta.is_tumour).sum())} normals, {expr.shape[0]} genes")
        if a.check:  # stratified subsample, 40 per covariate level (tumour + normal in TCGA)
            keep = meta.groupby("cov", group_keys=False).apply(lambda g: g.sample(min(len(g), 40), random_state=0)).index
            expr, meta = expr[keep], meta.loc[keep]
        meta = meta.loc[expr.columns]
        out = meta[["is_tumour"]].copy().rename_axis("sample")
        for seed in SEEDS:
            fit, info = fit_one(name, expr, meta, seed, a.threads, a.check)
            mag = fit.Magnitude.astype(float)
            tum = meta.is_tumour.to_numpy()
            out[f"magnitude_s{seed}"] = mag.reindex(out.index)
            out[f"z_s{seed}"] = (mag - mag[tum].mean()) / mag[tum].std(ddof=1)
            out[f"phase_s{seed}"] = fit.Phase.astype(float).reindex(out.index)
            runs.append(dict(cohort=name, seed=seed, n_samples=expr.shape[1], n_seed_genes_used=info["n_seed_genes_used"],
                             n_eigengenes=info["n_eigengenes"], runtime_s=info["runtime_s"]))
        out["magnitude"] = out[[f"z_s{s}" for s in SEEDS]].mean(1)
        out = out.drop(columns=[f"z_s{s}" for s in SEEDS]).reset_index().assign(cohort=name)
        t = out[out.is_tumour]
        for i, j in [(0, 1), (0, 2), (1, 2)]:
            log(f"{name} seed concordance magnitude s{i}-s{j} (tumours): Spearman "
                f"{spearmanr(t[f'magnitude_s{i}'], t[f'magnitude_s{j}'])[0]:.3f}")
        res.append(out)
        cols = ["cohort", "sample", "is_tumour"] + [f"magnitude_s{s}" for s in SEEDS] + ["magnitude"] + [f"phase_s{s}" for s in SEEDS]
        pd.concat(res)[cols].to_csv(OUT / f"cyclops_magnitude{sfx}.csv", index=False)  # written after each cohort
        pd.DataFrame(runs).to_csv(OUT / f"cyclops_runs{sfx}.csv", index=False)
    log("done")


if __name__ == "__main__":
    main()
