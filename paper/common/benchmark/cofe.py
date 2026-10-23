"""COFE (sparse cyclic PCA with cross-validated sparsity; Ananthasubramaniam et al. 2025) on a samples x genes table.

Run as in the COFE paper's TCGA analysis: standardised features, 3-repeated 5-fold cross-validation with 5 restarts
per fit. The sparsity grid spans 1 (one feature) to sqrt(p) (no sparsity). COFE's own validity signal is a clear
minimum of the cross-validation error; `interior_min` records whether the minimum is inside the grid.

The fits are COFE's (`sparse_cyclic_pca_masked`, `sparse_cyclic_pca`) and the folds, seeds, restart selection and
error are those of `COFE.cross_validate`, but its driver is replaced, for two reasons. (1) It keeps a full AnnData per
fit (n_s x K x repeats x restarts = 1,500 fits; > 100 GB for BRCA). (2) Every fit writes its scores into the AnnData
it is given, and all fits receive the same object: run serially every score is the last fit's, and in parallel the
fits joblib batches into one worker share one object (228 of 300 cross-validation scores wrong in a test with 4
workers). Here each fit gets its own AnnData and returns only what is needed.
"""

from __future__ import annotations

import warnings

import anndata as ad
import COFE
import numpy as np
import pandas as pd
from COFE.analyse import _shuffled_checkerboard
from COFE.scpca import sparse_cyclic_pca, sparse_cyclic_pca_masked
from joblib import Parallel, delayed

MAX_THREADS = 16

FIT = dict(tol=1e-3, tol_z=1e-6, max_iter=400, scale_by_features=False)  # COFE.cross_validate defaults


def _cv_fit(X, s, cv_ind, seed):
    """One masked fit: (training rss, error on the masked entries)."""
    mask = np.zeros(X.size, dtype=bool)
    mask[cv_ind] = True
    a = ad.AnnData(np.asarray(X))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc = sparse_cyclic_pca_masked(a, mask.reshape(X.shape), s=s, seed=seed, **FIT).uns["scpca"]
    return sc["rss"], sc["cv_err"]


def _full_fit(X, s, seed):
    """One unmasked fit (a restart of COFE's final fit): rss, cyclic components, loadings, scores."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = sparse_cyclic_pca(ad.AnnData(np.asarray(X)), s=s, seed=seed, tol=FIT["tol"], max_iter=FIT["max_iter"],
                              scale_by_features=FIT["scale_by_features"])
    return a.uns["scpca"]["rss"], a.obsm["X_scpca"], a.varm["PCs"], a.uns["scpca"]


def run_cofe(logexpr: pd.DataFrame, s_choices=None, n_s: int = 20, seed: int = 0, ncores: int = MAX_THREADS,
             K: int = 5, repeats: int = 3, restarts: int = 5) -> dict:
    """`logexpr`: samples x genes (log scale). Returns phase (h, arbitrary origin and direction), the
    cross-validation curve, the chosen sparsity, and the loadings of the two cyclic components."""
    a = ad.AnnData(logexpr.to_numpy(float), obs=pd.DataFrame(index=logexpr.index.astype(str)),
                   var=pd.DataFrame(index=logexpr.columns.astype(str)))
    a.obs["train"] = True
    a.uns["period"] = 24.0
    a = COFE.preprocess_data(a, scale=True)
    if s_choices is None:
        s_choices = np.linspace(1.0, np.sqrt(a.n_vars), n_s)
    s_choices = np.asarray(s_choices, float)

    rng = np.random.default_rng(seed=seed)
    seeder = lambda: rng.integers(0, 2**32 - 1)
    folds = _shuffled_checkerboard(a.shape, K, repeats, seed=seeder())
    jobs = [(s, f, seeder()) for s in s_choices for f in range(len(folds)) for _ in range(restarts)]
    out = Parallel(n_jobs=ncores, backend="loky", backend_kwargs=dict(inner_max_num_threads=1))(
        delayed(_cv_fit)(a.X, s, folds[f], sd) for s, f, sd in jobs)
    rss, cv_err = np.array(out, float).reshape(len(s_choices), repeats, K, restarts, 2).transpose(4, 0, 1, 2, 3)
    # per fold the restart with the lowest training rss; error per masked entry, per repeat (as COFE)
    err = np.take_along_axis(cv_err, rss.argmin(-1)[..., None], -1)[..., 0]
    size = np.array([len(f) for f in folds]).reshape(repeats, K)
    ok = ~np.isnan(err)
    per_repeat = np.where(ok, err, 0).sum(-1) / (ok * size).sum(-1)
    cv = pd.DataFrame({"cv_err": per_repeat.mean(1), "cv_sd": per_repeat.std(1), "s": s_choices})
    best = int(np.argmin(cv.cv_err))

    final = np.random.default_rng(seeder())  # COFE's _multi_start seeding
    fits = Parallel(n_jobs=min(ncores, restarts))(
        delayed(_full_fit)(a.X, s_choices[best], final.integers(0, 2**32 - 1)) for _ in range(restarts))
    _, a.obsm["X_scpca"], a.varm["PCs"], sc = min(fits, key=lambda f: f[0])
    a.uns["scpca"] = {**sc, "best_s": s_choices[best], "s_choices": s_choices,
                      "runs": cv[["cv_err", "cv_sd"]].to_numpy().tolist()}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = COFE.estimate_phase(a, train=True)
    return {
        "phase": a.obs["phase"].reindex(logexpr.index.astype(str)).set_axis(logexpr.index).rename("phase"),
        "cv": cv,
        "best_s": float(s_choices[best]),
        "interior_min": 0 < best < len(cv) - 1,
        "loadings": pd.DataFrame(a.varm["PCs"], index=a.var_names, columns=["v1", "v2"]),
    }
