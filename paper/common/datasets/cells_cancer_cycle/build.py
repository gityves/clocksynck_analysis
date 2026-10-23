"""CoPhaser cell-cycle phase of every cell of the single-cell cancer cohorts (input of the cycling checks C2, C2b, C3, C4 of the
2026-10-04 amendments of docs/cancer/ANALYSIS_PLAN.md).

Source:   the cohort's tumour and adjacent-normal samples (`cancer_common.studies`: same samples and majority chemistry as
          `cancer_sc`), at most min(3,000, 150,000 / samples) cells per sample (seeded subsample); one CoPhaser model per cohort
          (cell-cycle mode, packaged human cell-cycle genes, 2,000 HVGs as context, RPE1 Fourier warm start, trainer from
          `auto_hyperparameters(cycle="cell_cycle")`, all compartments together, no batch term; clock and clock-like genes removed
          from both gene sets; trained on <= 80,000 sample-stratified cells, every loaded cell then inferred). CoPhaser is used as
          installed, not modified.
Table (`$SCCC_DATA/cells/cancer_cycle/<study>.csv.gz`, one row per cell, index sample_barcode):
          sample, patient, kind (tumour / adjacent_normal), compartment (marker class), total_umi, `stratum` (chemistry),
          theta, cycling_logit (CoPhaser), S, G2M (mean log CP10k of the Tirosh lists), `clock:<gene>` and `pseudo0:<gene>`,
          `pseudo1:<gene>` (counts of the 14 clock genes and of two expression-matched pseudo-clock sets),
          and the annotation: `histone_frac` (replication-dependent histone fraction of the cell's UMI), `state_t1` (T1 score
          state of the epithelial cells: G2M_hi / S_hi / noncyc / blank; control-adjusted scores, thresholds over the epithelial
          cells of the table), `u` (phase oriented S -> G2M, 0 at the S-score peak) and `g2m_peak_u` (study constant).
Replaces: legacy/scripts/s3_cellcycle_cophaser.py (`cells`, `fit`, `orient`), s5_cc_diagnostics.py (`annotate`), s4_g2m_checks.py
          (`chemistry`: the stratum is stored with the cell).
Check:    `python build.py check STUDY`: refits one cohort and compares with the legacy table (needs the GPU; CoPhaser training is
          not bit-reproducible: agreement of theta (circular correlation), of the cycling call and of the states is reported).
          `from_legacy(study)` annotates the legacy fit instead of refitting (the annotation is deterministic).

    python build.py fit STUDY        # CoPhaser fit and annotation -> $SCCC_DATA/cells/cancer_cycle/STUDY.csv.gz
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(Path(__file__).resolve().parents[3])]
import clocksync as cs  # noqa: E402
from clocksync.genes import CLOCK_GENES, canonical, gene_columns  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common.paths import DATA  # noqa: E402

from cancer_common import studies as S  # noqa: E402

import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("cells_cancer_build", HERE.parents[1] / "cells_cancer" / "build.py")
CC = importlib.util.module_from_spec(_s)
_s.loader.exec_module(CC)

OUT = DATA / "cells" / "cancer_cycle"
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # paper/
from common.paths import ARCHIVE  # noqa: E402
LEGACY = ARCHIVE / "outputs" / "states" / "cc"
GENES = list(CLOCK_GENES)
MAX_TRAIN, MAX_PER_SAMPLE, MAX_CELLS = 80_000, 3000, 150_000
HISTONE = re.compile(r"^(H1-[0-6]$|H2AC\d|H2BC\d|H3C\d|H4C\d|HIST[123]H)")


def cells(study: str) -> ad.AnnData:
    us = [u for u in S.majority_chemistry(S.units(study)[0]) if u[2] in ("tumour", "adjacent_normal")]
    cap = min(MAX_PER_SAMPLE, MAX_CELLS // max(len(us), 1))
    parts = []
    for smp, pat, kind, ld, chem, _ in us:
        X, _, umi, comp, bc, sym = ld()
        if len(umi) > cap:
            k = np.sort(np.random.default_rng(0).choice(len(umi), cap, replace=False))
            X, umi, comp, bc = sp.csr_matrix(X)[k], umi[k], comp[k], np.asarray(bc)[k]
        obs = pd.DataFrame({"sample": smp, "patient": pat, "kind": kind, "compartment": comp, "total_umi": umi,
                            "stratum": chem}, index=[f"{smp}_{b}" for b in bc])
        a = ad.AnnData(sp.csr_matrix(X), obs=obs, var=pd.DataFrame(index=pd.Index([str(s).upper() for s in sym])))
        a.var_names_make_unique()
        parts.append(a)
    a = ad.concat(parts, join="inner")
    a.obs = a.obs.astype({c: str for c in ("sample", "patient", "kind", "compartment", "stratum")})
    a.layers["total"] = a.X
    return a


def mean_log(X, umi, cols):
    M = X[:, cols].multiply(1e4 / np.maximum(umi, 1)[:, None]).tocsr()
    M.data = np.log1p(M.data)
    return np.asarray(M.mean(1)).ravel()


def orient(o: pd.DataFrame):
    """u in [0, 2 pi): 0 at the S-score peak, increasing towards the G2M-score peak (and the G2M peak on u)."""
    cyc = o.cycling_logit > 0
    z = np.exp(1j * o.theta[cyc])
    ps = np.angle(np.sum(z * (o.S[cyc] - o.S[cyc].mean())))
    pg = np.angle(np.sum(z * (o.G2M[cyc] - o.G2M[cyc].mean())))
    sign = 1.0 if np.mod(pg - ps, 2 * np.pi) < np.pi else -1.0
    return np.mod(sign * (o.theta.to_numpy() - ps), 2 * np.pi), float(np.mod(sign * (pg - ps), 2 * np.pi))


def score_states(cells: pd.DataFrame) -> np.ndarray:
    """T1 states: G2M_hi (top quartile of G2M), S_hi (top quartile of S, not G2M_hi), noncyc (bottom half of both)."""
    g2m = cells.G2M >= cells.G2M.quantile(0.75)
    s = (cells.S >= cells.S.quantile(0.75)) & ~g2m
    non = (cells.G2M <= cells.G2M.quantile(0.5)) & (cells.S <= cells.S.quantile(0.5))
    return np.where(g2m, "G2M_hi", np.where(s, "S_hi", np.where(non, "noncyc", "")))


def annotate(study: str, o: pd.DataFrame, a: ad.AnnData | None = None) -> pd.DataFrame:
    """Histone fraction, T1 score state and oriented phase of the cells of a CoPhaser table (same cells as the fit)."""
    a = cells(study) if a is None else a
    a = a[a.obs_names.isin(o.index)]
    o = o.loc[a.obs_names].copy()
    if "stratum" not in o:
        o["stratum"] = a.obs.stratum.to_numpy()
    X = sp.csr_matrix(a.X)
    umi = o.total_umi.to_numpy(float)
    hist = np.flatnonzero([bool(HISTONE.match(g)) for g in a.var_names])
    o["histone_frac"] = np.asarray(X[:, hist].sum(1)).ravel() / umi
    epi = (o.compartment == "epithelial").to_numpy()
    M = X[epi].multiply(1e4 / np.maximum(umi[epi], 1)[:, None]).tocsr()
    M.data = np.log1p(M.data)
    sets = CC.control_sets(a.var_names, np.asarray(M.mean(0)).ravel(), np.random.default_rng(0),
                           {"S": CC.S_GENES, "G2M": CC.G2M_GENES})
    sc = CC.scores(X[epi], umi[epi], sets)
    t1 = pd.Series("", index=o.index)
    t1[epi] = score_states(sc)
    o["state_t1"] = t1
    o["u"], g2m_u = orient(o)
    o["g2m_peak_u"] = g2m_u
    return o


def fit(study: str, seed: int = 0) -> pd.DataFrame:
    import torch
    from importlib import resources

    from cophaser import CoPhaser, DecoderPrior, Loss, Trainer, auto_hyperparameters, gene_sets, utils
    _tensor = torch.tensor  # compat shim (not a CoPhaser change): torch 2.9 / pandas 3 cannot build a tensor from a Series

    def tensor(x, *args, **kw):
        return _tensor(x.to_numpy() if isinstance(x, pd.Series) else x, *args, **kw)
    torch.tensor = tensor
    torch.manual_seed(seed)
    np.random.seed(seed)
    a = cells(study)
    X = sp.csr_matrix(a.X)
    umi = a.obs.total_umi.to_numpy(float)
    names = pd.Index(a.var_names)
    lm = pd.Series(np.log1p(np.asarray(X[:: max(1, a.n_obs // 20000)].mean(0)).ravel()), index=names)
    wf = cs.gtex()
    sets = {"clock": GENES}
    for j, w in enumerate(pseudo_waveforms(lm, 2, np.random.default_rng(1000), wf)):
        sets[f"pseudo{j}"] = list(w.genes)
    keep_cc = [g for g in gene_sets.HUMAN_CELL_CYCLE_GENE_SET if g in names and canonical(g) not in CC.CLOCK_LIKE]
    hvg = utils.get_variable_genes(ad.AnnData(X=X, var=pd.DataFrame(index=names)), n_variable_genes=2000)
    hvg = [g for g in hvg if canonical(g) not in CC.CLOCK_LIKE]
    rng = np.random.default_rng(seed)
    per = MAX_TRAIN // a.obs["sample"].nunique()
    tr = np.sort(np.concatenate([rng.choice(ix, min(per, len(ix)), replace=False)
                                 for ix in a.obs.groupby("sample", observed=True).indices.values()]))
    genes = sorted(set(keep_cc) | set(hvg))
    b = a[:, genes].copy()
    model = CoPhaser(keep_cc, context_genes=hvg)
    model.load_anndata(b[tr].copy(), layer_to_use="total", library_size_field="total_umi")
    DecoderPrior.load_fourier_coefficients_prior(resources.files("cophaser") / "resources" / "fourier_coefficients_RPE.csv", model)
    hp = auto_hyperparameters(model, cycle="cell_cycle", verbose=False)
    t0 = time.time()
    Trainer(model, Loss.compute_loss, **hp["trainer"]).train_model(device="cuda", silent=True)
    model.load_anndata(b, layer_to_use="total", library_size_field="total_umi", keep_mean_genes=True)
    model.to("cpu")
    model.eval()
    th, lg = [], []
    with torch.no_grad():
        for i in range(0, model.variable_genes.shape[0], 20000):
            sl = slice(i, i + 20000)
            _, s = model.forward(model.variable_genes[sl], model.rhythmic_genes[sl], model.library_size[sl], use_max_posterior=True)
            mt = s["mu_theta"].numpy()
            th.append(np.arctan2(mt[:, 1], mt[:, 0]))
            lg.append(s["cycling_logits"].numpy().ravel())
    o = a.obs.copy()
    o["theta"], o["cycling_logit"] = np.concatenate(th), np.concatenate(lg)
    o["S"] = mean_log(X, umi, [names.get_loc(g) for g in CC.S_GENES if g in names])
    o["G2M"] = mean_log(X, umi, [names.get_loc(g) for g in CC.G2M_GENES if g in names])
    for k, gl in sets.items():
        cols = gene_columns(names, gl)
        C = X[:, cols[cols >= 0]].toarray()
        for g, c in zip([g for g, c in zip(GENES, cols) if c >= 0], C.T):
            o[f"{k}:{g}"] = c
    o.attrs["minutes"] = (time.time() - t0) / 60
    return annotate(study, o, a)


def build(study: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    o = fit(study)
    o.to_csv(OUT / f"{study}.csv.gz")
    json.dump(dict(study=study, n_cells=len(o), minutes=o.attrs.get("minutes")), open(OUT / f"{study}.json", "w"))
    return OUT / f"{study}.csv.gz"


def from_legacy(study: str) -> pd.DataFrame:
    """The annotated table of the legacy CoPhaser fit (`outputs/states/cc/<study>.csv.gz`)."""
    o = pd.read_csv(LEGACY / f"{study}.csv.gz", index_col=0, dtype=dict(sample=str, patient=str))
    return annotate(study, o)


def load(study: str) -> pd.DataFrame:
    """Stored table, else the annotated legacy fit, else a new fit."""
    f = OUT / f"{study}.csv.gz"
    if f.exists():
        return pd.read_csv(f, index_col=0, dtype=dict(sample=str, patient=str, stratum=str), keep_default_na=False)
    if (LEGACY / f"{study}.csv.gz").exists():
        return from_legacy(study)
    build(study)
    return load(study)


def check(study: str) -> None:
    """Refit one cohort and compare with the legacy fit: the angle theta is only defined up to rotation and reflection, so the
    comparison is on the oriented phase u and on the calls."""
    new = fit(study)
    old = annotate(study, pd.read_csv(LEGACY / f"{study}.csv.gz", index_col=0, dtype=dict(sample=str, patient=str)))
    i = new.index.intersection(old.index)
    d = new.loc[i, "u"].to_numpy() - old.loc[i, "u"].to_numpy()
    r = float(np.abs(np.mean(np.exp(1j * d))))
    cyc = new.loc[i, "cycling_logit"].to_numpy() > 0, old.loc[i, "cycling_logit"].to_numpy() > 0
    print(f"{study}: {len(new)} cells (legacy {len(old)}, common {len(i)}); oriented phase u: |mean e^i(u - u_legacy)| {r:.3f}; "
          f"cycling call agreement {np.mean(cyc[0] == cyc[1]):.3f} (cycling fraction {cyc[0].mean():.3f} vs {cyc[1].mean():.3f}); "
          f"S, G2M equal: {bool(np.allclose(new.loc[i, 'S'], old.loc[i, 'S']))}; T1 states equal: "
          f"{bool((new.loc[i, 'state_t1'] == old.loc[i, 'state_t1']).all())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=["fit", "check"])
    ap.add_argument("study")
    a = ap.parse_args()
    print(build(a.study)) if a.cmd == "fit" else check(a.study)
