"""A134 D134.3 (descriptive): per-cell in-plane dispersion D with the common GTEx gene phases vs patient-specific gene
phases fitted on one random half of the patient's cells and evaluated on the other half (cross-fit), malignant vs normal
epithelium, on the 64 patients with >= 100 cells in both groups.

Patient-specific plane: top two eigenvectors of the half's residual cross-gene covariance (diagonal removed), rescaled to the
mean squared loading of the waveform, entering `vonmises._split_projectors` in place of the waveform. D on the held-out half = `_plane` of its
residual covariance. Each D has a null from 10 permuted planes (common: permuted waveform; patient-specific: fitted
loadings with gene rows permuted, 5 per fit); D minus null is reported. 5 random half splits x both directions per group.
Cells: legacy dispersion tables, malignant = tumour epithelium with the stored copy-number call `malignant`; normal =
adjacent-normal epithelium. Interval: patient bootstrap (1,000) of cell-weighted means.

    python score_d134_3.py
Outputs: output/a134_d134_3_groups.csv (per patient and group), a134_d134_3.csv (summary).
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "common" / "datasets"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.cells import vonmises as VM  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import data  # noqa: E402
from cancer_common import studies as S  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
GENES = list(cs.gtex().genes)
N_SPLITS, N_REP, N_NULL_C, N_NULL_P = 30, 5, 10, 5
CW = cs.gtex().amp * np.exp(-1j * cs.gtex().psi)


def resid_cov(C, N):
    r, w, scale = VM._residuals(C, N, 1, None)
    return (w[:, None] * r).T @ r / len(r) * scale


def fitted_plane_c(M):
    A = (M + M.T) / 2
    np.fill_diagonal(A, 0.0)
    val, vec = np.linalg.eigh(A)
    o = np.argsort(-np.abs(val))[:2]
    X = vec[:, o]
    X = X * np.sqrt(np.mean(np.abs(CW) ** 2) / np.mean(np.sum(X ** 2, 1)))  # waveform scale: D comparable with the common-phase D
    return X[:, 0] - 1j * X[:, 1], X  # `_split_projectors` takes X = [c.real, -c.imag]


def proj(X, rng):
    """Split projectors of a fitted plane (balance threshold relaxed to 0.05 until 10 splits exist, else None; fitted planes can be elongated)."""
    from clocksync.cells.dispersion import _projector, splits
    for mc in (0.2, 0.1, 0.05):
        sp_ = splits(X, N_SPLITS, rng, min_cond=mc)
        if len(sp_) >= 10:
            return [(_projector(X[A], np.ones(len(A))), _projector(X[B], np.ones(len(B))), A, B) for A, B in sp_]
    return None  # plane too elongated for balanced well-conditioned gene halves: the fit is dropped


def main():
    limit_threads(2)
    warnings.filterwarnings("ignore")
    ok = pd.read_csv(OUT / "scout_d134_3.csv", dtype={"patient": str})
    ok = ok[ok.half_ok_scout]
    assert len(ok) == 64
    wf = cs.gtex()
    cw = wf.amp * np.exp(-1j * wf.psi)
    rng = np.random.default_rng(1343)
    LA = VM._split_projectors(cw, N_SPLITS, rng)
    LAn = [VM._split_projectors(p.amp * np.exp(-1j * p.psi), N_SPLITS, rng) for p in [wf.permuted(rng) for _ in range(N_NULL_C)]]
    rows = []
    for study, w in ok.groupby("study"):
        cells = pd.read_csv(Path(S.DATA.parent / "x") if False else Path("/home/maxine/Documents/paychere/scCircadianCancer/outputs/celldisp") / f"{study}_cells.csv.gz",
                            dtype=dict(resection=str, stratum=str, sample=str, barcode=str), keep_default_na=False)
        cells = cells[cells.resection.isin(set(w.patient))]
        calls = S.read_calls(study)
        epi = cells[cells.compartment == "tumour:epithelial"]
        cl = calls.reindex(pd.MultiIndex.from_arrays([epi["sample"], epi.barcode])).to_numpy()
        groups = {"malignant": epi[cl == "malignant"], "normal": cells[cells.compartment == "adjacent_normal:epithelial"]}
        for gname, g in groups.items():
            for pat, gg in g.groupby("resection"):
                if len(gg) < 100:
                    continue
                C, N = gg[GENES].to_numpy(float), gg.depth.to_numpy(float)
                acc = dict(common=[], common_null=[], pat=[], pat_null=[])
                fails = 0
                for _ in range(N_REP):
                    h = rng.permutation(len(gg)) % 2 == 0
                    for fit, ev in ((h, ~h), (~h, h)):
                        M1 = resid_cov(C[fit], N[fit])
                        M2 = resid_cov(C[ev], N[ev])
                        cp, X = fitted_plane_c(M1)
                        LP = proj(X, rng)
                        if LP is None:
                            fails += 1
                            continue
                        acc["common"].append(VM._plane(M2, LA)[0])
                        acc["common_null"].append(np.mean([VM._plane(M2, L)[0] for L in LAn]))
                        acc["pat"].append(VM._plane(M2, LP)[0])
                        nul = []
                        for _ in range(N_NULL_P):
                            LQ = proj(X[rng.permutation(len(X))], rng)
                            if LQ is not None:
                                nul.append(VM._plane(M2, LQ)[0])
                        acc["pat_null"].append(np.mean(nul) if nul else np.nan)
                rows.append(dict(study=study, patient=pat, group=gname, cells=len(gg), fits=len(acc["pat"]), fit_fail=fails,
                                 **{k: float(np.nanmean(v)) if len(v) else np.nan for k, v in acc.items()}))
        print(study, flush=True)
    d = pd.DataFrame(rows)
    d["D_common"] = d.common - d.common_null
    d["D_patient"] = d.pat - d.pat_null
    d.to_csv(OUT / "a134_d134_3_groups.csv", index=False)
    # summary: cell-weighted means over patients present in both groups, bootstrap over patients
    w = d.pivot_table(index=["study", "patient"], columns="group", values=["cells", "D_common", "D_patient", "common", "pat"]).dropna()
    n = len(w)

    def est(idx):
        s = w.iloc[idx]
        out = {}
        for g in ("malignant", "normal"):
            c = s[("cells", g)].to_numpy()
            for k in ("D_common", "D_patient", "common", "pat"):
                out[f"{k}_{g}"] = float(np.sum(c * s[(k, g)]) / c.sum())
        for k in ("D_common", "D_patient"):
            out[f"{k}_mal_minus_norm"] = out[f"{k}_malignant"] - out[f"{k}_normal"]
        for g in ("malignant", "normal"):
            out[f"crossfit_gain_{g}"] = out[f"D_patient_{g}"] - out[f"D_common_{g}"]
        out["crossfit_gain_mal_minus_norm"] = out["crossfit_gain_malignant"] - out["crossfit_gain_normal"]
        return out
    obs = est(np.arange(n))
    b = pd.DataFrame([est(rng.integers(0, n, n)) for _ in range(1000)])
    s = pd.DataFrame({k: dict(estimate=v, lo=b[k].quantile(.025), hi=b[k].quantile(.975)) for k, v in obs.items()}).T.assign(patients=n)
    s.to_csv(OUT / "a134_d134_3.csv")
    print(s.round(3).to_string())


if __name__ == "__main__":
    main()
