"""A214 scout, step 1: power of the pooled per-cell D contrast (malignant minus normal epithelium, same patient) and
identifiability of damping vs dephasing. No malignant-vs-normal D difference and no amplitude is computed.

Effect sizes (exact von Mises model, `clocksync.cells.vonmises.model_stats`, no data):
- rewiring: each gene's acrophase shifted by N(0, sigma^2) in the malignant cells (independent draws per tumour; 40 draws
  averaged), D relative to the intact clock, sigma 1-8 h and uniform phases;
- damping vs dephasing at a matched pseudobulk amplitude ratio a = rho_mal / rho_norm (rho = lam R): damping lam -> a lam
  (R kept), dephasing R -> a R (lam kept); D relative to the intact clock. Grid of the intact (lam, kappa).
D is the estimator of `compute_cancer_split` (cell-weighted group D minus its mean over 10 acrophase-permuted planes).

Measurement noise (real cells, blind to the contrast): per eligible group (patient x compartment, >= 100 cells), D on two
random halves of the cells; the half difference dh has mean 0 whatever the true D, and dh^2 / 4 estimates the sampling
variance of the full-group D (variance ~ 1 / cells). Only dh, never D, is written. MDE = 2.49 x SE of the pooled mean
difference, SE^2 = sum w^2 (s2_mal + s2_norm) / (sum w)^2, w = equal or min(cells); true between-patient heterogeneity is
not included, so the MDE is a lower bound.

Depth needed (`--depth-only` reruns this part from noise_groups.csv): median regression of log |dh| / 2 on log clock UMI
per cell and log cells; per-group SE at a given depth and 500 cells (median-based, x 1 and x 2 for heavier tails), and the
patients needed for MDE 1.0 (about the 2 h bound at D 3.9) with both compartments at that depth.

    python power_scout.py [--jobs 4] [--depth-only]
Outputs: ../output/power/{model_rewiring,model_damp_dephase,noise_groups,power,depth_scaling}.csv
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import optimize, special

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import inventory_disk as inv  # noqa: E402  (loaders, eligibility constants)
from clocksync import gtex  # noqa: E402
from clocksync.cells import vonmises as VM  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[0] / "output" / "power"
WF = gtex()
CW = WF.amp * np.exp(-1j * WF.psi)
GENES = list(WF.genes)
N_SPLITS, N_NULL, MIN_CELLS, MIN_UMI = 30, 10, 100, 2.0
MU = np.linspace(0, 2 * np.pi, 12, endpoint=False)
SEED = 214


def planes():
    """Split projectors of the waveform and of 10 permuted planes (one stream, as `compute_cancer_split`)."""
    rng = np.random.default_rng(0)
    LA = VM._split_projectors(CW, N_SPLITS, rng)
    return LA, [VM._split_projectors(p.amp * np.exp(-1j * p.psi), N_SPLITS, rng) for p in [WF.permuted(rng) for _ in range(N_NULL)]]


LA, LAP = planes()


def R_of(k):
    return float(special.ive(1, k) / special.ive(0, k))


def kappa_of(R):
    return float(optimize.brentq(lambda lk: R_of(np.exp(lk)) - R, np.log(1e-4), np.log(500)))


def d_model(lam, kappa, c):
    """Expected D minus its permuted-plane null under the model."""
    return VM.model_stats(lam, kappa, c, LA, MU)[0] - np.mean([VM.model_stats(lam, kappa, c, L, MU)[0] for L in LAP])


# ---- effect sizes -------------------------------------------------------------------------------------------------

def model_rewiring() -> pd.DataFrame:
    rows = []
    for R in (0.25, 0.45, 0.7):
        k = kappa_of(R)
        for lam in (1.5, 2.5):
            d0 = d_model(lam, k, CW)
            for sig in (1, 2, 3, 4, 6, 8, "uniform"):
                rng = np.random.default_rng(SEED)
                eps = [rng.uniform(0, 2 * np.pi, len(CW)) if sig == "uniform" else rng.normal(0, sig * 2 * np.pi / 24, len(CW)) for _ in range(40)]
                r = np.array([d_model(lam, k, CW * np.exp(-1j * e)) for e in eps]) / d0
                rows.append(dict(R=R, lam=lam, sigma_h=sig, D_ratio=r.mean(), D_ratio_sd_per_tumour=r.std(ddof=1)))
    return pd.DataFrame(rows)


def model_damp_dephase() -> pd.DataFrame:
    rows = []
    for R in (0.25, 0.45, 0.7):
        for lam in (1.5, 2.5):
            d0 = d_model(lam, kappa_of(R), CW)
            for a in (0.5, 0.7, 0.85):
                dd = d_model(a * lam, kappa_of(R), CW) / d0
                dp = d_model(lam, kappa_of(a * R), CW) / d0
                rows.append(dict(R=R, lam=lam, a=a, D_ratio_damping=dd, D_ratio_dephasing=dp, gap=dp - dd))
    return pd.DataFrame(rows)


# ---- measurement noise --------------------------------------------------------------------------------------------

def d_group(C, N, rng):
    return VM.group_stats(C, N, LA, rng=rng)["D"] - np.mean([VM.group_stats(C, N, L, rng=rng)["D"] for L in LAP])


def study_noise(study: str) -> list[dict]:
    limit_threads(1)
    warnings.filterwarnings("ignore")
    P, grp, _, _, cells = inv.cohort(study)
    if "barcode" not in cells or not len(P):
        return []
    rng = np.random.default_rng([SEED, sum(map(ord, study))])
    out = []
    for (pt, comp), g in grp.items():
        if comp not in (inv.MAL, inv.ADJ, inv.NMAL) or len(g) < MIN_CELLS:
            continue
        C, N = g[GENES].to_numpy(float), g.depth.to_numpy(float)
        h = rng.permutation(len(N)) % 2 == 0
        dh = d_group(C[h], N[h], rng) - d_group(C[~h], N[~h], rng)
        u = C.sum(1)
        out.append(dict(study=study, patient=pt, compartment=comp, n_cells=len(N), clock_umi_median=float(np.median(u)),
                        clock_umi_mean=float(u.mean()), dh=dh))
    return out


def power(noise: pd.DataFrame) -> pd.DataFrame:
    noise = noise.assign(s2=noise.dh ** 2 / 4)
    rows = []
    for design, norms in (("a_adjacent", [inv.ADJ]), ("b_in_tumour", [inv.NMAL]), ("either_a_first", [inv.ADJ, inv.NMAL])):
        for umi_rule in (False, True):
            ok = noise if not umi_rule else noise[noise.clock_umi_median >= MIN_UMI]
            g = ok.set_index(["study", "patient", "compartment"])
            pairs = []
            for (st, pt), _ in ok.groupby(["study", "patient"]):
                if (st, pt, inv.MAL) not in g.index:
                    continue
                for nc in norms:
                    if (st, pt, nc) in g.index:
                        m, n = g.loc[(st, pt, inv.MAL)], g.loc[(st, pt, nc)]
                        pairs.append(dict(s2=m.s2 + n.s2, w=min(m.n_cells, n.n_cells)))
                        break
            p = pd.DataFrame(pairs)
            if not len(p):
                rows.append(dict(design=design, umi_rule=umi_rule, patients=0))
                continue
            for wt in ("equal", "min_cells"):
                w = np.ones(len(p)) if wt == "equal" else p.w.to_numpy(float)
                se = float(np.sqrt(np.sum(w ** 2 * p.s2) / w.sum() ** 2))
                rows.append(dict(design=design, umi_rule=umi_rule, weights=wt, patients=len(p), studies=ok.study.nunique(),
                                 se=se, mde=2.49 * se, per_patient_sd=float(np.sqrt(p.s2.mean()))))
    return pd.DataFrame(rows)


def depth_scaling(noise: pd.DataFrame) -> pd.DataFrame:
    import statsmodels.api as sm
    n = noise[noise.dh != 0]
    X = sm.add_constant(np.column_stack([np.log(n.clock_umi_mean), np.log(n.n_cells)]))
    fit = sm.QuantReg(np.log(np.abs(n.dh) / 2), X).fit(q=0.5)
    b, se_b = np.asarray(fit.params), np.asarray(fit.bse)
    rows = []
    for umi in (1, 2, 4, 8, 16, 30):
        se = float(np.exp(b[0] + b[1] * np.log(umi) + b[2] * np.log(500)) / 0.6745)
        rows.append(dict(clock_umi_per_cell=umi, cells=500, group_se=se, slope_umi=b[1], slope_umi_se=se_b[1], slope_cells=b[2],
                         slope_cells_se=se_b[2], n_for_mde1=int((2.49 * np.sqrt(2) * se) ** 2),
                         n_for_mde1_tails_x2=int((2.49 * np.sqrt(2) * 2 * se) ** 2)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--depth-only", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    if args.depth_only:
        d = depth_scaling(pd.read_csv(OUT / "noise_groups.csv"))
        d.to_csv(OUT / "depth_scaling.csv", index=False)
        print(d.round(3).to_string(index=False))
        return
    model_rewiring().to_csv(OUT / "model_rewiring.csv", index=False)
    model_damp_dephase().to_csv(OUT / "model_damp_dephase.csv", index=False)
    print("model done", flush=True)
    studies = pd.read_csv(HERE.parents[0] / "inventory" / "inventory_studies.csv")
    studies = studies[(studies.source == "per-cell table") & (studies.n_eligible_a.fillna(0) + studies.n_eligible_b.fillna(0) > 0)].study
    res = Parallel(n_jobs=args.jobs)(delayed(study_noise)(s) for s in studies)
    noise = pd.DataFrame([r for rr in res for r in rr])
    noise.to_csv(OUT / "noise_groups.csv", index=False)
    p = power(noise)
    p.to_csv(OUT / "power.csv", index=False)
    depth_scaling(noise).to_csv(OUT / "depth_scaling.csv", index=False)
    pd.set_option("display.width", 200)
    print(p.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
