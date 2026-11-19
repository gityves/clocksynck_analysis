"""A322 (registered): chord model of conflicting drives on Manella 2021 (liver, WAT, lung; AL vs daytime feeding).

Per tissue, A297's method: one projector over both arms, per-arm complex amplitude (OLS on ZT harmonics), r_t =
A_DF / A_AL, library bootstrap within ZT (2000) for its 2 x 2 covariance. Model: r_t = 1 + b_t (exp(i Delta) - 1),
Delta shared, b_t in [0, 1]; fit by minimising the summed Mahalanobis distance Q (grid on Delta, closed-form b_t).
H322.1: Q against chi-square (2 df) and a parametric bootstrap of Q (2000).

    python score_a322.py
Outputs (`../output/`): ratios.csv, fit.csv, tests.csv.
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2  # noqa: E402

spec = importlib.util.spec_from_file_location("a322_a297", PAPER / "exploratory/A297_feeding_model/scripts/score_a297.py")
A297 = importlib.util.module_from_spec(spec)
sys.modules["a322_a297"] = A297
spec.loader.exec_module(A297)

OUT = HERE.parents[1] / "output"
W = 2 * np.pi / 24
TISSUES = ("Liver", "WAT", "Lung")
GRID = np.arange(0, 24, 0.01)
N_PAR = 2000


def vec(z):
    return np.stack([np.real(z), np.imag(z)], -1)


def fit(r, P):
    """r: complex (T,), P: precision matrices (T, 2, 2). Returns Q, Delta_h, b (T,)."""
    best = (np.inf, None, None)
    y = vec(r - 1)
    for d in GRID:
        u = vec(np.exp(1j * d * W) - 1)
        num = np.einsum("i,tij,tj->t", u, P, y)
        den = np.einsum("i,tij,j->t", u, P, u)
        b = np.clip(num / den, 0, 1)
        e = y - b[:, None] * u
        q = float(np.einsum("ti,tij,tj->", e, P, e))
        if q < best[0]:
            best = (q, d, b)
    return best


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(322)
    rows, r, R = [], [], []
    for t in TISSUES:
        s, genes = A297.liver_states("manella2021", t, "arm")
        a0, a1 = (s[s.arm == a] for a in ("AlbCre_AL", "AlbCre_DF"))
        A0, A1 = A297.amplitude(a0), A297.amplitude(a1)
        B0, B1 = A297.boot_amp(a0, rng), A297.boot_amp(a1, rng)
        rt, rb = A1 / A0, B1 / B0
        r.append(rt), R.append(rb)
        rows.append(dict(tissue=t, genes=len(genes), n_AL=len(a0), n_DF=len(a1), amp_ratio=abs(rt),
                         amp_lo90=float(np.quantile(abs(rb), 0.05)), amp_hi90=float(np.quantile(abs(rb), 0.95)),
                         shift_h=float(np.mod(np.angle(rt) / W, 24)), re=rt.real, im=rt.imag))
    r, R = np.array(r), np.array(R)
    cov = np.array([np.cov(vec(R[k]).T) for k in range(len(TISSUES))])
    P = np.linalg.inv(cov)
    Q, D, b = fit(r, P)
    model = 1 + b * (np.exp(1j * D * W) - 1)
    ratios = pd.DataFrame(rows).assign(b=b, model_amp=np.abs(model), model_shift_h=np.mod(np.angle(model) / W, 24),
                                       mahalanobis=[float(vec(r[k] - model[k]) @ P[k] @ vec(r[k] - model[k])) for k in range(3)])
    ratios.to_csv(OUT / "ratios.csv", index=False)
    # parametric bootstrap of Q at the fitted parameters
    Lc = np.linalg.cholesky(cov)
    qs = []
    for _ in range(N_PAR):
        e = np.einsum("tij,tj->ti", Lc, rng.standard_normal((3, 2)))
        rs = model + e[:, 0] + 1j * e[:, 1]
        qs.append(fit(rs, P)[0])
    qs = np.array(qs)
    # library-bootstrap intervals of Delta and b (refit on each replicate's ratios, 200 replicates)
    reps = [fit(R[:, j], P) for j in range(200)]
    dd = np.array([x[1] for x in reps])
    dc = np.angle(np.exp(1j * (dd - D) * W)) / W + D
    bb = np.array([x[2] for x in reps])
    fitrow = dict(Q=Q, delta_h=D, delta_lo90=float(np.quantile(dc, 0.05)), delta_hi90=float(np.quantile(dc, 0.95)),
                  **{f"b_{t}": float(b[k]) for k, t in enumerate(TISSUES)},
                  **{f"b_{t}_lo90": float(np.quantile(bb[:, k], 0.05)) for k, t in enumerate(TISSUES)},
                  **{f"b_{t}_hi90": float(np.quantile(bb[:, k], 0.95)) for k, t in enumerate(TISSUES)})
    pd.DataFrame([fitrow]).to_csv(OUT / "fit.csv", index=False)
    p_chi, p_par = float(chi2.sf(Q, 2)), float((1 + np.sum(qs >= Q)) / (1 + N_PAR))
    tests = pd.DataFrame([dict(test="H322.1", Q=Q, p_chi2_2df=p_chi, p_parametric=p_par,
                               retained=bool(p_chi >= 0.05 and p_par >= 0.05))])
    tests.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 220)
    print(ratios.round(3).to_string(index=False))
    print(pd.Series(fitrow).round(3).to_string())
    print(tests.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
