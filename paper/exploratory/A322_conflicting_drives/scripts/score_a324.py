"""A324 (registered): replication of A322 in Xin 2021 (DRF vs NRF; liver, heart, kidney, VAT).

A322's estimation and chord fit (functions of `score_a322.py`), seed 324. H324.1: chord rejected (parametric bootstrap
p < 0.05). H324.2: every tissue with Mahalanobis distance > 3.84 lies above the chord (|r| > |model|).

    python score_a324.py
Outputs (`../output/xin/`): ratios.csv, fit.csv, tests.csv.
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2  # noqa: E402

spec = importlib.util.spec_from_file_location("a324_a322", HERE.parent / "score_a322.py")
B = importlib.util.module_from_spec(spec)
sys.modules["a324_a322"] = B
spec.loader.exec_module(B)

OUT = HERE.parents[1] / "output" / "xin"
TISSUES = ("liver", "heart", "kidney", "VAT")
W = B.W


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(324)
    rows, r, R = [], [], []
    for t in TISSUES:
        s, genes = B.A297.liver_states("xin2021", t, "arm")
        a0, a1 = (s[s.arm == a] for a in ("NRF", "DRF"))
        A0, A1 = B.A297.amplitude(a0), B.A297.amplitude(a1)
        B0, B1 = B.A297.boot_amp(a0, rng), B.A297.boot_amp(a1, rng)
        rt, rb = A1 / A0, B1 / B0
        r.append(rt), R.append(rb)
        rows.append(dict(tissue=t, genes=len(genes), n_NRF=len(a0), n_DRF=len(a1), amp_ratio=abs(rt),
                         amp_lo90=float(np.quantile(abs(rb), 0.05)), amp_hi90=float(np.quantile(abs(rb), 0.95)),
                         shift_h=float(np.mod(np.angle(rt) / W, 24)), re=rt.real, im=rt.imag))
    r, R = np.array(r), np.array(R)
    cov = np.array([np.cov(B.vec(R[k]).T) for k in range(len(TISSUES))])
    P = np.linalg.inv(cov)
    Q, D, b = B.fit(r, P)
    model = 1 + b * (np.exp(1j * D * W) - 1)
    maha = [float(B.vec(r[k] - model[k]) @ P[k] @ B.vec(r[k] - model[k])) for k in range(len(TISSUES))]
    ratios = pd.DataFrame(rows).assign(b=b, model_amp=np.abs(model), model_shift_h=np.mod(np.angle(model) / W, 24),
                                       mahalanobis=maha, above=np.abs(r) > np.abs(model))
    ratios.to_csv(OUT / "ratios.csv", index=False)
    Lc = np.linalg.cholesky(cov)
    qs = []
    for _ in range(B.N_PAR):
        e = np.einsum("tij,tj->ti", Lc, rng.standard_normal((len(TISSUES), 2)))
        qs.append(B.fit(model + e[:, 0] + 1j * e[:, 1], P)[0])
    qs = np.array(qs)
    p_par = float((1 + np.sum(qs >= Q)) / (1 + B.N_PAR))
    far = ratios[ratios.mahalanobis > 3.84]
    h2 = bool(len(far) and far.above.all())
    fitrow = dict(Q=Q, delta_h=D, **{f"b_{t}": float(b[k]) for k, t in enumerate(TISSUES)})
    pd.DataFrame([fitrow]).to_csv(OUT / "fit.csv", index=False)
    tests = pd.DataFrame([dict(test="H324.1", Q=Q, df=2 * len(TISSUES) - len(TISSUES) - 1, p_chi2=float(chi2.sf(Q, len(TISSUES) - 1)),
                               p_parametric=p_par, rejected=p_par < 0.05),
                          dict(test="H324.2", tissues_far=";".join(far.tissue), all_above=h2,
                               passes=bool(p_par < 0.05 and h2))])
    tests.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 220)
    print(ratios.round(3).to_string(index=False))
    print(pd.Series(fitrow).round(3).to_string())
    print(tests.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
