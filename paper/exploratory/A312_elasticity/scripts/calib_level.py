"""Post hoc calibration of the A312 elasticity level (2026-11-24, after the user asked why nearly every tissue reads
above 1): what does `clocksync.elasticity` return when every tissue scales exactly in proportion with the person's
volume (true elasticity 1)?

Simulated cohort at the GTEx design: the 27 peripheral tissues of the A294b table with their shared variance S (gain
|z| = sqrt S, host of unit mean square) and departure D (psi = D / 2 per axis); persons with a volume A_p = exp(N(0,
0.352)) (the GTEx person spread of A331, rescaled to unit mean square) and a clock time uniform over the day; each
tissue present with probability 0.55; counts Poisson around the GTEx waveform with per-gene noise. Elasticity of six
target tissues spanning S (10 splits each), for three draws. No GTEx data is used beyond the per-tissue S and D.

    python calib_level.py
Output: ../output/calib_level.csv.
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import elasticity as E  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
SIGMA_A, N, PRESENCE, DEPTH, TAU = 0.352, 250, 0.55, 2e7, 0.1
TARGETS = ["Testis", "Skin - Sun Exposed (Lower leg)", "Lung", "Liver", "Artery - Tibial", "Thyroid"]


def cohort(t, wf, rng):
    a = np.exp(rng.normal(0, SIGMA_A, N))
    a /= np.sqrt(np.mean(a ** 2))
    theta = rng.uniform(0, 2 * np.pi, N)
    base = np.full(len(wf.genes), -10.0) if wf.mesor is None else np.asarray(wf.mesor, float)
    rows = []
    for k, r in t.iterrows():
        z = np.sqrt(r.S) * np.exp(1j * rng.uniform(0, 2 * np.pi))
        off = rng.normal(0, 0.5, len(wf.genes))
        for p in range(N):
            if rng.random() > PRESENCE:
                continue
            c = a[p] * z * np.exp(1j * theta[p])
            beta = np.array([c.real, c.imag]) + rng.normal(0, np.sqrt(r.own / 2), 2)
            lr = base + off + wf.X @ beta + rng.normal(0, TAU, len(wf.genes))
            rows.append([f"p{p:03d}", k, DEPTH, *rng.poisson(DEPTH * np.exp(lr))])
    return pd.DataFrame(rows, columns=["donor", "compartment", "depth", *wf.genes])


def main():
    warnings.filterwarnings("ignore")
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(HERE.parents[2] / "A294_floor/output/tissues.csv", index_col="tissue")[["S", "own"]]
    wf = cs.gtex()
    out = []
    for draw in range(3):
        rng = np.random.default_rng([312, 77, draw])
        pb = cohort(t, wf, rng)
        comps = sorted(pb.compartment.unique())
        ds = cs.Dataset(pb, anchors=comps, waveform=wf, name=f"calib{draw}")
        for k in TARGETS:
            e = E.amplitude_elasticity(ds, k, n_splits=10, seed=[312, draw])
            out.append(dict(draw=draw, tissue=k, S=t.S[k], elasticity=e["elasticity"], leverage=e["leverage"]))
            print(out[-1], flush=True)
    res = pd.DataFrame(out)
    res.to_csv(OUT / "calib_level.csv", index=False)
    print(res.groupby("tissue")[["S", "elasticity", "leverage"]].mean().sort_values("S").round(3).to_string())
    print("median elasticity", round(res.elasticity.median(), 3), "median leverage", round(res.leverage.median(), 3))


if __name__ == "__main__":
    main()
