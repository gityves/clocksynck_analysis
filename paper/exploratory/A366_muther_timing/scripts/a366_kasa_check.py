"""A366 post hoc check (2026-11-24, on the user's remark that a circle centre is hard to find on a narrow arc): how far
does the algebraic (Kasa) circle fit move the centre towards the true ring centre at MuTHER's phase spread and noise?

Synthetic only (no MuTHER value): persons on a ring of radius 1 (no person amplitude), phases N(0, sd), centred on the
cohort mean, complex noise at A331's skin clock fraction (0.408) or none; 680 persons, 200 draws per row (seed 1).
Reports the true centre's distance from the cohort mean (relative to the mean centred radius) and the fit's remaining
error (relative to that distance: 0 = centre found, 1 = fit stayed at the cohort mean).

    python a366_kasa_check.py
Output: ../output/kasa_check.csv.
"""
import numpy as np
import pandas as pd

from a366_muther_timing import F_SKIN, OUT, circle_centre

pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
N, SD_OBS = int(float(pw["persons"])), float(pw["phase_sd_rad"])


def main():
    rng = np.random.default_rng(1)
    rows = []
    for sd in (SD_OBS, 1.0, 1.5, 2.5):
        for noisy in (True, False):
            dist, err = [], []
            for _ in range(200):
                T = np.exp(1j * rng.normal(0, sd, N))
                mu = T.mean()
                T = T - mu
                if noisy:
                    s = np.sqrt(np.mean(np.abs(T) ** 2) * (1 - F_SKIN) / F_SKIN / 2)
                    T = T + s * (rng.normal(size=N) + 1j * rng.normal(size=N))
                dist.append(abs(mu) / np.abs(T).mean())
                err.append(abs(circle_centre(T) + mu) / abs(mu))
            rows.append(dict(phase_sd_rad=sd, phase_sd_h=sd * 24 / (2 * np.pi), noise="A331 skin" if noisy else "none",
                             true_centre_distance=float(np.median(dist)), fit_error=float(np.median(err))))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "kasa_check.csv", index=False)
    print(t.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
