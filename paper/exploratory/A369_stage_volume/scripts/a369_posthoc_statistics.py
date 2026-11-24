"""A369 post hoc (2026-11-24, synthetic only, no STAGE clock value): would a radius-based statistic decide at STAGE's
geometry and clock fractions? Pairwise raw leverage and one tissue vs the mean of the other two, timing-only null vs
GTEx-size and half-size volumes; power of the latter at 66, 200 and 600 persons (timing null only, no companions).
Run from paper/: python exploratory/A369_stage_volume/scripts/a369_posthoc_statistics.py (prints only)."""
import sys, itertools, numpy as np, pandas as pd
sys.path.insert(0, "exploratory/A369_stage_volume/scripts")
import a369_stage as S
from clocksync import elasticity as E
pw = pd.read_csv("exploratory/A369_stage_volume/output/power.csv", index_col=0)["value"]
f = {k: float(pw[f"fraction_{k}"]) for k in S.PRIMARY}; sd = float(pw["phase_sd_rad"]); n = 66
def lev3(Z):  # mean pairwise raw radius leverage (A331 statistic)
    return np.mean([E.radius_covariance(Z[a], Z[b])["leverage"] for a, b in itertools.combinations(S.PRIMARY, 2)])
def lev_split(Z):  # one tissue vs the mean of the other two (GTEx-like half host), averaged
    out = []
    for k in S.PRIMARY:
        o = [j for j in S.PRIMARY if j != k]
        out.append(E.radius_covariance(Z[k], (Z[o[0]] + Z[o[1]]) / 2)["leverage"])
    return np.mean(out)
rng = np.random.default_rng(5)
for name, fn in (("pairwise raw leverage", lev3), ("one vs other two", lev_split)):
    for k in (1.0, 1.25):
        nul = [fn(S.sim(n, k * sd, 0, f, rng)) for _ in range(400)]
        alt = [fn(S.sim(n, sd, S.M.SIGMA_A, f, rng)) for _ in range(400)]
        half = [fn(S.sim(n, sd, S.M.SIGMA_A / 2, f, rng)) for _ in range(400)]
        q = np.quantile(nul, .95)
        print(f"{name:22s} null sd x{k}: med {np.median(nul):.3f} q95 {q:.3f} | GTEx med {np.median(alt):.3f} power {np.mean(np.array(alt) > q):.2f} | half power {np.mean(np.array(half) > q):.2f}")
for n_ in (66, 200, 600):
    nul = [lev_split(S.sim(n_, sd, 0, f, rng)) for _ in range(300)]; alt = [lev_split(S.sim(n_, sd, S.M.SIGMA_A, f, rng)) for _ in range(300)]
    print(f"n {n_}: one-vs-two power (timing only, no companions) {np.mean(np.array(alt) > np.quantile(nul, .95)):.2f}")
