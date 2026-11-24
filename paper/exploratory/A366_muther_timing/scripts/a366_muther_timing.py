"""A366 (registered): is A331's living fat-skin loudness a person volume or the shared biopsy time?

MuTHER (A255 states, as A331). Shared biopsy time puts fat and skin at the same place on the clock; on a daytime arc the
centred radius then grows with the distance from the cohort's mean phase. Steps:
- geometry: the clock-state cloud of each tissue and of the combined state; elongation = major / minor eigenvalue of
  the 2 x 2 covariance (about 1 on a ring, large on a narrow arc);
- circle-centred leverage: each tissue's states recentred on its fitted circle (algebraic fit), so the radius is
  loudness and the angle is phase, then `radius_covariance` (radii residualised on two harmonics of both angles);
  a first variant (radii residualised on the other tissue's major-axis position) was dropped before registration: on
  synthetic data it removed a true volume (the position carries amplitude);
- timing-only null and power: persons on a ring with phases spread to reproduce the observed combined elongation, with
  no person amplitude (null) or a GTEx-size one (log SD 0.352), noise at A331's clock fractions;
- companions: the same corrected statistic on A331's 200 matched and 200 responsive sets.

    python a366_muther_timing.py --stage power   # geometry, simulation, companions; no clock radius statistic
    python a366_muther_timing.py --stage score   # registered test and the phase-radius relation
Outputs (`../output/`): power.csv, geometry.csv, companions.csv, sim.csv | tests.csv. Aggregates only.
"""
import argparse
import importlib.util
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync import elasticity as E  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
F_FAT, F_SKIN = 0.696, 0.408   # A331 clock fractions of state variance
SIGMA_A, N_SIM, N_SETS = 0.352, 500, 200
SD_GRID = np.linspace(0.05, 3.0, 60)   # phase SD (rad) grid for matching the elongation


def a255():
    spec = importlib.util.spec_from_file_location("a255", PAPER / "exploratory/A255_living_tissues/scripts/compute_a255.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def axes(z):
    """Eigenvalues (major, minor) and major-axis unit vector of the 2 x 2 covariance of a centred complex cloud."""
    C = np.cov(np.vstack([z.real, z.imag]))
    w, V = np.linalg.eigh(C)
    return w[1], w[0], V[:, 1]


def elongation(z):
    a, b, _ = axes(z)
    return float(a / b)


def tpos(z):
    _, _, v = axes(z)
    return z.real * v[0] + z.imag * v[1]


def circle_centre(z):
    """Algebraic (Kasa) circle fit: centre of the ring the states lie on (the cohort mean is not the centre on an arc)."""
    A = np.column_stack([z.real, z.imag, np.ones(len(z))])
    c = np.linalg.lstsq(A, z.real ** 2 + z.imag ** 2, rcond=None)[0]
    return c[0] / 2 + 1j * c[1] / 2


def corrected_leverage(a, b):
    """Each tissue's states recentred on its fitted circle, so that the radius is loudness and the angle is phase; then
    `radius_covariance` (radii residualised on two harmonics of both angles)."""
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    return float(E.radius_covariance(a - circle_centre(a), b - circle_centre(b))["leverage"])


def sim_pair(n, sd, sigma, rng):
    """Fat and skin states for persons on a ring: phases N(0, sd), log amplitude SD sigma, A331 noise fractions."""
    A = np.exp(rng.normal(0, sigma, n) - sigma ** 2 / 2)
    T = A * np.exp(1j * rng.normal(0, sd, n))
    T = T - T.mean()
    v = np.mean(np.abs(T) ** 2)
    out = []
    for f in (F_FAT, F_SKIN):
        s = np.sqrt(v * (1 - f) / f / 2)
        out.append(T + s * (rng.normal(size=n) + 1j * rng.normal(size=n)))
    return out


def combined(zf, zs):
    rot = np.angle(np.nansum(zf * np.conj(zs)))
    return zf / np.abs(zf).mean() + zs * np.exp(1j * rot) / np.abs(zs).mean()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    m = a255()
    c, persons = m.load_muther(False)
    both = set(persons[(persons.has_fat == "True") & (persons.has_skin == "True")].person)
    idx = np.array([i in both for i in c.ids])
    Z = c.Z()
    zf, zs = Z["fat"][idx], Z["skin"][idx]
    n = int(idx.sum())

    if a.stage == "power":
        geo = pd.DataFrame([dict(cloud=k, elongation=elongation(z)) for k, z in
                            (("fat", zf), ("skin", zs), ("combined", combined(zf, zs)))])
        geo.to_csv(OUT / "geometry.csv", index=False)
        target = float(geo.set_index("cloud").elongation["combined"])
        rng = np.random.default_rng([366, 1])
        el = np.array([np.median([elongation(combined(*sim_pair(n, sd, 0.0, rng))) for _ in range(20)]) for sd in SD_GRID])
        sd_hat = float(SD_GRID[np.argmin(np.abs(np.log(el) - np.log(target)))])
        null = np.array([corrected_leverage(*sim_pair(n, sd_hat, 0.0, rng)) for _ in range(N_SIM)])
        null_raw = np.array([E.radius_covariance(*sim_pair(n, sd_hat, 0.0, rng))["leverage"] for _ in range(N_SIM)])
        alt = np.array([corrected_leverage(*sim_pair(n, sd_hat, SIGMA_A, rng)) for _ in range(N_SIM)])
        alt_half = np.array([corrected_leverage(*sim_pair(n, sd_hat, SIGMA_A / 2, rng)) for _ in range(N_SIM)])
        rows = []
        for kind, tag in (("matched", 3), ("responsive", 4)):
            for j, w in enumerate(c.pseudo(kind, m.N_SETS, tag)):
                Zp = c.Z(w)
                rows.append(dict(kind=kind, set=j, corrected=corrected_leverage(Zp["fat"][idx], Zp["skin"][idx])))
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "companions.csv", index=False)
        mq95 = float(comp[comp.kind == "matched"].corrected.quantile(0.95))
        thr = max(mq95, float(np.quantile(null, 0.95)))
        res = dict(persons=n, elongation_fat=float(geo.elongation[0]), elongation_skin=float(geo.elongation[1]),
                   elongation_combined=target, phase_sd_rad=sd_hat, phase_sd_h=sd_hat * 24 / (2 * np.pi),
                   timing_only_raw_leverage_median=float(np.median(null_raw)),
                   timing_only_raw_leverage_q95=float(np.quantile(null_raw, 0.95)),
                   timing_only_corrected_median=float(np.median(null)), timing_only_corrected_q95=float(np.quantile(null, 0.95)),
                   matched_median=float(comp[comp.kind == "matched"].corrected.median()), matched_q95=mq95,
                   responsive_median=float(comp[comp.kind == "responsive"].corrected.median()),
                   responsive_q95=float(comp[comp.kind == "responsive"].corrected.quantile(0.95)),
                   threshold=thr, expected_corrected_gtex=float(np.median(alt)), power_gtex=float(np.mean(alt > thr)),
                   power_half=float(np.mean(alt_half > thr)))
        res["gate"] = bool(res["power_gtex"] >= 0.8)
        pd.DataFrame(dict(null=null, null_raw=null_raw, alt=alt, alt_half=alt_half)).to_csv(OUT / "sim.csv", index=False)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("A366 not powered: not decidable at this design; not scored")
    comp = pd.read_csv(OUT / "companions.csv")
    L = corrected_leverage(zf, zs)
    rng = np.random.default_rng([366, 2])
    nullp = np.array([corrected_leverage(zf, zs[rng.permutation(n)]) for _ in range(2000)])
    mt, rs = comp[comp.kind == "matched"].corrected, comp[comp.kind == "responsive"].corrected
    from scipy.stats import spearmanr
    rf, rk = np.abs(zf), np.abs(zs)
    row = dict(test="H366.1", persons=n, corrected_leverage=L, raw_leverage_a331=E.radius_covariance(zf, zs)["leverage"],
               p_perm=float((1 + np.sum(nullp >= L)) / 2001), threshold=float(pw["threshold"]),
               timing_only_corrected_q95=float(pw["timing_only_corrected_q95"]), matched_q95=float(mt.quantile(.95)),
               pct_matched=float((mt < L).mean() * 100), responsive_q95=float(rs.quantile(.95)),
               pct_responsive=float((rs < L).mean() * 100),
               rho_radius_fat_vs_abs_t_skin=float(spearmanr(rf, np.abs(tpos(zs))).statistic),
               rho_radius_skin_vs_abs_t_fat=float(spearmanr(rk, np.abs(tpos(zf))).statistic),
               centre_shift_fat=float(abs(circle_centre(zf)) / np.abs(zf).mean()),
               centre_shift_skin=float(abs(circle_centre(zs)) / np.abs(zs).mean()))
    row["passes"] = bool(L > row["threshold"] and row["p_perm"] < 0.05)
    row["reading"] = ("living fat and skin share loudness beyond the shared biopsy time and random genes" if row["passes"]
                      else "A331 not separated from the shared biopsy time")
    pd.DataFrame([row]).to_csv(OUT / "tests.csv", index=False)
    print(pd.Series(row).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
