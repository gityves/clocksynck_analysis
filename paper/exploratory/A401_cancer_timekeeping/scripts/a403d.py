"""A403d (POST HOC, DESCRIPTIVE; registered in paper/docs/registered_experiments.md, "A403d"): can shared biopsy time on
hospital-hours arcs produce the colder-louder slope of A403c?

A403c table (patients with immune-hot, their cohorts, strata, cancers). Stratified arc model with no loudness differences;
immune-hot driven by the true phase; through the A403c estimator (phase-free loudness, pooled within-cohort slope).
                                                     -> output/a403d/{geometry,sim}.csv
    python a403d.py [--jobs 12] [--genome-jobs 4]                        (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import a403_a405_score as S
import a403b as B
import a403c as C3
from a402_a405 import U, log, trim, centred, cancer_hosts, loudness_frame, load_slim, readouts
from clocksync.config import limit_threads

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "A366_muther_timing" / "scripts"))
import a366_muther_timing as M6  # noqa: E402

OUT = U.OUT / "a403d"
REL, N_REP = 0.70, 200
R_GRID = (0.3, 0.5, 0.7)
T0_GRID = np.arange(0, 24, 2)
SD_GRID = np.linspace(0.05, 3.0, 40)
H = 2 * np.pi / 24


def sim_states(st, sd, rng):
    """Stratified arc: common centre phase 0, stratum offset N(0, 1 h), patient phases N(offset, sd); unit amplitude;
    noise at reliability REL; centred within stratum. Returns (true phase, observed centred state)."""
    us = np.unique(st)
    off = dict(zip(us, rng.normal(0, H, len(us))))
    ph = np.array([off[s] for s in st]) + rng.normal(0, sd, len(st))
    T = np.exp(1j * ph)
    Tc = centred(T, pd.Series(st))
    v = np.mean(np.abs(Tc) ** 2)
    s = np.sqrt(v * (1 - REL) / REL / 2)
    z = centred(T + s * (rng.normal(size=len(st)) + 1j * rng.normal(size=len(st))), pd.Series(st))
    return ph, z


def slope_of(d, z, y):
    """A403c estimator: phase-free log radius per cancer, pooled within-cohort slope of y (no covariates)."""
    ok = np.abs(z) > 0
    lr = pd.Series(np.where(ok, np.log(np.where(ok, np.abs(z), 1.0)), np.nan), index=d.index)
    lr = pd.Series(centred(lr.to_numpy(), d.stratum), index=d.index)
    pf, _ = C3.phase_free(lr, np.angle(z), d.stratum.to_numpy(), d.cancer.to_numpy())
    dd = d.assign(y=y, lpf=pf.to_numpy())
    un = B.units(dd, "y", [], "lpf")
    return S.W.pooled_slope(un) if un else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--genome-jobs", type=int, default=4)
    args = ap.parse_args()
    args.smoke = False
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    frames, donors, f388 = U.pool_sources(args, OUT)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    p401 = pd.read_csv(U.OUT / "power" / "p401.csv")
    cancers = [c for c in p401[(p401.readout.astype(str) == "14") & p401.p401_pass].cancer if c in set(donors.cancer) and c != "unspecified"]
    donors = donors[donors.cancer.isin(cancers)].reset_index(drop=True)
    hosts = cancer_hosts(frames, donors, cancers)
    sources, f388 = list(frames), set(f388)
    frames.clear()
    trim()
    loud = loudness_frame(hosts)
    slim, errs = load_slim(sources, donors, f388, args)
    assert len(slim) == 51 and errs.empty
    ro = readouts(slim, loud)
    df = loud.merge(ro.drop(columns="cohort"), on="donor").merge(S.covariates(hosts), on="donor", how="left")
    df["cohort"] = df.donor.str.split("|").str[0]
    assert abs(B.slope(df, "immune_hot", B.BASECOV)[0] - (-0.10228)) < 5e-4, "A403 not reproduced"
    d = df.dropna(subset=["immune_hot"]).reset_index(drop=True)
    st = d.stratum.astype(str).to_numpy()
    zc = centred(d.z.to_numpy(), d.stratum)
    obs_el = M6.elongation(zc)
    # same estimator on the observed data without covariates (reference for the simulated slopes)
    obs_slope = slope_of(d, zc, d.immune_hot.to_numpy(float))
    rng = np.random.default_rng([4034])
    el = np.array([np.median([M6.elongation(sim_states(st, sd, rng)[1]) for _ in range(10)]) for sd in SD_GRID])
    sd_hat = float(SD_GRID[np.argmin(np.abs(el - obs_el))])
    geo = pd.DataFrame([dict(patients=len(d), strata=len(set(st)), elongation_observed=obs_el, sd_hat_rad=sd_hat,
                             sd_hat_h=sd_hat / H, observed_slope_no_covariates=obs_slope)])
    geo.to_csv(OUT / "geometry.csv", index=False)
    log("geometry\n" + geo.round(4).to_string(index=False))
    rows = []
    for r in R_GRID:
        for t in T0_GRID:
            vals = []
            for _ in range(N_REP):
                ph, z = sim_states(st, sd_hat, rng)
                y = r * np.cos(ph - t * H) / np.sqrt(0.5) * 1.0
                y = r * (y / y.std()) + np.sqrt(1 - r ** 2) * rng.normal(size=len(y))
                vals.append(slope_of(d, z, y))
            v = np.array(vals)
            rows.append(dict(r=r, t0_h=int(t), median=float(np.nanmedian(v)), q05=float(np.nanquantile(v, .05)),
                             q95=float(np.nanquantile(v, .95)), share_le_observed=float(np.nanmean(v <= -0.126))))
        log(f"r {r} done")
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "sim.csv", index=False)
    pd.set_option("display.width", 200)
    log("sim\n" + T.round(4).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
