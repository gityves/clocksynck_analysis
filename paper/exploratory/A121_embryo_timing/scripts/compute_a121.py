"""A121: the mouse embryo clock against the known timing of fetal clock emergence (Qiu 2024, E15.0-E18.75).

Registered in paper/docs/registered_experiments.md (A121, before scoring). Seed 121.
  H121.1  R(E17.0-E18.75) - R(E15.0-E16.75): R = |mean exp(i(phi - w t))| of the host phase of 8 embryos per half (host
          refitted per half on its cell types, batch intercepts), clock times t permuted within batch in both halves
          (10,000), one-sided p.  Secondary (lineages.csv): the same on the clock-plane state of Hepatocytes and
          Mesoderm (centred per batch), descriptive.
  H121.2  common offset of the late embryo host (host phase minus clock time, `Dataset.timed`) against the frame offset
          of the timed mouse datasets (A21 colon C-DC and R-DC, A48 ages, A72 skin): inverse-variance weighted circular
          mean of the per-dataset offsets, SE by bootstrap over datasets and samples; predicted ZT of the embryo clock
          peak for lights-on 06:00, 07:00, 08:00; pass if the 95 % interval under 07:00 overlaps ZT12-ZT14;
          uninformative if the calibration SE > 1.5 h or the datasets differ by > 3 h.

    python compute_a121.py [--n-boot B] [--n-perm P]
Outputs: ../output/{h121,frame_offsets,lineages}.csv
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "paper"))
sys.path.insert(0, str(REPO))
import clocksync as cs  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import data  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
SEED = 121
W = 2 * np.pi / 24
HALVES = {"early": (15.0, 16.75), "late": (17.0, 18.75)}
STROMAL_VASCULAR = ("fibroblast", "smooth_muscle", "pericyte", "endothelial", "lymphatic_EC")
PERIPH = ["adrenal", "heart", "kidney", "lung", "muscle"]


def _load(rel, name):
    spec = importlib.util.spec_from_file_location(name, REPO / "paper" / "results" / rel)
    m = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(Path(spec.origin).parent))
    spec.loader.exec_module(m)
    return m


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def cmean_h(h, w=None):
    return float(np.angle(np.average(np.exp(1j * W * np.asarray(h)), weights=w)) / W) % 24


def host_phase(ds) -> pd.Series:
    """Host phase (rad) per donor, as `Dataset.timed` computes it (clocksync's own host call)."""
    h = A.host(ds.data, ds.anchors)
    return pd.Series(np.angle(h.f[:, 0] + 1j * h.f[:, 1]), index=h.resections)


def offset_of(phi: pd.Series, t: pd.Series) -> float:
    tt = t.reindex(phi.index.astype(str).str.split("#").str[0]).to_numpy(float)
    ok = np.isfinite(tt)
    return cmean_h(np.angle(np.exp(1j * (phi.to_numpy()[ok] - W * tt[ok]))) / W)


def boot_offsets(ds, t: pd.Series, n_boot: int, rng) -> np.ndarray:
    out = []
    for _ in range(n_boot):
        try:
            out.append(offset_of(host_phase(ds.resample(rng)), t))
        except (ValueError, np.linalg.LinAlgError):
            pass
    return np.array(out)


def se_of(point: float, boots: np.ndarray) -> float:
    return float(np.sqrt(np.mean(wrap(boots - point) ** 2)))


# ---- frame offsets of the timed mouse datasets ----------------------------------------------------------------------

def frame_offsets(n_boot: int, rng):
    rows, boots = [], {}

    def add(src, name, point, b, n):
        rows.append(dict(source=src, dataset=name, offset_h=point, se_h=se_of(point, b), n=n, n_boot=len(b)))
        boots[name] = b

    # A21: colon, arms C-DC and R-DC (independent mice; R-PC shares the mice of R-DC), anchors as in A21
    pb, donors = data.load("colon_chrono")
    for arm in ("C-DC", "R-DC"):
        p = pb[pb.donor.isin(donors.index[donors.arm == arm])]
        anchors = [k for k in sorted(p.compartment.unique()) if k in STROMAL_VASCULAR]
        ds = cs.Dataset(p, anchors, donors=donors, name=arm)
        t = donors.zt.astype(float)
        r = ds.timed(t, n_perm=1, seed=SEED)
        add("A21", f"A21 {arm}", float(r.offset_h), boot_offsets(ds, t, n_boot, rng), int(r.n))
    # A48: ages, host on the five peripheral tissues, ZT = circadian time folded to 24 h
    M = _load("figure_4/scripts/resilience/compute_mouse_ageing.py", "a48")
    (pb, donors), mod = M.load_wolff()
    q = mod.tissue_qc(pb)
    pb = pb[pb.compartment.isin(q.index[q.kept])]
    pb, _ = mod.drop_units(pb, PERIPH)
    age = donors.age.reindex(sorted(pb.donor.unique()))
    ct = (donors.ct % 24).reindex(age.index).astype(float)
    for a in M.AGES:
        ds = cs.Dataset(pb[pb.donor.isin(age.index[age == a])], anchors=PERIPH, name=a)
        r = ds.timed(ct, n_perm=1, seed=SEED)
        add("A48", f"A48 {a}", float(r.offset_h), boot_offsets(ds, ct, n_boot, rng), int(r.n))
    # A72: skin datasets, one compartment each: the state itself is the host (offset = arg a / w of the
    # regression of the complex state on e^{iwt}, as in the A72 script); bootstrap within ZT
    K = _load("figure_1/scripts/species/compute_mouse_skin_timed.py", "a72")
    T, _ = data._builder("mouse_skin_timed").build()
    for name in K.PRIMARY:
        d = T[T.cohort == name].reset_index(drop=True)
        c, zt = K.dataset_state(d), d.zt.to_numpy(float)
        point = K.offset(c, zt) % 24
        groups = [np.flatnonzero(zt == z) for z in np.unique(zt)]
        b = np.array([K.offset(c[i], zt[i]) % 24 for i in
                      (np.concatenate([rng.choice(g, len(g)) for g in groups]) for _ in range(n_boot))])
        add("A72", f"A72 {name}", point, b, len(d))
    return pd.DataFrame(rows), boots


def calibrate(fo: pd.DataFrame, boots: dict, n_boot: int, rng):
    """Inverse-variance weighted circular mean of the dataset offsets; SE by bootstrap over datasets and samples."""
    w = 1 / fo.se_h.to_numpy() ** 2
    point = cmean_h(fo.offset_h.to_numpy(), w)
    names = fo.dataset.tolist()
    reps = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(names), len(names))
        o = [rng.choice(boots[names[j]]) for j in pick]
        reps.append(cmean_h(o, w[pick]))
    reps = np.array(reps)
    off = fo.offset_h.to_numpy()
    spread = float(np.max(np.abs(wrap(off[:, None] - off[None, :]))))
    return point, se_of(point, reps), reps, spread


# ---- H121.1 ----------------------------------------------------------------------------------------------------------

def half_host(pb, donors, lo, hi):
    win, don = S.window(pb, donors, lo, hi)
    ds = cs.Dataset(win, anchors=sorted(win.compartment.unique()), donors=donors, name=f"E{lo}-E{hi}")
    return ds, win, don


def r_stat(phi, t):
    return np.abs(np.exp(1j * (phi - W * t)).mean(-1))


def paired_perm(phis, ts, blocks, n_perm, rng):
    """R per half and R(late) - R(early) with clock times permuted within batch in both halves (same draws)."""
    obs = {k: float(r_stat(phis[k], ts[k])) for k in phis}
    null = {k: np.empty(n_perm) for k in phis}
    for k in phis:
        parts = [np.flatnonzero(blocks[k] == b) for b in np.unique(blocks[k])]
        for i in range(n_perm):
            t = ts[k].copy()
            for p in parts:
                t[p] = ts[k][rng.permutation(p)]
            null[k][i] = r_stat(phis[k], t)
    d_obs = obs["late"] - obs["early"]
    d_null = null["late"] - null["early"]
    pf = lambda x, o: float((1 + np.sum(x >= o)) / (1 + n_perm))  # noqa: E731
    return obs, {k: pf(null[k], obs[k]) for k in obs}, d_obs, pf(d_null, d_obs), null


def h121_1(pb, donors, n_perm, rng):
    phis, ts, blocks, meta = {}, {}, {}, {}
    for k, (lo, hi) in HALVES.items():
        ds, win, don = half_host(pb, donors, lo, hi)
        phi = host_phase(ds)
        tt = donors.t_clock.reindex(phi.index.astype(str)).to_numpy(float)
        phis[k], ts[k] = phi.to_numpy(), tt
        blocks[k] = donors.batch.reindex(phi.index.astype(str)).to_numpy()
        chk = ds.timed("t_clock", n_perm=10, seed=SEED)
        meta[k] = dict(n=len(phi), cell_types=win.compartment.nunique(), R_clocksync=float(chk.R),
                       offset_h=float(chk.offset_h))
    obs, p_half, d, p, null = paired_perm(phis, ts, blocks, n_perm, rng)
    return obs, p_half, d, p, meta


# ---- lineages (secondary) --------------------------------------------------------------------------------------------

def lineages(n_perm, rng, seed_rows=("Hepatocytes", "Mesoderm")):
    pb, donors = data.load("mouse_qiu", level="trajectory")
    out = []
    states = {}
    for k, (lo, hi) in HALVES.items():
        don = S.embryos(donors, lo, hi)
        t = S.eligible(pb, min_cells=50, donors=don)
        n = t.groupby("compartment").donor.nunique()
        t = t[t.compartment.isin(n.index[n >= 6])].reset_index(drop=True)
        d = S.dataset(t, donors).data
        states[k] = d
    for lin in seed_rows:
        phis, ts, blocks = {}, {}, {}
        for k, d in states.items():
            i = np.flatnonzero(d.compartment == lin)
            r = d.resection[i]
            y = d.beta[i, 0] + 1j * d.beta[i, 1]
            b = donors.batch.reindex(r).to_numpy()
            for g in np.unique(b):
                y[b == g] -= y[b == g].mean()
            phis[k], ts[k], blocks[k] = np.angle(y), donors.t_clock.reindex(r).to_numpy(float), b
        obs, p_half, dd, p, _ = paired_perm(phis, ts, blocks, n_perm, rng)
        out.append(dict(lineage=lin, n_early=len(phis["early"]), n_late=len(phis["late"]), R_early=obs["early"],
                        p_early=p_half["early"], R_late=obs["late"], p_late=p_half["late"], diff=dd, p_diff=p))
    return pd.DataFrame(out)


# ---- H121.2 ----------------------------------------------------------------------------------------------------------

def peak_rows(o_e, se_e, o_c, se_c, spread, lights=(6, 7, 8), ref=(12.0, 14.0)):
    """Predicted ZT of the embryo clock peak.

    Primary (`peak_zt`): the clock peak is host phase 0 of the shared waveform frame, which the embryo reaches at wall
    time -o_e and the adult mouse at ZT -o_c; embryo ZT = -o_e - L (the calibration offset gives the adult reference
    `adult_phase0_zt` = -o_c and cancels from the difference). `alt_D`: literal subtraction of the frame offset,
    o_c - o_e - L. `alt_B`: adult peak assumed at ZT12 and the embryo shifted by its lead over the adult,
    12 + o_c - o_e - L. Interval: centre +- 1.96 sqrt(se_e^2 + se_c^2) (registered), and embryo-only."""
    rows = []
    for L in lights:
        for nm, c in (("A_phase0", -o_e - L), ("D_literal", o_c - o_e - L), ("B_adult_zt12", 12 + o_c - o_e - L)):
            c = c % 24
            hw = 1.96 * np.hypot(se_e, se_c)
            hw_e = 1.96 * se_e
            ov = lambda h: bool(min(abs(wrap(c - z)) for z in np.linspace(ref[0], ref[1], 41)) <= h)  # noqa: E731
            inside = lambda h: any(abs(wrap(c - z)) <= h for z in np.linspace(ref[0], ref[1], 81))  # noqa: E731
            rows.append(dict(test="H121.2", definition=nm, lights_on=f"{L:02d}:00", embryo_offset_h=o_e,
                             embryo_offset_se_h=se_e, frame_offset_h=o_c, frame_offset_se_h=se_c,
                             adult_phase0_zt=(-o_c) % 24, predicted_peak_zt=c, lo=(c - hw) % 24, hi=(c + hw) % 24,
                             halfwidth_h=hw, overlaps_zt12_14=inside(hw), overlaps_embryo_only=inside(hw_e),
                             dataset_spread_h=spread))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--part", default="all", choices=["all", "h1", "h2"])
    a = ap.parse_args()
    limit_threads(4)
    OUT.mkdir(exist_ok=True, parents=True)
    rng = np.random.default_rng(SEED)
    pd.set_option("display.width", 250)
    pb, donors = data.load("mouse_qiu")

    # H121.1
    obs, p_half, d, p, meta = h121_1(pb, donors, a.n_perm, rng)
    h1 = pd.DataFrame([dict(test="H121.1", statistic="R(late) - R(early)", n_early=meta["early"]["n"],
                            n_late=meta["late"]["n"], cell_types_early=meta["early"]["cell_types"],
                            cell_types_late=meta["late"]["cell_types"], R_early=obs["early"], p_R_early=p_half["early"],
                            R_late=obs["late"], p_R_late=p_half["late"], diff=d, p=p,
                            offset_early_h=meta["early"]["offset_h"], offset_late_h=meta["late"]["offset_h"],
                            R_early_clocksync=meta["early"]["R_clocksync"], R_late_clocksync=meta["late"]["R_clocksync"],
                            n_perm=a.n_perm)])
    print(h1.round(4).T.to_string(), flush=True)
    h1.to_csv(OUT / "h121_part1.csv", index=False)

    lin = lineages(a.n_perm, rng)
    lin.to_csv(OUT / "lineages.csv", index=False)
    print(lin.round(4).to_string(index=False), flush=True)

    # H121.2
    fo, boots = frame_offsets(a.n_boot, rng)
    fo.to_csv(OUT / "frame_offsets.csv", index=False)
    print(fo.round(3).to_string(index=False), flush=True)
    o_c, se_c, _, spread = calibrate(fo, boots, a.n_boot, rng)
    ds, win, don = half_host(pb, donors, *S.LATE)
    ds = cs.Dataset(win, anchors=sorted(win.compartment.unique()), donors=donors, name="late")
    chk = ds.timed("t_clock", n_perm=10, seed=SEED)
    o_e = float(chk.offset_h)
    b_e = boot_offsets(ds, donors.t_clock, a.n_boot, rng)
    se_e = se_of(o_e, b_e)
    print(f"embryo late: n={int(chk.n)} R={chk.R:.3f} offset={o_e:.2f} h se={se_e:.2f}; calibrated frame {o_c:.2f} se {se_c:.2f}"
          f" spread {spread:.2f}", flush=True)
    h2 = peak_rows(o_e, se_e, o_c, se_c, spread)
    informative = bool(se_c <= 1.5 and spread <= 3)
    h2["informative"] = informative
    pri = h2[(h2.definition == "A_phase0") & (h2.lights_on == "07:00")].iloc[0]
    verdict = "uninformative" if not informative else ("pass" if pri.overlaps_zt12_14 else "fail")
    h2["verdict_primary"] = verdict
    print(h2.round(2).to_string(index=False), "\nverdict H121.2:", verdict, flush=True)

    # Holm over H121.1-2 (H121.2 is an interval criterion without p: Holm multiplier m = 2 on the smallest p)
    h1["p_holm"] = min(1.0, 2 * float(h1.p.iloc[0]))
    h1["verdict"] = "pass" if (h1.p_holm.iloc[0] < 0.05) else ("pass_unadjusted" if h1.p.iloc[0] < 0.05 else "fail")
    out = pd.concat([h1, h2], ignore_index=True)
    out.to_csv(OUT / "h121.csv", index=False)
    pd.DataFrame(dict(replicate=np.arange(len(b_e)), embryo_offset_h=b_e)).to_csv(OUT / "embryo_offset_boot.csv", index=False)


if __name__ == "__main__":
    main()
