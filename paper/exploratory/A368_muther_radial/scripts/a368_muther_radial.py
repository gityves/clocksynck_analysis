"""A368 (registered): at matched phase, do living women with high fat clock amplitude also have high skin amplitude?

MuTHER (A255 states, as A331 / A366). On a daytime arc, position along the arc reads phase and the perpendicular
direction reads amplitude. Skin is rotated onto fat; the combined cloud's major axis gives x (along the arc) and y
(perpendicular, oriented away from the ring centre) for both tissues; Q = correlation of y_fat and y_skin after each is
regressed on a quadratic of both arc positions. Timing-only nulls over phase spread and shape, GTEx-size and half-size
volumes, companions (A331 sets).

    python a368_muther_radial.py --stage power   # simulations and companions; no clock Q
    python a368_muther_radial.py --stage score   # registered test
Outputs (`../output/`): power.csv, sim.csv, companions.csv | tests.csv. Aggregates only.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "A366_muther_timing" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a366_muther_timing as M  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
A366_POWER = HERE.parents[2] / "A366_muther_timing" / "output" / "power.csv"
N_SIM, N_PERM = 500, 2000


def frame(zf, zs):
    """x (along the arc) and y (perpendicular, away from the ring centre) of both tissues in a common frame."""
    zs = zs * np.exp(1j * np.angle(np.nansum(zf * np.conj(zs))))
    c = M.combined(zf, zs)
    _, _, v = M.axes(c)
    u = np.conj(v[0] + 1j * v[1])
    wc, wf, ws = c * u, zf * u, zs * u
    X = np.column_stack([np.ones(len(wc)), wc.real, wc.real ** 2])
    if np.linalg.lstsq(X, wc.imag, rcond=None)[0][2] > 0:
        wf, ws = np.conj(wf), np.conj(ws)
    return wf.real, wf.imag, ws.real, ws.imag


def Q(zf, zs, perm=None):
    ok = np.isfinite(zf) & np.isfinite(zs)
    zf, zs = zf[ok], zs[ok]
    xf, yf, xs, ys = frame(zf, zs)
    X = np.column_stack([np.ones(len(xf)), xf, xs, xf ** 2, xs ** 2, xf * xs])
    res = lambda y: y - X @ np.linalg.lstsq(X, y, rcond=None)[0]  # noqa: E731
    rf, rs = res(yf), res(ys)
    if perm is not None:
        rs = rs[perm]
    return float(np.corrcoef(rf, rs)[0, 1])


def sim(n, sd, sigma, rng, shape="normal"):
    """A366 `sim_pair` with a choice of phase distribution (uniform with the same SD)."""
    A = np.exp(rng.normal(0, sigma, n) - sigma ** 2 / 2)
    ph = rng.normal(0, sd, n) if shape == "normal" else rng.uniform(-sd * np.sqrt(3), sd * np.sqrt(3), n)
    T = A * np.exp(1j * ph)
    T = T - T.mean()
    v = np.mean(np.abs(T) ** 2)
    out = []
    for f in (M.F_FAT, M.F_SKIN):
        s = np.sqrt(v * (1 - f) / f / 2)
        out.append(T + s * (rng.normal(size=n) + 1j * rng.normal(size=n)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    m = M.a255()
    c, persons = m.load_muther(False)
    both = set(persons[(persons.has_fat == "True") & (persons.has_skin == "True")].person)
    idx = np.array([i in both for i in c.ids])
    Z = c.Z()
    zf, zs = Z["fat"][idx], Z["skin"][idx]
    n = int(idx.sum())
    sd = float(pd.read_csv(A366_POWER, index_col=0)["value"]["phase_sd_rad"])

    if a.stage == "power":
        rng = np.random.default_rng([368, 1])
        cols = {}
        for k in (0.75, 1.0, 1.25):
            cols[f"null_normal_{k}"] = [Q(*sim(n, k * sd, 0.0, rng)) for _ in range(N_SIM)]
        cols["null_uniform_1.0"] = [Q(*sim(n, sd, 0.0, rng, "uniform")) for _ in range(N_SIM)]
        cols["alt_gtex"] = [Q(*sim(n, sd, M.SIGMA_A, rng)) for _ in range(N_SIM)]
        cols["alt_half"] = [Q(*sim(n, sd, M.SIGMA_A / 2, rng)) for _ in range(N_SIM)]
        s = pd.DataFrame(cols)
        s.to_csv(OUT / "sim.csv", index=False)
        rows = []
        for kind, tag in (("matched", 3), ("responsive", 4)):
            for j, w in enumerate(c.pseudo(kind, m.N_SETS, tag)):
                Zp = c.Z(w)
                rows.append(dict(kind=kind, set=j, Q=Q(Zp["fat"][idx], Zp["skin"][idx])))
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "companions.csv", index=False)
        nq = {k: float(s[k].quantile(0.95)) for k in s if k.startswith("null")}
        mq95 = float(comp[comp.kind == "matched"].Q.quantile(0.95))
        thr = max(mq95, max(nq.values()))
        res = dict(persons=n, phase_sd_rad=sd, **{f"{k}_median": float(s[k].median()) for k in s},
                   **{f"{k}_q95": v for k, v in nq.items()}, matched_median=float(comp[comp.kind == "matched"].Q.median()),
                   matched_q95=mq95, responsive_median=float(comp[comp.kind == "responsive"].Q.median()),
                   responsive_q95=float(comp[comp.kind == "responsive"].Q.quantile(0.95)), threshold=thr,
                   power_gtex=float(np.mean(s.alt_gtex > thr)), power_half=float(np.mean(s.alt_half > thr)))
        res["gate"] = bool(res["power_gtex"] >= 0.8)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("A368 not powered: not decidable at this design; not scored")
    comp = pd.read_csv(OUT / "companions.csv")
    s = pd.read_csv(OUT / "sim.csv")
    q = Q(zf, zs)
    rng = np.random.default_rng([368, 2])
    m_ok = int((np.isfinite(zf) & np.isfinite(zs)).sum())
    null = np.array([Q(zf, zs, rng.permutation(m_ok)) for _ in range(N_PERM)])
    mt, rs = comp[comp.kind == "matched"].Q, comp[comp.kind == "responsive"].Q
    row = dict(test="H368.1", persons=n, Q=q, p_perm=float((1 + np.sum(null >= q)) / (1 + N_PERM)),
               threshold=float(pw["threshold"]), pct_timing_only=float((s["null_normal_1.0"] < q).mean() * 100),
               pct_alt_half=float((s.alt_half < q).mean() * 100), pct_alt_gtex=float((s.alt_gtex < q).mean() * 100),
               matched_q95=float(mt.quantile(.95)), pct_matched=float((mt < q).mean() * 100),
               responsive_q95=float(rs.quantile(.95)), pct_responsive=float((rs < q).mean() * 100))
    row["passes"] = bool(q > row["threshold"] and row["p_perm"] < 0.05)
    pd.DataFrame([row]).to_csv(OUT / "tests.csv", index=False)
    print(pd.Series(row).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
