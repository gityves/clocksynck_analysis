"""A369 (registered): do living patients share a clock amplitude across liver, muscle and visceral fat at matched phase?

STAGE (GSE40231, RMA from CEL; `prep_stage.py`), 66 bypass patients. A255 `Cohort` states (GTEx waveform, centred within
tissue). Q per tissue pair = A368 `Q`; primary Q3 = mean over the liver / muscle / fat pairs. Timing-only nulls over phase
spread and shape at the observed geometry, volume simulations, companions (expression-matched, responsive).

    python a369_stage.py --stage power   # gates, pairings, fractions, geometry, simulations, companions; no clock Q
    python a369_stage.py --stage score   # registered test
Outputs (`../output/`): gates.csv, power.csv, sim.csv, companions.csv | tests.csv, drops.csv.
"""
import argparse
import itertools
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A366_muther_timing" / "scripts"))
sys.path.insert(0, str(EXP / "A368_muther_radial" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a366_muther_timing as M  # noqa: E402
from a368_muther_radial import Q  # noqa: E402
from clocksync import elasticity as E  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
PREP = OUT / "prep"
PRIMARY = ["liver", "muscle", "visceral_fat"]
ARTERY = ["aortic_wall", "mammary_artery"]
N_SIM, N_PERM, N_SETS = 500, 2000, 200
SD_GRID = np.linspace(0.05, 3.0, 60)
TRIO = ["PER1", "PER2", "NR1D1"]


def load():
    A = M.a255()
    x = pd.read_parquet(PREP / "stage_expr.parquet").astype(float) * np.log(2.0)
    s = pd.read_csv(PREP / "stage_sheet.csv", dtype=str).set_index("sample")
    t = {}
    for tis in PRIMARY + ARTERY:
        v = x[(s.tissue == tis).to_numpy()]
        v.index = s.person[v.index].to_numpy()
        t[tis] = v
    return A, A.Cohort("stage", t)


def Q3(Z, idx, tissues=PRIMARY):
    return float(np.mean([Q(Z[a][idx], Z[b][idx]) for a, b in itertools.combinations(tissues, 2)]))


def combined(Z, idx, tissues=PRIMARY):
    ref = Z[tissues[0]][idx]
    out = 0
    for k in tissues:
        z = Z[k][idx]
        out = out + z * np.exp(1j * np.angle(np.sum(ref * np.conj(z)))) / np.abs(z).mean()
    return out


def sim(n, sd, sigma, f, rng, shape="normal"):
    A_ = np.exp(rng.normal(0, sigma, n) - sigma ** 2 / 2)
    ph = rng.normal(0, sd, n) if shape == "normal" else rng.uniform(-sd * np.sqrt(3), sd * np.sqrt(3), n)
    T = A_ * np.exp(1j * ph)
    T = T - T.mean()
    v = np.mean(np.abs(T) ** 2)
    return {k: T + np.sqrt(v * (1 - fk) / fk / 2) * (rng.normal(size=n) + 1j * rng.normal(size=n)) for k, fk in f.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    A, c = load()
    Z = c.Z()
    have = {k: np.isfinite(Z[k]) for k in Z}
    idx = have["liver"] & have["muscle"] & have["visceral_fat"]
    idx5 = idx & have["aortic_wall"] & have["mammary_artery"]
    n, n5 = int(idx.sum()), int(idx5.sum())
    print(f"{n} patients with the primary triad, {n5} with five tissues, {len(c.wf.genes)} clock genes", flush=True)

    if a.stage == "power":
        g = pd.DataFrame([dict(tissue=k, persons=int(have[k].sum()), **c.gate(k, 20 + i)) for i, k in enumerate(PRIMARY + ARTERY)])
        g.to_csv(OUT / "gates.csv", index=False)
        print(g.round(4).to_string(index=False), flush=True)
        prim = [k for k in PRIMARY if float(g.set_index("tissue").p[k]) < 0.05]
        if len(prim) < 3:
            pd.Series(dict(gate="False", note=f"primary gates passed: {prim}")).to_csv(OUT / "power.csv", header=["value"])
            raise SystemExit(f"A369: a primary gate failed ({prim}); not decidable")
        r = {(x, y): A.r_pair(Z[x][idx], Z[y][idx])[0] for x, y in itertools.combinations(PRIMARY, 2)}
        rr = lambda x, y: r[(x, y)] if (x, y) in r else r[(y, x)]  # noqa: E731
        f = {}
        for k in PRIMARY:
            o = [j for j in PRIMARY if j != k]
            f[k] = float(min(rr(k, o[0]) * rr(k, o[1]) / rr(o[0], o[1]), 0.99))
        target = M.elongation(combined(Z, idx))
        rng = np.random.default_rng([369, 1])
        el = np.array([np.median([M.elongation(combined(s, np.ones(n, bool))) for s in (sim(n, sd, 0.0, f, rng) for _ in range(20))])
                       for sd in SD_GRID])
        sd_hat = float(SD_GRID[np.argmin(np.abs(np.log(el) - np.log(target)))])
        cols = {}
        all_ = np.ones(n, bool)
        for k in (0.75, 1.0, 1.25):
            cols[f"null_normal_{k}"] = [Q3(sim(n, k * sd_hat, 0.0, f, rng), all_) for _ in range(N_SIM)]
        cols["null_uniform_1.0"] = [Q3(sim(n, sd_hat, 0.0, f, rng, "uniform"), all_) for _ in range(N_SIM)]
        cols["alt_gtex"] = [Q3(sim(n, sd_hat, M.SIGMA_A, f, rng), all_) for _ in range(N_SIM)]
        cols["alt_half"] = [Q3(sim(n, sd_hat, M.SIGMA_A / 2, f, rng), all_) for _ in range(N_SIM)]
        s = pd.DataFrame(cols)
        s.to_csv(OUT / "sim.csv", index=False)
        rows = []
        for kind, tag in (("matched", 3), ("responsive", 4)):
            for j, w in enumerate(c.pseudo(kind, N_SETS, tag)):
                rows.append(dict(kind=kind, set=j, Q3=Q3(c.Z(w), idx)))
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "companions.csv", index=False)
        nq = {k: float(s[k].quantile(0.95)) for k in s if k.startswith("null")}
        mq95 = float(comp[comp.kind == "matched"].Q3.quantile(0.95))
        thr = max(mq95, max(nq.values()))
        res = dict(persons=n, persons_5=n5, **{f"pairing_{x}_{y}": v for (x, y), v in r.items()},
                   **{f"fraction_{k}": v for k, v in f.items()}, elongation_combined=target, phase_sd_rad=sd_hat,
                   phase_sd_h=sd_hat * 24 / (2 * np.pi), **{f"{k}_median": float(s[k].median()) for k in s},
                   **{f"{k}_q95": v for k, v in nq.items()}, matched_median=float(comp[comp.kind == "matched"].Q3.median()),
                   matched_q95=mq95, responsive_median=float(comp[comp.kind == "responsive"].Q3.median()),
                   responsive_q95=float(comp[comp.kind == "responsive"].Q3.quantile(0.95)), threshold=thr,
                   power_gtex=float(np.mean(s.alt_gtex > thr)), power_half=float(np.mean(s.alt_half > thr)))
        res["gate"] = bool(res["power_gtex"] >= 0.8)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("A369 not powered: not decidable at this design; not scored")
    comp = pd.read_csv(OUT / "companions.csv")
    s = pd.read_csv(OUT / "sim.csv")
    q = Q3(Z, idx)
    rng = np.random.default_rng([369, 2])
    sub = {k: Z[k][idx] for k in PRIMARY}
    null = np.array([Q3({k: v[rng.permutation(n)] for k, v in sub.items()}, np.ones(n, bool)) for _ in range(N_PERM)])
    mt, rs = comp[comp.kind == "matched"].Q3, comp[comp.kind == "responsive"].Q3
    row = dict(test="H369.1", persons=n, Q3=q, p_perm=float((1 + np.sum(null >= q)) / (1 + N_PERM)),
               threshold=float(pw["threshold"]), pct_timing_only=float((s["null_normal_1.0"] < q).mean() * 100),
               pct_alt_half=float((s.alt_half < q).mean() * 100), pct_alt_gtex=float((s.alt_gtex < q).mean() * 100),
               matched_q95=float(mt.quantile(.95)), pct_matched=float((mt < q).mean() * 100),
               responsive_q95=float(rs.quantile(.95)), pct_responsive=float((rs < q).mean() * 100))
    for x, y in itertools.combinations(PRIMARY, 2):
        row[f"Q_{x}_{y}"] = Q(Z[x][idx], Z[y][idx])
        row[f"leverage_raw_{x}_{y}"] = E.radius_covariance(Z[x][idx], Z[y][idx])["leverage"]
    row["passes"] = bool(q > row["threshold"] and row["p_perm"] < 0.05)
    row["Q10_five_tissues"] = float(np.mean([Q(Z[x][idx5], Z[y][idx5]) for x, y in itertools.combinations(PRIMARY + ARTERY, 2)]))
    core = c.wf.drop(TRIO)
    row["Q3_core"] = Q3(c.Z(core), idx)
    dr = pd.DataFrame([dict(drop="+".join(d), trio=list(d) == TRIO, Q3=Q3(c.Z(c.wf.drop(list(d))), idx))
                       for d in itertools.combinations(c.wf.genes, 3)])
    dr.to_csv(OUT / "drops.csv", index=False)
    row.update(drops_median=float(dr.Q3.median()), trio_drop_pct=float((dr.Q3 < row["Q3_core"]).mean() * 100))
    pd.DataFrame([row]).to_csv(OUT / "tests.csv", index=False)
    print(pd.Series(row).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
