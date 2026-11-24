"""A370 (registered): do living adults share a clock amplitude between subcutaneous fat and muscle at matched phase?

AAGMEx (GSE95674 / GSE95675; `prep_aagmex.py`), 251 persons with both tissues. A255 `Cohort` states (GTEx waveform on the
genes present, centred within tissue); Q = A368 `Q`. Timing-only nulls over phase spread, shape and clock fractions at
the observed geometry, volume simulations, companions (expression-matched, responsive).

    python a370_aagmex.py --stage power   # gates, pairing, fractions, geometry, simulations, companions; no clock Q
    python a370_aagmex.py --stage score   # registered test
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
sys.path.insert(0, str(EXP / "A369_stage_volume" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a366_muther_timing as M  # noqa: E402
from a368_muther_radial import Q  # noqa: E402
from a369_stage import sim  # noqa: E402  (three-tissue-general simulator: dict of fractions)
from clocksync import elasticity as E  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
PREP = OUT / "prep"
T2 = ["adipose", "muscle"]
N_SIM, N_PERM, N_SETS = 500, 2000, 200
SD_GRID = np.linspace(0.05, 3.0, 60)
TRIO = ["PER1", "PER2", "NR1D1"]


def load():
    A = M.a255()
    x = pd.read_parquet(PREP / "aagmex_expr.parquet").astype(float) * np.log(2.0)
    s = pd.read_csv(PREP / "aagmex_sheet.csv", dtype=str).set_index("sample")
    t = {}
    for tis in T2:
        v = x[(s.tissue == tis).to_numpy()]
        v.index = s.person[v.index].to_numpy()
        t[tis] = v
    return A, A.Cohort("aagmex", t)


def q2(Z, idx):
    return Q(Z["adipose"][idx], Z["muscle"][idx])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    A, c = load()
    Z = c.Z()
    idx = np.isfinite(Z["adipose"]) & np.isfinite(Z["muscle"])
    n = int(idx.sum())
    print(f"{n} persons with both tissues, {len(c.wf.genes)} clock genes", flush=True)

    if a.stage == "power":
        g = pd.DataFrame([dict(tissue=k, persons=int(np.isfinite(Z[k]).sum()), **c.gate(k, 30 + i)) for i, k in enumerate(T2)])
        g.to_csv(OUT / "gates.csv", index=False)
        print(g.round(4).to_string(index=False), flush=True)
        if (g.p >= 0.05).any():
            pd.Series(dict(gate="False", note="a tissue gate failed")).to_csv(OUT / "power.csv", header=["value"])
            raise SystemExit("A370: a gate failed; not decidable")
        r = A.r_pair(Z["adipose"][idx], Z["muscle"][idx])[0]
        sh = g.set_index("tissue").share
        ratio = float(sh["adipose"] / sh["muscle"])
        f = {"adipose": float(min(r * np.sqrt(ratio), 0.99)), "muscle": float(min(r / np.sqrt(ratio), 0.99))}
        target = M.elongation(M.combined(Z["adipose"][idx], Z["muscle"][idx]))
        rng = np.random.default_rng([370, 1])
        all_ = np.ones(n, bool)
        comb = lambda s_: M.combined(s_["adipose"], s_["muscle"])  # noqa: E731
        el = np.array([np.median([M.elongation(comb(sim(n, sd, 0.0, f, rng))) for _ in range(20)]) for sd in SD_GRID])
        sd_hat = float(SD_GRID[np.argmin(np.abs(np.log(el) - np.log(target)))])
        feq = {k: float(min(r, 0.99)) for k in T2}
        cols = {}
        for k in (0.75, 1.0, 1.25):
            cols[f"null_normal_{k}"] = [q2(sim(n, k * sd_hat, 0.0, f, rng), all_) for _ in range(N_SIM)]
        cols["null_uniform_1.0"] = [q2(sim(n, sd_hat, 0.0, f, rng, "uniform"), all_) for _ in range(N_SIM)]
        cols["null_equal_fractions_1.0"] = [q2(sim(n, sd_hat, 0.0, feq, rng), all_) for _ in range(N_SIM)]
        cols["alt_gtex"] = [q2(sim(n, sd_hat, M.SIGMA_A, f, rng), all_) for _ in range(N_SIM)]
        cols["alt_half"] = [q2(sim(n, sd_hat, M.SIGMA_A / 2, f, rng), all_) for _ in range(N_SIM)]
        s = pd.DataFrame(cols)
        s.to_csv(OUT / "sim.csv", index=False)
        rows = []
        for kind, tag in (("matched", 3), ("responsive", 4)):
            for j, w in enumerate(c.pseudo(kind, N_SETS, tag)):
                rows.append(dict(kind=kind, set=j, Q=q2(c.Z(w), idx)))
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "companions.csv", index=False)
        nq = {k: float(s[k].quantile(0.95)) for k in s if k.startswith("null")}
        mq95 = float(comp[comp.kind == "matched"].Q.quantile(0.95))
        thr = max(mq95, max(nq.values()))
        res = dict(persons=n, pairing=r, share_ratio=ratio, **{f"fraction_{k}": v for k, v in f.items()},
                   elongation_combined=target, phase_sd_rad=sd_hat, phase_sd_h=sd_hat * 24 / (2 * np.pi),
                   **{f"{k}_median": float(s[k].median()) for k in s}, **{f"{k}_q95": v for k, v in nq.items()},
                   matched_median=float(comp[comp.kind == "matched"].Q.median()), matched_q95=mq95,
                   responsive_median=float(comp[comp.kind == "responsive"].Q.median()),
                   responsive_q95=float(comp[comp.kind == "responsive"].Q.quantile(0.95)), threshold=thr,
                   power_gtex=float(np.mean(s.alt_gtex > thr)), power_half=float(np.mean(s.alt_half > thr)))
        res["gate"] = bool(res["power_gtex"] >= 0.8)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("A370 not powered: not decidable at this design; not scored")
    comp = pd.read_csv(OUT / "companions.csv")
    s = pd.read_csv(OUT / "sim.csv")
    q = q2(Z, idx)
    rng = np.random.default_rng([370, 2])
    za, zm = Z["adipose"][idx], Z["muscle"][idx]
    null = np.array([Q(za, zm[rng.permutation(n)]) for _ in range(N_PERM)])
    mt, rs = comp[comp.kind == "matched"].Q, comp[comp.kind == "responsive"].Q
    row = dict(test="H370.1", persons=n, Q=q, p_perm=float((1 + np.sum(null >= q)) / (1 + N_PERM)),
               threshold=float(pw["threshold"]), pct_timing_only=float((s["null_normal_1.0"] < q).mean() * 100),
               pct_alt_half=float((s.alt_half < q).mean() * 100), pct_alt_gtex=float((s.alt_gtex < q).mean() * 100),
               matched_q95=float(mt.quantile(.95)), pct_matched=float((mt < q).mean() * 100),
               responsive_q95=float(rs.quantile(.95)), pct_responsive=float((rs < q).mean() * 100),
               leverage_raw=E.radius_covariance(za, zm)["leverage"])
    row["passes"] = bool(q > row["threshold"] and row["p_perm"] < 0.05)
    core = c.wf.drop([g for g in TRIO if g in c.wf.genes])
    row["Q_core"] = q2(c.Z(core), idx)
    dr = pd.DataFrame([dict(drop="+".join(d), trio=sorted(d) == sorted(TRIO), Q=q2(c.Z(c.wf.drop(list(d))), idx))
                       for d in itertools.combinations(c.wf.genes, 3)])
    dr.to_csv(OUT / "drops.csv", index=False)
    row.update(drops_median=float(dr.Q.median()), trio_drop_pct=float((dr.Q < row["Q_core"]).mean() * 100))
    pd.DataFrame([row]).to_csv(OUT / "tests.csv", index=False)
    print(pd.Series(row).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
