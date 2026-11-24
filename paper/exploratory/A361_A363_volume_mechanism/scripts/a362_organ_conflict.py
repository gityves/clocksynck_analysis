"""A362: do organs that lose the most amplitude under conflicting drives in mouse follow the person's volume more steeply?

Mouse conflict sensitivity k_o = -log(amplitude kept under feeding inversion), A358 organ table (liver, adipose, lung,
heart, kidney; mean of A322 Manella and A324 Xin where both exist). Human: A312 elasticity (`clocksync.elasticity`, 20
halvings of the other tissues, seed [312, tissue index]) of the mapped GTEx tissues on the core loop (blind) and on the 14
genes (A312, not blind), organ value = mean over its tissues (adipose: subcutaneous and visceral; heart: atrial
appendage and left ventricle). Each tissue's elasticity is corrected by its mean under the calibration null (every
tissue in proportion, simulated at the GTEx design); statistic c = least-squares slope of the corrected organ elasticity
on k_o; null and power by simulation.

    python a362_organ_conflict.py --stage power [--jobs 4]   # simulations only (no GTEx expression)
    python a362_organ_conflict.py --stage score               # GTEx elasticities and the registered test
Outputs (`../output/`): a362_sims.csv, a362_power.csv | a362_organs.csv, a362_tests.csv. Aggregates only.
"""
import argparse
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import vm_common as V

sys.path.insert(0, str(V.PAPER / "results" / "figure_2" / "scripts" / "backbone"))
import clocksync as cs  # noqa: E402
from clocksync import elasticity as E  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

MAP = {"liver": ["Liver"], "adipose": ["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)"], "lung": ["Lung"],
       "heart": ["Heart - Atrial Appendage", "Heart - Left Ventricle"], "kidney": ["Kidney - Cortex"]}
TARGETS = [t for v in MAP.values() for t in v]
DESCRIPTIVE = ["Muscle - Skeletal"]
N_NULL, N_ALT, N_SPLITS = 40, 20, 20
A312 = V.EXP / "A312_elasticity" / "output" / "tissues.csv"


def organ_k() -> pd.Series:
    o = pd.read_csv(V.EXP / "A358_A360_gain_sources/output/a358_organs.csv", index_col=0)
    return -np.log(o.kept.reindex(list(MAP)))


def slope(e_t: pd.Series, k: pd.Series):
    """Organ means of the (corrected) tissue elasticities; OLS intercept a and slope c on k; conflict share at mean k."""
    e = pd.Series({o: e_t[ts].mean() for o, ts in MAP.items()})
    c, a = np.polyfit(k[e.index].to_numpy(), e.to_numpy(), 1)
    kb = float(k.mean())
    return e, float(a), float(c), float(c * kb / (a + c * kb)) if a + c * kb != 0 else np.nan


def expo_for(kind, k):
    kb = float(k.mean())
    if kind == "null":
        return {}
    share = 1.0 if kind == "pure" else 0.5
    return {t: (1 - share) + share * k[o] / kb for o, ts in MAP.items() for t in ts}


def sim_one(args):
    kind, draw = args
    V.quiet()
    limit_threads(1)
    k = organ_k()
    t = pd.read_csv(A312, index_col="tissue")
    table = V.tissue_table(t.index)
    presence = (t.n / 244).clip(upper=1.0)
    rng = np.random.default_rng([362, {"null": 0, "pure": 1, "half": 2}[kind], draw])
    pb = V.simulate(table, presence, V.CORE, rng, expo=expo_for(kind, k))
    comps = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb, anchors=comps, waveform=V.CORE, name=f"a362{kind}{draw}")
    out = dict(kind=kind, draw=draw)
    for t_ in TARGETS:
        out[t_] = E.amplitude_elasticity(ds, t_, n_splits=N_SPLITS, seed=[362, draw, comps.index(t_)])["elasticity"]
    return out


def main():
    V.quiet()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    V.OUT.mkdir(parents=True, exist_ok=True)
    k = organ_k()
    if a.stage == "power":
        jobs = [("null", d) for d in range(N_NULL)] + [(kd, d) for kd in ("pure", "half") for d in range(N_ALT)]
        rows = []
        with Pool(a.jobs) as pool:
            for r in pool.imap_unordered(sim_one, jobs):
                rows.append(r)
                print(f"{len(rows)}/{len(jobs)} simulated cohorts ({time.time() - t0:.0f}s)", flush=True)
        s = pd.DataFrame(rows).sort_values(["kind", "draw"])
        s.to_csv(V.OUT / "a362_sims.csv", index=False)
        bias = s[s.kind == "null"][TARGETS].mean()
        res = {}
        for kind in ("null", "pure", "half"):
            cs_ = np.array([slope(r[TARGETS] - bias + 1, k)[2] for _, r in s[s.kind == kind].iterrows()])
            res[f"c_{kind}_median"] = float(np.median(cs_))
            if kind == "null":
                q95 = float(np.quantile(cs_, 0.95))
                res["c_null_q95"], res["c_null_sd"] = q95, float(cs_.std(ddof=1))
            else:
                res[f"power_{kind}"] = float(np.mean(cs_ > q95))
        res.update({f"bias_{t_}": float(bias[t_]) for t_ in TARGETS})
        res["gate"] = bool(res["power_pure"] >= 0.8)
        pd.Series(res).to_csv(V.OUT / "a362_power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    import compute_oscillator_gtex as A104
    from common import gtex
    pw = pd.read_csv(V.OUT / "a362_power.csv", index_col=0)["value"]
    sims = pd.read_csv(V.OUT / "a362_sims.csv", keep_default_na=False)  # "null" is a label, not NaN
    null = sims[sims.kind == "null"]
    bias = null[TARGETS].mean()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS)
    pb = pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)
    comps = sorted(pb.compartment.unique())
    ref = pd.read_csv(A312, index_col="tissue").elasticity
    rows, est = [], {}
    for gs, w in (("14 genes", V.WF), ("core", V.CORE)):
        ds = cs.Dataset(pb, anchors=comps, waveform=w, name=f"a362_{gs}")
        e = {t_: E.amplitude_elasticity(ds, t_, n_splits=N_SPLITS, seed=[312, comps.index(t_)])["elasticity"]
             for t_ in TARGETS + DESCRIPTIVE}
        if gs == "14 genes" and not np.allclose([e[t_] for t_ in TARGETS], ref[TARGETS].to_numpy(), atol=1e-9):
            raise SystemExit("A312 elasticities not reproduced on the 14 genes: stop before the core loop")
        est[gs] = pd.Series(e)
    for gs in ("core", "14 genes"):
        corr = est[gs][TARGETS] - bias + 1          # null bias from the core-loop simulation (used for both)
        eo, a_, c, share = slope(corr, k)
        cnull = np.array([slope(r[TARGETS] - bias + 1, k)[2] for _, r in null.iterrows()])
        p = float((1 + np.sum(cnull >= c)) / (1 + len(cnull)))
        rows.append(dict(test="H362.1" if gs == "core" else "H362.1 (14 genes, reported, not blind)", intercept=a_,
                         slope=c, conflict_share=share, p_sim=p, c_null_q95=float(np.quantile(cnull, 0.95)),
                         spearman=float(spearmanr(k[eo.index], eo).statistic),
                         passes=bool(c > np.quantile(cnull, 0.95)) if gs == "core" else np.nan,
                         gate=str(pw["gate"])))
        for o in eo.index:
            rows.append(dict(test=f"organ {gs}", organ=o, k=float(k[o]), elasticity_corrected=float(eo[o])))
        for t_ in TARGETS + DESCRIPTIVE:
            rows.append(dict(test=f"tissue {gs}", organ=t_, elasticity_raw=float(est[gs][t_]),
                             null_mean=float(bias.get(t_, np.nan))))
    out = pd.DataFrame(rows)
    out.to_csv(V.OUT / "a362_tests.csv", index=False)
    print(out.round(3).to_string(index=False), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
