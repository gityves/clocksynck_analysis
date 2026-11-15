"""Scout for the clock-autonomous core of muscle (follow-up of A261): can a program of the core's size clear the A261
rhythm gate and resolve central / tissue clock / behaviour on GSE106800, judged on matched sets only?

The core set itself is never read on the arrays. Power is measured on stand-in programs: reference-rhythmic genes
(GSE108539, p < 0.01, on the array, clock genes excluded) outside the core, drawn to match the core's reference
amplitudes by decile. Per stand-in:
  - gate: aligned loading |mean z_al| above the q95 of 200 non-rhythmic sets matched on mean array intensity (A261 rule);
  - resolution under each truth (central, tissue clock, behaviour): the misaligned loadings are replaced by
    z_mis,i = r e^{i w s} mean(z_al) + (z_mis,i - mean z_mis), i.e. the hypothesised mean plus the stand-in's own
    misaligned residuals, and classified by the A261 Hotelling tests at level 0.01.
Without a core list (--sizes), stand-ins are drawn without amplitude matching at each size.

    python scout_core_power.py [--core genes.txt] [--sizes 20 50 100 200] [--draws 300]
Outputs (output/): power.csv (one row per size or core: gate rate, resolution per truth), power_draws.csv.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/scouts_2026-11-15/T1_nightshift"))
sys.path.insert(0, str(PAPER / "exploratory/A261_nightshift_programs/scripts"))
sys.path.insert(0, str(PAPER / "exploratory/A257_nightshift_muscle/scripts"))
from score_a261 import CLOCK, HYP, W, per_subject_z, read  # noqa: E402
from scout_t1 import arrays, reference  # noqa: E402

OUT = HERE / "output"
SEED, N_NULL = 1116, 200


def setup():
    ref = reference().set_index("gene")
    d = arrays()
    genes = [g for g in ref.index if g in d["ctl"].columns and g not in CLOCK]
    ref = ref.loc[genes]
    from compute_a257 import load
    _, y = load()
    mi = (y / np.log(2.0)).mean(0)
    mi = mi[~mi.index.duplicated()].reindex(genes)
    C = 2 * ref.amp * np.cos(W * (19 - ref.acro_h))
    S = 2 * ref.amp * np.sin(W * (19 - ref.acro_h))
    return ref, d, pd.qcut(mi, 10, labels=False), C, S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--core", type=Path, help="one human symbol per line (the core set, fixed beforehand)")
    ap.add_argument("--label", default=None)
    ap.add_argument("--sizes", type=int, nargs="*", default=[20, 50, 100, 200])
    ap.add_argument("--draws", type=int, default=300)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    ref, d, dec, C, S = setup()
    al, mis = d["ctl"], d["mis"]
    rhy = ref.index[ref.p < 0.01]
    non = ref.index[ref.p > 0.5]
    by_dec = {k: np.array(non[(dec.reindex(non) == k).to_numpy()]) for k in range(10)}

    def loads(gs):
        gs = list(gs)
        return (per_subject_z(al[gs].to_numpy(float), C[gs].to_numpy(), S[gs].to_numpy()),
                per_subject_z(mis[gs].to_numpy(float), C[gs].to_numpy(), S[gs].to_numpy()))

    def gate_q95(decs):
        v = []
        for _ in range(N_NULL):
            pick = [rng.choice(by_dec[int(k)]) for k in decs]
            v.append(abs(loads(pick)[0].mean()))
        return float(np.quantile(v, 0.95))

    if a.core:
        core = {ln.strip().upper() for ln in open(a.core) if ln.strip()}
        core_r = [g for g in rhy if g.upper() in core]
        pool = rhy[~rhy.isin(core_r)]
        amp_dec = pd.qcut(ref.amp.loc[rhy], 10, labels=False)
        jobs = [(a.label or a.core.stem, len(core_r), core_r)]
        print(f"core {a.core.name}: {len(core)} listed, {len(core_r)} reference-rhythmic on the array "
              f"(clock genes excluded); complement {len(pool)}", flush=True)
    else:
        pool, amp_dec = rhy, None
        jobs = [(f"size_{G}", G, None) for G in a.sizes]

    rows, drows = [], []
    for label, G, core_r in jobs:
        if G < 2:
            rows.append(dict(set=label, G=G))
            continue
        # gate null matched on intensity: the core's deciles if given (a level, not a time statistic), else random
        decs = dec.loc[core_r].to_numpy() if core_r else rng.choice(dec.dropna().to_numpy(), G)
        q95 = gate_q95(decs)
        for i in range(a.draws):
            if core_r:
                pick = [rng.choice(pool[(amp_dec.loc[pool] == amp_dec.loc[g]).to_numpy()]) for g in core_r]
                pick = list(dict.fromkeys(pick))
            else:
                pick = list(rng.choice(pool, G, replace=False))
            zA, zM = loads(pick)
            row = dict(set=label, draw=i, G=len(pick), z_al_abs=float(abs(zA.mean())), gated=bool(abs(zA.mean()) > q95),
                       empirical=read(zA, zM)["cls"])
            for h, (s, r) in HYP.items():
                zMs = r * np.exp(1j * W * s) * zA.mean() + (zM - zM.mean())
                row[f"truth_{h}"] = read(zA, zMs)["cls"]
            drows.append(row)
        D = pd.DataFrame([r for r in drows if r["set"] == label])
        g = D[D.gated]
        out = dict(set=label, G=G, gate_q95=q95, gate_rate=float(D.gated.mean()),
                   median_aligned_over_q95=float(D.z_al_abs.median() / q95))
        for h in HYP:
            out[f"resolved_{h}"] = float((g[f"truth_{h}"] == h).mean()) if len(g) else np.nan
            out[f"power_{h}"] = float(((D[f"truth_{h}"] == h) & D.gated).mean())
        out["empirical_behaviour"] = float((g.empirical == "behaviour").mean()) if len(g) else np.nan
        rows.append(out)
        print(pd.Series(out).to_string(), flush=True)
    tag = f"_{a.label or a.core.stem}" if a.core else ""
    pd.DataFrame(rows).to_csv(OUT / f"power{tag}.csv", index=False)
    pd.DataFrame(drows).to_csv(OUT / f"power_draws{tag}.csv", index=False)


if __name__ == "__main__":
    main()
