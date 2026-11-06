"""A223: is the tumour-vs-paired-normal clock-gene level pattern a damped clock read at daytime surgery? (registered in
paper/docs/registered_experiments.md, A223).

Per type: mean paired log2 LFC of the 14 clock genes regressed (OLS, intercept) on x_g = -amp_g C_g / ln 2, C_g = mean cos(w t -
psi_g) over t in [8, 18) h (GTEx waveform); B = mean slope over types. Null: all psi rotated by delta (23 whole-hour offsets + 977
uniform draws). CPTAC-3 (A186 cache) primary, without PER1 secondary, TCGA (A221 cache) secondary.

    python compute_a223.py --stage power|tests
Outputs: ../output/{lfc.csv (tests only), power.csv, tests.csv, per_type.csv, window_scan.csv}
"""
from __future__ import annotations

import argparse
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync.genes import gene_columns

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "output"
EXP = HERE.parents[1]
CPTAC = EXP / "A186_cptac_c1" / "output" / "cache"
TCGA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/gtex_recount3/cache_a221")
WF = cs.gtex()
GENES = list(WF.genes)
W = 2 * np.pi / 24
WINDOW = (8.0, 18.0)
MIN_PAIRS, NROT, K = 20, 1000, 2.49
SEED = zlib.crc32(b"A223") % (2 ** 31)


def C_of(psi, window=WINDOW):
    t = np.linspace(*window, 241)[:-1]
    return np.cos(W * t[:, None] - psi[None, :]).mean(0)


def x_of(genes, delta=0.0, window=WINDOW):
    i = [GENES.index(g) for g in genes]
    return -WF.amp[i] * C_of(WF.psi[i] + W * delta, window) / np.log(2)


def lfc_cptac():
    out = {}
    for f in sorted(CPTAC.glob("*.npz")):
        z = np.load(f, allow_pickle=False)
        if len(z["cases"]) < MIN_PAIRS:
            continue
        g = list(z["genes"])
        ct, cn = z["ct"].astype(float), z["cn"].astype(float)
        cols = gene_columns(g, GENES, missing="raise")
        lt = np.log2(ct[:, cols] / ct.sum(1, keepdims=True) * 1e6 + 1)
        ln = np.log2(cn[:, cols] / cn.sum(1, keepdims=True) * 1e6 + 1)
        out[f.stem] = (lt - ln).mean(0)
    return pd.DataFrame(out, index=GENES).T


def lfc_tcga():
    out = {}
    for f in sorted(TCGA.glob("*.npz")):
        z = np.load(f, allow_pickle=False)
        cols = gene_columns(list(z["genes"]), GENES, missing="raise")
        out[f.stem] = (z["tum"][z["pair_t"]][:, cols] - z["nor"][z["pair_n"]][:, cols]).mean(0)
    return pd.DataFrame(out, index=GENES).T


def slopes(L: pd.DataFrame, genes, delta=0.0, window=WINDOW):
    x = x_of(genes, delta, window)
    X = np.column_stack([np.ones_like(x), x])
    Y = L[genes].to_numpy().T
    b = np.linalg.lstsq(X, Y, rcond=None)[0][1]
    return pd.Series(b, index=L.index)


def deltas(rng):
    return np.r_[np.arange(1, 24), rng.uniform(0, 24, NROT - 23)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", required=True, choices=["power", "tests"])
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    D = deltas(rng)
    sets = {"H223.1": ("CPTAC", lfc_cptac(), GENES), "H223.2": ("CPTAC", None, [g for g in GENES if g != "PER1"]),
            "H223.3": ("TCGA", lfc_tcga(), GENES)}
    sets["H223.2"] = ("CPTAC", sets["H223.1"][1], sets["H223.2"][2])
    pw, res, per = [], [], []
    for h, (cohort, L, genes) in sets.items():
        null = np.array([slopes(L, genes, d).mean() for d in D])
        pw.append(dict(test=h, cohort=cohort, types=len(L), genes=len(genes), null_mean=null.mean(), null_sd=null.std(ddof=1),
                       mde=K * null.std(ddof=1), scored=K * null.std(ddof=1) <= 0.5))
        if a.stage == "tests":
            b = slopes(L, genes)
            B = b.mean()
            res.append(dict(test=h, cohort=cohort, B=B, p=(1 + (null >= B).sum()) / (1 + len(null)), pctile=100 * (null < B).mean()))
            x = x_of(genes)
            for tp, bb in b.items():
                y = L.loc[tp, genes].to_numpy()
                r2 = np.corrcoef(x, y)[0, 1] ** 2
                per.append(dict(test=h, type=tp, beta=bb, r2=r2))
    pw = pd.DataFrame(pw)
    pw.to_csv(OUT / "power.csv", index=False)
    print(pw.round(4).to_string(index=False))
    if a.stage == "power":
        return
    res = pd.DataFrame(res)
    order = np.argsort(res.p.to_numpy())
    res.loc[res.index[order], "holm"] = np.maximum.accumulate([min(1, (3 - k) * res.p.iloc[i]) for k, i in enumerate(order)])
    res = res.merge(pw[["test", "null_mean", "null_sd", "mde"]], on="test")
    res.to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(per).to_csv(OUT / "per_type.csv", index=False)
    pd.concat({c: L for c, L in (("CPTAC", sets["H223.1"][1]), ("TCGA", sets["H223.3"][1]))}).to_csv(OUT / "lfc.csv")
    scan = [dict(cohort=c, centre=h, B=slopes(L, GENES, window=(h - 5, h + 5)).mean())
            for c, L in (("CPTAC", sets["H223.1"][1]), ("TCGA", sets["H223.3"][1])) for h in range(24)]
    pd.DataFrame(scan).to_csv(OUT / "window_scan.csv", index=False)
    print(res.round(4).to_string(index=False))
    print(pd.DataFrame(per).round(3).to_string(index=False))
    s = pd.DataFrame(scan)
    print(s.loc[s.groupby("cohort").B.idxmax()].to_string(index=False))


if __name__ == "__main__":
    main()
