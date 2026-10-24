"""A134 scout: feasibility counts and the null spread of the H134.1 statistic (no observed gain).

Reads output/scout_cells_<study>.csv (per patient and group: cells, mean copy-number profile value at the 14 clock-gene loci)
and the A87 pseudobulks of the same 79 patients. Writes
  scout_cn_feasibility.csv   patients with a copy-number state per locus (cells >= 20 / 50 / 100), per cohort and overall;
  scout_cn_distribution.csv  per locus: malignant and normal-epithelium patient-level score quantiles and the share of
                             patients above the normal-epithelium 97.5 % quantile (gain) and below its 2.5 % (loss);
  scout_null_spread.csv      bootstrap SD of the coherence gain (adjusted minus raw, malignant minus normal) when the
                             copy-number rows are permuted across patients, and MDE = 2.49 SD (1,000 bootstraps).
The observed gain (true copy-number labels) is NOT computed.

    python scout_null.py [--reps 1000]
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_5" / "scripts" / "cancer"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from clocksync.decompose import components, logfrac  # noqa: E402
from common import data  # noqa: E402

import _shared as sh  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
GENES, WF = sh.GENES, sh.WF
MAL, NORM = sh.MAL, sh.NORM
MDE_K = 2.49  # (z_0.95 + z_0.80) / one-sided 5 %, 80 % power ... as registered (2.49 x SD)


def load_cn() -> pd.DataFrame:
    parts = [pd.read_csv(f, dtype={"patient": str}) for f in sorted(OUT.glob("scout_cells_*.csv"))]
    c = pd.concat(parts, ignore_index=True)
    c["donor"] = c.study + "|" + c.patient
    return c


def feasibility(c: pd.DataFrame, pts: pd.DataFrame) -> pd.DataFrame:
    loci = pd.concat([pd.read_csv(f) for f in sorted(OUT.glob("scout_loci_*.csv"))])
    rows = []
    for grp in ("malignant", "normal"):
        g = c[c.group == grp].set_index("donor")
        for gene in GENES:
            v = g[f"cn_{gene}_mean"]
            for m in (20, 50, 100):
                ok = v[(g.cells >= m) & v.notna()]
                rows.append(dict(group=grp, gene=gene, min_cells=m, patients=int(len(ok)), of=len(pts),
                                 studies=int(ok.index.str.split("|").str[0].nunique())))
    f = pd.DataFrame(rows)
    # cohort-level locus coverage: nearest expressed gene within 1 Mb / own gene expressed
    cov = loci.groupby("gene").agg(cohorts=("study", "nunique"), own_expressed=("own", "sum"),
                                   max_dist_bp=("dist_bp", "max"), min_window=("win", "min"))
    f = f.merge(cov, left_on="gene", right_index=True)
    return f


def distribution(c: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in GENES:
        k = f"cn_{gene}_mean"
        m = c[(c.group == "malignant") & (c.cells >= 50)][k].dropna()
        n = c[(c.group == "normal") & (c.cells >= 20)][k].dropna()
        lo, hi = n.quantile(0.025), n.quantile(0.975)
        rows.append(dict(gene=gene, n_malignant=len(m), n_normal=len(n), mal_median=m.median(), mal_sd=m.std(),
                         mal_q05=m.quantile(.05), mal_q95=m.quantile(.95), norm_median=n.median(), norm_sd=n.std(),
                         norm_q025=lo, norm_q975=hi, gain_share=float((m > hi).mean()), loss_share=float((m < lo).mean())))
    return pd.DataFrame(rows)


def centre(x: np.ndarray, g: np.ndarray) -> np.ndarray:
    return x - pd.DataFrame(x).groupby(g).transform("mean").to_numpy()


def adjusted_coherence(y: np.ndarray, pv, k, cn: np.ndarray, g: np.ndarray) -> float:
    """Coherence of the residuals of each gene's centred log fraction on its own (stratum-centred) copy-number score."""
    z = centre(cn, g)
    den = (z ** 2).sum(0)
    beta = np.where(den > 0, (y * z).sum(0) / np.where(den > 0, den, 1), 0.0)
    return components(y - z * beta, pv, k, WF)["coherence"]


def null_spread(tab, cn_m, cn_n, reps, seed, within):
    """Bootstrap (patients) of gain_malignant - gain_normal with the copy-number rows permuted across the resampled patients."""
    rng = np.random.default_rng(seed)
    case, ref = tab["case"], tab["ref"]
    n = len(case)
    strata = case.stratum.to_numpy()
    out = []
    for _ in range(reps):
        idx = rng.integers(0, n, n)
        ys = []
        for d in (case, ref):
            x = d.iloc[idx]
            ys.append(logfrac(x, GENES, x.stratum.to_numpy()))
        g = strata[idx]
        if within:  # permute inside strata
            pm = np.arange(n)
            for s in np.unique(g):
                j = np.flatnonzero(g == s)
                pm[j] = rng.permutation(j)
            pn = pm.copy()
        else:
            pm, pn = rng.permutation(n), rng.permutation(n)
        res = {}
        for name, (y, pv, k), cn, p in (("mal", ys[0], cn_m[idx], pm), ("norm", ys[1], cn_n[idx], pn)):
            raw = components(y, pv, k, WF)["coherence"]
            adj = adjusted_coherence(y, pv, k, cn[p], g)
            res[name] = (raw, adj - raw)
        # alternative: the patient's malignant copy number used for the normal epithelium
        y, pv, k = ys[1]
        adj_b = adjusted_coherence(y, pv, k, cn_m[idx][pm], g)
        out.append(dict(gain_mal=res["mal"][1], gain_norm=res["norm"][1], gain_norm_malcn=adj_b - res["norm"][0],
                        raw_mal=res["mal"][0], raw_norm=res["norm"][0]))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=134)
    a = ap.parse_args()
    limit_threads(2)
    warnings.filterwarnings("ignore")
    pts = pd.read_csv(OUT / "scout_patients79.csv")
    c = load_cn()
    c = c[c.donor.isin(pts.donor)]
    feasibility(c, pts).to_csv(OUT / "scout_cn_feasibility.csv", index=False)
    distribution(c).to_csv(OUT / "scout_cn_distribution.csv", index=False)

    pb, _ = data.load("cancer_sc")
    pb = pb[pb.donor.isin(pts.donor)]
    mal = pb[pb.compartment == MAL].drop_duplicates("donor").set_index("donor")
    nor = pb[pb.compartment == NORM].drop_duplicates("donor").set_index("donor")
    # patients with a copy-number state at all 14 loci in both groups (>= 50 malignant cells)
    cm = c[c.group == "malignant"].set_index("donor")
    cn_ = c[c.group == "normal"].set_index("donor")
    cols = [f"cn_{g}_mean" for g in GENES]
    ok = [d for d in mal.index if d in nor.index and d in cm.index and d in cn_.index and cm.loc[d, "cells"] >= 50
          and cm.loc[d, cols].notna().all() and cn_.loc[d, cols].notna().all()]
    print("patients with all 14 loci, malignant and normal:", len(ok), "of", len(mal.index.intersection(nor.index)))
    case, ref = mal.loc[ok], nor.loc[ok]
    case = case.assign(unit=case.index)
    tab = {"case": case, "ref": ref}
    cn_m, cn_n = cm.loc[ok, cols].to_numpy(float), cn_.loc[ok, cols].to_numpy(float)
    # reproduction check of the raw coherence (A87: 0.228 malignant, 0.782 normal on the 79 patients)
    chk = {}
    for nm, d in (("malignant", case), ("normal", ref)):
        y, pv, k = logfrac(d, GENES, d.stratum.to_numpy())
        chk[nm] = components(y, pv, k, WF)["coherence"]
    print("raw coherence on these patients:", {k: round(v, 3) for k, v in chk.items()})
    rows = []
    for within in (False, True):
        r = null_spread(tab, cn_m, cn_n, a.reps, a.seed, within)
        d = r.gain_mal - r.gain_norm
        d2 = r.gain_mal - r.gain_norm_malcn
        for nm, v in (("mal_minus_norm_ownCN", d), ("mal_minus_norm_malCN", d2), ("mal_gain", r.gain_mal),
                      ("norm_gain_ownCN", r.gain_norm)):
            rows.append(dict(permutation="within_stratum" if within else "across_all", statistic=nm,
                             null_mean=float(v.mean()), null_sd=float(v.std(ddof=1)), mde=MDE_K * float(v.std(ddof=1)),
                             reps=a.reps, patients=len(ok), raw_coh_mal=chk["malignant"], raw_coh_norm=chk["normal"]))
    pd.DataFrame(rows).to_csv(OUT / "scout_null_spread.csv", index=False)
    print(pd.DataFrame(rows).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
