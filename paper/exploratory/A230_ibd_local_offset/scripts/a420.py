"""A420 (POST HOC, DESCRIPTIVE; registered in paper/docs/registered_experiments.md, "A420"): within an endoscopy (biopsies
sharing the hour), does a biopsy's endothelial (and fibroblast) clock loudness fall with its own IFN-gamma response?

TAURUS genome pseudobulks (biopsy x class), A230 sample metadata. Loudness = log radius of the class's 14-gene state
(`clocksync.pairwise.states`, one projector per class), states centred within library batch. x = the class's IFN-gamma
response score (STAT1, IRF1, CXCL9, CXCL10; mean within-batch z of log2 CPM); inflammation score descriptive.
Within-endoscopy slope (both demeaned within patient x timepoint, pooled, standardised); permutations of x within
endoscopy; companions: 200 matched and 200 responsive pseudo-waveforms on the same rows.
                                                     -> output/a420/{tests,companions}.csv
    python -I a420.py                                                    (seconds to minutes)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import clocksync as cs
from clocksync import pairwise
from clocksync.config import limit_threads
from clocksync.genes import canonical
from clocksync.paired import pseudo_waveforms

HERE = Path(__file__).resolve()
DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/ibd/TAURUS/pseudobulk_genome")
META = HERE.parents[1] / "output" / "data" / "taurus" / "donors.csv"
OUT = HERE.parents[1] / "output" / "a420"
CLASSES = {"endothelial": "gut:endothelial", "fibroblast": "gut:fibroblast"}
IFN = [canonical(g) for g in ("STAT1", "IRF1", "CXCL9", "CXCL10")]
N_SET, N_PERM, SEED = 200, 10_000, 420
BASE = cs.gtex()


def load():
    z = np.load(DATA / "counts.npz", allow_pickle=False)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    genes = [canonical(g) for g in pd.read_csv(DATA / "genes.csv").gene.astype(str)]
    rows = pd.read_csv(DATA / "rows.csv")
    keep = ~pd.Index(genes).duplicated()
    X = X[:, np.flatnonzero(keep)]
    genes = list(pd.Index(genes)[keep])
    return X, genes, rows


def table(X, genes, rows, comp, meta, gene_list):
    i = np.flatnonzero(rows.compartment.to_numpy() == comp)
    cols = [genes.index(g) for g in gene_list if g in genes]
    t = pd.DataFrame(X[i][:, cols].toarray(), columns=[genes[c] for c in cols])
    t.insert(0, "depth", rows.depth.to_numpy()[i])
    t.insert(0, "n_cells", rows.n_cells.to_numpy()[i])
    t.insert(0, "compartment", comp)
    t.insert(0, "donor", rows.donor.to_numpy()[i])
    t = t.merge(meta[["donor", "patient", "timepoint", "Batch", "inflammation_score"]], on="donor", how="inner")
    t["stratum"] = t.Batch.astype(str)
    t["endoscopy"] = t.patient.astype(str) + "|" + t.timepoint.astype(str)
    return t


def centred(v, g):
    s = pd.Series(v)
    return (s - s.groupby(np.asarray(g)).transform("mean")).to_numpy()


def loudness(t, wf):
    z = centred(pairwise.states(t, wf), t.stratum)
    return np.log(np.abs(z))


def ifn_score(t):
    g = [x for x in IFN if x in t.columns]
    lc = np.log2(t[g].div(t.depth, axis=0) * 1e6 + 1)
    gr = t.stratum.to_numpy()
    z = (lc - lc.groupby(gr).transform("mean")) / lc.groupby(gr).transform(lambda x: x.std(ddof=0)).replace(0, np.nan)
    return z.mean(axis=1).to_numpy()


def within(y, x, grp):
    ok = np.isfinite(y) & np.isfinite(x)
    y, x, grp = y[ok], x[ok], grp[ok]
    n = pd.Series(grp).map(pd.Series(grp).value_counts()).to_numpy()
    m = n >= 2
    y, x, grp = y[m], x[m], grp[m]
    yr, xr = centred(y, grp), centred(x, grp)
    return float(np.corrcoef(yr, xr)[0, 1]), yr, xr, grp


def main():
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    X, genes, rows = load()
    meta = pd.read_csv(META, keep_default_na=False)
    meta["inflammation_score"] = pd.to_numeric(meta.inflammation_score, errors="coerce")
    st_rows = rows[rows.compartment.isin(CLASSES.values())]
    tot = np.asarray(X[st_rows.index.to_numpy()].sum(0)).ravel()
    lm = pd.Series(np.log(tot / tot.sum() + 1e-12), index=genes)
    mean_counts = pd.Series(tot / len(st_rows), index=genes)
    lm = lm[(mean_counts >= 1) | lm.index.isin(list(BASE.genes))]
    wfs = {fam: pseudo_waveforms(lm, N_SET, np.random.default_rng([SEED, code]), BASE, fam) for code, fam in enumerate(("matched", "responsive"))}
    need = sorted(set(BASE.genes) | set(IFN) | {g for ws in wfs.values() for w in ws for g in w.genes})
    rows_out, comp_rows = [], []
    for k, comp in CLASSES.items():
        t = table(X, genes, rows, comp, meta, need)
        x = ifn_score(t)
        y = loudness(t, BASE)
        grp = t.endoscopy.to_numpy()
        r, yr, xr, g = within(y, x, grp)
        rng = np.random.default_rng([SEED, 9, len(k)])
        nul = []
        for _ in range(N_PERM):
            xp = xr.copy()
            for u in np.unique(g):
                i = np.flatnonzero(g == u)
                xp[i] = xp[rng.permutation(i)]
            nul.append(np.corrcoef(yr, xp)[0, 1])
        nul = np.array(nul)
        ri, *_ = within(y, t.inflammation_score.to_numpy(float), grp)
        rx, *_ = within(x, t.inflammation_score.to_numpy(float), grp)
        row = dict(cls=k, biopsies=int(len(yr)), endoscopies=int(len(np.unique(g))), r_loud_ifn=r,
                   p_two_sided=float((np.sum(np.abs(nul) >= abs(r)) + 1) / (N_PERM + 1)), p_lower=float((np.sum(nul <= r) + 1) / (N_PERM + 1)),
                   r_loud_inflammation=ri, r_ifn_inflammation=rx)
        for fam, ws in wfs.items():
            v = []
            for j, w in enumerate(ws):
                try:
                    v.append(within(loudness(t, w), x, grp)[0])
                except Exception:  # noqa: BLE001
                    v.append(np.nan)
                comp_rows.append(dict(cls=k, family=fam, set=j, r=v[-1]))
            v = np.array(v)[np.isfinite(v)]
            row.update({f"{fam}_median": float(np.median(v)), f"{fam}_q05": float(np.quantile(v, .05)), f"{fam}_pct": float(100 * np.mean(v <= r))})
        rows_out.append(row)
    T = pd.DataFrame(rows_out)
    T.to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(comp_rows).to_csv(OUT / "companions.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(T.round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
