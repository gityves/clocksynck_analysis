"""Scout (not registered): timed clock amplitude of T cells vs fibroblasts in the mouse colon chrono-atlas.
Per arm/compartment: debiased cosinor amplitude over mice (log-CPM of 14 clock-gene pseudobulks), bootstrap of A_T/A_fib.
    python amplitude_ratio_colon.py
"""
import sys, warnings
from pathlib import Path
PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np, pandas as pd
import clocksync as cs
from clocksync.config import limit_threads
from common import data
limit_threads(4)
warnings.filterwarnings("ignore")
OUT = Path(__file__).parent / "output"
ARMS = ("C-DC", "R-DC", "R-PC")
CTX = ["T", "fibroblast", "B", "plasma", "macrophage", "endothelial", "smooth_muscle"]
ALL = ["ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "RORA", "TEF", "CIART"]
CORE = [g for g in ALL if g != "RORA"]
SETS = {"all14": ALL, "core13": CORE}
MIN_CELLS, MIN_MICE, NB = 30, 6, 2000


def logcpm(t, genes):
    return np.log2(1 + t[genes].to_numpy(float) / t.depth.to_numpy(float)[:, None] * 1e6)


def amp2(Y, zt):
    """debiased squared amplitude (sum over genes); Y mice x genes."""
    n = len(zt)
    X = np.column_stack([np.ones(n), np.cos(2 * np.pi * zt / 24), np.sin(2 * np.pi * zt / 24)])
    B, *_ = np.linalg.lstsq(X, Y, rcond=None)
    s2 = ((Y - X @ B) ** 2).sum(0) / max(n - 3, 1)
    raw = B[1] ** 2 + B[2] ** 2
    cov = np.linalg.inv(X.T @ X)  # var(coef) = s2 * cov
    bias = s2 * (cov[1, 1] + cov[2, 2])  # E[a^2+b^2] under no rhythm
    return max(float((raw - bias).sum()), 0.0)


def main():
    pb, donors = data.load("colon_chrono")
    donors = donors.set_axis(donors.index.astype(str))
    pb = pb.assign(arm=pb.donor.map(donors.arm), zt=pb.donor.map(donors.zt).astype(float))
    rng = np.random.default_rng(1)
    rows, summ = [], []
    for arm in ARMS:
        t = pb[(pb.arm == arm) & (pb.n_cells >= MIN_CELLS)]
        # all-cell counts per mouse for depth reporting (before filter)
        tall = pb[pb.arm == arm]
        # timed clocksync (context)
        try:
            ds = cs.Dataset(t.drop(columns="arm").drop(columns="zt"), anchors=None, donors=donors, min_donors=MIN_MICE)
            T = ds.timed_loadings("zt")
        except Exception as e:
            T = pd.DataFrame(); print("clocksync timed failed", arm, e)
        for cname, genes in SETS.items():
            boots, est = {}, {}
            for c in CTX:
                g = t[t.compartment == c].sort_values("donor")
                if g.donor.nunique() < MIN_MICE:
                    continue
                Y, zt = logcpm(g, genes), g.zt.to_numpy()
                est[c] = amp2(Y, zt)
                idx_by_zt = [np.flatnonzero(zt == z) for z in np.unique(zt)]
                b = np.empty(NB)
                for k in range(NB):
                    ii = np.concatenate([rng.choice(ix, len(ix)) for ix in idx_by_zt])
                    b[k] = amp2(Y[ii], zt[ii])
                boots[c] = b
                rows.append(dict(arm=arm, compartment=c, gene_set=cname, n_genes=len(genes), A=np.sqrt(est[c]),
                                 A_boot_median=np.sqrt(np.median(b)), A_boot_lo=np.sqrt(np.quantile(b, .025)), A_boot_hi=np.sqrt(np.quantile(b, .975)),
                                 n_mice=g.donor.nunique(), cells_median=g.n_cells.median(), cells_min=g.n_cells.min(),
                                 umis_median=g[genes].sum(1).median(), depth_median=g.depth.median(),
                                 clocksync_abs_amplitude=T.abs_amplitude.get(c, np.nan) if len(T) else np.nan,
                                 clocksync_coupling=T.coupling.get(c, np.nan) if len(T) else np.nan))
            if "fibroblast" in boots:
                for c in boots:
                    if c == "fibroblast":
                        continue
                    r = np.sqrt(boots[c]) / np.sqrt(np.maximum(boots["fibroblast"], 1e-12))
                    summ.append(dict(arm=arm, compartment=c, gene_set=cname, A_c=np.sqrt(est[c]), A_fib=np.sqrt(est["fibroblast"]),
                                     ratio_point=np.sqrt(est[c]) / np.sqrt(est["fibroblast"]),
                                     ratio_boot_median=np.median(r), ratio_q025=np.quantile(r, .025), ratio_q975=np.quantile(r, .975),
                                     ratio_q95=np.quantile(r, .95)))
    A = pd.DataFrame(rows); R = pd.DataFrame(summ)
    A.to_csv(OUT / "amplitudes.csv", index=False)
    R.to_csv(OUT / "ratios_vs_fibroblast.csv", index=False)
    M = A.merge(R.drop(columns=["A_c", "A_fib"]), on=["arm", "compartment", "gene_set"], how="left")
    M.to_csv(OUT / "amplitude_ratios.csv", index=False)
    pd.set_option("display.width", 250)
    print(A.round(3).to_string()); print(R.round(3).to_string())
    # per-mouse T cell counts
    tc = pb[pb.compartment.isin(["T", "fibroblast"])].pivot_table(index="donor", columns="compartment", values="n_cells")
    print(tc.to_string())

main()
