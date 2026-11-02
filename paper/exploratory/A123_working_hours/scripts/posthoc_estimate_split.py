"""A123c post hoc (2026-11-02, descriptive, not registered): reproducibility of the clocksync working-hours organ map itself.

The 111 window donors (fast deaths 08-18 h) are split into random halves balanced within tissue (`split_donors`, as the truth
split-half in diag_truth_iii.py); clocksync is refitted on each half with the A123 cross-fitted waveform (out-of-window donors),
and the two half maps are compared (offset circular r, as `compare`). Each half map is also scored against truth (i) and (iii).
Set next to the truth (iii) split-half r (median 0.544, `a123c/diag/split_half.csv`): the estimate's reliability needs no truth.
Only summary statistics are written; donor identities and splits stay in memory.

    python posthoc_estimate_split.py [--reps 40] [--jobs 8]
"""
import argparse
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import importlib.util  # noqa: E402

_s = importlib.util.spec_from_file_location("cm", HERE.parent / "compute_methods.py")
cm = importlib.util.module_from_spec(_s)
_s.loader.exec_module(cm)
import compute_truth_iii as C3  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output" / "a123c" / "diag"
SEED = 20261102
G = {}


def r_off(est, tru):
    te = cm.ev.tissue_errors(est, tru)
    return cm.ev.circ_corr(te.offset_h, te.offset_h_true) if len(te) > 3 else np.nan


def one(i):
    limit_threads(1)
    tab, wf, per = G["tab"], G["wf"], G["per"]
    half = C3.split_donors(per, SEED + i)
    est = {}
    for h in (1, 2):
        keep = {d for d, v in half.items() if v == h}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            est[h] = cm.B.clocksync_estimates(tab[tab.donor.isin(keep)], wf, 0, seed=cm.SEED)
    return dict(rep=i, r_half1_half2=r_off(est[1], est[2]),
                r_half1_truth_i=r_off(est[1], G["ti"]), r_half2_truth_i=r_off(est[2], G["ti"]),
                r_half1_truth_iii=r_off(est[1], G["tiii"]), r_half2_truth_iii=r_off(est[2], G["tiii"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reps", type=int, default=40)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    limit_threads(1)
    import anndata as ad
    x = ad.read_h5ad(cm.LABELLED, backed="r")
    o = x.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy()
    o["tissue"] = o.tissue.astype(str)
    x.file.close()
    tissues = list(cm.A.window_tissues(o).index)
    D = cm.load_clock_data(tissues, False)
    table = pd.concat(D["rows"], ignore_index=True)
    donors = sorted(table.donor[cm.in_window(table.tod)].unique())
    G["wf"] = cm.B.cross_fit_waveform(table, donors)
    G["tab"] = table[table.donor.isin(donors) & cm.in_window(table.tod)]
    G["per"] = {t: sorted(g.donor.unique()) for t, g in G["tab"].groupby("compartment")}
    G["ti"] = pd.read_csv(HERE.parents[1] / "output" / "truth_i.csv", index_col=0)
    t3 = pd.read_csv(HERE.parents[1] / "output" / "a123c" / "truth_iii.csv", index_col=0)
    G["tiii"] = t3.loc[t3.kept, ["amplitude", "offset_h"]].rename_axis("compartment")  # tissues kept in truth (iii)
    full = cm.B.clocksync_estimates(G["tab"], G["wf"], 0, seed=cm.SEED)
    print(f"full map vs truth i {r_off(full, G['ti']):.3f}, vs truth iii {r_off(full, G['tiii']):.3f} "
          f"({len(donors)} donors)", flush=True)
    with Pool(a.jobs) as p:
        res = pd.DataFrame(p.map(one, range(a.reps)))
    rows = []
    for c in res.columns[1:]:
        v = res[c].dropna()
        q = np.quantile(v, [.025, .5, .975])
        rows.append(dict(stat=c, n=len(v), median=q[1], q025=q[0], q975=q[2], mean=v.mean()))
    s = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    s.to_csv(OUT / "estimate_split_half.csv", index=False)
    print(s.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
