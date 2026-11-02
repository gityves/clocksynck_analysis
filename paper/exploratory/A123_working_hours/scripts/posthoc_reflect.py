"""A123b post hoc diagnostic (2026-11-02, not registered): A123b orderer organ maps in design W vs truth (i), allowing a global reflection and a
stretch of the offsets (orderers spread a 10 h window over 24 h). Tissue-level aggregates only."""
import sys
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
D0 = Path(__file__).resolve().parent
sys.path.insert(0, str(D0))
spec = importlib.util.spec_from_file_location("cm", D0 / "compute_methods.py")
cm = importlib.util.module_from_spec(spec); sys.argv = ["x"]; spec.loader.exec_module(cm)
import anndata as ad
a = ad.read_h5ad(cm.LABELLED, backed="r")
o = a.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy(); o["tissue"] = o.tissue.astype(str); a.file.close()
tissues = list(cm.A.window_tissues(o).index)
D = cm.load_clock_data(tissues, False)
tod_all = pd.concat(list(D["tod"].values())).groupby(level=0).first()
truth = cm.ev.csvd_tissues(pd.DataFrame(D["Ci"]), cm.REF)
OUT = D0.parent / "output" / "a123b"
runs = OUT / "runs"
table = pd.concat(D["rows"], ignore_index=True)
rows, maps = [], {}
def centred(x):  # hours, centred on circular mean, in (-12, 12]
    m = np.angle(np.mean(np.exp(1j * x * np.pi / 12))) * 12 / np.pi
    return (x - m + 12) % 24 - 12
for m in ["clocksync", "chiral", "cofe_clock", "cofe_genome", "cyclops"]:
    for d in ["W", "F"]:
        est = cm.clocksync_est(d, table, runs)[0] if m == "clocksync" else cm.assemble(m, d, tissues, D, runs, tod_all)[0]
        te = cm.ev.tissue_errors(est, truth).dropna(subset=["offset_h", "offset_h_true"])
        x, y = centred(te.offset_h.to_numpy()), centred(te.offset_h_true.to_numpy())
        r = cm.ev.circ_corr(te.offset_h, te.offset_h_true)
        rr = cm.ev.circ_corr(-te.offset_h, te.offset_h_true)
        sgn = 1 if r >= rr else -1
        xs = sgn * x
        k = float(np.sum(xs * y) / np.sum(y * y))  # least-squares stretch of estimate on truth
        rows.append(dict(method=m, design=d, n=len(te), circ_r=r, circ_r_reflected=rr,
                         spearman_best=pd.Series(xs).corr(pd.Series(y), method="spearman"),
                         pearson_best=float(np.corrcoef(xs, y)[0, 1]), stretch=k,
                         sd_est_h=float(np.std(x)), sd_truth_h=float(np.std(y))))
        maps[m, d] = pd.DataFrame({"est_h": xs, "truth_h": y}, index=te.index)
out = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(out.round(3).to_string(index=False))
out.to_csv(OUT / "posthoc_reflect.csv", index=False)
pd.concat(maps, names=["method", "design", "tissue"]).to_csv(OUT / "posthoc_reflect_maps.csv")
