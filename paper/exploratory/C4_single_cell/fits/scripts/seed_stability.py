"""Seed stability of the CoPhaser cell-cycle phase (QC). Per study: circular correlation (Jammalamadaka-SenGupta) and median
absolute angular difference between the kept seed and every other usable seed (main run seeds 0-4, extra seeds 5-9 where
run), over the cells the kept seed calls cycling (p > 0.5); both thetas are after the package's orientation. Also the
kept seed vs the pilot's cells_<study>.csv.gz theta (reproduction).   python -I seed_stability.py"""
import json
from pathlib import Path
import numpy as np, pandas as pd

OUT = Path(__file__).resolve().parents[1] / "output"
PILOT = OUT.parents[1] / "pilot/output"


def wrap(x):
    return np.arctan2(np.sin(x), np.cos(x))


def circcorr(a, b):
    sa, sb = np.sin(a - np.angle(np.exp(1j * a).sum())), np.sin(b - np.angle(np.exp(1j * b).sum()))
    return float((sa * sb).sum() / np.sqrt((sa ** 2).sum() * (sb ** 2).sum()))


rows = []
for qf in sorted(OUT.glob("qc_*.json")):
    q = json.load(open(qf)); st = q["study"]; ks = q["kept_seed"]
    main = np.load(OUT / f"seed_theta_{st}.npz")
    th = {("main", int(k.split("_")[1])): main[k] for k in main.files if k.startswith("theta_")}
    pc = {("main", int(k.split("_")[1])): main[k] for k in main.files if k.startswith("pcyc_")}
    cells = main["cells"]
    ex = OUT / f"seed_theta_{st}_extra.npz"
    if ex.exists():
        e = np.load(ex); assert (e["cells"] == cells).all()
        th.update({("extra", int(k.split("_")[1])): e[k] for k in e.files if k.startswith("theta_")})
        pc.update({("extra", int(k.split("_")[1])): e[k] for k in e.files if k.startswith("pcyc_")})
    if ("main", ks) not in th:  # every seed vetoed: no usable seed, nothing to compare
        rows.append(dict(study=st, kept_seed=ks, set="none_usable"))
        continue
    kth, kpc = th[("main", ks)], pc[("main", ks)]
    use = kpc > 0.5
    if use.sum() < 50:
        use = np.ones(len(kth), bool)
    base = dict(study=st, kept_seed=ks, n_cells_used=int(use.sum()), cycling_cells_only=bool((kpc > 0.5).sum() >= 50))
    pl = pd.read_csv(PILOT / f"cells_{st}.csv.gz", index_col=0) if (PILOT / f"cells_{st}.csv.gz").exists() else None
    if pl is not None:
        pl = pl.reindex(cells)
        ok = pl.theta.notna().to_numpy()
        pu = use & ok
        rows.append({**base, "set": "pilot_cells_csv", "seed": q["kept_seed"], "circ_corr": circcorr(kth[pu], pl.theta.to_numpy()[pu]),
                     "med_abs_diff_rad": float(np.median(np.abs(wrap(kth[pu] - pl.theta.to_numpy()[pu])))), "n_cells_compared": int(pu.sum()),
                     "frac_cells_in_pilot": float(ok.mean())})
    for (s, k), t in sorted(th.items()):
        if s == "main" and k == ks:
            continue
        rows.append({**base, "set": s, "seed": k, "circ_corr": circcorr(kth[use], t[use]),
                     "med_abs_diff_rad": float(np.median(np.abs(wrap(kth[use] - t[use])))), "n_cells_compared": int(use.sum())})
df = pd.DataFrame(rows)
df.to_csv(OUT / "seed_stability.csv", index=False)
pd.set_option("display.width", 200)
print(df.round(3).to_string())
