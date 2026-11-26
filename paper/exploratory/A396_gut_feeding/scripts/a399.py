"""A399 (descriptive; registered 2026-11-26): is the quiet sheep digestive tract lymphocyte dilution?

Per tissue: median lymphocyte and plasma-cell marker scores (A295 markers; log((sum + 0.5) / depth)); tissue log gain from
the A383 fits. OLS of log gain on the lymphocyte score (and on both scores) over the peripheral tissues, with and without
the lymphoid organs; digestive minus other mean log gain, raw and on the residuals. Sheep, and baboon as a comparison.

    python a399.py
Outputs (`../output/`): a399_tissues.csv, a399_summary.csv. Aggregates only.
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory" / "A383_species_organisation" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a383_a384 as B  # noqa: E402

OUT = HERE.parents[1] / "output"
LYMPH = ["CD3D", "CD3E", "CD2", "CD247", "MS4A1", "CD79A", "CD19"]
PLASMA = ["JCHAIN", "MZB1"]
DIGESTIVE = {"sheep": ["Rumen", "Reticulum", "Omasum", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum"],
             "baboon": ["OES", "STF", "DUO", "ILE", "CEC", "ASC", "DEC"]}
LYMPHOID = {"sheep": ["Lymph node", "Thymus", "Spleen", "Bone marrow"], "baboon": ["MEL", "SPL", "BOM"]}


def scores(species, per):
    pb = B.A.table(species, LYMPH + PLASMA)
    pb = pb[pb.compartment.isin(per)]
    out = {}
    for name, g in (("lymph", LYMPH), ("plasma", PLASMA)):
        v = np.log((pb[g].sum(1) + 0.5) / pb.depth)
        out[name] = v.groupby(pb.compartment).median()
    return pd.DataFrame(out)


def ols_resid(y, X):
    Z = np.column_stack([np.ones(len(y))] + [X])
    b = np.linalg.lstsq(Z, y, rcond=None)[0]
    return y - Z @ b, b[1:]


def summarise(species, t):
    dig = t.index.isin(DIGESTIVE[species])
    rows = [dict(species=species, model="raw", slopes="", digestive_minus_other=float(t.gain[dig].mean() - t.gain[~dig].mean()),
                 n_digestive=int(dig.sum()), n_tissues=len(t))]
    for lab, cols, fit_mask in (("lymph", ["lymph"], None), ("lymph + plasma", ["lymph", "plasma"], None),
                                ("lymph, fitted without lymphoid organs", ["lymph"], ~t.index.isin(LYMPHOID[species]))):
        X = t[cols].to_numpy(float)
        m = np.ones(len(t), bool) if fit_mask is None else np.asarray(fit_mask)
        Z = np.column_stack([np.ones(m.sum()), X[m]])
        b = np.linalg.lstsq(Z, t.gain.to_numpy(float)[m], rcond=None)[0]
        res = t.gain.to_numpy(float) - np.column_stack([np.ones(len(t)), X]) @ b
        rows.append(dict(species=species, model=lab, slopes=";".join(f"{c} {v:.3f}" for c, v in zip(cols, b[1:])),
                         digestive_minus_other=float(res[dig].mean() - res[~dig].mean()), n_digestive=int(dig.sum()),
                         n_tissues=len(t), r2=float(1 - np.var(res[m]) / np.var(t.gain.to_numpy(float)[m]))))
    return rows


def main():
    warnings.filterwarnings("ignore")
    F = B.fits()
    tabs, rows = [], []
    for sp in ("sheep", "baboon"):
        tab, per = F[sp]
        t = scores(sp, per).join(tab.log_amp_rel.rename("gain")).dropna()
        t["digestive"] = t.index.isin(DIGESTIVE[sp])
        t["lymphoid"] = t.index.isin(LYMPHOID[sp])
        tabs.append(t.assign(species=sp).rename_axis("tissue").reset_index())
        rows += summarise(sp, t)
    T, S = pd.concat(tabs, ignore_index=True), pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT / "a399_tissues.csv", index=False)
    S.to_csv(OUT / "a399_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(T.sort_values(["species", "lymph"], ascending=[True, False]).round(3).to_string(index=False))
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
