"""A323 (registered): do cell-class gains and composition predict the GTEx tissue gains?

Composition: per GTEx tissue of the Eraslan 2022 single-nucleus atlas (clock-gene pseudobulk per donor x tissue x cell
type, `$SCCC_DATA/atlas/gtex_eraslan/pseudobulk.h5ad`), each family's share of the summed 14 clock-gene counts.
Prediction G_t = sum_c p_c exp(g_c) with the A293 Part 3 pooled class gains (lymphocyte, other: 0). Observed sqrt(S_t),
S = z2 x V of the A104 row. H323.1: Spearman rho > 0, exact one-sided permutation p over the 8 tissues.

    python score_a323.py
Outputs (`../output/`): shares.csv, tissues.csv, tests.csv.
"""
import itertools
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, theilslopes

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
OUT = HERE.parents[1] / "output"
DATA = Path(os.environ.get("SCCC_DATA", "/mnt/max8tb/data/CoPhaser_data/cancer_circadian"))
SRC = DATA / "atlas" / "gtex_eraslan" / "pseudobulk.h5ad"
A104 = PAPER / "results/figure_2/output/backbone/a104_gtex_tissues.csv"
A293 = PAPER / "exploratory/A293_organisation_levels/output/part3/families.csv"

TISSUE = {"prostate gland": "Prostate", "anterior wall of left ventricle": "Heart - Left Ventricle",
          "mucosa": "Esophagus - Mucosa", "gastrocnemius": "Muscle - Skeletal",
          "esophagus muscularis mucosa": "Esophagus - Muscularis", "lingula of left lung": "Lung",
          "skin of leg": "Skin - Sun Exposed (Lower leg)", "breast": "Breast - Mammary Tissue"}
FAMILY = {
    "fibroblast_stromal": ["fibroblast", "skeletal muscle fibroblast", "skin fibroblast", "fibroblast of breast",
                           "fibroblast of cardiac tissue"],
    "endothelial": ["endothelial cell of vascular tree", "endothelial cell of lymphatic vessel", "cardiac endothelial cell"],
    "mural": ["enteric smooth muscle cell", "smooth muscle cell", "smooth muscle cell of prostate", "contractile cell"],
    "parenchyma": ["skeletal muscle fiber", "fast muscle cell", "slow muscle cell", "cardiac muscle cell",
                   "subcutaneous adipocyte", "adipocyte"],
    "myeloid": ["macrophage", "alveolar macrophage", "elicited macrophage", "professional antigen presenting cell",
                "dendritic cell", "myeloid dendritic cell", "mast cell", "neutrophil", "epidermal Langerhans cell"],
    "lymphocyte": ["T cell", "mature NK T cell", "B cell"],
    "other": ["Schwann cell", "skeletal muscle satellite cell", "melanocyte of skin", "interstitial cell of Cajal",
              "neuroendocrine cell", "neuron"],
}
EPI_WORDS = ("epithelial", "basal", "luminal", "keratinocyte", "epidermal", "sebaceous", "sweat gland",
             "glandular secretory", "squamous")


def family_of(ct):
    for f, names in FAMILY.items():
        if ct in names:
            return f
    if any(w in ct for w in EPI_WORDS):
        return "epithelial"
    raise KeyError(f"cell type without a family: {ct}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(SRC)
    X = np.asarray(a.X.todense() if hasattr(a.X, "todense") else a.X, float)
    o = a.obs.assign(tissue=a.obs.tissue.astype(str), cell_type=a.obs.cell_type.astype(str))
    o["family"] = o.cell_type.map(family_of)
    o["clock"], o["cells"] = X.sum(1), o.n_cells.astype(float)
    sh = o.groupby(["tissue", "family"]).agg(clock=("clock", "sum"), cells=("cells", "sum")).reset_index()
    sh["p_clock"] = sh.clock / sh.groupby("tissue").clock.transform("sum")
    sh["p_cells"] = sh.cells / sh.groupby("tissue").cells.transform("sum")
    sh.to_csv(OUT / "shares.csv", index=False)
    g = pd.read_csv(A293).set_index("family").pooled.to_dict()
    eg = {f: float(np.exp(g.get(f, 0.0))) for f in FAMILY if f not in ("lymphocyte", "other")}
    eg.update(endothelial=1.0, lymphocyte=0.0, other=0.0)
    gt = pd.read_csv(A104).set_index("tissue")
    rows = []
    for t, x in sh.groupby("tissue"):
        w = x.family.map(eg).to_numpy()
        r = gt.loc[TISSUE[t]]
        keep = ~x.family.isin(["other"]).to_numpy()
        rows.append(dict(tissue=TISSUE[t], G=float(np.sum(x.p_clock * w)),
                         G_other_renorm=float(np.sum(x.p_clock[keep] * w[keep]) / x.p_clock[keep].sum()),
                         G_cells=float(np.sum(x.p_cells * w)), sqrtS=float(np.sqrt(r.z2 * r.V)),
                         nuclei=float(x.cells.sum())))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "tissues.csv", index=False)

    def exact(col):
        obs = spearmanr(T[col], T.sqrtS).statistic
        null = [spearmanr(T[col].to_numpy()[list(p)], T.sqrtS).statistic for p in itertools.permutations(range(len(T)))]
        return float(obs), float(np.mean(np.array(null) >= obs - 1e-12))

    res = []
    for col in ("G", "G_other_renorm", "G_cells"):
        rho, p = exact(col)
        lg, ls = np.log(T[col]), np.log(T.sqrtS)
        sl = theilslopes(ls, lg)[0]
        r2 = float(np.corrcoef(lg, ls)[0, 1] ** 2)
        res.append(dict(prediction=col, rho=rho, p_exact=p, theil_sen_slope=float(sl), r2_log=r2,
                        range_pred=float(T[col].max() / T[col].min()), range_obs=float(T.sqrtS.max() / T.sqrtS.min()),
                        H323_1=(p < 0.05 and rho > 0) if col == "G" else None))
    R = pd.DataFrame(res)
    R.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 220)
    print(sh.pivot(index="tissue", columns="family", values="p_clock").round(3).to_string())
    print(T.round(3).to_string(index=False))
    print(R.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
