"""Scout for a within-tissue proliferation-loudness test (before registration; no clock quantity computed).

GTEx fast deaths, A104 periphery rows (the A295 rows). Proliferation score per sample: mean log CPM of the Hallmark G2M
genes (clock genes removed; A127 / A391 list). Per tissue: within-tissue SD and 10-90 % fold range of the score;
split-half reliability (G2M genes split at random into two halves, 20 splits, Spearman-Brown of the within-tissue
correlation); correlation with the A295 lymphocyte and stromal-vascular scores, age, RIN and ischaemic time; and the mean
score difference (in within-tissue SD) of samples whose pathology notes name a lesion, by keyword class.
Aggregates only (GTEx controlled access: no per-sample row is written).

    python scout.py
Outputs (`../output/scout/`): tissues.csv, pathology.csv.
"""
import importlib.util
import re
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_2" / "scripts" / "backbone"))
sys.path.insert(0, str(PAPER / "exploratory" / "A143_brain_artefact" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_brain_artefact_checks as A143  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output" / "scout"
LYMPH = ["CD3D", "CD3E", "CD2", "CD247", "MS4A1", "CD79A", "CD19"]
STROMA = ["DCN", "LUM", "COL1A2", "PECAM1", "VWF", "CDH5"]
LESIONS = {"inflammation": r"inflam|itis\b|gastritis|colitis|esophagitis|lymphocytic infiltrat",
           "hyperplasia_metaplasia": r"hyperplas|metaplas|barrett|dysplas",
           "neoplasm": r"adenoma|carcinoma|tumou?r|neoplas|malignan|polyp",
           "fibrosis_atrophy": r"fibros|atroph|cirrho|sclero"}


def g2m():
    spec = importlib.util.spec_from_file_location("fig5_shared", PAPER / "results/figure_5/scripts/cancer/_shared.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    clock = {canonical(g) for g in cs.gtex().genes}
    return [g for g in m.hallmark_g2m() if canonical(g) not in clock]


def score(pb, genes):
    g = [x for x in genes if x in pb.columns]
    return np.log1p(pb[g].div(pb.depth, axis=0) * 1e6).mean(axis=1)


def main():
    warnings.filterwarnings("ignore")
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    G = g2m()
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    var = set(a.var_names)
    a.file.close()
    G = [g for g in G if g in var]
    pb, donors = gtex.load(min_donors=A104.MIN_DONORS, genes=G + LYMPH + STROMA)
    pb = pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)
    pb["donor"] = pb.donor.astype(str)
    att = A143.sample_attributes(pb)
    notes = pd.read_csv(A143.ATTRIBUTES, sep="\t", comment="#", low_memory=False, usecols=["SAMPID", "SMPTHNTS"])
    pb["notes"] = att.SAMPID.map(notes.set_index("SAMPID").SMPTHNTS).fillna("").str.lower().to_numpy()
    pb["prolif"] = score(pb, G)
    pb["lymph"] = np.log((pb[LYMPH].sum(1) + 0.5) / pb.depth)
    pb["stroma"] = np.log((pb[STROMA].sum(1) + 0.5) / pb.depth)
    pb["rin"] = att.SMRIN.to_numpy(float)
    dn = donors.copy()
    dn.index = dn.index.astype(str)
    pb["age"], pb["isch"] = pb.donor.map(dn.age).astype(float), pb.donor.map(dn.ischaemic_h).astype(float)
    rng = np.random.default_rng(398)
    halves = []
    for _ in range(20):
        p = rng.permutation(G)
        halves.append((list(p[: len(p) // 2]), list(p[len(p) // 2:])))
    rows, path = [], []
    for k, t in pb.groupby("compartment"):
        s = t.prolif
        rel = []
        for h1, h2 in halves:
            r = np.corrcoef(score(t, h1), score(t, h2))[0, 1]
            rel.append(2 * r / (1 + r))
        sd = float(s.std(ddof=1))
        row = dict(tissue=k, n=len(t), mean=float(s.mean()), sd=sd,
                   fold_10_90=float(np.exp(np.quantile(s, 0.9) - np.quantile(s, 0.1))),
                   split_half_reliability=float(np.median(rel)))
        for c in ("lymph", "stroma", "age", "rin", "isch"):
            ok = t[c].notna()
            row[f"rho_{c}"] = float(spearmanr(s[ok], t[c][ok]).statistic) if ok.sum() >= 10 else np.nan
        rows.append(row)
        z = (s - s.mean()) / sd
        for lab, pat in LESIONS.items():
            m = t.notes.str.contains(pat, regex=True).to_numpy()
            if m.sum() >= 5:
                path.append(dict(tissue=k, lesion=lab, n_lesion=int(m.sum()), n_other=int((~m).sum()),
                                 diff_sd=float(z[m].mean() - z[~m].mean())))
    T = pd.DataFrame(rows).sort_values("sd", ascending=False)
    T["between_tissue_sd_of_means"] = float(T["mean"].std(ddof=1))
    P = pd.DataFrame(path)
    T.to_csv(OUT / "tissues.csv", index=False)
    P.to_csv(OUT / "pathology.csv", index=False)
    pd.set_option("display.width", 250)
    print(f"G2M genes {len(G)}; {pb.donor.nunique()} donors, {pb.compartment.nunique()} tissues")
    print(T.round(3).to_string(index=False))
    print(P.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
