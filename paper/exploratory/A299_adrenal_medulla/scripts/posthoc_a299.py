"""A299 post hoc (not registered; chosen after H299.1 was scored, at the user's question whether the cortex effect is the
converse of the medulla hypothesis with clearer markers): what is the non-cortex part of a GTEx adrenal sample?

Per adrenal sample (Hardy 1-4, the A299 rows): marker scores (medulla, cortex, fat, capsule / stroma, blood), their
internal consistency (mean pairwise correlation of the log marker fractions) and their correlations; then the A299
phase / gain model for the extra scores, and for medulla and cortex fitted jointly. Aggregates only.

    python posthoc_a299.py
Outputs (`../output/posthoc/`): scores.csv (consistency, spread), correlations.csv, effects.csv.
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a299 as M  # noqa: E402

OUT = HERE.parents[1] / "output" / "posthoc"
SCORES = dict(M.SCORES, fat=["ADIPOQ", "PLIN1", "LEP", "CIDEA", "CIDEC"], stroma=["DCN", "LUM", "COL1A2", "COL3A1"],
              blood=["HBB", "HBA1", "HBA2"])


def fit_joint(x, h, S, G):
    Z = np.column_stack([G, S, h[:, None], S * h[:, None]])
    coef, *_ = np.linalg.lstsq(Z, x, rcond=None)
    k = S.shape[1]
    a, b = coef[-k - 1], coef[-k:]
    r = b / a
    return np.real(r), np.imag(r) / M.A295.W


def main():
    warnings.filterwarnings("ignore")
    M.limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    clock = list(M.cs.gtex().genes)
    markers = sorted({g for v in SCORES.values() for g in v})
    pbm, donors = M.gtex.load(min_donors=M.A104.MIN_DONORS, genes=clock + markers, hardy=(1, 2, 3, 4))
    pbm = pbm[~pbm.compartment.str.startswith("Brain")].reset_index(drop=True)
    ad = pbm[pbm.compartment == M.TISSUE].set_index(pbm.donor[pbm.compartment == M.TISSUE].astype(str))
    _, _, _, D = M.A292.prepared(pbm[["donor", "compartment", "depth"] + clock], donors.tod, "periphery")
    e = D[M.TISSUE]
    ad = ad.reindex(e["donor"])
    lf = np.log((ad[markers] + 0.5).div(ad.depth, axis=0))
    rows, sc = [], {}
    for k, g in SCORES.items():
        c = lf[g].corr().to_numpy()
        v = np.log((ad[g].sum(1) + 0.5) / ad.depth)
        sc[k] = (v - v.mean()) / v.std()
        rows.append(dict(score=k, genes=len(g), mean_pairwise_r=float(c[np.triu_indices(len(g), 1)].mean()),
                         sd_log_fraction=float(v.std()), median_cpm=float(np.exp(v.median()) * 1e6)))
    pd.DataFrame(rows).to_csv(OUT / "scores.csv", index=False)
    sc = pd.DataFrame(sc)
    sc.corr().to_csv(OUT / "correlations.csv")
    x, h, G = e["x"], e["h"], e["G"]
    eff = []
    rng = np.random.default_rng(2990)
    sets = [[k] for k in SCORES] + [["medulla", "cortex"], ["cortex", "fat"], ["cortex", "medulla", "fat", "stroma"]]
    for s in sets:
        S = sc[s].to_numpy(float)
        g0, p0 = fit_joint(x, h, S, G)
        bs = [fit_joint(x[i], h[i], S[i], G[i]) for i in (rng.integers(0, len(x), len(x)) for _ in range(2000))]
        P = np.array([b[1] for b in bs])
        Gn = np.array([b[0] for b in bs])
        for j, k in enumerate(s):
            eff.append(dict(model="+".join(s), score=k, phase_h=float(p0[j]), phase_lo90=float(np.quantile(P[:, j], 0.05)),
                            phase_hi90=float(np.quantile(P[:, j], 0.95)), gain=float(g0[j]),
                            gain_lo90=float(np.quantile(Gn[:, j], 0.05)), gain_hi90=float(np.quantile(Gn[:, j], 0.95))))
    pd.DataFrame(eff).to_csv(OUT / "effects.csv", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    print(sc.corr().round(2).to_string())
    print(pd.DataFrame(eff).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
