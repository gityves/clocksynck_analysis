"""A301 (registered in `paper/docs/registered_experiments.md`): is the adrenal phase effect carried by ACTH activity or by
cortex content?

A299 rows and model (GTEx adrenal, leave-one-out host of the periphery, Hardy 1-4 primary, fast deaths descriptive) with
the ACTH / cAMP-response and cortex-identity scores fitted jointly; the steroidogenic score of A299 for reference.
H301.1: ACTH phase change < 0 (one-sided bootstrap p < 0.05); H301.2: identity phase change within +-0.5 h (90 %).
Aggregates only.

    python score_a301.py
Outputs (`../output/a301/`): scores.csv, correlations.csv, effects.csv, tests.csv.
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import posthoc_a299 as PH  # noqa: E402
import score_a299 as M  # noqa: E402

OUT = HERE.parents[1] / "output" / "a301"
SCORES = {"acth": ["NR4A1", "NR4A2", "NR4A3", "SIK1", "CREM", "MC2R", "MRAP"],
          "identity": ["NR5A1", "NR0B1", "KCNK3", "VSNL1", "CCN3"],
          "steroidogenic": ["STAR", "CYP11A1", "CYP11B1", "CYP17A1", "CYP21A2"]}
BOUND, N_BOOT = 0.5, 2000


def run(hardy, label, rng):
    clock = list(M.cs.gtex().genes)
    markers = sorted({g for v in SCORES.values() for g in v})
    pbm, donors = M.gtex.load(min_donors=M.A104.MIN_DONORS, genes=clock + markers, hardy=hardy)
    pbm = pbm[~pbm.compartment.str.startswith("Brain")].reset_index(drop=True)
    sel = pbm.compartment == M.TISSUE
    ad = pbm[sel].set_index(pbm.donor[sel].astype(str))
    _, _, _, D = M.A292.prepared(pbm[["donor", "compartment", "depth"] + clock], donors.tod, "periphery")
    e = D[M.TISSUE]
    ad = ad.reindex(e["donor"])
    lf = np.log((ad[markers] + 0.5).div(ad.depth, axis=0))
    info, sc = [], {}
    for k, g in SCORES.items():
        c = lf[g].corr().to_numpy()
        v = np.log((ad[g].sum(1) + 0.5) / ad.depth)
        sc[k] = (v - v.mean()) / v.std()
        info.append(dict(donors=label, score=k, mean_pairwise_r=float(c[np.triu_indices(len(g), 1)].mean()), sd_log=float(v.std())))
    sc = pd.DataFrame(sc)
    x, h, G = e["x"], e["h"], e["G"]
    eff = []
    for s in (["acth", "identity"], ["acth", "identity", "steroidogenic"], ["acth"], ["identity"]):
        S = sc[s].to_numpy(float)
        g0, p0 = PH.fit_joint(x, h, S, G)
        bs = [PH.fit_joint(x[i], h[i], S[i], G[i]) for i in (rng.integers(0, len(x), len(x)) for _ in range(N_BOOT))]
        P, Gn = np.array([b[1] for b in bs]), np.array([b[0] for b in bs])
        for j, k in enumerate(s):
            eff.append(dict(donors=label, model="+".join(s), score=k, n=len(x), phase_h=float(p0[j]),
                            phase_lo90=float(np.quantile(P[:, j], 0.05)), phase_hi90=float(np.quantile(P[:, j], 0.95)),
                            p_phase_negative=float((1 + np.sum(P[:, j] >= 0)) / (1 + N_BOOT)), gain=float(g0[j]),
                            gain_lo90=float(np.quantile(Gn[:, j], 0.05)), gain_hi90=float(np.quantile(Gn[:, j], 0.95))))
    return info, sc.corr().assign(donors=label), eff


def main():
    warnings.filterwarnings("ignore")
    M.limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(301)
    res = [run((1, 2, 3, 4), "Hardy 1-4 (primary)", rng), run((1, 2), "fast deaths (descriptive)", rng)]
    info = pd.DataFrame([r for x in res for r in x[0]])
    cor = pd.concat([x[1] for x in res])
    eff = pd.DataFrame([r for x in res for r in x[2]])
    info.to_csv(OUT / "scores.csv", index=False)
    cor.to_csv(OUT / "correlations.csv")
    eff.to_csv(OUT / "effects.csv", index=False)
    j = eff[(eff.donors == "Hardy 1-4 (primary)") & (eff.model == "acth+identity")].set_index("score")
    tests = pd.DataFrame([
        dict(test="H301.1", score="acth", phase_h=j.loc["acth", "phase_h"], p=j.loc["acth", "p_phase_negative"],
             passes=bool(j.loc["acth", "phase_h"] < 0 and j.loc["acth", "p_phase_negative"] < 0.05)),
        dict(test="H301.2", score="identity", phase_h=j.loc["identity", "phase_h"], lo90=j.loc["identity", "phase_lo90"],
             hi90=j.loc["identity", "phase_hi90"],
             passes=bool(-BOUND < j.loc["identity", "phase_lo90"] and j.loc["identity", "phase_hi90"] < BOUND))])
    tests.to_csv(OUT / "tests.csv", index=False)
    print(info.round(3).to_string(index=False))
    print(cor.round(2).to_string())
    print(eff.round(2).to_string(index=False))
    print(tests.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
