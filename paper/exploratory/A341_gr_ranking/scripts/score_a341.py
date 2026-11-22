"""A341 (registered in `paper/docs/registered_experiments.md`; power note: detectable rho 0.90 at 7 organs, so reported
descriptively): glucocorticoid responsiveness vs clock maturity across the fetal organs.

Units and rows: the organ readouts as scored (A244 adrenal cortex and medulla, A243 `all` lung, liver, heart, kidney,
A239 `brain_mean:prenatal`). Per organ, residualised as the readout (cohort and quadratic age; brain: log age and sex):
(i) NR3C1 median log CPM; (ii) between-fetus variance of the A307 GC-target score (mean residual log fraction of FKBP5,
TSC22D3, ZBTB16, KLF9) minus its counting variance (exact log-count noise, `lognoise`; brain: arrays, no correction);
(iii) median log CPM of the 11 core genes; (iv) ARNTL - NR1D1 median log CPM. Outcome: the organ's C (scored l1 tables).
(i), (iii), (iv) are not computed for the brain (arrays: levels not comparable with counts). Spearman of the outcome
with each measure, one-sided permutation p over all orders (descriptive).

    python score_a341.py
Outputs (`../output/`): measures.csv, spearman.csv.
"""
import sys
from itertools import permutations
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP.parent))
sys.path.insert(0, str(EXP / "A328_fetal_running_loop" / "scripts"))
sys.path.insert(0, str(EXP / "A239_fetal_brain" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import lognoise  # noqa: E402
from common import development as dev  # noqa: E402

DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/natexp")
OUT = HERE.parents[1] / "output"
GC = ["FKBP5", "TSC22D3", "ZBTB16", "KLF9"]
CORE = ["ARNTL", "CLOCK", "NPAS2", "PER3", "CRY1", "CRY2", "NR1D2", "DBP", "RORA", "TEF", "CIART"]
OUTCOME = {"adrenal cortex": ("A244_fetal_adrenal_classes", "cortex"), "adrenal medulla": ("A244_fetal_adrenal_classes", "medulla"),
           "lung": ("A243_fetal_organs_pooled", "lung:all"), "liver": ("A243_fetal_organs_pooled", "liver:all"),
           "heart": ("A243_fetal_organs_pooled", "heart:all"), "kidney": ("A243_fetal_organs_pooled", "kidney:all"),
           "brain": ("A239_fetal_brain", "brain_mean:prenatal")}


def count_measures(x: pd.DataFrame, obs: pd.DataFrame) -> dict:
    """x: fetus x gene counts (NaN = gene absent from a cohort, set to 0 only for the summed depth), obs: cohort, age, depth."""
    depth = obs.depth.to_numpy(float)
    Z = dev.design(obs, "age", 2, ("cohort",))
    lcpm = np.log((x.fillna(0).to_numpy(float) + 0.5) / depth[:, None] * 1e6)
    L = pd.DataFrame(lcpm, columns=x.columns)
    g = x[GC].fillna(0).to_numpy(float)
    y = dev.residual(np.log((g + 0.5) / depth[:, None]), Z)
    score = y.mean(1)
    dof = 1 - (1 + Z.shape[1]) / len(y)
    noise = dof * lognoise.var_log(lognoise.expected(g, depth, obs.cohort.to_numpy())).mean(0).sum() / len(GC) ** 2
    return dict(n=len(x), NR3C1=float(L["NR3C1"].median()), gc_var=float(score.var(ddof=1 + Z.shape[1]) - noise),
                gc_noise=float(noise), core_level=float(L[CORE].median().median()),
                arntl_minus_nr1d1=float(L["ARNTL"].median() - L["NR1D1"].median()))


def main():
    OUT.mkdir(exist_ok=True)
    rows = []
    s = pd.read_csv(DATA / "fetal_organs_pooled" / "samples.csv")
    s = s[s.eligible]
    pb = pd.read_parquet(DATA / "fetal_organs_pooled" / "pseudobulk.parquet")
    for organ in ("lung", "liver", "heart", "kidney"):
        o = s[s.organ == organ].rename(columns={"age_weeks": "age"})
        rows.append(dict(organ=organ, **count_measures(pb.loc[o.id], o.reset_index(drop=True))))
    a = pd.read_csv(DATA / "fetal_adrenal_open" / "a244_class_samples.csv")
    q = pd.read_parquet(DATA / "fetal_adrenal_open" / "a244_class_pseudobulk.parquet")
    for cls in ("cortex", "medulla"):
        o = a[(a.compartment == cls) & (a.cells >= 200)]  # the A244 >= 200-cell rule
        rows.append(dict(organ=f"adrenal {cls}", **count_measures(q.loc[o.id], o.reset_index(drop=True))))
    import score_a239 as b  # brain arrays (Kang 2011): GC variance only
    x, obs = b.load()
    lo, hi = b.GROUPS["prenatal"]
    m = (obs.compartment == "brain_mean") & obs.stage.between(lo, hi)
    xb, ob = x[m.to_numpy()], obs[m].reset_index(drop=True)
    Z = b.cov(2)(ob)
    gb = [g for g in GC if g in xb.columns]  # FKBP5 is not on the arrays
    y = dev.residual(np.log(xb[gb].to_numpy(float)), Z)
    rows.append(dict(organ="brain", n=len(xb), gc_var=float(y.mean(1).var(ddof=1 + Z.shape[1])),
                     gc_genes=";".join(gb)))
    t = pd.DataFrame(rows)
    t["C"] = [pd.read_csv(EXP / OUTCOME[o][0] / "output" / "l1.csv").set_index("test").loc[OUTCOME[o][1], "C"]
              for o in t.organ]
    t.to_csv(OUT / "measures.csv", index=False)
    sp = []
    for col in ("NR3C1", "gc_var", "core_level", "arntl_minus_nr1d1"):
        u = t.dropna(subset=[col])
        k = len(u)
        r = spearmanr(u.C, u[col])[0]
        null = np.array([spearmanr(np.arange(k), p)[0] for p in permutations(range(k))])
        sp.append(dict(measure=col, organs=k, rho=r, p_one_sided=float(np.mean(null >= r - 1e-12))))
    sp = pd.DataFrame(sp)
    sp.to_csv(OUT / "spearman.csv", index=False)
    pd.set_option("display.width", 200)
    print(t.round(3).to_string(index=False), sp.round(3).to_string(index=False), sep="\n")


if __name__ == "__main__":
    main()
