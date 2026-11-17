"""A289 (registered in `paper/docs/registered_experiments.md`): is the brain's lower coupling explained by the blood-brain
barrier, and do brain clocks couple more than their amplitude predicts? Baboon (primary) and sheep (replication).

Reads the A42 readout-2 table (`results/figure_2/output/organs/r2_couplings.csv`: coupling of every tissue to the
peripheral-anchor host, own amplitude); no refit. OUT = pineal and pituitary; IN = the A42 brain regions.
  H289.1 / H289.2: Delta = median coupling(OUT) - median coupling(IN), exact one-sided permutation of the OUT labels.
  H289.3: residuals of the brain regions from the peripheral line coupling = a + b log(own amplitude); median residual,
          one-sided sign test; extrapolation check (peripheral tissues within the brain amplitude range).

    python score_a289.py
Outputs (`../output/`): tissues.csv (every tissue: species, group, coupling, own amplitude, residual), tests.csv.
"""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
OUT = HERE.parents[1] / "output"
SRC = PAPER / "results/figure_2/output/organs/r2_couplings.csv"
GLANDS = {"baboon": ("PIN", "PIT"), "sheep": ("Pineal", "Pituitary")}
OCULAR = ("COR", "ONH", "RET", "RPE")


def barrier(t, glands, inside):
    """Delta and its exact one-sided permutation p over the assignment of the OUT labels."""
    pool = t.loc[list(glands) + list(inside), "coupling"]
    k = len(glands)
    delta = pool[list(glands)].median() - pool[list(inside)].median()
    null = [pool.iloc[list(c)].median() - pool.drop(pool.index[list(c)]).median() for c in combinations(range(len(pool)), k)]
    return float(delta), float(np.mean(np.asarray(null) >= delta - 1e-12)), len(null)


def amplitude_line(t):
    per = t[t.group == "periphery"]
    b, a = np.polyfit(np.log(per.own_amplitude), per.coupling, 1)
    res = t.coupling - (a + b * np.log(t.own_amplitude))
    return res, a, b


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(SRC)
    rows, tabs = [], []
    for sp in ("baboon", "sheep"):
        t = d[d.species == sp].set_index("compartment")
        t["group"] = np.where(t.index.isin(GLANDS[sp]), "outside barrier",
                              np.where(t.index.isin(OCULAR), "ocular", t.group))
        inside = list(t.index[t.group == "brain"])
        delta, p, n = barrier(t, GLANDS[sp], inside)
        rows.append(dict(test="H289.1" if sp == "baboon" else "H289.2", species=sp, statistic="Delta", value=delta, p=p,
                         n_perm=n, n_in=len(inside), median_out=t.loc[list(GLANDS[sp]), "coupling"].median(),
                         median_in=t.loc[inside, "coupling"].median(), passed=bool(delta > 0 and p < 0.05)))
        if sp == "baboon":
            oc = list(t.index[t.group == "ocular"])
            dd, pp, nn = barrier(t, GLANDS[sp], oc)
            rows.append(dict(test="descriptive", species=sp, statistic="Delta vs ocular", value=dd, p=pp, n_perm=nn,
                             n_in=len(oc), median_out=t.loc[list(GLANDS[sp]), "coupling"].median(),
                             median_in=t.loc[oc, "coupling"].median(), passed=np.nan))
        res, a, b = amplitude_line(t)
        t["residual"] = res
        r = res[inside]
        lo, hi = t.loc[inside, "own_amplitude"].min(), t.loc[inside, "own_amplitude"].max()
        per = t[t.group == "periphery"]
        within = int(per.own_amplitude.between(lo, hi).sum())
        sign_p = binomtest(int((r > 0).sum()), len(r), 0.5, alternative="greater").pvalue
        rows.append(dict(test="H289.3" if sp == "baboon" else "H289.3 sheep (reported)", species=sp,
                         statistic="median brain residual", value=float(r.median()), p=float(sign_p), n_in=len(r),
                         above=int((r > 0).sum()), slope=b, intercept=a, brain_amp_lo=lo, brain_amp_hi=hi,
                         periph_amp_lo=per.own_amplitude.min(), periph_amp_hi=per.own_amplitude.max(),
                         periph_within_brain_range=within, extrapolated=within < 3,
                         passed=bool(r.median() > 0 and sign_p < 0.05)))
        tabs.append(t.assign(species=sp)[["species", "group", "n", "coupling", "own_amplitude", "residual"]])
    pd.concat(tabs).rename_axis("tissue").to_csv(OUT / "tissues.csv")
    tests = pd.DataFrame(rows)
    tests.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 220)
    print(tests.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
