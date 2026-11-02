"""A190 post hoc checks (not registered; chosen after A190 was scored, 2026-11-02).

1. Ambient / doublet contamination. sci-RNA-seq3 profiles of shallow classes can carry RNA of the dominant class. Per
   fetus, the contamination fraction a of a class from the parenchyma is the fraction of parenchyma-specific marker
   counts in that class relative to the parenchyma (adrenal: STAR, CYP11A1, CYP17A1, CYP21A2, CYP11B1, NR5A1; heart, as a
   comparison: MYH6, MYH7, TNNT2, NPPA, MYL7). Expected parenchyma-derived counts of a clock gene in the class are
   a x class depth x parenchyma fraction of the gene; their share of the observed counts is reported. The within-adrenal
   gate is rerun on corrected counts (expected contamination subtracted, floored at 0, rounded), 10,000 permutations.
2. Receptor expression (descriptive): counts per million of MTNR1A, MTNR1B (melatonin), NR3C1 (glucocorticoid), MC2R
   (ACTH) per organ and per class, fetuses pooled (median over fetuses).

    python posthoc_checks.py [--n-perm 10000]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[3]), str(HERE.parents[2] / "A189_fetal_within" / "scripts")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from score_a189 import GENES, WF, gate, genome, table  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 1901
MARKERS = {"Adrenal": ["STAR", "CYP11A1", "CYP17A1", "CYP21A2", "CYP11B1", "NR5A1"],
           "Heart": ["MYH6", "MYH7", "TNNT2", "NPPA", "MYL7"]}
OTHERS = ["stromal", "endothelial"]
RECEPTORS = ["MTNR1A", "MTNR1B", "NR3C1", "MC2R"]


def contamination(organ):
    """Per fetus and shallow class: marker-based contamination fraction and the share of clock-gene counts it explains."""
    comps = [f"{organ}:{c}" for c in ["parenchyma"] + OTHERS]
    t = table("class", comps, GENES + MARKERS[organ])
    full = t.groupby("donor").compartment.nunique()
    t = t[t.donor.isin(full.index[full == len(comps)])]
    par = t[t.compartment == f"{organ}:parenchyma"].set_index("donor")
    rows, corrected = [], [par.reset_index()]
    for c in OTHERS:
        x = t[t.compartment == f"{organ}:{c}"].set_index("donor").loc[par.index]
        a = (x[MARKERS[organ]].sum(1) / x.depth) / (par[MARKERS[organ]].sum(1) / par.depth)
        exp = a.to_numpy()[:, None] * x.depth.to_numpy()[:, None] * (par[GENES] / par.depth.to_numpy()[:, None])
        share = exp.sum(0) / x[GENES].sum(0)
        rows.append(dict(organ=organ, compartment=c, n_fetuses=len(x), a_median=a.median(), a_min=a.min(),
                         a_max=a.max(), clock_share_median=share.median(), clock_share_max=share.max(),
                         clock_share_max_gene=share.idxmax()))
        y = x.copy()
        y[GENES] = np.round(np.clip(x[GENES].to_numpy() - exp.to_numpy(), 0, None))
        corrected.append(y.reset_index())
    return pd.DataFrame(rows), pd.concat(corrected, ignore_index=True)


def receptors():
    rows = []
    for level in ("organ", "class"):
        g = genome(level)
        X = g[:, RECEPTORS].X.toarray() / g.obs.depth.to_numpy()[:, None] * 1e6
        d = pd.DataFrame(X, columns=RECEPTORS).assign(compartment=g.obs.compartment.to_numpy(), level=level)
        rows.append(d.groupby(["level", "compartment"]).agg(n_fetuses=("MC2R", "size"),
                                                          **{r: (r, "median") for r in RECEPTORS}).reset_index())
    return pd.concat(rows, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-perm", type=int, default=10000)
    a = ap.parse_args()
    limit_threads(4)
    np.seterr(divide="ignore", invalid="ignore")
    OUT.mkdir(exist_ok=True)
    pd.set_option("display.width", 250)
    tabs = {o: contamination(o) for o in MARKERS}
    cont = pd.concat([v[0] for v in tabs.values()], ignore_index=True)
    cont.to_csv(OUT / "posthoc_contamination.csv", index=False)
    print(cont.round(4).to_string(index=False), flush=True)
    corr = tabs["Adrenal"][1][["donor", "compartment", "n_cells", "depth", "stratum"] + GENES]
    g = gate(corr, WF, a.n_perm, SEED)
    res = pd.DataFrame([dict(test="H190.1 gate, contamination-corrected (post hoc)", n_fetuses=corr.donor.nunique(),
                             p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share)])
    res.to_csv(OUT / "posthoc_gate_corrected.csv", index=False)
    print(res.round(4).to_string(index=False), flush=True)
    rec = receptors()
    rec.to_csv(OUT / "posthoc_receptors.csv", index=False)
    print(rec[rec.level == "organ"].round(2).to_string(index=False))
    print(rec[rec.compartment.str.match("(Adrenal|Heart|Kidney|Lung|Liver):")].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
