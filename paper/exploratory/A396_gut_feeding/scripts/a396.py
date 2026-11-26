"""A396 (descriptive; registered 2026-11-26): gut loudness and the species' feeding pattern; the intensive-care proxy.

Part a: loudness rank (log relative amplitude, descending) of the mouse gut tissues (Stomach, Duodenum, Jejunum, Ileum)
among the 16 peripheral tissues of the A394 clock fit; the same summary for the gut organs of human, baboon and sheep
(A384 gains, ranks among each species' mapped organs); per feeding arm, gut minus other-tissue mean log gain (one fit per
arm, A394 setup). Part b: from A349 (`part2_tissues.csv`, pooled design), digestive and liver tissues against the rest
after the loudness trend, against the A297 arrhythmic-feeding benchmark, if A349 covers them.

    python a396.py
Outputs (`../output/`): part_a.csv, part_a_arms.csv, part_b.csv. Aggregates only.
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory" / "A394_mouse_organisation" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a394 as M  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
EXP = PAPER / "exploratory"
GUT = {"mouse": ["Stomach", "Duodenum", "Jejunum", "Ileum"],
       "human": ["Stomach", "Colon - Transverse", "Colon - Sigmoid", "Esophagus - Muscularis"],
       "baboon": ["STF", "ASC", "DEC", "DUO", "ILE", "CEC", "OES"],
       "sheep": ["Abomasum", "Colon", "Duodenum", "Ileum", "Cecum"]}
DIGESTIVE_B = ["Stomach", "Colon - Transverse", "Colon - Sigmoid", "Esophagus - Gastroesophageal Junction",
               "Esophagus - Muscularis", "Esophagus - Mucosa", "Small Intestine - Terminal Ileum"]
BENCH = (0.65, 0.47, 0.84)  # A297 Greenwell: -2 ln 0.724 (0.657 to 0.792)


def summary(species, gains, gut):
    r = gains.rank(ascending=False)
    g = [x for x in gut if x in r.index]
    return dict(species=species, organs=int(gains.notna().sum()), gut=";".join(g), gut_ranks=";".join(f"{r[x]:.0f}" for x in g),
                gut_median_rank=float(r[g].median()), gut_mean_percentile=float(100 * (1 - (r[g] - 1) / (len(r) - 1)).mean()),
                gut_minus_other_log_gain=float(gains[g].mean() - gains.drop(g).mean()))


def part_a():
    tab = pd.read_csv(M.OUT / "mouse_fit.csv", index_col="tissue")
    pb, per, brain = M.mouse_table()
    rows = [summary("mouse", tab.log_amp_rel.reindex(per).dropna(), GUT["mouse"])]
    gt = pd.read_csv(EXP / "A383_species_organisation" / "output" / "gains_a384.csv")
    for sp in ("human", "baboon", "sheep"):
        g = gt.dropna(subset=[sp]).set_index(sp)[f"log_gain_{sp}"].dropna()
        rows.append(summary(sp, g, GUT[sp]))
    arms = []
    for arm in sorted(pb.stratum.unique()):
        t = M.fit(pb[pb.stratum == arm], M.WF, per, brain, 0, name=f"a396 {arm}")
        arms.append(dict(arm=arm, **summary("mouse", t.log_amp_rel.reindex(per).dropna(), GUT["mouse"])))
    return pd.DataFrame(rows), pd.DataFrame(arms)


def part_b():
    t = pd.read_csv(EXP / "A349_icu_limits" / "output" / "part2_tissues.csv")
    t = t[t.design == "pooled"]
    rows = []
    for gs, g in t.groupby("gene_set"):
        dig = [x for x in g.tissue if x in DIGESTIVE_B]
        rows.append(dict(gene_set=gs, tissues=len(g), digestive=";".join(dig), liver="Liver" in set(g.tissue),
                         informative=bool(len(dig) >= 2), benchmark=BENCH[0]))
    return pd.DataFrame(rows)


def main():
    warnings.filterwarnings("ignore")
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    a, arms = part_a()
    b = part_b()
    a.to_csv(OUT / "part_a.csv", index=False)
    arms.to_csv(OUT / "part_a_arms.csv", index=False)
    b.to_csv(OUT / "part_b.csv", index=False)
    pd.set_option("display.width", 250)
    print(a.round(3).to_string(index=False), arms.round(3).to_string(index=False), b.to_string(index=False), sep="\n")


if __name__ == "__main__":
    main()
