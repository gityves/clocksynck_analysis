"""A336 part 1, A275 (Figure 2b; per-gene sharing, descriptive): clock-class sharing on 14 genes, on the core loop and on
every three-gene drop, from the saved per-gene table (`A275.../output/genes.csv`; class value = median over genes per
atlas, mean over atlases). The registered reading rule (immediate-early vs glucocorticoid vs clock) is re-applied with
the core clock class. Trio contribution on 200 atlas bootstraps (rng [336, 1]).

    python p1_a275.py
Output: ../output/part1/a275.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(HERE.parent))
import a336_common as C  # noqa: E402

SRC = PAPER / "exploratory/A275_shared_gene_classes/output"


def klass_value(G, klass, genes, atlases=None):
    g = G[(G.klass == klass) & G.gene.isin(genes)]
    per = g.groupby("atlas").sharing.median()
    if atlases is not None:
        return float(np.mean([per[a] for a in atlases if a in per]))
    return float(per.mean())


def rule(clock, gc, ie):
    if ie >= clock and gc < clock:
        return "handling state"
    if gc >= clock and ie < clock:
        return "hormonal state"
    if ie >= clock and gc >= clock:
        return "both"
    return "neither"


def main():
    G = pd.read_csv(SRC / "genes.csv")
    S = pd.read_csv(SRC / "classes.csv").set_index(["klass", "gene_class"])
    clock = sorted(G[G.gene_class == "clock"].gene.unique())
    assert len(clock) == 14
    core = [g for g in clock if g not in C.TRIO]
    dl, tl = C.drops()
    wfs = dict(dl)
    rows = []
    rb = np.random.default_rng([336, 1])
    for klass in ("SV_X", "LY_X"):
        v14 = klass_value(G, klass, clock)
        stored = float(S.loc[(klass, "clock"), "sharing"])
        assert np.isclose(v14, stored, atol=1e-9), (klass, v14, stored)
        vc = klass_value(G, klass, core)
        gc, ie, mt = (float(S.loc[(klass, c), "sharing"]) for c in ("gc_target", "immediate_early", "matched"))
        atl = sorted(G[G.klass == klass].atlas.unique())
        bs = []
        for _ in range(200):
            pick = list(rb.choice(atl, len(atl)))
            bs.append(klass_value(G, klass, clock, pick) - klass_value(G, klass, core, pick))
        dv = {lab: klass_value(G, klass, list(w.genes)) for lab, w in dl}
        assert np.isclose(dv[tl], vc)
        pr = C.drop_percentiles(dv, tl, wfs, direction=1)
        r14, rc = rule(v14, gc, ie), rule(vc, gc, ie)
        status = "holds" if r14 == rc else "falls"
        rows.append(dict(id=f"A275_classes_{klass}", figure_panel="Fig 2b", statistic=f"clock-class sharing, {klass} pairs",
                         direction=1, value_14_stored=stored, value_14=v14, reproduced=True, value_core=vc,
                         core_criterion="registered reading rule (immediate-early / glucocorticoid vs clock) unchanged",
                         core_status=status, diff_14_minus_core=v14 - vc, diff_lo90=float(np.quantile(bs, 0.05)),
                         diff_hi90=float(np.quantile(bs, 0.95)), diff_draws="200 atlas bootstraps (rng [336, 1])", **pr,
                         trio_gate_passes="not computable", reading=C.reading(status, pr["trio_drop_below_q10"]),
                         blind="core blind (descriptive entry)",
                         notes=f"rule 14 genes: {r14}; core: {rc}; gc {gc:.3f}, immediate-early {ie:.3f}, matched {mt:.3f}; "
                               f"trio-alone class value {klass_value(G, klass, list(C.TRIO)):.3f} (no gate)",
                         script="exploratory/A336_two_readouts/scripts/p1_a275.py"))
    df = C.write_rows("a275", rows)
    print(df.drop(columns=["script"]).round(4).T.to_string())


if __name__ == "__main__":
    main()
