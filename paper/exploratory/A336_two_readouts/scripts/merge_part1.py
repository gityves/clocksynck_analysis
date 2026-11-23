"""A336: two_readouts.csv = the part-1 rows (`../output/part1/*.csv`), the part-4 rows (A317 H317.1, A318 H318.1 from
`braindeath_core.csv`) and the part-2 row (GTEx person loudness, `person_gtex.csv`, if scored), plus every inventory
statistic not run (tier 2 / 3) with its reason, so that no statistic is dropped.

    python merge_part1.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import a336_common as C  # noqa: E402

OUT = HERE.parents[1] / "output"

# Stored 14-gene verdicts, read from each entry's own test table (filled where checked; empty = not filled here). A
# core "fall" on a statistic that already fails on 14 genes is no change of verdict.
VERDICT_14 = {
    "A227_H227.1": "PASS", "A227_H227.2": "PASS", "A271_H271.1": "FAIL", "A271_H271.1:myeloid": "FAIL",
    "A271_H271.1:T_NK": "FAIL", "A271_H271.1:CD8neg_T": "FAIL", "A272_H272.1": "FAIL (above matched, not responsive q95)",
    "A272_H272.2": "PASS", "A273_V273.1": "PASS", "A273_V273.2": "PASS", "A274_H274.1": "PASS", "A274_H274.2": "PASS",
    "A274_H274.3": "PASS", "A279_V279.1": "PASS", "A279_H279.1": "FAIL", "A279_Bmal1KO": "descriptive (conditions met)",
    "A280_H280.1": "PASS", "A280_H280.2": "FAIL", "A282_living_D": "descriptive (both p_cal < 0.05)",
    "A282_living_T_SV": "descriptive (matched p_cal < 0.05 only)", "A282_postmortem_D": "descriptive (both p_cal < 0.05)",
    "A283_H283.1": "FAIL (matched p_cal 0.010, responsive 0.299)", "A283_H283.2": "PASS", "A283_H283.3": "PASS",
    "A283_H283.4": "PASS", "A283_H283.5": "PASS", "A283_V283": "PASS", "A293_H293.3": "PASS (equivalence, A335)",
    "A326_H326.1:A281": "FAIL (power < 0.8 at r = 0.25)", "A326_H326.1:A287": "FAIL (power < 0.8 at r = 0.25)",
    "A326_H326.1:A321": "FAIL (power < 0.8 at r = 0.25)", "A288_H288.1": "FAIL (responsive p_cal 0.133)",
    "A288_H288.2": "PASS", "A288_H288.3": "FAIL", "A295_H295.1": "FAIL", "A324_H324.2": "FAIL",
    **{f"A233_H233.c.1:{c}": "PASS (equivalent)" for c in ("htn", "ihd", "copd", "depression", "smoking")},
    **{f"A233_H233.c.2-3:{c}:{k}": "null kept" for c in ("htn", "ihd", "copd", "depression", "smoking")
       for k in ("Q", "module_d_star", "module_offset")},
}


def part4():
    b = pd.read_csv(OUT / "braindeath_core.csv")
    rows = []
    for adj, ident, stored in (("unadjusted (H317.1)", "A317_H317.1", 0.620925338252315),
                               ("ischaemic time and RIN adjusted (H318.1)", "A318_H318.1", 0.5841177196571063)):
        f, c = (b[(b.adjustment == adj) & (b.gene_set == g)].iloc[0] for g in ("14 genes", "core"))
        first = b.iloc[0]
        r = dict(id=ident, figure_panel="Fig 5b", statistic=f"brain-death order D_rho, {adj}", direction=1,
                 value_14_stored=stored, value_14=f.D_rho, reproduced=bool(np.isclose(f.D_rho, stored, atol=1e-9)),
                 value_core=c.D_rho, core_criterion="D_rho > 0, relabelling p < 0.05, > core matched and responsive q95",
                 core_status="holds" if c.passes else "falls", n_drops=c.n_drops, drop_pct_all=c.pct_all,
                 n_similar=c.n_similar, drop_pct_similar=c.pct_similar, drop_q10=c.drop_q10, drop_median=c.drop_median,
                 trio_drop_below_q10=c.below_q10, reading=c.reading, blind="core blind",
                 notes=f"core p {c.p_labels:.4f}, matched q95 {c.matched_q95:.3f}, responsive q95 {c.responsive_q95:.3f}",
                 script="exploratory/A336_two_readouts/scripts/a336_braindeath.py")
        if ident == "A317_H317.1":
            r.update(diff_14_minus_core=first.diff_14_minus_core, diff_lo90=first.diff_lo90, diff_hi90=first.diff_hi90,
                     diff_draws="200 shared donor bootstraps (seed [336, 4])", trio_gate_passes=bool(first.trio_gate_passes),
                     trio_alone_value=first.trio_alone_D_rho)
        rows.append(r)
    return rows


def part2():
    f = OUT / "person_gtex.csv"
    if not f.exists():
        return []
    p = pd.read_csv(f).set_index("gene_set")
    f14, c = p.loc["14 genes"], p.loc["core"]
    return [dict(id="A336_P2_gtex_person", figure_panel="Fig 4f", statistic="GTEx person leverage (half-host radii)",
                 direction=1, value_14_stored=np.nan, value_14=f14.leverage, reproduced=True, value_core=c.leverage,
                 core_criterion="core leverage > core matched q95 and persons-permuted p < 0.05",
                 core_status="holds" if c.holds else "falls", diff_14_minus_core=f14.diff_14_minus_core,
                 diff_lo90=f14.diff_lo90, diff_hi90=f14.diff_hi90, diff_draws="1,000 shared donor bootstraps",
                 n_drops=c.n_drops, drop_pct_all=c.drop_pct_all, n_similar=c.n_similar, drop_pct_similar=c.drop_pct_similar,
                 drop_q10=c.drop_q10, drop_median=c.drop_median, trio_drop_below_q10=c.below_q10,
                 trio_gate_passes=f14.trio_gate_passes, trio_alone_value=f14.trio_alone_leverage, reading=c.reading,
                 blind="14 genes not blind (A331 power stage); core blind",
                 notes="14-gene reproduction: A331 power-stage cv2 median", script="exploratory/A336_two_readouts/scripts/a336_person_gtex.py")]


def main():
    parts = [pd.read_csv(f) for f in sorted((OUT / "part1").glob("*.csv"))]
    parts = [p for p in parts if list(p.columns) == C.COLUMNS]  # row files only (not replicate / pair tables)
    extra = pd.DataFrame(part4() + part2())
    tab = pd.concat(parts + [extra], ignore_index=True)
    inv = pd.read_csv(OUT / "inventory.csv")
    done = set(tab.id)
    alias = {"A304_H304.3-4": "A304_H304.3", "A306_gain_departure": "A306_", "A270_clockR": "A270_clockR",
             "A275_classes": "A275_classes", "A235_H235.1": "A235_H235.1_", "A233_H233.c.1": "A233_H233.c.1:",
             "A233_H233.c.2-3": "A233_H233.c.2-3:", "A271_H271.1": "A271_H271.1",
             "A326_H326.1": "A326_H326.1:"}
    ran = lambda i: i in done or any(str(d).startswith(alias.get(i, "\0")) for d in done)  # noqa: E731
    rest = inv[~inv.id.map(ran)]
    nr = pd.DataFrame(dict(id=rest.id, figure_panel=rest.figure_panel, statistic=rest.statistic,
                           reading=np.where(rest.tier == 3, "not run (tier 3)", "not run yet (tier " + rest.tier.astype(str) + ")"),
                           notes=rest.plan_note, script=rest.script))
    tab = pd.concat([tab, nr], ignore_index=True)
    for col in C.COLUMNS:
        if col not in tab:
            tab[col] = np.nan
    tab["verdict_14"] = tab.id.map(VERDICT_14)
    tab[C.COLUMNS + ["verdict_14"]].to_csv(OUT / "two_readouts.csv", index=False)
    print(tab.reading.value_counts().to_string())


if __name__ == "__main__":
    main()
