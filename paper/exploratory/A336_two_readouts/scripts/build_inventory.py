"""A336 inventory: the two subagent inventories (Figures 1, 2, 4 | 3, 5, 6) merged, with the run tier and planned
readouts fixed at registration. Read-only (no statistic computed).

    python build_inventory.py
Output: ../output/inventory.csv
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]

# tier 1: run in this amendment (in the registered order); 2: registered, run after tier 1 under the same rules;
# 3: not run (reason). gene_drop: "364" (all), "100" (100 random drops, seed 336, trio added), "none" (with reason).
PLAN = {
    # Figure 1b (first)
    "A252_H252.1": (1, 1, "364", ""), "A252_H252.2": (1, 1, "364", ""), "A252_H252.3_adipose": (1, 1, "364", ""),
    "A252_H252.4": (3, 0, "none", "already a gene drop (PER1, NR1D1); superseded by the core readout of H252.1-2"),
    "A256_H256.2": (1, 1, "364", ""), "A256_H256.3": (1, 1, "364", ""),
    "A256_H256.1": (3, 0, "none", "not scored in A256 (13 units)"),
    # Figures 4, 5, 6 (frame)
    "A293_H293.1": (1, 2, "364", ""), "A293_H293.2": (1, 2, "364", ""),
    "A295_H295.1": (1, 2, "364", ""), "A295_H295.2": (1, 2, "364", ""),
    "A291_H291.1": (1, 2, "364", ""),
    "A294b_H294b.1": (1, 2, "364", "core exists (A329); drops clock-only"), "A294b_H294b.2": (1, 2, "364", "as H294b.1"),
    "A312_H312.1": (1, 2, "364", "no companion possible (A312)"),
    "A331_H331.1": (1, 2, "364", "core exists (A335 1c)"),
    "A315_H315.1": (1, 2, "364", "core exists (A335 1d)"), "A315_ICU": (1, 2, "364", "core exists (A335 1d)"),
    "A317_H317.1": (1, 2, "364", "part 4"), "A318_H318.1": (1, 2, "364", "part 4"),
    "A304_H304.1_gain": (1, 2, "364", ""), "A304_H304.1_departure": (1, 2, "364", ""),
    "A304_H304.3-4": (1, 2, "364", ""),
    "A306_gain_departure": (1, 2, "100", "power stage 37 min not rerun; clock statistic only for drops"),
    "A297_H297.1": (1, 2, "364", "mouse: Per1, Per2, Nr1d1 dropped"),
    "A322_H322.1": (1, 2, "364", "mouse"), "A324_H324.1": (1, 2, "364", "mouse"), "A324_H324.2": (1, 2, "364", "mouse"),
    "A164_H164.1": (1, 2, "100", "6,307 s run; drops: clock statistic only"), "A164_H164.2": (1, 2, "100", "as H164.1"),
    "A235_H235.1": (1, 2, "100", "drops: clock statistic only"),
    "A233_H233.c.1": (2, 2, "100", "1.8 h score + 3 h power per waveform"), "A233_H233.c.2-3": (2, 2, "100", "as H233.c.1"),
    "A138_H138.4": (1, 2, "100", "1,603 s run"), "A138_H138.8": (1, 2, "100", ""),
    "A138_H138.1": (3, 0, "none", "threshold calibrated by a 5.6 h 14-gene simulation; recalibration out of budget"),
    "A227_H227.1": (2, 2, "100", "about 45 min per family at 2 workers (memory)"), "A227_H227.2": (2, 2, "100", "as H227.1"),
    "A271_H271.1": (2, 2, "100", "as A227"), "A271_H271.2": (1, 2, "100", "2,793 s"),
    "A254_H254.1": (1, 2, "364", ""), "A254_H254.3": (1, 2, "364", ""),
    "A183_H183.1": (1, 2, "364", "core exists (A335 1b)"), "A183_H183.2": (1, 2, "364", ""),
    "A308_H308": (1, 2, "364", "core exists (A335 1b)"), "A115_D115.1": (1, 2, "364", ""), "A332_D1": (1, 2, "364", ""),
    # Figures 1c, 1d
    "A270_H270.1": (1, 3, "364", "baboon: timed reference"), "A270_H270.2": (1, 3, "364", ""), "A270_clockR": (1, 3, "364", ""),
    "A257_H257.1": (1, 3, "364", "core exists (A335 1a)"), "A257_H257.2": (1, 3, "364", "core exists; trio-alone shift blind"),
    "A257_H257.3": (1, 3, "364", "core exists; swap-test tie fixed"),
    "A261_H261.3": (1, 3, "none", "waveform enters only as constants: core constants (4.06 h, 0.239) swapped in; drops not meaningful"),
    "A261_H261.1": (3, 0, "none", "control, waveform-free"),
    # Figure 2 (tier 2: calibrated cell-type tests, placebo nulls redrawn per projector)
    **{k: (2, 4, "100", "calibrated placebo null redrawn per projector; 17-30 min per waveform") for k in (
        "A272_H272.1", "A272_H272.2", "A274_H274.1", "A274_H274.2", "A274_H274.3", "A283_H283.1", "A283_H283.2",
        "A283_H283.3", "A283_H283.4", "A283_H283.5", "A283_V283", "A280_H280.1", "A280_H280.2", "A282_living_D",
        "A282_living_T_SV", "A282_postmortem_D", "A288_H288.1", "A288_H288.2", "A288_H288.3", "A326_H326.1",
        "A279_V279.1", "A279_H279.1", "A279_Bmal1KO", "A273_V273.1", "A273_V273.2")},
    "A275_classes": (1, 4, "364", "per-gene statistic: core = median over the 11 core genes of the saved per-gene table"),
    "A321_H321.1": (3, 0, "none", "already a core-loop test (C11 beside C14); reported as is"),
    # Figure 3 (only existing core values)
    **{k: (3, 0, "none", "Figure 3: existing core / drop value reported, no rerun") for k in (
        "A264_H264.1", "A264_H264.2", "A264_H264.3", "A266_H266.1", "A269_H269.1_lung", "A269_H269.1_liver",
        "A269_H269.2")},
    # not run
    "A293_H293.3": (2, 4, "none", "chain: reads the A272 / A274 / A273 / A279 pair tables (tier 2 reruns first)"),
    "A293_H293.4": (3, 0, "none", "chain: refit of about 34 atlas datasets with frozen anchors (hours, unknown)"),
    "A323_H323.1": (3, 0, "none", "chain on A293 H293.4 and the A104 S table"),
    "A311_H311.1": (3, 0, "none", "the reported primary share (tissues_primary.csv) is written by no script in the repo"),
    "A311_H311.2": (3, 0, "none", "as H311.1"),
    "A313_H313.1": (3, 0, "none", "about 7 h per waveform at 4 workers"), "A313_H313.2": (3, 0, "none", "as H313.1"),
    "A319_H319.1": (3, 0, "none", "published explant phases, no genes"), "A319_H319.2": (3, 0, "none", "as H319.1"),
    "A265_H265.1": (3, 0, "none", "COFE is genome-wide; only the clocksync comparator depends on the clock genes; 0.75 h run"),
    "A265_H265.2": (3, 0, "none", "6.8 h simulation"), "A267_H267.1": (3, 0, "none", "6.0 h run"),
    "A267_H267.2": (3, 0, "none", "waveform-free (genome-wide COFE)"),
}


def main():
    a = pd.read_csv(ROOT / "inventory_part_a.csv")
    b = pd.read_csv(ROOT / "inventory_part_b.csv")
    a["registered_MDE"] = ""
    inv = pd.concat([a, b[a.columns]], ignore_index=True).drop_duplicates("id", keep="last")
    miss = sorted(set(inv.id) - set(PLAN))
    if miss:
        raise SystemExit(f"no plan for: {miss}")
    p = pd.DataFrame.from_dict(PLAN, orient="index", columns=["tier", "order_block", "gene_drop", "plan_note"])
    inv = inv.join(p, on="id")
    inv["readouts"] = inv.tier.map({1: "14 genes | core | trio contribution | gene-drop | trio alone",
                                    2: "14 genes | core | trio contribution | gene-drop | trio alone (after tier 1)",
                                    3: "not run"})
    inv = inv.sort_values(["tier", "order_block", "id"])
    (ROOT / "output").mkdir(exist_ok=True)
    inv.to_csv(ROOT / "output" / "inventory.csv", index=False)
    print(inv.groupby("tier").size().to_string())


if __name__ == "__main__":
    main()
