"""Slide figure (results, before birth, 2.1): before birth, some organs carry a clock-shaped state shared between
individuals; others, at adult power, do not.

Clock-gene covariation between individuals (C, signed coherence of the 14 clock genes against the adult waveform,
cohort / batch and age adjusted) per organ:
- human fetus (blue): adrenal cortex cells (A244, single-cell pseudobulks, 34 fetuses); whole adrenal, blind replication
  (A328 part a, del Valle 2023, 32 fetuses, 7-11.5 wpc, stage and sex removed; H328.1); lung, liver, heart, kidney (A243,
  pooled open cohorts, 4-20 wpc); adrenal medulla cells (A244, 30 fetuses); brain (A239, arrays, 20 fetuses);
- late mouse embryo (orange; A309, Qiu 2024, E15-E18.75): adrenocortical cells, hepatocytes, lung progenitors,
  ventricular cardiomyocytes, nephron progenitors.
Organs along x; bars: measured C, filled when significant and above the 95th percentile of 200 expression-matched and
200 stress-responsive gene sets, open when not. Open diamond: median C of adult (brain: postnatal) tissue subsampled to the same design with the full-size clock
(A243, A244, A309 power tables; A328 power stage; A240 array curve). Grey line across each bar: 95th percentile of the expression-matched
sets.

Status: A243, A244, A239 registered (blind); A309 registered, NOT BLIND for the earlier fetal results; A328 part a
registered, BLIND. Aggregates only.
Data: exploratory/A243_fetal_organs_pooled/output/{l1,verdicts,position_in_adult}.csv;
exploratory/A244_fetal_adrenal_classes/output/{l1,verdicts_a244,position}.csv;
exploratory/A328_fetal_running_loop/output/{score_a/l1,score_a/verdicts,power_a/passrates}.csv;
exploratory/A239_fetal_brain/output/{l1,verdicts}.csv; exploratory/A240_fetal_amplitude/output/curve_array.csv;
exploratory/A309_mouse_embryo_organs/output/{l1,power_passrates}.csv.
Output: presentation/figures/development/development.{svg,png} (half).

    python scripts/development.py
"""
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "development"
TAB = sns.color_palette().as_hex()
SPECIES = {"human": TAB[0], "mouse": TAB[1]}  # human fetus, mouse embryo
MATCHED = "#BBBBBB"
ORGANS = ["adrenal cortex", "whole adrenal, blind", "liver", "lung", "heart", "kidney", "adrenal medulla", "brain"]
DX, BAR = {"human": -0.19, "mouse": 0.19}, 0.33
LABEL = {"adrenal cortex": "adrenal\ncortex", "whole adrenal, blind": "whole\nadrenal,\nblind",
         "adrenal medulla": "adrenal\nmedulla"}


def covariation():
    """One row per organ x species: C, pass, matched q95, adult median C."""
    a243 = EXP / "A243_fetal_organs_pooled/output"
    l1 = pd.read_csv(a243 / "l1.csv").set_index("test")
    v = pd.read_csv(a243 / "verdicts.csv").set_index("test")
    adult = pd.read_csv(a243 / "position_in_adult.csv").query("setting == 'a1'").set_index("test")
    rows = [dict(organ=o, species="human", C=l1.loc[f"{o}:all", "C"], ok=v.loc[f"{o}:all", "verdict"] == "PASS",
                 q95=l1.loc[f"{o}:all", "matched_q95"], adult=adult.loc[f"{o}:all", "adult_median_C"])
            for o in ("liver", "lung", "heart", "kidney")]
    a244 = EXP / "A244_fetal_adrenal_classes/output"
    c = pd.read_csv(a244 / "l1.csv").set_index("test").loc["cortex"]
    rows.append(dict(organ="adrenal cortex", species="human", C=c.C, q95=c.matched_q95,
                     ok=pd.read_csv(a244 / "verdicts_a244.csv").set_index("test").loc["cortex", "verdict"] == "PASS",
                     adult=pd.read_csv(a244 / "position.csv").set_index("test").loc["cortex", "adult_median"]))
    med = pd.read_csv(a244 / "verdicts_a244.csv").set_index("test").loc["medulla"]
    rows.append(dict(organ="adrenal medulla", species="human", C=med.C, ok=med.verdict.upper() == "PASS",
                     q95=pd.read_csv(a244 / "l1.csv").set_index("test").loc["medulla", "matched_q95"],
                     adult=pd.read_csv(a244 / "position.csv").set_index("test").loc["medulla", "adult_median"]))
    a328 = EXP / "A328_fetal_running_loop/output"
    dv = pd.read_csv(a328 / "score_a/l1.csv").set_index("test").loc["dv:all"]
    rows.append(dict(organ="whole adrenal, blind", species="human", C=dv.C, q95=dv.matched_q95,
                     ok=pd.read_csv(a328 / "score_a/verdicts.csv").set_index("test").loc["dv:all", "verdict"] == "PASS",
                     adult=pd.read_csv(a328 / "power_a/passrates.csv").query("test == 'dv:all' and setting == 'a1'")
                     .median_C.iloc[0]))
    br = pd.read_csv(EXP / "A239_fetal_brain/output/l1.csv").set_index("test").loc["brain_mean:prenatal"]
    curve = pd.read_csv(EXP / "A240_fetal_amplitude/output/curve_array.csv")
    rows.append(dict(organ="brain", species="human", C=br.C, q95=br.matched_q95,
                     ok=pd.read_csv(EXP / "A239_fetal_brain/output/verdicts.csv").set_index("test")
                     .loc["brain_mean:prenatal", "verdict"] == "PASS",
                     adult=curve.query("label == 'brain_mean:postnatal@20' and setting == 1.0").median_stat.iloc[0]))
    a309 = EXP / "A309_mouse_embryo_organs/output"
    m = pd.read_csv(a309 / "l1.csv").set_index("readout")
    pw = pd.read_csv(a309 / "power_passrates.csv").query("setting == 'a1'").set_index("readout")
    for r, o in (("adrenal", "adrenal cortex"), ("liver", "liver"), ("lung", "lung"), ("heart", "heart"),
                 ("kidney", "kidney")):
        rows.append(dict(organ=o, species="mouse", C=m.loc[r, "C"], ok=m.loc[r, "verdict"] == "PASS",
                         q95=m.loc[r, "matched_q95"], adult=pw.loc[r, "median_C"]))
    return pd.DataFrame(rows)


def panel_a(ax, d):
    """Organs along x; per organ one bar per species, at fixed slots (mouse slot empty when not measured) (measured C; filled = clock-specific, open = not), the adult-size
    expectation (open diamond) and the matched-genes q95 (grey tick) at the bar's position."""
    for _, r in d.iterrows():
        x, w = ORGANS.index(r.organ) + DX[r.species], BAR
        c = SPECIES[r.species]
        ax.bar([x], [r.C], width=w, color=c if r.ok else "white", edgecolor=c, lw=1.4, zorder=2)
        ax.plot([x - w / 2, x + w / 2], [r.q95] * 2, color=MATCHED, lw=2.2, solid_capstyle="butt", zorder=3)
        ax.scatter([x], [r.adult], s=34, marker="D", facecolor="white", edgecolor=c, lw=1.1, zorder=4)
    ax.axhline(0, color=pf.GREY, lw=0.6, zorder=1)
    ax.set_xticks(range(len(ORGANS)), [LABEL.get(o, o) for o in ORGANS])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(ORGANS) - 0.4)
    ax.set_ylim(-0.6, 1.0)
    ax.set_ylabel("clock-gene covariation\nbetween individuals, C")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("A", kind="half", height=125)
    d = covariation()
    panel_a(lay.axes["A"], d)
    h = [Patch(facecolor=SPECIES["human"], edgecolor=SPECIES["human"], label="human fetus"),
         Patch(facecolor=SPECIES["mouse"], edgecolor=SPECIES["mouse"], label="mouse embryo"),
         Patch(facecolor="white", edgecolor="#222222", label="not clock-specific"),
         Line2D([], [], ls="", marker="D", ms=7, mfc="white", mec="#222222", label="expected at adult size"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct")]
    lay.shared_legend("A", h, [x.get_label() for x in h], side="bottom", ncol=3, fontsize=10,
                      handletextpad=0.3, columnspacing=1.0)
    fig = lay.finalize()
    print(d.round(3).to_string(index=False))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "development")


if __name__ == "__main__":
    main()
