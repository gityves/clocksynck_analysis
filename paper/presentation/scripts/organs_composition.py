"""Slide figure (results, R2.3): back to the organs. A sample's cell composition does not set its organ's phase, and it
explains only part of the organ's gain.

A. Peak time (clocksync offsets negated, as CHIRAL Fig. 1B; user, 2026-11-23). Top row: peak time of each GTEx
   peripheral tissue relative to its person's body time (organ map, A104 table; one dot per
   tissue, coloured by organ system). Below: phase change of a tissue's clock state per SD of a marker score of its
   samples (A295, registered, H295.2 pass): stromal-vascular and lymphocyte scores, one dot per tissue, black: pooled
   value with 90 % interval; adrenal cortex score (A299, descriptive specificity arm), estimate and 90 % interval.
B. Gain, fold on a log scale. Top row: each tissue's gain on the body's time (relative amplitude, A104 table). Below:
   gain change per SD of the marker scores (A295 descriptive for stromal-vascular, H295.1 fail for lymphocytes; A299
   cortex), as fold = 1 + change.
C. Tissues made mostly of one cell class (post hoc choice: squamous epithelium, oesophagus mucosa and skin; adipocytes,
   visceral and subcutaneous adipose and breast; smooth muscle, oesophagus muscularis and sigmoid colon): each tissue's
   gain on the body's time, one labelled dot per tissue; right, highest over lowest gain of the row.
Status: display post hoc; every value is from a scored table. Writes the plotted rows to
presentation/output/organs_composition/rows.csv. GTEx aggregates only (one value per tissue).

Data: presentation/output/organs/couplings.csv; exploratory/A295_composition/output/{tissues,tests}.csv;
exploratory/A299_adrenal_medulla/output/tests.csv. Output: presentation/figures/organs_composition/organs_composition.{svg,png}.

    python scripts/organs_composition.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import paperfig as pf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _systems import COLOR, gtex_system, legend_in  # noqa: E402

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "organs_composition"
TAB_OUT = PAPER / "presentation" / "output" / "organs_composition"
INK, MUTED = "#222222", pf.GREY
TAB = sns.color_palette().as_hex()
SCORE = {"stromal_vascular": ("+1 SD stromal-\nvascular cells", TAB[5]), "lymphocyte": ("+1 SD lymphocytes", TAB[9])}
CORTEX = ("adrenal: +1 SD\ncortex cells", COLOR["endocrine"])
CLASSES = {  # tissues made mostly of one cell class (post hoc): GTEx tissue -> short label
    "squamous\nepithelium": {"Esophagus - Mucosa": "oesophagus\nmucosa", "Skin - Sun Exposed (Lower leg)": "skin,\nsun",
                             "Skin - Not Sun Exposed (Suprapubic)": "skin,\nno sun"},
    "adipocyte": {"Adipose - Visceral (Omentum)": "visceral", "Adipose - Subcutaneous": "subcut.",
                  "Breast - Mammary Tissue": "breast"},
    "smooth muscle": {"Esophagus - Muscularis": "oesophagus\nmuscularis", "Colon - Sigmoid": "sigmoid\ncolon"}}


def rows():
    """Long table of the plotted values: panel, row, label, value, interval, colour key."""
    org = pd.read_csv(PAPER / "presentation/output/organs/couplings.csv")
    org = org[org.group == "periphery"].set_index("tissue")
    tis = pd.read_csv(EXP / "A295_composition/output/tissues.csv")
    tests = pd.read_csv(EXP / "A295_composition/output/tests.csv").set_index("score")
    cx = pd.read_csv(EXP / "A299_adrenal_medulla/output/tests.csv")
    cx = cx[(cx.score == "cortex") & cx.donors.str.startswith("Hardy")].iloc[0]
    R = []
    for t, r in org.iterrows():
        R.append(dict(panel="phase", row="organ", label=t, value=r.offset_h, lo=r.offset_lo, hi=r.offset_hi,
                      key=gtex_system(t)))
        R.append(dict(panel="gain", row="organ", label=t, value=np.exp(r.log_amp_rel), lo=np.exp(r.log_amp_lo),
                      hi=np.exp(r.log_amp_hi), key=gtex_system(t)))
    for cls, ts in CLASSES.items():
        for t in ts:
            R.append(dict(panel="class", row=cls, label=t, value=np.exp(org.log_amp_rel[t]), lo=np.exp(org.log_amp_lo[t]),
                          hi=np.exp(org.log_amp_hi[t]), key=gtex_system(t)))
    for s in SCORE:
        t = tis[tis.score == s]
        for r in t.itertuples():
            R.append(dict(panel="phase", row=s, label=r.tissue, value=r.phase_change_h, lo=np.nan, hi=np.nan, key=s))
            R.append(dict(panel="gain", row=s, label=r.tissue, value=1 + r.gain_change, lo=np.nan, hi=np.nan, key=s))
        p = tests.loc[s]
        R.append(dict(panel="phase", row=s, label="pooled", value=p.pooled_phase_change_h, lo=p.phase_lo90,
                      hi=p.phase_hi90, key=s))
        R.append(dict(panel="gain", row=s, label="pooled", value=1 + p.pooled_gain_change,
                      lo=1 + p.pooled_gain_change - 1.645 * p.gain_se, hi=1 + p.pooled_gain_change + 1.645 * p.gain_se,
                      key=s))
    R.append(dict(panel="phase", row="cortex", label="Adrenal Gland", value=cx.phase_change_h, lo=cx.phase_lo90,
                  hi=cx.phase_hi90, key="cortex"))
    R.append(dict(panel="gain", row="cortex", label="Adrenal Gland", value=1 + cx.gain_change, lo=1 + cx.gain_lo90,
                  hi=1 + cx.gain_hi90, key="cortex"))
    R = pd.DataFrame(R)
    ph = R.panel == "phase"  # peak time (user, 2026-11-23): clocksync offsets (> 0 advanced) negated, as CHIRAL Fig. 1B
    R.loc[ph, ["value", "lo", "hi"]] = -R.loc[ph, ["value", "hi", "lo"]].to_numpy()
    return R


def draw(ax, R, order, labels):
    for i, row in enumerate(order):
        d = R[R.row == row]
        if row == "organ":
            y = i + np.random.default_rng(i).uniform(-0.2, 0.2, len(d))
            ax.scatter(d.value, y, s=22, c=[COLOR[k] for k in d.key], edgecolor="white", lw=0.3, zorder=2)
            continue
        if row == "cortex":
            r = d.iloc[0]
            ax.plot([r.lo, r.hi], [i, i], color=CORTEX[1], lw=1.6, zorder=2)
            ax.scatter([r.value], [i], s=34, color=CORTEX[1], zorder=3)
            continue
        t, p = d[d.label != "pooled"], d[d.label == "pooled"].iloc[0]
        y = i + np.random.default_rng(i).uniform(-0.2, 0.2, len(t))
        ax.scatter(t.value, y, s=16, color=SCORE[row][1], alpha=0.6, edgecolor="none", zorder=1)
        ax.plot([p.lo, p.hi], [i, i], color=INK, lw=1.6, zorder=2)
        ax.scatter([p.value], [i], s=34, color=INK, zorder=3)
    ax.set_yticks(range(len(order)), labels)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(order) - 0.5, -0.5)


def same_class(ax, R):
    """Gain of tissues made mostly of one cell class: one row per class, one labelled dot per tissue."""
    for i, cls in enumerate(CLASSES):
        d = R[R.row == cls].sort_values("value")
        ax.plot([d.value.min(), d.value.max()], [i, i], color=MUTED, lw=0.8, zorder=1)
        ax.scatter(d.value, [i] * len(d), s=30, c=[COLOR[k] for k in d.key], edgecolor="white", lw=0.4, zorder=2)
        for j, r in enumerate(d.itertuples()):
            ax.annotate(CLASSES[cls][r.label], (r.value, i), xytext=(0, 10 if j % 2 == 0 else -14),
                        textcoords="offset points", ha="center", va="bottom" if j % 2 == 0 else "top", fontsize=7.5,
                        color=INK, linespacing=0.9, zorder=3,
                        bbox=dict(boxstyle="square,pad=0.1", facecolor="white", edgecolor="none"))
        ax.annotate(f"×{d.value.max() / d.value.min():.1f}", (d.value.max(), i), xytext=(8, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=8, color=MUTED)
    ax.set_yticks(range(len(CLASSES)), list(CLASSES))
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(CLASSES) - 0.4, -0.6)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    TAB_OUT.mkdir(parents=True, exist_ok=True)
    R = rows()
    R.to_csv(TAB_OUT / "rows.csv", index=False)
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1, 1, 0.85])
    a, b, c = lay.axes["A"], lay.axes["B"], lay.axes["C"]
    comp = [SCORE["stromal_vascular"][0], SCORE["lymphocyte"][0], CORTEX[0]]
    draw(a, R[R.panel == "phase"], ["organ", "stromal_vascular", "lymphocyte", "cortex"], ["organ peak time"] + comp)
    a.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    a.set_xlim(-5.5, 2.5)
    a.set_xlabel("peak time relative to the body [h]")
    ad = R[(R.panel == "phase") & (R.row == "organ") & (R.label == "Adrenal Gland")].iloc[0]
    y_ad = a.collections[0].get_offsets()[list(R[(R.panel == "phase") & (R.row == "organ")].label).index(
        "Adrenal Gland")][1]
    a.annotate("adrenal", (ad.value, y_ad), xytext=(6, 0), textcoords="offset points", ha="left", va="center",
               fontsize=7.5, color=COLOR["endocrine"])
    legend_in(a, fontsize=7.5, loc="lower left")
    draw(b, R[R.panel == "gain"], ["organ", "stromal_vascular", "lymphocyte", "cortex"], ["organ gain"] + comp)
    same_class(c, R[R.panel == "class"])
    for ax, lim in ((b, (0.2, 3.2)), (c, (0.6, 3.2))):
        ax.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
        ax.set_xscale("log")
        ax.set_xlim(*lim)
        ax.minorticks_off()
    b.set_xticks([0.25, 0.5, 1, 2], ["0.25", "0.5", "1", "2"])
    c.set_xticks([0.75, 1, 1.5, 2, 3], ["0.75", "1", "1.5", "2", "3"])
    b.set_xlabel("gain [fold]")
    c.set_xlabel("organ gain [fold]")
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_composition")


if __name__ == "__main__":
    main()
