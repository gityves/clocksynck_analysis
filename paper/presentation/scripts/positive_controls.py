"""Slide figure (results, R1.5): the readout moves when the clock is broken or shifted (mouse positive controls, A111/A129).

A. Hepatocyte Bmal1 knockout vs AlbCre control (Manella 2021, ad libitum, 24 vs 24 mice; H111.3 registered, pass): change in
   coupling of each tissue at the control's host spread (d*), point and 95 % interval from the permutations.
B. Daytime feeding, two datasets (H111.1 registered): per tissue, offset after the switch minus offset before, centred on
   the circular mean of the dataset's tissues (the shift of each tissue relative to the others; a whole-body shift is not
   identifiable untimed and is removed from both). x: timed (from the recorded ZT), y: untimed (clocksync); dashed:
   identity. Circles: Xin 2021, food at night (NRF, normal for mice) -> in the day (DRF), 4 tissues (pass, Fisher-Lee
   0.988, exact p 0.042, the minimum with 4 tissues). Squares: Manella 2021 AlbCre controls, ad libitum -> daytime
   feeding, 3 tissues (Fisher-Lee 0.986; FAIL by design: exact p 0.167 is the minimum with 3 tissues, a registration flaw
   noticed at scoring).
Colours: organ systems (`_systems.py`): liver digestive, fat adipose, heart cardiovascular, kidney and lung other.

Data: exploratory/A111_feeding_controls/output/a111_h111_3.csv (primary, ad libitum rows), a111_h111_1.csv (primary rows,
xin2021 and manella2021).
Output: presentation/figures/positive_controls/positive_controls.{svg,png} (half-slide figure).

    python scripts/positive_controls.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

import paperfig as pf
from matplotlib.lines import Line2D

from _systems import COLOR

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "exploratory/A111_feeding_controls/output"
OUT = PAPER / "presentation" / "figures" / "positive_controls"
INK, MUTED = "#222222", pf.GREY
NAMES = {"Liver": "liver", "WAT": "white fat", "Lung": "lung", "liver": "liver", "VAT": "visceral fat", "heart": "heart",
         "kidney": "kidney"}
SYSTEM = {"liver": "digestive", "white fat": "adipose, skin, breast", "visceral fat": "adipose, skin, breast",
          "heart": "cardiovascular", "kidney": "other", "lung": "other"}
COL = {k: COLOR[v] for k, v in SYSTEM.items()}
FEEDING = {"xin2021": ("o", ["liver", "VAT", "heart", "kidney"], "night → day feeding (Xin)"),
           "manella2021": ("s", ["Liver", "WAT", "Lung"], "ad libitum → day feeding (Manella)")}
# label offsets [points] per (dataset, tissue), placed by eye to clear neighbouring marks
NUDGE = {("manella2021", "Liver"): (-8, 6, "right", "bottom"), ("manella2021", "WAT"): (8, -4, "left", "top"),
         ("manella2021", "Lung"): (-6, 6, "right", "bottom")}


def knockout(ax, t):
    t = t.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(t))
    ax.axvline(0, color=MUTED, lw=1, ls="--", zorder=1)
    for yi, (_, r) in zip(y, t.iterrows()):
        c = COL[NAMES[r.target]]
        ax.plot([r.d_star_lo, r.d_star_hi], [yi, yi], color=c, lw=1.8, solid_capstyle="butt", zorder=2)
        ax.scatter([r.d_star], [yi], s=64, color=c, zorder=3)
    ax.set_yticks(y, [NAMES[k] for k in t.target])
    ax.set_ylim(-0.6, len(t) - 0.4)
    ax.set_xlim(-0.5, 0.2)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("change in coupling,\nknockout - control")
    ax.set_title("liver Bmal1 knockout", fontsize=13, color=INK)


def feeding(ax, f):
    ax.plot([-9, 11], [-9, 11], color=MUTED, lw=1, ls="--", zorder=1)
    for ds, (m, tis, _) in FEEDING.items():
        r = f.loc[ds]
        for k in tis:
            a, b, name = r[f"timed_{k}"], r[f"untimed_{k}"], NAMES[k]
            ax.scatter([a], [b], s=50, marker=m, color=COL[name], edgecolor="white", lw=0.5, zorder=3)
            dx, dy, ha, va = NUDGE.get((ds, k), (7, -2, "left", "top"))
            ax.annotate(name, (a, b), xytext=(dx, dy), textcoords="offset points", fontsize=10, color=COL[name], ha=ha,
                        va=va)
    ax.set_xlim(-9, 11)
    ax.set_ylim(-9, 11)
    ax.set_xticks([-6, -3, 0, 3, 6, 9])
    ax.set_yticks([-6, -3, 0, 3, 6, 9])
    ax.set_aspect("equal")
    ax.set_xlabel("phase offset timed [h]")
    ax.set_ylabel("phase offset clocksync [h]")
    ax.set_title("daytime feeding", fontsize=13, color=INK)
    return [Line2D([], [], ls="", marker=m, ms=7, color=INK, label=lab) for m, _, lab in FEEDING.values()]


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    k = pd.read_csv(SRC / "a111_h111_3.csv")
    k = k[(k.dataset == "manella2021") & (k.variant == "primary") & k.contrast.str.startswith("AL:")]
    f = pd.read_csv(SRC / "a111_h111_1.csv").query("variant == 'primary'").set_index("dataset")
    lay = pf.Layout.slide("AB", kind="half", height=110)
    knockout(lay.axes["A"], k)
    h = feeding(lay.axes["B"], f)
    lay.shared_legend("AB", h, [x.get_label() for x in h], side="bottom", ncol=2, fontsize=10)
    fig = lay.finalize()
    print(k[["target", "d_star", "d_star_lo", "d_star_hi"]].round(3).to_string(index=False))
    print(f[["contrast", "n", "r_fl", "p_exact", "verdict"]].to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "positive_controls")


if __name__ == "__main__":
    main()
