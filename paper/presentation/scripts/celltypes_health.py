"""Slide figure (results, R2.2): inside tissues, stromal-vascular cells and macrophages share the tissue's clock state;
lymphocytes do not.

Readout (anchor-free, `clocksync.pairwise`): for two cell classes of the same donors, |r| = the complex correlation over
donors of their clock-plane states; a class's sharing = the median |r| over its pairs with the non-immune classes of
the tissue (X). Companions: 200 expression-matched and 200 stress-responsive gene sets per atlas; calibrated p from
5000 placebo draws (A283 form).
A. Per atlas, sharing of stromal-vascular cells (SV_X), macrophages (T_MAC, A288) and lymphocytes (LY_X), connected per
   atlas: the ten healthy body atlases of A272 and A274 (descriptive per atlas; A283 per-atlas table).
B. The registered contrasts, each a pooled difference in sharing: stromal-vascular minus lymphocyte in three
   independent atlas sets (A272 = H283.2, A274 = H283.3, IBDverse = H280.1) and in living donors only (A282), and
   macrophage minus lymphocyte (A288 H288.2); grey: medians of the expression-matched (circle) and stress-responsive
   (square) companion sets. All pass against both families (p_cal <= 0.016).

Data: exploratory/A276_pairwise_calibration/output/a283/{per_atlas,tests}.csv; exploratory/A280_celltype_replication2/
output/calibrated/tests.csv; exploratory/A282_pooled_living/output/tests.csv; exploratory/A288_macrophage_pairwise/output/
{per_atlas,tests}.csv. Output: presentation/figures/celltypes_health/celltypes_health.{svg,png}.

    python scripts/celltypes_health.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "celltypes_health"
INK, MUTED = "#222222", pf.GREY
TAB = sns.color_palette().as_hex()
CELL = {"stromal-vascular": TAB[5], "macrophage": TAB[1], "lymphocyte": TAB[9]}
ATLAS = {"lung_hlca": "lung", "heart": "heart", "organ_breast": "breast", "organ_colon": "colon (Pelka)",
         "organ_muscle": "muscle", "kidney_sc": "kidney", "kong_ileum": "ileum (Kong)", "kong_colon": "colon (Kong)", "organ_eye_anterior": "anterior eye",
         "organ_optic_nerve": "optic nerve"}


def per_atlas(ax, pa, mac):
    pa = pa[pa.atlas.isin(ATLAS)]
    sv = pa[pa.klass == "SV_X"].set_index("atlas").clock
    ly = pa[pa.klass == "LY_X"].set_index("atlas").clock
    mc = mac[mac.statistic == "T_MAC"].set_index("atlas").clock
    order = sv.sort_values(ascending=False).index
    for i, a in enumerate(order):
        xs = [v for v in (sv.get(a), mc.get(a), ly.get(a)) if v is not None and np.isfinite(v)]
        ax.plot(xs, [i] * len(xs), color=MUTED, lw=0.8, zorder=1)
        for s_, c in zip((sv, mc, ly), CELL.values()):
            if a in s_.index:
                ax.scatter([s_[a]], [i], s=34, color=c, edgecolor="white", lw=0.5, zorder=2)
    ax.set_yticks(range(len(order)), [ATLAS[a] for a in order])
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.8)
    ax.set_xlabel("clock state shared with\nthe tissue's other cells")
    h = [Line2D([], [], ls="", marker="o", ms=6, color=c, label=lab) for lab, c in CELL.items()]
    ax.legend(handles=h, loc="lower right", frameon=False, fontsize=9, handletextpad=0.1)


def contrasts(ax, rows):
    for i, r in enumerate(rows):
        c = CELL["macrophage"] if r["kind"] == "MAC" else CELL["stromal-vascular"]
        ax.scatter([r["clock"]], [i], s=50, color=c, zorder=3)
        ax.scatter([r["matched"]], [i], s=22, marker="o", facecolor="white", edgecolor=MUTED, lw=1.2, zorder=2)
        ax.scatter([r["responsive"]], [i], s=22, marker="s", facecolor="white", edgecolor=MUTED, lw=1.2, zorder=2)
        ax.plot([min(r["matched"], r["responsive"]), r["clock"]], [i, i], color=MUTED, lw=0.8, zorder=1)
    ax.set_yticks(range(len(rows)), [r["label"] for r in rows])
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_xlim(-0.1, 0.5)
    ax.set_xlabel("difference in shared clock state\nvs lymphocytes")
    h = [Line2D([], [], ls="", marker="o", ms=5, mfc="white", mec=MUTED, label="expression-matched genes"),
         Line2D([], [], ls="", marker="s", ms=5, mfc="white", mec=MUTED, label="stress-responsive genes")]
    ax.legend(handles=h, loc="lower right", frameon=False, fontsize=8, handletextpad=0.1)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    pa = pd.read_csv(EXP / "A276_pairwise_calibration/output/a283/per_atlas.csv")
    t283 = pd.read_csv(EXP / "A276_pairwise_calibration/output/a283/tests.csv").set_index("test")
    t280 = pd.read_csv(EXP / "A280_celltype_replication2/output/calibrated/tests.csv").set_index("test")
    t282 = pd.read_csv(EXP / "A282_pooled_living/output/tests.csv").query("arm == 'living' and statistic == 'D'").iloc[0]
    mac = pd.read_csv(EXP / "A288_macrophage_pairwise/output/per_atlas.csv")
    t288 = pd.read_csv(EXP / "A288_macrophage_pairwise/output/tests.csv").set_index("test")

    def row(r, label, kind="SV"):
        return dict(label=label, kind=kind, clock=r.clock, matched=r.matched_median, responsive=r.responsive_median,
                    p=(r.matched_p_cal, r.responsive_p_cal))
    rows = [row(t283.loc["H283.2"], "stromal-vascular,\nfive atlases"),
            row(t283.loc["H283.3"], "stromal-vascular,\nreplication"),
            row(t280.loc["H280.1 D"], "stromal-vascular,\nIBD cohort"),
            row(t282, "stromal-vascular,\nliving donors"),
            row(t288.loc["H288.2 D_MAC"], "macrophage", "MAC")]
    lay = pf.Layout.slide("AB", kind="wide", col_ratios=[1, 1])
    per_atlas(lay.axes["A"], pa, mac)
    contrasts(lay.axes["B"], rows)
    fig = lay.finalize()
    for r in rows:
        print(r["label"].replace("\n", " "), round(r["clock"], 3), [round(p, 4) for p in r["p"]])
    print(pf.collisions(fig))
    pf.save(fig, OUT / "celltypes_health")


if __name__ == "__main__":
    main()
