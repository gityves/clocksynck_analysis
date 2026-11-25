"""Slide figure (results, drives, 4.2): a 12 h shift of behaviour moves the human muscle clock by only about 4 h and
costs amplitude, while its rhythmic outputs follow behaviour.

Human vastus lateralis after 3.5 days of 12 h-shifted sleep, light, activity and meals (Wefers 2018, GSE106800; 11
men; biopsies at 07:00 and 19:00 in the aligned and the misaligned condition). Clock faces: angle = peak time of the
misaligned 07:00-to-19:00 change against the aligned one (+ = later, clockwise; the A257 shift negated: the clock state
is further along its cycle, so it peaks earlier), radius = its amplitude relative to the aligned change.
A. Core loop (11 genes; A335 part 1a, filled), 14 clock genes (A257, open) and the rhythmic outputs (1,587
   reference-rhythmic genes, A261 H261.3, orange). Clock genes: excursion ratio |r| with r = sum dm conj(dc) / sum |dc|^2
   and the shift, 90 % subject bootstrap (`nightshift_compute.py`; the registered R projects on the expected chord and
   is not this ratio, see the A257 erratum). Outputs: A261's shift and excursion ratio with their 90 % intervals.
B. Core loop per subject (11 men): angle and |dm| / |dc|; the group value of A in blue.

Status: A257 registered (not blind), A335 part 1a registered (core loop blind), A261 registered; the excursion ratio
display is post hoc (user, 2026-11-25). Open data.
Data: presentation/output/nightshift/{summary,per_subject}.csv; exploratory/A261_nightshift_programs/output/programs.csv.
Output: presentation/figures/nightshift/nightshift.{svg,png} (wide).

    python scripts/nightshift.py
"""
from pathlib import Path

import pandas as pd
from matplotlib.lines import Line2D

import paperfig as pf
from holds_organisation import CS, OUTPUTS, REF, clockface, vector

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "presentation" / "output" / "nightshift"
OUT = PAPER / "presentation" / "figures" / "nightshift"
NAMES = {0: "aligned", 12: "shifted 12 h"}
LIGHT = "#9ecae1"


def clocks(ax, summ, prog):
    clockface(ax, 1.5, (0.5, 1, 1.5), "core clock and outputs", NAMES, ring_label=2.3)
    vector(ax, 1, 0, REF)
    for g, open_ in (("14 genes", True), ("core loop", False)):
        r = summ.loc[g]
        vector(ax, r.excursion_ratio, -r.shift_h, CS, None, (r.excursion_lo90, r.excursion_hi90),
               (-r.shift_hi90, -r.shift_lo90), open_=open_)
    vector(ax, prog.ratio, -prog.shift_h, OUTPUTS, None, (prog.ratio_lo, prog.ratio_hi), (-prog.shift_hi, -prog.shift_lo))


def subjects(ax, sub, summ):
    clockface(ax, 1.2, (0.5, 1), "core clock, per subject", NAMES, ring_label=2.3)
    vector(ax, 1, 0, REF)
    for r in sub.itertuples():
        vector(ax, r.excursion_ratio, -r.shift_h, LIGHT, size=34)
    r = summ.loc["core loop"]
    vector(ax, r.excursion_ratio, -r.shift_h, CS, size=44)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    summ = pd.read_csv(SRC / "summary.csv").set_index("gene_set")
    sub = pd.read_csv(SRC / "per_subject.csv")
    sub = sub[sub.gene_set == "core loop"]
    prog = pd.read_csv(PAPER / "exploratory/A261_nightshift_programs/output/programs.csv").set_index("program")
    lay = pf.Layout.slide("AB", kind="wide")
    clocks(lay.axes["A"], summ, prog.loc["ALL_RHYTHMIC"])
    subjects(lay.axes["B"], sub, summ)
    h = [Line2D([], [], ls="", marker="o", ms=8, color=REF, label="aligned (reference)"),
         Line2D([], [], ls="", marker="o", ms=8, color=CS, label="core loop"),
         Line2D([], [], ls="", marker="o", ms=8, mfc="white", mec=CS, mew=1.2, label="14 clock genes"),
         Line2D([], [], ls="", marker="o", ms=8, color=LIGHT, label="one subject (core loop)"),
         Line2D([], [], ls="", marker="o", ms=8, color=OUTPUTS, label="rhythmic outputs")]
    lay.shared_legend("AB", h, [x.get_label() for x in h], side="bottom", ncol=5, fontsize=10, handletextpad=0.2,
                      columnspacing=1.0)
    fig = lay.finalize()
    print(summ.round(3).to_string())
    print(sub[["shift_h", "excursion_ratio"]].describe().round(2).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "nightshift")


if __name__ == "__main__":
    main()
