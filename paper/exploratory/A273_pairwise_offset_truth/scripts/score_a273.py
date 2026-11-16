"""A273: are untimed pairwise offsets between cell types of one tissue sample unbiased? Mouse colon with known ZT.

Registered in paper/docs/registered_experiments.md (A273; a calibration of the A272 pairwise readout). Per arm of the
colon chrono-atlas (8 mice, ZT 1 / 7 / 13 / 19, two per ZT):
  timed    per compartment, `Dataset.timed_loadings("zt")`: phase against ZT (`clock_offset_h`) and coupling to the known
           time; timed offset of a pair = difference of the two phases
  untimed  arg of the complex correlation of the two compartments' clock-plane states over the mice (the A272 statistic,
           no time used)
  error    untimed minus timed (wrapped to +-12 h), leave-one-mouse-out jackknife SD

    python score_a273.py [--design]
Outputs: paper/exploratory/A273_pairwise_offset_truth/output/{design,compartments,pairs,tests}.csv
"""

from __future__ import annotations

import importlib.util
import itertools
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

spec = importlib.util.spec_from_file_location("a273_a272", PAPER / "exploratory/A272_backbone_pairwise/scripts/score_a272.py")
A272 = importlib.util.module_from_spec(spec)
sys.modules["a273_a272"] = A272
spec.loader.exec_module(A272)

AMENDMENT = "A273_pairwise_offset_truth"
ARMS, SENS = ("C-DC", "R-DC"), ("R-PC",)
MIN_CELLS, MIN_MICE, MIN_COUPLING = 30, 7, 0.5


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def arm_table(pb, donors, arm):
    d = donors.index[donors.arm == arm]
    t = pb[pb.donor.isin(d) & (pb.n_cells >= MIN_CELLS)].copy()
    t["stratum"] = arm
    n = t.groupby("compartment").donor.nunique()
    return t[t.compartment.isin(n.index[n >= MIN_MICE])].reset_index(drop=True)


def untimed_pairs(t, leave=None):
    """{(a, b): complex r} over the mice having both compartments (optionally without one mouse)."""
    if leave is not None:
        t = t[t.donor != leave].reset_index(drop=True)
    z, v = A272.atlas_states(t, cs.gtex())
    rows = {c: g.set_index("donor")["_i"] for c, g in t.assign(_i=np.arange(len(t))).groupby("compartment")}
    out = {}
    for a, b in itertools.combinations(sorted(rows), 2):
        d = rows[a].index.intersection(rows[b].index)
        if len(d) >= MIN_MICE - (leave is not None):
            pr = dict(ia=rows[a].loc[d].to_numpy(), ib=rows[b].loc[d].to_numpy(), s=np.zeros(len(d), int))
            out[(a, b)] = A272.pair_r(z, v, pr)
    return out


def timed(t, donors, leave=None):
    if leave is not None:
        t = t[t.donor != leave]
    ds = cs.Dataset(t, anchors=None, donors=donors, min_donors=MIN_MICE - (leave is not None))
    return ds.timed_loadings("zt")


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--design", action="store_true")
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(AMENDMENT)
    out.mkdir(parents=True, exist_ok=True)
    pb, donors = data.load("colon_chrono")
    donors = donors.set_axis(donors.index.astype(str))
    tabs = {a: arm_table(pb, donors, a) for a in ARMS + SENS}
    design = pd.DataFrame([dict(arm=a, mice=t.donor.nunique(), compartments=t.compartment.nunique(),
                                pairs=t.compartment.nunique() * (t.compartment.nunique() - 1) // 2,
                                zts=sorted(donors.zt[donors.arm == a].unique())) for a, t in tabs.items()])
    design.to_csv(out / "design.csv", index=False)
    print(design.to_string())
    if args.design:
        return

    comp_rows, pair_rows = [], []
    for arm, t in tabs.items():
        T = timed(t, donors)
        for c, r in T.iterrows():
            comp_rows.append(dict(arm=arm, compartment=c, n=r["n"], timed_phase_h=r["clock_offset_h"], timed_coupling=r["coupling"]))
        U = untimed_pairs(t)
        jk = {k: [] for k in U}
        for m in sorted(t.donor.unique()):
            Tm, Um = timed(t, donors, leave=m), untimed_pairs(t, leave=m)
            for (a, b) in U:
                if (a, b) in Um and a in Tm.index and b in Tm.index:
                    jk[(a, b)].append(wrap(np.angle(Um[(a, b)]) * HOURS - wrap(Tm.clock_offset_h[a] - Tm.clock_offset_h[b])))
        for (a, b), r in U.items():
            td = float(wrap(T.clock_offset_h[a] - T.clock_offset_h[b]))
            ud = float(np.angle(r) * HOURS)
            e = np.asarray(jk[(a, b)], float)
            n = len(e)
            pair_rows.append(dict(arm=arm, a=a, b=b, r_abs=abs(r), untimed_h=ud, timed_h=td, error_h=float(wrap(ud - td)),
                                  jackknife_sd_h=float(np.sqrt((n - 1) / n * np.sum((e - e.mean()) ** 2))) if n > 1 else np.nan,
                                  min_timed_coupling=float(min(T.coupling[a], T.coupling[b]))))
    comp, pairs = pd.DataFrame(comp_rows), pd.DataFrame(pair_rows)
    pairs["rhythmic"] = pairs.min_timed_coupling >= MIN_COUPLING
    comp.to_csv(out / "compartments.csv", index=False)
    pairs.to_csv(out / "pairs.csv", index=False)

    def summary(g, label):
        slope = float(np.sum(g.untimed_h * g.timed_h) / np.sum(g.timed_h ** 2)) if len(g) else np.nan
        return dict(scope=label, pairs=len(g), median_abs_error_h=float(g.error_h.abs().median()),
                    share_within_1h=float((g.error_h.abs() <= 1).mean()), slope_untimed_on_timed=slope,
                    pearson=float(np.corrcoef(g.untimed_h, g.timed_h)[0, 1]) if len(g) > 2 else np.nan,
                    median_abs_timed_h=float(g.timed_h.abs().median()), max_abs_timed_h=float(g.timed_h.abs().max()),
                    median_jackknife_sd_h=float(g.jackknife_sd_h.median()), median_r_abs=float(g.r_abs.median()))

    prim = pairs[pairs.arm.isin(ARMS) & pairs.rhythmic]
    rows = [summary(prim, "primary: C-DC + R-DC, both compartments timed-rhythmic")]
    rows += [summary(pairs[(pairs.arm == a) & pairs.rhythmic], f"{a}, rhythmic") for a in ARMS + SENS]
    rows += [summary(pairs[pairs.arm.isin(ARMS)], "C-DC + R-DC, all pairs")]
    tests = pd.DataFrame(rows)
    tests.loc[0, "V273_1_median_error_le_1h"] = bool(tests.median_abs_error_h[0] <= 1.0)
    tests.loc[0, "V273_2_slope_in_0.7_1.3"] = bool(0.7 <= tests.slope_untimed_on_timed[0] <= 1.3)
    tests.to_csv(out / "tests.csv", index=False)
    print(comp.round(3).to_string())
    print(tests.round(3).T.to_string())


if __name__ == "__main__":
    main()
