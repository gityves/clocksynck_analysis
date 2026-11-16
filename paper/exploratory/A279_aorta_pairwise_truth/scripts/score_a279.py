"""A279: untimed pairwise offsets against known ZT in the mouse aorta, and macrophage against stromal-vascular timed
coupling. Auerbach 2026 (pools of 2 aortas, ZT 0/6/12/18; A201 loader).

Registered in paper/docs/registered_experiments.md (A279). As A273, per arm (aligned primary; misaligned, Bmal1 KO
descriptive): library x cell type pseudobulks (SMC, fibroblast, EC, macrophage; >= 30 cells; compartments in >= 10
libraries of the arm, or all libraries when the arm has fewer than 10):
  timed    `Dataset.timed_loadings("zt")`: phase and coupling to known time per compartment
  untimed  arg of the A272 pairwise statistic over libraries (no time used); error = untimed - timed (wrapped)
  V279.1   median |error| <= 1 h over aligned pairs with both timed couplings >= 0.5
  H279.1   aligned: macrophage timed coupling minus the smallest stromal-vascular one < 0, jackknife 90 % excluding 0

    python score_a279.py [--design]
Outputs: paper/exploratory/A279_aorta_pairwise_truth/output/{design,compartments,pairs,tests}.csv
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
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


A272 = _load("a279_a272", PAPER / "exploratory/A272_backbone_pairwise/scripts/score_a272.py")
A201 = _load("a279_a201", PAPER / "exploratory/A201_aorta_jetlag/scripts/score_a201.py")

AMENDMENT = "A279_aorta_pairwise_truth"
CLASSES = ("SMC", "Fibroblast", "EC", "Macrophage")  # labels of the source (A201 MAIN + EXTRA)
SV = ("SMC", "Fibroblast", "EC")
ARMS = {"aligned": "ali", "misaligned": "mis", "bmal1_ko": "ko"}
MIN_CELLS, MIN_LIBS, MIN_COUPLING = 30, 10, 0.5


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def tables():
    o, X, genes = A201.load()
    keep = o.cell_type.astype(str).isin(CLASSES).to_numpy()
    lib, C = A201.pseudobulk(o[keep].reset_index(drop=True), X[keep])
    cnt = pd.DataFrame(C, columns=[g.upper() if g.upper() in cs.gtex().genes else g for g in genes])
    out = {}
    for arm, cond in ARMS.items():
        m = (lib.cond == cond).to_numpy() & (lib.cells >= MIN_CELLS).to_numpy()
        key = lib[m][["sample_id", "cell_type"]].astype(str).reset_index(drop=True)
        pb = cs.from_bulk(cnt[m].round().reset_index(drop=True), key, donor="sample_id", compartment="cell_type",
                          depth=pd.Series(lib.total.to_numpy()[m]))
        pb["stratum"] = arm
        pb["n_cells"] = lib.cells.to_numpy()[m]
        n_lib = lib[lib.cond == cond].sample_id.nunique()
        need = min(MIN_LIBS, n_lib)
        n = pb.groupby("compartment").donor.nunique()
        pb = pb[pb.compartment.isin(n.index[n >= need])].reset_index(drop=True)
        don = lib[lib.cond == cond].groupby("sample_id").agg(zt=("zt", "first"), sex=("sex", "first"))
        don.index = don.index.astype(str)
        out[arm] = (pb, don, need)
    return out


def untimed_pairs(t, need, leave=None):
    if leave is not None:
        t = t[t.donor != leave].reset_index(drop=True)
    z, v = A272.atlas_states(t, cs.gtex())
    rows = {c: g.set_index("donor")["_i"] for c, g in t.assign(_i=np.arange(len(t))).groupby("compartment")}
    out = {}
    for a, b in itertools.combinations(sorted(rows), 2):
        d = rows[a].index.intersection(rows[b].index)
        if len(d) >= need - (leave is not None):
            out[(a, b)] = A272.pair_r(z, v, dict(ia=rows[a].loc[d].to_numpy(), ib=rows[b].loc[d].to_numpy(),
                                                 s=np.zeros(len(d), int)))
    return out


def timed(t, don, need, leave=None):
    if leave is not None:
        t = t[t.donor != leave]
    return cs.Dataset(t, anchors=None, donors=don, min_donors=need - (leave is not None)).timed_loadings("zt")


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--design", action="store_true")
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(AMENDMENT)
    out.mkdir(parents=True, exist_ok=True)
    tabs = tables()
    design = pd.DataFrame([dict(arm=a, libraries=t.donor.nunique(), min_libraries=need,
                                compartments=";".join(sorted(t.compartment.unique())),
                                libs_per_compartment=t.groupby("compartment").donor.nunique().to_dict(),
                                zts=sorted(don.zt.unique())) for a, (t, don, need) in tabs.items()])
    design.to_csv(out / "design.csv", index=False)
    print(design.to_string())
    if args.design:
        return

    comp_rows, pair_rows, h1 = [], [], {}
    for arm, (t, don, need) in tabs.items():
        T = timed(t, don, need)
        for c, r in T.iterrows():
            comp_rows.append(dict(arm=arm, compartment=c, n=r["n"], timed_phase_h=r["clock_offset_h"], timed_coupling=r["coupling"]))
        U = untimed_pairs(t, need)
        jk, jh = {k: [] for k in U}, []
        for m in sorted(t.donor.unique()):
            Tm, Um = timed(t, don, need, leave=m), untimed_pairs(t, need, leave=m)
            for (a, b) in U:
                if (a, b) in Um and a in Tm.index and b in Tm.index:
                    jk[(a, b)].append(wrap(np.angle(Um[(a, b)]) * HOURS - wrap(Tm.clock_offset_h[a] - Tm.clock_offset_h[b])))
            if "Macrophage" in Tm.index and all(c in Tm.index for c in SV):
                jh.append(Tm.coupling["Macrophage"] - min(Tm.coupling[c] for c in SV))
        for (a, b), r in U.items():
            td, ud = float(wrap(T.clock_offset_h[a] - T.clock_offset_h[b])), float(np.angle(r) * HOURS)
            e = np.asarray(jk[(a, b)], float)
            n = len(e)
            pair_rows.append(dict(arm=arm, a=a, b=b, r_abs=abs(r), untimed_h=ud, timed_h=td, error_h=float(wrap(ud - td)),
                                  jackknife_sd_h=float(np.sqrt((n - 1) / n * np.sum((e - e.mean()) ** 2))) if n > 1 else np.nan,
                                  min_timed_coupling=float(min(T.coupling[a], T.coupling[b]))))
        if arm == "aligned" and "Macrophage" in T.index and all(c in T.index for c in SV):
            d = float(T.coupling["Macrophage"] - min(T.coupling[c] for c in SV))
            e = np.asarray(jh)
            n = len(e)
            sd = float(np.sqrt((n - 1) / n * np.sum((e - e.mean()) ** 2)))
            h1 = dict(H279_1_diff=d, H279_1_jackknife_sd=sd, H279_1_lo90=d - 1.645 * sd, H279_1_hi90=d + 1.645 * sd,
                      H279_1_pass=bool(d < 0 and d + 1.645 * sd < 0))
    comp, pairs = pd.DataFrame(comp_rows), pd.DataFrame(pair_rows)
    pairs["rhythmic"] = pairs.min_timed_coupling >= MIN_COUPLING
    comp.to_csv(out / "compartments.csv", index=False)
    pairs.to_csv(out / "pairs.csv", index=False)

    def summary(g, label):
        span = float(g.timed_h.max() - g.timed_h.min()) if len(g) else np.nan
        slope = float(np.sum(g.untimed_h * g.timed_h) / np.sum(g.timed_h ** 2)) if len(g) and span >= 1 else np.nan
        return dict(scope=label, pairs=len(g), median_abs_error_h=float(g.error_h.abs().median()) if len(g) else np.nan,
                    share_within_1h=float((g.error_h.abs() <= 1).mean()) if len(g) else np.nan, timed_span_h=span,
                    slope_untimed_on_timed=slope, median_abs_timed_h=float(g.timed_h.abs().median()) if len(g) else np.nan,
                    median_jackknife_sd_h=float(g.jackknife_sd_h.median()) if len(g) else np.nan,
                    median_r_abs=float(g.r_abs.median()) if len(g) else np.nan)

    rows = [summary(pairs[(pairs.arm == "aligned") & pairs.rhythmic], "primary: aligned, both timed-rhythmic")]
    rows += [summary(pairs[pairs.arm == a], f"{a}, all pairs") for a in ARMS]
    tests = pd.DataFrame(rows)
    tests.loc[0, "V279_1_median_error_le_1h"] = bool(tests.median_abs_error_h[0] <= 1.0) if tests.pairs[0] else None
    for k, v in h1.items():
        tests.loc[0, k] = v
    tests.to_csv(out / "tests.csv", index=False)
    print(comp.round(3).to_string())
    print(pairs.round(3).to_string())
    print(tests.round(3).T.to_string())


if __name__ == "__main__":
    main()
