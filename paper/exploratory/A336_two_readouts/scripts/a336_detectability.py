"""A336 part 3 (registered, power only): a detection limit beside every "kept" claim of the main figures.

No observed value enters any computation: limits come from registered MDEs, null / bootstrap spreads and one simulation
on the reference group; observed values and intervals are copied beside them from the entries' output tables. No verdict.

Simulation (time kept in intensive care, A293 H293.2 order rho): GTEx fast deaths only, 200 random halvings (seed
[336, 3]); in half A one organ's clock states rotated by x h (each organ in turn, median over organs), or every organ by
N(0, x^2) h; detection when rho(A, B) (A293 / A132 Fisher-Lee estimator, pairs with >= 20 donors in both halves) falls
below the 5th percentile of the unshifted halvings; MDE = the smallest x on a 0.25 h grid with detection >= 0.8.

    python a336_detectability.py
Outputs (`../output/`): detectability.csv, detectability_order_sim.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
for d in ("A109_coupling_nulls", "A132_brain_death", "A317_braindeath_confirm"):
    sys.path.insert(0, str(EXP / d / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a317 as S  # noqa: E402
from score_a317 import B, H  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
W = 2 * np.pi / 24
Z80 = 2.49  # one-sided alpha 0.05, power 0.80
GRID = np.arange(0.25, 12.01, 0.25)


def order_sim(n_half=200):
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    fast = T.group_donors("fast")
    Z, M = S.raw(T, "clock", 0, T.R)
    rng = np.random.default_rng([336, 3])
    halves = []
    for _ in range(n_half):
        p = rng.permutation(fast)
        a, b = p[: len(p) // 2], p[len(p) // 2:]
        halves.append((a, b, B.pair_mask(M, a, b)))
    rho = lambda Zx, a, b, m: S.rho(Zx, M, b, a, m)  # noqa: E731
    base = np.array([rho(Z, a, b, m) for a, b, m in halves])
    thr = np.quantile(base, 0.05)
    rows = []
    for x in GRID:
        rot = np.exp(1j * x * W)
        for o, t in enumerate(H.BODY):
            det = []
            for a, b, m in halves:
                Zs = Z.copy()
                Zs[a, o] *= rot
                det.append(rho(Zs, a, b, m) < thr)
            rows.append(dict(design="one organ shifted", organ=t, x_h=x, detection=float(np.mean(det))))
        det = []
        for a, b, m in halves:
            Zs = Z.copy()
            Zs[a, :] *= np.exp(1j * rng.normal(0, x * W, len(H.BODY)))[None, :]
            det.append(rho(Zs, a, b, m) < thr)
        rows.append(dict(design="organ phases dispersed (SD)", organ="all", x_h=x, detection=float(np.mean(det))))
    sim = pd.DataFrame(rows)
    sim.attrs.update(base_median=float(np.median(base)), threshold=float(thr), n_half=n_half,
                     pairs_median=float(np.median([m.sum() for *_, m in halves])))
    return sim


def mde_from(sim, design):
    s = sim[sim.design == design]
    out = {}
    for o, g in s.groupby("organ"):
        ok = g[g.detection >= 0.8]
        out[o] = float(ok.x_h.min()) if len(ok) else np.nan
    return out


def main():
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    sim = order_sim()
    sim.to_csv(OUT / "detectability_order_sim.csv", index=False)
    one = mde_from(sim, "one organ shifted")
    disp = mde_from(sim, "organ phases dispersed (SD)")["all"]
    one_med = float(np.nanmedian(list(one.values()))) if np.isfinite(list(one.values())).any() else np.nan
    n_never = int(np.sum(~np.isfinite(list(one.values()))))
    rd = lambda p: pd.read_csv(EXP / p)  # noqa: E731
    g293 = rd("A293_organisation_levels/output/part1/groups.csv").set_index("test").loc["H293.2"]
    a235 = rd("A235_metabolic_offset/output/power.csv")
    a235 = a235[(a235.selection == "module") & (a235.readout == "offset_h") & (a235.set == "fast") & (a235.host == "A164")]
    a233 = rd("A233_chronic_conditions/output/power.csv")
    a233o = a233[(a233.selection == "module") & (a233.readout == "offset_h") & (a233.set == "pooled")]
    a233b = a233[(a233.selection == "bodywide") & (a233.readout == "d_star") & (a233.set == "pooled")]
    a164 = rd("A164_gtex_metabolic_matched/output/a164_power.csv")
    a164 = a164[(a164.statistic == "module") & (a164.set == "pooled")].set_index("design")
    p138 = rd("A138_cancer_estimability/output/power/power_summary.csv").set_index("hypothesis").loc["H138.4"]
    p254 = rd("A254_paired_lesions/output/power.csv").iloc[0]
    e293 = rd("A335_review_batch/output/part3_equivalence.csv").iloc[0]
    t295 = rd("A295_composition/output/tests.csv").set_index(rd("A295_composition/output/tests.csv").columns[0])
    eq315 = rd("A335_review_batch/output/equivalence.csv").set_index(["gene_set", "donors"])
    p315 = rd("A315_braindeath_gain/output/power.csv").set_index(rd("A315_braindeath_gain/output/power.csv").columns[0])["value"]
    m304 = rd("A304_cancer_gain_departure/output/power/mde.csv").set_index("contrast").loc["H304.1"]
    p306 = rd("A306_disease_gain_departure/output/power.csv")
    p306s = p306[p306.scored]
    se = lambda lo, hi: (hi - lo) / (2 * 1.645)  # noqa: E731
    # pairing: a dispersion of SD sigma shrinks a complex pairing by exp(-sigma^2 / 2); lower bound of the MDE in hours
    # from a reference pairing of 1 losing the registered MDE (larger for any reference below 1)
    sig = lambda frac: float(np.sqrt(-2 * np.log(1 - frac)) / W) if 0 < frac < 1 else np.nan  # noqa: E731
    rows = [
        dict(claim="time kept in intensive care (order between organs)", entry="A293 H293.2", statistic="Fisher-Lee rho of pair phases, ventilated post-mortem vs fast",
             observed=f"rho {g293.rho_clock:.3f} (matched q95 {g293.matched_q95:.3f}); no interval", mde_native="none registered",
             mde_hours=f"one organ shifted: median {one_med:.2f} h over organs ({n_never} of 12 never detected up to 12 h); all organs dispersed: SD {disp:.2f} h",
             mde_fold="", method="simulation on fast-death halvings (no ventilated value used)",
             positive_control="A317 brain death D_rho 0.621 (core 0.465, A336 part 4); A201 aorta jet lag -2.86 h; A319 central clock lost",
             time_mde_h=disp, loud_mde_fold=np.nan),
        dict(claim="chronic disease: time (metabolic offset)", entry="A235 H235.1", statistic="module offset vs own GTEx tissues (h), fast deaths",
             observed="T2D +0.02 h (-0.25 to +0.30); obesity +0.20 h (-0.23 to +0.64) (A235 score_tests)",
             mde_native="; ".join(f"{r.design} {r.mde_80:.3f} h" for r in a235.itertuples()), mde_hours="as native", mde_fold="",
             method="registered (A235 power.csv, fast, host A164)", positive_control="A201 aorta jet lag -2.86 h; A257 night shift 4.06 h (core)",
             time_mde_h=float(a235.mde_80.max()), loud_mde_fold=np.nan),
        dict(claim="chronic disease: time (module offset, five conditions)", entry="A233 H233.c.2-3", statistic="module offset (h), pooled",
             observed="no condition detected (A233 score_tests)", mde_native="; ".join(f"{r.design} {r.mde_80:.2f} h" for r in a233o.itertuples()),
             mde_hours="as native", mde_fold="", method="registered (A233 power.csv)", positive_control="A201 -2.86 h; A257 4.06 h (core)",
             time_mde_h=float(a233o.mde_80.median()), loud_mde_fold=np.nan),
        dict(claim="chronic disease: coupling (body-wide and module)", entry="A233 H233.c.1; A164 H164.1-2", statistic="coupling difference d* at common spread",
             observed="A164 T2D +0.019 (-0.021 to +0.060), obesity +0.004 (-0.081 to +0.089); A233 all five within +-0.05",
             mde_native="A164 " + "; ".join(f"{k} {v:.3f}" for k, v in a164.mde_80.items()) + "; A233 body-wide " + f"{a233b.mde_80.min():.3f}-{a233b.mde_80.max():.3f}",
             mde_hours="not available (coupling)", mde_fold="not available (coupling)", method="registered",
             positive_control="none on this statistic", time_mde_h=np.nan, loud_mde_fold=np.nan),
        dict(claim="chronic disease: loudness and departure", entry="A306", statistic="log gain / log departure ratio vs companions",
             observed="all scored rows not detected (A306 score_tests)", mde_native="; ".join(f"{r.condition} {r.readout} {r.mde_80:.3f}" for r in p306s.itertuples()),
             mde_hours="", mde_fold="; ".join(f"{r.condition} {r.readout} x{np.exp(r.mde_80):.2f}" for r in p306s.itertuples()),
             method="registered (A306 power.csv), fold = exp(log MDE)", positive_control="A315 ICU loudness Q* +0.404 (amplitude about 0.82); A257 R 0.239",
             time_mde_h=np.nan, loud_mde_fold=float(np.exp(p306s[p306s.readout == "gain"].mde_80.max()))),
        dict(claim="tumours keep the patient's time", entry="A138 H138.4", statistic="tumour-normal complex pairing (pooled)",
             observed="pairing 0.365 (null q95 0.157; matched q95 0.261; responsive q95 0.310)", mde_native=f"pairing {p138.mde:.3f} (vs 0)",
             mde_hours=f">= {sig(p138.mde):.2f} h dispersion SD (lower bound: reference pairing 1 losing {p138.mde:.3f})",
             mde_fold="", method="exp(-sigma^2/2) shrinkage of a complex pairing", positive_control="A201 -2.86 h (constant shift: rotation-invariant pairing does not see it)",
             time_mde_h=sig(p138.mde), loud_mde_fold=np.nan),
        dict(claim="chronic lesions keep time", entry="A254 H254.1", statistic="pooled lesion offset vs own tissue (h)",
             observed="+0.01 h (A254 score_pooled)", mde_native=f"{p254.pooled_offset_mde_h:.3f} h", mde_hours=f"{p254.pooled_offset_mde_h:.2f} h", mde_fold="",
             method="registered (2.49 x pooled null SD)", positive_control="A308 psoriasis -1.45 h (core -1.00 h); A201 -2.86 h",
             time_mde_h=float(p254.pooled_offset_mde_h), loud_mde_fold=np.nan),
        dict(claim="cell classes share time", entry="A293 H293.3", statistic="pooled endothelium minus fibroblast offset (h)",
             observed=f"{e293.estimate:.2f} h ({e293.lo90:.2f} to {e293.hi90:.2f}); not blind",
             mde_native=f"{Z80 * se(e293.lo90, e293.hi90):.2f} h", mde_hours=f"{Z80 * se(e293.lo90, e293.hi90):.2f} h", mde_fold="",
             method="2.49 x SE from the atlas-bootstrap interval width", positive_control="mouse colon untimed offset 1.33 h (A273); A279 aorta",
             time_mde_h=Z80 * se(e293.lo90, e293.hi90), loud_mde_fold=np.nan),
        dict(claim="composition does not move time", entry="A295 H295.2", statistic="pooled phase change per SD of lymphocyte score (h)",
             observed=f"{t295.loc['lymphocyte', 'pooled_phase_change_h']:+.3f} h per SD ({t295.loc['lymphocyte', 'phase_lo90']:+.3f} to {t295.loc['lymphocyte', 'phase_hi90']:+.3f})"
             if "phase_lo90" in t295 else "see A295 tests.csv",
             mde_native=f"{Z80 * se(t295.loc['lymphocyte', 'phase_lo90'], t295.loc['lymphocyte', 'phase_hi90']):.3f} h per SD",
             mde_hours=f"{Z80 * se(t295.loc['lymphocyte', 'phase_lo90'], t295.loc['lymphocyte', 'phase_hi90']):.2f} h per SD", mde_fold="",
             method="2.49 x SE from the interval width", positive_control="none on this statistic",
             time_mde_h=Z80 * se(t295.loc["lymphocyte", "phase_lo90"], t295.loc["lymphocyte", "phase_hi90"]), loud_mde_fold=np.nan),
        dict(claim="composition: loudness per SD", entry="A295 H295.1", statistic="pooled relative-gain change per SD of lymphocyte score",
             observed=f"{t295.loc['lymphocyte', 'pooled_gain_change']:+.4f} (SE {t295.loc['lymphocyte', 'gain_se']:.4f})",
             mde_native=f"{Z80 * t295.loc['lymphocyte', 'gain_se']:.3f} per SD", mde_hours="",
             mde_fold=f"x{np.exp(Z80 * t295.loc['lymphocyte', 'gain_se']):.3f} per SD", method="2.49 x registered SE",
             positive_control="stromal-vascular score +0.041 per SD (A295, same statistic)", time_mde_h=np.nan,
             loud_mde_fold=float(np.exp(Z80 * t295.loc["lymphocyte", "gain_se"]))),
        dict(claim="brain-death loudness equal", entry="A315 H315.1 (A335 part 3)", statistic="Q* (clock minus median matched log variance ratio)",
             observed="; ".join(f"{g} {d}: {r.Q_star:+.3f} ({r.lo90:+.3f} to {r.hi90:+.3f})" for (g, d), r in eq315.iterrows()),
             mde_native=f"registered: power {float(p315['power']):.3f} at Q* {float(p315['target']):.3f}; 2.49 x SD {Z80 * float(p315['sd_used']):.3f}",
             mde_hours="", mde_fold=f"amplitude x{np.exp(-Z80 * float(p315['sd_used']) / 2):.2f} (loss) / x{np.exp(Z80 * float(p315['sd_used']) / 2):.2f} (gain)",
             method="fold = exp(MDE / 2) (variance to amplitude)", positive_control="A315 ICU loss Q* +0.404 (core +0.457), same statistic and cohort",
             time_mde_h=np.nan, loud_mde_fold=float(np.exp(Z80 * float(p315["sd_used"]) / 2))),
        dict(claim="the departure does not grow (tumours)", entry="A304 H304.1", statistic="log departure ratio vs own normal",
             observed="Q_R -0.190, p 0.944 (A304 tests)", mde_native=f"departure x{m304.mde_k_D} (power at x2 {m304.power_dep_at_2})",
             mde_hours="", mde_fold=f"x{m304.mde_k_D}", method="registered simulation grid", positive_control="A304 gain loss x0.5 detected (power 0.84); none for a departure rise",
             time_mde_h=np.nan, loud_mde_fold=float(m304.mde_k_D)),
    ]
    tab = pd.DataFrame(rows)
    tab["reading_rule"] = "an arrow inside its detection circle reads 'not detectable', not 'no change'"
    tab.to_csv(OUT / "detectability.csv", index=False)
    print(f"order sim: base median rho {sim.attrs['base_median']:.3f}, threshold {sim.attrs['threshold']:.3f}, "
          f"pairs {sim.attrs['pairs_median']:.0f}; one organ MDE per organ {one}; dispersion SD MDE {disp}", flush=True)
    pd.set_option("display.width", 250, "display.max_colwidth", 120)
    print(tab[["claim", "mde_native", "mde_hours", "mde_fold"]].to_string(index=False))


if __name__ == "__main__":
    main()
