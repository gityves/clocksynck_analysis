"""T1 / T2 (PAPER_PLAN_2026-11-19 section 6.2, tier 0): assemble the two-axis table and the polar-plot data.

No new statistic and no resampling: every number is read from a scored entry's output CSV and converted as in
PAPER_PLAN section 10. Aggregates only (GTEx controlled access: no per-donor row is read into the outputs).

Conversions (`derived` in the table):
- variance ratio to amplitude: sqrt(exp(log ratio));
- A304 gain: exp(Q_G), already on the amplitude scale; A306 gain: log z2 ratio, amplitude exp(delta / 2);
- A257 excursion = R / cos(s), s in radians;
- Delta^2 = (1 - R)^2 + 4 R sin^2(phi / 2); loudness share = (1 - R)^2 / Delta^2. Where the condition gives an order
  statistic and no shift, phi is taken as 0 and flagged (`phi_assumed_0`).

    python assemble_two_axis.py
Outputs (`output/`): two_axis_table.csv, polar_points.csv.
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EX = HERE.parent
RES = EX.parent / "results"
OUT = HERE / "output"


def rd(p):
    return pd.read_csv(p)


def rel(p):
    return str(Path(p).resolve().relative_to(EX.parents[1]))


def delta(R, h, assumed=False):
    if R is None or not np.isfinite(R):
        return np.nan, np.nan
    phi = 0.0 if (h is None or not np.isfinite(h)) else 2 * np.pi * h / 24
    loud, time = (1 - R) ** 2, 4 * R * np.sin(phi / 2) ** 2
    d2 = loud + time
    return float(np.sqrt(d2)), float(loud / d2) if d2 > 0 else np.nan


def row(condition, entry, src, time_kind="", time=np.nan, time_lo=np.nan, time_hi=np.nan, amp=np.nan, amp_lo=np.nan,
        amp_hi=np.nan, amp_kind="", departure="", status="", note="", species="human"):
    assumed = time_kind != "shift_h"
    d, share = delta(amp, time if time_kind == "shift_h" else np.nan)
    return dict(condition=condition, species=species, time_kind=time_kind, time=time, time_lo90=time_lo,
                time_hi90=time_hi, amplitude_ratio=amp, amplitude_lo90=amp_lo, amplitude_hi90=amp_hi,
                amplitude_kind=amp_kind, departure=departure, Delta=d, Delta_loudness_share=share,
                phi_assumed_0=bool(assumed and np.isfinite(d)), status=status, entry=entry, source_csv=src, note=note)


def table():
    rows = []
    # intensive care (fast deaths -> ventilated post-mortem, living brain): order A293 H293.2, loudness A315 descriptive
    g = rd(EX / "A293_organisation_levels/output/part1/groups.csv").set_index("group")
    d315 = rd(EX / "A315_braindeath_gain/output/descriptive.csv").set_index("contrast")
    qs = float(d315.loc["fast_vs_postmortem", "Q_star"])
    rows.append(row("intensive care, living brain (vs fast deaths)", "A293 H293.2; A315 descriptive",
                    f"{rel(EX / 'A293_organisation_levels/output/part1/groups.csv')}; {rel(EX / 'A315_braindeath_gain/output/descriptive.csv')}",
                    "order rho (Fisher-Lee) kept", float(g.loc["ventilated_postmortem", "rho_clock"]),
                    amp=float(np.exp(-qs / 2)), amp_kind="sqrt variance ratio, clock-specific (Q* = Q - matched median)",
                    status="registered re-test; not blind; loudness descriptive",
                    note=f"Q* fast over ventilated {qs:.3f}; no interval in the outputs"))
    t315 = rd(EX / "A315_braindeath_gain/output/tests.csv").iloc[0]
    t317 = rd(EX / "A317_braindeath_confirm/output/tests.csv").set_index("test")
    rows.append(row("brain death (organ donors vs ventilated post-mortem)", "A315 H315.1; A317 H317.1",
                    f"{rel(EX / 'A317_braindeath_confirm/output/tests.csv')}; {rel(EX / 'A315_braindeath_gain/output/tests.csv')}",
                    "order changed: D_rho", float(t317.loc["H317.1", "D_rho"]),
                    amp=float(np.sqrt(t315.ratio_star)), amp_kind="sqrt variance ratio, clock-specific (Q*)",
                    status="aliased with procurement; not blind",
                    note=f"Q* {t315.Q_star:.3f}; raw Q {t315.Q:.3f} (90% {t315.boot_lo90:.3f} to {t315.boot_hi90:.3f}, raw scale); "
                         "no shift: Delta is the loudness term only"))
    # A257 night shift
    s = rd(EX / "A257_nightshift_muscle/output/score_tests.csv").set_index("test")
    h2, h3 = s.loc["H257.2 shift of the muscle clock (h)"], s.loc["H257.3 damping rho_mis / rho_ctl"]
    R, sh = float(h3.R), float(h2.s)
    rows.append(row("human muscle, day 3 of a 12 h behavioural shift (vs aligned)", "A257 H257.2, H257.3",
                    rel(EX / "A257_nightshift_muscle/output/score_tests.csv"), "shift_h", sh, float(h2.lo90),
                    float(h2.hi90), amp=R / np.cos(2 * np.pi * sh / 24), amp_kind="excursion = R / cos s (derived)",
                    status="not blind; 11 men",
                    note=f"R (projection on aligned chord) {R:.3f}, 90% {h3.lo90:.3f} to {h3.hi90:.3f}; excursion interval not in outputs"))
    # mouse conflicting feeding
    for f, lab, ent in [("A322_conflicting_drives/output/ratios.csv", "Manella 2021, daytime vs ad libitum", "A322"),
                        ("A322_conflicting_drives/output/xin/ratios.csv", "Xin 2021, day- vs night-restricted", "A324")]:
        for _, r in rd(EX / f).iterrows():
            rows.append(row(f"mouse {r.tissue}, {lab} feeding", ent, rel(EX / f), "shift_h", float(r.shift_h),
                            amp=float(r.amp_ratio), amp_lo=float(r.amp_lo90), amp_hi=float(r.amp_hi90),
                            amp_kind="timed amplitude ratio", species="mouse",
                            status="A324 blind for heart, kidney, VAT" if ent == "A324" else "not blind",
                            note="shift interval not in ratios.csv"))
    t297 = rd(EX / "A297_feeding_model/output/tests.csv")
    r = t297[t297.contrast.str.startswith("observed Greenwell")].iloc[0]
    rows.append(row("mouse liver, feeding rhythm removed (arrhythmic vs ad libitum)", "A297",
                    rel(EX / "A297_feeding_model/output/tests.csv"), "shift_h", float(r.shift_h), float(r.shift_lo90),
                    float(r.shift_hi90), float(r.amp_ratio), float(r.amp_lo90), float(r.amp_hi90),
                    "timed amplitude ratio", species="mouse", status="not blind"))
    # tumours
    t304 = rd(EX / "A304_cancer_gain_departure/output/score/tests.csv")
    gain = t304[(t304.contrast == "H304.1") & (t304.readout == "gain")].iloc[0]
    dep = t304[(t304.contrast == "H304.1") & (t304.readout == "departure")].iloc[0]
    pw = rd(EX / "A304_cancer_gain_departure/output/power/power.csv")
    pdbl = float(pw[(pw.contrast == "H304.1") & (pw.scenario == "departure x2.0")].power_dep_05.iloc[0])
    rows.append(row("tumours, high vs low proliferation (own normal tissue)", "A138; A304 H304.1",
                    f"{rel(EX / 'A304_cancer_gain_departure/output/score/tests.csv')}; {rel(EX / 'A304_cancer_gain_departure/output/power/power.csv')}",
                    "time kept (tumours pair with their patient, A138)", amp=float(np.exp(gain.statistic)),
                    amp_kind="gain exp(Q_G), clock-specific",
                    departure=f"not raised: exp(Q_R) {np.exp(dep.statistic):.2f}, p {dep.p:.2f}; power {pdbl:.3f} at a doubling",
                    status="A304 not blind for A138", note=f"Q_G {gain.statistic:.3f}, p {gain.p:.4f}; no interval in the outputs"))
    # chronic lesions
    l = rd(EX / "A254_paired_lesions/output/score_pooled.csv").iloc[0]
    rows.append(row("chronic lesions (periodontitis, osteoarthritis, carotid plaque, fibroid, hidradenitis; own tissue)",
                    "A254", rel(EX / "A254_paired_lesions/output/score_pooled.csv"), "shift_h", float(l.pooled_offset_h),
                    amp=float(np.exp(l.pooled_r / 2)),
                    amp_lo=float(np.exp((l.pooled_r - 1.645 * l.pooled_r_se) / 2)),
                    amp_hi=float(np.exp((l.pooled_r + 1.645 * l.pooled_r_se) / 2)),
                    amp_kind="sqrt variance ratio; interval r +- 1.645 se (derived from the reported se)",
                    status="powered for time (MDE %.2f h)" % l.mde_h,
                    note=f"log variance ratio r {l.pooled_r:.3f} (se {l.pooled_r_se:.3f}); time interval not in outputs; clock-specific only pooled, at the line"))
    # psoriasis
    t183 = rd(EX / "A183_skin_gse154200/output/tests.csv").set_index("hypothesis")
    pp = rd(EX / "A183_skin_gse154200/output/h183_1_paired_pooled.csv").set_index("Unnamed: 0")["value"]
    rf, se = float(pp["fixed"]), float(pp["se_fixed"])
    rows.append(row("psoriasis lesion (own skin)", "A183 H183.1; A308",
                    f"{rel(EX / 'A183_skin_gse154200/output/tests.csv')}; {rel(EX / 'A183_skin_gse154200/output/h183_1_paired_pooled.csv')}; {rel(EX / 'A308_psoriasis_calibration/output/tests.csv')}",
                    "shift_h", float(t183.loc["H183.1", "statistic_h"]), amp=float(np.exp(rf / 2)),
                    amp_lo=float(np.exp((rf - 1.645 * se) / 2)), amp_hi=float(np.exp((rf + 1.645 * se) / 2)),
                    amp_kind="sqrt variance ratio, no companions (descriptive); interval from fixed-effect se (derived)",
                    status="not blind; clock-specific after A308 calibration; loudness: check D1",
                    note=f"negative = lesion later (A183 convention); r fixed {rf:.3f} (se {se:.3f})"))
    a115 = rd(EX / "A115_A116_psoriasis/output/a115_d115_1.csv").set_index("series").loc["pooled"]
    rows.append(row("psoriasis, uninvolved vs healthy skin", "A115", rel(EX / "A115_A116_psoriasis/output/a115_d115_1.csv"),
                    "shift_h", float(a115.D115_h), float(a115.ci_lo), float(a115.ci_hi), status="descriptive (trait)",
                    note="negative = uninvolved later; interval is the bootstrap CI of the output"))
    # chronic systemic disease
    t235 = rd(EX / "A235_metabolic_offset/output/score_tests.csv")
    t306 = rd(EX / "A306_disease_gain_departure/output/score_tests.csv").set_index(["condition", "readout"])
    for cond, key, g306 in [("T2D (metabolic organs)", "t2d", "t2d"), ("obesity (metabolic organs)", "obese", None)]:
        r = t235[(t235.condition == key) & (t235.host == "A164")].iloc[0]
        gd = t306.loc[(g306, "gain")] if g306 else None
        dd = t306.loc[(g306, "departure")] if g306 else None
        rows.append(row(cond, "A164; A235; A306" if g306 else "A164; A235",
                        f"{rel(EX / 'A235_metabolic_offset/output/score_tests.csv')}" + (f"; {rel(EX / 'A306_disease_gain_departure/output/score_tests.csv')}" if g306 else ""),
                        "shift_h", float(r.fast_offset_h), float(r.fast_lo90), float(r.fast_hi90),
                        amp=float(np.exp(gd.delta / 2)) if g306 else np.nan,
                        amp_kind="exp(delta log z2 / 2), module minus companion" if g306 else "",
                        departure=(f"not detected: delta log psi {dd.delta:+.3f} (MDE {dd.mde_80:.3f})" if g306 else ""),
                        status="powered nulls; equivalence +-0.15 on d* (A164)",
                        note=f"offset in fast deaths, MDE {r.fast_mde_80:.2f} h" + (f"; gain MDE {gd.mde_80:.3f} (log z2)" if g306 else "")))
    t233 = rd(EX / "A233_chronic_conditions/output/score_pooled.csv")
    r = t233[(t233.design == "htn") & (t233.selection == "bodywide") & (t233.readout == "offset_h") & (t233.set == "pooled")].iloc[0]
    gd, dd = t306.loc[("htn", "gain")], t306.loc[("htn", "departure")]
    rows.append(row("hypertension (arteries, heart)", "A233; A306",
                    f"{rel(EX / 'A233_chronic_conditions/output/score_pooled.csv')}; {rel(EX / 'A306_disease_gain_departure/output/score_tests.csv')}",
                    "shift_h", float(r.estimate), float(r.lo90), float(r.hi90), amp=float(np.exp(gd.delta / 2)),
                    amp_kind="exp(delta log z2 / 2), module minus companion",
                    departure=f"not detected: delta log psi {dd.delta:+.3f} (MDE {dd.mde_80:.3f})", status="powered nulls",
                    note=f"body-wide pooled offset; gain MDE {gd.mde_80:.3f} (log z2)"))
    gd, dd = t306.loc[("smoking", "gain")], t306.loc[("smoking", "departure")]
    rows.append(row("smoking (lung)", "A306", rel(EX / "A306_disease_gain_departure/output/score_tests.csv"),
                    amp=float(np.exp(gd.delta / 2)), amp_kind="exp(delta log z2 / 2); descriptive (not powered)",
                    departure=f"not detected: delta log psi {dd.delta:+.3f} (MDE {dd.mde_80:.3f})", status="gain not powered"))
    for lab, f, ent in [("sex (women vs men)", "A208_sex_coupling/output/tests.csv", "A208"),
                        ("age", "A231_age_coupling/output/score_tests.csv", "A231")]:
        t = rd(EX / f)
        t = t[t.iloc[:, 1] == "offset_h"].iloc[0]
        rows.append(row(lab, ent, rel(EX / f), "shift_h", float(t.estimate), float(t.lo90), float(t.hi90),
                        status="d* equivalent within +-0.15; offset not detected",
                        note="CHIRAL's age dampening shrinks under death-type matching (A231)" if ent == "A231" else ""))
    return pd.DataFrame(rows)


def polar():
    rows = []
    org = rd(RES / "figure_2/output/organs/tissues.csv")
    org = org[org.cohort == "fast"]  # fast deaths, as the A104 / A294b rows
    t294 = rd(EX / "A294_floor/output/tissues.csv").set_index("tissue")
    t291 = rd(EX / "A291_attenuation/output/tissues.csv").set_index("tissue")
    t312 = rd(EX / "A312_elasticity/output/tissues.csv").set_index("tissue")
    o = org.set_index("tissue")
    src_o = f"{rel(EX / 'A294_floor/output/tissues.csv')}; {rel(EX / 'A291_attenuation/output/tissues.csv')}; {rel(RES / 'figure_2/output/organs/tissues.csv')}; {rel(EX / 'A312_elasticity/output/tissues.csv')}"
    for t in t294.index:
        r = dict(kind="organ", unit=t, gain_sqrtS=float(np.sqrt(t294.loc[t, "S"])), S=float(t294.loc[t, "S"]),
                 coupling=float(t291.loc[t, "coupling"]), psi=float(t291.loc[t, "psi"]),
                 departure_own=float(t294.loc[t, "own"]), departure_log_ratio_vs_matched=float(t294.loc[t, "log_ratio"]))
        if t in o.index:
            r.update(offset_h=float(o.loc[t, "untimed_offset_h"]), offset_lo=float(o.loc[t, "untimed_offset_h_lo"]),
                     offset_hi=float(o.loc[t, "untimed_offset_h_hi"]))
        if t in t312.index:
            r.update(elasticity=float(t312.loc[t, "elasticity"]), elasticity_lo90=float(t312.loc[t, "lo90"]),
                     elasticity_hi90=float(t312.loc[t, "hi90"]), person_leverage=float(t312.loc[t, "leverage"]))
        r.update(entry="A294b; A291; A104 offsets; A312", source_csv=src_o)
        rows.append(r)
    fam = rd(EX / "A293_organisation_levels/output/part3/families.csv")
    for _, f in fam.iterrows():
        rows.append(dict(kind="class", unit=f.family, class_log_gain_vs_endothelium=float(f.pooled),
                         class_gain_vs_endothelium=float(np.exp(f.pooled)), class_se=float(f.se_pooled),
                         class_tau=float(f.tau), atlases=int(f.atlases), entry="A293 Part 3 (H293.4)",
                         source_csv=rel(EX / "A293_organisation_levels/output/part3/families.csv"),
                         note="pooled log gain relative to endothelium (scale as in the entry)"))
    t312s = rd(EX / "A312_elasticity/output/power.csv").set_index("Unnamed: 0")["value"]
    tt = rd(EX / "A312_elasticity/output/tests.csv").set_index("test")
    rows.append(dict(kind="person", unit="shared amplitude across persons",
                     person_leverage=float(t312s["median_leverage"]), person_leverage_min=float(t312s["min_leverage"]),
                     elasticity=float(tt.loc["descriptive", "median_elasticity"]), entry="A312",
                     source_csv=f"{rel(EX / 'A312_elasticity/output/power.csv')}; {rel(EX / 'A312_elasticity/output/tests.csv')}",
                     note="leverage = corr of two independent half-host radii over donors (share of radius variance that is person-level); no amplitude CV in the outputs"))
    t294t = rd(EX / "A294_floor/output/tests.csv").set_index("test")
    rows.append(dict(kind="departure", unit="clock departure over matched (all organs)",
                     departure_log_ratio_vs_matched=float(t294t.loc["H294.1", "median_log_ratio"]),
                     departure_lo90=float(t294t.loc["H294.1", "lo90"]), departure_hi90=float(t294t.loc["H294.1", "hi90"]),
                     departure_fold=float(np.exp(t294t.loc["H294.1", "median_log_ratio"])),
                     entry="A294b H294.1", source_csv=rel(EX / "A294_floor/output/tests.csv")))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t = table()
    t.to_csv(OUT / "two_axis_table.csv", index=False, float_format="%.4g")
    p = polar()
    p.to_csv(OUT / "polar_points.csv", index=False, float_format="%.4g")
    pd.set_option("display.width", 250)
    print(t[["condition", "time_kind", "time", "time_lo90", "time_hi90", "amplitude_ratio", "amplitude_lo90",
             "amplitude_hi90", "Delta", "Delta_loudness_share"]].round(3).to_string())
    print(p.drop(columns=["source_csv", "entry", "note"], errors="ignore").round(3).to_string())


if __name__ == "__main__":
    main()
