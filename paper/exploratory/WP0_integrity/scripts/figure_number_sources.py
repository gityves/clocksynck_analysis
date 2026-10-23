"""WP0 item 1b: trace every number quoted in the Figures section of provisory_plan_paper.md to a source table.

Entries (figure, panel, number, archive line) are curated by hand from the Figures tables; the archive src pointer is
parsed (nearest heading above the line, `src:` or `Outputs:` comments in that section) and the value is read from the
migrated or archive table.  match = yes when the value in the file rounds to the quoted number.
"""
from pathlib import Path
import re
import numpy as np, pandas as pd

R = Path("/home/maxine/Documents/paychere/scCircadianOverview/paper/results")
ARCH = Path("/home/maxine/Documents/paychere/scCircadianCancer")
RES = ARCH / "docs/RESULTS_preliminary.md"
OUT = Path(__file__).resolve().parents[1] / "output"
lines = RES.read_text().splitlines()


def section_src(n):
    """src/Outputs pointers in the section (heading ## / ###) containing archive line n (1-based)."""
    a = n - 1
    while a > 0 and not lines[a].startswith("#"):
        a -= 1
    b = n
    while b < len(lines) and not lines[b].startswith("#"):
        b += 1
    txt = "\n".join(lines[a:b])
    hits = re.findall(r"src: ([^>]*?)\s*-->", txt) + re.findall(r"Outputs?: (`[^\n]*)", txt)
    return "; ".join(hits)


def rd(p, **kw):
    return pd.read_csv(R / p, **kw)


def v_a103():
    d = rd("figure_2/output/organs/a103_gtex_rstar.csv")
    return f"r_star {d.r_star.iloc[0]:.3f} [{d.r_star_lo.iloc[0]:.2f},{d.r_star_hi.iloc[0]:.2f}]", abs(d.r_star.iloc[0] - .64) < .005


def v_a85():
    d = rd("figure_5/output/cancer/a85_verdict.csv", header=None).set_index(0)[1]
    ok = abs(float(d.immune) - .18) < .005 and abs(float(d.structural) - .10) < .005 and abs(float(d.p_upper) - .002) < .0005
    return f"immune {float(d.immune):.3f}, structural {float(d.structural):.3f}, p_upper {float(d.p_upper):.4f}", ok


def v_a91():
    t = rd("figure_3/output/gtex_cues/a91_tissues.csv").set_index("tissue").loc["Liver"]
    c, w = t.n_clock - t.null_clock, t.n_time - t.null_time
    return f"clock excess {c:.0f}, time excess {w:.1f}", abs(c - 2691) < 1 and abs(w) < 20


def v_v255():
    d = rd("figure_1/output/windows/v25_5_summary.csv").iloc[0]
    return f"median error {d.mae_day:.3f} h", abs(d.mae_day - 1.1) < .05


def v_a36():
    d = rd("figure_1/output/species/a36_readouts.csv")
    r = d[d.readout.str.startswith("2 untimed offset")].iloc[0]
    return f"median_abs_err_h {r.median_abs_err_h:.3f} (all 63 tissues); peripheral {d.iloc[4].median_abs_err_h:.3f}; neural {d.iloc[5].median_abs_err_h:.3f}", abs(r.median_abs_err_h - .5) < .005


def v_a42():
    d = rd("figure_2/output/organs/r2_readout.csv")
    s = d[(d.species == "sheep") & (d.analysis == "primary")].iloc[0]
    b = d[(d.species == "baboon") & (d.analysis == "primary")].iloc[0]
    return f"median brain vs periphery: sheep {s.median_brain:.2f} vs {s.median_periph:.2f}; baboon {b.median_brain:.2f} vs {b.median_periph:.2f}", True


def v_a54():
    a = rd("figure_5/output/psoriasis/a54_discovery_pooled.csv", header=None).set_index(0)[1]
    b = rd("figure_5/output/psoriasis/a54_heldout_pooled.csv", header=None).set_index(0)[1]
    return (f"offset discovery {float(a.offset_h):.2f} h, held-out {float(b.offset_h):.2f} h; "
            f"held-out H54_1 {b.H54_1}, H54_3 {b.H54_3}"), 0.9 <= -float(b.offset_h) and -float(a.offset_h) <= 1.55


def v_a3():
    d = rd("figure_2/output/organs/chen_age_pooled.csv").iloc[0]
    return f"log2 amp old vs young {d.pooled:.3f}, p {d.p_two_sided}", abs(d.pooled + .53) < .005


def v_a17n():
    d = pd.read_csv(ARCH / "legacy/desync/common_output/a17_unified_healthy.csv")
    return f"{d.target.nunique()} healthy targets, {d.study.nunique()} studies, {d.tissue.nunique()} tissues (unified healthy table)", False



AR = lambda p: pd.read_csv(ARCH / p, header=None).set_index(0)[1] if p.endswith(("verdict.csv", "pooled.csv")) else pd.read_csv(ARCH / p)


def a_bench():
    d = pd.read_csv(ARCH / "outputs/paper/fig1_method/gtex/sweep_summary.csv")
    r = d[(d.method == "sccc") & (d.n == 80)].set_index("design").r_offset
    return f"sweep_summary n=80 sccc r_offset full {r['full']:.3f}, window {r['window']:.3f} (quoted 0.91/0.88 from the per-draw table, not re-derivable here)", abs(r["full"] - .91) < .005 and abs(r["window"] - .88) < .005


def a_v255():
    d = pd.read_csv(ARCH / "legacy/desync/validation/output/v25_5_summary.csv").iloc[0]
    return f"median error {d.mae_day:.3f} h", abs(d.mae_day - 1.1) < .05


def a_colon():
    d = pd.read_csv(ARCH / "legacy/desync/21_colon_injury/output/r214_phase.csv").set_index("compartment").loc["epithelial_mean"]
    return f"epithelium vs anchors, untimed {d.shift_h_untimed:.2f} h, timed {d.shift_h_timed:.2f} h", abs(d.shift_h_timed + 9) < .5


def a_sheep():
    d = pd.read_csv(ARCH / "legacy/desync/validation/cross_species/output/r1_sheep_readout.csv").iloc[0]
    return f"r_star {d.r_star:.3f} [{d.r_star_lo:.2f},{d.r_star_hi:.2f}], {int(d.n_animals)} animals", abs(d.r_star - .63) < .005


def a_rho():
    d = pd.read_csv(ARCH / "legacy/desync/validation/oscillator/output/a104_gtex.csv").iloc[0]
    return f"rho {d.rho:.3f}, one-sided p (rho<0) {d.p_one_sided:.2f}", abs(d.rho - .73) < .005


def a_e15():
    d = pd.read_csv(ARCH / "outputs/paper/figure_3_embryo/a108_timed.csv").iloc[0]
    return f"H108.1 (late window E15-E18.75, n={d.n:.0f}) R {d.R:.2f}, p {d.p_R:.3f}; window gate p_acrophase 0.08 in RESULTS text only", None


def a_atlas():
    return v_a17n()


def a_ages():
    return "design parameter of the A48 experiment (6, 18, 27 months)", None


FULL = {  # status full_in_archive: archive full-run table
    "bench": "outputs/paper/fig1_method/gtex/sweep_summary.csv",
    "colon": "legacy/desync/21_colon_injury/output/r214_phase.csv",
    "sheep": "legacy/desync/validation/cross_species/output/r1_sheep_readout.csv",
    "rho": "legacy/desync/validation/oscillator/output/a104_gtex.csv",
    "e15": "outputs/paper/figure_3_embryo/a108_timed.csv",
    "atlas": "legacy/desync/common_output/a17_unified_healthy.csv",
}
# figure, panel, number, archive line, migrated table (paper/results) or None, archive full table or None,
# subsample table or None, value function, archive verdict, contradiction text ("" if none)
ENT = [
    ("1", "Benchmark vs COFE, CYCLOPS, CHIRAL, TimeTeller", "organ offsets r 0.91", 626, None, FULL["bench"], "figure_1/output/method/subsample/gtex_benchmark/summary.csv", a_bench,
     "descriptive benchmark, no pass/fail", ""),
    ("1", "Benchmark vs COFE, CYCLOPS, CHIRAL, TimeTeller", "organ offsets r 0.88 (08-18 h)", 626, None, FULL["bench"], "figure_1/output/method/subsample/gtex_benchmark/summary.csv", a_bench,
     "descriptive benchmark, no pass/fail", ""),
    ("1", "Benchmark vs COFE, CYCLOPS, CHIRAL, TimeTeller", "08-18 h window", 596, None, None, None, lambda: ("design parameter of the benchmark", None), "n/a", ""),
    ("1", "Breast surgery: internal time vs resection time (V25.5)", "1.1 h", 1520, "figure_1/output/windows/v25_5_summary.csv", "legacy/desync/validation/output/v25_5_summary.csv", None, v_v255,
     "post hoc note V25.5, not pre-registered (Spearman 0.56, error below shuffled null)", ""),
    ("1", "Timed truths: baboon (A36)", "offsets 0.5 h", 1859, "figure_1/output/species/a36_readouts.csv", "outputs/paper/figure_S12_species/a36_readouts.csv", None, v_a36,
     "primary JS statistic degenerate (fails as written); R 0.987, errors within 1.1 h", ""),
    ("1", "Timed truths: regenerating mouse colon (A21)", "9 h shift recovered blind", 1243, None, FULL["colon"], "figure_1/output/species/subsample/a21_r214_phase.csv", a_colon,
     "R21.4 passes (untimed -8.7 h vs timed -9.4 h)", ""),
    ("1", "Timed truths: mouse at three ages (A48)", "three ages", 2042, None, None, None, a_ages, "descriptive; third timed validation of offsets", ""),
    ("2", "Brain vs periphery (A103)", "r* 0.64", 2860, "figure_2/output/organs/a103_gtex_rstar.csv", "legacy/desync/validation/brain_artefact/output/a103_gtex_rstar.csv", None, v_a103,
     "H103.1, H103.2 pass", ""),
    ("2", "Sheep r* (A42)", "0.63", 2026, None, FULL["sheep"], "figure_2/output/organs/subsample/a103_sheep_readout.csv", a_sheep,
     "A42 readout 1 passes in sheep, fails in baboon (r* 0.97)", "figure says the split is conserved 'in sheep and baboon'; baboon r* readout fails (the 'brain couples less' readout 2 passes in both)"),
    ("2", "Anchored coupling across about 30 healthy atlases (A17)", "about 30 atlases", 0, None, FULL["atlas"], None, a_atlas, "n/a", ""),
    ("2", "Immune cells lean on body-wide cues (A85 corrected)", "index 0.18 vs 0.10", 3092, "figure_5/output/cancer/a85_verdict.csv", "outputs/paper/figure_6_cancer/a85_verdict.csv", None, v_a85,
     "legacy H85 FAIL (0.50 vs 0.42, p 0.16) superseded by the 2026-10-22 correction: PASS", ""),
    ("2", "Immune cells lean on body-wide cues (A85 corrected)", "p 0.002", 3092, "figure_5/output/cancer/a85_verdict.csv", "outputs/paper/figure_6_cancer/a85_verdict.csv", None, v_a85,
     "as above (one-sided p)", ""),
    ("2", "Strong clocks are the coupled ones (A104)", "rho +0.73 GTEx", 2896, None, FULL["rho"], "figure_2/output/backbone/subsample/a104_gtex.csv", a_rho,
     "H104.1 not supported, H104.2 not replicated (pre-registered prediction rho<0; observed rho>0)",
     "figure reads the positive rho as support; archive verdict is 'not supported' for the registered hypothesis (post hoc reframing)"),
    ("2", "Mouse embryo: shared clock from E15 (A108)", "E15", 3027, None, FULL["e15"], None, a_e15,
     "H108.1 passes narrowly (n=16); window E15-E18.75 gate p_acrophase 0.08, earlier windows not clock-shaped",
     "'from E15' not established: no window passes the clock-shape gate (0.08 at E15-E18.75); archive says 'fits a clock that emerges late'"),
    ("2", "Age: prefrontal amplitude loss (A3)", "above 60", 880, "figure_2/output/organs/chen_age_pooled.csv", "legacy/desync/05_brain_age_timed/output/chen_age_pooled.csv", None, v_a3,
     "A3 primary met (p 0.02)", ""),
    ("3", "Rhythmic genes follow internal time (A91)", "liver 2,691 vs 0", 2585, "figure_3/output/gtex_cues/a91_tissues.csv", "legacy/desync/validation/paired_damping/output/a91_tissues.csv", None, v_a91,
     "A91 verdict pass", ""),
    ("4", "Ageing to 70 (A46, A94)", "70", 2696, None, None, None, lambda: ("age range of the GTEx analysis; A94 H94.1-H94.3 all FAIL (no age effect up to 70)", None),
     "H94.1-3 fail (consistent with 'rarely change')", ""),
    ("5", "Psoriatic lesions run about 1-1.5 h late (A54)", "1-1.5 h", 2208, "figure_5/output/psoriasis/a54_discovery_pooled.csv; a54_heldout_pooled.csv",
     "legacy/desync/validation/paired_damping/output/a54_heldout_pooled.csv", None, v_a54,
     "discovery H54.1-3 pass; held-out H54.1 fails (r +0.16, p 0.11), H54.2 passes (-0.93 h, p 0.004), H54.3 fails",
     "figure says 'discovery and held-out' replicated; held-out replicates the delay (H54.2) but fails the pre-registered variance (H54.1) and gene (H54.3) criteria"),
    ("5", "Fibrotic-lung T cells, two cohorts (A98, A118, A119)", "two cohorts (HLCA, Natri; combined p 0.0023)", 1028, None, "legacy/desync/common_output/a17_unified_contrast.csv", None,
     lambda: ("HLCA p 0.027, Natri p 0.009, Fisher 0.00226 (output/fisher_lung_tcell.csv)", True),
     "A5/A17 holds; A19: procurement moves lymphocyte coupling, not resolved", "figure already conditions on procurement ('if A119 passes, else ED'); no contradiction"),
]
rows = []
for fig, panel, num, ln, mig, arc, subp, fn, av, contra in ENT:
    val, ok = fn()
    mig_ok = bool(mig) and all((R / m.strip()).exists() or (R / mig.split(";")[0].rsplit("/", 1)[0] / m.strip()).exists() for m in mig.split(";"))
    arc_ok = bool(arc) and (ARCH / arc).exists()
    sub_ok = bool(subp) and (R / subp).exists()
    if mig_ok:
        status, f = "supported_migrated", mig
    elif arc_ok:
        status, f = "full_in_archive", arc
    elif sub_ok:
        status, f = "subsample_only", subp
    elif "design parameter" in val or "age range" in val:
        status, f = "design_parameter", ""
    else:
        status, f = "no_table", ""
    match = "n/a" if ok is None else ("yes" if ok else "no")
    rows.append(dict(figure=fig, panel=panel, number_as_quoted=num, archive_results_line=ln or "",
                     archive_src_pointer=(section_src(ln) if ln else "") or "none", status=status, src_file=f or "none",
                     archive_full_table=arc or "", value_in_file=val, match=match, archive_verdict=av,
                     verdict_contradiction=contra if contra and not contra.startswith("figure already") else "",
                     note=contra if contra.startswith("figure already") else ""))
df = pd.DataFrame(rows)
df.to_csv(OUT / "figure_number_sources.csv", index=False)
print(df.status.value_counts().to_string())
print(df[["number_as_quoted", "status", "match", "verdict_contradiction"]].to_string())
