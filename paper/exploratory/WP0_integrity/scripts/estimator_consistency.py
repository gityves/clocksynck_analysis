"""WP0 item 4: estimator behind every host-referenced quantity in the Figures section (plan "## Figures" to the end).

Curated table (estimator read from the producing migrated script, evidence file:line relative to paper/results, and from
the archive amendment text).  flag = ok | n/a (no coupling) | FLAG: reason.
"""
import pandas as pd
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "output"
H = ["figure", "panel", "amendment", "estimator", "anchors", "a22_confirmed", "flag", "evidence"]
A17 = "anchored A17 (leave-one-out host)"
rows = [
 ("1", "Benchmark vs COFE/CYCLOPS/CHIRAL/TimeTeller", "benchmark", "host-referenced amplitude/offset; no coupling", "GTEx waveform", "n/a", "n/a (no coupling)", "figure_1/scripts/method/compute_gtex_benchmark.py"),
 ("1", "Timed truths: regenerating mouse colon (9 h shift)", "A21", A17, "stromal-vascular rule present in colon", "no: V17.2 fails, injured anchor (lymphatic EC) decouples",
  "FLAG: A21 R21.2 coupling contrast attenuated; an anchor decouples in R-DC (A22 V17.2 fail). Phase (R21.4) passes", "figure_1/scripts/species/compute_colon_injury.py:66-67"),
 ("1", "Timed truths: baboon, mouse three ages, living skin", "A36, A48, A37", "timed host phase / offsets; no coupling", "peripheral tissues", "n/a", "n/a (no coupling)", "figure_4/scripts/resilience/compute_mouse_ageing.py:125"),
 ("2", "Donors' shared departure from time of death; organs follow with unequal gains", "01_gtex_organs (pre-A1)", "timed per-tissue departure, factor analysis; not the A17 estimator", "none (timed)", "n/a",
  "FLAG: pre-A17 analysis; gains and departures are not anchored couplings (A103 r* builds on them)", "figure_2/scripts/organs/compute_organs_timed.py:1-8"),
 ("2", "Brain vs periphery r* 0.64", "A103", "two block hosts, every tissue of the block an anchor; disattenuated phase correlation", "all tissues of each block", "n/a",
  "FLAG: not the A17 anchor rule; a host-correlation statistic, not a coupling", "figure_2/scripts/organs/compute_brain_periphery_species.py:6,150-151"),
 ("2", "Sheep r* 0.63; brain couples less in sheep and baboon", "A42", "r*: block hosts (as A103). Coupling readout 2: peripheral-tissue host, brain read on it", "all peripheral tissues", "n/a",
  "FLAG: organ-level host, not the A17 rule; r* readout 1 fails in baboon", "figure_2/scripts/organs/compute_brain_periphery_species.py:199-200"),
 ("2", "Anchored coupling by cell family across healthy atlases", "A17 unified", A17 + ", learned and fixed rules, frozen", "learned / A17 frozen / fixed stromal-neuronal", "n/a (healthy atlases; A22 covers disease contrasts)", "ok",
  "figure_2/scripts/atlas_fits/compute_atlas_fits.py:5,65"),
 ("2", "Healthy mouse colon against known time", "A21 R21.1", A17, "stromal-vascular rule", "n/a", "ok (Spearman 0.95 vs timed coupling)", "figure_1/scripts/species/compute_colon_injury.py:66-67,88"),
 ("2", "Immune cells lean on body-wide cues (index 0.18 vs 0.10)", "A85 corrected", "c_k = a H_-k + b Q_-k on the mean of the other classes (reference-free); not anchored", "none (leave-class-out mean)", "n/a",
  "FLAG: not the A17 estimator; index is 0.5 when neither regressor carries signal (archive: uninformative before correction)", "figure_5/scripts/cancer/compute_cue_following.py:5-8"),
 ("2", "Strong clocks are the coupled ones (rho +0.73)", "A104", "c on a leave-one-out host of all peripheral tissues; own amplitude A", "all peripheral tissues", "n/a",
  "FLAG: organ-level host, not the A17 rule; registered prediction rho<0 not supported", "figure_2/scripts/backbone/compute_oscillator_gtex.py:72"),
 ("2", "Neurons carry the brain host, glia less; retina", "02, 04 (A1, A2), A18 3.5", "02/04: free fit (pre-A17, 2026-10-11); A18 3.5 anchored contrast on learned anchors", "02/04 none; A18 3.5 learned", "A18 3.5: no (not in A22 list)",
  "FLAG: 02/04 panels rest on pre-A17 free-fit couplings; only A18 3.5 is anchored", "figure_2/scripts/backbone/compute_followups_a18.py:56,61; extended_data/scripts/host_identity/compute_cortex_modules.py:51"),
 ("2", "Mouse embryo shared clock", "A108-A108c", "gate and timed phase R; no coupling", "n/a", "n/a", "n/a (no coupling)", "figure_2/scripts/embryo/compute_embryo_timed.py"),
 ("2", "Age: prefrontal amplitude loss; peripheral phase advance", "A3, A94", "timed amplitude (A3); GTEx radius/phase (A94); no coupling", "n/a", "n/a", "n/a (no coupling)", "figure_4/scripts/resilience/compute_gtex_age.py"),
 ("3", "Fast, intermediate, slow, ventilated deaths: inter-organ coherence", "A92, A93, A109", "pair coherence C of clock states, tissues as anchors; not the anchored coupling c", "all reference tissues", "n/a",
  "FLAG: coherence has no clock-specific null (A109 adds it); not the A17 coupling", "figure_3/scripts/cues_lost/compute_ventilator_tissues.py:4,139"),
 ("3", "A step, not a drift: ventilation days", "A46, A92, A95", "radius / wall-time R / coherence on all-tissue anchors; no coupling", "all tissues", "n/a", "n/a (no coupling)", "figure_3/scripts/gtex_cues/compute_ventilation_dose.py:147"),
 ("3", "Liver, muscle, adrenal loosen; brain leaves the body", "A33, A93", "coupling on leave-one-out host of the reference tissues (tissue-level)", "all reference tissues", "n/a",
  "FLAG: tissue-level host, not the A17 rule; A33 tissue-specific, not a class effect", "figure_2/scripts/organs/compute_organs_timed.py:55-56"),
 ("3", "Rhythmic genes follow internal time; drug targets", "A91, A100, A102", "gene-level excess over a null; no coupling", "n/a", "n/a", "n/a (no coupling)", "figure_3/scripts/gtex_cues/compute_internal_external.py"),
 ("4", "Forest plot of all disease contrasts with MDEs (d*)", "A17-A22", A17 + ", learned and fixed rules, contrasts at common spread d*", "learned / A17 frozen / fixed", "yes (22 A20 contrasts; anchors dropped where flagged)",
  "ok; the A17 table predates the A22 refit: ild_natri capillary EC, pd, SCZ_hbcc, BD_hbcc, endometriosis had an anchor dropped", "figure_4/scripts/resilience/compute_contrasts_a17.py:41,47"),
 ("4", "Hierarchical pooling by class and family: coupling d* (y1)", "A20, A22", "pooled y1 = d* of A17 contrasts, inputs on A22-confirmed anchors", "learned, group-confirmed", "yes", "ok (coupling readouts move <= 0.02 after A22)",
  "figure_4/scripts/spread_multiplicity/compute_spread_inputs.py:56,83"),
 ("4", "Hierarchical pooling: own-noise 'damping' (y2, y3)", "A20, A22", "log z2 ratio and log psi ratio on the same anchored host", "learned, group-confirmed", "yes",
  "FLAG: 'damping' withdrawn as a headline in A22; psi contrasts not recovered in A21 (R21.2 0.43 < 0.5), slope no longer excludes 0", "figure_4/scripts/spread_multiplicity/compute_spread_inputs.py:56"),
 ("4", "Paired same-time designs: tumours, inflamed gut, skin", "V52, A52, A53", "paired clock-plane variance, same-time unaffected sample; no coupling", "n/a", "n/a", "n/a (no coupling)", "figure_5/scripts/cancer/compute_tumour_gate.py"),
 ("4", "Ageing, sex, BMI, T2D, ventilation dose", "A46, A47, A94, A99", "A17 anchored age slopes (cell classes); GTEx radius/C elsewhere", "learned (fixed fallback, flagged)", "no (not in A22 list)", "ok for the cell-class age slopes; GTEx measures are not couplings", "figure_4/scripts/resilience/compute_age.py:90-95"),
 ("4", "COVID-19 decedents keep organ clocks and endothelial sharing", "A105, A105b", "organ coupling: every organ an anchor; lung arms: two fixed anchors (endothelial, fibroblast)", "all organs / 2 fixed", "no (not in A22 list)",
  "FLAG: organ-level host; lung contrast has 2 anchors (A17 asks for >= 3) and no A22 confirmation", "figure_4/scripts/natural_experiments/compute_covid_couplings.py:34,43,52"),
 ("4", "Fatty liver; inflamed gut stroma", "A43, A44, A96, A45", "contrast on fixed stromal anchors; A45 learned-rule sensitivity", "fixed rule / learned (A45)", "no (not in A22 list)", "ok (A17-style fixed rule); no A22 group confirmation", "figure_4/scripts/resilience/compute_masld_gut.py:68; compute_a45_learned.py:76"),
 ("5", "Malignant cells: clock-shaped share; decomposition amplitude/coherence", "A86, A87, A97", "gate clock share; relative amplitude |z| and coherence on anchored host (A87)", "learned anchors (cancer)", "no", "ok (A87 on A17 anchors); A86 share is a gate statistic", "figure_5/scripts/cancer/compute_anchored_cancer.py:5-10"),
 ("5", "Psoriatic delay, core clock vs cues, treatment", "A54-A81b", "paired phase offset arg(z_lesion/z_uninvolved); no coupling", "n/a", "n/a", "n/a (no coupling)", "figure_5/scripts/psoriasis/compute_replication.py"),
 ("5", "Fibrotic-lung T cells lose relative amplitude, two cohorts", "A98, A5, A17, A19", "anchored decomposition (|z|, psi, c) on learned anchors; contrast A17", "learned", "yes (A22: ild_natri capillary EC dropped)",
  "FLAG: original A5 coupling (Natri 0.45, HLCA 0.14) was free-fit; confirm the A17 contrast table equals the A22 refit; procurement (A19) unresolved", "figure_5/scripts/lung_breaks/compute_decomposition.py:106-118"),
 ("5", "Per-cell synchrony vs coherence in malignant cells", "WP0 item 2", "per-cell synchrony vs A86 coherence", "n/a", "n/a", "n/a (no coupling)", "figure_2/scripts/backbone/compute_per_cell.py"),
 ("ED", "Estimator validation: A17 simulations, small-n calibration, anchors, A22", "A17, A22", A17, "learned / fixed", "yes", "ok", "extended_data/scripts/calibration/compute_anchored_simulations.py"),
 ("ED", "Early immune couplings near 1 (A4, A5 free fits)", "A4, A5, A6", "free fit, small-sample inflation (lymphocytes 0.88 -> 0.35 corrected -> 0.09 anchored)", "none", "n/a",
  "FLAG: if cited anywhere it is pre-A17; no figure-table panel cites it", "archive RESULTS_preliminary.md:1026-1027"),
]
df = pd.DataFrame(rows, columns=H)
df.to_csv(OUT / "estimator_consistency.csv", index=False)
fl = df[df.flag.str.startswith("FLAG")]
print(len(df), "rows;", len(fl), "flagged;", (df.flag.str.startswith("ok")).sum(), "ok;", df.flag.str.startswith("n/a").sum(), "n/a")
