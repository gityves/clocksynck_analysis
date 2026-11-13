# A256 design addendum (draft, 2026-11-13; no expression value of these datasets read)

Prepared by `scripts/prep_a256.py` (outputs in `output/prep/`: `<GSE>_expr.parquet`, `<GSE>_sheet.csv`, `design.csv`,
`prep_summary.csv`). Data: `$SCCC_DATA/repeat_biopsy/<GSE>/` (each with PROVENANCE.md), platform tables in `.../platforms/`.
Unit = person x arm (x bout / session). Elapsed hours are from the unit's first biopsy.

## Retained datasets (>= 8 persons with a first and a later same-day biopsy, per-sample matrix public, 12-14 of 14 clock genes)

| dataset | tissue / platform | arm (class) | what happened between the biopsies | elapsed h | t0 clock h (source) | units |
|---|---|---|---|---|---|---|
| GSE40589 | subcutaneous fat, Agilent GPL6848 (log ratios) | INF_NACL (**nothing_done**) | saline 5 mL/h for 240 min, supine, fasted | 4.0 | 9.0 assumed | 7 |
| | | EU_NACL (intervention) | euglycemic hyperinsulinemic clamp + saline | 4.0 | | 8 |
| | | HC_NACL (intervention) | hyperglycemic (140 mg/dl) hyperinsulinemic clamp + saline, diazoxide | 4.0 | | 8 |
| | | INF_GIP / EU_GIP / HC_GIP (intervention) | GIP 2 pmol/kg/min for 240 min, alone / with EU / with HC clamp | 4.0 | | 7 / 8 / 8 |
| GSE9105 | vastus lateralis, HG-U133A (12 of 14 clock genes: PER1, CIART not on the platform annotation) | insulin_clamp (intervention) | 4-h euglycemic hyperinsulinemic clamp, 80 mU/m2/min | 1.5, 5.0 | 9.0 assumed | 12 |
| GSE59088 | vastus lateralis, HuGene 1.0 ST | control (**nothing_done**) | rested 08:30-10:30 (no exercise), then fasting rest; first biopsy here = 2.5 h time point | 2.5 | **13.0 stated (derived)** | 6 |
| | | endurance (intervention) | 2 h cycling ended 10:30, i.e. BEFORE the first biopsy; fasting rest between biopsies | 2.5 | 13.0 | 6 |
| | | resistance (intervention) | 30 min leg resistance exercise ended 10:30, BEFORE the first biopsy; fasting rest between | 2.5 | 13.0 | 6 |
| GSE220899 | vastus lateralis, Clariom S (GPL24324) | 30%1RM and 80%1RM bouts (intervention) | 4 sets squat + 4 sets leg extension to failure; one week apart, each with its own PRE | 4.0, 7.0 (bout of 1 h assumed + 3 / 6 h) | 9.0 assumed | 11 / 10 per bout (22 / 20) |
| GSE202295 | vastus lateralis, RNA-seq counts | NGT, T2D (intervention) | 30 min cycling at 85% HRmax BEFORE the first (post) biopsy; 3 h rest; contralateral leg | 3.0 | 9.0 assumed | 16 / 19 |
| GSE163356 | vastus lateralis, RNA-seq (non-integer deposited counts) | session1, session9 (intervention) | HIIT bout 10 x 4 min at 91% HRmax; pre and 3 h post | 4.0 (bout of 1 h assumed + 3 h) | 9.0 assumed | 10 / 9 |
| GSE126296 | vastus lateralis, HuGene 1.0 ST | sprint (intervention) | 3 x 30 s sprints, 20 min rest between; contralateral leg | 3.67 (220 min, the sample title) | 9.0 assumed | 14 |

Pooled nothing_done units: GSE40589 INF_NACL 7 + GSE59088 control 6 = 13 (each dataset alone < 10).
H256.4 (arms in the same study): GSE40589 (INF_NACL vs the five intervention arms, same persons) and GSE59088 (control vs two exercise arms).

### Quotes the classification and clock times rest on
- GSE40589 (Diabetologia 2015, PMID 25994074; GEO design text): "All studies were done in the morning in the overnight fasted state (>10h since last meal)."
  "Either the participants received only a GIP- or a saline- infusion (0.9% NaCl-isotonic solution) for 240 min." "From each participant two adipose tissue biopsies were taken before and after infusions (t = 0, 240 min), from contralateral sites." "The participants were kept in a supine position throughout the experiments." Crossover: "At different investigation days participants underwent euglycemic (EU)- and hyperglycemic (HC), hyperinsulinemic clamps ... Between examination days an intermission time of at least 7 days was maintained."
- GSE9105 (PMC3581328): "After a 10-h overnight fast, all subjects reported to the General Clinical Research Center (GCRC) at 0800 ... Biopsies of the vastus lateralis muscle were taken at -60, 30, and 240 min of the insulin infusion." (GEO titles say 0, 30, 240 min; the paper's offsets give 1.5 and 5.0 h from the first biopsy.) The paper's 3-person saline control and 5 validation subjects are not in GEO (36 samples = 12 x 3).
- GSE59088 (Sci Data, PMC4432635): "Subjects arrived after an overnight fast at 8 am on day 1 ... Non-exercise control subjects rested from 8.30 am to 10.30 am. ... After exercise, all subjects rested under fasting conditions (with water offered ad libitum) until 3.30 pm (corresponding to 5 h post exercise)." "The non-exercise control group was included, to control for potential non-exercise related effects of fasting conditions and the invasive procedures." The GEO 'Pre' array is the biopsy "at least three day before commencing this training phase" (>= 10 weeks earlier): NOT same-day, dropped. The used pair is 2.5 h vs 5 h post (13:00 and 15:30).
- GSE220899 (PMC9856538): "participants completed a resistance exercise bout in the morning hours (between 7:00 AM-10:00 AM) following an overnight fast ... 4 sets ... provided 5 min of rest between sets." Biopsies "prior to ... 3-hours following ... 6-hours following each bout"; bouts one week apart. Bout length is not stated (8 sets, 5 min rests: about 1 h assumed).
- GSE202295 (PMC9451165): "Participants returned to the clinic 5 to 10 days after the baseline measurement and performed an acute bout of exercise ... Directly after the exercise bout, additional ... skeletal muscle biopsies were taken. The participants then rested for 3 hours before the last blood samples and biopsy were taken from the contralateral leg." The 'rest' (basal) biopsy is from the baseline visit: not same-day, dropped. (GEO text gives both "30-minute ... 85% of maximal heart rate" and "35 min ... 80% heart rate reserve".)
- GSE163356 (J Appl Physiol 2022; full text not open, GEO + abstract only): "Before and 3 h after the 1st and 9th exercise bout, skeletal muscle biopsies were obtained" (10 x 4 min at 91% HRmax). No clock time, no bout duration stated.
- GSE126296 (PMC6812755): "The subjects reported to the laboratory in the morning after an overnight fast ... three all-out cycle sprints of 30-s duration ... with 20 min of rest between sprints ... The first biopsy sample was obtained randomly in either right or left leg before the first sprint. The second sample was obtained 2 h and 20 min after the third sprint from the contralateral leg." GEO titles: pre 0 min, post 220 min.

### Classification rule applied
nothing_done = no exercise, meal, insulin, glucose or drug between the biopsies (saline infusion, rest). Arms in which exercise
preceded the FIRST biopsy (GSE59088 endurance / resistance, GSE202295) are labelled intervention (post-exercise recovery, not
a clean rest interval); the lead may reclassify them. All intervals are fasted (limit already in A256).

## Rejected
| dataset | reason |
|---|---|
| GSE7146 | two protocols of 6 persons each (2 h clamp on GPL96, 3 h on GPL80), < 8 per protocol and platform, person overlap unknown; insulin clamp only |
| GSE164081 | no per-sample matrix (only TPM per time point pooled over persons; three series on three platforms) |
| GSE43219 | rest biopsy on a separate day ("a day separated from exercise tests by at least one week"), not a same-day pair |
