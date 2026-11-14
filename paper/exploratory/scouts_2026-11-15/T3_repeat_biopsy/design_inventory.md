# T3 design inventory: same-person, same-tissue, same-day biopsies (read-only; no expression statistic computed)

Sources: paper/exploratory/scouts_2026-11-12/repeat_biopsy_scout.md (S), A252_elapsed_time/output/design.csv (D252),
A256_elapsed_replication/output/design.csv (D256) and design_addendum_draft.md (ADD), registered_experiments.md A252/A256/A257,
series-matrix headers (GSE320593, GSE106800). "units" = person x arm/visit/leg/bout, as in the A252/A256 design tables.
"assumed" = clock time of the first biopsy is a protocol guess registered before scoring (not in the data). Inferred items are marked (inf).

## Studies in repeat_biopsy/ (18 folders; all 18 were used in A252 or A256)

| GSE | tissue; platform; data | persons with >= 2 same-day biopsies | clock time (t0) / elapsed h | between the biopsies, per arm (persons) | use | caveats |
|---|---|---|---|---|---|---|
| GSE118280 | subcut. adipose; Affy HuGene 2.0 ST; RMA log2 only (RAW not downloaded); genome-wide (25582 genes; 13 of 14 clock genes, CIART absent) | 15 men; 30 person-diet visits x 3 biopsies | 08:40 (stated), 12:20, 19:00 = 0, 3.67, 10.33 h; the only study with stated clock times | meals; HC/HF vs HF/HC 4-week diets, crossover, two visits (15 persons, 30 units) | A252 primary (H252.1) | meals between biopsies in all arms; 4-week diet is a between-visit perturbation; men only |
| GSE252357 | vastus lateralis; normalised counts, protein-coding only (14392 genes; 14 clock genes) | 13 (S) | first biopsy 1 h after breakfast, overnight fast; clock time not given, 09:00 assumed; 0, 1.75, 4.25, 9.25 h (24 h point = 25.25 h not used) | CTRL: 45 min rest, time-matched biopsies (5 persons; breakfast, and lunch after the 3 h biopsy, are standard in all arms); RE: 45 min leg press/extension (8 persons) | A252 pooled; rest arm descriptive (rho 2.89 on 5 persons, per A256 text) | only explicit no-intervention arm among counts; meals not withheld; 9.25 h point has 10 samples (A252 score table: 10 units) |
| GSE142401 | adipose; Illumina HT-12 v4, quantile-normalised log2 (20760 genes; 14) | 19 men (10 control, 9 MetS), 2 meals each = 38 units | 0 h fasting and +4 h; clock time not stated, 09:00 assumed | high-fat breakfast (dairy or soy), crossover on different days | A252 pooled | meal in every interval; men only |
| GSE231509 | skeletal muscle; raw counts (59653 genes; 14) | 21 (7 healthy, 7 obese, 7 obese+T2D) | arrival 09:00 after 12 h fast (biopsy time not given), 09:30 assumed; 0 and 1 h | mixed-meal tolerance test | A252 pooled | 1 h gap only (15 deg predicted); ids repeat across groups |
| GSE151066 | muscle; RSEM expected counts (60920; 14) | 19 older adults (10 active, 9 sedentary) | pre, immediately post, +3 h post = 0, 0.67, 3.67 h; 09:00 assumed | 40 min cycling | A252 pooled | exercise in interval |
| GSE273041 | muscle; non-integer normalised counts (25262; 14) | 20 athletes (10 high-carb, 10 ketogenic) | BL, after 3 h run, +2 h = 0, ~3, ~5 h; 09:00 assumed | 3 h treadmill run | A252 pooled | long exercise; diet groups |
| GSE303983 | abdominal subcut. adipose; raw gene counts (61228; 14) | 45 | pre, 1.5 h after exercise end = 2.25 h from first; 09:00 assumed | one exercise session (LOW 60 min / MOD 45 min / HIGH 10x1 min); arm per person not in GEO metadata | A252 pooled | elapsed time approximate (exercise duration unknown per person) |
| GSE120862 | muscle; DESeq-normalised counts (56609; 14) | 7 men; 28 leg x training units | baseline, +1 h, +4 h after exercise (~0, 2, 5 h); 09:00 assumed | 60 min one-leg knee extension; other leg = within-person no-exercise control; repeated after 2 months training | A252 pooled | n=7 only; baseline may be 48 h after previous exercise (S); control leg is a no-exercise leg but still fasted/meal state shared |
| GSE4386 | right atrium; HG-U133 Plus 2 RMA (22836; 13) | 20 (sevoflurane 10, propofol 10) | start vs end of off-pump CABG; no clock time, elapsed NA | surgery + anaesthesia | A252 descriptive only | not usable for T3 (time unknown) |
| GSE75432 | muscle; raw counts (64258; 14) | 10 | before vs >= 2 h after total knee arthroplasty with tourniquet; no clock time | surgery + ischaemia | A252 descriptive only | not usable for T3 |
| GSE151648 | liver graft; raw counts (60498; 14) | 40 (IRI+ 23, IRI- 17) | pre- vs post-transplant; elapsed NA | cold ischaemia, reperfusion | A252 descriptive only | not same body, not usable for T3 |
| GSE40589 | subcut. adipose, contralateral sites; Agilent GPL6848 log ratios | 52 persons-arms: INF_NACL 7, EU_NACL 8, HC_NACL 8, INF_GIP 7, EU_GIP 8, HC_GIP 8 (units; persons overlap across crossover days) | 0 and 240 min = 4.0 h; morning, overnight fasted (>10 h); 09:00 assumed | INF_NACL: saline 5 mL/h x 240 min, supine, fasted = nothing done (7); NACL + euglycaemic clamp (8); NACL + hyperglycaemic clamp with diazoxide (8); GIP alone (7); GIP + EU (8); GIP + HC (8) | A256 (replication) | crossover with >= 7 days between visits; log ratios, not counts |
| GSE9105 | vastus lateralis; HG-U133A (12 of 14 clock genes) | 12 | biopsies at -60, 30, 240 min of insulin infusion (paper) = 0, 1.5, 5.0 h; GCRC arrival 08:00 after 10 h fast; 09:00 assumed | 4 h euglycaemic hyperinsulinaemic clamp (12) | A256 | paper's 3-person saline control not in GEO (36 samples = 12 x 3) |
| GSE59088 | vastus lateralis; HuGene 1.0 ST (14) | 18 (6 per arm) | used pair: 2.5 h vs 5 h post (stated 13:00 and 15:30) = 2.5 h apart; t0 13.0 stated (derived) | control: rested 08:30-10:30, no exercise, fasting rest, "to control for non-exercise effects of fasting and invasive procedures" (6; nothing done); endurance: 2 h cycling ended 10:30, before the first biopsy (6); resistance: 30 min leg exercise ended 10:30 (6) | A256 | GEO "Pre" array is >= 10 weeks earlier (dropped). Best fasting rest design but 6 persons |
| GSE220899 | vastus lateralis; Clariom S | 11 and 10 per bout (22 / 20 units) | PRE, +3 h, +6 h after bout, elapsed 4.0 and 7.0 h (bout ~1 h assumed); 09:00 assumed (exercise 07:00-10:00) | resistance bouts 30%1RM and 80%1RM, one week apart, own PRE each | A256 | bout length not stated |
| GSE202295 | vastus lateralis; RNA-seq counts | NGT 16, T2D 19 | post-exercise biopsy vs 3 h recovery (contralateral leg) = 3.0 h; 09:00 assumed | 30 min cycling before the first (post) biopsy; 3 h rest | A256 | resting basal biopsy is from another day (dropped) |
| GSE163356 | vastus lateralis; non-integer deposited counts | 10 (session 1), 9 (session 9) | pre and +3 h post = 4.0 h (bout 1 h assumed); 09:00 assumed | HIIT 10 x 4 min at 91% HRmax | A256 | full text not open; no clock time |
| GSE126296 | vastus lateralis; HuGene 1.0 ST | 14 | pre and 220 min post = 3.67 h, contralateral leg; 09:00 assumed | 3 x 30 s sprints, 20 min rest between | A256 | |

Rejected for A256 (not in the downloaded set unless listed): GSE7146 (two 6-person clamp protocols), GSE164081 (no per-sample matrix), GSE43219 (rest biopsy another day).
Not downloaded / not usable: GSE260643 (LV pre/post ischaemia, AM vs PM, DEG tables only), GSE108643, GSE305038 (DE tables only), GSE15480, GSE87487 (liver transplant).

## Timed same-person pair outside repeat_biopsy/
- GSE106800 (cancer_circadian/geo/GSE106800; Wefers 2018): vastus lateralis, HuGene 2.1 ST (GPL17692), log2 series matrix. 14 men in the paper summary (series header), 12 subjects with 47 arrays in the sheet (scout), 11 with all four arrays (A257). Biopsies at 07:00 and 19:00 (stated "time of sampling" in GEO), after a 3-day aligned and a 3.5-day misaligned (12 h behaviour shift) protocol, so 12 h apart with stated clock times. "Evening (baseline)" sample is labelled first in the sheet; whether 19:00 and 07:00 fall on the same calendar day is not stated in the header (inf: A257 treats 07:00 -> 19:00 as 12 h elapsed). Not rest/fasting: behaviour and fasting state per protocol; 12 h shift condition is an intervention. Used in A124 (seen) and A257 (not blind). Same-person, exact clock times, the only 12 h span.

## living_multitissue/ (9 folders): no same-tissue same-day repeat biopsy design
- GSE24293 liver, GSE24294 omental fat, GSE24297 subcut. fat, GSE24335 (mixed tissues; MGH RYGB, Agilent two-colour, one sample per tissue per person): different tissues, one time point (inf from PROVENANCE).
- GSE24883 paired SAT/VAT, 32 women; GSE73439 SAT vs VAT same woman (pregnancy/labour); E-TABM-1140 MuTHER (Illumina HT-12, not repeated): different depots, not repeats of one tissue.
- GSE133979 rectus abdominis and SAT, PDAC cachexia vs controls: two tissues, one time.
- GSE320593 (liver and subcut. adipose RNA-seq counts; series title on perioperative dexamethasone): "paired non-surgical and surgical" biopsies of the same person and tissue, with and without DEX (titles such as A1_..._non-surgical_DEX / _surgical_DEX; batches A, B, C). Same-person, same-tissue pair, but clock time and interval are not in the header; surgical = general anaesthesia plus DEX acting within 0.5-2 h per the summary; not usable as a rest design. Not used in A252/A256 (not in either design table). Candidate only if the paper gives times (inf).

## Summary: arm type x (studies, persons, elapsed hours) among the 17 datasets scored or scorable in A252/A256 (exact times only for GSE118280 and GSE59088 t0, GSE106800)

| arm type | studies | persons / units (as in design tables) | typical elapsed h |
|---|---|---|---|
| nothing done (rest / saline, fasted) | GSE40589 INF_NACL; GSE59088 control; GSE252357 CTRL (rest, with standard meals, classed as rest arm in A252 text, not "nothing_done" in D256 because not in A256) | 7 + 6 + 5 = 18 persons | 4.0; 2.5; 1.75 / 4.25 / 9.25 |
| meal only | GSE118280, GSE142401, GSE231509 | 15 men (30 units); 19 men (38 units); 21 | 3.67 and 10.33; 4; 1 |
| insulin clamp / infusion | GSE9105; GSE40589 EU/HC +/- GIP, NACL clamps | 12; 8 + 8 (NACL) + 8 + 8 + 7 (GIP arms) | 1.5 and 5; 4 |
| exercise before or inside the interval | GSE252357 RE, GSE151066, GSE273041, GSE303983, GSE120862, GSE220899, GSE202295, GSE163356, GSE126296, GSE59088 endurance and resistance | 8, 19, 20, 45, 7, 11 + 10, 16 + 19, 10 + 9, 14, 6 + 6 | 0.67-9.25, mostly 2-5 |
| surgery / ischaemia (no time) | GSE4386, GSE75432, GSE151648 | 20, 10, 40 | NA |
| 12 h behavioural shift, stated clock times | GSE106800 | 11-12 men | 12 |

Persons in arms where nothing is done between the biopsies (fasting rest / saline control): A256 classification (nothing_done) = GSE40589 INF_NACL 7 + GSE59088 control 6 = 13 (D256; each study < 10). Adding GSE252357 CTRL = 5 persons (A252 text "only the rest arm of GSE252357 (5 persons) has none") gives 18, but the GSE252357 arm had a standard breakfast and lunch (lunch after the 3 h biopsy), so it is "rest with meals", not fasted. Persons are not shared between these arms (inf: different studies). Elapsed hours: 4.0 (GSE40589), 2.5 (GSE59088; t0 13:00 stated, so a different clock phase), 1.75-9.25 (GSE252357).

## Caveats for T3 (same-person, two biopsies, gene-level rotation)
- Only GSE118280 (3 biopsies, 08:40/12:20/19:00, with meals) and GSE106800 have clock times from the data; all others carry assumed 09:00 (GSE59088: stated 13:00 derived from protocol).
- Every interval except the 13 + 5 units above contains a meal, exercise, insulin, drug or surgery; "rest" biopsies still carry a wound response from the first biopsy.
- Genome-wide matrices exist for all; counts (not only normalised) for GSE231509, GSE151066, GSE303983, GSE151648, GSE75432, GSE202295, GSE252357 (normalised, protein-coding only), GSE273041 and GSE163356 (non-integer), GSE120862 (DESeq-normalised). Arrays log2 or log ratios: GSE118280, GSE142401, GSE4386, GSE40589, GSE9105, GSE59088, GSE220899, GSE126296, GSE106800.
- Numbers of persons are from S, D252, D256, ADD; the GSE252357 CTRL count 5 and RE 8 are from S.
