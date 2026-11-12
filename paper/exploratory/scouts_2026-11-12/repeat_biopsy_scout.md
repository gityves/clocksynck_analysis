# Repeat-biopsy scout, 2026-11-12

Goal: public human solid-tissue datasets with >= 2 biopsies of the same tissue per person within about 24-48 h and a known (or protocol-derived) elapsed time, to test whether the clock state rotates ~15 deg/h and whether interventions shift it. No clock statistic was computed; only sample sheets and detection summaries.

Sources: NCBI GEO E-utilities (gds) searches (62 queries, scripts/search_geo.py and search_geo2.py), GEO brief pages read for every candidate (scripts/geo_brief.py), GEO FTP listings (scripts/ftp_list.py). Downloaded data: `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/repeat_biopsy/<GSE>/` (96 MB total, each with PROVENANCE.md and sha256); platform annotations in `repeat_biopsy/platforms/`. Sample sheets: `output/repeat_<GSE>_samples.csv`; clock-gene detection: `output/repeat_<GSE>_clockgenes.csv`; summary: `output/repeat_summary.json`. Code: `scripts/build_sheets.py` (metadata + presence / median percentile only), `scripts/download.py`.

"Pairs" below = subjects with >= 2 samples of the same tissue in the same window (<= 48 h), counted from the samples actually present in the downloaded matrix. Elapsed hours are from metadata or the paper's protocol as stated per dataset; "NA" means GEO gives no interval.

## Ranked candidate table

Ranking criteria: n >= 10 subjects with pairs, elapsed time exact or protocol-bound, count data, clock time known, presence of a no-intervention arm.

| # | Accession | Tissue | n subjects (pairs) | Time points, elapsed h | Intervention between biopsies | Data type | Clock time of first biopsy known? | Downloaded |
|---|---|---|---|---|---|---|---|---|
| 1 | [GSE118280](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118280) | subcutaneous adipose | 15 men (15; 30 subject-diet visits x 3 biopsies) | 08:40, 12:20, 19:00 same day = 0, 3.67, 10.33 h | meals of HC/HF or HF/HC 4-wk diets (crossover, two visits) | RMA log2 array (Affy HuGene 2.0 ST, Entrez CDF) | yes, 08:40 | yes |
| 2 | [GSE252357](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252357) | skeletal muscle (vastus lat.) | 13 (RE 8, CTRL 5) | PRE, 30 min, 3 h, 8 h post (same day) = 0, 1.75, 4.25, 9.25 h; 24 h = 25.25 h | RE group: 45 min leg-press/extension; CTRL: 45 min rest, time-matched biopsies (no-intervention arm) | normalised counts, protein-coding only | no (1 h after breakfast, overnight fast) | yes |
| 3 | [GSE142401](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142401) | adipose | 19 men (10 control, 9 metabolic syndrome), 2 meals each (38 pairs) | 0 h fasting, 4 h post meal | high-fat breakfast (dairy or soy based), crossover | quantile-normalised log2 array (Illumina HT-12 v4) | no | yes |
| 4 | [GSE231509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231509) | skeletal muscle | 21 (7 healthy, 7 obese, 7 obese+T2D) | before, 1 h after meal = 0, 1 h | mixed-meal tolerance test | raw counts | arrival 09:00 after 12 h fast (biopsy time not given) | yes |
| 5 | [GSE4386](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4386) | right atrium | 20 (sevoflurane 10, propofol 10) | beginning vs end of off-pump CABG; elapsed NA | surgery + general anaesthesia | RMA log2 array (HG-U133 Plus 2) | no | yes |
| 6 | [GSE151648](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151648) | liver graft (TRANSPLANT) | 40 (IRI+ 23, IRI- 17) | pre-transplant vs post-transplant (post-reperfusion); elapsed NA | cold ischaemia, implantation, reperfusion | raw counts (HiSeq 3000) | no | yes |
| 7 | [GSE75432](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75432) | skeletal muscle | 10 | immediately before surgery vs 2 h after surgery with tourniquet (>= 2 h) | total knee arthroplasty + tourniquet ischaemia | raw counts per sample (RAW.tar) | no | yes |
| 8 | [GSE273041](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273041) | skeletal muscle | 20 (10 high-carb, 10 low-carb athletes) | BL, H0 (after 3 h run), H2 = 0, ~3, ~5 h | 3 h treadmill run | normalised (non-integer) counts | no | yes |
| 9 | [GSE151066](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151066) | skeletal muscle | 19 older adults (10 active, 9 sedentary) | pre, immediately post, 3 h post = 0, 0.67, 3.67 h | 40 min cycling | RSEM expected counts | no | yes |
| 10 | [GSE303983](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303983) | abdominal subcutaneous adipose | 45 (arm per subject not in GEO metadata) | pre, 1.5 h after exercise end | one exercise session (LOW 60 min / MOD 45 min / HIGH 10x1 min) | raw gene counts | no | yes |
| 11 | [GSE120862](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120862) | skeletal muscle | 7 only (4 leg x training units each, 28 units) | baseline, +1 h, +4 h after one-legged exercise (~0, 2, 5 h from baseline) | 60 min one-leg knee extension; other leg = within-person control; repeated after 2 months training | DESeq-normalised counts | no | yes (n < 10) |
| 12 | [GSE260643](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260643) | left ventricle | up to 73 patients (AM: 60 pre / 56 post samples; PM: 22 pre / 17 post) | pre-ischaemia vs post-ischaemia (~80 min cross-clamp) | aortic valve replacement, cardioplegia, cross-clamp | NOT public: only two DEG tables (AM_Pre vs PM_Pre, AM_Post vs PM_Post); raw reads in SRA only | yes: AM median 10:32 (8-12 h) vs afternoon | no (no per-sample matrix) |
| 13 | [GSE106800](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106800) | skeletal muscle | 12 subjects x 4 samples (47) | 19:00 vs 07:00 (12 h) after 3-day aligned and after 3.5-day misaligned protocol; TIMED design | 12 h behaviour-cycle shift | log2 array (13 MB matrix) | yes (07:00, 19:00) | no (budget; cheap to add) |
| 14 | [GSE15480](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15480) | liver graft (TRANSPLANT) | 12 pairs (IPC 6, standard 6) | end of cold storage vs 90 min after reperfusion | deceased-donor liver transplant | array (series matrix 3 MB) | no | no |
| 15 | [GSE87487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87487) | liver graft (TRANSPLANT) | 10 pairs (20 samples) | pre- vs post-reperfusion; elapsed NA | liver transplant | raw counts | no | no |
| 16 | [GSE9105](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9105) | skeletal muscle | 12 (36 samples) | baseline, 30 min, 240 min of 180-min hyperinsulinaemic clamp | euglycaemic insulin clamp | array | no | no |
| 17 | [GSE7146](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7146) | skeletal muscle | ~12 (24 samples, titles pre/post clamp) | pre vs post clamp | euglycaemic insulin clamp | array | no | no |
| 18 | [GSE40589](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40589) | subcutaneous adipose | 96 samples (0 vs 240 min, GIP vs NaCl, clamps) | 0, 240 min | GIP or saline infusion +/- clamps (NaCl arm is a no-drug 240 min control), morning fasted | array | morning, fasted | no |
| 19 | [GSE59088](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59088) | skeletal muscle | 18 (54 samples; endurance, strength, no-training control) | pre, 2.5 h, 5 h after last session of 12-wk training | exercise session | array | no | no |
| 20 | [GSE220899](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220899) | skeletal muscle | 64 samples (PRE, 3 h, 6 h, two loads one week apart) | 0, 3, 6 h | resistance exercise | array | no | no |
| 21 | [GSE202295](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202295) | skeletal muscle | 111 samples (NGT and T2D) | rest, post, 3 h recovery | 30 min cycling | counts | no | no |
| 22 | [GSE164081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164081) | skeletal muscle | 10 athletes, 50 samples | prior, 2 min, 1, 3, 6 h | 60 min intermittent cycling | CAGE/RNA-seq (no gene count table found in suppl list) | no | no |
| 23 | [GSE163356](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163356) | skeletal muscle | 11 men | pre vs 3 h post, session 1 and session 9 | HIIT bout | counts + logCPM | no | no |
| 24 | [GSE43219](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43219) | skeletal muscle | 14 | rest vs 30 min after exercise | single exercise bout | array | no | no |
| 25 | [GSE126296](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126296) | skeletal muscle | 14 | pre vs 2 h 20 min after sprint exercise (opposite leg) | 3x30 s sprints | array | no | no |
| 26 | [GSE108643](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108643) | skeletal muscle | 30 men (58 samples) | fasted pre, then breakfast, exercise 4 h later, biopsy immediately after (~4+ h) | breakfast + cycling | only Cuffdiff DE table public; no per-sample matrix | no | no (no matrix) |
| 27 | [GSE305038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305038) | skeletal muscle | 7 subjects (28 samples) | pre and 4 h after exercise, after active vs inactive day | exercise | only edgeR xlsx tables public | no | no (n < 10, no matrix) |
| 28 | [GSE111351](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111351) | right atrium | 15 samples, one per patient (AM vs PM surgery, different patients) | not paired | aortic valve replacement | array | surgery AM vs PM | no (not repeat biopsy) |
| 29 | [GSE108539](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108539) | skeletal muscle | timed series | around-the-clock (timed; already in `muscle_timed/`) | none | array | yes | already local |
| 30 | [GSE79934](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79934) | skeletal muscle | 5 samples (< 10 subjects; biopsies at 04:30, 08:30, ...) | timed | none | array | yes | no (n < 10) |

Not found / not available: no public dataset for design (e): anaesthesia without surgery. Search for start-vs-end-of-surgery adipose, liver (non-transplant), skin, gut mucosa and endometrium repeat sampling within hours returned no qualifying matrix. GSE17992 (dexamethasone, 4 days apart) and GSE124935 (GH vs saline in crossover, one biopsy per condition) do not fit.

## Summary of the downloaded datasets (from output/repeat_summary.json)

| GSE | samples in matrix | subjects | subjects with >= 2 same-window samples | genes in matrix | clock genes present | median percentile range of clock genes |
|---|---|---|---|---|---|---|
| GSE118280 | 90 | 15 | 15 | 25582 | 13 (CIART absent) | 51.3-93.6 |
| GSE252357 | 61 | 13 | 13 | 14392 | 14 | 38.6-95.2 (protein-coding only, so percentiles are against a smaller gene set) |
| GSE142401 | 76 | 19 | 19 (38 subject-meal units) | 20760 | 14 | 17.5-92.8 |
| GSE231509 | 42 | 21 | 21 | 59653 | 14 | 83.3-98.6 |
| GSE4386 | 40 | 20 | 20 | 22836 | 14 | 45.0-91.0 |
| GSE151648 | 80 | 40 | 40 | 60498 | 14 | 82.2-99.1 |
| GSE75432 | 20 | 10 | 10 | 64258 | 14 | 87.2-98.3 |
| GSE273041 | 60 | 20 | 20 | 25262 | 14 | 40.6-93.9 |
| GSE151066 | 57 | 19 | 19 | 60920 | 14 | 83.5-99.1 |
| GSE303983 | 90 | 45 | 45 | 61228 | 14 | 78.8-98.2 |
| GSE120862 | 84 | 7 | 7 (28 subject-leg-training units) | 56609 | 14 | 82.5-98.4 |

Percentile = share of genes in the matrix whose median (across samples) is <= the clock gene's median. Arrays are collapsed to gene symbol (highest-median probe); GSE118280 mapped by Entrez ID; BMAL1 and CHRONO are read as ARNTL and CIART in symbol-based tables. Counts-based NPAS2 and ARNTL are usually the lowest-ranked of the 14; in the arrays NPAS2 and NR1D1 are. Per-gene values are in `output/repeat_<GSE>_clockgenes.csv`.

## Per-dataset notes and problems

- **GSE118280 (best timed design).** Same person, three biopsies on one day at fixed clock times (08:40, 12:20, 19:00), two diet visits per subject (HF_HC and HC_HF, 45 samples each). The gaps are exact, and the clock time is known, so this separates elapsed time from time of day in a within-person design. Caveats: meals differ in composition between arms and occur between biopsies; the two visits are 4-week diets, so the diet itself is a between-visit perturbation; no age; male only; CIART (Entrez 400970) is not in the Brainarray CDF. Only the normalised matrix is on GEO (RAW.tar is 734 MB CEL files, not downloaded).
- **GSE252357 (only true no-intervention arm among the downloads).** CTRL (n = 5: BK, BM, BU, BX, BY; one has 4 samples) are time-point-matched biopsies without exercise; RE n = 8. Elapsed times are derived from the protocol (biopsy PRE at 1 h after breakfast; exercise or rest from +30 min for 45 min; later biopsies timed from exercise end), so 30 min = 1.75 h, 3 h = 4.25 h, 8 h = 9.25 h, 24 h = 25.25 h; all assumptions are stated in the sheet. Standardised meals (breakfast, lunch after the 3 h biopsy) occur in all arms. Counts are normalised and restricted to protein-coding genes (14392 genes). 24 h biopsy is the next morning, 8 h has 10 samples, 24 h has 12.
- **GSE142401.** Fasting 0 h and 4 h post high-fat meal, two meals per subject on different days, controls (n = 10) and metabolic syndrome (n = 9), men only; no age in GEO; clock time not stated. Probe level collapsed with GPL10558 annotation. The paper (PMID 32928582) should be checked for time of day.
- **GSE231509.** MMTT with biopsies before and 1 h after; arrival 09:00 after a 12 h fast; age 23-64 in the sheet (21 of 21 subjects have age). Individual numbers repeat across diagnosis groups, so subject id = group + individual. Only 1 h gap (15 degrees expected if clock follows).
- **GSE4386.** Beginning vs end of off-pump CABG, atrial tissue, 20 patients. Elapsed time NOT in GEO (operation duration unknown), so it can only test the sign of a shift or use the paper (PMID 17325502) for durations. Sex and age not parsed (GEO characteristics are sparse and misaligned).
- **GSE151648 / GSE87487 / GSE15480 (liver transplant).** Pre-implantation vs post-reperfusion biopsies of the same graft. Counts in GSE151648 (40 pairs, 23 IRI+). The recipient operation is at an arbitrary clock time, the donor liver spent hours in cold storage (clock stops or is not comparable), so this tests the effect of ischaemia and reperfusion on the clock readout rather than elapsed time. Elapsed time is NA in GEO (GSE15480 states 90 min after reperfusion). Flag as transplant.
- **GSE75432.** 10 TKA patients; the second biopsy is 2 h after surgery and tourniquet, so the true interval is longer than 2 h (surgery duration not given); the sheet records 2.0 as a lower bound. Count files per sample come from a 3.9 MB RAW.tar (parsed in memory, not extracted).
- **GSE273041.** BL rest, H0 immediately after a 3 h run, H2 two hours later (nominally 0, 3, 5 h, assuming BL immediately before the run); 20 elite athletes (10 high-carbohydrate, 10 ketogenic). The "raw_counts" file has non-integer values (normalised). No age or clock time in GEO.
- **GSE151066.** Pre, post, 3 h post a 40 min cycling bout; 19 older adults. The counts file has 78 columns, of which the 57 GEO samples were matched by column name; ids are Ensembl without version. No age/clock time in GEO.
- **GSE303983.** 45 regular exercisers, pre and 1.5 h after the end of one session; which of LOW/MOD/HIGH each subject did is not in the metadata (titles are all "3X-NN"), so the interval from the first biopsy is 1.5 h plus 20-60 min of exercise; the arm may be in the paper's supplement. Gene counts file ids are `ENSG|SYMBOL`; BMAL1 is the symbol used for ARNTL.
- **GSE120862.** Only 7 men (below the n >= 10 criterion) but a unique within-person no-exercise control (non-exercised leg), biopsies at baseline, +1 h, +4 h, in untrained and again after training; the baseline is described as 48 h after the previous exercise in the paper summary, so check the protocol before treating "baseline" as the pre-exercise sample of the same day. Counts are DESeq-normalised. Mapped from two GPL matrices (GPL16791 and GPL18573).
- **GSE260643 (high value, not usable as deposited).** Human left-ventricle pre- vs post-ischaemia biopsies with AM (60 pre, 56 post) and PM (22 pre, 17 post) surgeries, i.e. a natural timed experiment with the clock time of surgery known; but GEO holds only two DEG tables (AM_Pre vs PM_Pre, AM_Post vs PM_Post), no counts. Raw reads would have to come from SRA (outside the processed-matrices rule); the authors could be asked for counts (draft only, per CLAUDE.md "ask first").
- **GSE108643, GSE305038.** Count data not public (differential-expression tables only).
- **GSE111351.** Atrial biopsies AM vs PM from different patients, so it is not repeat sampling.
- **GSE106800.** Not downloaded; 12 men with 07:00 and 19:00 biopsies under aligned and misaligned conditions; timed design, array.
- **Overall gaps.** No dataset with a pure no-intervention arm, n >= 10 and same-day biopsies of gut, skin, endometrium or fat was found; GSE252357 CTRL (n = 5) is the only explicit no-intervention arm with counts. Clock time of the first biopsy is known only for GSE118280, GSE231509 (arrival only), GSE260643 and GSE106800.

## Files

- Scripts: `scripts/search_geo.py`, `search_geo2.py`, `geo_brief.py`, `ftp_list.py`, `gsm_meta.py`, `download.py`, `build_sheets.py`
- Output: `output/repeat_<GSE>_samples.csv`, `output/repeat_<GSE>_clockgenes.csv`, `output/repeat_summary.json` for the 11 downloaded GSEs
- Data and provenance: `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/repeat_biopsy/<GSE>/PROVENANCE.md`
