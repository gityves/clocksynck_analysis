# Prompt for a fresh lead agent: fetal clocks, pooled organs (A243) and adrenal cell classes (A244) (written 2026-11-11)

You lead two registered analyses for the paper "The human body keeps internal time" (Yves Paychere, EPFL, Naef lab):
**A243** (is a clock of the size read in the fetal adrenal present in the heart, kidney, lung or liver when open
cohorts are pooled?) and **A244** (is the fetal adrenal signal carried by the cortex, the medulla, or both?). Both
are registered in `paper/docs/registered_experiments.md` as *staged* amendments: hypotheses, statistic, rules,
families and readings are fixed; you write a **design addendum** for each (cohorts, classifier thresholds, adult
power), commit it, and only then compute a fetal clock statistic, once. Work autonomously; report to the user after
each addendum and after each scoring.

## Read first
1. `CLAUDE.md` (repository root): rules, budget, data policy, long-job procedure, git (local commits only).
2. Memory index (`MEMORY.md`), in particular: PI questions file, machine budget, blood has no coherent core clock,
   estimable bulk designs, CoPhaser single-cell, Phase 8 needs PI review.
3. `paper/docs/registered_experiments.md`: **A243 and A244 in full** (end of file), then A242 (the design they copy)
   and A237 (design points of the developmental series).
4. `paper/docs/results_preliminary.md`: A242, A240, A241, A237-A239, A190 with its post hoc checks, A189.
   `paper/docs/PI_QUESTIONS.md`: Q32 (the developmental section; these two amendments were registered on the user's
   instruction while Q32 is with the PI).
5. Code to reuse, not rebuild:
   - `paper/common/development.py` (`design`, `adjust`, `coherence(strata=)`, `surrogates`, `l1_test`, `shared`,
     `l2_prepare`, `l2_pseudo`, `l2_row`, `score`, `verdicts`, `main`) and `paper/common/tests/test_development.py`.
   - `paper/exploratory/A242_fetal_adrenal_replication/scripts/`: `prep_adrenal.py` (parsers for 10x, Microwell,
     DeepSAGE, ENCODE, cellxgene objects; fetus identity), `feasibility_a242.py` (adult stand-ins, thinning, caps,
     settings a1 / a0.75 / a0.5 / shuffled), `score_a242.py` (pooled tests on `development.main`), `position_a242.py`.
   - `paper/exploratory/A241_fetal_adrenal_host/scripts/score_a241.py` (`contamination`: marker-based foreign RNA
     and its subtraction), `paper/exploratory/A190_fetal_adrenal/scripts/posthoc_checks.py`.
   - Loaders: `common.data.load("fetal_cao", level="organ")`, the Cardoso-Moreira loader used by
     `paper/exploratory/A237_fetal_organs/scripts/`, `common.gtex.load([...], hardy=..., genes=..., depth=...)`.
   - `paper/exploratory/scouts_2026-11-10/FETAL_SCOUT.md` (what was already looked at, what is controlled access).

## Where the question stands (numbers from `results_preliminary.md`, A242)
- Fetal adrenal, nine open cohorts, 59 fetuses, none from Cao 2020: `late` (from 12 weeks, 19 fetuses) C 0.552, p
  0.022; `schroder` (13 bulk) C 0.651, Holm p 0.0048; `early` (before 12 weeks, 25 fetuses) C 0.732, p 0.0007; all
  clock-specific. `all` (46 fetuses, descriptive) C 0.750, p 0.0006. Against adult stand-ins the fetal signal lies
  between half and three quarters of the adult between-individual clock variation.
- The other organs are negative in single cohorts (A189, A237-A239), which excludes a clock of 0.75-1 times the
  adult one only (A240). So "the adrenal starts first" is not shown: that is what A243 can decide, organ by organ.
- Whether the adrenal signal is cortical, medullary or a matter of cell proportions is open (A190 was confounded by
  cortex RNA): that is A244.

## Order of work
**A244 first** (all data on disk: `$SCCC_DATA/natexp/{kameneva2021,wang2022,cheng2022,gse291218,kildisiute2021,
dong2020,han2020}/raw/`), **the A243 scout in parallel** through at most 2 Sonnet 5.5 subagents.

### A244
1. Prep (`paper/exploratory/A244_fetal_adrenal_classes/scripts/`): per cell, the marker counts of the registered
   lists and the total; fix one per-cell rule and its thresholds on the marker counts alone; compare with the
   published annotation where a source has one (the Kildisiute object does). Build the class pseudobulks (>= 200
   cells per fetus and class). The classifier must not read the clock genes or any pseudo-clock set.
2. Foreign-RNA fractions per fetus and class, design table, adult feasibility per class (thin GTEx adrenal to the
   class and cohort clock-gene counts); for H244.2 choose the second adult tissue and calibrate the gate (pass rate
   and false-positive rate with mis-paired donors). If the gate needs pairing permutations within cohort, add that to
   `development.py` (or to clocksync if it is of general use) with a test.
3. Write `### A244 design addendum` under A244, commit it with the prep and feasibility scripts. Then one scoring
   run, then the docs (below).

### A243
1. Scout (metadata only, no clock statistic): open human fetal heart, kidney, lung and liver data with gene-level
   counts and an age per sample. Search GEO, BioStudies / ArrayExpress, CELLxGENE, the HCA data portal, ENCODE
   (fetal tissue RNA-seq of the same programme as the A242 `encode_pa` cohort), descartes (Cao 2020), the Human Cell
   Landscape (Han 2020, GSE134355, on disk for the adrenal), Roost 2015 (GSE66302, on disk for the adrenal), and the
   reference lists of the fetal atlases. Apply the eligibility rules E1-E6 of A243 literally; keep a table of every
   source considered and the rule that excludes it. Verify what subagents report (accessions, sample counts, ages,
   file formats, sizes) against the source metadata before using it. Report download sizes to the user before
   downloading more than 20 GB in total; A243 may add at most 60 GB (cap 500 GB, about 408 GB used on 2026-11-11).
2. Download (processed matrices only; each source in its own new directory under `$SCCC_DATA/natexp/<source>/raw/`
   with a `PROVENANCE.md`: URL, date, sha256; parse with `python -I`, never execute), prep whole-organ pseudobulks
   and the cohort design and clock-detection tables (sums only), adult feasibility per organ and readout.
3. Write `### A243 design addendum` under A243, commit with the prep and feasibility scripts. Then one scoring run
   for the readouts that are powered, the position among the adult subsamples beside the adrenal, then the docs.

## After each scoring
- `paper/docs/results_preliminary.md`: one entry per amendment (headline with the verdicts, power first, each
  registered readout, descriptives, reading as registered, limits, code line with runtimes, `<!-- src: -->`).
- `paper/docs/VERDICT_TABLE.md`: one row per registered readout. `paper/docs/METHOD.md` when a tool or a calibration
  is new. `paper/docs/PI_QUESTIONS.md`: update Q32 with the outcome (a positive in another organ, an informative
  negative at half size, or a cortex negative all change the developmental section: say what changes and what does
  not, and give a recommendation). A result that needs the PI before anything else can follow: new question, pause
  the dependents.
- Commit your own files only (another session commits presentation work to the same branch), with the attribution
  line. Then report to the user: verdicts, adult power, what it means for "which organ starts first", what is not
  shown.

## Rules that bite here
- **Nothing is chosen after seeing a fetal clock statistic.** The addendum may fix only what the registration leaves
  to it. A readout whose adult pass rate at the full clock is < 0.8 is "not estimable" and is not computed. Any
  further cut (another organ, another window, another class) is a new amendment: propose it, do not run it.
- Numbers in the docs come from the output tables, never from memory or from a subagent's message.
- A negative is read with its adult pass rate: informative only at the size the power table reaches; never "no
  clock". The comparison with the adrenal follows the rule written in A243 (pass rate >= 0.8 at 0.5).
- Fetal liver is mostly blood cells at these ages; blood carries no coherent core clock in the adult. Do not explain
  a liver result away post hoc: report it with that limit.
- GTEx time of death and phenotypes are controlled access: aggregates only in outputs, logs and docs.
- Controlled-access fetal data (EGA, dbGaP, Synapse) and any contact with data owners: ask the user first; nothing
  of the kind is needed for A243 or A244 as registered.
- Long runs (> 20 min) only through `RUN_CAP_GB=<gb> paper/common/run.sh <script.py> [args]`, detached, with
  `paper/common/watch_jobs.sh` in the background; check `paper/common/jobs.sh` before reading any output (a run can
  write its tables and still fail). At most 16 workers and about 100 GB summed, 20 GB left free.
- Do not touch the uncommitted changes in the `clocksync/` repository that are not yours (`.gitignore`,
  `pyproject.toml`, `tests/test_additions.py`, `src/clocksync/utils.py`, `tests/__init__.py`).
- Do not start Phase 8 work (figures, naming, tutorials).

## Pitfalls met in A242
- `pd.read_csv` reads a setting named "null" as NaN: the no-clock setting is called `shuffled`.
- Companion sets are drawn among genes present in every cohort of a readout: check the gene lists in a dry run (no
  statistic) before scoring; the responsive sets draw from few genes (14-18 in A242).
- Adult libraries cannot be thinned upwards: `feasibility_a242.cap` lowers a gene's targets over a cohort; report the
  share of capped targets.
- If the adult pool of an organ holds fewer donors in the 08-18 h window than the design has fetuses (kidney cortex
  is small in GTEx), draw distinct donors within a cohort and say so in the addendum.
- Samples deposited in two series by one laboratory can be the same cells (GSE193672 / GSE291218): compare barcodes.
- The shell working directory resets between commands: use absolute paths. Python:
  `/home/maxine/miniconda3/envs/sccc/bin/python`.

## State at hand-over
A242 scored and documented (commit `5642801`); A243 and A244 registered, nothing prepared, no job running, no data
of the four organs downloaded beyond what earlier amendments use (Cao 2020 and Cardoso-Moreira 2019 are on disk;
the Roost 2015 table on disk was read for the adrenal only; of Han 2020 only the adrenal files are on disk).
