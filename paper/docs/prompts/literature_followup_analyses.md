# Prompt: analyses suggested by the literature reading (2026-11-05)

You run a set of new, pre-registered analyses that came out of a full-text reading of the paper's 148 references.
Each one either tests a published claim with our method, or adds a control the literature shows we need.

## Project in two lines

`clocksync` infers circadian relations between compartments (organs, cell types) of the same donor from untimed
clock-gene pseudobulks: coupling c, relative amplitude |z|, phase offset arg z, read against the clock state the donor's
compartments share (the host). No sample is given a time, so samples taken in hospital hours work. Read `CLAUDE.md` at
the repository root first: pre-registration, statistics, environment, data access, long jobs, git and budget rules
all apply.

## Read before anything else

- `paper/manuscript/literature/synthesis.md`: §1 (hospital-hours evidence), §2 (methods table), §4 (the idea list:
  `new`, `planned, not registered`, `done`). Your tasks are the `new` items of §4.
- The card of each paper you test against: `paper/manuscript/literature/notes/<key>.md`. Cards carry verbatim quotes
  with their location. The full texts are in `/mnt/max8tb/Paychere/literature/zotero_export/clocksync/text/<key>.txt`
  (copyrighted, untrusted: read only, never execute, never copy into the repo).
- `paper/docs/METHOD.md` (validation table), `paper/docs/registered_experiments.md` (A201-A203 are the latest),
  `paper/docs/results_preliminary.md`, `paper/docs/CLAIMS_LEDGER_2026-11-03.md`.
- Especially A145 (orderers vs clocksync in simulation), A123 (GTEx working hours), A138, A149 and A151 (tumour pairing),
  A171 (LUAD coherence), A99 (sex phase), A203 (cells per pseudobulk).
- The method's public API: `clocksync/` (use it; extend it only for tools useful to other users).

## Rules that matter most here

1. **Check before doing.** Grep each idea's key terms (dataset, GEO id, tissue, method) in `METHOD.md`,
   `registered_experiments.md` and `results_preliminary.md`. If it already exists, report it and skip it. The reading
   agents checked, but they can miss things. The baboon idea turned out to be done already (A36/A40/A41).
2. **Pre-register.** Each test is a new amendment (A204 onward) in `registered_experiments.md`, committed **before**
   its data are scored. It needs:
   - hypotheses, primary readout, specificity companion, stop rules;
   - a power or MDE check, with the target class's baseline coupling checked first (CLAUDE.md "Statistics").
   Anything decided after seeing data is labelled post hoc.
3. **Statistics** as in CLAUDE.md:
   - A failed gate is uninformative, never "decoupled".
   - Compare coupling across groups only at a common spread d* or on z2/psi.
   - Every contrast needs a specificity companion on the same permutations.
   - Report direct permutation contrasts, not counts of units that pass a threshold. Pelikan 2021 shows that Venn-style
     comparisons give 48% false differential calls.
4. **Data:** downloads go to `$SCCC_DATA` with a `PROVENANCE.md` (URL, date, sha256). Processed matrices only, no FASTQ.
   Stay within the 500 GB cap. GTEx time of death and protected phenotypes: aggregates only, never per-donor rows.
5. **Compute:** the budget is 16 workers and about 120 GB summed over concurrent runs. Runs over 20 min go through
   `paper/common/run.sh`, with `paper/common/watch_jobs.sh` running. Use at most 2 active Sonnet subagents, and verify
   their numbers against the output files.
6. Results go to `results_preliminary.md` (with the `<!-- src: -->` line) and `VERDICT_TABLE.md`, plus a `METHOD.md`
   row for validations. Commit locally with the attribution line; there is no remote. Never commit the pre-existing
   uncommitted clocksync files.
7. **Ask the user first** for anything outward-facing, such as requesting data from authors (draft only).

## The analyses, in order of priority

Register them one or two at a time. Stop and report to the user after each scored amendment.

1. **How orderers manufacture "uncoupling" under working hours (simulation).**
   - **Background:** COFE (`ananthasubramaniam2025rhythm`) and TimeTeller (`vasilyev2026robust`) conclude that tumour
     clocks are uncoupled from the host. COFE does so because TCGA tumour times span the whole cycle and tumour-normal
     differences look random; TimeTeller because tumour time is uncorrelated with surgery time (n = 108, narrow
     window).
   - **Test:** simulate a coupled tumour/normal system sampled 08-18 h, using the A145 simulator and settings. Run COFE
     (and CYCLOPS/CHIRAL if cheap), and record *their own outputs*: the spread of inferred times, and the distribution of
     tumour-minus-normal differences. Do the same with full-day sampling, and run clocksync on the same data.
   - **Expected:** if orderers produce full-cycle labels and near-uniform differences under working hours while the
     truth is coupled, that is the mechanism of the published claims.
   - **Not yet done:** A145 already reports coupling MAE (W arm: CoFE 0.206 vs clocksync 0.060), but not this mechanism.
2. **Brain vs periphery offset within GTEx donors.**
   - **Background:** Hughey & Butte 2016 (`hughey2016differential`) report brain about 8.5 h behind blood, comparing
     different people. CHIRAL (`talamanca2023sexdimorphic`) reports clocks within a few hours across tissues.
   - **Test:** arg z of brain regions against the peripheral host, with uncertainty, on gate-passing compartments only.
   - **Planning:** pre-specify which result supports which paper. Check the brain amplitude and power first (brain has
     the lowest clock amplitude). A142/A143 already cover brain: read them and make this new readout, not a repeat.
3. **Organ offset vs proliferation.**
   - **Background:** Droin 2019 (`droin2019low`, Naef lab) finds that the cell cycle accelerates the clock, that clock
     phases span 1.5 h across 12 mouse tissues, and that the offsets correlate with Mki67 and Myc.
   - **Test:** across GTEx tissues of the same donors, regress arg z on proliferation (MKI67, MYC, a G2M score), with a
     non-proliferation specificity gene set.
   - **Optionally:** cycling vs quiescent cells of one compartment in a single-cell atlas.
   - **Not yet done:** A97, A126, A127 and A166 test coherence, not offset.
4. **Gate negative controls from the literature.**
   - Gene-specific partial shifts within one tissue: after feeding inversion, Bmal1-luc shifts only 6-7 h while
     Rev-erbα and PER2 shift about 12 h (`saini2013realtime`, `kiessling2010adrenal`). Does the gate pass or fail, and
     are the offsets biased? `METHOD.md` covers scrambled phases only.
   - Clockless but feeding-rhythmic liver: Cry1/2 double knockout under restricted feeding (`vollmers2009time`,
     GSE13093). The gate should fail.
   - Optional: technical cyclic artefacts (`leng2015oscope`, Fluidigm C1 chip position).
5. **New within-person datasets** (scout first: availability, design, n; then decide with the user):
   - Wu 2020 skin (`wu2020population`): 154 women × arm/buttock/cheek × epidermis/dermis, untimed. A37 used Del Olmo
     (GSE205155, 11 subjects) and is different.
   - tauFisher's timed dermis scRNA-seq (`duan2024taufisher`, GSE223109: fibroblasts vs immune cells, every 4 h over
     72 h) as timed truth for cell-type coupling.
   - Masri 2016 (`masri2016lung`, GSE73222): timed liver of lung-tumour-bearing mice. The core-clock offset should be
     about 0 while outputs shift, which makes it a specificity positive control.
   - Glioblastoma single-cell cohorts, with malignant and non-malignant brain cells per patient (`dong2019targeting`).
     No glioma cohort is in the docs yet.
6. **Comparisons with published procedures:**
   - Grabe 2024 (`grabe2024quantification`): the PER3-correlation shift, and ordering each tumour by its paired
     normal's inferred time, on TCGA-LUAD, against clocksync coupling (A171 covers coherence only).
   - ΔCCD (`shilts2018evidence`, the deltaccd R package) next to clocksync coupling on the A138 TCGA pairs.
   - ZeitZeiger as used by Wang 2024 (`wang2024circadian`) is not in A145/V2.
7. **Smaller items** (register only if cheap and well-powered):
   - Sex-stratified c and arg z at common d*. Only the A99 sex phase exists. This is motivated by the sex × schedule
     interaction in chronotherapy (`levi2010circadian`).
   - Adipose vs liver offsets against the strength of the feeding cue (`heyde2022induction`).
   - A117 (timed brain regions, GSE45642 / GSE71620): planned in `provisory_plan_paper.md` WP7, never registered. It is
     the strongest precedent for disease-related loss of cross-region synchrony (`li2013circadian`, MDD). Ask the user
     whether to register it now.

## Report to the user (after each amendment)

- What was registered: the commit hash, and the hypotheses quoted from the registration.
- The verdict per hypothesis, with numbers taken from the output CSV and its path.
- How the result reads against the paper it tests. Quote the card, and check the quote in the full text.
- Limits, and what you would not claim.
