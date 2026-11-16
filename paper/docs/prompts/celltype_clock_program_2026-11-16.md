# Prompt for a lead agent: the cell-type organisation of the human clock (written 2026-11-16)

You lead the cell-type line of the paper "The human body keeps internal time" (Yves Paychere, EPFL, Naef lab; method
`clocksync`). Work autonomously; report as listed under "Reports and stops".

**Where it stands.** Without any host or anchor, the stromal-vascular cells of a tissue (fibroblasts, endothelium,
mural cells) share the donor's clock state with the tissue's other cells, and lymphocytes do not. This holds against
expression-matched and stress-responsive gene sets, replicates blind in fresh atlases, and matches timed mouse colon
(A272, A273, A274, A275). The cell classes of the cortex share a clock state too (A274). Cell types of one tissue
differ by well under an hour, so phase order gives no driver.

**What the paper still lacks, and what you are after.** (i) The result is not yet robust at the level of single
atlases. (ii) It says nothing yet about why lymphocytes are apart. (iii) It has not been used on disease, cancer or
the brain, where the anchored readout failed or was anchor-dependent. (iv) No human truth with known time exists at
cell-type level. Blocks A to D below address these in that order. A negative that removes a claim is a result: do not
tune a design to keep one.

## Read first
1. `CLAUDE.md` and the memory index (`MEMORY.md`): all entries apply.
2. `paper/docs/STATUS_2026-11-16.md` (state and findings; its section 6 lists the weak points this plan targets), then
   section 1 of `paper/docs/STATUS_ROADMAP_2026-11-14.md`.
3. RP (`paper/docs/results_preliminary.md`) entries A272, A273, A274, A275, A263 Part B, A271 Part A, A270, and their
   registrations in `paper/docs/registered_experiments.md` (A272 has a design addendum that fixes the statistic).
4. `paper/docs/PI_QUESTIONS.md`: Q35 (open: ordering of the paper), Q12 and Q6 (brain time of death, anchor-dependent
   brain contrasts).
5. `paper/docs/prompts/target_line_what_internal_time_controls_2026-11-15.md`: sections "Rules that bite here" and
   "Pitfalls" apply unchanged.

Before starting, run `git log` and `git status`. Another session works on the same branch (presentation figures,
benchmark runs): commit your own files only, never touch its uncommitted files or its jobs.

## The tool you inherit
- **Statistic** (A272 addendum): for two compartments of one atlas, the complex correlation over donors of their
  clock-plane states, centred within stratum and pooled over strata, **without** the Poisson correction (it
  over-corrects at pseudobulk counts). |r| is the shared state, arg r the offset. Class statistics: median |r| over
  pairs of a class per atlas, mean over atlases; T_SV (stromal-vascular with other families), T_LY, D = T_SV - T_LY.
- **Companions:** 200 expression-matched and 200 stress-responsive sets per atlas
  (`clocksync.paired.pseudo_waveforms(..., kind="matched" | "responsive")`), set j pooled across atlases as set j. A
  result is clock-specific only above the q95 of both families.
- **Code:** `exploratory/A272_backbone_pairwise/scripts/score_a272.py` (functions `pair_index`, `pair_r`,
  `atlas_states`, `pseudo_sets`, `placebo`), generalised in `exploratory/A274_celltype_map/scripts/score_a274.py`
  (any list of atlases, family grouping, eligibility by twin check, `--stage design | score`, `--placebo`). A whole
  run takes about 30 s with 6 workers.
- **Genome-wide tables** (needed for companions): `$SCCC_DATA/atlas_wg/<folder>/pseudobulk.h5ad`, loaders in
  `common/datasets/_atlas/genome.py` and `common/datasets/<name>_genome/build.py`, read with
  `results/figure_2/scripts/gene_rhythms/_shared.py::genome_table(name)`. New tables are built from CELLxGENE Census
  by `exploratory/celltype_map_builds/run.sh` (template; env `census_env`; about 20 to 60 min per atlas). An atlas
  enters an analysis only if its genome table reproduces every clock count of its clock-gene table (`twin.csv`).
  `islets` and `lung_natri_control` have no exact twin; do not approximate them.
- **Offsets:** unbiased against known time at the 1 h scale (A273). Report an offset only for a pair above both
  companions, with its donor-bootstrap interval.

## Block A. Make the result robust (do first; cheap, no new data)

- **A1. False-positive rate of the pairwise test.** The A274 placebo had one of six tests at its q95 and all six
  clock values above the companion medians. Run the placebo over many independent draws (at least 200) on the A272 and
  A274 atlases and report the rejection rate of each test at the q95 rule. Register as a calibration (V-type) with
  the acceptance band fixed in advance. If the rate exceeds the band, find the cause (the clock plane against
  companion planes under permuted donors) before any further test of this line is scored, and report.
- **A2. One pooled analysis.** A272 (five atlases, not blind for the anchored readout) and A274 (eight, blind) were
  scored separately. Register a descriptive pooled table over all thirteen: family by family matrix with both
  companions, per-atlas rows, and the split **living** (surgical or biopsy tissue) against **post-mortem** donors,
  taken from each atlas's metadata. The question fixed in advance: does the ranking hold in living tissue alone,
  where no agonal state is shared by the cells of a donor? No new test statistic; label post hoc what is.
- **A3. More lymphocyte atlases.** The lymphocyte side rests on few atlases. Scout which healthy atlases with
  lymphocytes and stromal-vascular cells from the same donors are still without a genome table (the A17 target list;
  `$SCCC_DATA/atlas_inventory`), build those that can have an exact twin, and score them as a second replication of
  H274.1 with the same statistic. Report sizes before building; the data cap is 500 GB (438 GB used).
- **A4. Mouse calibration beyond colon.** Find timed mouse single-cell datasets with several cell classes per animal
  in other tissues (liver, aorta, skin, lung are candidates; check what is on disk first: `mouse_skin_timed`,
  `aorta_jetlag_auerbach2026`, `mouse_ageing`). For each usable one repeat A273: timed coupling per cell class, and
  untimed against timed pairwise offsets. Registered question: is the order stromal-vascular, parenchymal, myeloid,
  lymphocyte the same in every tissue?

## Block B. Why are lymphocytes apart? (new biology, on data at hand)

Two readings are open after A275 (lymphocytes share glucocorticoid-target and immediate-early states with the tissue,
not the clock state): a weak clock, or a clock that keeps the time of where the cell came from.

- **B1. Amplitude against sharing.** Per compartment, the relative clock amplitude and its noise floor (the readouts
  of `Dataset.fit`; see METHOD on what |z| measures and the count floor from A271 Part B post hoc). Do lymphocytes
  have a clock-plane spread above their counting noise that is simply not shared (own time), or none (weak clock)?
  Lymphocyte pairs with each other are low as well (A274 family matrix: 0.214 against matched 0.255), which favours
  a weak clock: test it rather than assume it. Check power first on a class with a known clock at matched counts
  (downsample fibroblasts to lymphocyte counts and see what survives).
- **B2. Residency.** Split immune classes by residency where the atlas annotation allows it without new clustering:
  tissue-resident against circulating T cells, plasma cells against B cells, resident macrophages against monocytes.
  Mouse colon gives the expectation (macrophage 0.514, plasma 0.495, B 0.167, T 0.000). Registered hypothesis:
  sharing with the tissue increases with residency. Fix the marker or label rule per atlas before any statistic.
- State in advance what each outcome means for the text of the health chapter.

## Block C. Use the anchor-free readout where the anchored one failed

- **C1. Disease contrasts by cell class.** Inventory first, no statistic: for every atlas with a disease arm
  (Kong ileum and colon: Crohn's; kidney: chronic kidney disease; lung: fibrosis; PsychAD and SEA-AD: Alzheimer's,
  schizophrenia, bipolar disorder; liver fibrosis), list donors per arm, pairs per class, and what earlier amendments
  already scored with the anchored readout (A230 inflamed gut: no epithelial shift, MDE 1.7 h; A167 liver fibrosis:
  no change; A141 / Q6: brain contrasts anchor-dependent, moved to Extended Data). Then register only contrasts with
  power: does a named cell class lose its shared state with the stromal-vascular cells, or shift against them, in the
  disease arm? The project rules bite here: sharing scales with between-donor spread, so compare arms at a common
  spread or on offsets; every contrast needs its specificity companion on the same permutations; check the class's
  baseline in controls and the detectable difference before choosing it. The anchor-free readout removes the
  anchor dependence that sank the brain contrasts; say so only if the stability check (A141 design) confirms it.
- **C2. The tumour in the same language.** A227 pairs tumour compartments with normal tissue through the anchored
  readout; fibroblasts and endothelium clear the responsive control, immune cells do not (A271 Part A). With the
  pairwise readout inside the tumour sample: do malignant cells share the clock state of the tumour's own
  stromal-vascular cells, compared with the epithelium of the adjacent normal sample and its stromal-vascular cells?
  Genome-wide tumour tables exist for the A187 cohorts (`exploratory/C4_single_cell/A187/output/genome_*.csv.gz`);
  read Q24 and the cancer rules in the memory before designing. Depth was the limit of earlier single-cell tumour
  questions (roadmap 1.5): do the power calculation first and stop there if it fails.
- **C3. Cortex, third cohort.** PsychAD and SEA-AD disagree on the sign of the neuron-vascular offset, and each is
  below the matched q95 on its own. `psychad_mssm` and `psychad_radc` exist as clock-gene datasets
  (`common/datasets/_atlas/cortex.py`); build their genome tables and score them as a replication of H274.3 and as
  the tie-break of the offset sign (registered with the sign rule fixed in advance).

## Block D. What needs the PI or new access (prepare, do not send)

- **D1. Human truth with known time.** Brain cohorts with single-nucleus data and recorded time of death (the Rush
  cohorts behind `psychad_radc` and the ROSMAP single-nucleus sets; verify at the source what exists and what access
  it needs). A request draft already exists and is the PI's to send: Q12,
  `exploratory/A141_brain_stability/drafts/tod_request_draft.md`. Update the draft's scientific case with the
  cell-type result (what the times would test: the A273 calibration in humans, and the cortex offset sign), write it
  as a new question in `PI_QUESTIONS.md`, and stop there. Bulk brain cohorts with time of death are on disk
  (`$SCCC_DATA/brain_tod`): they cannot give cell classes; do not rescore them.
- **D2. Clock and cell cycle.** In gut epithelium and cycling T cells, is the distribution of cell-cycle phases of a
  donor's cells related to the donor's clock state (mouse: division is gated by time of day)? CoPhaser gives
  cell-cycle phase per cell; the memory entry on CoPhaser fixes how it may be used. Scout only: which atlases have
  the depth (median UMI above 1.5k), cells per donor, and what the power would be. Write the outcome as a PI
  question: this is a project of its own and uses the lab's unpublished tools.
- **D3. Liver zonation.** Do pericentral and periportal hepatocytes of a person share one clock phase, as in the
  mouse (Droin 2021, the lab's work)? Scout only: zonation score from genes fixed in the mouse literature, cells per
  donor and zone in `liver_xu` and other liver single-nucleus cohorts. PI question with the power.

## Later, each to be registered when blocks A and B are done
- **Outputs per cell type.** Rescore the per-cell-type gene rhythm maps (A50; `results/figure_2/scripts/gene_rhythms/`)
  with the stromal-vascular state as reference and both companion families: which genes of which cell class follow
  the donor's clock state. This is the cell-type version of the "what internal time controls" question; keep the
  masking caveat (following the state is not being driven by the local clock).
- **Cell type against organ.** Partition the offsets of all atlases into a cell-class term and an organ term.
- **Age.** PsychAD HBCC (145 healthy donors) and the eye atlases: does the sharing of any class change with age, at
  a common spread?

## Order and lanes
A1 first, alone: everything else depends on the test being calibrated. Then A2 and B1 (minutes each, data at hand)
while at most 2 Sonnet 5.5 subagents run clearly defined metadata or plumbing tasks (A3 and A4 scouts, C1 inventory,
genome loaders): verify what they return against the files. Then B2, C3, C1, C2. Block D scouts whenever a lane is
free. Machine budget: at most 21 threads and 83 GB of summed caps, at most 2 subagents.

## Registration and writing
- Every hypothesis is an amendment in `registered_experiments.md`, committed before its data are scored; next free
  number A276 (check the end of the file). Use staged registrations (hypotheses first, design addendum with atlases
  and pair counts before scoring) as A274 did. Placebo before every scoring run. State what has been seen.
- After each scoring: RP entry (headline with verdicts, design, placebo, each readout, per-atlas values, reading as
  registered, limits, code line, `<!-- src: -->`), rows in `VERDICT_TABLE.md`, `METHOD.md` when a tool or
  calibration is new. Numbers from output tables only. Commit with the attribution lines.
- The anchor-free pairwise readout is used by three amendments from script folders. If A1 confirms its calibration,
  move it into `clocksync` as a public function with tests (CLAUDE.md: "Modify clocksync if tools potentially useful
  for future users are missing"); the package is its own git repository, with foreign uncommitted changes that are
  not yours to touch.

## Reports and stops
- Report without waiting: after the first `git` check; after A1; after each scout or inventory with what you intend
  to register; after each scoring.
- Stop and wait: (1) A1 outside its band and the cause not found; (2) anything outward-facing (access requests,
  messages to authors or data owners: drafts only); (3) blocks A to C scored: write the tally as a new question in
  `PI_QUESTIONS.md` after Q35, with what the cell-type section can claim and your reading for the ordering of the
  paper, which stays the PI's choice.
- A report states verdicts with their numbers, the power behind every negative, what is not shown, and what comes
  next.

## Not in scope
Driver-follower readings from phase order (offsets are below the resolution, and A273 shows none in the mouse);
brain-death mechanism; further plan 1 scouts; Phase 8 (naming, figures, tutorials); pushing to the remote (the
agent's push is blocked; the user pushes); GTEx per-donor data in any output.

## State at hand-over (2026-11-16)
Last commits of this line: `943d7cb` (A274 scored) and the commit adding this prompt. Registered and unscored: none.
No background job of this line (three benchmark runs of the other session may be running or stuck: not yours).
