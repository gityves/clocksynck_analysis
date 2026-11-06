# Prompt for a fresh lead agent: pooled single-cell round on tumour clocks (written 2026-11-06)

You lead a new round of analyses for the paper "The human body keeps internal time" (Yves Paychere, EPFL, Naef lab).
You register, run and record analyses. You delegate well-defined tasks to Sonnet 5.5 subagents (at most 2 active), check
their numbers against the output files, and report to the user after each amendment. Work autonomously. Before
proposing anything, check that it was not already done (grep `METHOD.md`, `registered_experiments.md`,
`results_preliminary.md`, `VERDICT_TABLE.md`), and reuse the existing code and loaders instead of rebuilding them.

## Read first
1. `CLAUDE.md` (repository root): rules, budget (at most 16 workers, about 120 GB summed), data policy, long-job
   procedure (`paper/common/run.sh`, `watch_jobs.sh`), git (local commits only, never push).
2. Memory index (`MEMORY.md`): paper framing, brain-death framing, the proliferation-test history, check-before-proposing.
3. `paper/docs/NOVELTY_ASSESSMENT_2026-11-06.md`: where the paper stands, and section 6 (the user's decisions of
   2026-11-06).
4. `paper/docs/results_preliminary.md`: the entries A176 (with its 2026-11-06 correction), A209/A209b, A210, A211, A212,
   A213, and the single-cell entries WP0 item 2 (2026-10-24, second and third parts: per-cell D), A134, A138, A149, A187,
   A187b, A188, A195. `METHOD.md` section on D, coherent amplitude, lambda and R.
5. The code of these single-cell analyses: `paper/exploratory/A138_cancer_estimability/`, `A149_paired_tumour_normal/`,
   `C4_single_cell/` (A187 holds per-study malignant pseudobulks, inferred-genome and match files; `fits/` holds
   CoPhaser cell-cycle fits of malignant cells for 18 studies, with per-cell files, per-patient cycling fractions and
   `qc_rerun_all.csv`, on disk but not committed; `REGISTRATION_DRAFT.md`), `WP0_integrity/`
   (per-cell D and its calibrated simulation), and the clocksync API (`decompose`, `paired`, `simulate`).

## Where the question stands (all numbers from the files named in results_preliminary.md)
- Bulk: proliferating and ER− tumours carry less in-plane clock signal (A176, 7 of 7). "Swamped, not disordered" was
  withdrawn: H176.2 was not informative.
- A212 (plasmode): for the proliferation loss, swamping by non-clock variation is excluded (in-plane variance falls,
  off-plane flat). A damped or dephased clock and mildly rewired gene relations (about 2 h per-gene shifts) are not
  separable in bulk (AUC 0.89-0.91, bar 0.9).
- A213: subclonal genome fraction does not lower coherence at matched proliferation (TCGA, T +0.021, MDE 0.058). Copy-number
  burden does in TCGA only (−0.063, clock-specific), and the A210 MYC result survives burden matching (ratio 0.77).
- Single cell so far: malignant coupling to the patient's other cells c 0.076 [0.014, 0.139] over about 225 patients
  (A138 H138.1, under the bound). Per-cell D, malignant 3.89 vs normal epithelium 3.90 (79 patients). The calibrated
  simulation gives 4.70 coherent, 0.42 scrambled and 2.85 with 3 h jitter, so full scrambling is rejected but mild
  rewiring is not. Malignant cells carry about 4.3 clock UMI per cell. Non-cycling malignant clock signal is at the 42nd
  pseudo-clock percentile (A195); cell-cycle locking was not estimable (A187).

## The round: can pooled public single-cell tumour data separate the explanations?
Questions, in priority order:
1. **Rewired vs intact gene relations in malignant cells.** Per patient, per-cell D for malignant cells minus normal
   epithelium, pooled across studies, with an equivalence bound set to the D difference that 2 h per-gene shifts
   produce at the real depth.
2. **Weaker clock per cell vs cells out of phase.** The pseudobulk clock amplitude is per-cell amplitude times
   synchrony (lambda R = rho). Cells out of phase keep the within-sample spread D but lose pseudobulk amplitude; a weaker
   per-cell clock loses both. The descriptive WP0 pattern (D equal to normal, pseudobulk coherence lower) hints at
   dephasing but is untested. lambda and R rest on Assumption W (V11): state it.
3. **The same split by tumour proliferation** (the bulk loss is a proliferation effect).
4. **Clones (exploratory):** per-clone pseudobulks from copy-number-inferred clones, and whether clones in one tumour
   share a phase.

## Step 1: scout and power, no outcome read (do this first, then report)
- Inventory the single-cell tumour studies already on disk (A138/A149/A187 inputs) and the Curated Cancer Cell Atlas
  (3CA, Tirosh lab; processed matrices with malignant annotation and inferred copy number); add CELLxGENE Census
  tumour data if useful. Verify sizes and formats yourself; nothing about 3CA's contents was checked.
- Eligibility to propose: droplet (UMI) data, patients with >= 100 malignant and >= 100 normal epithelial cells of
  the same lineage, median clock UMI per cell >= 2 (the METHOD.md eligibility rule).
- Report per study: patients eligible, cells, clock UMI per cell, tissue, presence of normal epithelium, cell-cycle
  annotation possible, and download size. Downloads: processed matrices only, `$SCCC_DATA/<source>/` with a
  `PROVENANCE.md` (URL, date, sha256), parse with `python -I`, never execute; the 500 GB cap holds.
- Power and identifiability by simulation at the real design (clocksync `simulate` or the WP0 calibrated simulator):
  the MDE of the pooled D difference against 2 h rewiring; whether (D, pseudobulk amplitude) separates damping from
  dephasing at the real depth and patient number. No real malignant-vs-normal D difference or amplitude is computed
  in this step.
- Then bring the user the inventory and power table, and the proposed registration. Register (A214 onwards) only what
  is powered: a failed gate or an unpowered test is reported as not estimable, never as "no difference".

## Rules that bite here
- Pre-register each hypothesis in `registered_experiments.md` and commit before scoring; anything after seeing results
  is post hoc. Results to `results_preliminary.md` (end each entry with a `<!-- src: -->` line) and `VERDICT_TABLE.md`.
- Every contrast needs a specificity companion on the same permutations (expression-matched pseudo-clock sets, and
  responsive sets where relevant; companion rule at the 5th percentile).
- Compare coupling only at a common d* or on z2/psi. MDE = 2.49 x null SD.
- Numbers in text come from output tables. Verify subagent numbers against the files.
- Anything outward-facing (emails to data owners, access requests) is drafted and brought to the user first.

## State at handover
All work is committed (last commit dd2aece, then this prompt). No jobs or subagents are running. Open user-decision
items elsewhere (not part of this round): genome-burden test beyond A213, Masri GSE73222, Wu 2020 skin/glioblastoma,
ZeitZeiger, Heyde adipose vs liver feeding, Oscope artefact control, and a pure no-clock arm for the gate false-claim
rate. That arm is the A209b design with clock-plane noise psi = 0 in every compartment, so clock and pseudo genes are
exchangeable. A209b's exchangeability check failed because each compartment's psi is clock-shaped by construction, so no
unbiased false-claim rate exists yet for gate + matched companion. It needs a registration, and runs in about 15 min with
`A209_gate_companion/scripts/compute_a209.py`.
