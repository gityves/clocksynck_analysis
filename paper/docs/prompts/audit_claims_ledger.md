# Prompt: independent audit of the claims ledger

You are auditing a scientific paper's claims ledger against the project's recorded results. You did not write the ledger
and should not assume it is right. Your job is to find what is missing, wrong, overstated or mislabelled.

## Project in two lines

`clocksync` infers circadian clock relations between organs or cell types of the same individual (donor x compartment
pseudobulks) without sample times, so it works on samples taken in hospital hours. The paper aims at Cell; its draft
headline is that organs keep a shared internal time across development, ageing, loss of daily cues and chronic disease,
and that clock breaks are local (brain, psoriatic skin, proliferating tumours). Read `CLAUDE.md` at the repository root
first (rules, data access, environment).

## Files (repository root: `/home/maxine/Documents/paychere/scCircadianOverview`)

- **Under audit:** `paper/docs/CLAIMS_LEDGER_2026-11-03.md`.
- **Results, primary sources:**
  - `paper/docs/VERDICT_TABLE.md`: one row per registered test from A109 onward, with the output file of each.
  - `paper/docs/results_preliminary.md`: every scored amendment, with errata and post hoc notes. Each entry ends with a
    `<!-- src: ... -->` line naming its output tables.
  - `paper/docs/registered_experiments.md`: the analysis plans (what was registered as primary, secondary, stop rules).
- **Results, secondary sources** (summaries; check them against the primary ones, do not trust them):
  `paper/docs/FINDINGS_MAP_2026-11-03.md`, `paper/docs/STATE_2026-10-31.md`, `paper/docs/HANDOFF_2026-10-31.md`,
  `paper/docs/METHOD.md`, `paper/docs/provisory_plan_paper.md` (figure plan), `paper/docs/literature_positioning.md`,
  `paper/docs/PI_QUESTIONS.md`.
- **Archive (A1-A108c):** `/home/maxine/Documents/paychere/scCircadianCancer/docs/RESULTS_preliminary.md` and
  `PLAN_internal_desync.md`; `paper/docs/AMENDMENT_INDEX.md` maps archive amendments to migrated scripts. Archive results
  that the figures still use count as results.
- **Output tables:** `paper/exploratory/<amendment folder>/output/` and `paper/results/figure_*/output/`.

## Tasks

1. **Coverage.** List every scored result that bears on the paper: all rows of `VERDICT_TABLE.md`, every scored entry of
   `results_preliminary.md` (including post hoc and descriptive ones), and the archive results that the figure plan or
   `FINDINGS_MAP` cites. For each, say whether the ledger (a) uses it, (b) lists it under "Do not claim", or (c) omits
   it. For every (c), judge whether it matters: does it support, qualify or contradict a ledger row, or would it make a
   claim of its own? Pay particular attention to results that contradict or weaken a ledger row.
2. **Numbers.** For every number in the ledger, find it in the output table named by the `src` line or the verdict table,
   and report any mismatch (value, interval, p, Holm p, n, direction). Read the CSV, not the summary docs. Where the ledger
   rounds, accept rounding to the stated digits.
3. **Labels.** Check each row's provenance and status against the docs: registered primary vs secondary vs post hoc vs
   descriptive vs archive; failed gates; not estimable; "seen data"; errata that changed a verdict or a reading (search
   `results_preliminary.md` for "erratum", "superseded", "withdrawn", "post hoc"). Flag any row whose status (A / B / C) is
   higher than its evidence.
4. **Wording.** Check each "Sentence" against the ledger's own six rules (bounds not absence; scope; the definition of
   internal desynchrony as a phase offset between compartments; positive control next to nulls; provenance; convergence
   not refutation). Flag overreach and propose tighter wording. Also flag sentences that undersell a result the evidence
   supports.
5. **Objections.** For each section, add any strong reviewer objection the ledger does not anticipate, with whether the
   recorded results answer it (cite the result) or not.
6. **Headline.** Say whether the headline (section 0) follows from the A and B rows, which rows it actually depends on,
   and what single missing result would most change it.

## Rules

- **Read only.** Do not edit any existing file, do not commit, do not run analyses longer than a few minutes, do not
  download anything. Reading CSVs with short Python one-liners is fine
  (`$(conda info --base)/envs/sccc/bin/python`).
- **Controlled access.** GTEx time of death and GTEx protected phenotypes are controlled access: never print or copy
  per-donor rows; aggregates only.
- **Evidence for every finding.** Each finding cites a file and line (or CSV path and row). If you cannot verify
  something, say "unverified" instead of guessing. Do not infer results from file names.
- **No recommendations about journals, authorship or contacting other groups**; those are the PI's decisions.

## Output

Write one file, `paper/docs/AUDIT_CLAIMS_LEDGER_<YYYY-MM-DD>.md`, with:
1. a summary of at most 15 lines: the most important omissions, errors and overstatements, ranked by how much they
   would change the paper;
2. the coverage table (result id | what it shows, one line | in ledger: used / do-not-claim / omitted | matters? why);
3. number mismatches (ledger row | ledger value | source value | file);
4. label and status problems;
5. wording problems with proposed text;
6. missing objections;
7. the headline assessment.

Keep it concise and factual. Return the file path and the summary when done.
