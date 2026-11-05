# Handoff (2026-11-05): literature follow-up analyses (`prompts/literature_followup_analyses.md`)

## State
- **A204 registered** (commit 70e3564), code committed (e5e2473). Item 1 of the prompt (COFE / TimeTeller mechanism).
- **Running:** `compute_a204.py --jobs 8` through run.sh (runner PID 2394297, log `paper/results/logs/compute_a204*.log`),
  630 tasks of about 30 s; started 23:25, expected about 40-60 min. A detached waiter (PID 2396471) runs
  `score_a204.py` when it ends -> `A204_orderer_uncoupling/output/score.log` and `{gate,tests,cofe_cells,auc,tt_cells}.csv`.
- If the run failed (`paper/common/jobs.sh` shows `!!`), rerun the same command: tasks are cached and resume.
- **Not yet done:** read gate.csv first (if it fails, H204 is "not scored"); write results to `results_preliminary.md`
  (with `<!-- src: -->`), `VERDICT_TABLE.md`, a `METHOD.md` row if relevant; report to the user per the prompt
  (quote the cards: COFE "distributed throughout the circadian cycle", TimeTeller "slope -0.027, p = 0.80", both checked
  in the full texts on 2026-11-05).

## Checks made before registering (rule 1)
- A138 Q2 H138.3 already registered the Rayleigh readout on a genome-scale simulator; parked after its simulator gate
  failed (ROADMAP 4.8). A204 does not score it; it uses the A145 clock-gene simulator.
- COFE's paper projects matched healthy samples on the tumour-trained model (test samples); A204 does the same.
- clocksync needs >= 2 anchors for anchored fits; the tumour-normal readout is the A138 H138.4 pairing statistic.

## Next items (prompt order), none registered yet
2. GTEx brain vs periphery arg z (read A142/A143 first; check brain amplitude and power).
3. Organ offset vs proliferation (Droin 2019).
4. Gate negative controls (gene-specific partial shifts; Cry1/2 DKO GSE13093).
5. Dataset scouting (ask the user before downloads that matter). 6-7 as in the prompt; ask the user about A117.
