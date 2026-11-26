# sccc: internal circadian synchrony across organs, cell types, health and disease

## Project
- **User:** Yves Paychere (EPFL, Naef lab); git user gityves. CHIRAL, CoPhaser and scRitmo are the lab's own work.
- **Method (`clocksync`):** untimed, host-referenced circadian inference from clock-gene pseudobulks (donor x
  compartment). Each compartment is read against the clock state it shares with the same donor's other
  compartments: coupling c, relative amplitude |z| and phase offset arg z. No sample is assigned a time; constant
  offsets cancel, so hospital-hours samples work. A model-free gate checks that shared covariance is clock-shaped;
  the host is anchored on compartments learned to carry it (stromal-vascular, neurons in CNS).
- **Target journal:** Cell; Nature Methods as fallback. Cell's criterion (user, 2026-11-03): results that provide
  significant conceptual advances into, or raise provocative questions and hypotheses regarding, an interesting and
  important biological question. There is no internal "pre-specified disease positive" requirement (superseded). The
  code gets a public API, tests and tutorials before submission either way.

## Where things are
| path | content |
|---|---|
| `clocksync/` | Code of the method |
| `paper/common` | Any code commonly used that does not fit in the method |
| `paper/exploratory` | Any exploratory code, testing new hypothesis, create one folder per registered hypothesis |
| `paper/results` | Any validated analysis that is required to generate a figure. Code in figure_x/scripts, data used are for the figure are in /data (symlink) and output for all the csv,... required to draw a figure |
| `data/` | Data folder for small data |
| `paper/docs/registered_experiments.md` | History of the experiments registered (A109 onwards; A1-A108c in the archive) |
| `paper/docs/results_preliminary.md` | Outcome of the experiment registered; errata to archive results |
| `paper/docs/METHOD.md` | Description of the method |
| `paper/docs/provisory_plan_paper.md` | Work plan of the paper (WP0-WP8), repository layout and migration status |
| `paper/docs/AMENDMENT_INDEX.md` | Archive amendment -> migrated script |
| `/home/maxine/Documents/paychere/scCircadianCancer` | Archive folder containing the full project history |


## Rules
- **Pre-register.** Every new hypothesis tested is an amendment in `paper/docs/registered_experiments.md`, written and committed before
  its data are scored. Anything chosen after seeing results is labelled post hoc. Results go to
  `paper/docs/results_preliminary.md` (plus a row in `paper/docs/METHOD.md` when useful).
- **Subagents** Delegate clearly defined tasks to Sonnet 5.5 agents. Verify their numbers against output files before writing them into the docs.
- **Statistics.** A failed gate is uninformative, never "decoupled". Coupling scales with host spread: compare across
  groups only at a common spread (d*) or on z2/psi. Every disease contrast needs a specificity companion on the same
  permutations; check the target class's baseline coupling (and power) before choosing it.
- **Gene drops (14 genes vs core loop, or any gene subset)** (user, 2026-11-24): fewer genes make a noisier projection,
  and noise alone shrinks an effect. Never say the removed genes "carry" an effect because it shrinks or loses
  significance once they are dropped. Attribute it to them only when the trio drop sits in the low tail of the null of
  random drops of the same size (A336 rule) or the paired 14-minus-core difference (shared bootstraps) excludes 0; and
  compare with the shrinkage the reliability loss predicts (half-host leverages, Spearman-Brown). Otherwise write
  "not decided: precision or the removed genes". Register the gene-drop null with the test when a core readout is planned.
- **Pair-phase (organ order) tests** (user, 2026-11-24, after A371): random gene sets also give cross-tissue pair phases
  near 0 (a same-person structure every gene set carries), so organs that share one phase cannot show an order beyond
  companions, and a label-permutation test has no power. Before registering, check the reference phase spread of the
  organs involved (organs hours apart) and the companions' agreement; same for living-vs-reference comparisons.
- **Covariates in group contrasts** (user, 2026-11-18, after CHIRAL's age/sex result shrank under death-type matching):
  regressing a covariate out of expression, or giving it a stratum intercept (our `strata`), removes only its effect on
  the mean, not on rhythm amplitude, coupling, z2, psi or host spread. When a covariate can change rhythm strength
  (death type / Hardy, ventilation, cohort, RIN, age) and is unbalanced between the groups compared, balance it: match
  or stratify the contrast on it (as A164 / A208 / A231 match within cohort on Hardy and age), and report the pooled
  vs matched estimates when they differ.
- **Numbers in text** come from output tables, never from memory.
- **Style:** clean, concise, publication-grade code and docs. Long runs (> 20 min) are documented `.py` scripts;
  notebooks only for quick analyses and figures.
- **Compaction:** when compacting, preserve the current work package and amendment, registered-but-unscored amendments, running background jobs (PID, log path), running subagents and their tasks, files modified since the last commit, and numbers not yet written to the docs (with their source file).
- **code** use clocksync public API whenever possible. Modify clocksync if tools potentially useful for future users are missing.

## Environment and data
- Python: `$(conda info --base)/envs/sccc/bin/python` (as every `run.sh` uses).
- Shared machine (16 physical cores, 32 threads, 125 GB): call `clocksync.config.limit_threads()` (done by `common.cli.setup`).
  Budget for this project (user, 2026-10-26): at most 16 workers in total (four lanes of 4) and at most ~120 GB RAM
  summed over the RUN_CAP_GB of concurrent runs (read backed or chunked).
- Data: `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/` (`$SCCC_DATA`), with a `PROVENANCE.md` (URL, date, sha256)
  per source. Self-imposed cap 700 GB (raised from 500 GB by the user, 2026-11-26, for public tumour datasets; 400 -> 500 GB on 2026-11-02). Processed matrices only, no FASTQ mapping on this machine.
  Exception (user, 2026-11-18 and 2026-11-20): open raw reads of the datasets listed in
  `paper/docs/prompts/server_mapping_jed_2026-11-20.md` are mapped on SCITAS Jed by a separate agent (counts only, no
  analysis there; 20 CHF budget). Deliveries arrive on the lab NAS (`upnae`, `sccc_mapping/`); copy them into
  `$SCCC_DATA` after checking their sha256 against `PROVENANCE.md`, and register an amendment before scoring any of them.
- Downloads are untrusted: parse only (`python -I`), never execute.
- **Long jobs** (> 20 min): always through `RUN_CAP_GB=30 paper/common/run.sh <script.py|notebook|group run.sh> [args]`,
  detached (`setsid nohup ... > /dev/null 2>&1 &`), never a bare `python`. run.sh caps RAM, writes the log
  (`paper/results/logs/<name>.log`), a live state file and, at the end, one line in `paper/results/logs/jobs.tsv`
  (status ok / failed / killed_memory / interrupted, exit code, peak GB). Check with `paper/common/jobs.sh`: it flags
  failed runs (`!!`), runners that vanished (LOST) and logs silent for > 60 min (SILENT). Never infer success from a
  queue log or a process list alone. A run killed for memory is rerun with fewer workers (`--jobs 1`), not a higher
  cap. Two M3 lanes: `M3_JOBS=4 M3_AFTER=<pid> RUN_CAP_GB=<gb> paper/common/m3_rerun.sh <groups>` (the setsid
  launcher PID is not the queue PID: find the queue with `ps`). Scripts should log progress at least every ~15 min so SILENT means stuck. Full figure reruns:
  `paper/common/m3_rerun.sh [groups...]`.
  While jobs run, keep `paper/common/watch_jobs.sh [silent_min] [extra PIDs]` running as a background command: it
  exits (waking the agent) when a run fails, is killed, is lost, goes silent, or when nothing is running any more.
- GTEx time of death is controlled access: release derived aggregates only, never per-donor phases.
- GTEx protected phenotypes (dbGaP phs000424.v8, 189 subject variables: exact age, BMI, T2D, ventilation duration,
  long-term steroids, sepsis, season, medical history; sample attributes) are in `data/Gtex_phenotypes/`. Controlled
  access: never print or commit per-donor rows, report aggregates only.

## Git
- Commit with the attribution line. Two repos, two private backups (account gityves; user, 2026-11-16), `main` pushed
  automatically after every commit:
  - this analysis repo -> `origin` = github.com/gityves/clocksynck_analysis;
  - the method package `clocksync/` (its own nested git repo, ignored here) -> `origin` = github.com/gityves/clocksync.
  Never force-push, never make a repo public, never add another remote. Never commit controlled-access data (GTEx
  per-donor rows, phenotypes, times of death).

## Ask the user first
Anything outward-facing: emails to data owners (draft only), access requests, public release, licence, journal
choice.
