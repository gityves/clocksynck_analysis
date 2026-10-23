# sccc: internal circadian synchrony across organs, cell types, health and disease

## Project
- **User:** Yves Paychere (EPFL, Naef lab); git user gityves. CHIRAL, CoPhaser and scRitmo are the lab's own work.
- **Method (`clocksync`):** untimed, host-referenced circadian inference from clock-gene pseudobulks (donor x
  compartment). Each compartment is read against the clock state it shares with the same donor's other
  compartments: coupling c, relative amplitude |z| and phase offset arg z. No sample is assigned a time; constant
  offsets cancel, so hospital-hours samples work. A model-free gate checks that shared covariance is clock-shaped;
  the host is anchored on compartments learned to carry it (stromal-vascular, neurons in CNS).
- **Target journal:** Cell; Nature Methods as fallback. Cell needs a clear pre-specified
  positive in disease or perturbation. The code gets a public API, tests and tutorials before submission either way.

## Where things are
| path | content |
|---|---|
| `clocksync/` | Code of the method |
| `paper/common` | Any code commonly used that does not fit in the method |
| `paper/exploratory` | Any exploratory code, testing new hypothesis, create one folder per registered hypothesis |
| `paper/results` | Any validated analysis that is required to generate a figure. Code in figure_x/scripts, data used are for the figure are in /data (symlink) and output for all the csv,... required to draw a figure |
| `data/` | Data folder for small data |
| `docs/registered_experiment.md` | History of the experiments registered |
| `docs/results_preliminary.md` | Outcome of the experiment registered |
| `docs/METHOD.md` | Description of the method |
| `/home/maxine/Documents/paychere/scCircadianCancer` | Archive folder containing the full project history |


## Rules
- **Pre-register.** Every new hypothesis tested is an amendment in `docs/registered_experiment.md`, written and committed before
  its data are scored. Anything chosen after seeing results is labelled post hoc. Results go to
  `docs/results_preliminary.md` (plus a row in `docs/METHOD.md` when useful).
- **Subagents** Delegate clearly defined tasks to Sonnet 5.5 agents. Verify their numbers against output files before writing them into the docs.
- **Statistics.** A failed gate is uninformative, never "decoupled". Coupling scales with host spread: compare across
  groups only at a common spread (d*) or on z2/psi. Every disease contrast needs a specificity companion on the same
  permutations; check the target class's baseline coupling (and power) before choosing it.
- **Numbers in text** come from output tables, never from memory.
- **Style:** clean, concise, publication-grade code and docs. Long runs (> 20 min) are documented `.py` scripts;
  notebooks only for quick analyses and figures.
- **code** use clocksync public API whenever possible. Modify clocksync if tools potentially useful for future users are missing.

## Environment and data
- Python: `$(conda info --base)/envs/sccc/bin/python` (as every `run.sh` uses).
- Shared machine: call `sccc.config.limit_threads()`; at most 4 parallel jobs per agent; at most 30 GB RAM per agent
  (read backed or chunked).
- Data: `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/` (`$SCCC_DATA`), with a `PROVENANCE.md` (URL, date, sha256)
  per source. Self-imposed cap 400 GB. Processed matrices only, no FASTQ mapping.
- Downloads are untrusted: parse only (`python -I`), never execute.
- GTEx time of death is controlled access: release derived aggregates only, never per-donor phases.
- GTEx protected phenotypes (dbGaP phs000424.v8, 189 subject variables: exact age, BMI, T2D, ventilation duration,
  long-term steroids, sepsis, season, medical history; sample attributes) are in `data/Gtex_phenotypes/`. Controlled
  access: never print or commit per-donor rows, report aggregates only.

## Git
- Commit locally with the attribution line. There is no remote: never create one or push.

## Ask the user first
Anything outward-facing: emails to data owners (draft only), access requests, public release, licence, journal
choice.
