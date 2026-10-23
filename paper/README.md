# Paper code

Analyses of the paper, on the public API of [`clocksync`](../clocksync). Layout, migration status and working rules:
`docs/provisory_plan_paper.md`, section "Repository and migration".

```
common/        paths (read paths.py first), CLI, dataset loading, GTEx access, frozen A17 anchors, benchmark methods,
               run.sh (RAM-capped runner); datasets/<name>/build.py (builders -> $SCCC_DATA/pseudobulk/)
results/       <figure>/scripts/<group>/ (validated analyses), <figure>/output/<group>/ (tables), <figure>/data/
exploratory/   <A1xx_slug>/ per registered hypothesis not yet in a figure
data/          small data (controlled-access GTEx phenotypes are git-ignored)
docs/          plan, registered experiments, results, method, amendment index
```

Run a group: `results/<figure>/scripts/<group>/run.sh [--subsample N]` (full run by default). A heavy step:
`[RUN_CAP_GB=30] common/run.sh <script> [args...]`. Python: `$(conda info --base)/envs/sccc/bin/python`. Data roots:
`$SCCC_DATA` (default `/mnt/max8tb/data/CoPhaser_data/cancer_circadian`), archive `$SCCC_ARCHIVE` (default
`/home/maxine/Documents/paychere/scCircadianCancer`, read only).

GTEx time of death and protected phenotypes are controlled access: outputs hold aggregates only.
