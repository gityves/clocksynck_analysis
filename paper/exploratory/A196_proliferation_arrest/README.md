# A196: does arresting proliferation restore tumour clock coherence?

Registration: `paper/docs/registered_experiments.md`, A196 (commit a0dc0c5). Outcome: `paper/docs/results_preliminary.md`.
Data: GSE93204 (NeoPalAna; series matrix in `$SCCC_DATA/scout_endocrine`, provenance there) and the cached A168 tables
(POETIC, GSE20181, GSE87411; run `A168_endocrine/scripts/prep_a168.py` first).

    python -I scripts/prep_a196.py                   # GSE93204 -> output/neopal_{meta,expr}.parquet
    python -I scripts/score_a196.py [--power-only]   # about 40 s

Outputs (`output/`, not versioned): `power.csv`, `tests.csv`, `descriptive.csv`, `calibration_posthoc.csv`, `prep_log.csv`.
Title variants handled in prep: one `C1D15_Redo` (taken as C1D15, no duplicate) and one `S` (surgery).
