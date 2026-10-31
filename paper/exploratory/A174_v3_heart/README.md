# A174: V3 heart (GSE270258) port, builder only

Port of the archive validation V3 (`legacy/scripts/v3_heart.py`). This folder holds the data build; nothing is scored here
(gate, host fit and AM/PM tests wait for the committed re-scoring note).

- **Data:** Jia 2024, left ventricle (end-stage DCM) of 16 patients, one snRNA-seq library each, aortic cross-clamp at 06-10 h
  (AM) or 14-18 h (PM); arms NCR (normal rhythm, primary) and ACR (sleep apnoea, selected on clock genes), 4 patients per
  arm x window. Source `$SCCC_DATA/validation/GSE270258/` (raw-count h5ad, `samples.csv`, `PROVENANCE.md`). Only the window is
  known: no per-patient surgery time exists (`exact_crossclamp_time` is empty), so the window midpoints (8 h, 16 h) are the
  only time information.
- **Rule:** no author labels, and the COARSE rules match Census cell-type names, so cells are called by marker scores
  (the legacy rule: mean log1p CP10K over the markers of a class; best > 0.3 and second < 0.7 x best, else dropped).
  Markers, thresholds and `MIN_CELLS` are constants at the top of `scripts/build_v3.py`; no clock gene is a marker. Classes follow
  `common/datasets/_atlas/heart.py::CLASSES`, compartment `heart:<class>`. Differences from the legacy: pericyte + smooth muscle
  merged into `mural`; pseudobulks >= 20 cells (legacy 50; filter `n_cells`).
- **Run:** `RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A174_v3_heart/scripts/build_v3.py` (36 s, 3.7 GB).
- **Outputs** (`output/`, git-ignored): `v3_heart.parquet` (donor, compartment, n_cells, depth, 14 clock-gene counts; rule stamp
  `r1-...` in the parquet metadata under `common.data.RULES_KEY`), `v3_heart_donors.parquet` (group, window, window start/end/mid
  h, sex, age), `cells_per_class.csv` (cells per patient x class, incl. `unassigned`), `labels.csv.gz` (per-cell class).
  The table is not in the shared `PSEUDOBULK` folder (the stamp scheme needs a builder under `common/datasets/`); move the
  builder there when the dataset enters a figure.
