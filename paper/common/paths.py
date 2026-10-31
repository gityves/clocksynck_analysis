"""Paths of the paper code. Data roots come from environment variables with the lab machine's defaults.

Analyses are grouped (a group = one folder of compute scripts sharing helpers); each group sits in the figure whose
claim it supports: `results/<figure>/scripts/<group>/`, outputs in `results/<figure>/output/<group>/`.
"""

from __future__ import annotations

import os
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
REPO = PAPER.parent
RESULTS = PAPER / "results"
DATASETS = PAPER / "common" / "datasets"  # builders: datasets/<name>/build.py -> standard pseudobulk tables
PAPER_DATA = PAPER / "data"  # small data kept in the repository
GTEX_PHENOTYPES = PAPER_DATA / "Gtex_phenotypes"  # controlled access: aggregates only, never per-donor rows
DATA = Path(os.environ.get("SCCC_DATA", "/mnt/max8tb/data/CoPhaser_data/cancer_circadian"))  # raw and built data
PSEUDOBULK = DATA / "pseudobulk"  # standard pseudobulk tables
# Archive of the project history (read only): legacy tables that some builders still convert, legacy outputs that
# the checks compare against, the pre-registrations A1-A108c (docs/PLAN_internal_desync.md).
ARCHIVE = Path(os.environ.get("SCCC_ARCHIVE", "/home/maxine/Documents/paychere/scCircadianCancer"))
LEGACY = ARCHIVE / "legacy"

GROUPS = {  # group -> figure folder under results/
    "method": "figure_1", "windows": "figure_1", "species": "figure_1",
    "organs": "figure_2", "backbone": "figure_2", "embryo": "figure_2", "atlas_fits": "figure_2",
    "gene_rhythms": "figure_2",
    "cues_lost": "figure_3", "gtex_cues": "figure_3",
    "resilience": "figure_4", "spread_multiplicity": "figure_4", "natural_experiments": "figure_4",
    "cancer": "figure_5", "cancer_pipeline": "figure_5", "lung_breaks": "figure_5", "psoriasis": "figure_5",
    "psoriasis_details": "figure_5",
    "calibration": "extended_data", "gate": "extended_data", "host_identity": "extended_data", "spikeins": "extended_data",
    "simulations": "extended_data",
}


def figure_of(group: str) -> Path:
    if group not in GROUPS:
        raise KeyError(f"unknown analysis group {group!r}; add it to common.paths.GROUPS")
    return RESULTS / GROUPS[group]


def group_dir(group: str) -> Path:
    """Code folder of a group: results/<figure>/scripts/<group>/."""
    return figure_of(group) / "scripts" / group


def group_out(group: str) -> Path:
    """Output folder of a group (not created): results/<figure>/output/<group>/."""
    return figure_of(group) / "output" / group


def out_dir(group: str, subsample: bool = False) -> Path:
    """results/<figure>/output/<group>/ (or .../<group>/subsample/ for checks on subsamples); created."""
    p = group_out(group) / ("subsample" if subsample else "")
    p.mkdir(parents=True, exist_ok=True)
    return p


def source_data(group: str) -> Path:
    p = group_out(group) / "source_data"
    p.mkdir(parents=True, exist_ok=True)
    return p


def exploratory_dir(amendment: str, subsample: bool = False) -> Path:
    """Outputs of a registered hypothesis not yet validated for a figure: exploratory/<amendment>/output/; created.

    `amendment` is the folder name, e.g. "A109_coupling_nulls"."""
    p = PAPER / "exploratory" / amendment / "output" / ("subsample" if subsample else "")
    p.mkdir(parents=True, exist_ok=True)
    return p
