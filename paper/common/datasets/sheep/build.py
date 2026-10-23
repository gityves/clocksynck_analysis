"""Sheep, 343 animals x 45 tissues and cell types (SheepGTEx v0, untimed): tissues as compartments of one animal.

Source: `$SCCC_DATA/farm_gtex/sheep/raw/` (`SheepGTEx_v0.Gene.count.txt.gz`, 25,699 genes x 8,171 passing RNA
samples with NCBI symbols; `SheepGTEx_v0.MetaData.xlsx` Table S1 gives Individual, tissue / cell subtype, BioProject,
stage, sex, breed). The authors assign samples to animals (`Individual`). No sampling times exist.
Legacy builder replaced: `desync/validation/farm_gtex/scripts/build_sheep.py` (A34; checks against its
`output/sheep_clock.csv`).
Used by: A34 (fit_sheep, match to GTEx tissues), A41 (sheep with the baboon waveform), A42 (brain vs periphery,
couplings across species).

Choices (as legacy, fixed in farm_gtex/method.md before fitting): samples with Filter == "pass"; compartment = tissue /
cell subtype with skeletal-muscle subtypes merged into "Muscle"; PBMC, alveolar macrophage, gametes, embryo and
extra-embryonic tissues, milk, lymph, hair follicle and unspecified "Brain" are dropped, as is every prenatal sample;
replicate libraries of an animal x compartment are summed (`n_lib`); depth = genome-wide count (sum over all genes of the
matrix); stratum = BioProject (the pooled submissions of the same animals; one per animal x compartment, first
value); inclusion rules iterated to a fixed point: (compartment, BioProject) groups of at least 3 animals, animals
with at least 3 compartments, compartments with at least 20 animals (343 animals, 45 compartments, 38 BioProjects).
The 14 clock genes are present under their human symbols (ARNTL for BMAL1). Donors table: stage, sex, breed are
descriptive (animal-level).
Note: `stratum` is per (animal x compartment) as in the legacy fits (124 of 343 animals have compartments from
several BioProjects); `clocksync.validate` requires one stratum per donor, so the sheep fit needs per-row strata
(see INDEX_bulk.md).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import clocksync as cs
from common.paths import DATA

RAW = DATA / "farm_gtex" / "sheep" / "raw"
MIN_ANIMALS, MIN_CELL, MIN_TISSUES = 20, 3, 3
MERGE = {"Longissimus muscle": "Muscle", "Biceps muscle": "Muscle", "Infraspinatus muscle": "Muscle",
         "Skeletal muscle": "Muscle"}
NOT_TISSUE = {"PBMC", "Alveolar macrophage", "Oocyte", "Granulosa cell", "Embryo", "Milk", "Lymph", "Amnion",
              "Placentome", "Hair follicle", "Brain"}


def samples() -> pd.DataFrame:
    """One row per RNA sample: clock-gene counts, genome-wide depth and the Table S1 annotation."""
    genes = list(cs.gtex().genes)
    meta = pd.read_excel(RAW / "SheepGTEx_v0.MetaData.xlsx", sheet_name="Table S1", header=1)
    meta = meta[meta.Filter == "pass"].set_index("Sample")
    depth, clock = {}, {}
    for ch in pd.read_csv(RAW / "SheepGTEx_v0.Gene.count.txt.gz", sep="\t", index_col=0, chunksize=2000):
        ch = ch.astype(np.float64)
        s = ch.sum(0)
        for c in s.index:
            depth[c] = depth.get(c, 0.0) + s[c]
        for g in ch.index.intersection(genes):
            clock[g] = ch.loc[g]
    cnt = pd.DataFrame(clock)[genes]
    smp = cnt.join(pd.Series(depth, name="depth"))
    return smp.join(meta[["Individual", "BioProject", "Tissue/cell type", "Tissue/cell subtype",
                          "Developmental stage (corrected)", "Sex (corrected)", "Breed (corrected)"]], how="left")


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pb, donors): animal x compartment pseudobulk after the A34 inclusion rules; donors indexed by animal."""
    genes = list(cs.gtex().genes)
    smp = samples()
    s = smp[smp.Individual.notna()].copy()
    s["compartment"] = s["Tissue/cell subtype"].replace(MERGE)
    s = s[~s.compartment.isin(NOT_TISSUE) & (s["Developmental stage (corrected)"] != "Prenatal")]
    agg = {g: "sum" for g in genes} | {"depth": "sum", "BioProject": "first", "Developmental stage (corrected)": "first",
                                       "Sex (corrected)": "first", "Breed (corrected)": "first"}
    s["n_lib"] = 1
    agg["n_lib"] = "sum"
    t = s.groupby(["Individual", "compartment"]).agg(agg).reset_index()
    for _ in range(10):
        n0 = len(t)
        t = t[t.groupby(["compartment", "BioProject"]).Individual.transform("nunique") >= MIN_CELL]
        t = t[t.groupby("Individual").compartment.transform("nunique") >= MIN_TISSUES]
        t = t[t.groupby("compartment").Individual.transform("nunique") >= MIN_ANIMALS]
        if len(t) == n0:
            break
    t = t.rename(columns={"Individual": "donor", "BioProject": "stratum"})
    t["donor"], t["stratum"] = t.donor.astype(str), t.stratum.astype(str)
    donors = (t.drop_duplicates("donor").set_index("donor")
              [["Developmental stage (corrected)", "Sex (corrected)", "Breed (corrected)"]]
              .rename(columns={"Developmental stage (corrected)": "stage", "Sex (corrected)": "sex",
                               "Breed (corrected)": "breed"}))
    pb = t[["donor", "compartment", "stratum", "n_lib", "depth"] + genes].reset_index(drop=True)
    return pb, donors
