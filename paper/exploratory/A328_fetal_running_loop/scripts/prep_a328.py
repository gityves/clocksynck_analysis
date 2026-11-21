"""A328 part a preparation: del Valle 2023 fetal adrenal sheet and cohort clock-gene depth (`A328 note 1`).

Reads the delivered exon counts of the 32 adrenal runs and writes only what the power stage may read: one row per fetus
(stage, sex, library size = sum of exon counts) and the cohort median count of each clock gene. No per-fetus
clock-gene value is written or printed.

    python prep_a328.py
Outputs (`../output/prep/`): dv_sheet.csv, dv_clock_detection.csv.
"""
from pathlib import Path

import pandas as pd

from clocksync.genes import CLOCK_GENES, canonical

SRC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/sccc_mapping/delvalle2023_fetal_adrenal")
OUT = Path(__file__).resolve().parents[1] / "output" / "prep"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(SRC / "samples.tsv", sep="\t")
    s = s[s.sdrf_Characteristics_organism_part == "adrenal gland"]
    exon = pd.read_csv(SRC / "counts_exon.tsv.gz", sep="\t", index_col=0)
    sym = exon.pop("symbol").astype(str).map(canonical)
    runs = s.run_accession.tolist()
    sheet = pd.DataFrame(dict(donor=s.sample_title.str.replace(" ", "_").to_numpy(), run=runs, compartment="adrenal",
                              stage=s.sdrf_Characteristics_age.astype(float).to_numpy(),
                              sex=s.sdrf_Characteristics_sex.to_numpy(), depth=exon[runs].sum().to_numpy()))
    assert sheet.donor.is_unique and len(sheet) == 32
    clock = exon.loc[sym.isin(CLOCK_GENES).to_numpy(), runs]
    clock.index = sym[sym.isin(CLOCK_GENES)].to_numpy()
    clock = clock.groupby(level=0).sum().reindex(list(CLOCK_GENES))
    assert clock.notna().all().all()
    det = clock.median(axis=1).rename("median_count").to_frame().T
    det.index = ["delvalle"]
    sheet.to_csv(OUT / "dv_sheet.csv", index=False)
    det.to_csv(OUT / "dv_clock_detection.csv")
    print(sheet.groupby(["stage", "sex"]).size().to_string(), f"median library size {sheet.depth.median():.0f}",
          sep="\n")


if __name__ == "__main__":
    main()
