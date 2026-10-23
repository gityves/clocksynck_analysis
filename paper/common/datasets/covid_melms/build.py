"""Melms et al. 2021 (Nature): lung snRNA of 19 COVID-19 autopsies and 7 surgical controls.

Source:   GEO GSE171524 (per-sample raw count csv.gz, genes x cells; author metadata `GSE171524_lung_metaData.txt.gz`),
          `$SCCC_DATA/covid/GSE171524/{ext,raw}`.
Replaces: legacy/desync/validation/natural_experiments/scripts/build_melms.py (+ pbcommon.py), build_covid.md.
Used by:  A105, A105c (post hoc).

Choices (as the legacy builder)
- cells: only barcodes of the author metadata (post-QC nuclei); counts are the raw integer UMI counts of the
  per-sample files (the CellBender `processed_data` file is not needed);
- donor: the add-on sample `L04covaddon` is summed into donor `L04cov`; 19 COVID-19 donors, 7 controls;
- compartment: author `cell_type_main` (level "main", default) or `cell_type_intermediate` (level "intermediate");
- stratum: the assay (one value, 10x 3' v3 nuclei); depth: genome-wide count of the cells of the pseudobulk;
- every pseudobulk with at least one cell is kept (the analyses apply their own minimum of 20 or 50 cells).

Donor table: group (Control / COVID-19), age, sex, intubation_days (IMV), interval_death_symptoms_onset_days, pmi_h,
n_samples. Controls are living surgical patients: no PMI or IMV (confounded with autopsy vs surgery).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(Path(__file__).resolve().parents[3])]
import clocksync as cs  # noqa: E402
from clocksync.genes import gene_columns  # noqa: E402
from common.paths import DATA  # noqa: E402
from natexp_common import GENES, clock_adata, compare, read_legacy, standard  # noqa: E402

ROOT = DATA / "covid" / "GSE171524"
LEVELS = {"main": "cell_type_main", "intermediate": "cell_type_intermediate"}
ASSAY = "10x 3' v3 nuclei"
SAMPLE_COLS = ["group", "age", "sex", "intubation_days", "interval_death_symptoms_onset_days", "pmi_h"]
CHUNK = 2000  # genes per row chunk of the genes x cells csv


def _donor(biosample: str) -> str:
    return re.sub(r"addon$", "", biosample)


def metadata() -> pd.DataFrame:
    """Author cell metadata indexed by barcode, with the donor of every cell."""
    m = pd.read_csv(ROOT / "raw" / "GSE171524_lung_metaData.txt.gz", sep="\t", skiprows=[1], low_memory=False)
    m = m.set_index("NAME")
    m["donor"] = m.biosample_id.map(_donor)
    return m


def read_sample(path: Path, meta: pd.DataFrame) -> ad.AnnData:
    """One genes x cells csv.gz streamed in row chunks: the 14 clock-gene rows and the column sums (cell depth) of
    the cells present in the author metadata."""
    barcodes = pd.read_csv(path, nrows=0, index_col=0).columns.to_numpy()
    symbols = pd.read_csv(path, usecols=[0]).iloc[:, 0].astype(str).to_numpy()
    cols = gene_columns(symbols, GENES)
    if (cols < 0).any():
        raise KeyError(f"{path.name}: clock genes missing")
    depth = np.zeros(len(barcodes))
    X = np.zeros((len(barcodes), len(GENES)))
    r0 = 0
    for df in pd.read_csv(path, index_col=0, chunksize=CHUNK):
        v = df.to_numpy(np.float64)
        depth += v.sum(0)
        for j, c in enumerate(cols):
            if r0 <= c < r0 + len(v):
                X[:, j] = v[c - r0]
        r0 += len(v)
    keep = np.isin(barcodes, meta.index)
    obs = meta.loc[barcodes[keep], ["donor"] + list(LEVELS.values())].astype(str)
    return clock_adata(X[keep], GENES, depth[keep], obs)


def samples_in(ext: Path, n: int | None = None) -> list[Path]:
    files = sorted(ext.glob("GSM*_raw_counts.csv.gz"))
    return files if n is None else files[:n]


def donor_table(meta: pd.DataFrame) -> pd.DataFrame:
    smp = meta.drop_duplicates("biosample_id").set_index("biosample_id")[["donor"] + SAMPLE_COLS]
    don = smp.groupby("donor").agg(**{c: (c, "first") for c in SAMPLE_COLS}, n_samples=("group", "size"))
    don["samples"] = smp.reset_index().groupby("donor").biosample_id.agg("+".join)
    return don.rename_axis("donor")


def build(level: str = "main", donors: list | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pseudobulk, donors). `donors`: restrict to these donors (checks; only their sample files are read)."""
    meta = metadata()
    if donors is not None:
        meta = meta[meta.donor.isin(donors)]
    files = [f for f in samples_in(ROOT / "ext")
             if _donor(re.match(r"GSM\d+_(.+)_raw_counts", f.name).group(1)) in set(meta.donor)]
    cells = [read_sample(f, meta) for f in files]
    pb = cs.pseudobulk(cells, donor="donor", compartment=LEVELS[level], depth="depth", min_cells=1)
    pb["stratum"] = ASSAY
    return standard(pb), donor_table(meta)


def from_legacy(level: str = "main") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Standard table from the legacy `pseudobulk_<level>/pseudobulk.h5ad` and `samples.csv`."""
    d = ROOT / f"pseudobulk_{level}"
    pb = read_legacy(d / "pseudobulk.h5ad")
    pb["stratum"] = ASSAY
    don = pd.read_csv(d / "samples.csv").set_index("donor_id").rename_axis("donor")
    return standard(pb), don


if __name__ == "__main__":
    # check on a slice: three donors, both levels
    sl = ["C51ctr", "L01cov", "L04cov"]
    for lev in LEVELS:
        a, da = build(lev, sl)
        b, db = from_legacy(lev)
        print(lev, compare(a, b, sl), "| donors table equal:",
              da.loc[sl, ["group", "sex", "n_samples"]].astype(str).equals(db.loc[sl, ["group", "sex", "n_samples"]].astype(str)))
