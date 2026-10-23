"""Delorey et al. 2021 (Nature): snRNA of lung, heart, liver and kidney from 24 COVID-19 autopsy donors.

Source:   GEO GSE171668 (per-sample CellRanger `raw_feature_bc_matrix.h5.gz`, per-organ author metadata csv.gz),
          Nature Extended Data Table 1 (`suppl/41586_2021_3570_MOESM3_ESM.xlsx`, sheet `donor_metadata`),
          all in `$SCCC_DATA/covid/GSE171668/`.
Replaces: legacy/desync/validation/natural_experiments/scripts/build_delorey.py (+ pbcommon.py), build_covid.md.
Used by:  A105 (organs, lung classes), A105b, A105c (post hoc).

Choices (as the legacy builder)
- cells: barcodes of the author metadata (post-QC nuclei) of samples with `method == nuclei` (cryo and fresh
  preparations dropped), `doublet == TRUE` dropped; counts are the raw CellRanger UMI (all-droplet matrix, `GRCh38premrna_`
  prefix stripped from gene names);
- donor: leading `D<number>` of the author donor field; all nuclei samples of a donor x organ are summed;
- level "organ" (default): compartment = organ (lung, heart, liver, kidney); level "class": `organ:class`, the author
  annotation mapped to broad classes by `CLASSMAP` (doublet, technical and mixed labels dropped; fixed before any clock
  quantity was computed);
- stratum: the assay (10x 3' nuclei, one value); every pseudobulk with at least one nucleus is kept.

Donor table: Extended Data Table 1 (age range, sex, IMV days, symptoms-to-death days, PMI for sample processing in
hours, organ failure, `melms_id` for the lung of donors D25-D32 that are also in Melms) and `has_<organ>`, `n_organs`.
"""

from __future__ import annotations

import gzip
import io
import re
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(Path(__file__).resolve().parents[3])]
import clocksync as cs  # noqa: E402
from common.paths import DATA  # noqa: E402
from natexp_common import clock_adata, compare, read_legacy, standard  # noqa: E402

ROOT = DATA / "covid" / "GSE171668"
LEVELS = ("organ", "class")
ASSAY = "10x 3' nuclei"

ORGANS = ["lung", "heart", "liver", "kidney"]
AGE = {0: "0-2", 0.2: "2-5", 0.5: "5-10", 0.1: "10-17", 1.8: "18-29", 3: "30-39", 4: "40-49", 5: "50-59",
       6: "60-69", 7: "70-79", 8: ">=80"}

# author label (regex, first match wins) -> broad class; None = dropped. Applied to `manual_annotation`
# (lung: `manual_annotation_cluster`).
CLASSMAP = {
    "lung": [("^Myeloid$", "myeloid"), ("^Epithelial$", "epithelial"), ("^Endothelial$", "endothelial"),
             ("^Fibroblasts$", "fibroblast"), ("^T\\+NK$", "T_NK"), ("^Vascular Contractile$", "vascular_contractile"),
             ("^Ciliated$", "ciliated"), ("^B\\+Plasma$", "B_plasma"), ("^Secretory$", "secretory"),
             ("^MAST$", "mast"), ("^RBC$", "RBC"), ("^Mesothelial$", "mesothelial")],
    "heart": [(r"^\d+: (CM \+|Fib \+|EC \+)", None), (r"^\d+: (Mito )?CM( \d|$| \()", "cardiomyocyte"),
              (r"^\d+: (Low count )?Fib( \d|$)", "fibroblast"), (r"^\d+: (Mito )?EC( \d| \(|$)", "endothelial"),
              (r"^\d+: Pericyte", "pericyte"), (r"^\d+: Macrophage", "macrophage"), (r"^\d+: T cell", "T"),
              (r"^\d+: VSMC", "vSMC"), (r"^\d+: Neuronal", "neuronal"), (r"^\d+: Lymph EC", "lymphatic_EC"),
              (r"^\d+: Adipocyte", "adipocyte")],
    "kidney": [(r"Technical|Doublet|\(Mes\)|\(EC\)|\(Imm\)|Mito-rich (Distal|Int|\(Imm\))", None),
               (r"^(PCT|Mito-rich PCT|Stressed_PCT)", "PT"), ("^TAL$", "TAL"), ("^tDL$", "tDL"), ("^DCT$", "DCT"),
               ("^CNT$", "CNT"), ("^CD-PC$", "CD_PC"), ("^CD-IC", "CD_IC"), ("^Podocyte$", "podocyte"),
               ("^PEC$", "PEC"), (r"EC$", "endothelial"), ("^Fibroblast", "fibroblast"),
               (r"Pericyte", "pericyte"), ("^(Macrophages|Myeloid)$", "myeloid"), ("^(T-cells|Proliferating T)", "T"),
               ("^NK$", "NK"), ("^Bcells$", "B"), ("^Plasma$", "plasma")],
    "liver": [(r"Technial|Doublet|\((Hep|Imm|EC|Mes)\)|Mito-rich (doublet|\(Mes\))|T-cell/Doublet", None),
              (r"(?i)hepatocyte", "hepatocyte"), ("^Cholangiocytes$", "cholangiocyte"),
              (r"^(LSEC|VWF\+ANO2|Mito-rich EC|Lymphatic EC)", "endothelial"),
              (r"^(Kuppfer cells|CD163\+Macrophage|SPP1\+Macrophage|Macrophage-)", "macrophage"),
              ("^(T/NKT|MAIT)$", "T_NKT"), ("^NK$", "NK"), (r"^(RSPO3\+ fibro|Fibroblast)$", "fibroblast"),
              ("^Stellate", "stellate"), ("^Monocyte$", "monocyte"), ("^Erythroid", "erythroid"),
              ("^Plasma cells$", "plasma"), ("^Cycling$", "cycling"), ("^Myocyte$", None)],
}


def classify(organ, label):
    for pat, cls in CLASSMAP[organ]:
        if re.search(pat, label):
            return cls
    return None


def clinical(xlsx):
    d = pd.read_excel(xlsx, sheet_name="donor_metadata", header=None).dropna(how="all")
    head = d.index[d[0] == "Donor"][0]
    t = d.loc[head + 1:].copy()
    t.columns = [str(c).strip() for c in d.loc[head]]
    t = t[t["Donor"].astype(str).str.match(r"^D\d+$")]

    def hours(s):
        m = re.fullmatch(r"\s*(?:(\d+)\s*h)?\s*(?:(\d+)\s*min)?\s*", str(s))
        return np.nan if m is None else float(m.group(1) or 0) + float(m.group(2) or 0) / 60

    out = pd.DataFrame(dict(
        donor=t.Donor.to_numpy(), broad_id=t["Broad ID"].to_numpy(),
        age_range=[AGE.get(float(v), str(v)) for v in t["Age (range)"]],
        sex=np.where(t.Sex.astype(float) == 1, "female", "male"),
        imv_days=pd.to_numeric(t.IMV, errors="coerce").to_numpy(),
        symptoms_to_death_days=pd.to_numeric(t["Number of days: Symptoms to death"], errors="coerce").to_numpy(),
        pmi_h=[hours(v) for v in t["Post-Mortem Interval (PMI) for sample processing"]],
        covid=t["COVID-19 status"].to_numpy(), analysis=t["Analysis type"].to_numpy(),
        organ_failure=t["Organ failure at time of death"].to_numpy(),
        biobanked=t["Biobanked tissue"].to_numpy(), melms_id=t["Melms et al. ID for lung"].to_numpy()))
    return out.set_index("donor")




def read_h5(gz: Path):
    """10x raw_feature_bc_matrix.h5.gz -> (csc genes x barcodes, barcodes, gene symbols); decompressed in memory."""
    with gzip.open(gz, "rb") as f:
        buf = io.BytesIO(f.read())
    with h5py.File(buf, "r") as h:
        m = h["matrix"]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:]))
        bc = np.array([b.decode().removesuffix("-1") for b in m["barcodes"][:]])
        genes = [g.decode().removeprefix("GRCh38premrna_") for g in m["features/name"][:]]
    return X, bc, genes


def cell_metadata() -> pd.DataFrame:
    """Author metadata of the four organs: donor, organ, sample, class, with the nuclei / non-doublet selection."""
    frames = []
    for o in ORGANS:
        d = pd.read_csv(ROOT / "raw" / f"GSE171668_{o}_metadata.csv.gz", encoding="utf-8-sig")
        d["organ"] = o
        d["donor"] = d.donor.str.extract(r"^(D\d+)")[0]
        d["sample"] = d.barcodes.str.rsplit("-", n=1).str[0]
        label = d["manual_annotation_cluster" if o == "lung" else "manual_annotation"].astype(str)
        d["cls"] = [classify(o, s) for s in label]
        frames.append(d)
    m = pd.concat(frames, ignore_index=True)
    m = m[(m.method == "nuclei") & (~m.doublet.astype(bool))].copy()
    m["organ_key"] = m.organ
    m["class_key"] = np.where(m.cls.notna(), m.organ + ":" + m.cls.astype(str), None)
    return m


def read_sample(path: Path, g: pd.DataFrame, sample: str) -> ad.AnnData:
    X, bc, genes = read_h5(path)
    col = pd.Series(np.arange(len(bc)), index=bc)
    idx = col.reindex(g.barcodes.str[len(sample) + 1:]).to_numpy()
    ok = ~np.isnan(idx)
    idx = idx[ok].astype(int)
    Xs = X[:, idx].tocsr()  # genes x selected nuclei
    obs = g.loc[ok, ["donor", "organ_key", "class_key"]].copy()
    obs.index = g.barcodes[ok].to_numpy()
    return clock_adata(Xs.T, genes, np.asarray(Xs.sum(0)).ravel(), obs)


def donor_table(m: pd.DataFrame) -> pd.DataFrame:
    clin = clinical(next((ROOT / "suppl").glob("*MOESM3_ESM.xlsx")))
    don = clin.copy()
    have = m.groupby(["donor", "organ"]).size()
    for o in ORGANS:
        don[f"has_{o}"] = don.index.isin(have.xs(o, level="organ").index)
    don["n_organs"] = don[[f"has_{o}" for o in ORGANS]].sum(axis=1)
    return don.rename_axis("donor")


def build(level: str = "organ", donors: list | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pseudobulk, donors). `donors`: restrict to these donors (checks; only their samples are read)."""
    m = cell_metadata()
    don = donor_table(m)
    if donors is not None:
        m = m[m.donor.isin(donors)]
    cells = []
    for f in sorted((ROOT / "ext").glob("GSM*_raw_feature_bc_matrix.h5.gz")):
        sample = re.match(r"GSM\d+_(.+)_raw_feature_bc_matrix", f.name).group(1)
        g = m[m["sample"] == sample]
        if len(g):
            cells.append(read_sample(f, g, sample))
    pb = cs.pseudobulk(cells, donor="donor", compartment="organ_key" if level == "organ" else "class_key",
                       depth="depth", min_cells=1)
    pb["stratum"] = ASSAY
    return standard(pb), don if donors is None else don.loc[don.index.isin(donors)]


def from_legacy(level: str = "organ") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Standard table from the legacy `pseudobulk_<level>/pseudobulk.h5ad` and `samples.csv` (class level:
    compartments `organ:class`)."""
    d = ROOT / f"pseudobulk_{level}"
    pb = read_legacy(d / "pseudobulk.h5ad")
    pb["stratum"] = ASSAY
    don = pd.read_csv(d / "samples.csv").set_index("donor_id").rename_axis("donor")
    return standard(pb), don


if __name__ == "__main__":
    sl = ["D1", "D3"]  # D3: two samples of heart and liver summed
    for lev in LEVELS:
        a, da = build(lev, sl)
        b, db = from_legacy(lev)
        print(lev, compare(a, b, sl), "| donors table equal:",
              da.loc[sl, ["sex", "imv_days", "n_organs"]].astype(str).equals(db.loc[sl, ["sex", "imv_days", "n_organs"]].astype(str)))
