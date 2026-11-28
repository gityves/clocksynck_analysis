"""Shared helpers of the A426 builds: donor x compartment genome-wide count pseudobulks (A401 format), the A149 rules and
the provenance files. COUNTS ONLY: no clock-gene statistic, loudness or immune score is computed in any build.

Output of a source (`$SCCC_DATA/<dir>/`): `pseudobulk.parquet` (one row per donor x sample x compartment with >= 20 cells:
donor = `<source>|<patient>`, sample, compartment = `tumour:<class>`, stratum = `<source>|<chemistry>`, n_cells, depth,
median_umi of the group, units, sample_median_umi, then one count column per gene symbol), `samples.csv`, `donors.csv`,
`PROVENANCE.md` (URL, date, sha256 of each downloaded file).
"""
from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "build_geo"))
import geobuild as G  # noqa: E402  (pseudobulk, read_h5_10x, marker_comp)

DATA = G.OUT.parent
MIN_UMI = 1500  # A149: sample median UMI
CLASSES = ("fibroblast", "endothelial", "mural", "epithelial", "malignant", "T_NK", "B_plasma", "myeloid", "mast")


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def sample_pseudobulk(X, sym, umi, comp, donor, sample, stratum, units="counts"):
    """Pseudobulk of one sample (cells x genes counts `X`; `comp` = class per cell, `unassigned` / None dropped)."""
    comp = np.asarray(comp, dtype=object)
    ok = ~pd.isna(comp) & (comp != "unassigned") & (comp != "drop")
    lab = np.where(ok, "tumour:" + comp.astype(str), "drop")
    obs = pd.DataFrame({"donor": donor, "sample": sample, "compartment": lab, "stratum": stratum})
    pb = G.pseudobulk(X, sym, obs, units)
    pb["sample_median_umi"] = float(np.median(umi))
    return pb


def concat_frames(parts):
    """Row-wise union of pseudobulk frames with different gene sets (missing genes 0)."""
    meta = ["donor", "sample", "compartment", "stratum", "n_cells", "depth", "median_umi", "units", "sample_median_umi"]
    genes = pd.Index(pd.unique(np.concatenate([[c for c in p.columns if c not in meta] for p in parts])))
    out = []
    for p in parts:
        g = p.drop(columns=meta).reindex(columns=genes, fill_value=0).astype(np.float32)
        out.append(pd.concat([p[meta].reset_index(drop=True), g.reset_index(drop=True)], axis=1))
    return pd.concat(out, ignore_index=True)


def write_source(directory, pb, samples, donors, provenance: str):
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    pb.to_parquet(d / "pseudobulk.parquet", index=False)
    samples.to_csv(d / "samples.csv", index=False)
    donors.to_csv(d / "donors.csv", index=False)
    (d / "PROVENANCE.md").write_text(provenance)
    log(f"wrote {len(pb)} pseudobulks to {d}")


def donor_counts(pb: pd.DataFrame) -> pd.DataFrame:
    """Per donor: cells per class (pooled over samples; classes with >= 50 cells) and the A426 unit counts."""
    t = pb[pb.compartment.str.startswith("tumour:")]
    t = t.groupby(["donor", "compartment"], as_index=False).n_cells.sum().query("n_cells >= 50")  # pooled over samples, >= 50 cells
    n = t.pivot_table(index="donor", columns=t.compartment.str.split(":").str[1], values="n_cells", aggfunc="sum").fillna(0)
    for c in CLASSES:
        if c not in n:
            n[c] = 0
    n["n_stroma_classes"] = (n[["fibroblast", "endothelial", "mural"]] >= 50).sum(1)
    n["n_immune_cells"] = n[["T_NK", "B_plasma", "myeloid", "mast"]].sum(1)
    return n.reset_index()
