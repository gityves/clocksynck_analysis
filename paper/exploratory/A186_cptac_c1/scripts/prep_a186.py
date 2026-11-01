"""A186 prep: per CPTAC-3 type, tumour / normal gene-count matrices of the paired cases (unstranded STAR counts, aliquots of one
sample type in one case summed, gene symbols with duplicates summed). Cache: output/cache/<TYPE>.npz; log: output/pairs.csv.
Run with `python -I` (the downloads are untrusted data; this script only parses them).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
TYPES = ["LUAD", "CCRCC", "LSCC", "SKCM", "HNSCC", "PDAC", "UCEC", "STAD"]


def read_counts(path: Path) -> pd.Series:
    t = pd.read_csv(path, sep="\t", comment="#", usecols=["gene_id", "gene_name", "unstranded"])
    t = t[t.gene_id.str.startswith("ENSG")]
    return t.groupby("gene_name").unstranded.sum()


def build(tp: str):
    d = DATA / f"cptac_{tp}"
    s = pd.read_csv(d / f"{tp}_samples.tsv", sep="\t")
    files = {p.name.split(".")[0]: p for p in d.glob("*.rna_seq.augmented_star_gene_counts.tsv")}
    miss = set(s.file_id) - set(files)
    assert not miss, f"{tp}: {len(miss)} files missing"
    ok = s.groupby("case").sample_type.agg(lambda x: {"Primary Tumor", "Solid Tissue Normal"} <= set(x))
    cases = sorted(ok.index[ok])
    out, info = {}, {}
    for kind in ("Primary Tumor", "Solid Tissue Normal"):
        rows = []
        for c in cases:
            f = s[(s.case == c) & (s.sample_type == kind)].file_id
            x = sum(read_counts(files[i]) for i in f)
            rows.append(x)
            info[(c, kind)] = len(f)
        out[kind] = pd.concat(rows, axis=1, keys=cases).T.fillna(0)
    genes = sorted(set(out["Primary Tumor"].columns) & set(out["Solid Tissue Normal"].columns))
    return cases, out["Primary Tumor"][genes], out["Solid Tissue Normal"][genes], info


def main():
    (HERE / "output" / "cache").mkdir(exist_ok=True, parents=True)
    rows = []
    for tp in TYPES:
        cases, ct, cn, info = build(tp)
        np.savez(HERE / "output" / "cache" / f"{tp}.npz", cases=np.array(cases), genes=np.array(ct.columns, dtype=str),
                 ct=ct.to_numpy(np.float32), cn=cn.to_numpy(np.float32))
        rows.append(dict(type=tp, n_pairs=len(cases), n_genes=ct.shape[1],
                         cases_multi_tumour_aliquots=sum(v > 1 for (c, k), v in info.items() if k == "Primary Tumor"),
                         cases_multi_normal_aliquots=sum(v > 1 for (c, k), v in info.items() if k != "Primary Tumor"),
                         median_depth_tumour=float(ct.sum(axis=1).median()), median_depth_normal=float(cn.sum(axis=1).median())))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(HERE / "output" / "pairs.csv", index=False)


if __name__ == "__main__":
    main()
