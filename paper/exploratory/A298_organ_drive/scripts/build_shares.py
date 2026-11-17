"""A298 data plumbing: share of clock-gene mRNA contributed by each cell family in GTEx tissues (composition only).

Method (no clock statistic is computed or read):
1. Each atlas compartment is assigned a cell family with `clocksync.anchors.family`.
2. Family profile = counts summed over all healthy pseudobulk rows of the family / summed depth x 1e6 (CPM).
   Families with < 2,000 summed cells (column `n_cells`) are dropped.
3. Markers per family and atlas: among genes with family CPM >= 10 and present in GTEx, the 50 with the largest
   log2((CPM_f + 1) / (max CPM over the other families + 1)); the 14 clock genes are excluded.
4. GTEx profile per tissue: fast deaths (Hardy 1-2), mean over samples of counts / depth x 1e6, on markers + clock genes.
5. NNLS of the tissue profile on the family profiles over the union of the atlas markers (raw CPM, each gene weighted by
   1 / (mean family CPM + 1)); coefficients normalised to sum 1 = mRNA fraction per family.
6. Clock-gene mRNA share of a family = fraction x sum of the family's CPM over the 14 clock genes, normalised to sum 1.
GTEx is controlled access: only per-tissue aggregates are written, never per-donor rows.
Outputs: output/shares.csv, output/design.csv.
"""
from __future__ import annotations

import gc
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_2/scripts/gene_rhythms"))
import clocksync as cs  # noqa: E402
from clocksync.anchors import family  # noqa: E402
from common import gtex  # noqa: E402
import _shared as S  # noqa: E402

OUT = HERE.parent / "output"
CLOCK = list(cs.gtex().genes)
MIN_CELLS, MIN_CPM, N_MARK, MIN_DONORS = 2000, 10.0, 50, 40
MAP = {
    "heart_global_genome": ["Heart - Left Ventricle", "Heart - Atrial Appendage"],
    "lung_hlca_genome": ["Lung"],
    "organ_breast_genome": ["Breast - Mammary Tissue"],
    "organ_muscle_genome": ["Muscle - Skeletal"],
    "organ_adipose_genome": ["Adipose - Subcutaneous"],
    "organ_colon_genome": ["Colon - Transverse"],
    "skin_ganier_genome": ["Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"],
    "kidney_sc_genome": ["Kidney - Cortex"],
    "liver_xu_genome": ["Liver"],
}


def family_profiles(pb: pd.DataFrame, genes: list) -> tuple[pd.DataFrame, list, pd.Series]:
    """CPM profile (genes x families) of the kept families, dropped families, summed cells per family."""
    fam = pb.compartment.map(family)
    cells = pb.groupby(fam).n_cells.sum()
    dropped = sorted(cells.index[cells < MIN_CELLS])
    kept = [f for f in sorted(cells.index) if f not in dropped]
    cols = {}
    for f in kept:
        m = (fam == f).to_numpy()
        cols[f] = pb.loc[m, genes].to_numpy(float).sum(0) / pb.loc[m, "depth"].sum() * 1e6
    return pd.DataFrame(cols, index=genes), dropped, cells


def markers(prof: pd.DataFrame, allowed: set) -> dict:
    out = {}
    for f in prof.columns:
        other = prof.drop(columns=f).max(axis=1) if prof.shape[1] > 1 else pd.Series(0.0, index=prof.index)
        lfc = np.log2((prof[f] + 1) / (other + 1))
        ok = (prof[f] >= MIN_CPM) & prof.index.isin(allowed) & ~prof.index.isin(CLOCK)
        out[f] = list(lfc[ok].sort_values(ascending=False).index[:N_MARK])
    return out


def gtex_vars() -> set:
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        return set(map(str, a.var_names))
    finally:
        a.file.close()


def run_atlas(name: str, tissues: list, allowed: set):
    pb, donors = S.genome_table(name)
    donors_used = len(donors)
    if "healthy" in donors.columns:
        h = set(donors.index[donors.healthy.astype(bool)])
        pb = pb[pb.donor.isin(h)].reset_index(drop=True)
        donors_used = len(h & set(pb.donor))
    clock_in = [g for g in CLOCK if g in pb.columns]
    genes = S.gene_columns(pb)
    prof, dropped, cells = family_profiles(pb, genes)
    del pb
    gc.collect()
    mk = markers(prof, allowed)
    union = sorted({g for v in mk.values() for g in v})
    clock_cpm = prof.loc[clock_in].sum()
    A = prof.loc[union]
    g, _ = gtex.load(tissues=tissues, genes=union + CLOCK, min_donors=MIN_DONORS)
    rows, design = [], []
    for t in tissues:
        sub = g[g.compartment == t]
        if sub.empty:
            design.append(dict(tissue=t, atlas=name, note="tissue absent or < 40 donors"))
            continue
        cpm = sub[union].to_numpy(float) / sub[["depth"]].to_numpy(float) * 1e6
        y = cpm.mean(0)
        w = 1.0 / (A.mean(axis=1).to_numpy() + 1.0)
        coef, _ = nnls(A.to_numpy() * w[:, None], y * w)
        frac = coef / coef.sum()
        fit = A.to_numpy() @ coef
        r2 = 1 - ((y - fit) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        r2w = 1 - (((y - fit) * w) ** 2).sum() / (((y * w) - (y * w).mean()) ** 2).sum()
        cs_ = frac * clock_cpm.to_numpy()
        cs_ = cs_ / cs_.sum()
        for i, f in enumerate(prof.columns):
            rows.append(dict(tissue=t, atlas=name, family=f, mrna_fraction=frac[i], clock_share=cs_[i],
                             n_markers=len(union), fit_r2=r2, fit_r2_weighted=r2w))
        design.append(dict(tissue=t, atlas=name, families_kept=";".join(prof.columns), families_dropped=";".join(dropped),
                           donors_atlas=donors_used, gtex_samples=len(sub), note=""))
    return rows, design


def main():
    t0 = time.time()
    OUT.mkdir(exist_ok=True)
    allowed = gtex_vars()
    rows, design, failed = [], [], []
    for name, tissues in MAP.items():
        try:
            r, d = run_atlas(name, tissues, allowed)
            rows += r
            design += d
        except Exception as e:  # record and continue
            failed.append((name, repr(e)))
            design += [dict(tissue=t, atlas=name, note=f"FAILED: {e!r}") for t in tissues]
        gc.collect()
        print(name, "done", round(time.time() - t0), "s", flush=True)
    sh, de = pd.DataFrame(rows), pd.DataFrame(design)
    sh.to_csv(OUT / "shares.csv", index=False)
    de.to_csv(OUT / "design.csv", index=False)
    pd.set_option("display.width", 250, "display.max_rows", 500, "display.max_colwidth", 80)
    print(sh.round(4).to_string(index=False))
    print(de.to_string(index=False))
    print("failed:", failed)
    print("peak GB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, "runtime s", round(time.time() - t0))


if __name__ == "__main__":
    main()
