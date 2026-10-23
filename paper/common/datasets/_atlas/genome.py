"""Genome-wide donor x compartment tables of the atlases used for cell-type gene rhythms (A27, A31b, A39, A50):
shared builder of the `<atlas>_genome` datasets. Same rows as the clock-gene tables (same groups, classes, majority assay,
minimum cells and donors), every gene instead of the 14 clock genes: `Dataset.gene_rhythms` input.

Raw inputs: `$SCCC_DATA/atlas_wg/<dataset>/pseudobulk.h5ad` (Census group pseudobulks with all genes, sparse; built by
desync/common/build_pseudobulk.py without a gene panel). Legacy: desync/validation/gonogo/scripts/sc_predict.py::load (heart,
muscle, breast, colon) and gonogo_new/scripts/atlases.py::table (lung_hlca, breast_hbca, liver_xu, skin_ganier, adipose).
Their covariates (`cov_age`, `cov_male`) are the donor columns `age` (years from the Census development stage; the legacy
used the stage's first number, 80+ = 80, missing -> median) and `male`. Legacy tables were cut to the GTEx-significant
genes; here genes are canonical symbols (BMAL1 -> ARNTL), first column of a duplicated symbol, genes without a count in the
table dropped.

Choices per dataset: see `SPECS`. All donors are `healthy` (these atlases are used as healthy references); lung_hlca keeps the
A4 healthy donors (Healthy, Healthy (tumor adjacent)) with stratum = assay | majority study; classes with >= 20 donors.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import census as C
import heart as H
import lung as LU
import organ as O
import organs as OR

NEW = {  # atlas folder: (organ prefix, region map, keep normal groups)
    "breast_hbca": ("breast", None), "liver_xu": ("liver", None), "skin_ganier": ("skin", None),
}
SPECS = {
    "heart_global_genome": "heart_global",
    "organ_breast_genome": "breast_kumar", "organ_muscle_genome": "muscle", "organ_colon_genome": "colon_crc",
    "organ_adipose_genome": "adipose_emont", "breast_hbca_genome": "breast_hbca", "liver_xu_genome": "liver_xu",
    "skin_ganier_genome": "skin_ganier", "lung_hlca_genome": "lung_hlca",
}
CLOCK_TWIN = {"heart_global_genome": "heart_global", "organ_breast_genome": "organ_breast",
              "organ_muscle_genome": "organ_muscle", "organ_colon_genome": "organ_colon",
              "organ_adipose_genome": "organ_adipose", "lung_hlca_genome": "lung_hlca"}


def _normal(o):
    return o.disease.astype(str).to_numpy() == "normal"


def spec(name: str) -> O.Spec:
    ds = SPECS[name]
    if name == "heart_global_genome":
        s = H.spec_global()
        s.path = C.ATLAS_WG / ds
        return s
    if name in ("organ_breast_genome", "organ_muscle_genome", "organ_colon_genome", "organ_adipose_genome"):
        organ = {"organ_breast_genome": "breast", "organ_muscle_genome": "muscle", "organ_colon_genome": "colon",
                 "organ_adipose_genome": "adipose"}[name]
        s = OR.spec(organ)
        if organ == "adipose":  # A31b: subcutaneous depot only
            s.region = {"subcutaneous adipose tissue": "subcutaneous"}
        s.path = C.ATLAS_WG / ds
        return s
    organ, region = NEW[ds]
    return O.Spec(C.ATLAS_WG / ds, ds, organ, O.COARSE, keep=_normal, region=region)


def build(name: str, donors=None, genes="all"):
    if genes != "all":
        raise ValueError("genome tables hold every gene (use the clock-gene dataset for the 14-gene panel)")
    ds = SPECS[name]
    path = C.ATLAS_WG / ds
    o, X, names = C.read_sparse(path, donors)
    dummy = pd.DataFrame(index=range(len(o)))
    full = donors is None
    if name == "lung_hlca_genome":
        ld = LU.hlca_donors()
        cls = o.ann_level_3.astype(str).map(LU.L3)
        cls[o.ann_level_2.astype(str) == "Smooth muscle"] = "smooth_muscle"
        res = "hlca_full:" + o.donor_id.astype(str)
        healthy = ld.index[~ld.lung_condition.str.contains(";") & ld.lung_condition.isin(LU.HEALTHY)]
        keep = (o.tissue_sampling_method.astype(str).isin(LU.HLCA_TISSUE) & (o.tissue.astype(str) != "nose")).to_numpy() \
            & cls.notna().to_numpy() & res.isin(healthy).to_numpy()
        F = C.keep_majority(C.frame(o, dummy, "hlca_full", "lung:" + cls.astype(str), o.assay, keep=keep))
        block = ld.study_major
        pb = C.sum_groups_sparse(F, X, names, F.resection.unique(), 20, "batch", 20 if full else None)
        pb["stratum"] = pb.stratum + "|" + pb.resection.map(block).astype(str)
        d = LU._hlca_cov()
        cols = ("age", "sex", "disease", "lung_condition", "study_major")
    else:
        s = spec(name)
        F = O.groups(s, donors, read=(o, dummy))
        pb = C.sum_groups_sparse(F, X, names, F.resection.unique(), 20, "batch", 20 if full else None)
        d = C.census_donors(s.path, s.prefix)
        cols = ("age", "sex", "disease")
    pb, out = O.finish(pb, d, {"healthy": O.flag(pb.resection.unique())}, cols=cols)
    out["male"] = (out.get("sex", pd.Series(index=out.index, dtype=object)).astype(str) == "male").astype(float)
    return pb, out


def from_legacy(name: str):
    raise NotImplementedError(
        "the legacy genome-wide tables are not stored (sc_predict.load rebuilds them); the clock-gene columns equal "
        f"those of `{CLOCK_TWIN.get(name, 'none')}` (see check_clock_twin)")


def check_clock_twin(name: str, n: int = 8, seed: int = 0) -> dict:
    """Clock-gene columns, n_cells, depth and strata of a slice of the genome table against the clock-gene dataset it
    shares its rows with (counts must be equal; rows of compartments dropped by that dataset's rules are ignored)."""
    import importlib.util
    from pathlib import Path

    twin = CLOCK_TWIN[name]
    f = Path(__file__).resolve().parents[1] / twin / "build.py"
    spec_ = importlib.util.spec_from_file_location(f"build_{twin}", f)
    mod = importlib.util.module_from_spec(spec_)
    spec_.loader.exec_module(mod)
    pb_t, _ = mod.from_legacy()
    rng = np.random.default_rng(seed)
    ids = sorted(rng.permutation(sorted(pb_t.donor.unique()))[:n])
    pb_g, _ = build(name, donors=[i.split(":", 1)[1] for i in ids])
    a = pb_t[pb_t.donor.isin(ids)].set_index(["donor", "compartment"]).sort_index()
    b = pb_g[pb_g.donor.isin(ids)].set_index(["donor", "compartment"]).sort_index()
    common = a.index.intersection(b.index)
    genes = list(C.GENES)
    diff = np.abs(a.loc[common, ["n_cells", "depth"] + genes].to_numpy(float)
                  - b.loc[common, ["n_cells", "depth"] + genes].to_numpy(float)).max()
    return dict(dataset=name, twin=twin, slice_donors=len(ids), rows_twin=len(a), rows_genome=len(b), common=len(common),
                only_twin=len(a.index.difference(b.index)), genes=pb_g.shape[1] - 5,
                max_abs_diff_clock=float(diff),
                strata_equal=bool((a.loc[common, "stratum"] == b.loc[common, "stratum"]).all()))
