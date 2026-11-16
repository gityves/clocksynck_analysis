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
import cortex as CX
import heart as H
import kong as K
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
    "organ_optic_nerve_genome": "optic_nerve", "organ_eye_anterior_genome": "eye_anterior",  # A17 / A27 extension
    "psychad_hbcc_genome": "psychad_hbcc", "seaad_dlpfc_genome": "seaad_dlpfc", "seaad_mtg_genome": "seaad_mtg",
    "psychad_mssm_genome": "psychad_mssm", "psychad_radc_genome": "psychad_radc",  # A280
    "kidney_sc_genome": "kpmp_sc", "islets_genome": "islets_t2d",
    "kong_colon_genome": "kong_co", "kong_ileum_genome": "kong_ti",  # three folders each (`kong.SEGMENTS`)
}
KONG = {"kong_colon_genome": "colon", "kong_ileum_genome": "ileum"}
CORTEX = ("psychad_hbcc_genome", "seaad_dlpfc_genome", "seaad_mtg_genome", "psychad_mssm_genome", "psychad_radc_genome")  # rules of `cortex.build`
DISEASE = {  # name: (label column rule, healthy min donors, block column or None, donor columns, extra donor column)
    "kidney_sc_genome": (["normal"], ["chronic kidney disease"], 20, "dataSource", "ckd_sc"),
    "islets_genome": (["normal"], ["type 2 diabetes mellitus"], 15, None, "islets_t2d"),
}
CLOCK_TWIN = {"heart_global_genome": "heart_global", "organ_breast_genome": "organ_breast",
              "organ_muscle_genome": "organ_muscle", "organ_colon_genome": "organ_colon",
              "organ_adipose_genome": "organ_adipose", "lung_hlca_genome": "lung_hlca",
              "organ_optic_nerve_genome": "organ_optic_nerve", "organ_eye_anterior_genome": "organ_eye_anterior",
              "psychad_hbcc_genome": "psychad_hbcc", "seaad_dlpfc_genome": "seaad_dlpfc",
              "seaad_mtg_genome": "seaad_mtg", "psychad_mssm_genome": "psychad_mssm",
              "psychad_radc_genome": "psychad_radc", "kidney_sc_genome": "kidney_sc", "islets_genome": "islets",
              "kong_colon_genome": "kong_colon", "kong_ileum_genome": "kong_ileum"}


ORGAN = {"organ_breast_genome": "breast", "organ_muscle_genome": "muscle", "organ_colon_genome": "colon",
         "organ_adipose_genome": "adipose", "organ_optic_nerve_genome": "optic_nerve",
         "organ_eye_anterior_genome": "eye_anterior"}

LEGACY_TWIN = tuple(CLOCK_TWIN)[:6]  # checked against the legacy A17 pickles (first version of the check)


def _normal(o):
    return o.disease.astype(str).to_numpy() == "normal"


def spec(name: str) -> O.Spec:
    ds = SPECS[name]
    if name == "heart_global_genome":
        s = H.spec_global()
        s.path = C.ATLAS_WG / ds
        return s
    if name in ORGAN:
        organ = ORGAN[name]
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
    if name in KONG:  # as kong.build: three files under the donor key kong:, majority assay per file, >= 15 donors per class
        import scipy.sparse as sp

        seg, Fs, Xs, off = KONG[name], [], [], 0
        for s in K.specs(seg):
            s.path = C.ATLAS_WG / s.path.name
            o, X, names = C.read_sparse(s.path, donors)
            F = O.groups(s, donors, read=(o, pd.DataFrame(index=range(len(o)))))
            F["_row"] += off
            Fs.append(F), Xs.append(X)
            off += X.shape[0]
            if len(Xs) > 1 and not np.array_equal(names, names0):
                raise ValueError(f"{s.path.name}: genes differ from the first file")
            names0 = names
        pb = C.sum_groups_sparse(pd.concat(Fs, ignore_index=True), sp.vstack(Xs).tocsr(), names0,
                                 pd.concat(Fs).resection.unique(), 20, "batch", 15 if donors is None else None)
        pb, out = O.finish(pb, K._donors(seg), {"healthy": O.flag(pb.resection.unique())})
        out["male"] = (out.sex.astype(str) == "male").astype(float)
        return pb, out
    o, X, names = C.read_sparse(path, donors)
    dummy = pd.DataFrame(index=range(len(o)))
    full = donors is None
    if name in CORTEX:  # as cortex.build: 10x groups, BRAIN classes (case-sensitive), >= 50 nuclei, stratum <ds>|<assay>
        cell = C.classify(o.cell_type.astype(str), CX.BRAIN, ignore_case=False)
        umi = o.assay.astype(str).str.startswith("10x").to_numpy()
        F = C.keep_majority(C.frame(o, dummy, ds, o.tissue_general.astype(str) + ":" + cell.astype(str), o.assay,
                                    keep=umi & cell.notna().to_numpy()))
        pb = C.sum_groups_sparse(F, X, names, None, CX.MIN_CELLS, "batch", None)
        pb["stratum"] = ds + "|" + pb.stratum
        d = CX._donors(ds)
        extra = {c: CX.labels(c, d).rename(c) for c, v in CX.CONTRASTS.items() if v[0] == ds}
        if full and ds in CX.HEALTHY_FIT:
            extra["healthy"] = pd.Series(True, index=CX.healthy_donors(pb, d))
        pb = C.standard(pb)
        out = C.donors_frame(d, pb.donor.unique(), extra)
        if "healthy" in out:
            out["healthy"] = out.healthy.fillna(False).astype(bool)
        out["male"] = (out.get("sex", pd.Series(index=out.index, dtype=object)).astype(str) == "male").astype(float)
        return pb, out
    if name in DISEASE:  # as kidney_sc / islets builders: control + case donors, classes of their CLASSES
        ctl, case, min_h, block_col, tag = DISEASE[name]
        mod = _clock_module(CLOCK_TWIN[name])
        s = O.Spec(path, ds, mod.SPEC.organ, mod.SPEC.classes)
        d = C.census_donors(path, ds)
        lab = O.labels_from(d, "disease", ctl, case)
        if donors is not None:
            lab = lab[lab.index.isin([f"{ds}:{x}" for x in donors])]
        F = O.groups(s, donors, read=(o, dummy))
        pb = C.sum_groups_sparse(F, X, names, lab.index, 20, "batch", None)
        if block_col:
            pb["stratum"] = pb.stratum + "|" + pb.resection.map(d[block_col].astype(str)).astype(str)
        extra = {tag: lab.rename(tag)}
        if full:
            extra["healthy"] = O.flag(O.healthy_members(pb, lab.index[lab == "control"], min_h))
        pb, out = O.finish(pb, d, extra, cols=("age", "sex", "disease") + ((block_col,) if block_col else ()))
        out["male"] = (out.sex.astype(str) == "male").astype(float)
        return pb, out
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


def _clock_module(twin: str):
    import importlib.util
    from pathlib import Path

    f = Path(__file__).resolve().parents[1] / twin / "build.py"
    spec_ = importlib.util.spec_from_file_location(f"build_{twin}", f)
    mod = importlib.util.module_from_spec(spec_)
    spec_.loader.exec_module(mod)
    return mod


def from_legacy(name: str):
    raise NotImplementedError(
        "the legacy genome-wide tables are not stored (sc_predict.load rebuilds them); the clock-gene columns equal "
        f"those of `{CLOCK_TWIN.get(name, 'none')}` (see check_clock_twin)")


def check_clock_twin(name: str, n: int = 8, seed: int = 0) -> dict:
    """Clock-gene columns, n_cells, depth and strata of a slice of the genome table against the clock-gene dataset it
    shares its rows with (counts must be equal; rows of compartments dropped by that dataset's rules are ignored)."""
    twin = CLOCK_TWIN[name]
    mod = _clock_module(twin)  # rebuilt clock-gene table (the legacy pickles for the first atlases)
    pb_t, _ = mod.from_legacy() if name in LEGACY_TWIN else mod.build()
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
