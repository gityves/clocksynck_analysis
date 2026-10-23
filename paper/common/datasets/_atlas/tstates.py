"""Lung T-cell states (A6, A7): the T class split by per-cell states, in the HLCA (healthy and fibrotic donors) and in Natri et al. 2024
(controls and cases): shared builder of `hlca_tstates` and `natri_tstates`.

Legacy: desync/08_lung_t_residency/scripts/residency.py (`t_compartment`, `groups`, `load`, `donor_sets`); per-cell states from
desync/common/cell_states.py (score = mean log1p(1e4 count / raw_sum) of RESIDENT minus CIRCULATING markers, tertile within block, cycling =
>= 3 of 12 genes incl. MKI67 / TOP2A; folded into the group pseudobulks as obs columns `t_state`, `t_state_core`, `t_cycling`, `t_phase`).
Raw inputs (group pseudobulks): $SCCC_DATA/atlas/hlca_tcell and natri2024_tcell; labels `inputs/lung_donors.csv`.

`build(name, mode=...)`: the T class is replaced by the compartments of one mode (the legacy `T_MODES`): "T" (one class `lung:T`), "state"
(`lung:T_resident_like`, `lung:T_circulating_like`; intermediate dropped), "state_core", "cd8_state" (`lung:T_CD8_resident_like`,
`..._circulating_like`), "subset" (`lung:T_CD4`, `lung:T_CD8`), "noncycling" (`lung:T_noncycling`), "g1" (`lung:T_G1`, Natri only). The other
classes are those of `lung.py` (HLCA: A4 classes, lung tissue, majority assay, stratum = study | assay; Natri: A5 classes, all cells,
stratum = majority processing site). Pseudobulks >= `t_min` cells for T compartments (20; 50 as sensitivity) and >= 50 (HLCA) / 20
(Natri) for the other classes; no minimum of donors per class (A6 used 20 in healthy sets: apply it per donor set). Donor column
`tstate_set`: HLCA healthy / fibrotic donors (single lung condition; `HEALTHY`, `FIBROTIC`), Natri control / case (primary donors).
No legacy table was stored: checked against twins (`check_twins`): mode "T" equals `hlca_fibrotic`, mode "noncycling" in Natri controls
equals `lung_natri_control`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import census as C
import lung as LU
import organ as O

PARAM = {"hlca_tstates": dict(dir="hlca_tcell", ref="lung:macrophage", min_cells=50, subset="ann_level_4", cd4="CD4 T cells",
                              cd8="CD8 T cells"),
         "natri_tstates": dict(dir="natri2024_tcell", ref="lung:monocyte", min_cells=20, subset="manual_annotation_1", cd4="CD4",
                               cd8="CD8/NKT")}
MODES = ("T", "state", "state_core", "cd8_state", "subset", "noncycling", "g1")


def t_compartment(o: pd.DataFrame, p: dict, mode: str) -> pd.Series:
    sub = o[p["subset"]].astype(str).map({p["cd4"]: "CD4", p["cd8"]: "CD8"})
    state = o.t_state.astype(str).where(o.t_state.astype(str) != "intermediate")
    if mode == "T":
        return pd.Series("T", index=o.index)
    if mode == "state":
        return "T_" + state
    if mode == "state_core":
        s = o.t_state_core.astype(str)
        return "T_" + s.where(s != "intermediate")
    if mode == "cd8_state":
        return ("T_CD8_" + state).where(sub == "CD8")
    if mode == "subset":
        return "T_" + sub
    if mode == "noncycling":
        return pd.Series("T_noncycling", index=o.index).where(o.t_cycling.astype(str) == "noncycling")
    if mode == "g1":
        if "t_phase" not in o:
            raise ValueError("mode 'g1' needs the author cell-cycle phase (Natri only)")
        return pd.Series("T_G1", index=o.index).where(o.t_phase.astype(str) == "G1")
    raise ValueError(mode)


def build(name: str, donors=None, genes="clock", mode: str = "state", t_min: int = 20):
    p = PARAM[name]
    path = C.ATLAS / p["dir"]
    o, X = C.read(path, genes, donors)
    if name == "hlca_tstates":
        cls = o.ann_level_3.astype(str).map(LU.L3)
        cls[o.ann_level_2.astype(str) == "Smooth muscle"] = "smooth_muscle"
        keep = (o.tissue_sampling_method.astype(str).isin(LU.HLCA_TISSUE) & (o.tissue.astype(str) != "nose")).to_numpy()
        batch, prefix = o.assay, "hlca_full"
    else:
        cls = LU._natri_cls(o)
        keep, batch, prefix = np.ones(len(o), bool), o.Processing_site, "natri2024_ild"
    is_t = cls == "T"
    cls = cls.where(~is_t, t_compartment(o, p, mode))
    F = C.frame(o, X, prefix, "lung:" + cls.astype(str), batch, keep=keep & cls.notna().to_numpy())
    if name == "hlca_tstates":
        ld = LU.hlca_donors()
        ld = ld[~ld.lung_condition.str.contains(";")]
        F = C.keep_majority(F)
        sets = pd.concat([pd.Series("healthy", index=ld.index[ld.lung_condition.isin(LU.HEALTHY)]),
                          pd.Series("fibrotic", index=ld.index[ld.lung_condition.isin(LU.FIBROTIC)])])
        stratum = F.resection.map(ld.study_major) + "|" + F.batch
        F = F.assign(batch=stratum.to_numpy())
        d = LU._hlca_cov()
        cols = ("age", "sex", "disease", "lung_condition", "study_major")
    else:
        d = C.census_donors(C.ATLAS / "natri2024_ild", prefix)
        sets = LU._natri_labels(d).rename({"control": "control", "case": "case"})
        stratum_by = C.majority(F)
        F = F.assign(batch=F.resection.map(stratum_by).to_numpy())
        cols = ("age", "sex", "disease", "Diagnosis", "Status")
    F = F[F.resection.isin(sets.index)]
    pb = C.donor_table(F, None, 1, stratum="batch")
    tmin = np.where(pb.compartment.str.startswith("lung:T_"), t_min, p["min_cells"])
    pb = pb[pb.n_cells >= tmin]
    return O.finish(pb, d, {"tstate_set": sets.rename("tstate_set")}, cols=cols)


def check_twins(n: int = 8, seed: int = 0) -> list:
    """Slices of mode "T" (HLCA) and mode "noncycling" (Natri controls) against `hlca_fibrotic` and `lung_natri_control`."""
    import importlib.util

    out = []
    for name, mode, twin, only in (("hlca_tstates", "T", "hlca_fibrotic", None), ("natri_tstates", "noncycling", "lung_natri_control", "control")):
        f = C.HERE.parent / twin / "build.py"
        spec = importlib.util.spec_from_file_location(f"b_{twin}", f)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        pb_t, _ = mod.from_legacy()
        rng = np.random.default_rng(seed)
        ids = sorted(rng.permutation(sorted(pb_t.donor.unique()))[:n])
        pb_b, d_b = build(name, donors=[i.split(":", 1)[1] for i in ids], mode=mode)
        if only:
            pb_b = pb_b[pb_b.donor.isin(d_b.index[d_b.tstate_set == only])]
        a = pb_t[pb_t.donor.isin(ids)].set_index(["donor", "compartment"]).sort_index()
        b = pb_b[pb_b.donor.isin(ids)].set_index(["donor", "compartment"]).sort_index()
        c = a.index.intersection(b.index)
        num = ["n_cells", "depth"] + list(C.GENES)
        out.append(dict(dataset=name, mode=mode, twin=twin, slice_donors=len(ids), rows_twin=len(a), rows_build=len(b),
                        only_twin=len(a.index.difference(b.index)), only_build=len(b.index.difference(a.index)),
                        max_abs_diff=float(np.abs(a.loc[c, num].to_numpy(float) - b.loc[c, num].to_numpy(float)).max()),
                        strata_equal=bool((a.loc[c, "stratum"].astype(str) == b.loc[c, "stratum"].astype(str)).all())))
    return out


def from_legacy(name: str):
    raise NotImplementedError("no legacy table was stored for the T-state compartments; see check_twins")
