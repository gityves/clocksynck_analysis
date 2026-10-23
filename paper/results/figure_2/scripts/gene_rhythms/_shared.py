"""Helpers shared by the scripts of gene_rhythms: genome-wide atlas tables, GTEx per-gene truths, untimed gene
rhythms of cell types on the host (`Dataset.gene_rhythms`), their recombination into a bulk rhythm and the scores against
the timed truth. Clock quantities come from `clocksync`; what is local is marked "local: not in clocksync" and listed in the
README. GTEx time of death is controlled access: callers write aggregates (per gene and tissue), never per-donor rows.
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from common import data  # noqa: E402
from common.anchors import anchors as a17_anchors  # noqa: E402
from common.gtex import LABELLED, FAST  # noqa: E402


HERE = Path(__file__).resolve().parent
W = 2 * np.pi / 24
T0 = 13.0  # rise / fall is read at 13 h of the host frame (A31b)
MIN_MEAN = 5.0
META = ("donor", "compartment", "stratum", "n_cells", "depth")
FRAME = {"ARNTL", "BMAL1", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "RORA", "TEF",
         "CIART"}  # the clock-plane genes (legacy compare.FRAME)
FOLLOWERS = ("ADGRA2 ADIRF-AS1 ARHGAP42 ARNTL2 BHLHE40 BHLHE41 CDC25B CDK18 CHST7 CIPC DDIT4L EPB41L1 FAM76A FOXO1 "
             "GUCY1A1 HLF JAM2 KLF10 LDB1 LONRF3 LPCAT3 MAMDC2 MPP2 NAMPT NFIL3 PPARA RAB3IL1 RORC SUN2 TSPAN4 USP2 "
             "WEE1").split()  # the 32 frozen follower genes of D10 (A28)
assert len(FOLLOWERS) == 32

# A27 atlases (clock-gene dataset names are the A17 targets) and A31b atlases: genome table, anchors, GTEx tissues
ATLASES = {
    "heart": dict(genome="heart_global_genome", anchors="heart", tissues=["Heart - Left Ventricle", "Heart - Atrial Appendage"]),
    "muscle": dict(genome="organ_muscle_genome", anchors="organ_muscle", tissues=["Muscle - Skeletal"]),
    "breast": dict(genome="organ_breast_genome", anchors="organ_breast", tissues=["Breast - Mammary Tissue"]),
    "colon": dict(genome="organ_colon_genome", anchors="organ_colon", tissues=["Colon - Transverse", "Colon - Sigmoid"]),
    "lung_hlca": dict(genome="lung_hlca_genome", anchors="lung_hlca", tissues=["Lung"], new=True),
    "breast_hbca": dict(genome="breast_hbca_genome", anchors="breast_hbca", tissues=["Breast - Mammary Tissue"], new=True),
    "liver_xu": dict(genome="liver_xu_genome", anchors="liver_xu", tissues=["Liver"], new=True),
    "skin_ganier": dict(genome="skin_ganier_genome", anchors="skin_ganier",
                        tissues=["Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"], new=True),
    "adipose_emont": dict(genome="organ_adipose_genome", anchors="adipose_emont", tissues=["Adipose - Subcutaneous"], new=True),
}
A27 = ["heart", "muscle", "breast", "colon"]
A31B = ["lung_hlca", "breast_hbca", "liver_xu", "skin_ganier", "adipose_emont"]
POOL_TISSUE = {k: v["tissues"][0] for k, v in ATLASES.items()}
TRUTH_TISSUES = sorted({t for v in ATLASES.values() for t in v["tissues"]})
BRAIN_TISSUES = ["Brain - Frontal Cortex (BA9)", "Brain - Cortex"]
CLOCK = list(cs.gtex().genes)


# ---- genome-wide atlas tables -------------------------------------------------------------------------------------

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS, PAPER_DATA  # noqa: E402
def _patch_census():
    """Workaround for the pandas 3 read-only array in `data/_atlas/census.sum_groups_sparse` (an in-place `&=` on a view)."""
    sys.path.insert(0, str(DATASETS / "_atlas"))
    import census as C

    if getattr(C, "_patched", False):
        return
    src = inspect.getsource(C.sum_groups_sparse).replace("ok = (pb.n_cells >= min_cells).to_numpy()",
                                                         "ok = (pb.n_cells >= min_cells).to_numpy().copy()")
    ns = dict(C.__dict__)
    exec(src, ns)  # noqa: S102
    C.sum_groups_sparse, C._patched = ns["sum_groups_sparse"], True


def genome_table(name: str, subsample: int | None = None, seed: int = 0):
    """(pseudobulk of every gene, donors) of a `*_genome` / `gene_map_*` dataset, built from the public group pseudobulks
    (a few seconds to a minute; no cache). `subsample`: N random donors."""
    _patch_census()
    pb, donors = data._builder(name).build()
    pb["donor"] = pb["donor"].astype(str)
    donors.index = donors.index.astype(str)
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb, donors = pb[pb.donor.isin(keep)].reset_index(drop=True), donors.loc[donors.index.isin(keep)]
    return pb, donors


def class_rule(pb: pd.DataFrame, min_donors: int = 20) -> pd.DataFrame:
    """Compartments seen in at least `min_donors` donors (the class rule of the atlases, applied within an analysis set)."""
    n = pb.groupby("compartment").donor.nunique()
    return pb[pb.compartment.isin(n.index[n >= min_donors])].reset_index(drop=True)


def gene_columns(pb: pd.DataFrame) -> list:
    return [c for c in pb.columns if c not in META]


def a31b_anchors(atlas: str, comps) -> list:
    """Anchors of a new A31b atlas: the learned rule of A17 re-run on it (frozen copy `anchors_a31b.csv`), fixed rule if < 2."""
    a = pd.read_csv(HERE / "anchors_a31b.csv")
    a = a[a.atlas == atlas]
    got = [c for c in a.compartment[a.anchor] if c in set(comps)]
    return got if len(got) >= 2 else [c for c in a.compartment[a.anchor_fixed] if c in set(comps)]


def atlas_anchors(atlas: str, comps) -> list:
    spec = ATLASES[atlas]
    if spec.get("new") and atlas not in ("lung_hlca", "adipose_emont"):
        return a31b_anchors(atlas, comps)
    if atlas == "adipose_emont":
        return [c for c in a31b_anchors(atlas, comps)] or [c for c in a17_anchors("organ_adipose") if c in set(comps)]
    got = [c for c in a17_anchors(spec["anchors"]) if c in set(comps)]
    return got if len(got) >= 2 else [c for c in a17_anchors(spec["anchors"], "fixed") if c in set(comps)]


def make_dataset(pb: pd.DataFrame, anchors: list, name: str = "") -> cs.Dataset:
    """Clock-gene `Dataset` of a genome table (only the 14 clock columns enter the host)."""
    return cs.Dataset(pb[[c for c in pb.columns if c in META or c in CLOCK]], anchors, name=name)


# ---- GTEx per-gene truth (local: not in clocksync) ------------------------------------------------------------

def _dip() -> pd.Series:
    d = pd.read_csv(PAPER_DATA / "chiral_DIP.csv", sep=";", decimal=",")
    d.index = "GTEX-" + d["SUBJ.ID"].astype(str)
    return d.hour


def gtex_truth(tissue: str, version: str = "v1", n: int | None = None, seed: int = 0, hardy=None) -> pd.DataFrame:
    """Per-gene timed cosinor of one GTEx tissue: log(CPM + 0.5) ~ 1 [+ covariates] + cos wt + sin wt (w = 2 pi / 24).

    v1: fast deaths (Hardy 1-2) with a time of death. v2: v1 with RIN > 6 and covariates (ischaemic time and its square,
    sex, age bracket, Hardy class). v3: all death types, the CHIRAL internal phase (DIP) replaces the time of death, as v2.
    Genes with median CPM >= 5. Returns gene, n, amp, acro_h, b_cos, b_sin, se, p (2-df F), q (BH). Aggregates only.
    `n`: use n random samples (checks)."""
    import anndata as ad

    a = ad.read_h5ad(LABELLED, backed="r")
    try:
        o = a.obs
        hardy = list(FAST) if hardy is None else list(hardy)
        sel = (o.tissue.astype(str) == tissue).to_numpy().copy()
        if version == "v3":
            tod = o.subject.astype(str).map(_dip()).to_numpy(float)
            sel &= o.DTHHRDY.notna().to_numpy() & ~np.isnan(tod) & (o.rin > 6).to_numpy() & o.ischemic_time.notna().to_numpy()
        else:
            tod = o.tod.to_numpy(float).copy()
            sel &= o.DTHHRDY.isin(hardy).to_numpy() & ~np.isnan(tod)
            if version == "v2":
                sel &= (o.rin > 6).to_numpy() & o.ischemic_time.notna().to_numpy()
        idx = np.flatnonzero(sel)
        if n is not None and n < len(idx):
            idx = np.sort(np.random.default_rng(seed).permutation(idx)[:n])
        X = a.X[idx]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        ob, genes, t = o.iloc[idx], np.asarray(a.var_names), tod[idx] % 24
    finally:
        a.file.close()
    cpm = X / ob.library_size.to_numpy(float)[:, None] * 1e6
    ok = np.median(cpm, 0) >= 5
    Y = np.log(cpm[:, ok] + 0.5)
    cols = [np.ones(len(t)), np.cos(W * t), np.sin(W * t)]
    if version in ("v2", "v3"):
        Z = pd.get_dummies(pd.DataFrame(dict(sex=ob.sex.astype(str), age=ob.age_bracket.astype(str),
                                             hardy=ob.DTHHRDY.astype(str))), drop_first=True).to_numpy(float)
        isch = ob.ischemic_time.to_numpy(float)
        isch = (isch - isch.mean()) / isch.std()
        cols += [isch, isch ** 2] + list(Z.T)
    D = np.column_stack(cols)
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    res = Y - D @ B
    dfr = len(t) - D.shape[1]
    s2 = (res ** 2).sum(0) / dfr
    D0 = np.delete(D, [1, 2], axis=1)
    B0, *_ = np.linalg.lstsq(D0, Y, rcond=None)
    F = ((((Y - D0 @ B0) ** 2).sum(0) - (res ** 2).sum(0)) / 2) / s2
    p = sps.f.sf(F, 2, dfr)
    out = pd.DataFrame(dict(tissue=tissue, gene=genes[ok], n=len(t), amp=np.hypot(B[1], B[2]),
                            acro_h=(np.arctan2(B[2], B[1]) / W) % 24, b_cos=B[1], b_sin=B[2],
                            se=np.sqrt(s2 * np.linalg.inv(D.T @ D)[1, 1]), p=p))
    out["q"] = cs.stats.bh(out.p.to_numpy())
    return out


def truth_c(g: pd.DataFrame) -> pd.DataFrame:
    g = g.drop_duplicates("gene").set_index("gene")
    g["c"] = g.b_cos + 1j * g.b_sin
    return g


# ---- recombination of cell-type loadings into a bulk rhythm (A27; local: not in clocksync) ---------------------------

def shares(pb: pd.DataFrame, comps, genes) -> pd.DataFrame:
    """Share of each gene's tissue RNA that comes from each class: mean cell fraction x mean count per cell (rows sum to 1)."""
    nc = pb.pivot_table(index="donor", columns="compartment", values="n_cells", aggfunc="sum").fillna(0)
    frac = nc.div(nc.sum(1), axis=0)
    sh = {}
    for k in comps:
        q = pb[pb.compartment == k]
        sh[k] = frac[k].mean() * q[genes].mean().to_numpy() / q.n_cells.mean()
    sh = pd.DataFrame(sh, index=list(genes))
    return sh.div(sh.sum(1), axis=0)


def predict(rhy: pd.DataFrame, sh: pd.DataFrame) -> pd.Series:
    """Predicted bulk loading per gene (complex): share-weighted mean of the class loadings of `gene_rhythms`; genes with
    a covered share <= 0.5 dropped."""
    c = (rhy.pivot(index="gene", columns="compartment", values="c_re")
         + 1j * rhy.pivot(index="gene", columns="compartment", values="c_im"))
    G = c.index
    wv = sh.reindex(G)[c.columns].where(c.notna()).fillna(0.0)
    cov = wv.sum(1)
    pe = (wv * c.fillna(0)).sum(1) / cov
    return pe[cov > 0.5]


def ccorr(a, b) -> float:
    return float(cs.stats.js_corr(a, b))


def st(a, b):
    """(R, offset_h, rise / fall agreement at 13 h) of predicted vs true phases (radians)."""
    R = np.exp(1j * (a - b)).mean()
    sign = (np.sign(np.sin(a - W * T0)) == np.sign(np.sin(b - W * T0))).mean()
    return float(abs(R)), float(np.angle(R) / W % 24), float(sign)


def score(pred_c, true_c_) -> dict:
    a, b = np.angle(pred_c), np.angle(true_c_)
    rot = np.angle(np.exp(1j * (b - a)).mean())
    err = np.abs(np.angle(np.exp(1j * (a + rot - b)))) / W
    R, off, sign = st(a, b)
    return dict(n=len(a), R_err=R, R_rot_h=off, ccorr=ccorr(a, b), rot_h=rot / W, med_err_h=float(np.median(err)),
                within3h=float(np.mean(err <= 3)), amp_spearman=float(sps.spearmanr(np.abs(pred_c), np.abs(true_c_))[0]),
                sign=sign)


def gene_set(g: pd.DataFrame, which: str):
    if which == "followers":
        return g.index[g.index.isin(FOLLOWERS) & (g.p < 0.05)]
    if which == "p<0.01 non-frame":
        return g.index[(g.p < 0.01) & ~g.index.isin(FRAME)]
    if which == "frame genes":
        return g.index[g.index.isin(FRAME)]
    raise KeyError(which)


def score_vs_null(pred: pd.Series, nulls: list, g: pd.DataFrame, genes) -> dict:
    """Score of the prediction on `genes` of truth `g` (index gene, complex c) and the host-permutation null (list of
    predictions): R with its null 95 % quantile and p = (1 + #{null >= R}) / (1 + B); same for the rise / fall agreement."""
    gs = [x for x in genes if x in pred.index]
    row = dict(n=len(gs))
    if len(gs) < 5:
        return row
    row.update(score(pred.loc[gs].to_numpy(), g.loc[gs, "c"].to_numpy()))
    if nulls:
        ns = np.array([st(np.angle(q.reindex(gs).dropna().to_numpy()), np.angle(g.loc[q.reindex(gs).dropna().index, "c"].to_numpy()))
                        for q in nulls])
        row.update(null_R_95=float(np.quantile(ns[:, 0], .95)), null_R_median=float(np.median(ns[:, 0])),
                   p_R=(1 + (ns[:, 0] >= row["R_err"]).sum()) / (len(ns) + 1), null_sign_95=float(np.quantile(ns[:, 2], .95)),
                   p_sign=(1 + (ns[:, 2] >= row["sign"]).sum()) / (len(ns) + 1), n_null=len(ns))
    return row


def fit_atlas(pb: pd.DataFrame, anchors: list, genes: list, n_null: int, seed: int = 0, name: str = ""):
    """Untimed gene rhythms of every class of a genome table, recombined into a bulk prediction.

    Returns (prediction, [null predictions], per-class table, shares). `Dataset.gene_rhythms` does the per-class fits (host
    states from the clock genes, one intercept per stratum); the recombination `predict` is local."""
    ds = make_dataset(pb, anchors, name)
    rhy = ds.gene_rhythms(pb[list(META) + [g for g in genes if g in pb.columns]], genes=[g for g in genes if g in pb.columns],
                          min_mean=MIN_MEAN, n_null=n_null, seed=seed)
    sh = shares(pb, ds.compartments, [g for g in genes if g in pb.columns])
    main = rhy[rhy["null"] == -1]
    nulls = [predict(rhy[rhy["null"] == b], sh) for b in range(n_null)]
    return predict(main, sh), nulls, main, sh, ds
