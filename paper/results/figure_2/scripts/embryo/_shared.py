"""Helpers shared by embryo and natural_experiments (gate rows, Holm over gates, eligibility
rules, the mouse-embryo windows). Method calls are clocksync's; everything marked "local: not in clocksync" is a
small helper for a statistic the package has no public call for (listed in the folder READMEs)."""

from __future__ import annotations

import copy
import dataclasses
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.timed import W  # noqa: E402
from common.paths import DATA  # noqa: E402

PASS = 0.05


# ---- eligibility and gates ---------------------------------------------------------------------------------------

def eligible(pb: pd.DataFrame, min_cells: int = 0, min_donors: int = 0, donors=None, compartments=None) -> pd.DataFrame:
    """Rows of `donors` / `compartments` with >= `min_cells` cells; then compartments with >= `min_donors` donors
    (the order of the legacy `organ_fit.donor_table`: cells first, donors second)."""
    t = pb
    if donors is not None:
        t = t[t.donor.isin(list(donors))]
    if compartments is not None:
        t = t[t.compartment.isin(list(compartments))]
    if min_cells:
        t = t[t.n_cells >= min_cells]
    if min_donors:
        n = t.groupby("compartment").donor.nunique()
        t = t[t.compartment.isin(n.index[n >= min_donors])]
    return t.reset_index(drop=True)


def dataset(pb: pd.DataFrame, donors: pd.DataFrame | None = None, anchors=None, name: str = "") -> cs.Dataset:
    """Dataset whose host is fitted on `anchors` (default: every compartment, as the legacy gates and the baboon
    tests; the gate itself does not depend on the anchors)."""
    return cs.Dataset(pb, anchors=sorted(pb.compartment.unique()) if anchors is None else list(anchors),
                      donors=donors, name=name)


def gate_row(g: cs.Gate, pb: pd.DataFrame | None = None, **extra) -> dict:
    """One output row of a gate; pass = p_pairing and p_acrophase < 0.05 (amendment rule)."""
    row = dict(extra)
    if pb is not None:
        per = pb.groupby("donor").compartment.nunique()
        row.update(donors=pb.donor.nunique(), donors_3plus=int((per >= 3).sum()),
                   compartments=";".join(sorted(pb.compartment.unique())),
                   strata=";".join(sorted(pb.stratum.unique())) if "stratum" in pb else "")
    row.update(p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share,
               share_null_mean=g.share_null_mean, clock_cov=g.clock_cov, n_pairs=g.n_pairs, passed=g.passes)
    return row


def run_gate(pb: pd.DataFrame, test: str, seed: int, n_perm: int, **extra) -> dict:
    """Gate over all compartments of `pb`."""
    return gate_row(dataset(pb).gate(n_perm=n_perm, seed=seed), pb, test=test, **extra)


def holm_max(t: pd.DataFrame, mask=None) -> pd.DataFrame:
    """p_max = max(p_pairing, p_acrophase) and its Holm adjustment over the rows of `mask` (A105b, A107b)."""
    t = t.copy()
    t["p_max"] = t[["p_pairing", "p_acrophase"]].max(axis=1)
    m = np.ones(len(t), bool) if mask is None else np.asarray(mask, bool)
    t["p_holm"] = np.nan
    t.loc[m, "p_holm"] = cs.stats.holm(t.loc[m, "p_max"].to_numpy())
    return t


# ---- mouse embryo (Qiu 2024) ---------------------------------------------------------------------------------------

WINDOWS = {"E8.0-E10.0": (8.0, 10.0), "E10.25-E12.75": (10.25, 12.75), "E13.0-E14.75": (13.0, 14.75),
           "E15.0-E18.75": (15.0, 18.75)}
LATE = WINDOWS["E15.0-E18.75"]
BLOOD_LIVER = {"Hepatocytes", "Megakaryocytes", "Definitive_erythroid", "Primitive_erythroid", "White_blood_cells",
               "B_cells", "T_cells", "Mast_cells"}
CNS = {"CNS_neurons", "Neuroectoderm_and_glia", "Intermediate_neuronal_progenitors", "Oligodendrocytes",
       "Ependymal_cells"}
QIU = DATA / "natexp" / "qiu2024"


def embryos(donors: pd.DataFrame, lo: float, hi: float) -> pd.Index:
    return donors.index[(donors.day >= lo - 1e-9) & (donors.day <= hi + 1e-9)]


def window(pb: pd.DataFrame, donors: pd.DataFrame, lo: float, hi: float, min_cells: int = 50, frac: float = 0.75):
    """Embryos of the window; cell types with >= `min_cells` cells present in >= `frac` of the window's embryos
    (A108). The stratum (batch) is already a column of the table: it only has an effect when the window holds
    several batches, as pre-registered."""
    don = embryos(donors, lo, hi)
    t = eligible(pb, min_cells=min_cells, donors=don)
    n = t.groupby("compartment").donor.nunique()
    t = t[t.compartment.isin(n.index[n >= np.ceil(frac * len(don))])].reset_index(drop=True)
    return t, don


def trajectory_of() -> pd.Series:
    """cell type -> trajectory (lineage), from the legacy pseudobulk annotation. local: not in clocksync (the
    standard table keeps no cell-type to trajectory map)."""
    import anndata as ad
    obs = ad.read_h5ad(QIU / "pseudobulk_celltype" / "pseudobulk.h5ad", backed="r").obs
    return obs.drop_duplicates("cell_type").set_index("cell_type").trajectory.astype(str)


def qiu_panel(level: str = "celltype") -> pd.DataFrame:
    """Embryo x compartment pseudobulk with all 67 panel genes (14 clock genes and 53 others) of the legacy run.
    local: not in clocksync / the standard table (the builder keeps only the 14 clock genes); used by the
    pseudo-clock controls and the glucocorticoid score."""
    import anndata as ad
    import scipy.sparse as sp
    a = ad.read_h5ad(QIU / f"pseudobulk_{level}" / "pseudobulk.h5ad")
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    t = pd.DataFrame(X.astype(float), columns=[canonical(g) for g in a.var_names])
    t.insert(0, "depth", a.obs.depth.to_numpy(float))
    t.insert(0, "n_cells", a.obs.n_cells.to_numpy(int))
    t.insert(0, "stratum", a.obs.batch.astype(str).to_numpy())
    t.insert(0, "compartment", a.obs.cell_type.astype(str).to_numpy())
    t.insert(0, "donor", a.obs.donor_id.astype(str).to_numpy())
    return t


def residualised(ds: cs.Dataset, x: pd.Series) -> cs.Dataset:
    """Dataset whose clock-plane states are residualised on a donor covariate, linearly per compartment, before the
    host is fitted (removes a smooth drift, A108 sensitivity). local: not in clocksync (it replaces the states of
    `ds.data`)."""
    d = ds.data
    b = d.beta.copy()
    xv = x.reindex(d.resection).to_numpy()
    for k in np.unique(d.compartment):
        i = np.flatnonzero(d.compartment == k)
        M = np.column_stack([np.ones(len(i)), xv[i] - xv[i].mean()])
        B, *_ = np.linalg.lstsq(M, b[i], rcond=None)
        b[i] = b[i] - M[:, 1:] @ B[1:]
    new = copy.copy(ds)
    new.data = dataclasses.replace(d, beta=b)
    return new


# clock genes and strong clock- or glucocorticoid-driven genes, never drawn as pseudo-clock genes (as clocksync's list)
CLOCK_LIKE = set(cs.gtex().genes) | {
    "ARNTL2", "BHLHE40", "BHLHE41", "NFIL3", "HLF", "RORB", "RORC", "CIPC", "TIMELESS", "CSNK1D", "CSNK1E",
    "USP2", "TSPAN4", "WEE1", "TSC22D3", "KLF9", "KLF15", "ZBTB16", "FKBP5", "PDK4", "TXNIP", "SGK1", "DDIT4"}


def matched_sets(lm: pd.Series, genes, n_sets: int, rng) -> list[list[str]]:
    """Pseudo-clock sets from the panel's non-clock genes: each clock gene (random order) takes the unused non-clock
    gene closest in log mean expression (the 67-gene panel is too small for the quantile bins of
    `pseudo_clock_genes`; A108 deviation, documented in the amendment). local: not in clocksync."""
    pool0 = lm[[canonical(g) not in CLOCK_LIKE for g in lm.index]].dropna()
    out = []
    for _ in range(n_sets):
        pool, s = pool0.copy(), {}
        for g in rng.permutation(list(genes)):
            c = (pool - lm[g]).abs().idxmin()
            s[g] = c
            pool = pool.drop(c)
        out.append([s[g] for g in genes])
    return out


def complex_cosinor(y, hours, covariates=None, groups=None) -> complex:
    """Complex coefficient a of y = mu(+ group) + b x + a exp(i w t) for a complex clock-plane state y = beta_x +
    i beta_y (the rotating-wave component; amplitude |a|, acrophase arg a). local: not in clocksync
    (`timed.cosinor` is real-valued, i.e. one coordinate at a time)."""
    y = np.asarray(y, complex)
    n = len(y)
    D = [np.ones((n, 1))] if groups is None else [pd.get_dummies(pd.Series(groups).astype(str)).to_numpy(float)]
    if covariates is not None:
        c = np.asarray(covariates, float).reshape(n, -1)
        D.append(c - c.mean(0))
    M = np.column_stack(D + [np.exp(1j * W * np.asarray(hours, float))])
    return complex(np.linalg.lstsq(M, y, rcond=None)[0][-1])


def timed_wave(y, hours, n_perm: int, rng, covariates=None, groups=None, blocks=None) -> dict:
    """|a| of the complex cosinor with a permutation p (hours permuted over samples, within `blocks` when given;
    the add-one p of `clocksync.stats.perm_p`)."""
    hours = np.asarray(hours, float)
    a = complex_cosinor(y, hours, covariates, groups)
    parts = [np.arange(len(hours))] if blocks is None else \
        [np.flatnonzero(np.asarray(blocks) == g) for g in np.unique(np.asarray(blocks))]
    null = []
    for _ in range(n_perm):
        h = hours.copy()
        for i in parts:
            h[i] = hours[rng.permutation(i)]
        null.append(abs(complex_cosinor(y, h, covariates, groups)))
    return dict(amp=abs(a), phase_h=(np.angle(a) / W) % 24, p=cs.stats.perm_p(null, abs(a), "greater", tie=0.0))


# ---- placenta (Arutyunyan 2023): maternal soup and the ambient check -------------------------------------------------

MAT = ["mat:dS", "mat:dNK", "mat:myeloid", "mat:T", "mat:endothelial", "mat:epithelial", "mat:perivascular"]
ANCH_PLACENTA = ["mat:dS", "mat:endothelial", "mat:perivascular"]
FET = ["fet:VCT", "fet:SCT", "fet:EVT", "fet:Hofbauer", "fet:fibroblast"]
PRIMARY_FET = ["fet:VCT", "fet:SCT"]


def placenta_table(pb: pd.DataFrame, min_cells: int = 20, min_donors: int = 8) -> pd.DataFrame:
    """Maternal host and fetal compartments with >= 20 cells and >= 8 pregnancies (A106)."""
    return eligible(pb, min_cells=min_cells, min_donors=min_donors, compartments=MAT + FET)


def soup_share(pb: pd.DataFrame, donors: pd.DataFrame, min_cells: int = 20) -> pd.Series:
    """Maternal soup share per fetal compartment: XIST per count in a fetal compartment of male fetuses (which cannot
    express it) relative to the pregnancy's maternal median (A106, build_fetal.md). Median over pregnancies.
    local: not in clocksync."""
    male = donors.index[donors.fetal_sex_inferred == "male"]
    a = pb[pb.donor.isin(male) & (pb.n_cells >= min_cells) & pb.compartment.str.match(r"^(mat|fet):")].copy()
    a["f"] = a.XIST / a.depth
    mat = a[a.compartment.str.startswith("mat:")].groupby("donor").f.median()
    fet = a[a.compartment.str.startswith("fet:")].copy()
    fet["share"] = fet.f / fet.donor.map(mat)
    return fet.groupby("compartment").share.median()


def synthetic_ambient(pb: pd.DataFrame, targets, share: pd.Series, factor: float = 10.0) -> pd.DataFrame:
    """Rows of `targets` replaced by (1 - a) x the across-pregnancy mean profile of the compartment + a x the
    pregnancy's maternal soup profile (sum of its maternal compartments), a = factor x soup share, both scaled to the
    pseudobulk's depth and rounded: fetal cells without a clock of their own (A106 ambient check). local: not in
    clocksync."""
    g = list(cs.gtex().genes)
    t = pb.copy()
    mat = t[t.compartment.isin(MAT)]
    soup = mat.groupby("donor")[g].sum().div(mat.groupby("donor").depth.sum(), axis=0)
    for k in targets:
        i = t.compartment == k
        prof = t.loc[i, g].sum() / t.loc[i, "depth"].sum()
        a = min(factor * share.get(k, np.nan), 1.0)
        s = soup.reindex(t.loc[i, "donor"]).to_numpy()
        s = np.where(np.isnan(s), prof.to_numpy()[None, :], s)
        mix = (1 - a) * prof.to_numpy()[None, :] + a * s
        t.loc[i, g] = np.round(mix * t.loc[i, "depth"].to_numpy()[:, None])
    return t
