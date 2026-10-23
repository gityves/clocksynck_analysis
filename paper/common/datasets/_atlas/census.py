"""Shared steps of the single-cell atlas builders (paper/common/datasets/<atlas>/build.py).

Input of every builder: the group-level pseudobulk of an atlas, `$SCCC_DATA/<root>/<dataset>/pseudobulk.h5ad`
(one row per donor x tissue x cell type x assay [x source label], summed raw counts of a gene panel, `n_cells`,
`depth` = genome-wide count of the group) with its `donors.csv`. These files are the output of step 1, which the
builders do not repeat (it streams the CELLxGENE Census over the network, or a GEO matrix, in the legacy
`desync/common/build_pseudobulk.py` and `desync/14_skin|16_ibd|19_endometrium/scripts/build_*.py`).
Step 2, here, is the donor x compartment table of the clock genes:

1. each group gets a class (first matching regular expression of an ordered tuple on a label; unmatched groups are
   dropped), the compartment is `<organ>:<class>` (or `<organ>:<region>_<class>`);
2. each donor keeps the groups of its majority assay (the assay with most cells), the stratum;
3. groups are summed per donor x compartment (x stratum); tables keep pseudobulks of at least `min_cells` cells and,
   optionally, compartments with at least `min_donors` donors.

Standard table (`clocksync` input): donor, compartment, stratum, n_cells, depth, clock-gene counts. Donor ids keep
the dataset prefix of the legacy tables (`<dataset>:<donor_id>`), so that the legacy and the rebuilt tables are
comparable row by row. `donors`: one row per donor (index `donor`).

`slice`: the builders accept `donors=[donor_id, ...]` (Census donor ids without prefix) to read only these donors;
used to check a builder against the legacy table without rebuilding the dataset. Per-donor steps (majority assay,
sums) are unaffected by the slice; the `min_donors` rule is not applied to a slice.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # paper/code

import clocksync as cs  # noqa: E402
from clocksync.genes import gene_columns  # noqa: E402
from common.paths import DATA  # noqa: E402

ATLAS, ATLAS_WG, GEO = DATA / "atlas", DATA / "atlas_wg", DATA / "geo"
GENES = list(cs.gtex().genes)
TARGETS = HERE / "a17_targets.csv"


# ---- reading ---------------------------------------------------------------------------------------------------
SCORE = (HERE / "inputs" / "score_all.txt").read_text().split()  # Hallmark TNFA / IFNG score genes (A51)


def scored_path(path: Path) -> Path:
    """The rebuild of a group pseudobulk with the score genes added (A51): `$SCCC_DATA/immune_damping/<same subpath>`."""
    p = Path(path).relative_to(DATA)
    return Path(path) if p.parts[0] == "immune_damping" else DATA / "immune_damping" / p


def read(path: Path, genes="clock", donors=None, donor_col: str = "donor_id") -> tuple[pd.DataFrame, pd.DataFrame]:
    """(obs, counts): groups of `<path>/pseudobulk.h5ad`; counts columns = the 14 clock genes, with `genes="all"` every
    var (by `feature_name`, canonical symbols), with `genes="a7"` every gene of the `<dataset>_a7` rebuild (A6/A7 panel), with `genes="scored"` the clock genes followed by the score genes found in the file (read from
    the immune_damping rebuild of the same dataset), or a list of genes. `donors`: keep only these donor ids."""
    if isinstance(genes, str) and genes == "scored":
        path = scored_path(path)
    elif isinstance(genes, str) and genes == "a7":  # the rebuild with the A6/A7 panel (clock, glucocorticoid, activation genes)
        path = Path(path) if str(path).endswith("_a7") else Path(str(path) + "_a7")
        genes = "all"
    a = ad.read_h5ad(Path(path) / "pseudobulk.h5ad")
    o = a.obs
    if donors is not None:
        sel = np.flatnonzero(o[donor_col].astype(str).isin([str(d) for d in donors]).to_numpy())
        a, o = a[sel], a.obs.iloc[sel]
    names = a.var.feature_name.astype(str) if "feature_name" in a.var else a.var_names.astype(str)
    M = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    M = M.astype(float, copy=False)
    if isinstance(genes, str) and genes == "all":
        from clocksync.genes import canonical
        return o, pd.DataFrame(M, columns=[canonical(g) for g in names])
    if isinstance(genes, str) and genes == "scored":
        have = gene_columns(names, SCORE) >= 0
        genes = GENES + [g for g, h in zip(SCORE, have) if h and g not in GENES]
    elif not isinstance(genes, str):
        genes = list(genes)
    else:
        genes = GENES
    cols = gene_columns(names, genes)
    if (cols < 0).any():
        raise KeyError(f"genes missing: {np.array(genes)[cols < 0]}")
    return o, pd.DataFrame(M[:, cols], columns=genes)


def classify(labels: pd.Series, classes, ignore_case: bool = True) -> pd.Series:
    """First matching (name, regex) of `classes` on each label; NaN if none. `ignore_case=False` reproduces the
    case-sensitive `str.contains` of `run_atlas.coarse` (cortex, retina, Census-wide COARSE of `run_atlas.load`)."""
    flags = re.I if ignore_case else 0
    m = {}
    for lab in pd.Series(labels.unique()):
        m[lab] = next((name for name, pat in classes if re.search(pat, lab, flags=flags)), np.nan)
    return labels.map(m)


# ---- groups -> donor x compartment ------------------------------------------------------------------------------
def frame(o: pd.DataFrame, X: pd.DataFrame, prefix: str, compartment: pd.Series, batch: pd.Series,
          keep=None, extra: dict | None = None) -> pd.DataFrame:
    """Rows = groups kept (compartment not NaN and `keep`): resection, compartment, batch, n_cells, depth, counts."""
    ok = compartment.notna().to_numpy() & (np.ones(len(o), bool) if keep is None else np.asarray(keep))
    F = X[ok].reset_index(drop=True)
    F["_row"] = np.flatnonzero(ok)  # row of the group in the input (used to sum sparse genome-wide matrices)
    F.insert(0, "resection", (f"{prefix}:" + o.donor_id.astype(str)[ok]).to_numpy())
    F.insert(1, "compartment", compartment[ok].astype(str).to_numpy())
    F.insert(2, "batch", pd.Series(batch).astype(str)[ok].to_numpy())
    F["n_cells"], F["depth"] = o.n_cells.to_numpy()[ok], o.depth.to_numpy(float)[ok]
    for k, v in (extra or {}).items():
        F[k] = pd.Series(np.asarray(v))[ok].to_numpy()
    return F


def majority(F: pd.DataFrame, col: str = "batch") -> pd.Series:
    """Per resection, the value of `col` with most cells (ties: last in sort order, as the legacy code)."""
    return (F.groupby(["resection", col]).n_cells.sum().reset_index().sort_values("n_cells")
            .drop_duplicates("resection", keep="last").set_index("resection")[col])


def keep_majority(F: pd.DataFrame, col: str = "batch") -> pd.DataFrame:
    return F[F[col].to_numpy() == majority(F, col).reindex(F.resection).to_numpy()]


def donor_table(F: pd.DataFrame, keep=None, min_cells: int = 20, stratum: pd.Series | str | None = "batch",
                min_donors: int | None = None, sum_over: str | None = None) -> pd.DataFrame:
    """Groups summed per resection x compartment (x batch). `stratum`: "batch" (the column), a Series
    resection -> label (replaces it), or None. `keep`: resections to use. Counts and n_cells and depth are summed."""
    if keep is not None:
        F = F[F.resection.isin(keep)]
    by_batch = isinstance(stratum, str) and stratum == "batch"
    cols = ["resection", "compartment"] + (["batch"] if by_batch else [])
    num = [c for c in F.columns if c not in ("resection", "compartment", "batch") and not c.startswith("_") and
           pd.api.types.is_numeric_dtype(F[c])]
    pb = F.groupby(cols, as_index=False)[num].sum()
    if by_batch:
        pb["stratum"] = pb.pop("batch")
    elif stratum is not None:
        pb["stratum"] = pb.resection.map(stratum)
    pb = pb[pb.n_cells >= min_cells]
    if min_donors:
        n = pb.groupby("compartment").resection.nunique()
        pb = pb[pb.compartment.isin(n[n >= min_donors].index)]
    return pb.reset_index(drop=True)


def standard(pb: pd.DataFrame) -> pd.DataFrame:
    """Standard column order and names: donor, compartment, [stratum], n_cells, depth, clock genes."""
    pb = pb.rename(columns={"resection": "donor"})
    head = ["donor", "compartment"] + (["stratum"] if "stratum" in pb else []) + ["n_cells", "depth"]
    genes = [c for c in pb.columns if c not in head]
    out = pb[head + genes].copy()
    out["n_cells"] = out.n_cells.astype(int)
    if "stratum" in out:
        out["stratum"] = out.stratum.astype(str)
    return out.sort_values(["donor", "compartment"]).reset_index(drop=True)


# ---- donor tables ------------------------------------------------------------------------------------------------
ORD = {w: i + 1 for i, w in enumerate(("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth",
                                       "ninth", "tenth"))}
SEAAD_BIN = {"90+ years old": 92.0, "78 to 89 years old": 83.5, "65 to 77 years old": 71.0,
             "Less than 65 years old": 60.0}


def stage_age(s, seaad_bin: str | None = None) -> float:
    """Age in years from a Census `development_stage` (A19 rule, `desync/20_age/scripts/ages.py`)."""
    s = str(s)
    if m := re.match(r"(\d+)-year-old", s):
        return float(m.group(1))
    if s.startswith("80 year-old and over"):
        return SEAAD_BIN.get(seaad_bin, 85.0) if seaad_bin else 85.0
    if s.startswith("90 year-old and over"):
        return 92.0
    if m := re.match(r"(\w+) decade", s):
        n = ORD.get(m.group(1))
        return (n - 1) * 10 + 4.5 if n else np.nan
    return np.nan


def census_donors(path: Path, prefix: str) -> pd.DataFrame:
    """`donors.csv` of a Census build: index `<prefix>:<donor_id>`, columns age (years, A19 rule), sex, disease and
    the file's other donor-level fields under their own names."""
    d = pd.read_csv(Path(path) / "donors.csv", dtype={"donor_id": str})
    binc = d["Age at death"] if "Age at death" in d else pd.Series(None, index=d.index, dtype=object)
    d["age"] = [stage_age(s, b if isinstance(b, str) else None) for s, b in
                zip(d.development_stage if "development_stage" in d else [""] * len(d), binc)]
    d.index = f"{prefix}:" + d.donor_id
    d.index.name = "donor"
    return d


def donors_frame(d: pd.DataFrame, donors, extra: dict | None = None, cols=("age", "sex", "disease")) -> pd.DataFrame:
    """Donor table for the donors of `donors` (index): `cols` of `d` and `extra` columns (Series by donor)."""
    out = d.reindex(pd.Index(donors, name="donor"))[[c for c in cols if c in d]].copy()
    for k, v in (extra or {}).items():
        out[k] = v.reindex(out.index)
    return out


# ---- A17 target metadata ------------------------------------------------------------------------------------------
def targets(dataset: str) -> pd.DataFrame:
    """A17 targets (healthy fits and contrasts of the re-analysis) that read this dataset: name, kind, study,
    tissue, label (amendment), ref (reference compartment of the generation-1 fit), groups, within (permutation
    blocks: 'stratum', 'block:<column of donors>' or empty), primary (measures of the pre-specified test),
    donors, compartments, cases, controls."""
    t = pd.read_csv(TARGETS, keep_default_na=False)
    return t[t.dataset == dataset].reset_index(drop=True)


def read_sparse(path: Path, donors=None, donor_col: str = "donor_id"):
    """(obs, X, names) of a genome-wide group pseudobulk: X sparse (groups x genes, float64), names = canonical gene
    symbols (`feature_name`); duplicated symbols keep their first column."""
    from clocksync.genes import canonical
    import scipy.sparse as sp

    a = ad.read_h5ad(Path(path) / "pseudobulk.h5ad")
    o, X = a.obs, a.X
    if donors is not None:
        sel = np.flatnonzero(o[donor_col].astype(str).isin([str(d) for d in donors]).to_numpy())
        o, X = o.iloc[sel], X[sel]
    names = np.array([canonical(g) for g in a.var.feature_name.astype(str)])
    first = ~pd.Index(names).duplicated()
    X = sp.csr_matrix(X).astype(np.float64)[:, np.flatnonzero(first)]
    return o, X, names[first]


def sum_groups_sparse(F: pd.DataFrame, X, names, keep=None, min_cells: int = 20, stratum="batch",
                      min_donors: int | None = None) -> pd.DataFrame:
    """`donor_table` for a sparse genome-wide matrix: F = group rows with `_row` (as `frame`), X = the input matrix.
    Genes without a count in the table are dropped (the 14 clock genes stay)."""
    import scipy.sparse as sp

    if keep is not None:
        F = F[F.resection.isin(keep)]
    F = F.reset_index(drop=True)
    by_batch = isinstance(stratum, str) and stratum == "batch"
    cols = ["resection", "compartment"] + (["batch"] if by_batch else [])
    codes, uniq = pd.MultiIndex.from_frame(F[cols]).factorize(sort=True)
    S = sp.csr_matrix((np.ones(len(F)), (codes, np.arange(len(F)))), shape=(len(uniq), len(F)))
    G = (S @ X[F["_row"].to_numpy()]).tocsr()
    pb = pd.DataFrame(list(uniq), columns=cols)
    pb["n_cells"] = np.bincount(codes, weights=F.n_cells.to_numpy(float), minlength=len(uniq)).astype(int)
    pb["depth"] = np.bincount(codes, weights=F.depth.to_numpy(float), minlength=len(uniq))
    if by_batch:
        pb["stratum"] = pb.pop("batch")
    elif stratum is not None:
        pb["stratum"] = pb.resection.map(stratum)
    ok = (pb.n_cells >= min_cells).to_numpy().copy()
    if min_donors:
        n = pb[ok].groupby("compartment").resection.nunique()
        ok &= pb.compartment.isin(n[n >= min_donors].index).to_numpy()
    G = G[np.flatnonzero(ok)]
    pb = pb[ok].reset_index(drop=True)
    nz = np.asarray(G.sum(0)).ravel() > 0
    used = np.flatnonzero(nz | np.isin(np.asarray(names), GENES))  # the clock genes stay even without a count
    counts = pd.DataFrame(G[:, used].toarray(), columns=np.asarray(names)[used])
    return pd.concat([pb, counts], axis=1)


def canon(g: str) -> str:
    """Canonical gene symbol (clock-gene aliases such as BMAL1 -> ARNTL)."""
    from clocksync.genes import canonical
    return canonical(g)
