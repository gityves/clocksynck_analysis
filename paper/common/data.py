"""Standard pseudobulk tables, one per dataset: `PSEUDOBULK/<name>.parquet` (one row per donor x compartment:
donor, compartment, [stratum], [n_cells], depth, clock-gene counts) and `PSEUDOBULK/<name>_donors.parquet` (one
row per donor: groups, covariates, sampling times; controlled-access columns never leave the machine).

Each dataset has a builder `paper/common/datasets/<name>/build.py` with `build() -> (pb, donors)` from the public raw
files. Until the data are rebuilt, a builder may also provide `from_legacy() -> (pb, donors)`, which converts the
table the legacy code built; `load` uses it when the standard table is missing.
"""

from __future__ import annotations

import importlib.util

import numpy as np
import pandas as pd

from common.paths import DATASETS, PSEUDOBULK  # noqa: E402
def _builder(name: str):
    f = DATASETS / name / "build.py"
    if not f.exists():
        raise FileNotFoundError(f"no builder for dataset {name!r} ({f})")
    spec = importlib.util.spec_from_file_location(f"build_{name}", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _key(name: str, kw: dict) -> str:
    return name + "".join(f"__{k}-{v}" for k, v in sorted(kw.items()))


def save(name: str, pb: pd.DataFrame, donors: pd.DataFrame | None = None, **kw) -> None:
    """Write a standard table; keyword options of the builder (e.g. level="class") are part of the file name."""
    PSEUDOBULK.mkdir(parents=True, exist_ok=True)
    pb.to_parquet(PSEUDOBULK / f"{_key(name, kw)}.parquet", index=False)
    if donors is not None:
        donors.rename_axis("donor").to_parquet(PSEUDOBULK / f"{_key(name, kw)}_donors.parquet")


def load(name: str, subsample: int | None = None, seed: int = 0, **kw) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """(pseudobulk, donors) of a dataset; `subsample`: keep that many random donors (checks); other keywords go to
    the builder (e.g. level="class")."""
    f, fd = PSEUDOBULK / f"{_key(name, kw)}.parquet", PSEUDOBULK / f"{_key(name, kw)}_donors.parquet"
    if f.exists():
        pb = pd.read_parquet(f)
        donors = pd.read_parquet(fd) if fd.exists() else None
    else:
        mod = _builder(name)
        if not hasattr(mod, "from_legacy"):
            raise FileNotFoundError(f"{f} missing: run paper/common/datasets/{name}/build.py")
        pb, donors = mod.from_legacy(**kw)
    if donors is not None:
        donors.index = donors.index.astype(str)
    pb["donor"] = pb["donor"].astype(str)
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
        donors = donors.loc[donors.index.isin(keep)] if donors is not None else None
    return pb, donors
