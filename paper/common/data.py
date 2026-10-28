"""Standard pseudobulk tables, one per dataset: `PSEUDOBULK/<name>.parquet` (one row per donor x compartment:
donor, compartment, [stratum], [n_cells], depth, clock-gene counts) and `PSEUDOBULK/<name>_donors.parquet` (one
row per donor: groups, covariates, sampling times; controlled-access columns never leave the machine).

Each dataset has a builder `paper/common/datasets/<name>/build.py` with `build() -> (pb, donors)` from the public raw
files. Until the data are rebuilt, a builder may also provide `from_legacy() -> (pb, donors)`, which converts the
table the legacy code built; `load` uses it when the standard table is missing (with a warning: the legacy class
rules were never audited against the current ones).

Rule version. `save` stamps each pseudobulk parquet with `RULES_KEY` = `rules_hash(name)`, the fingerprint of the class
rules the builder used; `load` raises `StaleTableError` when the stored stamp differs from the current one (the rules
changed after the table was built: rebuild it, or pass `allow_stale=True` to read it anyway). A table without a stamp
(saved before the scheme existed) only warns. Scheme (`r1-<16 hex>`, sha256 over): (1) the module-level constants
(UPPER_CASE names holding str / number / tuple / list / dict / set values, e.g. COARSE, CLASSES, BRAIN, L1, REGIONS) of
the builder module and, recursively, of every module under `datasets/` that it imports; (2) the builder source,
normalised through `ast` (comments, formatting and docstrings do not count). Edits to the functions of shared modules
(`_atlas/*.py`) and to input CSVs are not covered: change a rule constant, or rebuild, when they alter the numbers.
"""

from __future__ import annotations

import ast
import functools
import hashlib
import importlib.util
import re
import types
import warnings

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from common.paths import DATASETS, PSEUDOBULK  # noqa: E402

RULES_KEY = b"sccc_rules_hash"
_SCHEME = "r1"


class StaleTableError(RuntimeError):
    """The stored table was built under class rules that differ from the current builder's."""


def _builder(name: str):
    f = DATASETS / name / "build.py"
    if not f.exists():
        raise FileNotFoundError(f"no builder for dataset {name!r} ({f})")
    spec = importlib.util.spec_from_file_location(f"build_{name}", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _plain(v):
    """Order-stable plain form of a constant (None when it is not a rule-like value)."""
    if isinstance(v, (str, bytes, int, float, bool)) or v is None:
        return v
    if isinstance(v, (tuple, list)):
        out = [_plain(x) for x in v]
        return None if any(x is None and y is not None for x, y in zip(out, v)) else out
    if isinstance(v, (set, frozenset)):
        return sorted(map(repr, v))
    if isinstance(v, dict):
        return sorted((repr(k), repr(_plain(x))) for k, x in v.items())
    return None


def _constants(mod, seen: set) -> list:
    """Rule constants of a module and of the `datasets/` modules it imports, as sorted (module, name, value)."""
    out = []
    root = str(DATASETS)
    for k, v in sorted(vars(mod).items()):
        if isinstance(v, types.ModuleType) and str(getattr(v, "__file__", "")).startswith(root) and v.__name__ not in seen:
            seen.add(v.__name__)
            out += _constants(v, seen)
        elif re.fullmatch(r"[A-Z][A-Z0-9_]*", k) and (p := _plain(v)) is not None:
            out.append((mod.__name__.removeprefix("build_"), k, repr(p)))
    return out


def _normalised_source(path) -> str:
    tree = ast.parse(path.read_text())
    for n in ast.walk(tree):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and n.body \
                and isinstance(n.body[0], ast.Expr) and isinstance(getattr(n.body[0], "value", None), ast.Constant) \
                and isinstance(n.body[0].value.value, str):
            n.body = n.body[1:] or [ast.Pass()]
    return ast.unparse(tree)


@functools.lru_cache(maxsize=None)
def rules_hash(name: str) -> str:
    """Fingerprint of the class rules the builder of `name` applies now (scheme in the module docstring)."""
    mod = _builder(name)
    h = hashlib.sha256()
    h.update(repr(_constants(mod, {mod.__name__})).encode())
    h.update(_normalised_source(DATASETS / name / "build.py").encode())
    return f"{_SCHEME}-{h.hexdigest()[:16]}"


def _stored_hash(f) -> str | None:
    md = pq.read_schema(f).metadata or {}
    v = md.get(RULES_KEY)
    return v.decode() if v else None


def _key(name: str, kw: dict) -> str:
    return name + "".join(f"__{k}-{v}" for k, v in sorted(kw.items()))


def save(name: str, pb: pd.DataFrame, donors: pd.DataFrame | None = None, **kw) -> None:
    """Write a standard table; keyword options of the builder (e.g. level="class") are part of the file name."""
    PSEUDOBULK.mkdir(parents=True, exist_ok=True)
    stamp = {RULES_KEY: rules_hash(name).encode()} if (DATASETS / name / "build.py").exists() else {}

    def write(df, f, index):
        t = pa.Table.from_pandas(df, preserve_index=index)
        pq.write_table(t.replace_schema_metadata({**(t.schema.metadata or {}), **stamp}), f)

    write(pb, PSEUDOBULK / f"{_key(name, kw)}.parquet", False)
    if donors is not None:
        write(donors.rename_axis("donor"), PSEUDOBULK / f"{_key(name, kw)}_donors.parquet", True)


def _check_rules(name: str, f, allow_stale: bool) -> None:
    if not (DATASETS / name / "build.py").exists():
        return
    stored, now = _stored_hash(f), rules_hash(name)
    if stored is None:
        warnings.warn(f"{f.name}: no rule stamp (saved before rule versioning); current rules are {now}", UserWarning, stacklevel=3)
    elif stored != now:
        msg = f"{f.name}: built under class rules {stored}, the builder of {name!r} now has {now}; rebuild the table"
        if not allow_stale:
            raise StaleTableError(msg + " (or load(..., allow_stale=True))")
        warnings.warn(msg, UserWarning, stacklevel=3)


def load(name: str, subsample: int | None = None, seed: int = 0, allow_stale: bool = False, **kw) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """(pseudobulk, donors) of a dataset; `subsample`: keep that many random donors (checks); other keywords go to
    the builder (e.g. level="class"). A stored table whose rule stamp differs from the current builder's raises
    `StaleTableError` unless `allow_stale` (then it warns); a legacy fallback and an unstamped table warn."""
    f, fd = PSEUDOBULK / f"{_key(name, kw)}.parquet", PSEUDOBULK / f"{_key(name, kw)}_donors.parquet"
    if f.exists():
        _check_rules(name, f, allow_stale)
        pb = pd.read_parquet(f)
        donors = pd.read_parquet(fd) if fd.exists() else None
    else:
        mod = _builder(name)
        if not hasattr(mod, "from_legacy"):
            raise FileNotFoundError(f"{f} missing: run paper/common/datasets/{name}/build.py")
        warnings.warn(f"{name}: no standard table, loading the legacy conversion (from_legacy); its class rules are not "
                      f"versioned or audited. Build and save it: paper/common/datasets/{name}/build.py", UserWarning, stacklevel=2)
        pb, donors = mod.from_legacy(**kw)
    if donors is not None:
        donors.index = donors.index.astype(str)
    pb["donor"] = pb["donor"].astype(str)
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
        donors = donors.loc[donors.index.isin(keep)] if donors is not None else None
    return pb, donors
