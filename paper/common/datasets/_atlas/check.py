"""Check a builder against the legacy table on a slice of donors (no full rebuild, nothing saved).

    python check.py <dataset> [<dataset> ...] [--n 8] [--seed 0] [--twins]

For each dataset `paper/common/datasets/<dataset>/build.py`: `from_legacy()` is run on the full legacy table (rows, donors,
compartments, strata); `build(donors=[...])` is run on a slice of donors (random donors plus one case and one control
of every contrast column) and compared with the legacy rows of the same donors: rows present in only one table,
largest absolute difference of n_cells, depth and counts, strata, and the donor-table columns the two share.
`--twins`: the genome-wide datasets (`*_genome`, no legacy table) against the clock-gene dataset that has the same rows, and the
T-state datasets against `hlca_fibrotic` / `lung_natri_control`.
Rows only in the rebuilt slice must belong to compartments the legacy table dropped by its `min_donors` rule (that
rule needs the whole dataset and is not applied to a slice).
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from census import GENES  # noqa: E402


def load_builder(name: str):
    f = HERE.parent / name / "build.py"
    spec = importlib.util.spec_from_file_location(f"build_{name}", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pick(pb: pd.DataFrame, donors: pd.DataFrame, n: int, seed: int) -> list:
    rng = np.random.default_rng(seed)
    ids = list(rng.permutation(sorted(pb.donor.unique()))[:n])
    for c in donors.columns:
        if donors[c].dropna().isin(["case", "control"]).all() and donors[c].notna().any():
            for v in ("case", "control"):
                s = donors.index[donors[c] == v]
                if len(s):
                    ids.append(rng.choice(s))
    return sorted(set(ids))


def _maxdiff(x: pd.DataFrame, y: pd.DataFrame) -> float:
    """Largest absolute difference; NaN cells must be NaN in both (else inf)."""
    x, y = x.to_numpy(float), y.to_numpy(float)
    if (np.isnan(x) != np.isnan(y)).any():
        return float("inf")
    return float(np.nanmax(np.abs(x - y))) if x.size else 0.0


def compare(name: str, n: int = 8, seed: int = 0, rebuild: bool = True, **kw) -> dict:
    mod = load_builder(name)
    pb_l, don_l = mod.from_legacy()
    g = [c for c in GENES]
    res = dict(dataset=name, rows=len(pb_l), donors=pb_l.donor.nunique(), compartments=pb_l.compartment.nunique(),
               strata=pb_l.stratum.nunique() if "stratum" in pb_l else 0,
               genes_ok=all(c in pb_l for c in g), n_gene_columns=pb_l.shape[1] - 5, dup=int(pb_l.duplicated(["donor", "compartment"]).sum()),
               donor_cols=list(don_l.columns) if don_l is not None else [])
    hook = getattr(mod, "check_legacy", None) or getattr(getattr(mod, "M", None), "check_legacy", None)
    if hook:
        res.update(hook(getattr(mod, "DATASET", name), pb_l, don_l))
    if not rebuild:
        return res
    ids = pick(pb_l, don_l, n, seed)
    raw = [i.split(":", 1)[1] if ":" in i else i for i in ids]
    pb_b, don_b = mod.build(donors=raw, **kw)
    a = pb_l[pb_l.donor.isin(ids)].set_index(["donor", "compartment"]).sort_index()
    b = pb_b[pb_b.donor.isin(ids)].set_index(["donor", "compartment"]).sort_index()
    common = a.index.intersection(b.index)
    num = ["n_cells", "depth"] + [c for c in a.columns if c in b.columns and c not in ("stratum", "n_cells", "depth")]  # every shared count column
    only_l, only_b = a.index.difference(b.index), b.index.difference(a.index)
    extra_comp = set(b.loc[only_b].index.get_level_values(1))
    explained = extra_comp.isdisjoint(set(pb_l.compartment))  # extra rows only in compartments the legacy dropped
    res.update(slice_donors=len(ids), slice_rows_legacy=len(a), slice_rows_build=len(b), only_legacy=len(only_l),
               only_build=len(only_b), only_build_explained=explained,
               max_abs_diff=_maxdiff(a.loc[common, num], b.loc[common, num]) if len(common) else np.nan,
               strata_equal=bool((a.loc[common, "stratum"].astype(str) == b.loc[common, "stratum"].astype(str)).all())
               if "stratum" in a and len(common) else None)
    if don_l is not None and don_b is not None:
        cols = [c for c in don_b.columns if c in don_l.columns]
        bad = []
        for c in cols:
            x = don_l.loc[[i for i in ids if i in don_l.index], c]
            y = don_b.reindex(x.index)[c]
            same = (x.astype(str).to_numpy() == y.astype(str).to_numpy()) | (x.isna().to_numpy() & y.isna().to_numpy())
            if not same.all():
                if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
                    same = np.isclose(x.astype(float), y.astype(float), equal_nan=True)
            if not same.all():
                bad.append(c)
        miss = sorted(set(don_l.columns) - set(don_b.columns) - {"healthy"} - {c for c in don_l if c.startswith("block_")})
        res.update(donor_cols_compared=len(cols), donor_cols_differing=bad, donor_cols_missing_in_build=miss,
                   donor_cols_only_in_build=sorted(set(don_b.columns) - set(don_l.columns)))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("datasets", nargs="*")
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--legacy-only", action="store_true")
    ap.add_argument("--twins", action="store_true", help="also run the twin checks of the genome-wide and T-state datasets")
    args = ap.parse_args()
    if args.twins:
        import genome
        import tstates
        for name in genome.CLOCK_TWIN:
            print(genome.check_clock_twin(name, args.n, args.seed), flush=True)
        for r in tstates.check_twins(args.n, args.seed):
            print(r, flush=True)
        if not args.datasets:
            return
    for name in args.datasets:
        try:
            r = compare(name, args.n, args.seed, rebuild=not args.legacy_only)
        except Exception as e:  # reported, the others continue
            import traceback
            traceback.print_exc()
            r = dict(dataset=name, error=f"{type(e).__name__}: {e}")
        print(r, flush=True)


if __name__ == "__main__":
    main()
