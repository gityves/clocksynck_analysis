"""Rebuild every A17 dataset table under the current class rules, compare it with the table currently loaded (stored
parquet or legacy conversion) and save it with the rule stamp (`data.save`); roadmap 2.1 (2026-10-28).

Datasets: the A17 registry (`datasets/_atlas/a17_targets.csv`) plus the scored variants of the class-fix rebuild.
A dataset whose builder fails is reported and left as it is. Per compartment: rows, donors, cells and depth before and
after. The previous files go to `$SCCC_DATA/pseudobulk/pre_stamp/` before anything is overwritten.

    python rebuild_and_stamp.py [--jobs 4] [names ...]
Output: paper/results/logs/rebuild_and_stamp.md and rebuild_and_stamp.csv (changed compartments only).
"""

from __future__ import annotations

import shutil
import sys
import traceback
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PAPER))
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from common import cli, data  # noqa: E402
from common.paths import DATASETS, PSEUDOBULK  # noqa: E402

SCORED = ["ia_crohn_scored", "heart_dcm_scored", "heart_myocarditis_scored", "kong_cd_scored", "periodontitis_scored"]
LOG = PAPER / "results" / "logs"


def summary(pb: pd.DataFrame) -> pd.DataFrame:
    agg = dict(rows=("donor", "size"), donors=("donor", "nunique"), depth=("depth", "sum"))
    if "n_cells" in pb:
        agg["cells"] = ("n_cells", "sum")
    return pb.groupby("compartment").agg(**agg)


def one(name: str) -> dict:
    from clocksync.config import limit_threads
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        prev = PSEUDOBULK / "pre_stamp" / f"{name}.parquet"  # a rerun compares with the table before the first stamp
        old = pd.read_parquet(prev) if prev.exists() else data.load(name, allow_stale=True)[0]
        pb, donors = data._builder(name).build()
        a, b = summary(old), summary(pb)
        d = a.join(b, how="outer", lsuffix="_old", rsuffix="_new").fillna(0)
        cols = [c for c in ("rows", "donors", "cells", "depth") if f"{c}_old" in d and f"{c}_new" in d]
        changed = d[(d[[f"{c}_old" for c in cols]].to_numpy() != d[[f"{c}_new" for c in cols]].to_numpy()).any(1)]
        bak = PSEUDOBULK / "pre_stamp"
        bak.mkdir(exist_ok=True)
        for f in (PSEUDOBULK / f"{name}.parquet", PSEUDOBULK / f"{name}_donors.parquet"):
            if f.exists() and not (bak / f.name).exists():
                shutil.copy2(f, bak / f.name)
        data.save(name, pb, donors)
        return dict(name=name, status="ok", hash=data.rules_hash(name), n_changed=len(changed),
                    changed=changed.reset_index().assign(dataset=name))
    except Exception as e:  # noqa: BLE001  a failed builder is a row, never silently skipped
        return dict(name=name, status=f"error: {e!r}"[:300], trace=traceback.format_exc()[-1500:], n_changed=None)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("names", nargs="*")
    args = ap.parse_args(argv)
    cli.setup(args)
    names = args.names or sorted(set(pd.read_csv(DATASETS / "_atlas" / "a17_targets.csv").dataset) | set(SCORED))
    print(f"{len(names)} datasets", flush=True)
    res = Parallel(n_jobs=args.jobs, verbose=5)(delayed(one)(n) for n in names)
    ch = [r.pop("changed") for r in res if "changed" in r]
    C = pd.concat([c for c in ch if len(c)], ignore_index=True) if any(len(c) for c in ch) else pd.DataFrame()
    C.to_csv(LOG / "rebuild_and_stamp.csv", index=False)
    S = pd.DataFrame([{k: v for k, v in r.items() if k != "trace"} for r in res])
    with open(LOG / "rebuild_and_stamp.md", "w") as fh:
        fh.write("# Rebuild and stamp (roadmap 2.1)\n\n```\n" + S.to_string(index=False) + "\n```\n\n")
        fh.write("## Changed compartments\n\n" + ("```\n" + C.to_string(index=False) + "\n```" if len(C) else "none") + "\n\n")
        for r in res:
            if "trace" in r:
                fh.write(f"## {r['name']} failed\n\n```\n{r['trace']}\n```\n\n")
    print(S.to_string(), flush=True)
    print(C.to_string() if len(C) else "no compartment changed", flush=True)


if __name__ == "__main__":
    main()
