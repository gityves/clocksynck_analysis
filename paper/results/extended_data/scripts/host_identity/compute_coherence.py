"""A17 steps 0-1: gate and model-free clock-plane coherence of every healthy fit and every contrast table.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17, step 1 as revised by the note of 2026-10-12). Per target (36 healthy
and contrast tables of the 70 A17 targets; the 21 cancer studies belong to cancer): the gate (1,000
permutations; healthy fits only: the learning set is the gate-passing healthy fits) and the coherence of every
compartment pair, rho_kl = |complex correlation| of the clock coefficients over the donors the pair shares (centred within
strata), its excess over the shuffled-pairing null (200 shuffles), and per compartment the mean excess over its partners
with a donor-bootstrap SE (200; healthy fits only). Contrast tables: pair coherence on the pooled donors, label-free
(no bootstrap), for the confirmation step of the anchor rule only.

    python compute_coherence.py [--subsample N] [--jobs 4] [--only NAME ...]
Outputs: paper/results/extended_data/output/host_identity/a17_gates.csv, a17_coherence_classes.csv, a17_coherence_pairs.csv
(column names of the legacy files, `name` = target).
Legacy: legacy/desync/common/anchors.py (step coherence).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("atlas_fits")))
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from _shared import FOLDER_S04, Target, registry  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402


def one(name: str, i: int, subsample, seed: int, n_perm_gate: int, n_perm: int, n_boot: int) -> dict:
    """Seed = position of the target in the registry (legacy: the seed of `anchors.one`)."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        T = Target(name, subsample, seed, new_classes=True)
        ds = T.dataset(anchors=None)
        out = dict(name=name)
        if T.kind == "healthy":
            g = ds.gate(n_perm=n_perm_gate, seed=i)
            out["gate"] = dict(name=name, clock_cov=g.clock_cov, p_pairing=g.p_pairing, clock_share=g.clock_share,
                               share_null_mean=g.share_null_mean, p_acrophase=g.p_acrophase, n_pairs=g.n_pairs,
                               n_resections=g.n_donors, n_compartments=g.n_compartments, passed=bool(g.passes))
            pairs, classes = ds.coherence(n_perm=n_perm, n_boot=n_boot, seed=i)
            out["classes"] = classes.assign(name=name)
        else:
            pairs, _ = ds.coherence(n_perm=n_perm, n_boot=0, seed=i)
        out["pairs"] = pairs.assign(name=name)
        return out
    except Exception as e:  # reported; the other targets continue
        return dict(error=f"{name}: {type(e).__name__}: {e}")


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER_S04, run["quick"])
    reg = registry()
    todo = [(n, i) for i, n in enumerate(reg.name) if args.only is None or n in args.only]
    n_pg, n_p, n_b = (100, 50, 20) if run["quick"] else (1000, 200, 200)
    res = Parallel(n_jobs=args.jobs)(delayed(one)(n, i, args.subsample, args.seed, n_pg, n_p, n_b) for n, i in todo)
    for r in res:
        if "error" in r:
            print("FAILED", r["error"], flush=True)
    res = [r for r in res if "error" not in r]
    pd.DataFrame([r["gate"] for r in res if "gate" in r]).to_csv(out / "a17_gates.csv", index=False)
    pd.concat([r["classes"] for r in res if "classes" in r]).to_csv(out / "a17_coherence_classes.csv", index=False)
    pd.concat([r["pairs"] for r in res]).to_csv(out / "a17_coherence_pairs.csv", index=False)
    print(f"{len(res)} targets", flush=True)


if __name__ == "__main__":
    main()
