"""A17 step 3 identity diagnostics: is the host the same state whatever the anchors, and is it stable over donors?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17, "Identity diagnostics"). Per healthy fit with a passing gate and
>= 3 anchors, for the learned and the fixed rule: the host is read as f1 + i f2 and similarity = |complex correlation|
of centred (within stratum) host states over common donors (invariant to the rotation and scale the model leaves free).
- agreement: free host (every compartment an anchor) vs anchored host;
- stability: share of donor-bootstrap refits (200) of the anchored host whose state has similarity > 0.7 with the
  full-data host (first copy of each distinct donor drawn); < 0.8 is labelled "unstable";
- learned vs fixed: similarity of the two anchored hosts (fits that have both).
Local helpers (not in clocksync): `_shared.similarity`, `stability`, `agreement_free` (built on `Dataset.fit(...).host`,
`Dataset.resample`).

    python compute_host_identity.py [--subsample N] [--jobs 4] [--only NAME ...] [--n-stab B]
Outputs: paper/results/extended_data/output/host_identity/a17_identity.csv (target, rule, n_anchors, agreement_free, stability,
median_similarity, unstable, similarity_learned_fixed). Gate: the gate of `compute_coherence.py` in the same folder.
Legacy: legacy/desync/common/run_anchored.py (step healthy: agreement_free, stability, unstable).
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

from _shared import (FOLDER_S04, Target, agreement_free, donor_strata, frozen, host_of, registry,  # noqa: E402
                     similarity, stability)  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402


def anchor_sets(name: str) -> dict:
    a = frozen()
    a = a[a.target == name]
    out = {}
    for rule, col in (("learned", "anchor"), ("fixed", "anchor_fixed")):
        comps = a.compartment[a[col].astype(bool)].tolist()
        ok = len(comps) >= 3 and not (rule == "learned" and str(a.status.iloc[0]).startswith("no anchor"))
        out[rule] = comps if ok else None
    return out


def one(name: str, passed: bool, args, n_stab: int) -> list:
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        T = Target(name, args.subsample, args.seed)
        rows, hosts = [], {}
        for rule, anc in anchor_sets(name).items():
            anc = [a for a in anc or [] if a in set(T.pb.compartment)]
            if not passed or len(anc) < 3:
                rows.append(dict(target=name, rule=rule, status="gate failed" if not passed else "no anchor"))
                continue
            ds = T.dataset(anchors=anc)
            h = host_of(ds)
            hosts[rule] = (h, donor_strata(ds))
            st = stability(ds, n_boot=n_stab, seed=args.seed)
            rows.append(dict(target=name, rule=rule, status="ok", n_anchors=len(anc),
                             agreement_free=agreement_free(ds, h), stability=st["stability"],
                             median_similarity=st["median_similarity"], n_refits=st["n"],
                             unstable=st["stability"] < 0.8))
        if len(hosts) == 2:
            s = similarity(hosts["learned"][0], hosts["fixed"][0], hosts["learned"][1])
            for r in rows:
                r["similarity_learned_fixed"] = s
        return rows
    except Exception as e:
        return [dict(target=name, status=f"error: {type(e).__name__}: {e}")]


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--n-stab", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_stab = args.n_stab if args.n_stab is not None else (10 if run["quick"] else 200)
    out = out_dir(FOLDER_S04, run["quick"])
    gates = pd.read_csv(out / "a17_gates.csv").set_index("name").passed
    reg = registry()
    names = [n for n in reg.name[reg.kind == "healthy"] if n in gates.index and (args.only is None or n in args.only)]
    res = Parallel(n_jobs=args.jobs)(delayed(one)(n, bool(gates[n]), args, n_stab) for n in names)
    t = pd.DataFrame([r for x in res for r in x])
    t.to_csv(out / "a17_identity.csv", index=False)
    print(t.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
