"""A17 (steps 0, 3, 8) on the healthy atlases: gate, anchored couplings of every compartment with the frozen learned and
fixed anchor sets, host-identity diagnostics, and the cohort table. One row per healthy fit x rule x compartment.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17; origin A1-A16, A20). Per healthy dataset (34 fits): the model-free
gate (1,000 permutations; nothing is estimated when it fails); anchors from `common/anchors_a17.csv` (learned leave one
study out, or the fixed stromal-vascular / neuronal rule), each anchor on the host of the others; coupling c_k, relative
amplitude and offset with donor-bootstrap intervals (host refitted, default 500). Host identity (agreement, stability) is
in host_identity. Descriptive. Exploratory-labelled cohorts (Tabula Sapiens: confirmation only) are fitted as in the legacy run.

    python compute_atlas_fits.py [--subsample N] [--jobs 4] [--only NAME ...] [--n-boot B]
Outputs: paper/results/figure_2/output/atlas_fits/a17_healthy_fits.csv (legacy a17_healthy.csv), a17_healthy_gates.csv (legacy
a17_gates.csv), cohort_table.csv, a17_healthy_hosts_summary.csv (host spread per fit; per-donor states are not written).
Legacy: legacy/desync/common/run_anchored.py (step healthy), anchors.py (gates).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from _shared import FOLDER_5, Target, frozen, registry  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

N_PERM_GATE, N_BOOT = 1000, 500


def anchor_sets(name: str) -> dict:
    """Frozen A17 anchors: {rule: list or None} (None: fewer than 3, or "no anchor" for the learned rule)."""
    a = frozen()
    a = a[a.target == name]
    out = {}
    for rule, col in (("learned", "anchor"), ("fixed", "anchor_fixed")):
        comps = a.compartment[a[col].astype(bool)].tolist()
        ok = len(comps) >= 3 and not (rule == "learned" and str(a.status.iloc[0]).startswith("no anchor"))
        out[rule] = comps if ok else None
    return out


def fit_target(name: str, i: int, args, n_boot: int, n_perm: int) -> dict:
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        T = Target(name, args.subsample, args.seed)
        ds_free = T.dataset(anchors=None)
        g = ds_free.gate(n_perm=n_perm, seed=i)
        gate = dict(target=name, **{k: getattr(g, k) for k in ("clock_cov", "clock_share", "share_null_mean",
                                                            "p_pairing", "p_acrophase", "n_pairs", "n_donors",
                                                            "n_compartments")}, passed=bool(g.passes))
        rows, hosts = [], []
        for rule, anc in anchor_sets(name).items():
            head = dict(target=name, rule=rule, gate_p_pairing=g.p_pairing, gate_p_acrophase=g.p_acrophase,
                        gate_passed=bool(g.passes))
            anc = [a for a in anc or [] if a in set(T.pb.compartment)] or None
            if not g.passes or anc is None or len(anc) < 3:
                rows.append(pd.DataFrame([head | dict(status="gate failed" if not g.passes else "no anchor")]))
                continue
            ds = T.dataset(anchors=anc)
            fit = ds.fit(n_boot=n_boot, seed=args.seed, gate=g)
            tab = fit.table.reset_index().rename(columns={"rel_amplitude": "amplitude"})
            tab = tab.assign(**head, status="ok")
            rows.append(tab)
            hosts.append(dict(target=name, rule=rule, n_donors=len(fit.host), host_spread=float(fit.host[["x", "y"]]
                              .var().sum()), median_radius=float(fit.host.radius.median()), n_anchors=len(anc)))
        s = T.pb.groupby("compartment").agg(n=("donor", "nunique"), cells=("n_cells", "median"),
                                            depth=("depth", "median"))
        r = registry().loc[name]
        cohort = dict(target=name, dataset=r.dataset, kind=r.kind, study=r.study, tissue=r.tissue,
                      donors=T.pb.donor.nunique(), compartments=len(s), min_donors_comp=int(s.n.min()),
                      median_donors_comp=float(s.n.median()), median_cells=float(s.cells.median()),
                      median_depth=float(s.depth.median()))
        return dict(rows=pd.concat(rows, ignore_index=True), gate=gate, cohort=cohort, hosts=hosts)
    except Exception as e:  # reported; the other targets continue
        return dict(error=f"{name}: {type(e).__name__}: {e}")


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--n-perm-gate", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_boot = args.n_boot if args.n_boot is not None else (20 if run["quick"] else N_BOOT)
    n_perm = args.n_perm_gate if args.n_perm_gate is not None else (100 if run["quick"] else N_PERM_GATE)
    out = out_dir(FOLDER_5, run["quick"])
    reg = registry()
    names = [n for n in reg.name[reg.kind == "healthy"] if args.only is None or n in args.only]
    res = Parallel(n_jobs=args.jobs)(delayed(fit_target)(n, i, args, n_boot, n_perm)
                                     for i, n in enumerate(names))
    for r in res:
        if "error" in r:
            print("FAILED", r["error"], flush=True)
    res = [r for r in res if "error" not in r]
    fits = pd.concat([r["rows"] for r in res], ignore_index=True)
    fits.to_csv(out / "a17_healthy_fits.csv", index=False)
    pd.DataFrame([r["gate"] for r in res]).to_csv(out / "a17_healthy_gates.csv", index=False)
    pd.DataFrame([r["cohort"] for r in res]).to_csv(out / "cohort_table.csv", index=False)
    pd.DataFrame([h for r in res for h in r["hosts"]]).to_csv(out / "a17_healthy_hosts_summary.csv", index=False)
    print(fits.groupby(["rule", "status"]).target.nunique(), flush=True)


if __name__ == "__main__":
    main()
