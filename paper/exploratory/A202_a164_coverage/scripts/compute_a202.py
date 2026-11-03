"""A202: is the A164 equivalence interval valid at its own design? Known coupling changes injected into the real GTEx
counts by binomial thinning (`clocksync.power.damped`), labels exchangeable within the A164 matched sets.

Registered in paper/docs/registered_experiments.md (A202). Choices not in the registration: ../README.md.
Arm P: point estimates (cells 0, -0.075, -0.15, +0.15; 1,000 replicates). Arm B: the A164 matched-set bootstrap per
replicate (cells 0, -0.15; 100 replicates, 200 bootstrap replicates). Controlled access: only aggregates are written
(per replicate: module estimates and truths per cohort); matched sets and labels stay in memory.

    python compute_a202.py --arm P [--jobs 16] [--reps 1000] [--subsample N]
    python compute_a202.py --arm B [--jobs 16] [--reps 100] [--n-boot 200]
Outputs: ../output/raw/<arm>/<design>_<cell>_<chunk>.csv (resumable), ../output/design.csv (per cohort x design x
tissue: donors, pooled c0, damping s and realised change per cell).
"""

from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
A164 = PAPER / "exploratory/A164_gtex_metabolic_matched"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_3/scripts/gtex_cues"))
sys.path.insert(0, str(A164 / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import _covariates as K  # noqa: E402
import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
import compute_a164 as C  # noqa: E402
from clocksync import power  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

DESIGNS = ("t2d", "obese")
CELLS = {"P": (0.0, -0.075, -0.15, 0.15), "B": (0.0, -0.15)}
REPS = {"P": 1000, "B": 100}
MIN_GAP = 0.02  # a tissue with c0 < |delta| + MIN_GAP is damped to s = 0 (realised change c0)
OUT = HERE.parent / "output"
say = C.say


def damping(c0: float, delta: float) -> tuple[float, float]:
    """(s, realised |change|) with c0 - c'(s) = |delta|, c'(s) = s^2 c0 / (s^2 c0 + 1 - c0) (`power.mde`)."""
    if delta == 0 or not np.isfinite(c0) or c0 <= 0:
        return 1.0, 0.0
    if c0 < abs(delta) + MIN_GAP:
        return 0.0, c0
    c1 = c0 - abs(delta)
    return float(np.sqrt(c1 * (1 - c0) / (c0 * (1 - c1)))), abs(delta)


class Cohort:
    """One GTEx cohort at one A164 design: the Dataset, matched sets per module tissue, each tissue's host on the
    labelled donors (unchanged by relabelling, so fitted once) and its pooled coupling c0."""

    def __init__(self, name: str, pb: pd.DataFrame, plan: dict):
        self.name, self.pb = name, pb
        self.anchors = [a for a in S.TISSUES if a in set(pb.compartment)]
        self.ds = cs.Dataset(pb, anchors=self.anchors)
        self.plan = {t: s for t, s in plan.items() if t in C.MET}
        self.donors = sorted({d for s in self.plan.values() for x in s for d in x})
        self.cases = sorted({s[0] for v in self.plan.values() for s in v})
        d = A.Data.of(self.ds.study)
        self.host, self.c0 = {}, {}
        for t, sets in self.plan.items():
            ids = [x for s in sets for x in s]
            dt = d.take(np.flatnonzero(np.isin(d.resection, ids)))
            self.host[t] = A._hosts(dt, self.anchors, [t])[t]
            self.c0[t] = float(A.on_host(dt, t, self.host[t], ids)["coupling"])

    def thinned(self, t: str, lab: pd.Series, delta: float, rng) -> tuple[pd.DataFrame, float]:
        """Long table with tissue t damped in the case (delta < 0) or control (delta > 0) donors, and the true change
        of d* = c_case - c_control; (None, NaN) when the arm's regression is rank-deficient. Only t is thinned, so no host of t contains a damped compartment."""
        s, ch = damping(self.c0[t], delta)
        if ch == 0:
            return self.pb, 0.0
        arm = "case" if delta < 0 else "control"
        tab = power.damped(self.ds.study, lab, self.anchors, t, s, rng, arm, h=self.host[t])
        if tab is None:  # rank-deficient arm (`power.damped`): the tissue is left out of this replicate
            return None, np.nan
        return tab.rename(columns={"resection": "donor"}), float(np.sign(delta) * ch)


def keyed(sets, u):
    """Sets reordered so that the permuted case (largest key) comes first, keeping the set's original case as its id."""
    out = []
    for s in sets:
        k = int(np.argmax([u[d] for d in s]))
        out.append((s[0], [s[k]] + [x for i, x in enumerate(s) if i != k]))
    return out


def boot_se(co: Cohort, tabs: dict, ksets: dict, n_boot: int, rng) -> float:
    """A164 matched-set bootstrap of the module mean: original cases resampled with replacement jointly over tissues,
    each draw brings its set (with the permuted labels) in every tissue; donors copied, Dataset and host refitted."""
    by_case = {t: {cid: s for cid, s in ksets[t]} for t in ksets}
    vals = []
    for _ in range(n_boot):
        draw = rng.choice(len(co.cases), len(co.cases))
        row = []
        for t in ksets:
            sets = [(by_case[t][co.cases[i]], j) for j, i in enumerate(draw) if co.cases[i] in by_case[t]]
            mp = pd.DataFrame(sorted({(d, f"{d}@{j}") for s, j in sets for d in s}), columns=["orig", "copy"])
            rep = tabs[t].merge(mp, left_on="donor", right_on="orig")
            rep["donor"] = rep["copy"]
            rep = rep.drop(columns=["orig", "copy"]).reset_index(drop=True)
            ds = cs.Dataset(rep, anchors=[a for a in co.anchors if a in set(rep.compartment)])
            lab = pd.concat([C.labels([s], copy=j) for s, j in sets])
            row.append(C.dstar(ds, lab, t))
        row = np.asarray(row, float)
        vals.append(np.nanmean(row) if np.isfinite(row).any() else np.nan)
    return float(np.nanstd(vals, ddof=1))


def replicate(co: Cohort, di: int, ci: int, delta: float, rep: int, n_boot: int, seed: int) -> dict:
    ss = np.random.SeedSequence([202, di, ci, rep, seed, ["vent", "fast"].index(co.name)])
    rng = np.random.default_rng(ss)
    u = dict(zip(co.donors, rng.random(len(co.donors))))
    est, tru, tabs, ksets = {}, {}, {}, {}
    for t, sets in co.plan.items():
        ksets[t] = keyed(sets, u)
        lab = C.labels([s for _, s in ksets[t]])
        tabs[t], tru[t] = co.thinned(t, lab, delta, rng)
        if tabs[t] is None:
            est[t] = np.nan
            continue
        ds = co.ds if tabs[t] is co.pb else cs.Dataset(tabs[t], anchors=co.anchors)
        est[t] = C.dstar(ds, lab, t)
    ok = [t for t in est if np.isfinite(est[t])]
    out = dict(cohort=co.name, rep=rep, n_tissues=len(ok), n_thin_failed=sum(v is None for v in tabs.values()),
               D=float(np.mean([est[t] for t in ok])) if ok else np.nan,
               truth=float(np.mean([tru[t] for t in ok])) if ok else np.nan)
    if n_boot:
        out["se_boot"] = boot_se(co, {t: tabs[t] for t in ok}, {t: ksets[t] for t in ok}, n_boot, rng)
    return out


def chunk(cohorts, di, ci, delta, reps, n_boot, seed):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    return [replicate(co, di, ci, delta, r, n_boot, seed) for r in reps for co in cohorts]


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--arm", choices=("P", "B"), required=True)
    ap.add_argument("--reps", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=200)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    warnings.filterwarnings("ignore")
    quick = run["quick"]
    reps = args.reps or (4 if quick else REPS[args.arm])
    n_boot = 0 if args.arm == "P" else (5 if quick else args.n_boot)
    out = OUT / "subsample" if quick else OUT
    raw = out / "raw" / args.arm
    raw.mkdir(parents=True, exist_ok=True)

    ph = K.protected()
    cohorts, drow = {}, []
    for c, h in C.COH.items():
        pb, don = K.cohort(h, subsample=args.subsample, seed=args.seed)
        pb = pb[pb.compartment.isin(set(S.TISSUES) | set(C.TISS))].reset_index(drop=True)
        q = don[["hardy", "male", "ischaemic_h"]].join(ph[["age", "bmi", "t2d", "steroid"]]).dropna(
            subset=["age", "ischaemic_h"])
        td = {t: set(pb.donor[pb.compartment == t]) for t in C.TISS}
        for dsg in DESIGNS:
            plan, _, _ = C.plan_of(q, td, dsg, 3 if quick else 10)
            if not any(t in C.MET for t in plan):
                continue
            co = Cohort(c, pb, plan)
            cohorts[dsg, c] = co
            for t in co.plan:
                for delta in CELLS["P"]:
                    s, ch = damping(co.c0[t], delta)
                    drow.append(dict(cohort=c, design=dsg, tissue=t, n_sets=len(co.plan[t]),
                                     n_donors=sum(len(x) for x in co.plan[t]), c0=co.c0[t], delta=delta, s=s,
                                     realised=np.sign(delta) * ch))
    D = pd.DataFrame(drow)
    D.to_csv(out / "design.csv", index=False)
    say("design written\n", D[D.delta == -0.15].round(3).to_string(index=False))

    t0 = time.time()
    for di, dsg in enumerate(DESIGNS):
        cos = [cohorts[k] for k in cohorts if k[0] == dsg]
        for ci, delta in enumerate(CELLS["P"]):
            if delta not in CELLS[args.arm]:
                continue
            size = max(1, int(np.ceil(reps / (args.jobs * 4))))
            blocks = [list(range(i, min(i + size, reps))) for i in range(0, reps, size)]
            todo = [(k, b) for k, b in enumerate(blocks) if not (raw / f"{dsg}_{ci}_{k:04d}.csv").exists()]
            gen = Parallel(n_jobs=args.jobs, return_as="generator")(
                delayed(chunk)(cos, di, ci, delta, b, n_boot, args.seed) for _, b in todo)
            for n, ((k, _), res) in enumerate(zip(todo, gen), 1):
                pd.DataFrame(res).assign(design=dsg, cell=ci, delta=delta).to_csv(raw / f"{dsg}_{ci}_{k:04d}.csv",
                                                                                 index=False)
                if n % max(1, len(todo) // 8) == 0 or n == len(todo):
                    say(f"arm {args.arm} {dsg} delta {delta:+.3f}: {n}/{len(todo)} chunks, {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
