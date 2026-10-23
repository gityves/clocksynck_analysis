"""A33, glucocorticoid-adjusted coupling done as in the legacy: target-specific adjustment on [1, GC of the target, GC of the
host tissues].

Pre-registered (post hoc-informed) in docs/PLAN_internal_desync.md (archive) (A33, "with GC" readout). Same design, estimand and
prediction as `compute_ventilator_tissues.py` (loss = coupling(fast) - coupling(ventilated); D = mean loss of food-entrained
minus vascular tissues; tissue-label permutation over all splits). What changes is the adjustment. The legacy residualised the
target's clock-plane state and the host on [1, GC score of the target, mean GC score of the donor's host tissues] (Frisch-Waugh
after the host fit). Here, for every target T, the pseudobulk gets a donor column `gcH` (mean GC score of the donor's reference
tissues other than T, the rows of T excluded) and `Dataset.adjust(["gc", "gcH"])` residualises the states of all compartments
(target and anchors alike, per compartment) on both before the host and the coupling are fitted. The target part equals the
legacy; the host is built from adjusted anchors instead of projected afterwards (a deviation, see the README). Run on the
clocksync `Dataset.adjust` (not `fit(covariates=)`, which cannot take a target-specific covariate).
Estimand/null: as A33. Interval: donor bootstrap (`Dataset.resample` of the adjusted data, default 100 replicates here; 1000 in
the main script). Outputs: a33_gc_couplings.csv, a33_gc_summary.csv. Legacy: ventilator_tissues/scripts/a33_ventilator_tissues.py.

    python compute_ventilator_gc.py [--subsample N] [--n-boot 100] [--jobs 4]
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402
import compute_ventilator_tissues as V  # noqa: E402

FOLDER = "cues_lost"


def target_dataset(sub: pd.DataFrame, T: str, present) -> cs.Dataset:
    """Dataset of one cohort whose states are adjusted on the target's [gc, gcH] (gcH: mean GC of the donor's other reference
    tissues, T left out)."""
    host_t = [t for t in S.V2_TISSUES if t in present and t != T]
    per = sub[sub.compartment.isin(host_t)].groupby("donor").gc.mean()
    tab = sub.assign(gcH=sub.donor.map(per)).reset_index(drop=True)
    ds = cs.Dataset(tab[["donor", "compartment", "depth", *S.CLOCK]], anchors=[t for t in S.V2_TISSUES if t in present],
                    min_donors=2)
    ds.table = ds.table.assign(gc=tab.gc.to_numpy(), gcH=tab.gcH.to_numpy())
    return ds.adjust(["gc", "gcH"])


def coupling(ds: cs.Dataset, T: str) -> dict:
    t = ds.fit(n_boot=0, gate=False, targets=[T]).table.loc[T]
    return dict(c_gc=t.coupling, psi_gc=t.psi, V_gc=t.V, n=t.n, boundary_gc=t.boundary)


def score(dss: dict, targets) -> pd.DataFrame:
    return pd.DataFrame([dict(tissue=T, cohort=c, **coupling(dss[c][T], T)) for c in dss for T in targets])


def loss(tab: pd.DataFrame) -> pd.Series:
    w = tab.pivot(index="tissue", columns="cohort", values="c_gc")
    return w["fast"] - w["vent"]


def boot_one(b, dss, targets):
    limit_threads(1)
    try:
        rep = {c: {T: d.resample(np.random.default_rng([0, b])) for T, d in dd.items()} for c, dd in dss.items()}
        return loss(score(rep, targets))
    except (ValueError, np.linalg.LinAlgError):
        return None


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = args.n_boot if args.n_boot is not None else (10 if run["quick"] else 100)
    min_donors = 5 if run["quick"] else 20
    tissues = V.FOOD + V.VASC
    pb, donors = V.load(tissues, args.subsample, args.seed)
    present = set(pb.compartment)
    targets = [t for t in tissues if t in present]
    for t in list(targets):
        n = [pb[(pb.compartment == t) & pb.donor.isin(donors.index[donors.hardy.isin(h)])].donor.nunique()
             for h in V.COHORTS.values()]
        if min(n) < min_donors:
            targets.remove(t)
    food, vasc = [t for t in V.FOOD if t in targets], [t for t in V.VASC if t in targets]
    dss = {}
    for c, hardy in V.COHORTS.items():
        sub = pb[pb.donor.isin(donors.index[donors.hardy.isin(hardy)])].reset_index(drop=True)
        dss[c] = {T: target_dataset(sub, T, present) for T in targets}
    tab = score(dss, targets)
    L = loss(tab)
    obs = V.gdiff(L, food, vasc)
    perm = np.array([L.loc[list(F)].mean() - L.loc[[t for t in targets if t not in F]].mean()
                     for F in itertools.combinations(targets, len(food))])
    reps = [r for r in Parallel(n_jobs=args.jobs)(delayed(boot_one)(b, dss, targets) for b in range(n_boot)) if r is not None]
    Db = np.array([V.gdiff(r, food, vasc) for r in reps])
    Db = Db[np.isfinite(Db)]
    tab.assign(loss_gc=tab.tissue.map(L)).to_csv(out / "a33_gc_couplings.csv", index=False)
    summ = pd.DataFrame([dict(readout="loss_gc_target_specific", D=obs, lo=np.quantile(Db, 0.025) if len(Db) else np.nan,
                              hi=np.quantile(Db, 0.975) if len(Db) else np.nan,
                              p_boot=(1 + np.sum(Db <= 0)) / (1 + len(Db)), n_boot_ok=len(Db),
                              p_perm=np.mean(perm >= obs - 1e-4), n_perm=len(perm), n_food=len(food), n_vasc=len(vasc))])
    summ.to_csv(out / "a33_gc_summary.csv", index=False)
    print(tab.assign(loss=tab.tissue.map(L)).round(3).to_string(index=False), "\n", summ.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
