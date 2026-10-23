"""A33 (stretch, post hoc-informed): ventilated donors, food-entrained vs systemic tissues.

Pre-registered (post hoc-informed) in docs/PLAN_internal_desync.md (archive) (A33). Coupling of a tissue to the host fitted on the
clock genes of the other 14 reference tissues (the anchors; an anchor is read on its leave-one-out host), in ventilated
donors (Hardy 0) and in fast deaths (Hardy 1-2); loss = coupling(fast) - coupling(ventilated). Prediction: a larger
loss in food-entrained tissues (liver, adipose x2, muscle, colon, pancreas, stomach) than in vascular / stromal ones
(arteries x3, tibial nerve, skin x2, thyroid): D = mean loss(food) - mean loss(vascular) > 0, one-sided donor bootstrap.
Companion: tissue-label permutation of D over all splits of the tested tissues into groups of the same sizes.

Two versions: `loss0` (no covariate; exact legacy readout) and `loss_gc` (glucocorticoid score as a covariate). In the
legacy `loss_gc` the target's states and the host were residualised on [1, GC of the target, GC of the host tissues]
after the host fit (Frisch-Waugh on the internal `on_host`). `clocksync` has no covariate-adjusted coupling, so here the
clock-plane state of every tissue is residualised on its own GC score (within tissue, per cohort) before the host and
the couplings are refitted (local helper `adjust_gc`; a deviation, see the README).

    python compute_ventilator_tissues.py [--subsample N] [--jobs 4] [--n-boot 1000]
Outputs: paper/results/figure_2/output/organs/a33_{tissues,summary,couplings}.csv  (aggregates only)
Legacy: legacy/desync/validation/ventilator_tissues/scripts/a33_ventilator_tissues.py
"""

from __future__ import annotations

import copy
import dataclasses
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402

FOLDER = "cues_lost"
FOOD = ["Liver", "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Muscle - Skeletal", "Colon - Transverse",
        "Small Intestine - Terminal Ileum", "Pancreas", "Stomach"]
VASC = ["Artery - Tibial", "Artery - Aorta", "Artery - Coronary", "Nerve - Tibial", "Skin - Sun Exposed (Lower leg)",
        "Skin - Not Sun Exposed (Suprapubic)", "Thyroid"]
GC = ["FKBP5", "TSC22D3", "ZBTB16", "KLF9", "KLF15", "DDIT4"]
COHORTS = {"vent": (0,), "fast": (1, 2)}


def load(targets, subsample, seed):
    """pb (clock genes + GC counts, depth = library size) of Hardy 0-2 donors in the reference tissues and targets, the
    per-row GC score (mean z of log CPM within tissue over all three Hardy classes, A27), and the donor table."""
    pb, donors = gtex.load(tissues=sorted(set(S.V2_TISSUES) | set(targets)), hardy=(0, 1, 2),
                           genes=list(S.CLOCK) + GC, depth="library_size", require_tod=False)
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb = pb[pb.donor.isin(keep)]
    lg = np.log(pb[GC].to_numpy(float) / pb.depth.to_numpy(float)[:, None] * 1e6 + 0.5)
    z = pd.DataFrame(lg, columns=GC).groupby(pb.compartment.to_numpy()).transform(lambda x: (x - x.mean()) / x.std())
    pb = pb.assign(gc=z.mean(1).to_numpy()).drop(columns=GC).reset_index(drop=True)
    return pb, donors


def adjust_gc(ds: cs.Dataset) -> cs.Dataset:
    """local: not in clocksync. Dataset whose clock-plane states are residualised on the row's own GC score (linearly,
    per tissue, mean kept) before the host and the couplings are fitted."""
    d = ds.data
    gc = ds.table.gc.to_numpy(float)
    b = d.beta.copy()
    for k in np.unique(d.compartment):
        i = np.flatnonzero(d.compartment == k)
        M = np.column_stack([np.ones(len(i)), gc[i] - gc[i].mean()])
        B, *_ = np.linalg.lstsq(M, b[i], rcond=None)
        b[i] = b[i] - M[:, 1:] @ B[1:]
    new = copy.copy(ds)
    new.data = dataclasses.replace(d, beta=b)
    return new


def couplings(ds: cs.Dataset, targets) -> pd.DataFrame:
    """c0 (no covariate) and c_gc per target on the host of the other reference tissues."""
    out = {}
    for lab, d in (("0", ds), ("gc", adjust_gc(ds))):
        f = d.fit(n_boot=0, gate=False, targets=targets).table
        out[lab] = f[["coupling", "psi", "V", "n", "boundary"]].rename(
            columns={c: f"{c}{'0' if lab == '0' else '_gc'}" for c in ("coupling", "psi", "V", "boundary")})
    t = out["0"].join(out["gc"].drop(columns="n"))
    return t.rename(columns={"coupling0": "c0", "coupling_gc": "c_gc"})


def score(dss: dict, targets) -> pd.DataFrame:
    rows = []
    for c, ds in dss.items():
        t = couplings(ds, targets).reset_index().rename(columns={"compartment": "tissue"}).assign(cohort=c)
        rows.append(t)
    return pd.concat(rows, ignore_index=True)


def losses(tab: pd.DataFrame) -> pd.DataFrame:
    w = tab.pivot(index="tissue", columns="cohort", values=["c0", "c_gc"])
    return pd.DataFrame({"loss0": w[("c0", "fast")] - w[("c0", "vent")],
                         "loss_gc": w[("c_gc", "fast")] - w[("c_gc", "vent")]})


def gdiff(L: pd.Series, food, vasc) -> float:
    return L.loc[food].mean() - L.loc[vasc].mean()


def boot_one(b: int, dss: dict, targets):
    limit_threads(1)
    rng = np.random.default_rng([0, b])
    try:
        rep = {c: ds.resample(rng) for c, ds in dss.items()}
        tab = score(rep, targets)
    except (ValueError, np.linalg.LinAlgError):
        return None
    return tab.assign(replicate=b)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None, help="donor bootstrap replicates (default 1000; 20 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = args.n_boot if args.n_boot is not None else (20 if run["quick"] else 1000)
    min_donors = 5 if run["quick"] else 20
    tissues = FOOD + VASC
    pb, donors = load(tissues, args.subsample, args.seed)
    present = set(pb.compartment)
    skipped = [t for t in tissues if t not in present]
    print("skipped (absent):", skipped, flush=True)
    food, vasc = [t for t in FOOD if t in present], [t for t in VASC if t in present]
    targets = food + vasc
    dss = {}
    for c, hardy in COHORTS.items():
        keep = donors.index[donors.hardy.isin(hardy)]
        sub = pb[pb.donor.isin(keep)].reset_index(drop=True)
        ds = cs.Dataset(sub[["donor", "compartment", "depth", *S.CLOCK]], anchors=[t for t in S.V2_TISSUES if t in present],
                        min_donors=2)
        ds.table = ds.table.assign(gc=sub.gc.to_numpy())
        dss[c] = ds
    for t in list(targets):
        n = {c: len(set(ds.table.resection[ds.table.compartment == t])) for c, ds in dss.items()}
        if min(n.values()) < min_donors:
            print(f"skipped (< {min_donors} donors): {t} {n}", flush=True)
            targets.remove(t)
    food, vasc = [t for t in food if t in targets], [t for t in vasc if t in targets]

    tab = score(dss, targets)
    L = losses(tab)
    obs = {k: gdiff(L[k], food, vasc) for k in ("loss_gc", "loss0")}
    print(tab.round(3).to_string(index=False), "\n", L.round(3).to_string(), "\nD:", {k: round(v, 4) for k, v in obs.items()},
          flush=True)

    reps = Parallel(n_jobs=args.jobs)(delayed(boot_one)(b, dss, targets) for b in range(n_boot))
    n_fail = sum(r is None for r in reps)
    bt = pd.concat([r for r in reps if r is not None], ignore_index=True)
    BL = {b: losses(x) for b, x in bt.groupby("replicate")}
    w = tab.pivot(index="tissue", columns="cohort")
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.join(L)
    w["group"] = ["food" if t in food else "vascular" for t in w.index]
    for k in ("loss_gc", "loss0"):
        M = pd.DataFrame({b: x[k] for b, x in BL.items()})
        w[f"{k}_lo"], w[f"{k}_hi"] = M.quantile(0.025, axis=1), M.quantile(0.975, axis=1)
    w.loc[food + vasc].to_csv(out / "a33_tissues.csv")
    tab.to_csv(out / "a33_couplings.csv", index=False)

    summ = []
    for k in ("loss_gc", "loss0"):
        Db = np.array([gdiff(x[k], food, vasc) for x in BL.values()])
        okb = np.isfinite(Db)
        Db = Db[okb]
        Lk = L[k].loc[targets]
        perm = np.array([Lk.loc[list(F)].mean() - Lk.loc[[t for t in targets if t not in F]].mean()
                         for F in itertools.combinations(targets, len(food))])
        summ.append(dict(readout=k, D=obs[k], lo=np.quantile(Db, 0.025), hi=np.quantile(Db, 0.975),
                         p_boot=(1 + np.sum(Db <= 0)) / (1 + len(Db)), n_boot_ok=int(len(Db)),
                         n_boot_failed=n_fail + int((~okb).sum()), p_perm=np.mean(perm >= obs[k] - 1e-4),
                         n_perm=len(perm), mean_loss_food=Lk.loc[food].mean(), mean_loss_vasc=Lk.loc[vasc].mean(),
                         n_food=len(food), n_vasc=len(vasc), skipped=";".join(skipped)))
    summ = pd.DataFrame(summ)
    summ.to_csv(out / "a33_summary.csv", index=False)
    print(summ.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
