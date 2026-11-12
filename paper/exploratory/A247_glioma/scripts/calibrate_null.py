"""A247 note 1: which pairing statistic and null are calibrated when two classes of a patient share their depth?

Healthy-brain calibration only (control donors of PsychAD; microglia and oligodendrocytes thinned to the clock-gene counts of
the P247.0 design; one thinning probability per sample, see `thin_total`). No statistic is computed on GBmap pairings: GBmap supplies the
thinning targets and the strata. For each subsample, three arms (true pairs, donors mis-paired, true pairs at the atlas's own
depth) and every combination of
  statistic  raw     |sum za conj(zb)| / sqrt(sum |za|^2 sum |zb|^2), states centred within stratum (registered)
             phase   the same on unit-modulus states
             resid   states residualised within stratum on the log clock-gene count of their own pseudobulk, then raw
             resid2  residualised on the log clock-gene counts of both pseudobulks
  null       perm    rows of the second class permuted within stratum (registered)
             flip    random sign per patient on the first class (moduli stay paired)
             block   permutation within stratum restricted to blocks of 4 neighbours in the second class's clock-gene count
the share of subsamples with pooled p < 0.05 -> output/calibration_null.csv.

`--design H247.1 | H247.2`: the same for the design of a malignant test, with the atlas oligodendrocytes thinned to the malignant
clock-gene counts and the atlas microglia to those of the other class (`ROLE`): the pass rate then answers "would malignant cells
tied to the patient's glia as healthy glia are tied to each other be detected at this depth?". `--assay nucleus | cell`
restricts the design to those strata.

    python calibrate_null.py [--design P247.0] [--assay nucleus] [--jobs 3] [--n-sub 200] [--n-perm 500]
"""

from __future__ import annotations

import sys
from itertools import product
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compute_a247 as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

STATS, NULLS = ("raw", "phase", "resid", "resid2"), ("perm", "flip", "block")
BLOCK = 4
# design -> GBmap classes whose counts are the thinning targets of (atlas microglia, atlas oligodendrocytes)
ROLE = {"P247.0": ("myeloid", "oligodendrocyte"), "H247.1": ("myeloid", "malignant"), "H247.2": ("oligodendrocyte", "malignant")}


def thin_total(x: pd.DataFrame, target: np.ndarray, rng) -> pd.DataFrame:
    """Binomial thinning with ONE probability per sample (total clock-gene count of the target over that of the donor, at most
    1), so the donor keeps its own gene pattern. `compute_a247.thin` thins gene by gene to the target's counts, which copies
    the GBmap patient's clock-gene pattern into the atlas sample: its pass and mis-paired rates then reflect GBmap, not the atlas."""
    c = x[C.CLOCK].to_numpy(float)
    p = np.minimum(1.0, target.sum(1) / np.maximum(c.sum(1), 1))[:, None]
    new = rng.binomial(c.astype(np.int64), np.broadcast_to(p, c.shape)).astype(float)
    out = x[C.CLOCK].copy()
    out[:] = new
    out["depth"] = x.depth.to_numpy(float) * p[:, 0]
    return out


def _resid(z, X):
    X = np.column_stack([np.ones(len(z))] + list(X))
    if len(z) <= X.shape[1]:
        return np.zeros_like(z)
    return z - X @ np.linalg.lstsq(X, z, rcond=None)[0]


def transform(za, zb, la, lb, how):
    if how == "resid":
        za, zb = _resid(za, [la]), _resid(zb, [lb])
    elif how == "resid2":
        za, zb = _resid(za, [la, lb]), _resid(zb, [la, lb])
    za, zb = za - za.mean(), zb - zb.mean()
    if how == "phase":
        za, zb = za / np.maximum(abs(za), 1e-12), zb / np.maximum(abs(zb), 1e-12)
        za, zb = za - za.mean(), zb - zb.mean()
    return za, zb


def stat_null(za, zb, lb, null, n_perm, rng):
    """(observed, null draws) of the raw statistic for already transformed states."""
    den = np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zb) ** 2))
    if den == 0:
        return 0.0, np.zeros(n_perm)
    obs = abs(np.sum(za * np.conj(zb))) / den
    n = len(za)
    if null == "flip":
        s = rng.choice([-1.0, 1.0], size=(n_perm, n))
        return obs, np.abs((s * (za * np.conj(zb))[None, :]).sum(1)) / den
    if null == "perm":
        perms = np.argsort(rng.random((n_perm, n)), axis=1)
    else:
        order = np.argsort(lb)
        perms = np.tile(np.arange(n), (n_perm, 1))
        for b in range(0, n, BLOCK):
            i = order[b:b + BLOCK] if n - b >= 2 * BLOCK or b + BLOCK >= n else order[b:]   # last block takes the remainder
            perms[:, i] = i[np.argsort(rng.random((n_perm, len(i))), axis=1)]
            if len(i) > BLOCK:
                break
    return obs, np.abs((za[None, :] * np.conj(zb[perms])).sum(1)) / den


def one(task):
    ds, k, n_perm = task
    mg, ol = C.G["atlas"][ds]
    rng = np.random.default_rng([C.SEED, 71, k, sum(map(ord, ds))])
    slots = C.G["slots"]
    pick = rng.permutation(len(mg))[:len(slots)]
    pick2 = pick[rng.permutation(len(pick))]
    out = {}
    for arm, po in (("pass", pick), ("mispaired", pick2), ("native", pick)):
        big = arm == "native"
        m = thin_total(mg.iloc[pick], np.full_like(C.G["tm"], 1e12) if big else C.G["tm"], rng)
        o = thin_total(ol.iloc[po], np.full_like(C.G["to"], 1e12) if big else C.G["to"], rng)
        df = pd.concat([x.reset_index(drop=True).assign(donor=slots.slot.to_numpy(), compartment=cls, stratum=slots.stratum.to_numpy())
                        for cls, x in (("myeloid", m), ("oligodendrocyte", o))], ignore_index=True)
        tab = C.Tab(df)
        designs = [d for d in C.make_designs("P247.0", tab, pd.DataFrame(index=slots.slot)) if d["eligible"]]  # atlas classes carry the P247.0 labels
        units = []
        for d in designs:
            s, pts = d["stratum"], d["pts"]
            la, lb = (np.log(tab.rows(c, s).loc[pts, C.CLOCK].to_numpy(float).sum(1) + 1) for c in ("myeloid", "oligodendrocyte"))
            units.append((len(pts), C.state(tab, "myeloid", s, pts), C.state(tab, "oligodendrocyte", s, pts), la, lb))
        w = np.array([u[0] for u in units], float)
        for st, nu in product(STATS, NULLS):
            r = [stat_null(*transform(za, zb, la, lb, st), lb, nu, n_perm, rng) for _, za, zb, la, lb in units]
            obs = np.sum(w * [x[0] for x in r]) / w.sum()
            nul = (w[:, None] * np.array([x[1] for x in r])).sum(0) / w.sum()
            out[(arm, st, nu)] = (1 + np.sum(nul >= obs)) / (1 + n_perm)
    return ds, k, out


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--design", choices=list(ROLE), default="P247.0")
    ap.add_argument("--assay", choices=["all", "cell", "nucleus"], default="all")
    args = ap.parse_args(argv)
    cli.setup(args)
    df, pat = C.load_table()
    tab = C.Tab(df)
    rows = [(d["stratum"], p) for d in C.make_designs(args.design, tab, pat) if d["eligible"] for p in d["pts"]
            if args.assay == "all" or d["stratum"].endswith("|" + args.assay)]
    slots = pd.DataFrame(rows, columns=["stratum", "patient"]).assign(slot=lambda t: [f"slot{i}" for i in range(len(t))])
    C.G["slots"] = slots
    for key, cls in zip(("tm", "to"), ROLE[args.design]):
        r = df[df.compartment == cls].set_index(["stratum", "donor"])
        C.G[key] = r.loc[list(zip(slots.stratum, slots.patient)), C.CLOCK].to_numpy(float)
    C.G["atlas"] = {ds: x for ds in ("psychad_hbcc", "psychad_mssm") if len((x := C.atlas_pairs(ds))[0]) >= len(slots)}
    print(f"{args.design} / {args.assay}: {len(slots)} patients in {slots.stratum.nunique()} strata; atlases {list(C.G['atlas'])}", flush=True)
    tasks = [(ds, k, args.n_perm) for ds in C.G["atlas"] for k in range(args.n_sub)]
    with Pool(args.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        res = pool.map(one, tasks, chunksize=4)
    R = pd.DataFrame([dict(atlas=ds, k=k, arm=a, statistic=s, null=n, p=p) for ds, k, o in res for (a, s, n), p in o.items()])
    T = (R.assign(hit=R.p < .05).groupby(["atlas", "statistic", "null", "arm"]).hit.mean().unstack("arm")
         .rename(columns={"pass": "pass_rate", "mispaired": "false_positive_rate", "native": "pass_rate_native_depth"}).reset_index())
    T["design"], T["assay"], T["patients"], T["strata"] = args.design, args.assay, len(slots), slots.stratum.nunique()
    T["subsamples"], T["n_perm"] = args.n_sub, args.n_perm
    T.to_csv(C.OUT / f"calibration_null_{args.design}_{args.assay}.csv", index=False)
    print(T.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
