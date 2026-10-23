"""Decomposition of the nominal coupling breaks (A98, descriptive): is a loss one of amplitude, own noise, coherence,
off-plane variance, clock-shaped share or alignment with the donor's host?

Pre-registered as A98 (2026-10-18; docs/PLAN_internal_desync.md (archive), method addendum in the legacy folder). Contrasts
(target class, specificity classes fixed before the run): Natri ILD (lung:T; lung:fibroblast), HLCA fibrotic (lung:T;
lung:fibroblast), Parkinson DMV Braak LB 5-6 vs 0 (brain:DMV_neuron_other; DMV oligodendrocyte, astrocyte, microglia),
T2D islets (islet:beta; alpha, stellate), atopic-dermatitis skin (skin:T; skin:fibroblast), with the A17 data, labels and
learned anchors. Per class and group, on the pooled anchored host: relative amplitude |z|, own noise psi, coupling
(`Dataset.fit`), log plane and log off-plane variance and gene coherence (`clocksync.decompose.components`, floored at
log 1e-9), their unfloored plane and off-plane variances, clock-shaped share (plane / total variance of the
Poisson-corrected gene covariance) and alignment (|Poisson-corrected complex correlation| of the class state with the
mean state of the other anchors, `stats.ccorr`; undefined when a corrected variance is <= 0). Estimand: case minus
control of every quantity. Interval: donor bootstrap within group (`--n-boot`, 200 by default; the pre-registration
asks for 1,000), percentile 95 %. Reading rule of A87: a quantity contributes when its interval excludes no change in
the direction of loss (amplitude, coupling, log plane, coherence, share, alignment down; psi and log off-plane up);
"floor-limited" when a floored log variance sits at its floor in either group or in > 5 % of the draws.

    python compute_decomposition.py [--subsample N] [--seed S] [--jobs 4] [--only NAME ...] [--n-boot N]
Outputs (paper/results/figure_5/output/lung_breaks/): decomposition_groups.csv, decomposition_differences.csv
Legacy: legacy/desync/validation/breaks_decomposition/scripts/a98_breaks.py -> output/a98_{groups,differences}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.decompose import components, logfrac  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "lung_breaks"
WF = cs.gtex()
GENES = list(WF.genes)
CONTRASTS = {
    "ild_natri": ("Natri ILD", "lung:T", ["lung:fibroblast"]),
    "hlca_fibrotic": ("HLCA fibrotic", "lung:T", ["lung:fibroblast"]),
    "pd": ("PD DMV Braak 5-6 vs 0", "brain:DMV_neuron_other",
           ["brain:DMV_oligodendrocyte", "brain:DMV_astrocyte", "brain:DMV_microglia"]),
    "islets_t2d": ("T2D islets", "islet:beta", ["islet:alpha", "islet:stellate"]),
    "skin_ad": ("AD skin lesional", "skin:T", ["skin:fibroblast"]),
}
LOSS = dict(amp=-1, coupling=-1, psi=+1, log_plane=-1, log_off=+1, plane=-1, off=+1, coherence=-1, share=-1, align=-1)
FLOORED = ("log_plane", "log_off")
LOG_FLOOR = float(np.log(1e-9))
QUANT = list(LOSS)
MIN_DONORS = 5


def _plane() -> np.ndarray:
    """local: not in clocksync. Orthogonal projector on the clock plane of the gene-centred design."""
    X = WF.X - WF.X.mean(0)
    return X @ np.linalg.solve(X.T @ X, X.T)


PL = _plane()


def centre_by(b: np.ndarray, strata: np.ndarray) -> np.ndarray:
    b = np.array(b, float, copy=True)
    for q in np.unique(strata):
        i = strata == q
        b[i] -= b[i].mean(0)
    return b


def share(y, pv, k) -> float:
    """local: not in clocksync (A86 clock-shaped share): ||P S P||^2 / ||S||^2 of the Poisson-corrected gene covariance."""
    S = y.T @ y / (len(y) - k) - np.diag(pv)
    return float(np.sum((PL @ S @ PL) ** 2) / np.sum(S ** 2))


def anchor_host(ds: cs.Dataset, anc: list, labelled) -> pd.DataFrame:
    """local: not in clocksync (A54 / A98 alignment host). Per donor: mean state of the anchors, each centred within
    its compartment x stratum over the labelled donors; variance sum(trace V) / n^2. Reads `Dataset.data` (clock-plane
    states of every pseudobulk)."""
    d = ds.data
    tr = np.trace(d.V, axis1=1, axis2=2)
    parts = []
    for a in anc:
        i = np.flatnonzero((d.compartment == a) & np.isin(d.resection, list(labelled)))
        if len(i):
            b = centre_by(d.beta[i], d.stratum[i])
            parts.append(pd.DataFrame(dict(donor=d.resection[i], x=b[:, 0], y=b[:, 1], v=tr[i])))
    if not parts:
        return None
    g = pd.concat(parts).groupby("donor")
    n = g.size()
    return pd.DataFrame(dict(x=g.x.mean(), y=g.y.mean(), v=g.v.sum() / n ** 2))


def class_stats(pb: pd.DataFrame, lab: pd.Series, anc: list, k: str) -> dict | None:
    """Quantities of class k per group on the pooled host of the other anchors. The two groups' pseudobulks of k are two
    non-anchor compartments of one Dataset ("k@control", "k@case"), so both are read on the same host fitted on all
    labelled donors. Donors without any other anchor are dropped."""
    anc = [a for a in anc if a != k]
    if len(anc) < 3:  # `Dataset.fit` needs >= 3 anchors (the legacy leave-one-out host also ran on 2)
        return None
    t = pb[pb.donor.isin(lab.index)]
    have = set(t[t.compartment.isin(anc)].donor)
    rk = t[(t.compartment == k) & t.donor.isin(have)]
    grp = rk.donor.map(lab)
    if min((grp == g).sum() for g in sh.GROUPS) < MIN_DONORS:
        return None
    tt = pd.concat([t[t.compartment.isin(anc)], rk.assign(compartment=k + "@" + grp.to_numpy())], ignore_index=True)
    ds = cs.Dataset(tt, anchors=anc, name=k)
    names = {g: f"{k}@{g}" for g in sh.GROUPS}
    f = ds.fit(targets=list(names.values()), n_boot=0, gate=False).table
    hd = anchor_host(ds, anc, lab.index)
    d = ds.data
    tr = np.trace(d.V, axis1=1, axis2=2)
    out = {}
    for g, nm in names.items():
        r = f.loc[nm]
        rows = rk[grp.to_numpy() == g]
        y, pv, kk = logfrac(rows, GENES, rows["stratum"].to_numpy() if "stratum" in rows else None)
        comp = components(y, pv, kk, WF)
        Sg = y.T @ y / (len(y) - kk) - np.diag(pv)
        plane = float(np.trace(PL @ Sg @ PL))
        al = np.nan
        i = np.flatnonzero(d.compartment == nm)
        m = np.isin(d.resection[i], hd.index) if hd is not None else np.zeros(len(i), bool)
        if m.sum() >= 5:
            i = i[m]
            q = d.stratum[i]
            bl = centre_by(d.beta[i], q)
            hh = hd.loc[d.resection[i]]
            bn = centre_by(hh[["x", "y"]].to_numpy(), q)
            ev = lambda b, v: np.sum(b[:, 0] ** 2 + b[:, 1] ** 2) / (len(b) - 1) - np.mean(v)
            if ev(bl, tr[i]) > 0 and ev(bn, hh.v.to_numpy()) > 0:
                al = abs(cs.stats.ccorr(bl, bn, tr[i], hh.v.to_numpy()))
        out[g] = dict(amp=float(np.sqrt(r.z2)), psi=float(r.psi), coupling=float(r.coupling), boundary=bool(r.boundary),
                      **comp, plane=plane, off=float(np.trace(Sg) - plane), share=share(y, pv, kk), align=al,
                      n=int(r.n))
    return out


def resample(pb: pd.DataFrame, lab: pd.Series, rng) -> tuple[pd.DataFrame, pd.Series]:
    """Donors drawn with replacement within each group; repeats relabelled `<donor>#<j>` (all their rows are copied)."""
    parts, labs, j = [], {}, 0
    by = {d: x for d, x in pb.groupby("donor")}
    for g in sh.GROUPS:
        ids = lab.index[lab == g].to_numpy()
        for d in rng.choice(ids, len(ids)):
            new = f"{d}#{j}"
            parts.append(by[d].assign(donor=new))
            labs[new] = g
            j += 1
    return pd.concat(parts, ignore_index=True), pd.Series(labs)


def run_contrast(name: str, n_boot: int, subsample, seed: int):
    limit_threads(1)
    label, target, spec = CONTRASTS[name]
    labels = None
    if name == "pd":
        def labels(idx):
            st = sh.braak_stage(idx)
            return st[st.isin(["0", "5-6"])].map({"0": "control", "5-6": "case"})
    pb, donors, lab, blocks, row = sh.load_target(name, subsample, seed, labels=labels)
    anc = [a for a in sh.anchor_set(name, "learned") if a in set(pb.compartment)]
    rng = np.random.default_rng(98 + list(CONTRASTS).index(name) + 1000 * seed)
    est, dif = [], []
    for role, k in [("target", target)] + [("specificity", s) for s in spec]:
        if k not in set(pb.compartment):
            continue
        obs = class_stats(pb, lab, anc, k)
        if obs is None:
            print(name, k, "not computed (< 5 donors in a group or < 3 other anchors)", flush=True)
            continue
        boots = {g: [] for g in sh.GROUPS}
        for _ in range(n_boot):
            b, bl = resample(pb, lab, rng)
            try:
                r = class_stats(b, bl, anc, k)
            except (ValueError, np.linalg.LinAlgError):
                continue
            if r is not None:
                for g in sh.GROUPS:
                    boots[g].append(r[g])
        boots = {g: pd.DataFrame(v) for g, v in boots.items()}
        nb = len(boots["case"])
        for g in sh.GROUPS:
            est.append(dict(contrast=name, role=role, compartment=k, group=g, donors=obs[g]["n"],
                            boundary=obs[g]["boundary"], **{q: obs[g][q] for q in QUANT}))
        for q in QUANT:
            b = (boots["case"][q] - boots["control"][q]).to_numpy(float) if nb else np.array([])
            b = b[np.isfinite(b)]
            d0 = obs["case"][q] - obs["control"][q]
            lo, hi = (np.quantile(b, 0.025), np.quantile(b, 0.975)) if len(b) >= 10 else (np.nan, np.nan)
            loss = LOSS[q]
            flo = np.nan
            if q in FLOORED and nb:
                bq = np.concatenate([boots[g][q].to_numpy(float) for g in sh.GROUPS])
                flo = float(np.mean(bq <= LOG_FLOOR + 1e-6))
            at_floor = q in FLOORED and (min(obs["case"][q], obs["control"][q]) <= LOG_FLOOR + 1e-6 or flo > 0.05)
            reading = ("floor-limited" if at_floor else "n/a" if not np.isfinite(lo) else
                       "contributes" if (lo > 0 if loss > 0 else hi < 0) else
                       "opposite" if (hi < 0 if loss > 0 else lo > 0) else "no change")
            dif.append(dict(contrast=name, role=role, compartment=k, quantity=q, control=obs["control"][q],
                            case=obs["case"][q], diff=d0, lo=lo, hi=hi, n_boot_ok=len(b), floor_frac=flo,
                            loss_direction="rise" if loss > 0 else "fall", reading=reading,
                            donors_control=obs["control"]["n"], donors_case=obs["case"]["n"]))
        print(name, k, flush=True)
    return pd.DataFrame(est), pd.DataFrame(dif)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="+", default=list(CONTRASTS))
    ap.add_argument("--n-boot", type=int, default=None, help="bootstrap draws (default 200; 20 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    nb = args.n_boot if args.n_boot is not None else run["n_boot"]
    res = Parallel(n_jobs=args.jobs)(delayed(run_contrast)(n, nb, args.subsample, args.seed) for n in args.only)
    est = pd.concat([e for e, _ in res], ignore_index=True)
    dif = pd.concat([d for _, d in res], ignore_index=True)
    est.to_csv(out / "decomposition_groups.csv", index=False)
    dif.to_csv(out / "decomposition_differences.csv", index=False)
    pd.set_option("display.width", 220)
    print(dif[dif.role == "target"][["contrast", "quantity", "control", "case", "diff", "lo", "hi", "reading"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
