"""A287: structured clock-gene covariance in tissue lymphocytes without a waveform: a different clock, or none?

Registered in paper/docs/registered_experiments.md (A287). The nine A281 pairs (lymphocyte, fibroblast-stromal on >= 15
shared donors). Per compartment: y = log((count + 0.5) / depth) of the clock genes, centred per gene within stratum and
per row over the genes; U = (lambda1 + lambda2) / trace of the gene correlation matrix of y. References: 200 matched sets
(A281 sets), per-gene permutations within strata, 200 responsive sets (reported). Structure: pi >= 0.95 and p_perm < 0.05.
Power: fibroblast-stromal rows thinned to the lymphocytes' median clock-gene count (20 replicates).

  H287.1  in powered atlases (thinned structure in >= 80 % of replicates), lymphocytes carry structure in >= half
  H287.2  (read only if H287.1 fails) lymphocytes do not in >= 2/3

    python score_a287.py [--design] [--jobs 6]
Outputs (../output/): design.csv, compartments.csv, thinned.csv, phases.csv, tests.csv
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync import pairwise as pw  # noqa: E402
from common import cli  # noqa: E402

spec = importlib.util.spec_from_file_location("a287_a281", PAPER / "exploratory/A281_lymphocyte_clock/scripts/score_a281.py")
A281 = importlib.util.module_from_spec(spec)
sys.modules["a287_a281"] = A281
spec.loader.exec_module(A281)
import cloudpickle  # noqa: E402

cloudpickle.register_pickle_by_value(A281)

SEED, N_SETS, N_PERM, N_THIN, N_PERM_THIN = 287, 200, 1000, 20, 200
CLOCK, WF = A281.CLOCK, A281.WF
META = ["donor", "compartment", "stratum", "n_cells", "depth"]


def prep(t, genes):
    """Clock-gene-like matrix of a set: log fractions centred per gene within stratum, then per row over the genes."""
    y = A281.y_of(t, genes)
    return y - y.mean(1, keepdims=True)


def U(y):
    """Share of the correlation matrix's variance in its two leading eigenvalues (nan without variance)."""
    sd = y.std(0)
    ok = sd > 0
    if ok.sum() < 3:
        return np.nan
    R = np.corrcoef(y[:, ok], rowvar=False)
    ev = np.linalg.eigvalsh(R)
    return float(ev[-2:].sum() / ev.sum())


def permuted(y, groups, rng):
    z = y.copy()
    for i in groups:
        for j in range(y.shape[1]):
            z[i, j] = y[rng.permutation(i), j]
    return z


def readout(t, sets, n_perm, rng, resp=None):
    """U of the clock, pi among the sets, p_perm (and responsive pi)."""
    t = t.reset_index(drop=True)
    y = prep(t, CLOCK)
    u = U(y)
    comp = np.array([U(prep(t, s)) for s in sets])
    groups = A281.groups_of(t)
    null = np.array([U(permuted(y, groups, rng)) for _ in range(n_perm)])
    out = dict(U=u, pi=float(np.mean(comp[np.isfinite(comp)] < u)), matched_median=float(np.nanmedian(comp)),
               p_perm=float((1 + np.sum(null >= u)) / (1 + n_perm)))
    if resp is not None:
        cr = np.array([U(prep(t, s)) for s in resp])
        out.update(responsive_pi=float(np.mean(cr[np.isfinite(cr)] < u)), responsive_median=float(np.nanmedian(cr)))
    out["structure"] = bool(out["pi"] >= 0.95 and out["p_perm"] < 0.05)
    return out, y


def phases(y):
    """Angle of each gene's loading on the two leading eigenvectors of the correlation matrix."""
    R = np.corrcoef(y, rowvar=False)
    w, v = np.linalg.eigh(R)
    return np.angle(v[:, -1] + 1j * v[:, -2])


def circ_corr(a, b):
    """|complex correlation| of exp(i a) with exp(i b) (rotation and reflection invariant up to sign)."""
    za, zb = np.exp(1j * a), np.exp(1j * b)
    za, zb = za - za.mean(), zb - zb.mean()
    plus = abs(np.sum(za * np.conj(zb))) / np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zb) ** 2))
    zbc = np.conj(zb)
    minus = abs(np.sum(za * np.conj(zbc))) / np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zbc) ** 2))
    return float(max(plus, minus))


def pair_task(n, ly, sv, tl, ts, matched, resp, seed):
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(seed)
    d = sorted(set(tl.donor) & set(ts.donor))
    tl = tl[tl.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    ts = ts[ts.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    target = float(tl[CLOCK].sum(1).median())
    rows, ph = [], []
    for role, comp, t in (("lymphocyte", ly, tl), ("fibroblast_stromal", sv, ts)):
        r, y = readout(t, matched, N_PERM, rng, resp)
        rows.append(dict(atlas=n, role=role, compartment=comp, donors=len(d),
                         clock_count_median=float(t[CLOCK].sum(1).median()), **r))
        if r["structure"]:
            p = phases(y)
            ph.append(dict(atlas=n, role=role, compartment=comp,
                           circ_corr_waveform=circ_corr(p, np.asarray(WF.psi, float)),
                           **{g: float(x) for g, x in zip(CLOCK, p)}))
    genes = sorted({g for s in matched for g in s} | set(CLOCK))
    thin = []
    for rep in range(N_THIN):
        rr = np.random.default_rng([SEED, seed[1], 100 + rep])
        tt = ts.copy()
        s = tt[CLOCK].sum(1).to_numpy(float)
        p = np.minimum(1.0, target / np.maximum(s, 1))
        tt[genes] = rr.binomial(tt[genes].to_numpy(np.int64), p[:, None])
        tt["depth"] = tt.depth.to_numpy(float) * p
        r, _ = readout(tt, matched, N_PERM_THIN, rr)
        thin.append(dict(atlas=n, rep=rep, **r))
    return rows, ph, thin


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--design", action="store_true")
    args = ap.parse_args()
    cli.setup(args)
    out = HERE.parents[1] / "output"
    out.mkdir(parents=True, exist_ok=True)
    prs = pd.read_csv(PAPER / "exploratory/A281_lymphocyte_clock/output/pairs.csv")
    A281.A150.ATLASES.update(A281.ATLASES)
    keys = list(A281.ATLASES)
    tasks, design = [], []
    for r in prs.itertuples():
        k = keys.index(r.atlas)
        A = A281.A150.load_atlas(r.atlas, None, 0)
        pb = A["pb"].reset_index(drop=True)
        matched = [list(w.genes) for w in A281.A272.pseudo_sets(A, k, N_SETS, "matched")]
        sub = pb[pb.compartment.isin([r.lymphocyte, r.stromal])].reset_index(drop=True)
        genes = [c for c in sub.columns if c not in META]
        resp = [list(w.genes) for w in pw.companions(sub[META + genes], N_SETS, "responsive",
                                                     np.random.default_rng([SEED, k, 1]))]
        cols = sorted({g for s in matched + resp for g in s} | set(CLOCK))
        sub = sub[META + cols]
        design.append(dict(atlas=r.atlas, lymphocyte=r.lymphocyte, stromal=r.stromal, shared_donors=r.shared_donors))
        tasks.append(delayed(pair_task)(r.atlas, r.lymphocyte, r.stromal, sub[sub.compartment == r.lymphocyte],
                                        sub[sub.compartment == r.stromal], matched, resp, [SEED, k]))
    pd.DataFrame(design).to_csv(out / "design.csv", index=False)
    print(pd.DataFrame(design).to_string(), flush=True)
    if args.design:
        return
    res = Parallel(n_jobs=args.jobs)(tasks)
    comp = pd.DataFrame([x for r in res for x in r[0]])
    ph = pd.DataFrame([x for r in res for x in r[1]])
    thin = pd.DataFrame([x for r in res for x in r[2]])
    rate = thin.groupby("atlas").structure.mean().rename("thinned_rate")
    ly = comp[comp.role == "lymphocyte"].set_index("atlas").join(rate)
    powered = ly[ly.thinned_rate >= 0.8]
    n = len(powered)
    share = float(powered.structure.mean()) if n else np.nan
    h1 = bool(n >= 3 and share >= 0.5)
    h2 = (bool(1 - share >= 2 / 3) if n >= 3 and not h1 else None)
    tests = pd.DataFrame([dict(atlases=len(ly), powered=n, lymph_structure=int(powered.structure.sum()) if n else 0,
                               share=share, H287_1_pass=h1 if n >= 3 else None, H287_2_pass=h2)])
    comp.to_csv(out / "compartments.csv", index=False)
    ph.to_csv(out / "phases.csv", index=False)
    thin.to_csv(out / "thinned.csv", index=False)
    tests.to_csv(out / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(comp.round(3).to_string(), rate.round(3).to_string(), ph.round(3).to_string(), tests.T.to_string(), sep="\n")


if __name__ == "__main__":
    main()
