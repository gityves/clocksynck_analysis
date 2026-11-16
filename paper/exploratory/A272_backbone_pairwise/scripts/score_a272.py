"""A272: is the stromal-vascular backbone of the healthy atlases seen without a host? Pairwise, anchor-free readout.

Registered in paper/docs/registered_experiments.md (A272). Atlases, rows and strata: those of A150 (the gate-passing
A17 healthy atlases with a genome-wide twin). No host and no anchors: for every pair of compartments (a, b) of an atlas,
the Poisson-corrected complex correlation of their clock-plane states over the donors that have both, pooled over
strata (states centred within stratum):

    r_ab = sum_s sum_d a_d conj(b_d) / sqrt(E_a E_b),   E_x = sum_s [sum_d |x_d|^2 - (n_s - 1) mean_s(v_x)]

|r_ab| is the shared state of the two compartments, arg r_ab their offset (h, a minus b). One projector per atlas and
gene set. The same statistic on 200 expression-matched and 200 stress-responsive pseudo-clock sets (same donors).

  T_SV    mean over atlases of the median |r| over pairs (stromal-vascular, X), X = any other non-lymphocyte family
  T_LY    the same for pairs (lymphocyte, X)
  H272.1  T_SV of the clock above the q95 of both companion families (set j pooled across atlases as set j)
  H272.2  D = T_SV - T_LY of the clock above the q95 of D in both companion families (atlases with lymphocytes)
  D272.*  family x family matrix, per-compartment centrality, offsets with donor-bootstrap intervals for pairs above
          both companions, triangle closure, ambient over-correction (clock), uncorrected correlation

    python score_a272.py --stage design                      # pairs, donors, degrees of freedom; no expression read into a statistic
    python score_a272.py --stage score [--placebo] [--n-pseudo 200] [--n-boot 500] [--jobs 4] [--only NAME ...]
`--placebo`: the donors of every compartment are permuted within stratum before pairing (code test; output/placebo/).
Outputs: paper/exploratory/A272_backbone_pairwise/output/{design,pairs,tests,family_matrix,centrality,offsets,
closure,sets}.csv (aggregates; public atlases).
"""

from __future__ import annotations

import importlib.util
import itertools
import sys
import time
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


A150 = _load("a272_a150", PAPER / "exploratory/A150_atlas_specificity/scripts/compute_a150.py")

AMENDMENT = "A272_backbone_pairwise"
SEED = 272
ATLASES = list(A150.ATLASES)
META = A150.META
CLOCK = A150.CLOCK
SV = ("fibroblast_stromal", "endothelial", "mural")
LY = "lymphocyte"
KINDS = {"matched": 1, "responsive": 2}
MIN_STRATUM, MIN_DF = 3, 15
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.0f} s]", *a, flush=True)


def group(family: str) -> str:
    return "SV" if family in SV else ("LY" if family == LY else "X")


# ---- design ----------------------------------------------------------------------------------------------------

def pair_index(pb: pd.DataFrame) -> list:
    """Per pair of compartments: row indices (a rows, b rows, stratum codes) of the donors having both, restricted to
    strata with >= MIN_STRATUM shared donors; kept when the pooled degrees of freedom reach MIN_DF."""
    pb = pb.reset_index(drop=True)
    rows = {c: g.set_index("donor") for c, g in pb.assign(_i=np.arange(len(pb))).groupby("compartment")}
    out = []
    for a, b in itertools.combinations(sorted(rows), 2):
        d = rows[a].index.intersection(rows[b].index)
        if len(d) == 0:
            continue
        ia, ib, st = rows[a].loc[d, "_i"].to_numpy(), rows[b].loc[d, "_i"].to_numpy(), rows[a].loc[d, "stratum"].to_numpy()
        codes, counts = np.unique(st, return_inverse=True)[1], None
        counts = np.bincount(codes)
        keep = counts[codes] >= MIN_STRATUM
        if not keep.any():
            continue
        codes = np.unique(codes[keep], return_inverse=True)[1]
        df = int(keep.sum() - codes.max() - 1)
        if df >= MIN_DF:
            out.append(dict(a=a, b=b, ia=ia[keep], ib=ib[keep], s=codes, n=int(keep.sum()), df=df,
                            n_strata=int(codes.max() + 1)))
    return out


# ---- statistic -------------------------------------------------------------------------------------------------

def _centre(z, s):
    m = np.zeros(s.max() + 1, complex)
    np.add.at(m, s, z)
    return z - (m / np.bincount(s))[s]


def pair_r(z, v, pr, corrected=True) -> complex:
    """Pooled complex correlation of a pair; NaN when a corrected variance is not positive."""
    s = pr["s"]
    a, b = _centre(z[pr["ia"]], s), _centre(z[pr["ib"]], s)
    ea, eb = np.sum(np.abs(a) ** 2), np.sum(np.abs(b) ** 2)
    if corrected:
        ns = np.bincount(s)
        for x, idx in (("a", pr["ia"]), ("b", pr["ib"])):
            mv = np.bincount(s, weights=v[idx]) / ns
            corr = np.sum((ns - 1) * mv)
            if x == "a":
                ea -= corr
            else:
                eb -= corr
    if ea <= 0 or eb <= 0:
        return complex(np.nan, np.nan)
    return np.sum(a * np.conj(b)) / np.sqrt(ea * eb)


def atlas_states(pb: pd.DataFrame, wf) -> tuple:
    g = list(wf.genes)
    L = P.projector(pb[g], wf, counts=True)
    beta, v = P.states(pb, L, wf, depth=pb.depth.to_numpy(float))
    return beta[:, 0] + 1j * beta[:, 1], v


def all_pairs(pb, wf, pairs, corrected=True) -> np.ndarray:
    z, v = atlas_states(pb, wf)
    return np.array([pair_r(z, v, pr, corrected) for pr in pairs])


def _chunk(pb, wfs, pairs):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    return np.array([all_pairs(pb[META + list(w.genes)], w, pairs) for w in wfs])


def pseudo_sets(A: dict, k: int, n: int, kind: str) -> list:
    pb, genes = A["pb"], A["genes"]
    X = pb[genes].to_numpy(float)
    frac = (X / pb.depth.to_numpy(float)[:, None]).mean(0)
    keep = (X.mean(0) >= 1) | np.isin(genes, CLOCK)
    lm = pd.Series(np.log(frac + 1e-12), index=genes)[keep]
    return P.pseudo_waveforms(lm, n, np.random.default_rng([SEED, k, KINDS[kind]]), cs.gtex(), kind)


def placebo(pb: pd.DataFrame, rng) -> pd.DataFrame:
    """Donor labels permuted within compartment x stratum: pairing between compartments is broken."""
    pb = pb.copy()
    for _, idx in pb.groupby(["compartment", "stratum"]).indices.items():
        pb.iloc[idx, pb.columns.get_loc("donor")] = rng.permutation(pb.donor.to_numpy()[idx])
    return pb


# ---- summaries -------------------------------------------------------------------------------------------------

def klass(ga: str, gb: str) -> str:
    s = {ga, gb}
    return "SV_X" if s == {"SV", "X"} else ("LY_X" if s == {"LY", "X"} else ("SV_LY" if s == {"SV", "LY"} else ga + "_" + gb))


def t_stat(absr: np.ndarray, meta: pd.DataFrame, cls: str, atlases) -> float:
    """Mean over atlases of the median |r| over the pairs of a class (NaN pairs left out)."""
    v = [np.nanmedian(absr[((meta.atlas == n) & (meta.klass == cls)).to_numpy()]) for n in atlases]
    return float(np.mean(v))


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("design", "score"), required=True)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--chunk", type=int, default=10)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--placebo", action="store_true")
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(AMENDMENT)
    if args.placebo:
        out = out / "placebo"
    out.mkdir(parents=True, exist_ok=True)
    names = [n for n in ATLASES if args.only is None or n in args.only]
    rng = np.random.default_rng([SEED, 99])

    atlas, pairs, meta = {}, {}, []
    for n in names:
        A = A150.load_atlas(n, None, 0)
        A["pb"] = A["pb"].reset_index(drop=True)
        if args.placebo:
            A["pb"] = placebo(A["pb"], rng)
        atlas[n], pairs[n] = A, pair_index(A["pb"])
        for pr in pairs[n]:
            fa, fb = A150.FAM.family(pr["a"]), A150.FAM.family(pr["b"])
            meta.append(dict(atlas=n, a=pr["a"], b=pr["b"], family_a=fa, family_b=fb, klass=klass(group(fa), group(fb)),
                             n=pr["n"], df=pr["df"], n_strata=pr["n_strata"]))
        log(n, "rows", len(A["pb"]), "donors", A["pb"].donor.nunique(), "pairs", len(pairs[n]))
    meta = pd.DataFrame(meta)
    design = (meta.groupby(["atlas", "klass"]).agg(pairs=("n", "size"), n_median=("n", "median"), df_min=("df", "min"))
              .reset_index())
    design.to_csv(out / "design.csv", index=False)
    print(design.to_string())
    if args.stage == "design":
        return

    # clock (corrected, uncorrected, ambient over-correction)
    r_clock = np.concatenate([all_pairs(atlas[n]["pb"], cs.gtex(), pairs[n]) for n in names])
    r_unc = np.concatenate([all_pairs(atlas[n]["pb"], cs.gtex(), pairs[n], corrected=False) for n in names])
    r_amb = []
    for n in names:
        corr, _, _ = A150.ambient_correct(atlas[n])
        r_amb.append(all_pairs(corr.reset_index(drop=True), cs.gtex(), pairs[n]))
    r_amb = np.concatenate(r_amb)

    # companions
    R = {}
    for kind in KINDS:
        tasks = []
        for k, n in enumerate(names):
            wfs = pseudo_sets(atlas[n], ATLASES.index(n), args.n_pseudo, kind)
            for i in range(0, args.n_pseudo, args.chunk):
                w = wfs[i:i + args.chunk]
                cols = sorted({g for x in w for g in x.genes})
                tasks.append((n, i, atlas[n]["pb"][META + cols], w, pairs[n]))
        res = Parallel(n_jobs=args.jobs, verbose=0)(delayed(_chunk)(*t[2:]) for t in tasks)
        mats = {n: np.full((args.n_pseudo, len(pairs[n])), np.nan, complex) for n in names}
        for t, m in zip(tasks, res):
            mats[t[0]][t[1]:t[1] + len(m)] = m
        R[kind] = np.concatenate([mats[n] for n in names], axis=1)
        log(kind, "sets done")

    # pair table
    absc = np.abs(r_clock)
    pt = meta.copy()
    pt["r_abs"], pt["offset_h"] = absc, np.angle(r_clock) * HOURS
    pt["r_abs_uncorrected"], pt["r_abs_ambient"] = np.abs(r_unc), np.abs(r_amb)
    for kind in KINDS:
        a = np.abs(R[kind])
        pt[f"{kind}_median"], pt[f"{kind}_q95"] = np.nanmedian(a, 0), np.nanquantile(a, 0.95, 0)
        pt[f"{kind}_valid"] = np.isfinite(a).sum(0)
    pt["above_both"] = (pt.r_abs > pt.matched_q95) & (pt.r_abs > pt.responsive_q95)

    # H272.1, H272.2
    with_ly = [n for n in names if ((meta.atlas == n) & (meta.klass == "LY_X")).any()]
    tests, sets = [], []
    for kind in KINDS:
        a = np.abs(R[kind])
        tsv = np.array([t_stat(a[j], meta, "SV_X", names) for j in range(len(a))])
        dj = np.array([t_stat(a[j], meta, "SV_X", with_ly) - t_stat(a[j], meta, "LY_X", with_ly) for j in range(len(a))])
        c_sv = t_stat(absc, meta, "SV_X", names)
        c_d = t_stat(absc, meta, "SV_X", with_ly) - t_stat(absc, meta, "LY_X", with_ly)
        for h, c, v in (("H272.1 T_SV", c_sv, tsv), ("H272.2 D", c_d, dj)):
            ok = np.isfinite(v)
            tests.append(dict(test=h, family=kind, clock=c, pseudo_median=float(np.nanmedian(v)),
                              pseudo_q95=float(np.nanquantile(v, 0.95)), n_sets=int(ok.sum()),
                              p=float((1 + np.sum(v[ok] >= c)) / (1 + ok.sum())), above_q95=bool(c > np.nanquantile(v, 0.95))))
        sets.append(pd.DataFrame(dict(family=kind, set=range(len(a)), T_SV=tsv, D=dj)))
    tests = pd.DataFrame(tests)
    for h in tests.test.unique():
        tests.loc[tests.test == h, "passes"] = bool(tests.loc[tests.test == h, "above_q95"].all())
    for lab, r in (("uncorrected", r_unc), ("ambient", r_amb)):
        a = np.abs(r)
        tests = pd.concat([tests, pd.DataFrame([dict(test=f"descriptive T_SV {lab}", clock=t_stat(a, meta, "SV_X", names)),
                                                dict(test=f"descriptive D {lab}", clock=t_stat(a, meta, "SV_X", with_ly)
                                                     - t_stat(a, meta, "LY_X", with_ly))])], ignore_index=True)
    tests.to_csv(out / "tests.csv", index=False)
    pd.concat(sets).to_csv(out / "sets.csv", index=False)

    # family x family matrix (median over pairs, pooled over atlases) and per-compartment centrality
    fam = pt.assign(f1=np.minimum(pt.family_a, pt.family_b), f2=np.maximum(pt.family_a, pt.family_b))
    fm = fam.groupby(["f1", "f2"]).agg(pairs=("r_abs", "size"), atlases=("atlas", "nunique"), clock=("r_abs", "median"),
                                       matched=("matched_median", "median"), responsive=("responsive_median", "median"),
                                       share_above_both=("above_both", "mean")).reset_index()
    fm.to_csv(out / "family_matrix.csv", index=False)
    cen = []
    for (n, c), _ in pd.concat([pt[["atlas", "a"]].rename(columns={"a": "c"}), pt[["atlas", "b"]].rename(columns={"b": "c"})]).groupby(["atlas", "c"]):
        m = ((pt.atlas == n) & ((pt.a == c) | (pt.b == c))).to_numpy()
        other = np.where(pt.a[m] == c, pt.family_b[m], pt.family_a[m])
        m2 = m.copy()
        m2[np.where(m)[0][other == A150.FAM.family(c)]] = False  # partners of another family only
        if not m2.any():
            continue
        row = dict(atlas=n, compartment=c, family=A150.FAM.family(c), partners=int(m2.sum()), clock=float(np.nanmedian(absc[m2])))
        for kind in KINDS:
            v = np.nanmedian(np.abs(R[kind])[:, m2], 1)
            row[f"{kind}_median"], row[f"{kind}_q95"] = float(np.nanmedian(v)), float(np.nanquantile(v, 0.95))
        row["excess_matched"], row["excess_responsive"] = row["clock"] - row["matched_median"], row["clock"] - row["responsive_median"]
        cen.append(row)
    pd.DataFrame(cen).to_csv(out / "centrality.csv", index=False)

    # offsets of the pairs above both companions: donor bootstrap within stratum
    brng = np.random.default_rng([SEED, 7])
    lo, hi = np.full(len(pt), np.nan), np.full(len(pt), np.nan)
    off = 0
    for n in names:
        z, v = atlas_states(atlas[n]["pb"], cs.gtex())
        for k, pr in enumerate(pairs[n]):
            i = off + k
            if not pt.above_both.iloc[i]:
                continue
            ang = []
            for _ in range(args.n_boot):
                idx = np.concatenate([brng.choice(np.where(pr["s"] == s)[0], (pr["s"] == s).sum()) for s in range(pr["s"].max() + 1)])
                r = pair_r(z, v, dict(ia=pr["ia"][idx], ib=pr["ib"][idx], s=pr["s"][idx]))
                ang.append(np.angle(r * np.conj(r_clock[i])) * HOURS if np.isfinite(r) else np.nan)
            q = np.nanquantile(ang, [0.05, 0.95])
            lo[i], hi[i] = pt.offset_h.iloc[i] + q[0], pt.offset_h.iloc[i] + q[1]
        off += len(pairs[n])
    pt["offset_lo90"], pt["offset_hi90"] = lo, hi
    pt.to_csv(out / "pairs.csv", index=False)

    # triangle closure among pairs above both companions
    clo = []
    for n in names:
        g = pt[(pt.atlas == n) & pt.above_both]
        o = {(a, b): h for a, b, h in zip(g.a, g.b, g.offset_h)}
        o.update({(b, a): -h for (a, b), h in list(o.items())})
        comps = sorted({c for k in o for c in k})
        for a, b, c in itertools.combinations(comps, 3):
            if (a, b) in o and (b, c) in o and (a, c) in o:
                e = (o[(a, b)] + o[(b, c)] - o[(a, c)] + 12) % 24 - 12
                clo.append(dict(atlas=n, a=a, b=b, c=c, closure_h=e))
    clo = pd.DataFrame(clo, columns=["atlas", "a", "b", "c", "closure_h"])
    clo.to_csv(out / "closure.csv", index=False)
    log("done")
    print(tests.to_string())
    print(fm.to_string())
    if len(clo):
        print("closure |h| median", float(clo.closure_h.abs().median()), "triangles", len(clo))


if __name__ == "__main__":
    main()
