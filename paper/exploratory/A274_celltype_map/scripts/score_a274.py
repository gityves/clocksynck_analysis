"""A274: replication of the A272 cell-type ranking in atlases A272 did not use, and the shared clock state between
cortical cell classes. Anchor-free pairwise readout (the A272 statistic after its design addendum).

Registered in paper/docs/registered_experiments.md (A274, staged; design addendum before scoring). Atlases: the A17
healthy targets whose genome-wide twin was built for this amendment; an atlas enters only if the genome table
reproduces the clock-gene table (twin check, `--stage design`).

  H274.1  D = T_SV - T_LY over the body atlases with >= 3 pairs (SV, X) and >= 3 pairs (LY, X): above the q95 of both
          companion families
  H274.2  T_SV over the body atlases with >= 3 pairs (SV, X): above the q95 of the matched family (responsive reported)
  H274.3  T_CX = mean over the cortex atlases of the median |r| over pairs of compartments of different families:
          above the q95 of both companion families

    python score_a274.py --stage design
    python score_a274.py --stage score [--placebo] [--n-pseudo 200] [--n-boot 500] [--jobs 4] [--only NAME ...]
Outputs: paper/exploratory/A274_celltype_map/output/{design,twin,pairs,tests,family_matrix,sets}.csv
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

spec = importlib.util.spec_from_file_location("a274_a272", PAPER / "exploratory/A272_backbone_pairwise/scripts/score_a272.py")
A272 = importlib.util.module_from_spec(spec)
sys.modules["a274_a272"] = A272
spec.loader.exec_module(A272)
A150 = A272.A150
import cloudpickle  # noqa: E402

for _m in (A272, A150):  # file-loaded modules: ship by value to joblib workers
    cloudpickle.register_pickle_by_value(_m)

AMENDMENT = "A274_celltype_map"
SEED = 274
BODY = {"organ_optic_nerve": "organ_optic_nerve_genome", "organ_eye_anterior": "organ_eye_anterior_genome",
        "kidney_sc": "kidney_sc_genome", "islets": "islets_genome", "kong_ileum": "kong_ileum_genome",
        "kong_colon": "kong_colon_genome", "lung_natri_control": "lung_natri_control_genome"}
CORTEX = {"psychad_hbcc": "psychad_hbcc_genome", "seaad_dlpfc": "seaad_dlpfc_genome", "seaad_mtg": "seaad_mtg_genome"}
ALL = {**BODY, **CORTEX}
KINDS = A272.KINDS
META = A272.META
MIN_PAIRS = 3
log = A272.log


def pseudo_sets(A, k, n, kind):
    pb, genes = A["pb"], A["genes"]
    X = pb[genes].to_numpy(float)
    frac = (X / pb.depth.to_numpy(float)[:, None]).mean(0)
    keep = (X.mean(0) >= 1) | np.isin(genes, A272.CLOCK)
    lm = pd.Series(np.log(frac + 1e-12), index=genes)[keep]
    return P.pseudo_waveforms(lm, n, np.random.default_rng([SEED, k, KINDS[kind]]), cs.gtex(), kind)


def _chunk(pb, wfs, pairs):
    return A272._chunk(pb, wfs, pairs)


def t_class(absr, meta, cls, atlases):
    v = [np.nanmedian(absr[((meta.atlas == n) & (meta.klass == cls)).to_numpy()]) for n in atlases]
    return float(np.mean(v)) if len(v) else np.nan


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
    A150.ATLASES.update(ALL)
    rng = np.random.default_rng([SEED, 99])

    atlas, pairs, meta, twin = {}, {}, [], []
    for n in ALL:
        if args.only is not None and n not in args.only:
            continue
        try:
            A = A150.load_atlas(n, None, 0)
        except Exception as e:  # noqa: BLE001  (not built, or rows missing: left out and named)
            twin.append(dict(atlas=n, eligible=False, reason=f"{type(e).__name__}: {str(e)[:150]}"))
            log(n, "not eligible:", twin[-1]["reason"])
            continue
        rep = A["rep"]
        ok = rep["rows_with_count_diff"] == 0
        twin.append(dict(atlas=n, eligible=ok, reason="" if ok else "clock counts differ from the clock-gene table", **rep,
                         donors=A["pb"].donor.nunique(), compartments=A["pb"].compartment.nunique()))
        if not ok:
            continue
        A["pb"] = A["pb"].reset_index(drop=True)
        if args.placebo:
            A["pb"] = A272.placebo(A["pb"], rng)
        atlas[n], pairs[n] = A, A272.pair_index(A["pb"])
        cortex = n in CORTEX
        for pr in pairs[n]:
            fa, fb = A150.FAM.family(pr["a"]), A150.FAM.family(pr["b"])
            kl = ("CX_diff" if fa != fb else "CX_same") if cortex else A272.klass(A272.group(fa), A272.group(fb))
            meta.append(dict(atlas=n, a=pr["a"], b=pr["b"], family_a=fa, family_b=fb, klass=kl, n=pr["n"], df=pr["df"],
                             n_strata=pr["n_strata"]))
        log(n, "rows", len(A["pb"]), "donors", A["pb"].donor.nunique(), "pairs", len(pairs[n]))
    pd.DataFrame(twin).to_csv(out / "twin.csv", index=False)
    meta = pd.DataFrame(meta)
    design = (meta.groupby(["atlas", "klass"]).agg(pairs=("n", "size"), n_median=("n", "median"), df_min=("df", "min"))
              .reset_index())
    design.to_csv(out / "design.csv", index=False)
    print(pd.DataFrame(twin).to_string())
    print(design.to_string())
    fams = meta.groupby("atlas").apply(lambda g: sorted(set(g.family_a) | set(g.family_b)), include_groups=False)
    print(fams.to_string())
    if args.stage == "design":
        return

    names = list(atlas)
    cnt = design.pivot(index="atlas", columns="klass", values="pairs").reindex(names).fillna(0)
    body = [n for n in names if n in BODY]
    sv = [n for n in body if cnt.get("SV_X", pd.Series(0, index=names))[n] >= MIN_PAIRS]
    both = [n for n in sv if cnt.get("LY_X", pd.Series(0, index=names))[n] >= MIN_PAIRS]
    cx = [n for n in names if n in CORTEX and cnt.get("CX_diff", pd.Series(0, index=names))[n] >= MIN_PAIRS]
    log("H274.1 atlases", both, "| H274.2", sv, "| H274.3", cx)

    r_clock = np.concatenate([A272.all_pairs(atlas[n]["pb"], cs.gtex(), pairs[n]) for n in names])
    R = {}
    for kind in KINDS:
        tasks = []
        for n in names:
            wfs = pseudo_sets(atlas[n], list(ALL).index(n), args.n_pseudo, kind)
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

    absc = np.abs(r_clock)
    pt = meta.copy()
    pt["r_abs"], pt["offset_h"] = absc, np.angle(r_clock) * HOURS
    for kind in KINDS:
        a = np.abs(R[kind])
        pt[f"{kind}_median"], pt[f"{kind}_q95"] = np.nanmedian(a, 0), np.nanquantile(a, 0.95, 0)
    pt["above_both"] = (pt.r_abs > pt.matched_q95) & (pt.r_abs > pt.responsive_q95)

    def stat(a, which):
        if which == "H274.1 D":
            return t_class(a, meta, "SV_X", both) - t_class(a, meta, "LY_X", both)
        if which == "H274.2 T_SV":
            return t_class(a, meta, "SV_X", sv)
        return t_class(a, meta, "CX_diff", cx)

    tests, sets = [], []
    for h, used in (("H274.1 D", both), ("H274.2 T_SV", sv), ("H274.3 T_CX", cx)):
        for kind in KINDS:
            if not used:
                tests.append(dict(test=h, family=kind, atlases="", note="no eligible atlas"))
                continue
            a = np.abs(R[kind])
            v = np.array([stat(a[j], h) for j in range(len(a))])
            c = stat(absc, h)
            ok = np.isfinite(v)
            tests.append(dict(test=h, family=kind, atlases=";".join(used), clock=c, pseudo_median=float(np.nanmedian(v)),
                              pseudo_q95=float(np.nanquantile(v, 0.95)), n_sets=int(ok.sum()),
                              p=float((1 + np.sum(v[ok] >= c)) / (1 + ok.sum())), above_q95=bool(c > np.nanquantile(v, 0.95))))
            sets.append(pd.DataFrame(dict(test=h, family=kind, set=range(len(a)), value=v)))
    tests = pd.DataFrame(tests)
    tests.to_csv(out / "tests.csv", index=False)
    pd.concat(sets).to_csv(out / "sets.csv", index=False)

    # per-atlas statistics (descriptive)
    rows = []
    for n in names:
        for cls in sorted(meta.klass[meta.atlas == n].unique()):
            m = ((meta.atlas == n) & (meta.klass == cls)).to_numpy()
            row = dict(atlas=n, klass=cls, pairs=int(m.sum()), clock=float(np.nanmedian(absc[m])),
                       share_above_both=float(pt.above_both[m].mean()))
            for kind in KINDS:
                v = np.nanmedian(np.abs(R[kind])[:, m], 1)
                row[f"{kind}_median"], row[f"{kind}_q95"] = float(np.nanmedian(v)), float(np.nanquantile(v, 0.95))
            rows.append(row)
    pd.DataFrame(rows).to_csv(out / "per_atlas.csv", index=False)
    fam = pt.assign(f1=np.minimum(pt.family_a, pt.family_b), f2=np.maximum(pt.family_a, pt.family_b),
                    part=np.where(pt.atlas.isin(list(CORTEX)), "cortex", "body"))
    fm = fam.groupby(["part", "f1", "f2"]).agg(pairs=("r_abs", "size"), atlases=("atlas", "nunique"), clock=("r_abs", "median"),
                                               matched=("matched_median", "median"), responsive=("responsive_median", "median"),
                                               share_above_both=("above_both", "mean")).reset_index()
    fm.to_csv(out / "family_matrix.csv", index=False)

    # offsets of the pairs above both companions (donor bootstrap within stratum)
    brng = np.random.default_rng([SEED, 7])
    lo, hi = np.full(len(pt), np.nan), np.full(len(pt), np.nan)
    off = 0
    for n in names:
        z, v = A272.atlas_states(atlas[n]["pb"], cs.gtex())
        for k, pr in enumerate(pairs[n]):
            i = off + k
            if not pt.above_both.iloc[i]:
                continue
            ang = []
            for _ in range(args.n_boot):
                idx = np.concatenate([brng.choice(np.where(pr["s"] == s)[0], (pr["s"] == s).sum()) for s in range(pr["s"].max() + 1)])
                r = A272.pair_r(z, v, dict(ia=pr["ia"][idx], ib=pr["ib"][idx], s=pr["s"][idx]))
                ang.append(np.angle(r * np.conj(r_clock[i])) * HOURS)
            q = np.nanquantile(ang, [0.05, 0.95])
            lo[i], hi[i] = pt.offset_h.iloc[i] + q[0], pt.offset_h.iloc[i] + q[1]
        off += len(pairs[n])
    pt["offset_lo90"], pt["offset_hi90"] = lo, hi
    pt.to_csv(out / "pairs.csv", index=False)
    log("done")
    warnings.filterwarnings("ignore")
    print(tests.to_string())
    print(pd.DataFrame(rows).round(3).to_string())


if __name__ == "__main__":
    main()
