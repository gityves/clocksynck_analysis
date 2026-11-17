"""A284: does sharing of the tissue's clock state increase with residency? Resident against circulating T cells, plasma
against B cells, macrophages against monocytes.

Registered in paper/docs/registered_experiments.md (A284, with its design addendum). Per atlas and contrast, donors
with both classes: the deeper of the two immune rows of a donor thinned binomially (all genes) to the other's
clock-gene count (10 replicates); T_res / T_circ = median clock |r| (`clocksync.pairwise`, A272 statistic) of each class
with the tissue classes X (non-immune families); Delta = T_res - T_circ, averaged over replicates; pooled = mean over
atlases. Null: the two rows of each donor swapped at random (2000 draws). Companions: 200 matched + 200 responsive
sets on the same thinned tables, the clock's pi against the swap distribution of pi (A283 form).

    python score_a284.py [--design] [--draws 2000] [--jobs 6]
Outputs (../output/): design.csv, per_atlas.csv, tests.csv
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import build_a284 as B  # noqa: E402
from clocksync import pairwise as pw  # noqa: E402
from common import cli  # noqa: E402

SEED, N_REP, N_SETS, MIN_DONORS = 284, 10, 200, 15
IMMUNE = ("lymphocyte", "myeloid")
CONTRASTS = {  # hypothesis: {atlas: (resident compartment, circulating compartment)}
    "H284.1 T": {"ibd_ileum": ("ileum:T_resident", "ileum:T_circulating"),
                 "ibd_rectum": ("rectum:T_resident", "rectum:T_circulating"),
                 "organ_breast": ("breast:T_resident", "breast:T_circulating"),
                 "kong_ileum": ("ileum:T_resident", "ileum:T_circulating")},
    "H284.2 plasma-B": {"organ_breast": ("breast:plasma", "breast:B"), "organ_colon": ("colon:plasma", "colon:B"),
                        "kong_ileum": ("ileum:plasma", "ileum:B"), "kong_colon": ("colon:plasma", "colon:B"),
                        "ibd_ileum": ("ileum:plasma", "ileum:B"), "ibd_rectum": ("rectum:plasma", "rectum:B")},
    "H284.3 macrophage-monocyte": {"lung_hlca": ("lung:macrophage", "lung:monocyte"),
                                   "organ_breast": ("breast:macrophage_DC", "breast:monocyte")},
}
ATLAS_GENOME = {"lung_hlca": "lung_hlca_genome", "organ_breast": "organ_breast_genome",
                "organ_colon": "organ_colon_genome", "kong_ileum": "kong_ileum_genome",
                "kong_colon": "kong_colon_genome", "ibd_ileum": "ibd_ileum_genome", "ibd_rectum": "ibd_rectum_genome"}


def _a280():
    spec = importlib.util.spec_from_file_location("a284_a280", PAPER / "exploratory/A280_celltype_replication2/scripts/score_a280.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["a284_a280"] = m
    spec.loader.exec_module(m)
    return m.A274.A150  # A150 with the anchor fallback (IBD targets have no learned anchors)


def atlas_table(n, A150):
    """Healthy rows of the atlas (A17 targets, genome-wide, A272/A274 strata); for the T atlases the A284 table's
    T_resident / T_circulating rows added (and Kong ileum's epithelial rows without the intraepithelial lymphocytes)."""
    A150.ATLASES.update(ATLAS_GENOME)
    base = A150.load_atlas(n, None, 0)["pb"].reset_index(drop=True)
    if n not in CONTRASTS["H284.1 T"]:
        return base, {}
    t, check = B.table(n)
    st = base.drop_duplicates("donor").set_index("donor").stratum
    add = t[t.compartment.str.endswith(("T_resident", "T_circulating")) & t.donor.isin(st.index)].copy()
    if n == "kong_ileum":  # IEL were classed as epithelial in the source table (design addendum)
        base = base[base.compartment != "ileum:epithelial"]
        add = pd.concat([add, t[(t.compartment == "ileum:epithelial") & t.donor.isin(st.index)]])
    add["stratum"] = add.donor.map(st)
    cols = list(base.columns)
    return pd.concat([base, add.reindex(columns=cols).fillna(0.0)], ignore_index=True), check


def unit(n, res, circ, pb, family, k):
    """Thinned replicates of one atlas x contrast: per replicate the |r| of every set on every (class, X) pair."""
    from clocksync.config import limit_threads
    limit_threads(1)
    warnings.filterwarnings("ignore")
    fam = {c: family(c) for c in pb.compartment.unique()}
    X = [c for c in fam if fam[c] not in IMMUNE and c not in (res, circ)]
    both = sorted(set(pb.donor[pb.compartment == res]) & set(pb.donor[pb.compartment == circ]))
    t = pb[pb.donor.isin(both) & pb.compartment.isin([res, circ] + X)].reset_index(drop=True)
    meta = ["donor", "compartment", "stratum", "n_cells", "depth"]
    genes = [c for c in t.columns if c not in meta]
    import clocksync as cs
    clock = list(cs.gtex().genes)
    wfs = [cs.gtex()]
    for j, kind in enumerate(("matched", "responsive")):
        wfs += pw.companions(t[meta + genes], N_SETS, kind, np.random.default_rng([SEED, k, j]))
    use = sorted({g for w in wfs for g in w.genes})
    t = t[meta + use]
    ir = t.index[t.compartment == res].to_numpy()
    ic = t.index[t.compartment == circ].to_numpy()
    r_d, c_d = t.donor[ir].to_numpy(), t.donor[ic].to_numpy()
    ic = ic[pd.Index(c_d).get_indexer(r_d)]
    reps = []
    for rep in range(N_REP):
        rng = np.random.default_rng([SEED, k, 100 + rep])
        tt = t.copy()
        cnt = tt[use].to_numpy(np.int64)
        s = tt[clock].sum(1).to_numpy(float)
        for a, b in zip(ir, ic):  # thin the deeper row of the pair to the other's clock-gene count
            hi, lo = (a, b) if s[a] > s[b] else (b, a)
            p = s[lo] / max(s[hi], 1.0)
            cnt[hi] = rng.binomial(cnt[hi], p)
            tt.loc[hi, "depth"] = tt.depth[hi] * p
        tt[use] = cnt
        Z = np.array([pw.states(tt, w) for w in wfs])
        reps.append(Z)
    return dict(atlas=n, res=res, circ=circ, X=X, key=t[["donor", "compartment", "stratum"]], ir=ir, ic=ic, Z=reps,
                donors=len(both))


def delta(u, swap=None):
    """Delta (sets) averaged over replicates; `swap`: boolean per donor pair, resident and circulating rows exchanged."""
    key = u["key"].copy()
    if swap is not None:
        comp = key.compartment.to_numpy().copy()
        comp[u["ir"][swap]], comp[u["ic"][swap]] = u["circ"], u["res"]
        key["compartment"] = comp
    pairs = pw.pair_index(key)
    sel_r = [i for i, p in enumerate(pairs) if u["res"] in (p["a"], p["b"])]
    sel_c = [i for i, p in enumerate(pairs) if u["circ"] in (p["a"], p["b"])]
    sel_r = [i for i in sel_r if u["circ"] not in (pairs[i]["a"], pairs[i]["b"])]
    sel_c = [i for i in sel_c if u["res"] not in (pairs[i]["a"], pairs[i]["b"])]
    out = []
    for Z in u["Z"]:
        A = pw.abs_r(Z, pairs)
        out.append(np.nanmedian(A[:, sel_r], 1) - np.nanmedian(A[:, sel_c], 1))
    return np.mean(out, 0), len(sel_r), len(sel_c)


def draws(us, ds):
    """Pooled Delta under within-donor swaps for the draws `ds` (seeds [284, 9, d, atlas])."""
    from clocksync.config import limit_threads
    limit_threads(1)
    warnings.filterwarnings("ignore")
    return np.array([np.mean([delta(u, np.random.default_rng([SEED, 9, d, i]).random(len(u["ir"])) < 0.5)[0]
                              for i, u in enumerate(us)], 0) for d in ds])


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--design", action="store_true")
    ap.add_argument("--draws", type=int, default=2000)
    args = ap.parse_args()
    cli.setup(args)
    out = HERE.parents[1] / "output"
    out.mkdir(parents=True, exist_ok=True)
    A150 = _a280()
    family = A150.FAM.family
    atlases = sorted({a for c in CONTRASTS.values() for a in c})
    tabs, checks = {}, []
    for n in atlases:
        tabs[n], chk = atlas_table(n, A150)
        if chk:
            checks.append(chk)
    design = []
    for h, per in CONTRASTS.items():
        for n, (res, circ) in per.items():
            pb = tabs[n]
            both = set(pb.donor[pb.compartment == res]) & set(pb.donor[pb.compartment == circ])
            design.append(dict(hypothesis=h, atlas=n, resident=res, circulating=circ, donors_both=len(both),
                               eligible=len(both) >= MIN_DONORS))
    D = pd.DataFrame(design)
    D.to_csv(out / "design.csv", index=False)
    pd.DataFrame(checks).to_csv(out / "table_checks.csv", index=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(checks).to_string(), D.to_string(), sep="\n", flush=True)
    if args.design:
        return

    jobs = [(h, r.atlas, r.resident, r.circulating) for h in CONTRASTS for r in D[(D.hypothesis == h) & D.eligible].itertuples()]
    units = Parallel(n_jobs=args.jobs)(delayed(unit)(n, res, circ, tabs[n], family, k)
                                       for k, (h, n, res, circ) in enumerate(jobs))
    rows, tests = [], []
    fam = {"matched": slice(1, 1 + N_SETS), "responsive": slice(1 + N_SETS, 1 + 2 * N_SETS)}
    for h in CONTRASTS:
        us = [u for (hh, *_), u in zip(jobs, units) if hh == h]
        if not us:
            continue
        obs = []
        for u in us:
            d, nr, nc = delta(u)
            obs.append(d)
            rows.append(dict(hypothesis=h, atlas=u["atlas"], donors=u["donors"], pairs_res=nr, pairs_circ=nc,
                             delta_clock=float(d[0]), **{f"delta_{f}_median": float(np.nanmedian(d[s])) for f, s in fam.items()}))
        pooled = np.mean(obs, 0)

        chunks = np.array_split(np.arange(args.draws), args.jobs)
        null = np.concatenate(Parallel(n_jobs=args.jobs)(delayed(draws)(us, c) for c in chunks))
        row = dict(hypothesis=h, atlases=len(us), delta=float(pooled[0]),
                   p_swap=float((1 + np.sum(null[:, 0] >= pooled[0])) / (1 + len(null))))
        for f, s in fam.items():
            o = pw.percentile(pooled, s)[0]
            row.update({f"{f}_median": float(np.nanmedian(pooled[s])), f"{f}_pi": float(o),
                        f"{f}_p_cal": float(pw.calibrated_p(o, pw.percentile(null, s))[0])})
        tests.append(row)
        print(h, "done", flush=True)
    T = pd.DataFrame(tests)
    order = np.argsort(T.p_swap.to_numpy())
    holm, run = np.empty(len(T)), 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(T) - r) * T.p_swap[i]))
        holm[i] = run
    T["p_holm"] = holm
    T["passes"] = (T.delta > 0) & (T.p_holm < 0.05)
    T.to_csv(out / "tests.csv", index=False)
    pd.DataFrame(rows).to_csv(out / "per_atlas.csv", index=False)
    print(pd.DataFrame(rows).round(4).to_string(), T.round(4).to_string(), sep="\n")


if __name__ == "__main__":
    main()
