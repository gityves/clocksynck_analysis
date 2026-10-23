"""A31: untimed per-cell-type gene phases against true times (mouse colon chrono-atlas, known ZT).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A31; run note and the post hoc window check E1). Hypotheses: A31 primary
(pooled over cell types, the R of the frozen follower list that is timed-rhythmic in its class exceeds its host-permutation
null 95 % quantile), also reported: all expressed followers, all rhythmic non-clock genes (timed p < 0.01).
Estimand: R = |mean exp(i (arg c_untimed - arg c_timed))| over (gene, class) pairs, after one rotation per class fixed on
the 14 clock genes (the untimed frame is arbitrary); offset (h), A = mean cos, share within 3 h, rise/fall agreement at the
four ZTs. Untimed: `Dataset.gene_rhythms` (gene ~ host state of the other classes, anchors fibroblast, smooth muscle,
endothelium; an anchor on the host of the other two). Timed: gene ~ stratum + cos/sin(ZT) (F test, local).
Null: animals permuted over the host (within arm), 200 permutations, p = (1 + #null >= obs) / (1 + n).
`--window`: E1 (post hoc, exploratory): host refitted on animals of 3 or 2 consecutive ZTs, timed truth from the full design.

    python compute_celltype_truth.py [--subsample N] [--arms C-DC,R-DC] [--tag primary] [--window]
Outputs: paper/results/figure_1/output/method/celltype_truth_{scores,hosts,genes}_<tag>.csv[.gz], celltype_truth_window_E1.csv
Legacy: legacy/desync/validation/celltype_truth/scripts/{score,window}.py (output/scores_primary.csv, window_E1.csv)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "method"
W = 2 * np.pi / 24
MIN_CELLS, MIN_PER_ZT, MIN_MEAN = 30, 3, 5.0
ANCHORS = ["fibroblast", "smooth_muscle", "endothelial"]
CLOCK = list(cs.gtex().genes)
FOLLOWERS = ["ADGRA2", "ADIRF-AS1", "ARHGAP42", "ARNTL2", "BHLHE40", "BHLHE41", "CDC25B", "CDK18", "CHST7", "CIPC",
             "DDIT4L", "EPB41L1", "FAM76A", "FOXO1", "GUCY1A1", "HLF", "JAM2", "KLF10", "LDB1", "LONRF3", "LPCAT3",
             "MAMDC2", "MPP2", "NAMPT", "NFIL3", "PPARA", "RAB3IL1", "RORC", "SUN2", "TSPAN4", "USP2", "WEE1"]
ZTS = np.array([1.0, 7.0, 13.0, 19.0])


def ols(D, Y):
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    return B, ((Y - D @ B) ** 2).sum(0)


def timed_cosinor(q: pd.DataFrame, Y):
    """local: not in clocksync. Gene ~ arm + cos/sin(ZT): complex coefficient and F-test p (2 df) per gene."""
    Dm = pd.get_dummies(q.arm).to_numpy(float)
    t = q.zt.to_numpy() * W
    B, rss1 = ols(np.c_[Dm, np.cos(t), np.sin(t)], Y)
    _, rss0 = ols(Dm, Y)
    df2 = len(q) - Dm.shape[1] - 2
    F = ((rss0 - rss1) / 2) / (rss1 / df2)
    return B[-2] + 1j * B[-1], stats.f.sf(F, 2, df2)


def prepare(pb, donors, kept, genes):
    """Per class: rows, expressed genes, timed fit (full design) and gene sets."""
    P = {}
    for k in kept:
        q = pb[pb.compartment == k].reset_index(drop=True)
        q["arm"] = donors.arm.reindex(q.donor).to_numpy()
        q["zt"] = donors.zt.reindex(q.donor).to_numpy()
        C = q[genes].to_numpy(float)
        ok = C.mean(0) >= MIN_MEAN
        g = np.asarray(genes)[ok]
        Y = np.log((C[:, ok] + 0.5) / q.depth.to_numpy(float)[:, None])
        ct, pt = timed_cosinor(q, Y)
        clock, fol = np.isin(g, CLOCK), np.isin(g, FOLLOWERS)
        P[k] = dict(genes=g, ct=ct, pt=pt, clock=clock,
                    sets={"a1_followers_rhythmic": fol & (pt < 0.05), "a2_followers_all": fol,
                          "b_rhythmic_p01": ~clock & (pt < 0.01)})
    return P


def scores(cu, ct):
    if len(ct) == 0:
        return dict(n=0)
    d = np.angle(cu) - np.angle(ct)
    m = np.exp(1j * d).mean()
    a, b = np.angle(cu)[:, None], np.angle(ct)[:, None]
    rf = (np.sign(np.sin(a - W * ZTS)) == np.sign(np.sin(b - W * ZTS))).mean()
    return dict(n=len(ct), R=abs(m), offset_h=float(np.angle(m) / W), A=float(np.cos(d).mean()),
                within3h=float((np.abs(np.angle(np.exp(1j * d))) <= 3 * W).mean()), risefall=float(rf),
                amp_spearman=stats.spearmanr(np.abs(cu), np.abs(ct))[0] if len(ct) > 2 else np.nan)


def untimed(ds, genome, genes):
    """{class: (complex untimed loading per expressed gene in the frame of the timed fit, rotation)} from
    `Dataset.gene_rhythms`; genes are matched to the timed fit's gene list."""
    r = ds.gene_rhythms(genome, genes=genes, min_mean=MIN_MEAN)
    return {k: g.set_index("gene") for k, g in r.groupby("compartment")}


def loadings(P, R):
    out = {}
    for k, p in P.items():
        g = R[k].reindex(p["genes"])
        cu = (g.c_re + 1j * g.c_im).to_numpy()
        cl = p["clock"] & np.isfinite(cu)
        rot = np.angle(np.exp(1j * (np.angle(cu[cl]) - np.angle(p["ct"][cl]))).mean())
        out[k] = (cu * np.exp(-1j * rot), rot)
    return out


def evaluate(P, L):
    rows, pooled = [], {}
    for k, p in P.items():
        cu, rot = L[k]
        for s, sel in list(p["sets"].items()) + [("clock14", p["clock"])]:
            sel = sel & np.isfinite(cu)
            rows.append(dict(compartment=k, set=s, rot_h=rot / W, **scores(cu[sel], p["ct"][sel])))
            if s != "clock14":
                acc = pooled.setdefault(s, ([], []))
                acc[0].append(cu[sel]); acc[1].append(p["ct"][sel])
    for s, (a, b) in pooled.items():
        rows.append(dict(compartment="POOLED", set=s, **scores(np.concatenate(a), np.concatenate(b))))
    return pd.DataFrame(rows)


def build_dataset(pb, donors, kept):
    clock = pb[pb.compartment.isin(kept)][["donor", "compartment", "depth", *CLOCK]].copy()
    clock["stratum"] = donors.arm.reindex(clock.donor).to_numpy()
    return cs.Dataset(clock, anchors=ANCHORS, donors=donors)


def permuted_genome(pb, donors, rng):
    """Animals permuted over the host within arm: the genome rows keep their own arm as stratum."""
    g = pb.copy()
    g["stratum"] = donors.arm.reindex(g.donor).to_numpy()
    perm = {}
    for _, ids in donors.loc[donors.index.isin(pb.donor)].groupby("arm").groups.items():
        ids = np.array(sorted(ids))
        perm.update(dict(zip(ids, rng.permutation(ids))))
    g["donor"] = g.donor.map(perm)
    return g


def run_null(ds, pb, donors, P, genes, n_perm, seed):
    rng = np.random.default_rng(seed)
    out = []
    for b in range(n_perm):
        g = permuted_genome(pb, donors, rng)
        out.append(evaluate(P, loadings(P, untimed(ds, g, genes))).assign(perm=b))
    return pd.concat(out)


def summarise(obs, null, n_perm):
    key = ["compartment", "set"]
    q = null.groupby(key).agg(R_null95=("R", lambda x: np.nanquantile(x, 0.95)), R_null_med=("R", "median"),
                              A_null95=("A", lambda x: np.nanquantile(x, 0.95)),
                              rf_null95=("risefall", lambda x: np.nanquantile(x, 0.95))).reset_index()
    res = obs.merge(q, on=key, how="left")
    nm = null.set_index(key)
    for col in ("R", "A", "risefall"):
        res[f"p_{col}"] = [(1 + (nm.loc[(r.compartment, r.set), col].to_numpy() >= getattr(r, col) - 1e-12).sum())
                           / (1 + n_perm) if np.isfinite(getattr(r, col, np.nan)) else np.nan
                           for r in res.itertuples()]
    return res


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--arms", default="C-DC,R-DC")
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--window", action="store_true", help="E1: hosts from 2-3 consecutive ZTs (post hoc)")
    args = ap.parse_args()
    opt = cli.setup(args)
    n_perm = args.n_perm or (10 if opt["quick"] else 100 if args.window else 200)
    out = out_dir(FOLDER, opt["quick"])
    min_per_zt = 1 if opt["quick"] else MIN_PER_ZT  # a donor subsample cannot hold 3 animals per ZT

    pb, donors = data.load("colon_chrono_genome", subsample=args.subsample, seed=args.seed)
    arms = args.arms.split(",")
    donors = donors[donors.arm.isin(arms)]
    pb = pb[pb.donor.isin(donors.index) & (pb.n_cells >= MIN_CELLS)].reset_index(drop=True)
    per_zt = (pb.assign(zt=donors.zt.reindex(pb.donor).to_numpy()).groupby(["compartment", "zt"]).donor.nunique()
              .unstack().reindex(columns=ZTS).fillna(0))
    kept = sorted(per_zt.index[(per_zt >= min_per_zt).all(1)])
    print(f"animals per ZT (>= {MIN_CELLS} cells):\n{per_zt.astype(int).to_string()}\nkept: {kept}")
    if len(kept) < 3 or not set(ANCHORS) <= set(kept):
        print("INFEASIBLE: fewer than 3 classes with enough animals per ZT")
        return
    pb = pb[pb.compartment.isin(kept)].reset_index(drop=True)
    genes = [c for c in pb.columns if c not in ("donor", "compartment", "n_cells", "depth")]
    ds = build_dataset(pb, donors, kept)
    P = prepare(pb, donors, kept, genes)

    if args.window:
        rows = []
        zt = donors.zt
        for m in (3, 2):
            for i in range(4):
                pick = [ZTS[(i + j) % 4] for j in range(m)]
                keep = zt.index[zt.isin(pick)]
                sub = pb[pb.donor.isin(keep)].reset_index(drop=True)
                dss = build_dataset(sub, donors, kept)
                Ps = {k: dict(p) for k, p in P.items()}  # timed fit and sets from the full design
                obs = evaluate(Ps, loadings(Ps, untimed(dss, sub.assign(stratum=donors.arm.reindex(sub.donor).to_numpy()), genes)))
                null = run_null(dss, sub, donors, Ps, genes, n_perm, 100 * m + i)
                o = summarise(obs, null, n_perm)
                o = o[o.compartment == "POOLED"].assign(zts="/".join(f"{p:g}" for p in pick), n_animals=len(keep))
                rows.append(o)
                print(o.round(3).to_string(index=False), flush=True)
        pd.concat(rows).to_csv(out / "celltype_truth_window_E1.csv", index=False)
        return

    g0 = pb.assign(stratum=donors.arm.reindex(pb.donor).to_numpy())
    R = untimed(ds, g0, genes)
    L = loadings(P, R)
    obs = evaluate(P, L)
    null = run_null(ds, pb, donors, P, genes, n_perm, 31)
    res = summarise(obs, null, n_perm).assign(data=args.tag)
    res.to_csv(out / f"celltype_truth_scores_{args.tag}.csv", index=False)
    pd.concat([pd.DataFrame(dict(compartment=k, gene=p["genes"], c_timed=p["ct"], p_timed=p["pt"],
                                 c_untimed=L[k][0], clock=p["clock"])) for k, p in P.items()]
              ).to_csv(out / f"celltype_truth_genes_{args.tag}.csv.gz", index=False)
    # host diagnostic: untimed host phase (all anchors) against the recorded ZT, per arm-centred donor
    t = ds.timed(donors.zt, n_perm=1000)
    t.to_frame("value").to_csv(out / f"celltype_truth_hosts_{args.tag}.csv")
    pd.set_option("display.width", 250)
    print(res.drop(columns=["data"]).round(3).to_string(index=False))
    pr = res[(res.compartment == "POOLED") & (res.set == "a1_followers_rhythmic")]
    if len(pr):
        pr = pr.iloc[0]
        print(f"PRIMARY ({args.tag}): pooled R(a1) = {pr.R:.3f}, null95 = {pr.R_null95:.3f}, p = {pr.p_R:.3f}, "
              f"n = {pr.n}, offset = {pr.offset_h:.1f} h, A = {pr.A:.2f} -> {'PASS' if pr.R > pr.R_null95 else 'FAIL'}")
    print("host vs ZT:", t.round(3).to_dict())


if __name__ == "__main__":
    main()
