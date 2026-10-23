"""A31: untimed per-cell-type gene phases against true times (mouse colon chrono-atlas, known ZT).

Pre-registered (A31, written before running; design disclosed in its run note). Mouse colon chrono-atlas (Carmona-Alcocer 2026):
16 mice sampled at ZT 1, 7, 13, 19 (4 per ZT), arms C-DC and R-DC pooled (stratum = arm), classes with >= 30 cells per
pseudobulk and >= 3 animals at every ZT (12 classes). Per class: the host of the three stromal anchors (fibroblast,
smooth_muscle, endothelial; leave-one-out for an anchor), untimed loadings of every gene on the host state
(`Dataset.gene_rhythms`) against the timed cosinor on ZT (arm intercepts, 2-df F). The untimed loadings are put in the ZT frame
by the circular mean over the 14 clock genes of their phase difference to the timed ones (one rotation per class). Scores per
set (a1 followers timed-rhythmic p < 0.05, a2 all expressed followers, b all rhythmic non-clock genes p < 0.01, clock14): R,
offset, A = mean cos(d), within 3 h, rise / fall at ZT {1, 7, 13, 19}, amplitude Spearman. Null: host states permuted within
arm (`n_null`), all rescored. Primary: pooled R of the 31 mouse orthologs of the 32 D10 followers (a1) above the null 95 % quantile.
Feasibility: >= 3 classes. Sensitivity: arms C-DC + R-PC. Exploratory E1 (post hoc): the host refitted on animals of
consecutive ZTs (3 or 2 ZTs), scored against the full-design timed truth.

    python compute_celltype_truth.py [--subsample N] [--arms C-DC R-DC] [--tag primary] [--no-window]
Outputs: paper/results/figure_2/output/gene_rhythms/a31_{scores,hosts,window_E1}_<tag>.csv, a31_genes_<tag>.csv.gz
Legacy: celltype_truth/scripts/{build,score,window}.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as sps  # noqa: E402

import clocksync as cs  # noqa: E402
import _shared as S  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "gene_rhythms"
W = S.W
MIN_CELLS, MIN_PER_ZT, MIN_MEAN = 30, 3, 5.0
ANCHORS = ["fibroblast", "smooth_muscle", "endothelial"]
ZTS = np.array([1.0, 7.0, 13.0, 19.0])


def prepare_table(pb, donors, arms, quick_donors=None):
    pb = pb[pb.donor.isin(donors.index[donors.arm.isin(arms)]) & (pb.n_cells >= MIN_CELLS)].copy()
    pb["zt"] = pb.donor.map(donors.zt)
    pb["arm"] = pb.donor.map(donors.arm)
    pb["stratum"] = pb.arm
    per_zt = pb.groupby(["compartment", "zt"]).donor.nunique().unstack().reindex(columns=ZTS).fillna(0)
    kept = sorted(per_zt.index[(per_zt >= MIN_PER_ZT).all(axis=1)])
    return pb[pb.compartment.isin(kept)].reset_index(drop=True), per_zt, kept


def timed_fit(q, Y):
    """Cosinor on ZT with arm intercepts, 2-df F p (local: not in clocksync)."""
    Dm = pd.get_dummies(q.arm).to_numpy(float)
    t = q.zt.to_numpy() * W
    D1 = np.c_[Dm, np.cos(t), np.sin(t)]
    B = np.linalg.lstsq(D1, Y, rcond=None)[0]
    rss1 = ((Y - D1 @ B) ** 2).sum(0)
    rss0 = ((Y - Dm @ np.linalg.lstsq(Dm, Y, rcond=None)[0]) ** 2).sum(0)
    df2 = len(q) - Dm.shape[1] - 2
    return B[-2] + 1j * B[-1], sps.f.sf(((rss0 - rss1) / 2) / (rss1 / df2), 2, df2)


def prepare_truth(pb, kept, rh_main):
    """Per class: genes (those of the untimed fit), timed loading and p, gene sets."""
    P = {}
    for k in kept:
        g = rh_main[rh_main.compartment == k].gene.to_numpy()
        q = pb[pb.compartment == k]
        C = q[g].to_numpy(float) if set(g) <= set(q.columns) else None
        Y = np.log((C + 0.5) / q.depth.to_numpy(float)[:, None])
        ct, pt = timed_fit(q, Y)
        clock = np.isin(g, S.CLOCK)
        fol = np.isin(g, S.FOLLOWERS)
        P[k] = dict(genes=g, ct=ct, pt=pt, clock=clock,
                    sets={"a1_followers_rhythmic": fol & (pt < 0.05), "a2_followers_all": fol, "b_rhythmic_p01": ~clock & (pt < 0.01)})
    return P


def loadings(P, rh, null_id):
    """Untimed loadings in the ZT frame per class (rotation from the 14 clock genes)."""
    r = rh[rh["null"] == null_id]
    out = {}
    for k, p in P.items():
        x = r[r.compartment == k].drop_duplicates("gene").set_index("gene").reindex(p["genes"])
        cu = (x.c_re + 1j * x.c_im).to_numpy()
        cl = p["clock"] & np.isfinite(cu)
        rot = np.angle(np.exp(1j * (np.angle(cu[cl]) - np.angle(p["ct"][cl]))).mean())
        out[k] = (cu * np.exp(-1j * rot), rot)
    return out


def scores(cu, ct, scale=None):
    ok = np.isfinite(cu)
    cu, ct = cu[ok], ct[ok]
    if len(ct) == 0:
        return dict(n=0)
    d = np.angle(cu) - np.angle(ct)
    m = np.exp(1j * d).mean()
    a, b = np.angle(cu)[:, None], np.angle(ct)[:, None]
    rf = (np.sign(np.sin(a - W * ZTS)) == np.sign(np.sin(b - W * ZTS))).mean()
    amp = np.abs(cu) * (1.0 if scale is None else scale[ok])
    return dict(n=len(ct), R=abs(m), offset_h=float(np.angle(m) / W), A=float(np.cos(d).mean()),
                within3h=float((np.abs(np.angle(np.exp(1j * d))) <= 3 * W).mean()), risefall=float(rf),
                amp_spearman=sps.spearmanr(amp, abs(ct))[0] if len(ct) > 2 else np.nan)


def evaluate(P, L):
    rows, pooled = [], {}
    for k, p in P.items():
        cu, rot = L[k]
        cl = p["clock"] & np.isfinite(cu)
        scale = np.median(np.abs(p["ct"][cl]) / np.abs(cu[cl])) if cl.any() else np.nan
        for s, sel in p["sets"].items():
            rows.append(dict(compartment=k, set=s, rot_h=rot / W, **scores(cu[sel], p["ct"][sel])))
            acc = pooled.setdefault(s, ([], [], []))
            acc[0].append(cu[sel]); acc[1].append(p["ct"][sel]); acc[2].append(np.full(sel.sum(), scale))
        rows.append(dict(compartment=k, set="clock14", rot_h=rot / W, **scores(cu[p["clock"]], p["ct"][p["clock"]])))
    for s, (a, b, sc) in pooled.items():
        rows.append(dict(compartment="POOLED", set=s, **scores(np.concatenate(a), np.concatenate(b), np.concatenate(sc))))
    return pd.DataFrame(rows)


def fit_hosts(pb, kept, genes, n_null, seed):
    clock = pb[list(S.META) + [c for c in S.CLOCK if c in pb.columns]]
    ds = cs.Dataset(clock, ANCHORS, name="colon_chrono")
    rh = ds.gene_rhythms(pb[list(S.META) + genes], genes=genes, min_mean=MIN_MEAN, n_null=n_null, seed=seed)
    zt = pb.groupby("donor").zt.first()
    return ds, rh, zt


def compare_null(obs, null, n_null):
    key = ["compartment", "set"]
    nq = null.groupby(key).agg(R_null95=("R", lambda x: np.nanquantile(x, .95)), R_null_med=("R", "median"),
                              A_null95=("A", lambda x: np.nanquantile(x, .95)),
                              rf_null95=("risefall", lambda x: np.nanquantile(x, .95))).reset_index()
    res = obs.merge(nq, on=key, how="left")
    nm = null.set_index(key)
    for col in ("R", "A", "risefall"):
        res[f"p_{col}"] = [(1 + (nm.loc[(r.compartment, r.set), col].to_numpy() >= getattr(r, col) - 1e-12).sum()) / (1 + n_null)
                           if np.isfinite(getattr(r, col, np.nan)) else np.nan for r in res.itertuples()]
    return res


def run_design(pb, kept, genes, n_null, seed, P=None):
    ds, rh, zt = fit_hosts(pb, kept, genes, n_null, seed)
    main = rh[rh["null"] == -1]
    P = P or prepare_truth(pb, kept, main)
    obs = evaluate(P, loadings(P, rh, -1))
    null = pd.concat([evaluate(P, loadings(P, rh, b)).assign(perm=b) for b in range(n_null)])
    return ds, rh, P, compare_null(obs, null, n_null), zt


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--arms", nargs="+", default=["C-DC", "R-DC"])
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--no-window", action="store_true")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_null = run["n_null"]
    pb_all, donors = S.genome_table("colon_chrono_genome")
    pb, per_zt, kept = prepare_table(pb_all, donors, args.arms)
    print("animals per ZT:\n", per_zt.astype(int).to_string(), "\nkept:", kept)
    if len(kept) < 3:
        print("INFEASIBLE: fewer than 3 classes with >= 3 animals per ZT")
        return
    genes = [c for c in pb.columns if c not in list(S.META) + ["zt", "arm"]]
    ds, rh, P, res, zt = run_design(pb, kept, genes, n_null, 31)
    res.insert(0, "data", args.tag)
    res.to_csv(out / f"a31_scores_{args.tag}.csv", index=False)
    arm = pb.groupby("donor").arm.first()
    diag = pd.DataFrame([dict(compartment="host", **ds.timed(zt, n_perm=1000 if run["quick"] else 10000).to_dict())])
    diag.to_csv(out / f"a31_hosts_{args.tag}.csv", index=False)
    pd.concat([pd.DataFrame(dict(compartment=k, gene=p["genes"], c_timed=p["ct"], p_timed=p["pt"], clock=p["clock"]))
               for k, p in P.items()]).to_csv(out / f"a31_genes_{args.tag}.csv.gz", index=False)
    pd.set_option("display.width", 250)
    pr = res[(res.compartment == "POOLED") & (res.set == "a1_followers_rhythmic")].iloc[0]
    print(res[res.compartment == "POOLED"].round(3).to_string(index=False))
    print(f"PRIMARY ({args.tag}): pooled R(a1) = {pr.R:.3f}, null95 = {pr.R_null95:.3f}, p = {pr.p_R:.3f}, n = {pr.n:.0f}, "
          f"offset = {pr.offset_h:.1f} h, A = {pr.A:.2f} -> {'PASS' if pr.R > pr.R_null95 else 'FAIL'}")
    print(diag.round(3).to_string(index=False))
    if args.no_window:
        return
    rows = []
    zt_a = pb.groupby("donor").zt.first()
    for m in (3, 2):
        for i in range(4):
            pick = [ZTS[(i + j) % 4] for j in range(m)]
            keep = zt_a.index[zt_a.isin(pick)]
            sub = pb[pb.donor.isin(keep)].reset_index(drop=True)
            _, _, _, rw, _ = run_design(sub, kept, genes, n_null, 100 * m + i, P=P)
            rw = rw[rw.compartment == "POOLED"].assign(zts="/".join(f"{p:g}" for p in pick), n_animals=len(keep))
            rows.append(rw)
            print(rw[["set", "R", "R_null95", "offset_h", "A", "zts", "n_animals"]].round(3).to_string(index=False), flush=True)
    w = pd.concat(rows)
    w.to_csv(out / f"a31_window_E1_{args.tag}.csv", index=False)
    print(w.groupby(["n_animals", "set"])[["R", "R_null95", "offset_h", "A", "within3h"]].median().round(3).to_string())


if __name__ == "__main__":
    main()
