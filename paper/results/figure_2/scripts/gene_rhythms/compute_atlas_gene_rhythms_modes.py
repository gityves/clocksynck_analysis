"""A27 part 1 and A31b with covariate modes, the composition term and the single-cell z > 2 sets.

A27 part 1 is post hoc and exploratory (no acceptance criterion); the A31b grid is the one written in A31b (A27 pipeline unchanged). The
amendments leave the covariate modes and the composition term to `sc_predict.py`; they are implemented here as there. Modes of the per-class gene
regression on the host state (one intercept per stratum, always): `none`; `basic` = donor age (median fill), male, log depth, log cells of the
pseudobulk; `pcs` = basic + the first 3 principal components of the stratum-centred, standardised log fractions of the class's genes (computed per
class). Composition term: the log fraction of the class among the modelled classes of the donor (log(frac + 1e-4)) is regressed on the same design;
its loading c_comp_k is added to every gene's loading of class k before the share-weighted recombination (`pred_total`; `pred_expr` is the
expression term alone). Single-cell z = |sum_k s_gk c_gk| / sqrt(sum_k (s_gk se_gk)^2); the set "p<0.01 non-frame, sc z>2" keeps the genes with z > 2
(in the fit, and in each permutation fit for its own null). Truth, null and estimands as `compute_atlas_gene_rhythms.py` (R, offset, rise / fall at 13 h;
host states permuted within stratum, 50 permutations); sets: followers (GTEx p < 0.05), p<0.01 non-frame, p<0.01 non-frame with sc z > 2, q<0.05
non-frame, clock-plane genes. The legacy GTEx split-half ceiling (`compare.py`) is not computed.

    python compute_atlas_gene_rhythms_modes.py [--subsample N] [--atlases heart ...] [--modes none basic pcs] [--n-null 50]
Reads gtex_cosinor{,_v2,_v3}.csv.gz (compute_gtex_truth.py). Outputs: paper/results/figure_2/output/gene_rhythms/modes_peratlas.csv,
modes_composition.csv, pred_<atlas>_<mode>.csv.gz (columns c_re, c_im: expression term; t_re, t_im: with composition; z)
Legacy: gonogo/scripts/{sc_predict,compare}.py (COV=none|basic|pcs), gonogo_new/scripts/{predict,compare_new}.py
-> output/gonogo_<atlas>[_basic|_pcs][_v2|_v3].csv, composition_<atlas>[_basic|_pcs].csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _shared as S  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "gene_rhythms"
warnings.filterwarnings("ignore")
W = S.W
N_PC = 3
COMP = "__COMP__"  # gene symbols are upper-cased by gene_rhythms
SETS = ("followers", "p<0.01 non-frame", "p<0.01 non-frame, sc z>2", "q<0.05 non-frame", "frame genes")
VERSIONS = {"v1": "", "v2": "_v2", "v3": "_v3"}


def covariate_table(pb: pd.DataFrame, donors: pd.DataFrame, genes: list, mode: str) -> tuple[pd.DataFrame, list]:
    """local: not in clocksync. Genome table with the per-row covariates of the mode and the composition pseudo-gene."""
    t = pb.copy()
    age = pd.to_numeric(donors.age, errors="coerce").reindex(t.donor)
    t["cov_age"] = age.fillna(age.median()).to_numpy()
    t["cov_male"] = (donors.sex.astype(str).reindex(t.donor) == "male").astype(float).to_numpy()
    t["log_depth"], t["log_cells"] = np.log(t.depth.to_numpy(float)), np.log(t.n_cells.to_numpy(float))
    nc = t.pivot_table(index="donor", columns="compartment", values="n_cells", aggfunc="sum").fillna(0)
    frac = nc.div(nc.sum(1), axis=0)
    fr = np.array([frac.at[d, k] for d, k in zip(t.donor, t.compartment)])
    t[COMP] = t.depth.to_numpy(float) * (fr + 1e-4) - 0.5  # log((count + 0.5) / depth) = log(frac + 1e-4)
    covs = []
    if mode in ("basic", "pcs"):
        covs = [c for c in ("cov_age", "cov_male", "log_depth", "log_cells") if np.isfinite(t[c]).all() and t[c].std() > 0]
        # the legacy standardised columns within each compartment; constant-within-compartment columns drop out of the fit
    if mode == "pcs":
        for c in range(N_PC):
            t[f"pc{c + 1}"] = 0.0
        for k, q in t.groupby("compartment"):
            C = q[genes].to_numpy(float)
            ok = C.mean(0) >= S.MIN_MEAN
            Y = np.log((C[:, ok] + 0.5) / q.depth.to_numpy(float)[:, None])
            Dm = pd.get_dummies(q.stratum.astype(str)).to_numpy(float)
            Yc = Y - Dm @ np.linalg.lstsq(Dm, Y, rcond=None)[0]
            sd = Yc.std(0)
            Yc = Yc[:, sd > 0] / sd[sd > 0]
            u = np.linalg.svd(Yc, full_matrices=False)[0][:, :N_PC]
            t.loc[q.index, [f"pc{c + 1}" for c in range(u.shape[1])]] = u
        covs += [f"pc{c + 1}" for c in range(N_PC)]
    return t, covs


def predict_all(rhy: pd.DataFrame, sh: pd.DataFrame, comp: pd.Series) -> pd.DataFrame:
    """local: not in clocksync. pred_expr, pred_total and the single-cell z per gene (genes with covered share > 0.5)."""
    c = (rhy.pivot(index="gene", columns="compartment", values="c_re") + 1j * rhy.pivot(index="gene", columns="compartment", values="c_im"))
    se = rhy.pivot(index="gene", columns="compartment", values="se")
    wv = sh.reindex(c.index)[c.columns].where(c.notna()).fillna(0.0)
    cov = wv.sum(1)
    pe = (wv * c.fillna(0)).sum(1)
    pc = (wv * comp.reindex(c.columns)).sum(1)
    z = pe.abs() / np.sqrt(((wv * se.fillna(0)) ** 2).sum(1))
    keep = cov > 0.5
    return pd.DataFrame(dict(pred_expr=pe / cov, pred_total=(pe + pc) / cov, z=z))[keep]


def gene_sets(g: pd.DataFrame, which: str, z: pd.Series | None):
    if which == "followers":
        return g.index[g.index.isin(S.FOLLOWERS) & (g.p < 0.05)]
    if which == "p<0.01 non-frame":
        return g.index[(g.p < 0.01) & ~g.index.isin(S.FRAME)]
    if which == "p<0.01 non-frame, sc z>2":
        ok = g.index[(g.p < 0.01) & ~g.index.isin(S.FRAME)]
        return ok.intersection(z.index[z > 2]) if z is not None else ok[:0]
    if which == "q<0.05 non-frame":
        return g.index[(g.q < 0.05) & ~g.index.isin(S.FRAME)]
    return g.index[g.index.isin(S.FRAME)]


def score_kind(P: pd.DataFrame, nulls: list, g: pd.DataFrame, which: str, kind: str) -> dict:
    gs = [x for x in gene_sets(g, which, P.z) if x in P.index]
    if len(gs) < 5:
        return dict(n=len(gs))
    row = S.score(P.loc[gs, kind].to_numpy(), g.loc[gs, "c"].to_numpy())
    ns = []
    for q in nulls:
        gq = [x for x in gene_sets(g, which, q.z) if x in q.index]
        if len(gq) >= 2:
            ns.append(S.st(np.angle(q.loc[gq, kind].to_numpy()), np.angle(g.loc[gq, "c"].to_numpy())))
    if ns:
        ns = np.array(ns)
        nc = np.array([S.ccorr(np.angle(q.loc[[x for x in gene_sets(g, which, q.z) if x in q.index], kind].to_numpy()),
                               np.angle(g.loc[[x for x in gene_sets(g, which, q.z) if x in q.index], "c"].to_numpy())) for q in nulls
                       if len([x for x in gene_sets(g, which, q.z) if x in q.index]) >= 3])
        row.update(null_R_95=float(np.quantile(ns[:, 0], .95)), p_R=(1 + (ns[:, 0] >= row["R_err"]).sum()) / (len(ns) + 1),
                   null_sign_95=float(np.quantile(ns[:, 2], .95)), p_sign=(1 + (ns[:, 2] >= row["sign"]).sum()) / (len(ns) + 1),
                   null_ccorr_95=float(np.quantile(nc, .95)) if len(nc) else np.nan,
                   p_perm=(1 + (nc >= row["ccorr"]).sum()) / (len(nc) + 1) if len(nc) else np.nan, n_null=len(ns))
    return row


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--atlases", nargs="+", default=S.A27 + S.A31B)
    ap.add_argument("--modes", nargs="+", default=["none", "basic", "pcs"], choices=["none", "basic", "pcs"])
    ap.add_argument("--n-null", type=int, default=None)
    ap.add_argument("--tag", default="", help="suffix of the two table names (to run atlas groups in parallel)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_null = args.n_null or (10 if run["quick"] else 50)
    out = out_dir(FOLDER, run["quick"])
    truths = {}
    for v, suf in VERSIONS.items():
        t = pd.read_csv(out / f"gtex_cosinor{suf}.csv.gz")
        truths[v] = {k: S.truth_c(d) for k, d in t.groupby("tissue")}
    genes = sorted(set().union(*[set(g.index) for g in truths["v1"].values()]))
    rows, comps = [], []
    for name in args.atlases:
        spec = S.ATLASES[name]
        pb, donors = S.genome_table(spec["genome"], args.subsample, args.seed)
        pb = S.class_rule(pb, 8 if run["quick"] else 20)
        cl = sorted(pb.compartment.unique())
        anchors = S.atlas_anchors(name, cl)
        if len(anchors) < 2 or len(cl) < 2:
            print(name, "infeasible", flush=True)
            continue
        g_in = [g for g in genes if g in pb.columns]
        ds = S.make_dataset(pb, anchors, name)
        sh = S.shares(pb, ds.compartments, g_in)
        for mode in args.modes:
            t, covs = covariate_table(pb, donors, g_in, mode)
            rhy = ds.gene_rhythms(t[list(S.META) + covs + g_in + [COMP]], genes=g_in + [COMP], min_mean=S.MIN_MEAN, n_null=n_null,
                                  seed=args.seed, covariates=covs or None)
            P, Ns, cc = None, [], None
            for b in range(-1, n_null):
                r = rhy[rhy["null"] == b]
                cr = r[r.gene == COMP].set_index("compartment")
                p = predict_all(r[r.gene != COMP], sh, cr.c_re + 1j * cr.c_im)
                if b < 0:
                    P, cc = p, cr
                else:
                    Ns.append(p)
            comps.append(pd.DataFrame(dict(atlas=name, mode=mode, compartment=cc.index, amp=np.hypot(cc.c_re, cc.c_im).to_numpy(),
                                           acro_h=(np.arctan2(cc.c_im, cc.c_re) / W % 24).to_numpy())))
            P.assign(c_re=P.pred_expr.to_numpy().real, c_im=P.pred_expr.to_numpy().imag, t_re=P.pred_total.to_numpy().real,
                     t_im=P.pred_total.to_numpy().imag)[["c_re", "c_im", "t_re", "t_im", "z"]].rename_axis("gene").to_csv(
                out / f"pred_{name}_{mode}.csv.gz")
            for v, tr in truths.items():
                for tis in spec["tissues"]:
                    if tis not in tr:
                        continue
                    for which in SETS:
                        for kind in ("pred_expr", "pred_total"):
                            rows.append(dict(atlas=name, tissue=tis, truth=v, mode=mode, set=which, pred=kind,
                                             group="A31b" if name in S.A31B else "A27",
                                             **score_kind(P, Ns, tr[tis], which, kind)))
            print(name, mode, "done", flush=True)
            pd.DataFrame(rows).to_csv(out / f"modes_peratlas{args.tag}.csv", index=False)
            pd.concat(comps).to_csv(out / f"modes_composition{args.tag}.csv", index=False)
    res = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    print(res[(res.truth == "v1") & (res.set == "p<0.01 non-frame")][["atlas", "tissue", "mode", "pred", "n", "R_err", "R_rot_h", "null_R_95", "p_R"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
