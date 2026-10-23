"""A27 part 1 and A31b: can untimed single-cell atlases give the rhythms of non-clock genes? Go / no-go against GTEx.

A27 part 1 (post hoc, exploratory; no acceptance criterion): heart (nuclei), skeletal muscle, breast, colon atlases.
A31b (pre-registered, confirmatory): HLCA lung, HBCA breast, Xu liver, Ganier skin, Emont adipose; selection fixed before
scoring. Per atlas: the host from the clock genes of the A17 anchors (leave-one-out for an anchor), then per class and gene
the loading on the host state (`Dataset.gene_rhythms`, untimed, one intercept per stratum), recombined into a bulk
rhythm with RNA shares (local). Truth: GTEx timed cosinor of the matched tissue (`compute_gtex_truth.py`, v1 time of death,
v2 with RIN and covariates, v3 CHIRAL DIP). Null: host states permuted within stratum.
Estimand R = |mean exp(i(pred - true))| (offset-free concentration), the offset, rise / fall agreement at 13 h of the host
frame, per gene set: (a) the 32 frozen D10 follower genes with GTEx p < 0.05, (b) GTEx p < 0.01 non-clock-plane genes,
(c) the clock-plane genes. A31b primary (pooled over the five atlases, one tissue each, v1): pooled follower R above the
95 % quantile of the pooled per-permutation null AND rise / fall agreement > 0.5 with p < 0.05.

    python compute_atlas_gene_rhythms.py [--subsample N] [--atlases heart lung_hlca ...] [--n-null 50]
Outputs: paper/results/figure_2/output/gene_rhythms/atlas_{peratlas,pooled,primary,classes}.csv, pred_<atlas>.csv.gz,
null_<atlas>.csv.gz (all genes / genes with GTEx p < 0.05, loadings complex as c_re, c_im).
Legacy: gonogo/scripts/{sc_predict,compare}.py (A27), gonogo_new/scripts/{atlases,predict,score}.py (A31b).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _shared as S  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "gene_rhythms"
VERSIONS = {"v1": "", "v2": "_v2", "v3": "_v3"}
SETS = ("followers", "p<0.01 non-frame", "frame genes")


def load_truths(out):
    """{version: {tissue: truth with complex c}} from compute_gtex_truth.py."""
    res = {}
    for v, suf in VERSIONS.items():
        f = out / f"gtex_cosinor{suf}.csv.gz"
        if f.exists():
            t = pd.read_csv(f)
            res[v] = {k: S.truth_c(d) for k, d in t.groupby("tissue")}
    return res


def run_atlas(name, args, run, genes, quick):
    spec = S.ATLASES[name]
    pb, _ = S.genome_table(spec["genome"], args.subsample, args.seed)
    pb = S.class_rule(pb, 8 if quick else 20)
    comps = sorted(pb.compartment.unique())
    anchors = S.atlas_anchors(name, comps)
    info = dict(atlas=name, donors=pb.donor.nunique(), classes=len(comps), anchors=";".join(anchors),
                strata=pb.stratum.nunique())
    print(info, flush=True)
    if len(anchors) < 2 or len(comps) < 2:
        return None, dict(info, status=f"infeasible ({len(anchors)} anchors, {len(comps)} classes)")
    pred, nulls, main, sh, ds = S.fit_atlas(pb, anchors, genes, args.n_null or run["n_null"], seed=args.seed, name=name)
    info.update(status="ok", n_genes_pred=len(pred), n_null=len(nulls))
    return (pred, nulls, main, sh), info


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--atlases", nargs="+", default=S.A27 + S.A31B)
    ap.add_argument("--n-null", type=int, default=None, help="host permutations (default 50 for a full run, 10 for a check)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    run["n_null"] = 10 if run["quick"] else 50
    out = out_dir(FOLDER, run["quick"])
    truths = load_truths(out)
    genes = sorted(set().union(*[set(g.index) for g in truths["v1"].values()]))
    fits, infos = {}, []
    for name in args.atlases:
        res, info = run_atlas(name, args, run, genes, run["quick"])
        infos.append(info)
        if res is None:
            continue
        pred, nulls, main_, sh = res
        fits[name] = res
        sig = set().union(*[set(t.index[t.p < 0.05]) for tr_ in truths.values() for k_, t in tr_.items()
                            if k_ in S.ATLASES[name]["tissues"]])
        keep = [g for g in pred.index if g in sig or g in S.FRAME or g in S.FOLLOWERS]
        pd.DataFrame(dict(c_re=pred.values.real, c_im=pred.values.imag), index=pred.index.rename("gene")) \
            .to_csv(out / f"pred_{name}.csv.gz")  # all genes (A30 reads it); nulls only for the scored genes
        pd.concat([pd.DataFrame(dict(perm=b, gene=q.index, c_re=q.values.real, c_im=q.values.imag)).query("gene in @keep")
                   for b, q in enumerate(nulls)]).to_csv(out / f"null_{name}.csv.gz", index=False)
        cl = main_.assign(rhythmic_p05=main_.p < 0.05).query("gene not in @S.FRAME").groupby("compartment") \
            .agg(donors=("n", "first"), genes=("gene", "size"), rhythmic_p05=("rhythmic_p05", "sum")).reset_index().assign(atlas=name)
        infos[-1]["classes_table"] = cl
    pd.concat([i.pop("classes_table") for i in infos if "classes_table" in i]).to_csv(out / "atlas_classes.csv", index=False)
    pd.DataFrame(infos).to_csv(out / "atlas_fits.csv", index=False)

    per, pooled = [], []
    for v, tr in truths.items():
        for which in SETS:
            d, nd = [], []
            for name, (pred, nulls, *_r) in fits.items():
                for tis in S.ATLASES[name]["tissues"]:
                    if tis not in tr:
                        continue
                    g = tr[tis]
                    gs = [x for x in S.gene_set(g, which) if x in pred.index]
                    row = dict(truth=v, mode="none", atlas=name, tissue=tis, set=which, group="A31b" if name in S.A31B else "A27")
                    row.update(S.score_vs_null(pred, nulls, g, gs))
                    if which != "p<0.01 non-frame":
                        row["genes"] = ";".join(gs)
                    per.append(row)
                    if tis == S.POOL_TISSUE[name] and name in S.A31B and len(gs):
                        d.append(pd.DataFrame(dict(atlas=name, gene=gs, a=np.angle(pred.loc[gs].to_numpy()),
                                                   b=np.angle(g.loc[gs, "c"].to_numpy()))))
                        nd.append(pd.DataFrame([dict(atlas=name, perm=b, a=np.angle(q.reindex(gs).dropna().to_numpy()),
                                                     b=np.angle(g.loc[q.reindex(gs).dropna().index, "c"].to_numpy()))
                                                for b, q in enumerate(nulls)]))
            if not d:
                continue
            D = pd.concat(d)
            R, rot, sign = S.st(D.a.to_numpy(), D.b.to_numpy())
            ND = pd.concat(nd)
            ns = np.array([S.st(np.concatenate(q.a.values), np.concatenate(q.b.values)) for _, q in ND.groupby("perm")])
            pooled.append(dict(truth=v, mode="none", set=which, atlases=";".join(D.atlas.unique()), n=len(D), R=R, rot_h=rot,
                               null_R_95=np.quantile(ns[:, 0], .95), p_R=(1 + (ns[:, 0] >= R).sum()) / (len(ns) + 1), sign=sign,
                               null_sign_med=np.median(ns[:, 2]), null_sign_95=np.quantile(ns[:, 2], .95),
                               p_sign=(1 + (ns[:, 2] >= sign).sum()) / (len(ns) + 1), n_null=len(ns),
                               per_atlas=str({k: f"R={S.st(x.a.to_numpy(), x.b.to_numpy())[0]:.2f},"
                                                 f"sign={S.st(x.a.to_numpy(), x.b.to_numpy())[2]:.2f},n={len(x)}"
                                              for k, x in D.groupby("atlas")})))
    per, pooled = pd.DataFrame(per), pd.DataFrame(pooled)
    per.to_csv(out / "atlas_peratlas.csv", index=False)
    pooled.to_csv(out / "atlas_pooled.csv", index=False)
    if len(pooled):
        p = pooled[(pooled.truth == "v1") & (pooled.set == "followers")]
        if len(p):
            p = p.iloc[0]
            prim = dict(n=p.n, atlases=p.atlases, R=p.R, null_R_95=p.null_R_95, p_R=p.p_R, R_pass=bool(p.R > p.null_R_95),
                        sign=p.sign, p_sign=p.p_sign, sign_pass=bool(p.sign > 0.5 and p.p_sign < 0.05))
            prim["primary_pass"] = prim["R_pass"] and prim["sign_pass"]
            pd.Series(prim).to_csv(out / "atlas_primary.csv", header=False)
            print(pd.Series(prim).to_string())
    pd.set_option("display.width", 250)
    if len(per):
        print(per.query("truth == 'v1' and set == 'p<0.01 non-frame'")[["atlas", "tissue", "n", "R_err", "R_rot_h", "null_R_95", "p_R"]]
              .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
