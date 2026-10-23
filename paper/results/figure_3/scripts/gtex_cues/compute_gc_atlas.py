"""A27 part 2, D7: the glucocorticoid score in the single-cell atlases (post hoc, exploratory).

Post hoc, exploratory (docs/PLAN_internal_desync.md (archive), A27 part 2, D7; no hypothesis ids). Atlases heart (nuclei), breast, skeletal
muscle, colon (the A27 part 1 atlases; genome-wide pseudobulk tables of `paper/common/datasets/*_genome`). Donor GC score = mean over
FKBP5, TSC22D3, ZBTB16, KLF9, KLF15, DDIT4 of log CPM z-scored within compartment x stratum, averaged over the donor's compartments
(composition-robust). Estimands:
  (1) R^2 of the GC score on the host state f (clock genes of the A17 anchors, no times; both centred by stratum);
  (2) the GC effect on the 14 clock genes (within compartment x stratum, pooled; the compartment's own GC and the donor's GC):
      fraction inside the GTEx waveform plane span{1, u, v} and the clock-time phase and magnitude it pulls towards;
  (3) the D3 regression (atlas pseudobulk = sum of the donor's compartments; every gene with mean count >= 10 on the host state,
      stratum dummies, age, sex, log depth) rerun "as before", with the GC score as a covariate, and with f residualised on GC
      + the GC covariate: agreement with the GTEx timed cosinor of the matched tissue (fast deaths; genes p < 0.01 outside the
      clock plane): R = |mean exp(i (arg c_est - arg c_true))|, offset (h), and the same for the clock-plane genes.
Does the GC covariate improve the non-clock agreement or remove the 12 h offset? (Legacy reading: no.)
Controlled access: GTEx time of death in memory only; aggregates are written.

    python compute_gc_atlas.py [--atlases heart breast muscle colon] [--subsample N]
Outputs (paper/results/figure_3/output/gtex_cues/): diag7_<atlas>.csv (variant scores per GTEx tissue), diag7_<atlas>_gc.csv
(GC summary), diag7_summary.csv
Legacy: legacy/desync/validation/glucocorticoid/scripts/diag7_gc_atlas.py (+ gonogo/scripts/diag3_atlas_pb.py). Local: the per-gene
regressions. Depends on gene_rhythms/_shared.py (genome tables, class rule, anchors) read-only.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _genes as G  # noqa: E402
import _shared as S  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("gene_rhythms")))
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location("s08_shared", group_dir("gene_rhythms") / "_shared.py")
S8 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S8)

W = S.W
FOLDER = G.FOLDER
CK = list(cs.gtex().genes)
FRAME = set(CK) | {"BMAL1"}
U = np.asarray(cs.gtex().X)
PB = np.c_[np.ones(len(CK)), U]
Q, _ = np.linalg.qr(PB)
A27_ATLASES = ["heart", "breast", "muscle", "colon"]


def gtex_tissue(tis):
    """Legacy diag3 truth inputs: log(CPM + 0.5) of genes with median CPM >= 5, fast deaths with a time of death, covariates sex,
    age bracket, Hardy class dummies and ischaemic time (z)."""
    Y, genes, ob = S.expression(tis, hardy=(1, 2), tod=True, median_cpm=5.0)
    Z = pd.get_dummies(pd.DataFrame(dict(sex=ob.sex, age=ob.age_bracket, hardy=ob.hardy.astype(str))),
                       drop_first=True).to_numpy(float)
    isch = ob.ischemic_min.fillna(ob.ischemic_min.median()).to_numpy(float)
    return Y, genes, ob.tod.to_numpy(float), np.c_[Z, (isch - isch.mean()) / isch.std()]


def scores(c_est, c_true):
    a, b = np.angle(c_est), np.angle(c_true)
    R = np.exp(1j * (a - b)).mean()
    return abs(R), float(np.angle(R) / W % 24)


def centre_by(x, strata):
    return x - pd.DataFrame(x).groupby(np.asarray(strata)).transform("mean").to_numpy().reshape(x.shape)


def one_atlas(name, subsample, seed):
    spec = S8.ATLASES[name]
    pb, donors = S8.genome_table(spec["genome"], subsample, seed)
    pb = S8.class_rule(pb, 8 if subsample else 20)
    comps = sorted(pb.compartment.unique())
    anchors = S8.atlas_anchors(name, comps)
    ds = S8.make_dataset(pb, anchors, name=name)
    host = ds.fit(n_boot=0, gate=False, targets=[anchors[0]]).host
    gcs = [g for g in S.GC if g in pb.columns]
    lg = np.log((pb[gcs].to_numpy(float) + 0.5) / pb.depth.to_numpy(float)[:, None] * 1e6)
    key = pb.compartment.astype(str) + "|" + pb.stratum.astype(str)
    z = pd.DataFrame(lg).groupby(key.to_numpy()).transform(lambda x: (x - x.mean()) / (x.std() + 1e-9))
    pb = pb.assign(gc_c=z.mean(axis=1).to_numpy())
    gc = pb.groupby("donor").gc_c.mean()
    don = [d for d in gc.index if d in set(host.index) and np.isfinite(gc[d])]
    strat = pb.groupby("donor").stratum.agg(lambda s: s.mode().iloc[0]).loc[don].to_numpy()
    f = centre_by(host.loc[don, ["x", "y"]].to_numpy(), strat)
    g = gc.loc[don].to_numpy()
    g = g - pd.Series(g).groupby(strat).transform("mean").to_numpy()
    bf = np.linalg.lstsq(g[:, None], f, rcond=None)[0][0]
    out = dict(atlas=name, donors=len(don), gc_sd=float(np.std(gc.loc[don])), r2_gc_on_f=float(1 - ((f - np.outer(g, bf)) ** 2).sum() / (f ** 2).sum()))
    # (2) GC effect on the clock genes within compartment x stratum
    cl = pb.set_index(["donor", "compartment"])
    lc = np.log((cl[CK].to_numpy(float) + 0.5) / cl.depth.to_numpy(float)[:, None] * 1e6)
    kk = (cl.index.get_level_values(1).astype(str) + "|" + cl.stratum.astype(str)).to_numpy()
    Y = lc - pd.DataFrame(lc).groupby(kk).transform("mean").to_numpy()
    gown = cl.gc_c.to_numpy()
    gdon = gc.reindex(cl.index.get_level_values(0)).to_numpy()
    for lab, gg in (("own", gown), ("donor", gdon)):
        ok = np.isfinite(gg) & np.isfinite(Y).all(1)
        x = gg[ok] - pd.Series(gg[ok]).groupby(kk[ok]).transform("mean").to_numpy()
        b = (x @ Y[ok]) / (x @ x)
        c = np.linalg.lstsq(PB, b, rcond=None)[0]
        out[f"{lab}_inplane"] = float(((Q @ (Q.T @ b)) ** 2).sum() / (b ** 2).sum())
        out[f"{lab}_pull_h"], out[f"{lab}_pull_mag"] = float(np.arctan2(c[2], c[1]) / W % 24), float(np.hypot(c[1], c[2]))
        for gname in ("PER1", "ARNTL", "NR1D1", "NR1D2", "DBP", "PER3", "TEF"):
            out[f"{lab}_{gname}"] = float(b[CK.index(gname)])
    # (3) D3 regression
    genes = [c for c in S8.gene_columns(pb)]
    agg = pb.groupby("donor")[genes + ["depth", "n_cells"]].sum()
    C = agg.loc[don, genes].to_numpy(float)
    okg = C.mean(0) >= 10
    Gn = np.asarray(genes)[okg]
    Yg = np.log((C[:, okg] + 0.5) / agg.loc[don, "depth"].to_numpy(float)[:, None])
    cov = donors.reindex(don)[["age", "male"]].astype(float)
    zc = np.c_[cov.fillna(cov.mean()).to_numpy(), np.log(agg.loc[don, "depth"].to_numpy(float))]
    zc = zc[:, zc.std(0) > 0]
    Zs = np.c_[pd.get_dummies(pd.Series(strat)).to_numpy(float)[:, 1:], (zc - zc.mean(0)) / zc.std(0)]
    f_res = f - np.outer(g, bf)

    def est(fm, extra):
        D = np.c_[np.ones(len(fm)), Zs, extra, fm]
        B = np.linalg.lstsq(D, Yg, rcond=None)[0]
        return pd.Series(B[-2] + 1j * B[-1], index=Gn)

    variants = {"as before": est(f, np.empty((len(f), 0))), "GC covariate": est(f, g[:, None]),
                "f residualised on GC + GC covariate": est(f_res, g[:, None])}
    rows = []
    for tis in spec["tissues"]:
        Yt, Gt, t, Zt = gtex_tissue(tis)
        c_true, p = S.timed_cosinor(Yt, t, Zt)
        truth = pd.DataFrame(dict(c=c_true, p=p), index=Gt)
        sel = truth[(truth.p < 0.01) & ~truth.index.isin(FRAME)].index
        gf = truth.index[truth.index.isin(FRAME)]
        for lab, ce in variants.items():
            gs = sel.intersection(ce.index)
            gfr = gf.intersection(ce.index)
            R, rot = scores(ce[gs].to_numpy(), truth.c[gs].to_numpy())
            fR, frot = scores(ce[gfr].to_numpy(), truth.c[gfr].to_numpy())
            rows.append(dict(atlas=name, tissue=tis, variant=lab, genes=len(gs), R=R, rot_h=rot, frame_R=fR, frame_rot_h=frot))
    return pd.DataFrame(rows), pd.Series(out)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--atlases", nargs="*", default=A27_ATLASES, choices=A27_ATLASES)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    summ = []
    for name in args.atlases:
        res, gcs = one_atlas(name, args.subsample, args.seed)
        res.to_csv(out / f"diag7_{name}.csv", index=False)
        gcs.to_csv(out / f"diag7_{name}_gc.csv")
        summ.append(gcs)
        pd.set_option("display.width", 220)
        print(gcs.to_string(), "\n", res.round(2).to_string(index=False), flush=True)
    pd.DataFrame(summ).to_csv(out / "diag7_summary.csv", index=False)


if __name__ == "__main__":
    main()
