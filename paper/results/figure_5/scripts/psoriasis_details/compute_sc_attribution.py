"""A60 and A61: which cells could carry the bulk psoriatic lesion delay, and is it a keratinocyte state mix? Descriptive.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A60, A61; the A60 scaled-NNLS variant is post hoc, labelled; A60
correction of 2026-10-17: ARNTL / BMAL1 collapsed). Nothing here is a test.
A60: (1) coupling map: the A58 status-split fit (builder `psoriasis_sc`), coupling, |z| and offset per class x status,
per-class lesional-minus-uninvolved offset, unit-bootstrap intervals; (2) bulk weights: NNLS deconvolution of the four
A54 bulk count series (GSE197056, GSE268958, GSE121212, GSE283313) with whole-skin class signatures (30 markers per class,
>= 20 CPM) -> class mRNA fractions m_k, class share w_kg of every clock gene; also the single-cell proportions and
(post hoc) the scaled NNLS; (3) mixing (`_attribution.py`): predicted bulk delay arg(z_L / z_N), composition-only,
within-class only, per-class contributions; (4) gene fingerprint: Spearman over genes of the A54 per-gene bulk shift
(`a54_*_genes.csv`, `compute_replication.py`) against each class's lesional share.
A61: the same with compartments endothelial, fibroblast, mural and the keratinocyte states kc_basal / kc_cycling /
kc_differentiated (builder `psoriasis_sc_kc_states`); mixing over the states (state mRNA share x state CPM per clock gene);
reading rule: state mix if composition < -0.5 h in >= 80 % of bootstraps, dysregulation if within < 0 in >= 80 % and >= 2
states negative. Without GSE278330 as a sensitivity; unit bootstrap (500 / 200).

    python compute_sc_attribution.py [--subsample N] [--n-boot 500] [--jobs 4] [--only a60|a61]
Outputs (paper/results/figure_5/output/psoriasis_details/): a60_{weights,gene_shares,fingerprint}.csv,
a60_{all,no_GSE278330}_{coupling,prediction}.csv, a61_{all,no_GSE278330}_{coupling,state_offsets,prediction,state_shares,
verdict}.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a60_attribution.py, a61_kc_states.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import _attribution as M  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
import _singlecell as C  # noqa: E402
from common import cli  # noqa: E402

SKC = ["kc_basal", "kc_cycling", "kc_differentiated"]


def one(t, ks, states, gw, seed=None):
    tt = C.boot(t, np.random.default_rng(seed)) if seed is not None else t
    z = M.fit_z(tt, ks)
    return M.summary(z, ks), {name: M.predict(z, g, states) for name, g in gw.items()}


def intervals(s0, bs, name):
    allb = pd.concat([b[0].assign(rep=i) for i, b in enumerate(bs)])
    rows = []
    for _, r in s0.iterrows():
        bb = allb[(allb.cls == r.cls) & (allb.status == r.status)]
        lo, hi = C.interval(bb.offset_h, r.offset_h)
        rows.append(dict(r, coupling_lo=bb.coupling.quantile(0.025), coupling_hi=bb.coupling.quantile(0.975),
                         offset_lo=lo, offset_hi=hi))
    return pd.DataFrame(rows), allb


def predictions(p0, bs, names):
    pr = []
    for name in names:
        for key, v in p0[name].items():
            bv = np.array([b[1][name].get(key, np.nan) for b in bs], float)
            bv = bv[np.isfinite(bv)]
            lo, hi = C.interval(bv, v)
            pr.append(dict(weights=name, quantity=key, h=v, lo=lo, hi=hi, frac_negative=float(np.mean(bv < 0)),
                           n_boot=len(bv)))
    return pd.DataFrame(pr)


def run_a60(t, ks, gw, tag, n_boot, jobs, out, seed=60):
    s0, p0 = one(t, ks, ks, gw)
    seeds = np.random.default_rng(seed).integers(0, 2**31, n_boot)
    bs = [b for b in Parallel(n_jobs=jobs)(delayed(C.safe)(one, t, ks, ks, gw, s) for s in seeds) if b]
    cm, _ = intervals(s0, bs, tag)
    cm.to_csv(out / f"a60_{tag}_coupling.csv", index=False)
    pr = predictions(p0, bs, list(gw))
    pr.to_csv(out / f"a60_{tag}_prediction.csv", index=False)
    print(tag, "\n", cm.round(2).to_string(), "\n", pr.round(2).to_string(), flush=True)


def a60(args, run, out, n):
    pb, _ = C.load("psoriasis_sc", args.subsample, args.seed)
    t = C.status_split(pb)
    ks = C.classes(t)
    sig, prop = M.signatures(M.S_allgenes("pseudobulk_all"), ks)
    b = M.bulk()
    mk = M.markers(sig["pooled"], b.columns)
    m, ms = M.fractions(b, sig["pooled"], mk), M.fractions(b, sig["pooled"], mk, scaled=True)
    L = M.bulk_projector(b)
    w = {s: M.shares(m[s], sig[s]) for s in "LN"}
    wp = {s: M.shares(prop[s], sig[s]) for s in "LN"}
    ws = {s: M.shares(ms[s], sig[s]) for s in "LN"}
    gw = {"deconvolution": {s: M.gains(w[s], L) for s in "LN"}, "sc_proportion": {s: M.gains(wp[s], L) for s in "LN"},
          "deconvolution_scaled_posthoc": {s: M.gains(ws[s], L) for s in "LN"}}
    wt = pd.DataFrame({"m_L": m["L"], "m_N": m["N"], "m_scaled_L": ms["L"], "m_scaled_N": ms["N"],
                       "sc_prop_L": prop["L"], "sc_prop_N": prop["N"], "gain_L": gw["deconvolution"]["L"].abs(),
                       "gain_N": gw["deconvolution"]["N"].abs()})
    wt.index.name = "cls"
    wt.to_csv(out / "a60_weights.csv")
    pd.concat({f"{nm}_{s}": x[s] for nm, x in (("deconvolution", w), ("sc_proportion", wp), ("scaled_posthoc", ws))
               for s in "LN"}, axis=1).to_csv(out / "a60_gene_shares.csv")
    print(f"{len(mk)} markers\n", wt.round(3).to_string(), flush=True)
    src = H.out_dir(H.FOLDER, run["quick"])  # A54 per-gene shifts written by compute_replication.py
    if not (src / "a54_discovery_genes.csv").exists():
        src = H.out_dir(H.FOLDER)
    sh = pd.concat([pd.read_csv(src / f"a54_{x}_genes.csv", index_col=0).gene_shift_h for x in ("discovery", "heldout")],
                   axis=1).mean(1)
    fp = pd.DataFrame([dict(cls=k, spearman=spearmanr(w["L"][k], sh.reindex(w["L"].index)).statistic,
                            p_descriptive=spearmanr(w["L"][k], sh.reindex(w["L"].index)).pvalue) for k in ks])
    fp.to_csv(out / "a60_fingerprint.csv", index=False)
    print(fp.round(2).to_string(), flush=True)
    run_a60(t, ks, gw, "all", n["n_boot"], args.jobs, out)
    t2 = t[t.stratum != "GSE278330"].reset_index(drop=True)
    run_a60(t2, C.classes(t2), gw, "no_GSE278330", max(n["n_boot"] * 2 // 5, 10), args.jobs, out)


def state_weights(studies):
    a = M.S_allgenes("pseudobulk_kc_states_all")
    a = a[a.obs.study.isin(studies)]
    ob = a.obs.reset_index(drop=True)
    cpm, share = {}, {}
    for st, s in C.STATUS.items():
        m = (ob.status == st).to_numpy()
        cpm[s] = M.class_cpm(a, {k: m & (ob.compartment == k).to_numpy() for k in SKC})
        share[s] = ob[m].groupby("compartment").depth.sum().reindex(SKC).pipe(lambda v: v / v.sum())
    return cpm, pd.DataFrame(share)


def run_a61(t, tag, L, jobs, n_boot, out, seed=61):
    ks = C.classes(t)
    cpm, share = state_weights(t.stratum.unique())
    states = [k for k in SKC if k in ks]
    g = {s: M.gains(M.shares(share[s].reindex(states), cpm[s][states]), L) for s in "LN"}
    gw = {"states": g}
    s0, p0 = one(t, ks, states, gw)
    seeds = np.random.default_rng(seed).integers(0, 2**31, n_boot)
    bs = [b for b in Parallel(n_jobs=jobs)(delayed(C.safe)(one, t, ks, states, gw, s) for s in seeds) if b]
    cm, allb = intervals(s0, bs, tag)
    diff = []
    p = s0.set_index(["cls", "status"]).offset_h
    bv = allb.pivot_table(index="rep", columns=["cls", "status"], values="offset_h")
    for k in ks:
        if (k, "L") not in p or (k, "N") not in p:
            continue
        v = float(C.wrap(p[(k, "L")] - p[(k, "N")]))
        bd = np.asarray(C.wrap((bv[(k, "L")] - bv[(k, "N")]).dropna().to_numpy()))
        lo, hi = C.interval(bd, v)
        diff.append(dict(cls=k, l_minus_n_h=v, lo=lo, hi=hi, frac_negative=float(np.mean(bd < 0))))
    diff = pd.DataFrame(diff)
    pr = predictions(p0, bs, ["states"]).drop(columns="weights").set_index("quantity")
    neg = int((diff[diff.cls.isin(SKC)].l_minus_n_h < 0).sum())
    mix = pr.h["composition"] < -0.5 and pr.frac_negative["composition"] >= 0.8
    dys = pr.frac_negative["within"] >= 0.8 and neg >= 2
    reading = "both" if mix and dys else "state mix" if mix else "dysregulation" if dys else "undetermined"
    v = pd.Series(dict(units=t.donor.nunique(), compartments=";".join(ks), n_boot=len(bs), kc_states_negative=neg,
                       reading=reading))
    cm.to_csv(out / f"a61_{tag}_coupling.csv", index=False)
    diff.to_csv(out / f"a61_{tag}_state_offsets.csv", index=False)
    pr.to_csv(out / f"a61_{tag}_prediction.csv")
    share.to_csv(out / f"a61_{tag}_state_shares.csv")
    v.to_csv(out / f"a61_{tag}_verdict.csv", header=False)
    print(tag, "\n", share.round(3).to_string(), "\n", diff.round(2).to_string(), "\n", pr.round(2).to_string(), "\n",
          v.to_string(), flush=True)


def a61(args, run, out, n):
    sv, _ = C.load("psoriasis_sc", args.subsample, args.seed)
    kc, _ = C.load("psoriasis_sc_kc_states", args.subsample, args.seed)
    pb = pd.concat([sv[sv.compartment.isin(C.DERMAL)], kc], ignore_index=True)
    t = C.status_split(pb)
    L = M.bulk_projector(M.bulk())
    run_a61(t, "all", L, args.jobs, n["n_boot"], out)
    run_a61(t[t.stratum != "GSE278330"].reset_index(drop=True), "no_GSE278330", L, args.jobs,
            max(n["n_boot"] * 2 // 5, 10), out)


def main(argv=None):
    ap = H.parser(__doc__)
    ap.add_argument("--only", choices=("a60", "a61"), default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    n = C.sizes(args, run, 2000, 500)
    if args.only in (None, "a60"):
        a60(args, run, out, n)
    if args.only in (None, "a61"):
        a61(args, run, out, n)


if __name__ == "__main__":
    main()
