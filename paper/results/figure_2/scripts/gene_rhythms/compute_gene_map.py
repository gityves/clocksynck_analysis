"""A50: pooled gene-rhythm maps per cell type, lung then prefrontal cortex (H50.1-H50.4).

Pre-registered (A50, written before running). Normal adult lung (HLCA 200 donors + 172 donors of 21 further Census datasets,
stratum = assay | collection) and DLPFC (452 donors; PsychAD 382, other collections 70). Pipeline of A31b / A39 unchanged: host
from the A17 anchors (lung: the lung_hlca target; DLPFC: psychad_mssm), per class and gene the loading on the host state
(`Dataset.gene_rhythms`), bulk recombination, 50 host permutations, GTEx truth (Lung; BA9 primary, Cortex secondary), score
set b (GTEx p < 0.01 non-clock-plane genes), R with p = (1 + #{null >= R}) / 51.
- H50.1: lung without HLCA predicts GTEx Lung: R above the 95 % null, p < 0.05 (sensitivity without disease-study controls).
- H50.2: pooled lung R >= HLCA alone (A31b, 0.56) and above its null; donor bootstrap interval (`--n-boot`, resampling donors
  within stratum, host and gene fits redone; 200 for the full run).
- H50.3: per class, gene phases of HLCA alone vs the other donors, over genes rhythmic in both (Wald p < 0.05 each, non-clock
  plane): circular correlation > 0 with p < 0.05 against a gene-pairing null (2000, seed 5003); classes need >= 20 such genes;
  the test passes if >= 50 % of eligible classes replicate. DLPFC: PsychAD alone vs the other collections.
- H50.4: DLPFC pooled predicts GTEx BA9: R above the 95 % null, p < 0.05; split replication as H50.3.
- Secondary: the map (per class and gene: amplitude, phase, se, p, BH q) and per-class rhythmic gene counts.

    python compute_gene_map.py [--subsample N] [--n-boot B] [--jobs J]    (N: donors per analysis set, checks)
Outputs: paper/results/figure_2/output/gene_rhythms/h50_{scores,boot,h503_<atlas>,h503_<atlas>_summary}.csv, map_<fit>.csv.gz,
map_counts.csv, h50_sets.csv
Legacy: gene_map/scripts/{pool,fit,score_gm,boot,h503,gtex_truth_gm}.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import _shared as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.anchors import anchors as a17  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "gene_rhythms"
SEED_H503, NPERM_H503 = 5003, 2000
SETS = {
    "lung": {"others": lambda m: m.source == "other", "pooled": lambda m: m.source.isin(["hlca", "other"]),
             "others_nondisease": lambda m: (m.source == "other") & ~m.disease_study.astype(bool),
             "pooled_nondisease": lambda m: (m.source == "hlca") | ((m.source == "other") & ~m.disease_study.astype(bool))},
    "dlpfc": {"pooled": lambda m: m.source == "dlpfc", "psychad": lambda m: m.group == "PsychAD",
              "other": lambda m: m.group == "other"},
}
ANCHOR_TARGET = {"lung": "lung_hlca", "dlpfc": "psychad_mssm"}
TRUTH_TISSUES = {"lung": ["Lung"], "dlpfc": S.BRAIN_TISSUES}


def anchors_for(atlas, comps):
    got = [c for c in a17(ANCHOR_TARGET[atlas]) if c in set(comps)]
    return got if len(got) >= 2 else [c for c in a17(ANCHOR_TARGET[atlas], "fixed") if c in set(comps)]


def fit_set(pb, atlas, genes, n_null, min_donors=20, name=""):
    pb = S.class_rule(pb, min_donors)
    comps = sorted(pb.compartment.unique())
    anchors = anchors_for(atlas, comps)
    info = dict(donors=pb.donor.nunique(), classes=len(comps), strata=pb.stratum.nunique(), anchors=";".join(anchors))
    if len(anchors) < 2 or len(comps) < 2:
        return None, dict(info, status=f"infeasible ({len(anchors)} anchors, {len(comps)} classes)")
    pred, nulls, main, sh, _ = S.fit_atlas(pb, anchors, genes, n_null, name=name)
    return (pred, nulls, main), dict(info, status="ok")


def score_row(fit, truths, tissue):
    pred, nulls, _ = fit
    g = truths[tissue]
    s = S.score_vs_null(pred, nulls, g, S.gene_set(g, "p<0.01 non-frame"))
    return dict(n=s.get("n"), R=s.get("R_err"), offset_h=s.get("R_rot_h"), ccorr=s.get("ccorr"), null95=s.get("null_R_95"),
                null_median=s.get("null_R_median"), n_perm=s.get("n_null"), p=s.get("p_R"))


def bh_by_class(x):
    return x.groupby("compartment").p.transform(lambda s: cs_bh(s.to_numpy()))


def cs_bh(p):
    import clocksync as cs
    return cs.stats.bh(p)


def class_table(main):
    x = main.copy()
    x["c"] = x.c_re + 1j * x.c_im
    x["phase_h"] = (np.angle(x.c) / S.W) % 24
    x["q"] = bh_by_class(x)
    x["clock_plane"] = x.gene.isin(S.FRAME)
    return x


def replicate(A, B, na, nb):
    """H50.3: per class, circular correlation of the gene phases of two halves over genes rhythmic in both."""
    rng = np.random.default_rng(SEED_H503)
    rows = []
    for k in sorted(set(A.compartment) & set(B.compartment)):
        a, b = A[A.compartment == k].drop_duplicates("gene").set_index("gene"), B[B.compartment == k].drop_duplicates("gene").set_index("gene")
        both = a.index.intersection(b.index)
        both = both[~both.isin(S.FRAME)]
        sel = both[(a.p[both] < 0.05) & (b.p[both] < 0.05)]
        row = dict(compartment=k, donors_a=int(a.n.iloc[0]), donors_b=int(b.n.iloc[0]), genes_both=len(both),
                   rhythmic_a=int((a.p[both] < 0.05).sum()), rhythmic_b=int((b.p[both] < 0.05).sum()), n_set=len(sel),
                   expected_by_chance=0.05 * 0.05 * len(both))
        if len(sel) >= 3:
            pa, pb_ = np.angle(a.c[sel].to_numpy()), np.angle(b.c[sel].to_numpy())
            obs = S.ccorr(pa, pb_)
            nul = np.array([S.ccorr(pa, pb_[rng.permutation(len(sel))]) for _ in range(NPERM_H503)])
            row.update(ccorr=obs, null95=float(np.quantile(nul, .95)), p=(1 + np.sum(nul >= obs)) / (NPERM_H503 + 1),
                       offset_h=float(np.angle(np.exp(1j * (pa - pb_)).mean()) / S.W % 24))
        row["eligible"] = len(sel) >= 20
        row["replicates"] = bool(row["eligible"] and row.get("ccorr", -1) > 0 and row.get("p", 1) < 0.05)
        rows.append(row)
    r = pd.DataFrame(rows)
    el = r[r.eligible] if len(r) else r
    s = dict(halves=f"{na} vs {nb}", classes_shared=len(r), classes_eligible=len(el),
             classes_replicating=int(el.replicates.sum()) if len(el) else 0,
             fraction=float(el.replicates.mean()) if len(el) else np.nan,
             passes=bool(len(el) > 0 and el.replicates.mean() >= 0.5))
    return r, pd.Series(s)


def boot_one(pb, genes, b):
    """H50.2 replicate: donors of the pooled lung table resampled within stratum (seed [50, b]), fits redone."""
    limit_threads(1)
    rng = np.random.default_rng([50, b])
    strat = pb.groupby("donor").stratum.first()
    pick = []
    for _, ds in strat.groupby(strat):
        ds = ds.index.to_numpy()
        pick += list(rng.choice(ds, len(ds), replace=True))
    rows_of = pb.groupby("donor").indices
    idx, ids = [], []
    for j, d in enumerate(pick):
        idx += list(rows_of[d])
        ids += [f"{d}__{j}"] * len(rows_of[d])
    pbb = pb.iloc[idx].copy()
    pbb["donor"] = ids
    pbb = pbb.reset_index(drop=True)
    comps = sorted(pbb.compartment.unique())
    pred, *_ = S.fit_atlas(pbb, anchors_for("lung", comps), genes, 0, name="boot")
    return b, pred, len(set(pick))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None, help="H50.2 bootstrap replicates (default 200; 2 for a check)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick, out = run["quick"], out_dir(FOLDER, run["quick"])
    n_null, md = (10 if quick else 50), (8 if quick else 20)
    n_boot = args.n_boot if args.n_boot is not None else (2 if quick else 200)
    tr = pd.read_csv(out / "gtex_cosinor.csv.gz")
    truths = {t: S.truth_c(d) for t, d in tr.groupby("tissue")}
    genes = sorted(set(tr.gene))
    fits, infos, scores = {}, [], []
    for atlas in ("lung", "dlpfc"):
        pb_all, meta = S.genome_table(f"gene_map_{atlas}")
        pb_all = pb_all[list(S.META) + [g for g in genes if g in pb_all.columns]]
        for name, rule in SETS[atlas].items():
            donors = meta.index[rule(meta).to_numpy()]
            if quick:
                donors = pd.Index(np.random.default_rng(args.seed).permutation(donors)[:args.subsample or 100])
            pb = pb_all[pb_all.donor.isin(donors)].reset_index(drop=True)
            tag = f"{atlas}_{name}"
            fit, info = fit_set(pb, atlas, genes, n_null, md, tag)
            infos.append(dict(fit=tag, **info))
            print(tag, info, flush=True)
            if fit is None:
                continue
            fits[tag] = fit
            for t in TRUTH_TISSUES[atlas]:
                scores.append(dict(fit=tag, tissue=t, status="ok", **score_row(fit, truths, t)))
        if atlas == "lung":
            lung_pooled_pb = pb_all[pb_all.donor.isin(meta.index[SETS["lung"]["pooled"](meta).to_numpy()])]
            if quick:
                lung_pooled_pb = lung_pooled_pb[lung_pooled_pb.donor.isin(
                    np.random.default_rng(args.seed).permutation(lung_pooled_pb.donor.unique())[:args.subsample or 100])]
            lung_pooled_pb = lung_pooled_pb.reset_index(drop=True)
    # HLCA alone (the A31b fit), the reference of H50.2 and one half of H50.3
    hl, _ = S.genome_table("lung_hlca_genome", args.subsample if quick else None, args.seed)
    hl = hl[list(S.META) + [g for g in genes if g in hl.columns]]
    fit, info = fit_set(hl, "lung", genes, n_null, md, "lung_hlca_alone")
    infos.append(dict(fit="lung_hlca_alone", **info))
    if fit is not None:
        fits["lung_hlca_alone"] = fit
        scores.append(dict(fit="lung_hlca_alone", tissue="Lung", status="ok", **score_row(fit, truths, "Lung")))
    pd.DataFrame(infos).to_csv(out / "h50_sets.csv", index=False)
    sc = pd.DataFrame(scores)
    sc.to_csv(out / "h50_scores.csv", index=False)
    pd.set_option("display.width", 250)
    print(sc.round(4).to_string(index=False))

    # maps and H50.3
    cnt, tabs = [], {}
    for tag, (pred, nulls, main_) in fits.items():
        x = class_table(main_)
        tabs[tag] = x
        x[["compartment", "gene", "n", "amplitude", "phase_h", "se", "p", "q", "clock_plane"]].to_csv(out / f"map_{tag}.csv.gz", index=False)
        g = x[~x.clock_plane].groupby("compartment").agg(donors=("n", "first"), genes=("gene", "size"),
                                                         p05=("p", lambda s: int((s < 0.05).sum())), q10=("q", lambda s: int((s < 0.1).sum())))
        cnt.append(g.reset_index().assign(fit=tag))
    pd.concat(cnt).to_csv(out / "map_counts.csv", index=False)
    for atlas, (ha, hb) in {"lung": ("lung_hlca_alone", "lung_others"), "dlpfc": ("dlpfc_psychad", "dlpfc_other")}.items():
        if ha in tabs and hb in tabs:
            r, s = replicate(tabs[ha], tabs[hb], ha, hb)
            r.to_csv(out / f"h50_h503_{atlas}.csv", index=False)
            s.to_csv(out / f"h50_h503_{atlas}_summary.csv", header=False)
            print(atlas, s.to_dict())

    # H50.2 bootstrap
    if n_boot and "lung_pooled" in fits:
        g = truths["Lung"]
        gs = S.gene_set(g, "p<0.01 non-frame")
        res = Parallel(n_jobs=args.jobs)(delayed(boot_one)(lung_pooled_pb, genes, b) for b in range(n_boot))
        rows = []
        for b, pred, nd in res:
            m = [x for x in gs if x in pred.index]
            R, off, _ = S.st(np.angle(pred.loc[m].to_numpy()), np.angle(g.loc[m, "c"].to_numpy()))
            rows.append(dict(rep=b, n_genes=len(m), R=R, offset_h=off, donors_unique=nd))
        bt = pd.DataFrame(rows)
        bt.to_csv(out / "h50_boot.csv", index=False)
        ref = float(sc.set_index("fit").loc["lung_hlca_alone", "R"]) if "lung_hlca_alone" in set(sc.fit) else np.nan
        s = dict(reps=len(bt), R_pooled=float(sc.set_index("fit").loc["lung_pooled", "R"]), R_hlca_alone=ref, boot_median=bt.R.median(),
                 lo95=bt.R.quantile(.025), hi95=bt.R.quantile(.975), frac_ge_hlca=float((bt.R >= ref).mean()))
        pd.Series(s).to_csv(out / "h50_boot_summary.csv", header=False)
        print(s)


if __name__ == "__main__":
    main()
