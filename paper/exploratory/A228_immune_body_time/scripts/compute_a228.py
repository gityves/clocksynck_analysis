"""A228: does the tumour immune clock state carry body time, and is the T/NK pairing a composition effect?

Registered in paper/docs/registered_experiments.md (A228). The 34 A227 cohorts, tables, eligibility, statistic, permutation null,
pooling and pseudo-clock sets of A227 (`compute_a227`); normal-side rows are permuted within cohort.

Part 1 (cross-compartment, reference = the patient's adjacent-normal fibroblast): H228.1 tumour T_NK, H228.2 tumour myeloid;
  descriptive: the same against adjacent-normal endothelium, tumour endothelium vs normal fibroblast (positive control),
  adjacent-normal T_NK / myeloid vs adjacent-normal fibroblast.
Part 2 (T_NK subsets, per-cell marker rule on raw counts; CD3 = CD3D+CD3E+CD3G, CD8 = CD8A+CD8B):
  CD8 T: CD3 >= 1 & CD8 >= 1 & CD4 = 0 (H228.3); CD8-negative T: CD3 >= 1 & CD8 = 0 (H228.4);
  NK: CD3 = 0 & KLRF1+NCAM1+FCGR3A+GNLY+NKG7 >= 2 (H228.5); other cells dropped. Subset pseudobulks per patient and sample kind
  (>= 50 cells, sample median UMI >= 1,500 as the A227 units); genome-wide, cached in `$SCCC_DATA/a228/<cohort>.parquet`.
  The cohort-level compartment-presence rule is not applied to subsets (it concerns the A149 compartments).

  power  patients / cohorts per test and the within-cohort permutation null (SHUFFLED pairings only: the identity pairing is never
         evaluated, no observed statistic is computed, stored or printed), pooled with patient-count weights: null mean, SD,
         MDE = 2.49 x SD, scored (>= 40 patients, MDE <= 0.25); builds the subset caches; subset-fraction descriptives.
                                  -> output/{eligibility,power,subset_status,subset_fractions,subset_fraction_summary}.csv
  tests  the registered statistics (IMPLEMENTED, NOT RUN before the power note is committed): observed pooled uncorrected
         |complex correlation|, permutation p, 200 pseudo-clock sets, Holm over scored Part 1 primaries and, separately, scored
         Part 2 primaries; descriptive per cancer type (>= 12 patients).   -> output/{tests,pooled,per_cohort,pseudo,descriptive}.csv

    python compute_a228.py --stage power|tests [--jobs 2] [--n-perm 10000] [--n-pseudo 200] [--force] [--smoke --out DIR]
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
for p in ("A227_tme_compartment_pairing", "A149_paired_tumour_normal", "A151_fibroblast_pairing", "A226_stroma_pairing_confirm"):
    sys.path.insert(0, str(PAPER / "exploratory" / p / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
import compute_a226 as K  # noqa: E402
import compute_a227 as A  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
CACHE = L.DATA / "a228"
SEED = 228
MIN_CELLS = 50
T, N = "tumour", "adjacent_normal"
# test name -> (hypothesis, role, part, tumour-side label, normal-side label)
TESTS = {
    "T_NK_vs_fibroblast": ("H228.1", "primary", 1, f"{T}:T_NK", f"{N}:fibroblast"),
    "myeloid_vs_fibroblast": ("H228.2", "primary", 1, f"{T}:myeloid", f"{N}:fibroblast"),
    "T_NK_vs_endothelial": (None, "descriptive", 1, f"{T}:T_NK", f"{N}:endothelial"),
    "myeloid_vs_endothelial": (None, "descriptive", 1, f"{T}:myeloid", f"{N}:endothelial"),
    "endothelial_vs_fibroblast": (None, "control", 1, f"{T}:endothelial", f"{N}:fibroblast"),
    "normalT_NK_vs_fibroblast": (None, "descriptive", 1, f"{N}:T_NK", f"{N}:fibroblast"),
    "normalmyeloid_vs_fibroblast": (None, "descriptive", 1, f"{N}:myeloid", f"{N}:fibroblast"),
    "CD8T": ("H228.3", "primary", 2, f"{T}:CD8T", f"{N}:CD8T"),
    "CD8neg_T": ("H228.4", "primary", 2, f"{T}:CD8neg_T", f"{N}:CD8neg_T"),
    "NK": ("H228.5", "primary", 2, f"{T}:NK", f"{N}:NK"),
}
KINDS = list(TESTS)
SUBSETS = ["CD8T", "CD8neg_T", "NK"]
HYP = {k: (v[0], v[1]) for k, v in TESTS.items()}
MIN_PATIENTS, MDE_K, MAX_MDE, MIN_CANCER = A.MIN_PATIENTS, A.MDE_K, A.MAX_MDE, A.MIN_CANCER
MARK = {"CD3": ["CD3D", "CD3E", "CD3G"], "CD8": ["CD8A", "CD8B"], "CD4": ["CD4"],
        "NK": ["KLRF1", "NCAM1", "FCGR3A", "GNLY", "NKG7"]}
META = ["donor", "compartment", "stratum", "n_cells", "depth"]


# ---- T_NK subsets ----------------------------------------------------------------------------------------------

def marker_cols(sym):
    """Column indices per marker group (aliases via `canonical`, duplicated symbols summed) and the missing genes."""
    can = np.array([L.canonical(str(g)) for g in sym])
    cols, missing = {}, []
    for grp, genes in MARK.items():
        cols[grp] = []
        for g in genes:
            i = np.flatnonzero(can == L.canonical(g))
            cols[grp].append(i)
            if not len(i):
                missing.append(g)
    return cols, missing


def classify(X, cols):
    """Subset of each cell (rows of the counts matrix X): CD8T / CD8neg_T / NK / dropped."""
    def tot(grp):
        idx = np.concatenate(cols[grp]) if len(cols[grp]) else np.array([], int)
        return np.asarray(X[:, idx].sum(1)).ravel() if len(idx) else np.zeros(X.shape[0])
    cd3, cd8, cd4, nk = tot("CD3"), tot("CD8"), tot("CD4"), tot("NK")
    lab = np.full(X.shape[0], "dropped", dtype=object)
    lab[(cd3 >= 1) & (cd8 >= 1) & (cd4 == 0)] = "CD8T"
    lab[(cd3 >= 1) & (cd8 == 0)] = "CD8neg_T"
    lab[(cd3 == 0) & (nk >= 2)] = "NK"
    return lab


def subset_file(study):
    return CACHE / f"{study}.parquet"


def build_subsets(study):
    """Genome-wide subset pseudobulks of the T_NK compartment (rows: patient x sample kind x subset; columns META + genes) and the
    per-sample subset fractions. Cached. Returns (study, status, fractions)."""
    from clocksync.config import limit_threads
    limit_threads(1)
    K.install()
    ff = CACHE / f"{study}_fractions.csv"
    if subset_file(study).exists() and ff.exists():
        return study, "ok", pd.read_csv(ff, dtype={"patient": str, "sample": str})
    us, _ = L.ST.units(study)
    us = [u for u in L.ST.majority_chemistry(us)]
    sums, ncell, depth, strat, genes, frac, missing = {}, {}, {}, {}, None, [], set()
    for sample, patient, kind, load, chem, _ in us:
        X, _g, umi, comp, barcodes, sym = load()
        lab = L.B._labels(comp, kind, sample, barcodes, None)
        i = np.flatnonzero(lab == f"{kind}:T_NK")
        if not len(i):
            continue
        cols, miss = marker_cols(sym)
        missing |= set(miss)
        if not len(np.concatenate(cols["CD3"])) or not len(np.concatenate(cols["CD8"])) or not len(cols["CD4"][0]):
            return study, "no CD3 / CD8 / CD4 gene in the count matrix", None
        genes = pd.Index([L.canonical(str(g)) for g in sym]) if genes is None else genes
        Xt = sp.csr_matrix(X[i])
        sub = classify(Xt, cols)
        frac.append(dict(study=study, patient=patient, sample=sample, kind=kind, n_T_NK=len(i),
                         **{f"n_{s}": int((sub == s).sum()) for s in SUBSETS + ["dropped"]}))
        for s in SUBSETS:
            j = np.flatnonzero(sub == s)
            if len(j):
                k = (patient, f"{kind}:{s}")
                sums[k] = sums.get(k, 0) + np.asarray(Xt[j].sum(0)).ravel()
                ncell[k] = ncell.get(k, 0) + len(j)
                depth[k] = depth.get(k, 0) + float(np.asarray(umi)[i][j].sum())
                strat[k] = chem
    if not sums:
        return study, "no T_NK cells", None
    keys = sorted(sums)
    G = pd.DataFrame(np.vstack([sums[k] for k in keys]), columns=genes)
    G = G.T.groupby(level=0).sum().T
    keepg = (G.sum(0) > 0) | G.columns.isin(list(L.q.WF.genes))
    G = G.loc[:, keepg].astype(np.float32)
    meta = pd.DataFrame({"donor": [k[0] for k in keys], "compartment": [k[1] for k in keys], "stratum": [strat[k] for k in keys],
                         "n_cells": [ncell[k] for k in keys], "depth": [depth[k] for k in keys]})
    CACHE.mkdir(parents=True, exist_ok=True)
    F_ = pd.DataFrame(frac)
    F_.to_csv(ff, index=False)
    pd.concat([meta, G.reset_index(drop=True)], axis=1).to_parquet(subset_file(study))
    return study, ("ok" if not missing else "ok; genes absent: " + ",".join(sorted(missing))), F_


def subset_table(study, columns=None):
    """Subset pseudobulks with >= 50 cells (None if the cohort has no cache)."""
    f = subset_file(study)
    if not f.exists():
        return None
    t = pd.read_parquet(f, columns=columns)
    return t[t.n_cells >= MIN_CELLS].reset_index(drop=True)


def wf_columns():
    return META + list(L.q.WF.genes)


# ---- units and nulls -------------------------------------------------------------------------------------------

def seed_of(kind, study):
    return [SEED, KINDS.index(kind), sum(map(ord, study))]


def pb_of(study, kind):
    """Table for the real states of a test: A227 cohort table (Part 1) or subset table restricted to the clock genes (Part 2)."""
    if TESTS[kind][2] == 1:
        return A.table(study)
    return subset_table(study, wf_columns())


def units_of(study, kind):
    pb = pb_of(study, kind)
    if pb is None:
        return None
    a, b = TESTS[kind][3:]
    if a not in set(pb.compartment) or b not in set(pb.compartment):
        return None
    return L.pair_states(pb, a, b)


def eligibility(study):
    from clocksync.config import limit_threads
    limit_threads(1)
    rows = []
    for kind, (_, _, part, a, b) in TESTS.items():
        pb = pb_of(study, kind)
        if pb is None:
            rows.append(dict(study=study, cancer=A.cancer(study), kind=kind, patients_both=0, eligible=False))
            continue
        A_, B_ = pb[pb.compartment == a], pb[pb.compartment == b]
        both = sorted(set(A_.donor) & set(B_.donor))
        sel = pd.concat([A_[A_.donor.isin(both)], B_[B_.donor.isin(both)]])
        rows.append(dict(study=study, cancer=A.cancer(study), kind=kind, patients_tumour_side=A_.donor.nunique(),
                         patients_normal_side=B_.donor.nunique(), patients_both=len(both), eligible=len(both) >= L.MIN_PAIRED,
                         median_cells_per_pseudobulk=float(sel.n_cells.median()) if len(sel) else np.nan))
    return rows


def power_one(study, n_perm):
    from clocksync.config import limit_threads
    limit_threads(1)
    out = {}
    for kind in KINDS:
        u = units_of(study, kind)
        if u:
            out[kind] = (u["n"], A.perm_null_raw(u, n_perm, np.random.default_rng(seed_of(kind, study))))
    return study, out


def fraction_tables(F_, args):
    F_.to_csv(args.out / "subset_fractions.csv", index=False)
    d = F_.copy()
    for s in SUBSETS + ["dropped"]:
        d[f"f_{s}"] = d[f"n_{s}"] / d.n_T_NK
    rows = []
    for (kind, ), g in d.groupby(["kind"]):
        r = dict(kind=kind, samples=len(g), patients=g.groupby("study").patient.nunique().sum(), cohorts=g.study.nunique(),
                 median_T_NK_cells=float(g.n_T_NK.median()))
        for s in SUBSETS + ["dropped"]:
            r[f"median_fraction_{s}"] = float(g[f"f_{s}"].median())
            r[f"cell_weighted_fraction_{s}"] = float(g[f"n_{s}"].sum() / g.n_T_NK.sum())
        rows.append(r)
    pd.DataFrame(rows).to_csv(args.out / "subset_fraction_summary.csv", index=False)


def stage_power(args):
    cohorts = A.cohort_list(args.smoke)
    sub = Parallel(n_jobs=args.jobs)(delayed(build_subsets)(s) for s in cohorts)
    status = pd.DataFrame([dict(study=s, cancer=A.cancer(s), subset_status=st) for s, st, _ in sub])
    status.to_csv(args.out / "subset_status.csv", index=False)
    fr = [f for _, _, f in sub if f is not None]
    if fr:
        fraction_tables(pd.concat(fr, ignore_index=True), args)
    elig = pd.DataFrame([r for x in Parallel(n_jobs=args.jobs)(delayed(eligibility)(s) for s in cohorts) for r in x])
    elig.to_csv(args.out / "eligibility.csv", index=False)
    res = dict(Parallel(n_jobs=args.jobs)(delayed(power_one)(s, args.n_perm) for s in cohorts))
    rows = []
    for kind in KINDS:
        cs = [s for s in cohorts if kind in res[s]]
        ns = [res[s][kind][0] for s in cs]
        nul = A.pool_null(ns, [res[s][kind][1] for s in cs]) if cs else np.array([np.nan])
        sd, n = float(nul.std()), int(sum(ns))
        rows.append(dict(kind=kind, hypothesis=HYP[kind][0], role=HYP[kind][1], part=TESTS[kind][2], tumour_side=TESTS[kind][3],
                         normal_side=TESTS[kind][4], patients=n, cohorts=len(cs), cohort_list=";".join(cs),
                         null_mean=float(nul.mean()), null_sd=sd, mde=MDE_K * sd, n_perm=args.n_perm,
                         scored=bool(n >= MIN_PATIENTS and MDE_K * sd <= MAX_MDE)))
    P = pd.DataFrame(rows)
    P.to_csv(args.out / "power.csv", index=False)
    print(P.drop(columns="cohort_list").round(4).to_string(), flush=True)
    print(status.subset_status.value_counts().to_string(), flush=True)


# ---- tests (registered statistics; not run before the power note) ----------------------------------------------

def pseudo_one(study, kind, n_pseudo):
    """Expression-matched pseudo-clock sets (A149 note 2 (iv), as A227): reference log-mean from the normal-side compartment."""
    from clocksync.config import limit_threads
    from clocksync.paired import pseudo_waveforms
    limit_threads(1)
    a, b = TESTS[kind][3:]
    u = units_of(study, kind)
    if TESTS[kind][2] == 1:
        G = A.genome_pair(study, (a, b), u["pts"])
    else:
        G = subset_table(study)
        G = G[G.donor.isin(u["pts"]) & G.compartment.isin([a, b])].reset_index(drop=True)
    ref = G[G.compartment == b]
    genes = [c for c in G.columns if c not in META]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(L.q.WF.genes))]
    rng = np.random.default_rng(seed_of(kind, study) + [44])
    vals = [A.raw_stat(p["za"], p["zb"]) if (p := L.pair_states(G, a, b, wf=w, pts=u["pts"])) else np.nan
            for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF)]
    return study, kind, np.array(vals)


def stage_tests(args):
    pw = pd.read_csv(args.out / "power.csv").set_index("kind")
    cohorts = A.cohort_list(args.smoke)
    kinds = [k for k in KINDS if HYP[k][1] != "primary" or bool(pw.loc[k, "scored"]) or args.force]
    for k in KINDS:
        if k not in kinds:
            print(f"{k}: not scored (power note); use --force to run anyway", flush=True)
    us = {(s, k): units_of(s, k) for k in kinds for s in cohorts}
    us = {key: u for key, u in us.items() if u}
    ps = Parallel(n_jobs=args.jobs)(delayed(pseudo_one)(s, k, args.n_pseudo) for (s, k) in us)
    pseudo = {(s, k): v for s, k, v in ps}
    pd.DataFrame({f"{s}|{k}": v for (s, k), v in pseudo.items()}).to_csv(args.out / "pseudo.csv", index=False)
    per, pooled, tests = [], [], []

    def pool(sel, kind):
        w = np.array([us[x]["n"] for x in sel], float)
        obs = float(np.sum(w * [A.raw_stat(us[x]["za"], us[x]["zb"]) for x in sel]) / w.sum())
        nul = A.pool_null(w, [A.perm_null_raw(us[x], args.n_perm, np.random.default_rng(seed_of(kind, x[0]) + [1])) for x in sel])
        return w, obs, nul, A.pooled_pseudo(np.array([pseudo[x] for x in sel]), w)

    for kind in kinds:
        sel = [x for x in us if x[1] == kind]
        if not sel:
            continue
        w, obs, nul, pse = pool(sel, kind)
        p1, p2 = (1 + np.sum(nul >= obs)) / (1 + len(nul)), (1 + np.sum(pse >= obs)) / (1 + len(pse))
        tests.append(dict(kind=kind, hypothesis=HYP[kind][0], role=HYP[kind][1], part=TESTS[kind][2], patients=int(w.sum()),
                          cohorts=len(sel), statistic=obs, null_mean=float(nul.mean()), null_sd=float(nul.std()), p_pairing=float(p1),
                          pseudo_median=float(np.median(pse)), pseudo_q95=float(np.quantile(pse, .95)), n_sets=len(pse),
                          p_clock=float(p2), clock_minus_pseudo_median=obs - float(np.median(pse))))
        pooled.append(tests[-1].copy())
        for x in sel:
            per.append(dict(kind=kind, cohort=x[0], cancer=A.cancer(x[0]), n=us[x]["n"],
                            statistic=A.raw_stat(us[x]["za"], us[x]["zb"])))
        for ca in sorted({A.cancer(x[0]) for x in sel}):          # descriptive: cancer types with >= 12 patients
            sc = [x for x in sel if A.cancer(x[0]) == ca]
            if sum(us[x]["n"] for x in sc) >= MIN_CANCER:
                wc, oc, _, pc = pool(sc, kind)
                pooled.append(dict(kind=kind, cancer=ca, part=TESTS[kind][2], patients=int(wc.sum()), cohorts=len(sc), statistic=oc,
                                   pseudo_median=float(np.median(pc)), clock_minus_pseudo_median=oc - float(np.median(pc))))
    Tt = pd.DataFrame(tests)
    Tt["p_pairing_holm"] = np.nan
    for part in (1, 2):                      # Holm over the scored primaries of each part
        m = (Tt.role == "primary") & (Tt.part == part)
        if m.any():
            Tt.loc[m, "p_pairing_holm"] = A.holm(Tt.loc[m, "p_pairing"])
    Tt["pass"] = (Tt.p_pairing_holm < 0.05) & (Tt.statistic > Tt.pseudo_q95)
    Tt.to_csv(args.out / "tests.csv", index=False)
    D = pd.DataFrame(pooled)
    D.to_csv(args.out / "pooled.csv", index=False)
    D.sort_values("clock_minus_pseudo_median", ascending=False).to_csv(args.out / "descriptive.csv", index=False)
    pd.DataFrame(per).to_csv(args.out / "per_cohort.csv", index=False)
    print(Tt.round(4).to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "tests"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--force", action="store_true", help="tests: also run primaries flagged not scored")
    ap.add_argument("--smoke", action="store_true", help="3 cohorts, 200 permutations, 5 pseudo sets")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    assert args.jobs <= 8, "at most 8 workers"
    cli.setup(args)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.smoke:
        args.n_perm, args.n_pseudo = 200, 5
    {"power": stage_power, "tests": stage_tests}[args.stage](args)


if __name__ == "__main__":
    main()
