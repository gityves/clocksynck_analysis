"""A227: which tumour-microenvironment compartments keep the same patient's normal-tissue clock state?

Registered in paper/docs/registered_experiments.md (A227). A149 cohorts (the 27 with estimable S1) plus the A226 cohorts; tables,
labels, states and pseudo-clock sets as A149 / A151 / A226 (`a149_lib`, `compute_a151.genome_pair`, `compute_a226.table`).
Compartments: myeloid, T_NK (primary), B_plasma, mural, mast (secondary); fibroblast, endothelial are descriptive anchors.

  power  patients / cohorts per compartment, and the within-cohort permutation null of the uncorrected |complex correlation| on
         the real states (SHUFFLED pairings only: the identity pairing is never evaluated, no observed value is computed, stored or
         printed), pooled with patient-count weights: null mean, SD, MDE = 2.49 x SD, scored flag (>= 40 patients, MDE <= 0.25).
                                                                    -> output/{eligibility,eligibility_by_cancer,power}.csv
  tests  the registered statistic (IMPLEMENTED, NOT RUN before the power note is committed): observed pooled uncorrected
         |complex correlation|, one-sided permutation p, 200 expression-matched pseudo-clock sets (A149 note 2 (iv)), Holm over scored
         primaries and, separately, scored secondaries; descriptive clock - pseudo median overall, in cancer types with >= 12
         patients, and for the anchors.                             -> output/{tests,pooled,per_cohort,pseudo,descriptive}.csv
         Compartments flagged "not scored" in power.csv are skipped unless --force; anchors are always descriptive.

    python compute_a227.py --stage power|tests [--jobs J] [--n-perm 10000] [--n-pseudo 200] [--force] [--smoke --out DIR]
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
for p in ("A149_paired_tumour_normal", "A151_fibroblast_pairing", "A226_stroma_pairing_confirm"):
    sys.path.insert(0, str(PAPER / "exploratory" / p / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
import compute_a151 as F  # noqa: E402
import compute_a226 as K  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
SEED = 227
HYP = {"myeloid": ("H227.1", "primary"), "T_NK": ("H227.2", "primary"), "B_plasma": ("H227.3", "secondary"),
       "mural": ("H227.4", "secondary"), "mast": ("H227.5", "secondary"),
       "fibroblast": (None, "anchor"), "endothelial": (None, "anchor")}
KINDS = list(HYP)
PAIRS = {k: (f"tumour:{k}", f"adjacent_normal:{k}") for k in KINDS}
MIN_PATIENTS, MDE_K, MAX_MDE, MIN_CANCER = 40, 2.49, 0.25, 12
A149_SCORE = PAPER / "exploratory" / "A149_paired_tumour_normal" / "output" / "score" / "per_cohort.csv"


def cohort_list(smoke=False):
    """The 27 A149 cohorts (estimable S1 in A149 `per_cohort.csv`) and the 7 A226 cohorts."""
    a = pd.read_csv(A149_SCORE)
    c = a[(a.kind == "s1") & a.estimable].cohort.tolist() + list(K.COHORTS)
    return c[:2] + c[-1:] if smoke else c


def table(study):
    """Pseudobulk table of a cohort. Existing A149 cohorts: stored `cancer_sc` tables (they carry all compartments); new A149
    cohorts: cached A149 builds (`output/pb`); A226 cohorts: cached A226 builds. All A149 rules (>= 50 cells, median UMI >= 1,500,
    compartment kept if in >= 60 % of patients or >= 6)."""
    if study in K.COHORTS:
        return K.table(study)
    return L.table(study)[0]


def cancer(study):
    """Cancer type (A149 rule; `a149_lib.cancer_of` is overridden by `compute_a226.install()`, so repeated here)."""
    if study in K.COHORTS:
        return K.COHORTS[study]
    c = L.NEW.get(study) or (L.ST.GEO_CANCER.get(study, "") if study.startswith("GSE") else L.B._inventory_cancer(study))
    return L.CANCER_ALIAS.get(c, c)


# ---- eligibility -----------------------------------------------------------------------------------------------

def eligibility(study):
    from clocksync.config import limit_threads
    limit_threads(1)
    pb = table(study)
    rows = []
    for kind, (a, b) in PAIRS.items():
        A_, B_ = pb[pb.compartment == a], pb[pb.compartment == b]
        both = sorted(set(A_.donor) & set(B_.donor))
        sel = pd.concat([A_[A_.donor.isin(both)], B_[B_.donor.isin(both)]])
        rows.append(dict(study=study, cancer=cancer(study), kind=kind, patients_total=pb.donor.nunique(),
                         patients_tumour=A_.donor.nunique(), patients_normal=B_.donor.nunique(), patients_both=len(both),
                         eligible=len(both) >= L.MIN_PAIRED,
                         median_cells_per_pseudobulk=float(sel.n_cells.median()) if len(sel) else np.nan))
    return rows


# ---- statistic and nulls ---------------------------------------------------------------------------------------

def raw_stat(za, zb):
    """Uncorrected |complex correlation| of centred states (A226 `raw_stat`). Used only by the tests stage."""
    za, zb = za - za.mean(), zb - zb.mean()
    return float(abs(np.sum(za * np.conj(zb))) / np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zb) ** 2)))


def perm_null_raw(u, n_perm, rng):
    """Null of the statistic: rows of the normal states shuffled across patients (n_perm uniform permutations), vectorised."""
    za, zb = u["za"] - u["za"].mean(), u["zb"] - u["zb"].mean()
    perms = np.argsort(rng.random((n_perm, len(za))), axis=1)
    den = np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zb) ** 2))
    return np.abs((za[None, :] * np.conj(zb[perms])).sum(1)) / den


def seed_of(kind, study):
    return [SEED, KINDS.index(kind), sum(map(ord, study))]


def units_of(study, kind):
    pb = table(study)
    a, b = PAIRS[kind]
    if a not in set(pb.compartment) or b not in set(pb.compartment):
        return None
    return L.pair_states(pb, a, b)


def power_one(study, n_perm):
    """Permutation nulls (shuffled pairings only) of every compartment of one cohort: patients and null draws."""
    from clocksync.config import limit_threads
    limit_threads(1)
    out = {}
    for kind in KINDS:
        u = units_of(study, kind)
        if u:
            out[kind] = (u["n"], perm_null_raw(u, n_perm, np.random.default_rng(seed_of(kind, study))))
    return study, out


def pool_null(ns, nulls):
    w = np.asarray(ns, float)
    return (w[:, None] * np.asarray(nulls)).sum(0) / w.sum()


def stage_power(args):
    OUT.mkdir(parents=True, exist_ok=True)
    cohorts = cohort_list(args.smoke)
    elig = pd.DataFrame([r for x in Parallel(n_jobs=args.jobs)(delayed(eligibility)(s) for s in cohorts) for r in x])
    elig.to_csv(args.out / "eligibility.csv", index=False)
    e = elig[elig.eligible]
    by = e.groupby(["cancer", "kind"]).agg(patients=("patients_both", "sum"), cohorts=("study", "nunique")).reset_index()
    by.to_csv(args.out / "eligibility_by_cancer.csv", index=False)
    res = dict(Parallel(n_jobs=args.jobs)(delayed(power_one)(s, args.n_perm) for s in cohorts))
    rows = []
    for kind in KINDS:
        cs = [s for s in cohorts if kind in res[s]]
        ns = [res[s][kind][0] for s in cs]
        nul = pool_null(ns, [res[s][kind][1] for s in cs]) if cs else np.array([np.nan])
        sd = float(nul.std())
        n = int(sum(ns))
        rows.append(dict(kind=kind, hypothesis=HYP[kind][0], role=HYP[kind][1], patients=n, cohorts=len(cs),
                         cohort_list=";".join(cs), null_mean=float(nul.mean()), null_sd=sd, mde=MDE_K * sd,
                         n_perm=args.n_perm, scored=bool(n >= MIN_PATIENTS and MDE_K * sd <= MAX_MDE)))
    P = pd.DataFrame(rows)
    P.to_csv(args.out / "power.csv", index=False)
    print(P.drop(columns="cohort_list").round(4).to_string(), flush=True)


# ---- tests (registered statistic; not run before the power note) ----------------------------------------------

def genome_pair(study, comps, pts):
    K.install()
    return F.genome_pair(study, comps, pts)


def pseudo_one(study, kind, n_pseudo):
    """Expression-matched pseudo-clock sets of the cohort (A149 note 2 (iv), as A226), uncorrected statistic."""
    from clocksync.config import limit_threads
    from clocksync.paired import pseudo_waveforms
    limit_threads(1)
    a, b = PAIRS[kind]
    u = L.pair_states(table(study), a, b)
    G = genome_pair(study, (a, b), u["pts"])
    ref = G[G.compartment == b]
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth")]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(L.q.WF.genes))]
    rng = np.random.default_rng(seed_of(kind, study) + [44])
    vals = [raw_stat(p["za"], p["zb"]) if (p := L.pair_states(G, a, b, wf=w, pts=u["pts"])) else np.nan
            for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF)]
    return study, kind, np.array(vals)


def holm(p):
    p = np.asarray(p, float)
    h, run = np.empty(len(p)), 0.0
    for r, i in enumerate(np.argsort(p)):
        run = max(run, min(1.0, (len(p) - r) * p[i]))
        h[i] = run
    return h


def pooled_pseudo(P, w):
    Wm = np.where(np.isfinite(P), w[:, None], 0.0)
    pse = np.nansum(P * Wm, 0) / Wm.sum(0)
    return pse[np.isfinite(pse)]


def stage_tests(args):
    pw = pd.read_csv(args.out / "power.csv").set_index("kind")
    cohorts = cohort_list(args.smoke)
    kinds = [k for k in KINDS if HYP[k][1] == "anchor" or bool(pw.loc[k, "scored"]) or args.force]
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
        obs = float(np.sum(w * [raw_stat(us[x]["za"], us[x]["zb"]) for x in sel]) / w.sum())
        nul = pool_null(w, [perm_null_raw(us[x], args.n_perm, np.random.default_rng(seed_of(kind, x[0]) + [1])) for x in sel])
        pse = pooled_pseudo(np.array([pseudo[x] for x in sel]), w)
        return w, obs, nul, pse

    for kind in kinds:
        sel = [x for x in us if x[1] == kind]
        w, obs, nul, pse = pool(sel, kind)
        p1, p2 = (1 + np.sum(nul >= obs)) / (1 + len(nul)), (1 + np.sum(pse >= obs)) / (1 + len(pse))
        tests.append(dict(kind=kind, hypothesis=HYP[kind][0], role=HYP[kind][1], patients=int(w.sum()), cohorts=len(sel), statistic=obs,
                          null_mean=float(nul.mean()), null_sd=float(nul.std()), p_pairing=float(p1), pseudo_median=float(np.median(pse)),
                          pseudo_q95=float(np.quantile(pse, .95)), n_sets=len(pse), p_clock=float(p2), clock_minus_pseudo_median=obs - float(np.median(pse))))
        pooled.append(tests[-1].copy())
        for x in sel:
            per.append(dict(kind=kind, cohort=x[0], cancer=cancer(x[0]), n=us[x]["n"], statistic=raw_stat(us[x]["za"], us[x]["zb"])))
        for ca in sorted({cancer(x[0]) for x in sel}):          # descriptive: cancer types with >= 12 patients
            sc = [x for x in sel if cancer(x[0]) == ca]
            if sum(us[x]["n"] for x in sc) >= MIN_CANCER:
                wc, oc, _, pc = pool(sc, kind)
                pooled.append(dict(kind=kind, cancer=ca, patients=int(wc.sum()), cohorts=len(sc), statistic=oc,
                                   pseudo_median=float(np.median(pc)), clock_minus_pseudo_median=oc - float(np.median(pc))))
    T = pd.DataFrame(tests, columns=["kind", "hypothesis", "role", "statistic", "p_pairing", "pseudo_q95"] if not tests else None)
    T["p_pairing_holm"] = np.nan
    for role in ("primary", "secondary"):
        m = T.role == role
        if m.any():
            T.loc[m, "p_pairing_holm"] = holm(T.loc[m, "p_pairing"])
    T["pass"] = (T.p_pairing_holm < 0.05) & (T.statistic > T.pseudo_q95)
    T.to_csv(args.out / "tests.csv", index=False)
    D = pd.DataFrame(pooled)
    D.to_csv(args.out / "pooled.csv", index=False)
    D.sort_values("clock_minus_pseudo_median", ascending=False).to_csv(args.out / "descriptive.csv", index=False)
    pd.DataFrame(per).to_csv(args.out / "per_cohort.csv", index=False)
    print(T.round(4).to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "tests"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--force", action="store_true", help="tests: also run compartments flagged not scored")
    ap.add_argument("--smoke", action="store_true", help="3 cohorts, 200 permutations")
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
