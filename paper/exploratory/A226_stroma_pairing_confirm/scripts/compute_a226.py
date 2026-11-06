"""A226: does tumour stroma follow the same patient's normal-tissue stroma clock-specifically, in cohorts not used before?

Registered in paper/docs/registered_experiments.md (A226). New paired tumour / adjacent-normal cohorts (`cohorts.csv`, fixed in
the power note), built by `build_a226.py` ($SCCC_DATA/a226/<study>/<study>.h5ad), pseudobulks with the A149 / `cancer_sc` rules
(marker compartments; the authors' labels are recorded in the h5ad but not used, as in A149), states and pseudo-clock sets as
A151 (`a149_lib.pair_states`, `compute_a151.genome_pair`).

  power  patients / depth per compartment and cohort, and the simulated null SD, MDE and scored flag. NO state, pairing
         statistic or permutation of observed states is computed.        -> output/{eligibility,power}.csv
  tests  the registered statistic (IMPLEMENTED, NOT RUN before the power note is committed): uncorrected |complex correlation|
         |sum za conj(zb)| / sqrt(sum |za|^2 sum |zb|^2) of the centred states, pooled with patient-count weights, within-cohort
         permutation null (normal states shuffled, 10,000), 200 expression-matched pseudo-clock sets (A149 note 2 (iv)); Holm over
         H226.1-2 (fibroblasts) and H226.3-4 (endothelial) separately.   -> output/{tests,pooled,per_cohort,pseudo}.csv
         A compartment flagged "not scored" in output/power.csv is skipped unless --force.

    python compute_a226.py --stage power|tests [--jobs J] [--n-perm 10000] [--n-pseudo 200] [--force]
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER / "exploratory" / "A149_paired_tumour_normal" / "scripts"))
sys.path.insert(0, str(PAPER / "exploratory" / "A151_fibroblast_pairing" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
import compute_a151 as F  # noqa: E402
from common import cli  # noqa: E402

ST = L.ST
A226 = ST.DATA / "a226"
OUT = HERE.parent / "output"
SEED = 226
PAIRS = F.PAIRS
# cohort: (cancer type, source description); fixed in the power note (scout/SCOUT.md)
COHORTS = {"GSE308103": "lung", "GSE183219": "lung", "GSE155109": "breast", "OPSCC_b6b5ea88": "head and neck",
           "RCC_5af90777": "kidney", "Qian2020_CRC": "colorectal", "Qian2020_LUAD": "lung"}
POSTHOC = PAPER / "exploratory" / "A151_fibroblast_pairing" / "output" / "posthoc_pooled.csv"
MIN_PATIENTS, MDE_K = 30, 2.49


# ---- cohort tables (A149 conventions, a226 root) ---------------------------------------------------------------------

def _units(study, where=None):
    """`cancer_common.studies.units` for the a226 h5ad files (same filters: patient known, not treated, not CD45+-sorted, sample
    median UMI >= 1,500; sample types tumour / adjacent_normal)."""
    from clocksync.genes import read_counts
    a = read_counts(A226 / study / f"{study}.h5ad")
    o = a.obs
    ok = (o.sample_type.isin(ST.KEEP) & (o.treatment != "treated") & ~o.sorting.str.contains(r"CD45\+|CD45p", regex=True)
          & o.patient.astype(str).ne(""))
    a = a[ok.to_numpy()].copy()
    depth = a.obs.groupby("sample", observed=True).total_umi.median()
    a = a[a.obs["sample"].isin(depth[depth >= (ST.MIN_UMI if SENS[0] is None else SENS[0])].index).to_numpy()].copy()
    sym = a.var_names.astype(str).to_numpy()
    comp = ST.marker_compartments(a.X, sym, a.obs.total_umi.to_numpy(float))
    chem = a.obs["chemistry"].astype(str) if "chemistry" in a.obs else pd.Series("", index=a.obs.index)
    out = []
    for smp, g in a.obs.groupby("sample", observed=True):
        i = np.flatnonzero(a.obs["sample"].to_numpy() == smp)
        out.append((str(smp), str(g.patient.iloc[0]), str(g.sample_type.iloc[0]),
                    (lambda i=i: (a.X[i], a.var_names, a.obs.total_umi.to_numpy(float)[i], comp[i], a.obs.index[i], sym)),
                    chem.iloc[i[0]], len(i)))
    return out, COHORTS[study]


_orig_units = ST.units
SENS = [None]       # sample median-UMI limit override (sensitivity only; None = the A149 rule, 1,500)


def install():
    """Route `cancer_common.studies.units` of the a226 cohorts to `_units` (idempotent; called in every worker process)."""
    ST.units = lambda study, where=None: _units(study, where) if study in COHORTS else _orig_units(study, where)  # noqa: E731
    L.cancer_of = lambda s: COHORTS[s]  # noqa: E731


install()


def table(study, sens=False):
    """Pseudobulk table with the A149 rules (>= 50 cells, sample median UMI >= 1,500, compartment kept if in >= 60 % of patients or
    >= 6); `sens=True`: sensitivity table without the median-UMI limit and without the compartment presence rule."""
    install()
    f = OUT / "pb" / f"{study}{'_sens' if sens else ''}.parquet"
    if f.exists():
        return pd.read_parquet(f)
    SENS[0] = 0 if sens else None
    try:
        pb, _ = L.B.build(study, cnv=False, compartment_rule=not sens)
    finally:
        SENS[0] = None
    f.parent.mkdir(parents=True, exist_ok=True)
    pb.to_parquet(f)
    return pb


# ---- power note (no pairing statistic) ---------------------------------------------------------------------------

def eligibility(study):
    """Patients per compartment (>= 50 cells per pseudobulk, sample median UMI >= 1,500, A149 rules) and depth. Counts only."""
    pb, sb = table(study), table(study, sens=True)
    rows = []
    for kind, (a, b) in PAIRS.items():
        A_, B_ = pb[pb.compartment == a], pb[pb.compartment == b]
        both = sorted(set(A_.donor) & set(B_.donor))
        sboth = len(set(sb[sb.compartment == a].donor) & set(sb[sb.compartment == b].donor))
        sel = pd.concat([A_[A_.donor.isin(both)], B_[B_.donor.isin(both)]])
        rows.append(dict(study=study, cancer=COHORTS[study], kind=kind, patients_total=pb.donor.nunique(),
                         patients_tumour=A_.donor.nunique(), patients_normal=B_.donor.nunique(), patients_both=len(both),
                         patients_both_sens=sboth, patients_total_sens=sb.donor.nunique(),
                         eligible=len(both) >= L.MIN_PAIRED, median_cells_per_pseudobulk=float(sel.n_cells.median()) if len(sel) else np.nan,
                         median_depth_per_pseudobulk=float(sel.depth.median()) if len(sel) else np.nan))
    return rows


def raw_stat(za, zb):
    za, zb = za - za.mean(), zb - zb.mean()
    return float(abs(np.sum(za * np.conj(zb))) / np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zb) ** 2)))


def sim_null(ns, n_sim, rng):
    """Pooled uncorrected |complex correlation| of independent isotropic complex Gaussian states, cohorts of sizes `ns`
    (states centred within cohort), weights = patients: its null distribution under no pairing."""
    out = np.empty(n_sim)
    w = np.asarray(ns, float)
    for s in range(n_sim):
        out[s] = np.sum(w * [raw_stat(rng.standard_normal(n) + 1j * rng.standard_normal(n),
                                      rng.standard_normal(n) + 1j * rng.standard_normal(n)) for n in ns]) / w.sum()
    return out


def power(elig, n_sim=20_000):
    """Null SD from simulated independent states with the cohorts' real patient counts, rescaled by the A151 calibration: the same
    simulation at the A151 cohort sizes (`A151 per_cohort.csv`) is compared with the spread of the A151 pseudo-clock sets
    (`posthoc_pooled.csv`: (q95 - median) / 1.645, the empirical null of the same statistic on the same patients). Only counts and
    published A151 summaries enter; no state of the new cohorts is computed."""
    ref = pd.read_csv(POSTHOC).set_index("kind")
    a151 = pd.read_csv(F.OUT / "per_cohort.csv")
    rows = []
    for kind in PAIRS:
        e = elig[(elig.kind == kind) & elig.eligible]
        ns = e.patients_both.astype(int).tolist()
        ns151 = a151[(a151.kind == kind) & a151.estimable].n.astype(int).tolist()
        rng = np.random.default_rng([SEED, len(kind)])
        s151, snew = sim_null(ns151, n_sim, rng), (sim_null(ns, n_sim, rng) if ns else np.array([np.nan]))
        sd_ref = (ref.loc[kind, "pseudo_q95"] - ref.loc[kind, "pseudo_median"]) / 1.645
        cal = sd_ref / s151.std()
        diff = ref.loc[kind, "clock_raw"] - ref.loc[kind, "pseudo_median"]
        sd_raw = float(snew.std())
        sd = sd_raw * cal
        es = elig[(elig.kind == kind) & (elig.patients_both_sens >= L.MIN_PAIRED)]     # sensitivity: no depth / presence rule
        nss = es.patients_both_sens.astype(int).tolist()
        sd_s = float(sim_null(nss, n_sim, rng).std() * cal) if nss else np.nan
        rows.append(dict(kind=kind, patients=int(sum(ns)), cohorts=len(ns), cohort_list=";".join(e.study), sim_null_mean=float(np.mean(snew)),
                         sim_null_sd_raw=sd_raw, a151_sim_sd=float(s151.std()), a151_sim_mean=float(s151.mean()),
                         a151_pseudo_median=float(ref.loc[kind, "pseudo_median"]), a151_pseudo_sd=float(sd_ref), calibration=float(cal),
                         null_sd=float(sd), mde=float(MDE_K * sd), a151_clock=float(ref.loc[kind, "clock_raw"]), a151_pseudo_q95=float(ref.loc[kind, "pseudo_q95"]),
                         a151_diff_clock_minus_pseudo_median=float(diff), half_diff=float(diff / 2),
                         mde_if_calibration_1=float(MDE_K * sd_raw), n_sim=n_sim,
                         sens_patients=int(sum(nss)), sens_cohorts=len(nss), sens_cohort_list=";".join(es.study), sens_null_sd=sd_s, sens_mde=float(MDE_K * sd_s) if nss else np.nan,
                         scored=bool(sum(ns) >= MIN_PATIENTS and MDE_K * sd <= diff / 2)))
    return pd.DataFrame(rows)


def stage_power(args):
    OUT.mkdir(parents=True, exist_ok=True)
    res = Parallel(n_jobs=args.jobs)(delayed(eligibility)(s) for s in COHORTS)
    elig = pd.DataFrame([r for x in res for r in x])
    elig.to_csv(OUT / "eligibility.csv", index=False)
    print(elig.round(1).to_string(), flush=True)
    P = power(elig)
    P.to_csv(OUT / "power.csv", index=False)
    print(P.round(4).T.to_string(), flush=True)


# ---- tests (registered statistic; not run before the power note) ---------------------------------------------------

def units_of(study, kind):
    pb = table(study)
    a, b = PAIRS[kind]
    if a not in set(pb.compartment) or b not in set(pb.compartment):
        return None
    return L.pair_states(pb, a, b)


def perm_null_raw(u, n_perm, rng):
    za, zb = u["za"] - u["za"].mean(), u["zb"] - u["zb"].mean()
    n = len(za)
    out = np.empty(n_perm)
    den = np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zb) ** 2))
    for k in range(n_perm):
        out[k] = abs(np.sum(za * np.conj(zb[rng.permutation(n)]))) / den
    return out


def pseudo_one(study, kind, n_pseudo):
    """200 expression-matched pseudo-clock sets of the cohort (A149 note 2 (iv), A151 `pseudo_one`), uncorrected statistic."""
    from clocksync.config import limit_threads
    from clocksync.paired import pseudo_waveforms
    limit_threads(1)
    a, b = PAIRS[kind]
    u = L.pair_states(table(study), a, b)
    G = F.genome_pair(study, (a, b), u["pts"])
    ref = G[G.compartment == b]
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth")]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(L.q.WF.genes))]
    rng = np.random.default_rng([SEED, 44, len(kind), sum(map(ord, study))])
    vals = [raw_stat(p["za"], p["zb"]) if (p := L.pair_states(G, a, b, wf=w, pts=u["pts"])) else np.nan
            for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF)]
    return study, np.array(vals)


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    h = np.empty(len(p))
    run = 0.0
    for r, i in enumerate(o):
        run = max(run, min(1.0, (len(p) - r) * p[i]))
        h[i] = run
    return h


def stage_tests(args):
    pw = pd.read_csv(OUT / "power.csv").set_index("kind")
    tests, pooled, per = [], [], []
    for (kind, hyps) in (("fibroblast", ("H226.1", "H226.2")), ("endothelial", ("H226.3", "H226.4"))):
        if not bool(pw.loc[kind, "scored"]) and not args.force:
            print(f"{kind}: not scored (power note); use --force to run anyway", flush=True)
            continue
        cohorts = [s for s in pw.loc[kind, "cohort_list"].split(";") if s]
        us = {s: units_of(s, kind) for s in cohorts}
        us = {s: u for s, u in us.items() if u}
        w = np.array([u["n"] for u in us.values()], float)
        obs_c = np.array([raw_stat(u["za"], u["zb"]) for u in us.values()])
        nulls = np.array([perm_null_raw(u, args.n_perm, np.random.default_rng([SEED, len(kind), sum(map(ord, s))])) for s, u in us.items()])
        nul = (w[:, None] * nulls).sum(0) / w.sum()
        obs = float((w * obs_c).sum() / w.sum())
        p1 = float((1 + np.sum(nul >= obs)) / (1 + len(nul)))
        ps = dict(Parallel(n_jobs=args.jobs)(delayed(pseudo_one)(s, kind, args.n_pseudo) for s in us))
        P = np.array([ps[s] for s in us])
        Wm = np.where(np.isfinite(P), w[:, None], 0.0)
        pse = np.nansum(P * Wm, 0) / Wm.sum(0)
        ok = np.isfinite(pse)
        p2 = float((1 + np.sum(pse[ok] >= obs)) / (1 + ok.sum()))
        h = holm([p1, p2])
        tests += [dict(hypothesis=hyps[0], kind=kind, statistic=obs, null_mean=float(nul.mean()), null_sd=float(nul.std()), p=p1, p_holm=h[0]),
                  dict(hypothesis=hyps[1], kind=kind, statistic=obs, pseudo_median=float(np.median(pse[ok])), pseudo_q95=float(np.quantile(pse[ok], .95)),
                       n_sets=int(ok.sum()), p=p2, p_holm=h[1])]
        pooled.append(dict(kind=kind, n=int(w.sum()), n_cohorts=len(us), statistic=obs))
        per += [dict(kind=kind, cohort=s, n=u["n"], statistic=float(o), null_mean=float(nl.mean()))
                for (s, u), o, nl in zip(us.items(), obs_c, nulls)]
        pd.DataFrame(ps).to_csv(OUT / f"pseudo_{kind}.csv", index=False)
    T = pd.DataFrame(tests)
    if len(T):
        T["pass"] = T.p_holm < 0.05
        T.to_csv(OUT / "tests.csv", index=False)
        pd.DataFrame(pooled).to_csv(OUT / "pooled.csv", index=False)
        pd.DataFrame(per).to_csv(OUT / "per_cohort.csv", index=False)
        print(T.round(4).to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "tests"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--force", action="store_true", help="tests: also run compartments flagged not scored")
    args = ap.parse_args(argv)
    cli.setup(args)
    {"power": stage_power, "tests": stage_tests}[args.stage](args)


if __name__ == "__main__":
    main()
