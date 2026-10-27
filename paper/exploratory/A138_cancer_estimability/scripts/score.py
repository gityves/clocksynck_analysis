"""A138 scoring of Q3 and Q1: H138.4, H138.5 (descriptive), H138.1 and D138.1.

Registered in paper/docs/registered_experiments.md (A138; notes 3, 4, 5). Uses the functions, data loading, weights,
cohorts and seeds of the power step (`power.py`, `q3_lib.py`): the permutation nulls are regenerated with the same
random streams and checked against the committed `output/power/*_null.csv` (null mean, SD, q95).
- H138.4: pooled (weights n_pairs, gated cancers) |complex correlation| of tumour and own-normal states; p =
  (1 + #{null >= obs}) / (10,000 + 1); companion: 200 expression-matched pseudo-clock sets (A52 rule), pooled 95th pct.
- H138.5 (descriptive, out of the Holm family): host-referenced malignant vs normal epithelium states; cohorts with
  >= 6 paired patients; note 4 (ii) data rule; companion: tumour vs adjacent fibroblasts (own permutations).
- H138.1: anchored coupling of the malignant compartment per Q1 cohort (the 17 cohorts of `h138_1_designs.csv`; learned
  anchors without tumour epithelium as the power simulation), 1,000 patient bootstraps, DerSimonian-Laird pooled.
  Note 5 rule: pass if pooled > 0.083; p = (1 + k) / 31 over the 30 c = 0 replicates of `h138_1_sim_cohorts.csv`
  (pooled estimates recomputed as `power.py`); the registered lower-bound rule is reported next to it.
- D138.1: per cohort variance split of malignant, normal epithelium and each stromal compartment of the same patients.
- Holm: over H138.1 and H138.4 (provisional: H138.2, H138.3, H138.6, H138.7 not yet scored).

    python score.py [--subsample N] [--jobs J] [--seed S]
Outputs (paper/exploratory/A138_cancer_estimability/output/score/, `subsample/` for checks): tests.csv, per_unit.csv,
holm.csv, d138_1.csv
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
import q3_lib as q  # noqa: E402
from common import cli  # noqa: E402
from power import SEED  # noqa: E402  (138: the stream of the power step)

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "output" / "score"
POWER = HERE / "output" / "power"
GROUPS = q.REPO / "paper" / "results" / "figure_2" / "output" / "backbone" / "a18_cancer_groups.csv"
T0 = time.time()
THRESH_H138_1 = 0.083  # note 5 (i)
N_BOOT = 1000
N_PSEUDO = 200


def log(*a):
    print(f"[{(time.time() - T0) / 60:7.1f} min]", *a, flush=True)


def perm_p(null: np.ndarray, obs: float) -> float:
    x = null[np.isfinite(null)]
    return float((1 + np.sum(x >= obs)) / (len(x) + 1))


def check_null(null: np.ndarray, ref: pd.Series, what: str, quick: bool):
    """The regenerated null must equal the committed power null (summary level)."""
    s = q.summarise(null)
    if quick:
        return
    for k in ("null_mean", "null_sd", "null_q95"):
        if not np.isclose(s[k], ref[k], rtol=1e-9, atol=1e-12):
            raise RuntimeError(f"{what}: regenerated null differs from the committed one ({k}: {s[k]} vs {ref[k]})")


# ---- H138.4 ------------------------------------------------------------------------------------------------------

def pseudo_pooled(cancers, st, samples, rng_seed: int) -> np.ndarray:
    """Pooled (weights n_pairs) statistic of `N_PSEUDO` expression-matched pseudo-clock gene sets (A52 rule: genes of
    similar mean log fraction in the normals, carrying the waveform's amplitudes and phases; `pseudo_waveforms`)."""
    from clocksync.paired import cohort_states, pseudo_waveforms
    per = []
    for i, c in enumerate(cancers):
        rng = np.random.default_rng([SEED, 44, rng_seed, i])
        case, _ = q.sh.pair_frames(samples, c)
        ct, cn = q.sh.tcga_full_counts(c, case.index)
        lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
        dt, dn = ct.sum(axis=1), cn.sum(axis=1)
        r = []
        for w in pseudo_waveforms(lm, N_PSEUDO, rng, q.WF):
            g = list(w.genes)
            bt, vt, bn, vn = cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)
            r.append(q.corr_stat(q.centred_complex(bt), vt, q.centred_complex(bn), vn))
        per.append(r)
        log(f"pseudo sets {c}")
    return q.pool(np.array(per), [st[c]["n"] for c in cancers])


def score_h138_4(args, quick: bool, n_perm: int, n_pseudo_on: bool):
    st = q.tcga_states(args.subsample, args.seed)
    ref = pd.read_csv(POWER / ("subsample/" if quick else "") / "h138_4_null.csv").set_index("scope")
    obs, nulls, rows = {}, {}, []
    for i, (c, u) in enumerate(st.items()):
        obs[c] = q.corr_stat(u["za"], u["vA"], u["zb"], u["vB"])
        nulls[c] = q.perm_null(u["za"], u["vA"], u["zb"], u["vB"], n_perm, np.random.default_rng([SEED, 4, args.seed, i]))
        check_null(nulls[c], ref.loc[c], f"H138.4 {c}", quick)
        rows.append(dict(hypothesis="H138.4", unit=c, n=u["n"], gate_p=u["gate_p"], gated=u["gated"], statistic=obs[c],
                         null_mean=float(np.nanmean(nulls[c])), p=perm_p(nulls[c], obs[c]) if np.isfinite(obs[c]) else np.nan))
    keep = [c for c in st if np.isfinite(nulls[c]).any() and st[c]["gated"]]
    w = [st[c]["n"] for c in keep]
    pooled_obs = float(q.pool(np.array([[obs[c]] for c in keep]), w)[0])
    pooled_null = q.pool(np.array([nulls[c] for c in keep]), w)
    check_null(pooled_null, ref.loc["pooled_gated"], "H138.4 pooled", quick)
    p = perm_p(pooled_null, pooled_obs)
    t = dict(hypothesis="H138.4", n_units=sum(w), n_cohorts=len(keep), statistic=pooled_obs, p=p, null_mean=float(pooled_null.mean()),
             null_sd=float(pooled_null.std(ddof=1)), null_q95=float(np.quantile(pooled_null, 0.95)), n_perm=len(pooled_null),
             verdict_rule="one-sided p < 0.05 (permutation)", pass_rule=bool(p < 0.05),
             per_unit=";".join(f"{c}={obs[c]:.4f}" for c in keep))
    if n_pseudo_on:
        samples = q.sh.load_tcga(args.subsample, args.seed)
        ps = pseudo_pooled(keep, st, samples, args.seed)
        t |= dict(companion="pseudo-clock pooled statistic (200 sets, A52 rule; descriptive)", companion_q95=float(np.nanquantile(ps, 0.95)),
                  companion_median=float(np.nanmedian(ps)), companion_share_ge_obs=float(np.mean(ps >= pooled_obs)))
    return t, rows


# ---- H138.5 ------------------------------------------------------------------------------------------------------

def one_study_obs(study: str, n_perm: int, i: int, seed: int):
    from clocksync.config import limit_threads
    limit_threads(1)
    pb, mal = q.study_table(study)
    u = q.sc_units(study, pb, mal)
    res = {}
    for part, x in u.items():
        if x is None:
            continue
        j = 0 if part == "primary" else 1
        null = q.perm_null(x["za"], x["vA"], x["zb"], x["vB"], n_perm, np.random.default_rng([SEED, 5, seed, i, j]))
        res[part] = dict(est=q.estimable(x["za"], x["vA"], x["zb"], x["vB"], null), n=x["n"], null=null,
                         obs=q.corr_stat(x["za"], x["vA"], x["zb"], x["vB"]), malignant=mal)
    return study, res


def score_h138_5(args, quick: bool, n_perm: int):
    studies = list(q.sh.sc_builder().STUDIES)
    if args.subsample:
        studies = studies[:6]
    gate = q.a17_gate()
    ref = pd.read_csv(POWER / ("subsample/" if quick else "") / "h138_5_null.csv")
    out = Parallel(n_jobs=args.jobs)(delayed(one_study_obs)(s, n_perm, i, args.seed) for i, s in enumerate(studies))
    rows, tests, keep = [], [], {"primary": [], "fibroblast": []}
    for study, res in out:
        for part, r in res.items():
            rr = ref[(ref.scope == study) & (ref.kind == part)]
            if len(rr) and not quick:
                check_null(r["null"], rr.iloc[0], f"H138.5 {study} {part}", quick)
                assert bool(rr.iloc[0].estimable) == r["est"]
            keep[part].append((study, r))
            rows.append(dict(hypothesis="H138.5" if part == "primary" else "H138.5 fibroblast", unit=study, n=r["n"], gate_p=np.nan,
                             gated=gate.get(study, False), estimable=r["est"], statistic=r["obs"], null_mean=float(np.nanmean(r["null"])),
                             p=perm_p(r["null"], r["obs"]) if np.isfinite(r["obs"]) else np.nan, malignant=r["malignant"] if part == "primary" else ""))
    for part, name in (("primary", "H138.5"), ("fibroblast", "H138.5 fibroblast companion")):
        lst = [x for x in keep[part] if x[1]["est"]]  # note 4 (ii)
        for scope, sel in (("pooled_all", lst), ("pooled_gate_passing", [x for x in lst if gate.get(x[0], False)])):
            if not sel:
                continue
            w = [r["n"] for _, r in sel]
            obs = float(q.pool(np.array([[r["obs"]] for _, r in sel]), w)[0])
            null = q.pool(np.array([r["null"] for _, r in sel]), w)
            rr = ref[(ref.scope == scope) & (ref.kind == part + "_pooled")]
            if len(rr) and not quick:
                check_null(null, rr.iloc[0], f"H138.5 {part} {scope}", quick)
            p = perm_p(null, obs)
            tests.append(dict(hypothesis=name + ("" if scope == "pooled_all" else " (gate-passing cohorts)"), n_units=sum(w), n_cohorts=len(sel),
                              statistic=obs, p=p, null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                              null_q95=float(np.quantile(null, 0.95)), n_perm=len(null), verdict_rule="descriptive (note 5 ii)",
                              pass_rule=np.nan, per_unit=";".join(f"{s}={r['obs']:.4f}" for s, r in sel),
                              detection_limit_mde=float(np.quantile(null, 0.95) + 0.84 * null.std(ddof=1)),
                              excluded_not_estimable=";".join(s for s, r in keep[part] if not r["est"])))
    return tests, rows


# ---- H138.1 and D138.1 --------------------------------------------------------------------------------------------

def stroma_comps(pb: pd.DataFrame) -> list:
    return sorted(set(pb.compartment[q.cls(pb.compartment).isin(q.STROMA)]))


def one_cohort(d: dict, ci: int, n_boot: int, seed: int):
    """Observed anchored fit of the Q1 cohort `d` on the patients with the malignant compartment. Fit 1 (H138.1): target
    = malignant, `n_boot` patient bootstraps (host refitted). Fit 2 (D138.1, point estimates): malignant, normal epithelium
    and the stromal compartments as targets; Poisson part = mean trace of the measurement covariance of the compartment."""
    from clocksync.config import limit_threads
    limit_threads(1)
    pb, mal = q.study_table(d["study"])
    pts = set(pb[pb.compartment == mal].donor)
    ss = int(np.random.SeedSequence([SEED, 6, seed, ci]).generate_state(1)[0])
    t1 = pb[pb.compartment.isin(d["anchors"] + [mal]) & pb.donor.isin(pts)]
    ds = q.cs.Dataset(t1, anchors=d["anchors"])
    f = ds.fit(n_boot=n_boot, seed=ss, targets=[mal], gate=False)
    rep = f.replicates[mal].to_numpy(float)
    rep = rep[np.isfinite(rep)]
    tab = f.table.loc[mal]
    h1 = dict(study=d["study"], malignant=mal, n_patients=int(tab.n), c=float(tab.coupling), se=float(np.std(rep, ddof=1)) if len(rep) > 2 else np.nan,
              lo95=float(tab.coupling_lo), hi95=float(tab.coupling_hi), n_boot=int(tab.n_boot), boundary=bool(tab.boundary), z2=float(tab.z2),
              psi=float(tab.psi), V=float(tab.V), n_anchors=len(d["anchors"]))
    comps = sorted(set(d["anchors"]) | {mal, q.NORM} | set(stroma_comps(pb)))
    t2 = pb[pb.compartment.isin(comps) & pb.donor.isin(pts)]
    ds2 = q.cs.Dataset(t2, anchors=d["anchors"])
    f2 = ds2.fit(n_boot=0, seed=ss, targets=[c for c in comps if c in set(t2.compartment)], gate=False)
    dd = ds2.data
    rows = []
    for k, r in f2.table.iterrows():
        i = np.flatnonzero(np.asarray(dd.compartment) == k)
        pois = float(np.mean(np.trace(np.asarray(dd.V)[i], axis1=1, axis2=2))) if len(i) else np.nan
        cl = q.cls(pd.Series([k])).iloc[0]
        role = "malignant" if k == mal else "normal_epithelium" if k == q.NORM else "stroma" if cl in q.STROMA else "other"
        rows.append(dict(study=d["study"], compartment=k, role=role, n_patients=int(r.n), c=r.coupling, host_coupled=r.z2 * r.V, patient_specific=2 * r.psi,
                         poisson=pois, boundary=bool(r.boundary), anchor=bool(r.anchor), z2=r.z2, psi=r.psi, V=r.V))
    return h1, rows


def between_cell(dd: pd.DataFrame, groups: pd.DataFrame, pts_by_study: dict) -> pd.DataFrame:
    """Between-cell within-patient dispersion on the clock plane: D of the A18 cell groups (>= 100 cells per patient x
    compartment, `a18_cancer_groups.csv`; D = trace of the within-group in-plane second moment, log scale, free of
    counting noise), cell-weighted over the groups of the same patients, raw and minus its acrophase-permuted null."""
    out = []
    for (s, k), g in dd.groupby(["study", "compartment"]):
        x = groups[(groups.study == s) & (groups.compartment == k) & groups.resection.isin(pts_by_study[s])]
        out.append(dict(study=s, compartment=k, between_cell_groups=len(x),
                        between_cell_raw=np.average(x.D, weights=x.cells) if len(x) else np.nan,
                        between_cell=np.average(x.D - x.D_null, weights=x.cells) if len(x) else np.nan))
    return pd.DataFrame(out)


def score_h138_1(args, quick: bool):
    gate = q.a17_gate()
    designs = [dd for s in sorted(s for s, g in gate.items() if g) if (dd := q.q1_design(s))]
    if args.subsample:
        designs = designs[:3]
    ref_d = pd.read_csv(POWER / ("subsample/" if quick else "") / "h138_1_designs.csv")
    if not quick:
        assert [d["study"] for d in designs] == list(ref_d.study), "H138.1 cohorts differ from the power step"
    n_boot = 10 if quick else N_BOOT
    res = Parallel(n_jobs=args.jobs)(delayed(one_cohort)(d, ci, n_boot, args.seed) for ci, d in enumerate(designs))
    coh = pd.DataFrame([r[0] for r in res])
    dec = pd.DataFrame([x for r in res for x in r[1]])
    m = q.S.meta(coh.c, coh.se)
    # null of the pooled estimate: the 30 c = 0 replicates of the power simulation, pooled as power.py does
    sim = pd.read_csv(POWER / ("subsample/" if quick else "") / "h138_1_sim_cohorts.csv")
    nul = np.array([q.S.meta(x.est, x.se)["random"] for _, x in sim[sim.c == 0].groupby("rep")])
    nul = nul[np.isfinite(nul)]
    if not quick:
        sm = pd.read_csv(POWER / "power_summary.csv").set_index("hypothesis").loc["H138.1"]
        assert np.isclose(nul.mean(), sm.null_mean) and np.isclose(nul.std(ddof=1), sm.null_sd) and np.isclose(np.quantile(nul, 0.95), sm.null_q95)
        assert len(nul) == 30
    k = int(np.sum(nul >= m["random"]))
    p = float((1 + k) / (len(nul) + 1))
    t = dict(hypothesis="H138.1", n_units=int(coh.n_patients.sum()), n_cohorts=int(m["k"]), statistic=m["random"], p=p, null_mean=float(nul.mean()),
             null_sd=float(nul.std(ddof=1)), null_q95=float(np.quantile(nul, 0.95)), n_perm=len(nul), verdict_rule=f"note 5: pooled > {THRESH_H138_1}; p = (1+k)/(n+1)",
             pass_rule=bool(m["random"] > THRESH_H138_1), per_unit=";".join(f"{r.study}={r.c:.3f}" for r in coh.itertuples()),
             pooled_lo95=m["lo"], pooled_hi95=m["hi"], pooled_se=m["se_random"], tau2=m["tau2"], pooled_fixed=m["fixed"], null_k_ge_obs=k,
             lower_bound_rule_pass=bool(m["lo"] > 0), n_boot=n_boot)
    # D138.1 table: variance split per cohort and compartment (+ between-cell dispersion where the cell groups exist)
    groups = pd.read_csv(GROUPS)
    pts_by_study = {}
    for d in designs:
        pb, mal = q.study_table(d["study"])
        pts_by_study[d["study"]] = set(pb[pb.compartment == mal].donor)
    dec = dec.merge(between_cell(dec, groups, pts_by_study), on=["study", "compartment"], how="left")
    return t, coh, dec


def d138_table(dec: pd.DataFrame) -> pd.DataFrame:
    """Shares of the host-coupled, patient-specific, between-cell and Poisson parts. Total A = host + patient + Poisson
    (the across-patient observed clock-plane variance); total B = A + between-cell (a within-patient part, on the
    per-cell scale). Stroma rows are aggregated over the stromal compartments (mean of parts)."""
    d = dec[dec.role.isin(["malignant", "normal_epithelium", "stroma"])].copy()
    s = d[d.role == "stroma"]
    cols = ["host_coupled", "patient_specific", "poisson", "between_cell", "between_cell_raw", "c"]
    agg = s.groupby("study")[cols + ["n_patients"]].mean().reset_index().assign(compartment="stroma (mean of compartments)", role="stroma_mean")
    d = pd.concat([d[d.role != "stroma"], s.assign(role="stroma"), agg], ignore_index=True)
    d["total_A"] = d.host_coupled + d.patient_specific + d.poisson
    d["total_B"] = d.total_A + d.between_cell.clip(lower=0)
    for c in ("host_coupled", "patient_specific", "poisson"):
        d[f"share_{c}_A"] = d[c] / d.total_A
        d[f"share_{c}_B"] = d[c] / d.total_B
    d["share_between_cell_B"] = d.between_cell.clip(lower=0) / d.total_B
    return d.sort_values(["study", "role", "compartment"]).reset_index(drop=True)


def contrast_patient_share(d: pd.DataFrame) -> pd.DataFrame:
    """Estimand of interest: patient-specific share, malignant minus normal epithelium, per cohort (both shares defined)."""
    w = d[d.role.isin(["malignant", "normal_epithelium"])].pivot(index="study", columns="role", values=["share_patient_specific_A", "share_patient_specific_B"])
    out = pd.DataFrame(index=w.index)
    for v in ("A", "B"):
        out[f"patient_share_{v}_malignant"] = w[(f"share_patient_specific_{v}", "malignant")]
        out[f"patient_share_{v}_normal"] = w[(f"share_patient_specific_{v}", "normal_epithelium")]
        out[f"patient_share_{v}_diff"] = out[f"patient_share_{v}_malignant"] - out[f"patient_share_{v}_normal"]
    return out.reset_index()


# ---- Holm --------------------------------------------------------------------------------------------------------

def holm(tests: list) -> pd.DataFrame:
    fam = ["H138.1", "H138.4"]
    p = pd.Series({t["hypothesis"]: t["p"] for t in tests if t["hypothesis"] in fam}).sort_values()
    rows, run_max, run_max6 = [], 0.0, 0.0
    for r, (h, v) in enumerate(p.items()):
        run_max = max(run_max, min(1.0, v * (len(p) - r)))
        run_max6 = max(run_max6, min(1.0, v * (6 - r)))  # six registered members, the four unscored taken as not smaller
        rows.append(dict(hypothesis=h, p=v, rank=r + 1, m_scored=len(p), holm_adj_scored=run_max, reject_scored_0_05=bool(run_max < 0.05),
                         m_registered_family=6, holm_adj_family6=run_max6, reject_family6_0_05=bool(run_max6 < 0.05),
                         provisional="H138.2, H138.3, H138.6, H138.7 (Q2) not yet scored: adjustment provisional; H138.5 not in the family (note 5 ii)"))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    out = OUT / "subsample" if quick else OUT
    out.mkdir(parents=True, exist_ok=True)
    n_perm = 50 if quick else q.N_PERM
    t4, rows4 = score_h138_4(args, quick, n_perm, True)
    log("H138.4:", {k: t4[k] for k in ("statistic", "p", "null_mean", "companion_q95")})
    t5, rows5 = score_h138_5(args, quick, n_perm)
    log("H138.5:", [(t["hypothesis"], round(t["statistic"], 4), t["p"]) for t in t5])
    t1, coh, dec = score_h138_1(args, quick)
    log("H138.1:", {k: t1[k] for k in ("statistic", "pooled_lo95", "p", "pass_rule", "lower_bound_rule_pass")})
    tests = [t1, t4, *t5]
    pd.DataFrame(tests).to_csv(out / "tests.csv", index=False)
    pu = pd.concat([pd.DataFrame(rows4 + rows5), coh.assign(hypothesis="H138.1", unit=coh.study, statistic=coh.c, n=coh.n_patients)], ignore_index=True)
    pu.to_csv(out / "per_unit.csv", index=False)
    holm([t1, t4]).to_csv(out / "holm.csv", index=False)
    d = d138_table(dec)
    d.to_csv(out / "d138_1.csv", index=False)
    contrast_patient_share(d).to_csv(out / "d138_1_contrast.csv", index=False)
    log("done")


if __name__ == "__main__":
    main()
