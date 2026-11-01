"""A188 (C5): do malignant cells follow their own tumour stroma less in high-proliferation (and ER-) tumours?

Registered in paper/docs/registered_experiments.md (A188). Cohorts and patients: those of H138.1 (A17 learned stromal
anchors, malignant = called malignant cells else tumour epithelium; `A138.../q3_lib.q1_design`). Proliferation = A97 G2M
score of the patient's malignant pseudobulk (mean z of log2 CPM + 1 of the Hallmark G2M genes without clock genes, z over
the patients of the cohort; counts from `prep_g2m_counts.py`), tertiles by rank within cohort, high vs low.

Stage 1 (per cohort, cached in output/cohorts): `Dataset.contrast(levels=(low, high), targets=[malignant] + companion,
require_host_gate_each_arm=True, anchor_readout="loo", n_perm=1000, blocks="none")` and the power curve of the malignant
compartment (`clocksync.power.mde`, the call `contrast(mde=True)` makes). Stage 2: pooling (inverse-variance on
se_d_star, same weights over the per-cohort permutation draws), the branch rule on the pooled MDE = 2.49 x pooled se
(<= 0.3 test, else descriptive only), leave-one-cohort-out. H188.2 (ER) is checked for eligible cohorts.

    python score_a188.py [--jobs 4] [--only STUDY ...] [--stage2-only]
Outputs (../output): prep_patients.csv, prep.csv, gates.csv, per_cohort.csv, power.csv, tests.csv or descriptive.csv, loco.csv
"""
from __future__ import annotations

import argparse, pickle, sys, threading, time, zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "output"
G2M = OUT / "g2m"
sys.path.insert(0, str(HERE.parents[1] / "A138_cancer_estimability" / "scripts"))
import numpy as np, pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import q3_lib as q  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync import power  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

N_PERM, N_BOOT, N_MDE, MIN_ARM, MDE_LIMIT, Z = 1000, 200, 100, 5, 0.3, 2.49
T0 = time.time()
TUMOUR_ANCHORS = "--tumour-anchors" in sys.argv


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def seed_of(name: str) -> int:
    return zlib.crc32(name.encode()) % (2 ** 31)


def designs() -> list:
    gate = q.a17_gate()
    ds = [d for s in sorted(s for s, g in gate.items() if g) if (d := q.q1_design(s))]
    ref = pd.read_csv(HERE.parents[1] / "A138_cancer_estimability/output/power/h138_1_designs.csv")
    assert [d["study"] for d in ds] == list(ref.study), "cohorts differ from H138.1"
    if TUMOUR_ANCHORS:
        for d in ds:
            d["anchors"] = [a for a in d["anchors"] if a.startswith("tumour:")]
        ds = [d for d in ds if len(d["anchors"]) >= 2]
    return ds


# ---- prep -----------------------------------------------------------------------------------------------------------

def g2m_scores(study: str, pts) -> pd.Series:
    """A97 score: mean z (ddof 0, over the cohort's patients) of log2(CPM + 1) of the G2M genes of the malignant pseudobulk."""
    t = pd.read_csv(G2M / f"g2m_{study}.csv", dtype={"patient": str}).set_index("patient")
    t = t.loc[sorted(set(t.index) & set(pts))]
    g = [c for c in t.columns if c not in ("n_cells", "depth")]
    lc = np.log2(t[g].div(t.depth, axis=0) * 1e6 + 1)
    z = (lc - lc.mean()) / lc.std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1).rename("g2m")


def tertiles(x: pd.Series) -> pd.Series:
    return pd.Series(pd.qcut(x.rank(method="first"), 3, labels=["low", "mid", "high"]).astype(str), index=x.index)


def companions(pb, d, mal, lab) -> list:
    """Registered rule: myeloid, else the non-anchor non-malignant tumour compartments (adjacent-normal compartments and the
    non-malignant epithelium of the tumour sample are left out: CHOICE). Each needs >= MIN_ARM patients in each arm."""
    arms = {g: set(lab.index[lab == g]) for g in ("low", "high")}
    ok = lambda c: all(len(set(pb[pb.compartment == c].donor) & a) >= MIN_ARM for a in arms.values())
    comps = set(pb.compartment)
    if "tumour:myeloid" in comps and ok("tumour:myeloid"):
        return ["tumour:myeloid"]
    skip = set(d["anchors"]) | {mal, q.EPI, "tumour:epithelial_nonmalignant", q.NORM}
    return sorted(c for c in comps if c.startswith("tumour:") and c not in skip and ok(c))


def prep(ds: list):
    pats, rows = [], []
    for d in ds:
        pb, mal = q.study_table(d["study"])
        pts = sorted(set(pb[pb.compartment == mal].donor))
        sc = g2m_scores(d["study"], pts)
        t = tertiles(sc)
        lab = t[t != "mid"]
        pats.append(pd.DataFrame(dict(study=d["study"], patient=sc.index, g2m=sc.to_numpy(), tertile=t.to_numpy())))
        comp = companions(pb, d, mal, lab)
        nst = pb[pb.compartment == mal].stratum.astype(str).nunique() if "stratum" in pb else 1
        rows.append(dict(cohort=d["study"], malignant=mal, n_patients=len(pts), n_g2m=len(sc), n_low=int((t == "low").sum()),
                         n_high=int((t == "high").sum()), anchors=";".join(d["anchors"]), n_anchors=len(d["anchors"]),
                         companion=";".join(comp), companion_rule="myeloid" if comp == ["tumour:myeloid"] else "mean of non-anchor tumour compartments" if comp else "none",
                         n_strata=nst, eligible_ge5_per_arm=bool((t == "low").sum() >= MIN_ARM and (t == "high").sum() >= MIN_ARM),
                         er_labels=False))
    pats = pd.concat(pats, ignore_index=True)
    pats.to_csv(OUT / "prep_patients.csv", index=False)
    p = pd.DataFrame(rows)
    p.to_csv(OUT / "prep.csv", index=False)
    return pats, p


# ---- stage 1 --------------------------------------------------------------------------------------------------------

def one_cohort(d: dict, lab: pd.Series, comp: list):
    limit_threads(1)
    f = OUT / "cohorts" / f"{d['study']}.pkl"
    if f.exists():
        return pickle.loads(f.read_bytes())
    t0, sd = time.time(), seed_of(d["study"])
    pb, mal = q.study_table(d["study"])
    t = pb[pb.compartment.isin(d["anchors"] + [mal] + comp) & pb.donor.isin(lab.index)]
    ds = cs.Dataset(t, anchors=d["anchors"])
    targets = [mal] + comp
    c = ds.contrast(lab, levels=("low", "high"), targets=targets, companion=None, blocks="none", n_perm=N_PERM, n_boot=N_BOOT,
                    seed=sd, require_host_gate_each_arm=True, n_gate=N_PERM, anchor_readout="loo")
    res = dict(study=d["study"], mal=mal, comp=comp, seed=sd, table=c.table, gate=c.gate, null=c.null, observed=c.observed, mde=None)
    if c.null is not None:  # both gates passed
        m = power.mde(ds.study, lab, ds.anchors, mal, n_rep=N_MDE, n_perm=200, seed=sd, case="high", control="low", within_stratum=False)
        res["mde"] = {k: v for k, v in m.items() if k != "curve"} | dict(curve=m["curve"])
    res["seconds"] = round(time.time() - t0)
    f.write_bytes(pickle.dumps(res))
    return res


def stage1(ds, pats, prepd, only, jobs):
    (OUT / "cohorts").mkdir(exist_ok=True)
    todo = []
    for d in ds:
        r = prepd[prepd.cohort == d["study"]].iloc[0]
        if not r.eligible_ge5_per_arm or (only and d["study"] not in only):
            continue
        p = pats[pats.study == d["study"]]
        lab = p[p.tertile != "mid"].set_index("patient").tertile
        todo.append((d, lab, r.companion.split(";") if isinstance(r.companion, str) and r.companion else []))
    say(f"stage 1: {len(todo)} cohorts, jobs {jobs}")
    for r in Parallel(n_jobs=jobs, return_as="generator_unordered")(delayed(one_cohort)(*a) for a in todo):
        say(f"  {r['study']}: gates {[(g, v.get('passes')) for g, v in r['gate'].items()]} ({r['seconds']} s)")


# ---- stage 2 --------------------------------------------------------------------------------------------------------

def load_results(ds):
    out = {}
    for d in ds:
        f = OUT / "cohorts" / f"{d['study']}.pkl"
        if f.exists():
            out[d["study"]] = pickle.loads(f.read_bytes())
    return out


def pooled(rows, nulls, w):
    """Inverse-variance mean of per-cohort values and, per permutation draw, of the per-cohort draws (same weights)."""
    w = np.asarray(w, float)
    obs = float(np.sum(w * np.asarray(rows)) / w.sum())
    nul = np.nansum(nulls * w[:, None], axis=0) / np.sum(np.where(np.isfinite(nulls), w[:, None], 0.0), axis=0)
    return obs, nul


def comp_series(r, scale):
    """Observed and permutation draws of the companion (mean over the companion compartments) and malignant d*."""
    nl, ob = r["null"][scale], r["observed"][scale]
    if not r["comp"]:
        return None
    return float(np.nanmean([ob[k] for k in r["comp"]])), nl[r["comp"]].mean(axis=1, skipna=True).to_numpy(float)


def summarize(res, label, studies, strict_branch=True):
    """Pooled d*, null-based p, specificity and MDE over `studies` (cohorts that passed both gates). Returns dict."""
    rows, nulls, cobs, cnull, w, se = [], [], [], [], [], []
    for s in studies:
        r = res[s]
        tb = r["table"].loc[r["mal"]]
        rows.append(tb.d_star); se.append(tb.se_d_star); w.append(1 / tb.se_d_star ** 2)
        nulls.append(r["null"]["d_star"][r["mal"]].to_numpy(float))
        cs_ = comp_series(r, "d_star")
        cobs.append(cs_[0] if cs_ else np.nan); cnull.append(cs_[1] if cs_ else np.full(N_PERM, np.nan))
    w = np.array(w); nulls = np.array(nulls)
    obs, nul = pooled(rows, nulls, w)
    pse = float(np.sum(w) ** -0.5)
    out = dict(analysis=label, n_cohorts=len(studies), cohorts=";".join(studies), pooled_d_star=obs, pooled_se=pse, mde=Z * pse,
               ci_lo=obs - 1.96 * pse, ci_hi=obs + 1.96 * pse, null_sd=float(np.nanstd(nul, ddof=1)), n_perm=int(np.isfinite(nul).sum()))
    out["p_one_sided_less"] = S.perm_p(nul, obs, "less")
    ok = [i for i, x in enumerate(cobs) if np.isfinite(x)]
    if ok and len(ok) == len(studies):
        cobs_p, cnul = pooled([cobs[i] for i in ok], np.array([cnull[i] for i in ok]), w[ok])
        _, nmal = pooled([rows[i] for i in ok], nulls[ok], w[ok])
        o_mal, _ = pooled([rows[i] for i in ok], nulls[ok], w[ok])
        diff, dnull = o_mal - cobs_p, nmal - cnul
        sd = float(np.nanstd(dnull, ddof=1))
        out.update(pooled_d_star_companion=cobs_p, spec_diff=diff, spec_p_one_sided_less=S.perm_p(dnull, diff, "less"),
                   spec_null_sd=sd, spec_ci_lo=diff - 1.96 * sd, spec_ci_hi=diff + 1.96 * sd)
    return out


def stage2(ds, prepd):
    res = load_results(ds)
    gates, per, pw = [], [], []
    for s, r in res.items():
        for arm, g in r["gate"].items():
            gates.append(dict(cohort=s, arm=arm, **{k: g.get(k) for k in ("p_pairing", "p_acrophase", "clock_share", "passes", "error")}))
        t = r["table"].reset_index()
        t.insert(0, "cohort", s)
        t["role"] = np.where(t.compartment == r["mal"], "malignant", "companion")
        both = all(g.get("passes", False) for g in r["gate"].values())
        t["both_gates_pass"] = both
        t["n_low_pats"], t["n_high_pats"] = t.get("donors_low"), t.get("donors_high")
        per.append(t)
        m = r["mde"]
        mdd = pd.Series(m).drop("curve") if m else pd.Series(dtype=float)
        tb = r["table"].loc[r["mal"]] if both else None
        pw.append(dict(cohort=s, both_gates_pass=both, n_low=None if tb is None else tb.get("donors_low"), n_high=None if tb is None else tb.get("donors_high"),
                       se_d_star=np.nan if tb is None else tb.se_d_star, mde_d_star_se=np.nan if tb is None else Z * tb.se_d_star,
                       **{f"power_{k}": v for k, v in mdd.items()}))
    gates = pd.DataFrame(gates); per = pd.concat(per, ignore_index=True); pw = pd.DataFrame(pw)
    gates.to_csv(OUT / "gates.csv", index=False)
    per.to_csv(OUT / "per_cohort.csv", index=False)
    mal = per[per.role == "malignant"]
    enter = [s for s in res if res[s]["null"] is not None and np.isfinite(res[s]["table"].loc[res[s]["mal"]].se_d_star)
             and res[s]["table"].loc[res[s]["mal"]].get("donors_low", 0) >= MIN_ARM and res[s]["table"].loc[res[s]["mal"]].get("donors_high", 0) >= MIN_ARM]
    pw["enters_pool"] = pw.cohort.isin(enter)
    # pooled power first
    se = float(np.sum([1 / pw.set_index("cohort").loc[s].se_d_star ** 2 for s in enter]) ** -0.5) if enter else np.nan
    pooled_row = dict(cohort="POOLED (inverse variance)", both_gates_pass=True, enters_pool=True, se_d_star=se, mde_d_star_se=Z * se)
    pw = pd.concat([pw, pd.DataFrame([pooled_row])], ignore_index=True)
    pw.to_csv(OUT / "power.csv", index=False)
    say("power written: pooled se", se, "MDE", Z * se, "cohorts", enter)
    if not enter:
        return
    h1 = summarize(res, "H188.1 proliferation (high - low)", enter)
    test = h1["mde"] <= MDE_LIMIT
    h1["branch"] = "test" if test else "descriptive only: not estimable at this n"
    if not test:
        for k in ("p_one_sided_less", "spec_p_one_sided_less"):
            h1.pop(k, None)
    loco = []
    for s in enter:
        rest = [x for x in enter if x != s]
        if rest:
            l = summarize(res, f"H188.1 without {s}", rest)
            if not test:
                for k in ("p_one_sided_less", "spec_p_one_sided_less"):
                    l.pop(k, None)
            loco.append(l | dict(left_out=s))
    if not loco:
        loco = [dict(analysis="leave-one-cohort-out", note=f"not defined: {len(enter)} cohort(s) entered the pool")]
    pd.DataFrame(loco).to_csv(OUT / "loco.csv", index=False)
    h2 = dict(analysis="H188.2 ER- vs ER+ (breast cohorts with ER labels)", n_cohorts=0, branch="not estimable: no eligible cohort",
              note="Wu (GSE176078) and Pal (GSE161529) are not H138.1 cohorts: no A17 learned anchors and no stromal pseudobulk table; "
                   "the one breast H138.1 cohort (SRP273482, 7 patients) has no ER labels and < 5 patients per arm")
    pd.DataFrame([h1, h2]).to_csv(OUT / ("tests.csv" if test else "descriptive.csv"), index=False)
    say("H188.1:", h1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--stage2-only", action="store_true")
    ap.add_argument("--tumour-anchors", action="store_true",
                    help="POST HOC sensitivity: host = the learned anchors that are tumour-sample compartments only (own tumour "
                         "stroma); outputs in output/post_hoc_tumour_anchors")
    a = ap.parse_args()
    if a.tumour_anchors:
        OUT = OUT / "post_hoc_tumour_anchors"
        OUT.mkdir(exist_ok=True)
    threading.Thread(target=lambda: [(time.sleep(600), say("alive")) for _ in iter(int, 1)], daemon=True).start()
    ds = designs()
    pats, prepd = prep(ds)
    say("prep written\n" + prepd.drop(columns=["anchors"]).to_string())
    if not a.stage2_only:
        stage1(ds, pats, prepd, a.only, a.jobs)
    stage2(ds, prepd)
    say("done")
