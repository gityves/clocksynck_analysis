"""A328 part b and A338 scoring (registered in `paper/docs/registered_experiments.md`: A328, A338, A328 note 2 with its
addendum, and the part b power note with the user's decision of 2026-11-22). One run.

Units (one row per fetus): del Valle whole adrenal and Cardoso-Moreira prenatal liver (scored, Holm over the two for
the gate); 10x adrenal cortex (A244 cells; and `cortex_minus1`), heart (Cardoso + Knight-Schrijver) and kidney
(Cardoso + Hochane) (descriptive); Cardoso-Moreira postnatal liver, heart, kidney (reference). Spliced = exon counts /
`spliced` layer, unspliced = intron (gene body minus exon) / `unspliced` layer; direct introns as sensitivity (bulk
units without 10x members). Log fractions residualised gene by gene on cohort and quadratic age (stage indicator in
del Valle), centred within the fit.

Estimators (power note): `matched` (primary): spliced counts thinned gene by gene to the unspliced count level within
each cohort, states averaged over 20 seeded thinning replicates, one spliced projector; `pergene` (companion):
per-gene regression of spliced and unspliced log fractions on the leave-gene-out spliced state, lead = arg sum_g
b_s,g conj(b_u,g); `plain` (registered original, never read alone). Gate: 10,000 permutations of the unspliced rows
within cohort; 90 % interval from 2,000 bootstrap resamples of fetuses within cohort. Companions: 200 expression-
matched and 200 responsive pseudo-clock sets (matched estimator, 5 thinning replicates, 1,000 permutations).

H328.3 (del Valle, liver; 14 genes and core): matched gate (Holm) < 0.05, matched interval within (0, 6 h), pergene
lead > 0; reported with and without the responsive condition (matched delta above the q95 of the responsive sets whose
own gate passes). A338 (liver, 14 genes, matched): phi_f = matched lead; phi_ref = postnatal liver matched lead if its
gate passes, else 1.27 h (A337); r = median over waveform genes of (E/I) prenatal / (E/I) postnatal (genes with mean
intron count >= 10 in both); T_eff = 24 h x r x tan(phi_ref) / tan(phi_f) with the joint bootstrap.

`--check` replaces every unspliced (and direct) count by binomially thinned spliced counts before any statistic (no
real unspliced value is used; expected lead about 0) to test the code.

    python score_a328b.py [--jobs 4] [--check]
Outputs (`../output/score_b[_check]/`): units.csv, compartments.csv, companions.csv, genes.csv, verdicts.csv, a338.csv.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[2] / "A337_timed_lead" / "scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, pseudo_waveforms  # noqa: E402
from common import development as dev  # noqa: E402
from common.mapping import load as mload  # noqa: E402
from power_a328b import CS_WPC, LOW_UNSPLICED, depth_match  # noqa: E402
from score_a337 import boot, interval, kappa, lead, within, wrap  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, N_PERM, N_BOOT, N_PSEUDO, N_PERM_PSEUDO, REPS, REPS_PSEUDO = 3283, 10000, 2000, 200, 1000, 20, 5
HOURS = 24 / (2 * np.pi)
LEAD_A337 = 1.27
WF = cs.gtex()
WAVEFORMS = {"all": WF, "core": WF.drop(["PER1", "PER2", "NR1D1"])}
SCORED = ("delvalle", "liver")
S = {}


# ---- units ---------------------------------------------------------------------------------------------------------

def _bulk(ds, cohort, comp, direct):
    (ce, se), (ci, si) = mload.bulk(ds, "exon"), mload.bulk(ds, "intron")
    m = (se.compartment == comp).to_numpy()
    lib = ce.index[m]
    se = se.loc[lib]
    if "stage" in se and se.stage.notna().all() and ds.endswith("prenatal"):
        age = [CS_WPC[v] if v in CS_WPC else float(v.rstrip("w")) for v in se.stage.astype(str)]
    else:
        age = se.age.to_numpy(float)
    obs = pd.DataFrame(dict(donor=se.donor.astype(str).to_numpy(), cohort=cohort, age=np.asarray(age, float),
                            depth_s=se.depth_exon.to_numpy(float), depth_u=si.loc[lib].depth_intron.to_numpy(float)))
    d = None
    if direct:
        cd, sd = mload.bulk(ds, "intron_direct")
        obs["depth_d"] = sd.loc[lib].depth_intron_direct.to_numpy(float)
        d = cd.loc[lib].reset_index(drop=True)
    return obs, ce.loc[lib].reset_index(drop=True), ci.loc[lib].reset_index(drop=True), d


def _tenx(ds, cohort, cortex):
    out = []
    for layer in ("spliced", "unspliced"):
        c, t = mload.pseudobulk(ds, layer, group="a244_class" if cortex else None)
        if cortex:
            c, t = c[t.compartment == "cortex"], t[t.compartment == "cortex"]
        out.append((c, t))
    (s, t), (u, tu) = out
    assert (t.index == tu.index).all()
    obs = pd.DataFrame(dict(donor=t.donor.to_numpy(), cohort=cohort, age=t.age.to_numpy(float),
                            depth_s=s.sum(axis=1).to_numpy(float), depth_u=u.sum(axis=1).to_numpy(float)))
    return obs, s.reset_index(drop=True), u.reset_index(drop=True), None


def _pool(name, members, primary, role):
    genes = members[0][1].columns
    for m in members[1:]:
        genes = genes.intersection(m[1].columns)
    obs = pd.concat([m[0] for m in members], ignore_index=True)
    obs["stage"] = obs.age.map(lambda a: f"{a:g}")
    assert obs.age.notna().all() and np.isfinite(obs[["depth_s", "depth_u"]].to_numpy()).all(), name
    cat = lambda k: pd.concat([m[k][genes] for m in members], ignore_index=True)  # noqa: E731
    d = cat(3) if all(m[3] is not None for m in members) else None
    return dict(name=name, primary=primary, role=role, obs=obs, s=cat(1), u=cat(2), d=d)


def load_units() -> list[dict]:
    cortex = [_tenx("kameneva2021_fetal_adrenal", "kameneva", True), _tenx("wang2022_fetal_adrenal", "wang", True),
              _tenx("cheng2022_fetal_adrenal", "cheng", True), _tenx("dong2020_fetal_adrenal", "dong", True)]
    units = [
        _pool("delvalle", [_bulk("delvalle2023_fetal_adrenal", "delvalle", "adrenal gland", True)], True, "fetal"),
        _pool("liver", [_bulk("cardosomoreira2019_prenatal", "cardoso", "liver", True)], True, "fetal"),
        _pool("cortex", cortex, False, "fetal"),
        _pool("heart", [_bulk("cardosomoreira2019_prenatal", "cardoso", "heart", False),
                        _tenx("knightschrijver2023_fetal_heart", "knightschrijver", False)], False, "fetal"),
        _pool("kidney", [_bulk("cardosomoreira2019_prenatal", "cardoso", "kidney", False),
                         _tenx("hochane2019_fetal_kidney", "hochane", False)], False, "fetal"),
    ]
    for organ in ("liver", "heart", "kidney"):
        units.append(_pool(f"postnatal_{organ}", [_bulk("cardosomoreira2019_postnatal", "postnatal", organ, True)],
                           False, "reference"))
    c = units[2]
    keep = ~c["obs"].donor.isin(LOW_UNSPLICED).to_numpy()
    assert (~keep).sum() == len(LOW_UNSPLICED)
    units.append(dict(c, name="cortex_minus1", obs=c["obs"][keep].reset_index(drop=True),
                      s=c["s"][keep].reset_index(drop=True), u=c["u"][keep].reset_index(drop=True), d=None))
    return units


def covariates(obs) -> np.ndarray:
    if obs.cohort.iloc[0] == "delvalle":
        return dev.design(obs, age=None, factors=("stage",))
    return dev.design(obs, age="age", df=2, factors=("cohort",))


# ---- estimators ----------------------------------------------------------------------------------------------------

def logfrac(x, depth, Z):
    y = dev.residual(np.log((x + 0.5) / np.broadcast_to(depth, x.shape)), Z)
    return y - y.mean(0)


def plane_states(s, u, ds, du, Z, groups, wf_genes, wf, estimator, rng, reps):
    """(zs, zu) for `plain` and `matched` (matched: zs averaged over `reps` thinning replicates)."""
    L = projector(pd.DataFrame(s, columns=wf_genes), wf)
    to_z = lambda y: (y @ L.T)[:, 0] + 1j * (y @ L.T)[:, 1]  # noqa: E731
    zu = to_z(logfrac(u, du[:, None], Z))
    if estimator == "plain":
        return to_z(logfrac(s, ds[:, None], Z)), zu
    zs = np.zeros(len(s), complex)
    for _ in range(reps):
        sm, dsm = depth_match(s, u, ds, du, groups, rng)
        zs += to_z(logfrac(sm, dsm, Z))
    return zs / reps, zu


def plane_stats(zs, zu, groups, rng, n_perm, n_boot) -> dict:
    d, k = lead(zu, zs), kappa(zu, zs)
    P = within(groups, rng, n_perm)
    kp = np.abs((zu[P] * np.conj(zs)[None]).sum(1)) / np.sqrt(np.sum(np.abs(zu) ** 2) * np.sum(np.abs(zs) ** 2))
    out = dict(kappa=k, gate_p=float((1 + np.sum(kp >= k)) / (1 + n_perm)), delta=d)
    if n_boot:
        B = boot(groups, rng, n_boot)
        draws = np.angle((zu[B] * np.conj(zs[B])).sum(1)) * HOURS
        out.update(zip(("lo", "hi"), interval(d, draws)), draws=wrap(draws - d) + d)
    return out


def pergene_stats(Ys, Yu, s, wf_genes, wf, groups, rng, n_perm, n_boot) -> tuple[dict, pd.DataFrame]:
    P, B = within(groups, rng, n_perm), boot(groups, rng, n_boot)
    c = cp = cb = 0
    ss, su, sup, rows = 0.0, 0.0, np.zeros(n_perm), []
    for j, g in enumerate(wf_genes):
        o = [k for k in range(len(wf_genes)) if k != j]
        Lj = projector(pd.DataFrame(s[:, o], columns=[wf_genes[k] for k in o]), wf.drop([g]))
        X = Ys[:, o] @ Lj.T
        H = np.linalg.pinv(X)
        bs, bu = H @ Ys[:, j], H @ Yu[:, j]
        bs, bu = bs[0] + 1j * bs[1], bu[0] + 1j * bu[1]
        c, ss, su = c + bs * np.conj(bu), ss + abs(bs) ** 2, su + abs(bu) ** 2
        bp = np.einsum("in,pn->pi", H, Yu[P, j])
        bp = bp[:, 0] + 1j * bp[:, 1]
        cp, sup = cp + bs * np.conj(bp), sup + np.abs(bp) ** 2
        Xb = X[B]
        XtX = np.einsum("bni,bnj->bij", Xb, Xb)
        f = lambda Y: np.linalg.solve(XtX, np.einsum("bni,bn->bi", Xb, Y[B, j])[..., None])[..., 0]  # noqa: E731
        bbs, bbu = f(Ys), f(Yu)
        cb = cb + (bbs[:, 0] + 1j * bbs[:, 1]) * np.conj(bbu[:, 0] + 1j * bbu[:, 1])
        rows.append(dict(gene=g, amp_s=abs(bs), amp_u=abs(bu), amp_ratio=abs(bu) / abs(bs),
                         lag_h=float(np.angle(bs * np.conj(bu)) * HOURS)))
    k = abs(c) / np.sqrt(ss * su)
    d = float(np.angle(c) * HOURS)
    draws = np.angle(cb) * HOURS
    lo, hi = interval(d, draws)
    return (dict(kappa=float(k), gate_p=float((1 + np.sum(np.abs(cp) / np.sqrt(ss * sup) >= k)) / (1 + n_perm)),
                 delta=d, lo=lo, hi=hi, draws=wrap(draws - d) + d), pd.DataFrame(rows))


# ---- tasks ---------------------------------------------------------------------------------------------------------

def counts(unit, genes, which):
    s = unit["s"][genes].to_numpy(np.int64)
    u = unit[which][genes].to_numpy(np.int64)
    du = unit["obs"]["depth_u" if which == "u" else "depth_d"].to_numpy(float)
    return s, u, unit["obs"].depth_s.to_numpy(float), du


def task_main(k):
    unit, wname, which = S["tasks"][k]
    wf, obs = WAVEFORMS[wname], unit["obs"]
    genes = list(wf.genes)
    rng = np.random.default_rng([SEED, k])
    s, u, ds, du = counts(unit, genes, which)
    Z, groups = covariates(obs), obs.cohort.to_numpy()
    base = dict(unit=unit["name"], role=unit["role"], primary=unit["primary"], waveform=wname,
                unspliced="intron" if which == "u" else "direct", n=len(obs))
    rows, draws = [], {}
    for est in ("plain", "matched"):
        zs, zu = plane_states(s, u, ds, du, Z, groups, genes, wf, est, rng, REPS)
        r = plane_stats(zs, zu, groups, rng, N_PERM, N_BOOT)
        draws[est] = r.pop("draws")
        rows.append(dict(base, estimator=est, **r))
    r, g = pergene_stats(logfrac(s, ds[:, None], Z), logfrac(u, du[:, None], Z), s, genes, wf, groups, rng, N_PERM,
                         N_BOOT)
    draws["pergene"] = r.pop("draws")
    rows.append(dict(base, estimator="pergene", **r))
    g = g.assign(unit=unit["name"], waveform=wname, unspliced=base["unspliced"],
                 mean_spliced=s.mean(0), mean_unspliced=u.mean(0))
    return rows, g, (unit["name"], wname, base["unspliced"], draws)


def task_pseudo(k):
    unit, wname = S["ptasks"][k]
    wf, obs = WAVEFORMS[wname], unit["obs"]
    sp = unit["s"]
    lm = dev.log_mean(sp, obs.depth_s.to_numpy(float), min_count=0.0)
    lm = lm[(sp.mean().reindex(lm.index) >= 1).to_numpy() | lm.index.isin(wf.genes)]
    Z, groups = covariates(obs), obs.cohort.to_numpy()
    rows = []
    for j, kind in enumerate(("matched", "responsive"), 1):
        sets = pseudo_waveforms(lm, N_PSEUDO, np.random.default_rng([SEED, k, j]), wf, kind)
        for i, w in enumerate(sets):
            genes = list(w.genes)
            s, u, ds, du = counts(unit, genes, "u")
            rng = np.random.default_rng([SEED, k, j, i])
            zs, zu = plane_states(s, u, ds, du, Z, groups, genes, w, "matched", rng, REPS_PSEUDO)
            r = plane_stats(zs, zu, groups, rng, N_PERM_PSEUDO, 0)
            rows.append(dict(unit=unit["name"], waveform=wname, kind=kind, set=i, **r))
    return pd.DataFrame(rows)


# ---- readings ------------------------------------------------------------------------------------------------------

def verdicts(t: pd.DataFrame, comp: pd.DataFrame) -> pd.DataFrame:
    q = (comp[comp.gate_p < .05].groupby(["unit", "waveform", "kind"]).delta.quantile(.95)
         .unstack("kind").rename(columns=lambda c: f"{c}_delta_q95"))
    share = comp.groupby(["unit", "waveform", "kind"]).gate_p.apply(lambda p: (p < .05).mean()).unstack("kind")
    m = t[(t.unspliced == "intron")].pivot_table(index=["unit", "waveform"], columns="estimator",
                                                values=["gate_p", "delta", "lo", "hi", "kappa"])
    m.columns = [f"{a}_{b}" for a, b in m.columns]
    v = m.join(q).join(share.rename(columns=lambda c: f"{c}_gate_share")).reset_index()
    v["scored"] = v.unit.isin(SCORED)
    v["gate_holm"] = np.nan
    for w in WAVEFORMS:
        sel = v.scored & (v.waveform == w)
        v.loc[sel, "gate_holm"] = dev.holm(v.loc[sel, "gate_p_matched"].to_numpy())
    v["H328_3"] = v.scored & (v.gate_holm < .05) & (v.lo_matched > 0) & (v.hi_matched < 6) & (v.delta_pergene > 0)
    v["H328_3_responsive"] = v.H328_3 & (v.delta_matched > v.responsive_delta_q95)
    for h in ("H328_3", "H328_3_responsive"):
        v[h] = v[h].astype(object).where(v.scored, None)
    return v


def a338(units, draws, t) -> pd.DataFrame:
    """Liver, 14 genes, matched: T_eff with the postnatal (if its gate passes) or A337 phi_ref and the lifetime ratio."""
    rng = np.random.default_rng([SEED, 338])
    get = lambda name: next(u for u in units if u["name"] == name)  # noqa: E731
    pre, post = get("liver"), get("postnatal_liver")
    row = t[(t.unit == "liver") & (t.waveform == "all") & (t.unspliced == "intron") & (t.estimator == "matched")].iloc[0]
    ref = t[(t.unit == "postnatal_liver") & (t.waveform == "all") & (t.unspliced == "intron")
            & (t.estimator == "matched")].iloc[0]
    phi_f = draws[("liver", "all", "intron")]["matched"] / HOURS
    if ref.gate_p < .05:
        phi_ref, ref_src = draws[("postnatal_liver", "all", "intron")]["matched"] / HOURS, "postnatal liver"
    else:
        phi_ref, ref_src = np.full(len(phi_f), LEAD_A337 / HOURS), "A337 (postnatal gate not passed)"
    genes = list(WF.genes)
    E1, I1 = pre["s"][genes].to_numpy(float), pre["u"][genes].to_numpy(float)
    E0, I0 = post["s"][genes].to_numpy(float), post["u"][genes].to_numpy(float)
    ok = (I1.mean(0) >= 10) & (I0.mean(0) >= 10)

    def ratio(i1, i0):
        return float(np.median((E1[i1][:, ok].sum(0) / I1[i1][:, ok].sum(0))
                               / (E0[i0][:, ok].sum(0) / I0[i0][:, ok].sum(0))))

    r = ratio(np.arange(len(E1)), np.arange(len(E0)))
    rb = np.array([ratio(rng.choice(len(E1), len(E1)), rng.choice(len(E0), len(E0))) for _ in range(len(phi_f))])

    def teff(pf, pr, rr):
        pf = np.asarray(pf, float)
        return np.where(pf <= 0, np.inf, np.where(pf >= np.pi / 2, 0.0, 24 * rr * np.tan(pr) / np.tan(
            np.clip(pf, 1e-9, np.pi / 2 - 1e-9))))

    T = teff(phi_f, rng.permutation(phi_ref), rb)
    T_hat = float(teff(row.delta / HOURS, np.median(phi_ref), r))
    Tf = np.where(np.isinf(T), 1e9, T)  # infinite T_eff (phi_f <= 0) kept as a very large value for the quantiles
    lo, hi = float(np.quantile(Tf, .05)), float(np.quantile(Tf, .95))
    h1 = (row.gate_p < .05) & (row.lo > 0) & (lo > 6)
    if h1 and lo <= 24 <= hi and hi < 96:
        reading = "daily timescale"
    elif h1 and lo > 96:
        reading = "slow input (days)"
    elif hi < 6:
        reading = "fast, pulsatile input"
    elif row.lo <= 0 <= row.hi:
        reading = "no daily-scale input of the reference size"
    else:
        reading = "timescale not resolved at this design"
    return pd.DataFrame([dict(unit="liver", waveform="all", estimator="matched", phi_f_h=row.delta,
                              phi_f_lo=row.lo, phi_f_hi=row.hi, gate_p=row.gate_p, phi_ref_source=ref_src,
                              phi_ref_h=float(np.median(phi_ref) * HOURS), r=r, r_lo=float(np.quantile(rb, .05)),
                              r_hi=float(np.quantile(rb, .95)), genes_r=int(ok.sum()), T_eff=T_hat, T_lo=lo, T_hi=hi,
                              H338_1=bool(h1), reading=reading)])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    global N_PERM, N_BOOT, N_PSEUDO, N_PERM_PSEUDO
    limit_threads(1)
    t0 = time.time()
    out = OUT / ("score_b_check" if a.check else "score_b")
    out.mkdir(parents=True, exist_ok=True)
    units = load_units()
    if a.check:  # no real unspliced count: unspliced := thinned spliced
        N_PERM, N_BOOT, N_PSEUDO, N_PERM_PSEUDO = 500, 300, 10, 200
        rng = np.random.default_rng(0)
        for u in units:
            th = lambda x: pd.DataFrame(rng.binomial(x.to_numpy(np.int64), 0.15), columns=x.columns)  # noqa: E731
            u["u"] = th(u["s"])
            u["obs"]["depth_u"] = u["u"].sum(axis=1).to_numpy(float)
            if u["d"] is not None:
                u["d"] = th(u["s"])
                u["obs"]["depth_d"] = u["d"].sum(axis=1).to_numpy(float)
    info = pd.DataFrame([dict(unit=u["name"], role=u["role"], n=len(u["obs"]), cohorts=u["obs"].cohort.nunique(),
                              genes=u["s"].shape[1], direct=u["d"] is not None,
                              clock_present=sum(g in u["s"] for g in WF.genes)) for u in units])
    info.to_csv(out / "units.csv", index=False)
    pd.set_option("display.width", 250)
    print(info.to_string(index=False), f"load {time.time() - t0:.0f}s", sep="\n", flush=True)
    S["tasks"] = [(u, w, which) for u in units for w in WAVEFORMS
                  for which in (("u", "d") if u["d"] is not None else ("u",))]
    S["ptasks"] = [(u, w) for u in units if u["role"] == "fetal" for w in WAVEFORMS]
    rows, genes, draws, comps = [], [], {}, []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, (r, g, (key0, key1, key2, d)) in enumerate(pool.imap_unordered(task_main, range(len(S["tasks"])))):
            rows += r
            genes.append(g)
            draws[(key0, key1, key2)] = d
            print(f"main {k + 1}/{len(S['tasks'])} {time.time() - t0:.0f}s", flush=True)
        for k, c in enumerate(pool.imap_unordered(task_pseudo, range(len(S["ptasks"])))):
            comps.append(c)
            print(f"companions {k + 1}/{len(S['ptasks'])} {time.time() - t0:.0f}s", flush=True)
    t = pd.DataFrame(rows).sort_values(["role", "unit", "waveform", "unspliced", "estimator"])
    t.to_csv(out / "compartments.csv", index=False)
    comp = pd.concat(comps)
    comp.to_csv(out / "companions.csv", index=False)
    pd.concat(genes).to_csv(out / "genes.csv", index=False)
    v = verdicts(t, comp)
    v.to_csv(out / "verdicts.csv", index=False)
    x = a338(units, {k[:3]: d for k, d in draws.items()}, t)
    x.to_csv(out / "a338.csv", index=False)
    print(t.drop(columns=["role", "primary"]).round(3).to_string(index=False),
          v.round(3).to_string(index=False), x.round(3).T.to_string(), f"runtime {time.time() - t0:.0f}s",
          sep="\n", flush=True)


if __name__ == "__main__":
    main()
