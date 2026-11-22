"""A328 part b and A338 power stage (registered in `paper/docs/registered_experiments.md`: A328, A338, A328 note 2 and its
addendum), before any fetal or postnatal unspliced statistic.

Read from the fetal data (`prep`, in memory): per compartment the design (cohort, age or stage, genome-wide spliced and
unspliced depths per fetus) and, per cohort, the median count of each clock gene in each count type. These are the
only quantities read; no per-fetus clock-gene value is written or printed.

Stand-ins (note 2): per compartment, n distinct GTEx samples of the matching tissue (deaths 08-18 h, any Hardy class).
Spliced counts = the GTEx clock-gene counts (clock at relative amplitude a, `simulate.scale_clock`), thinned binomially
to the compartment's spliced depth (target = cohort median count x fetal depth / cohort median depth, lowered by one
common factor per gene where the adult libraries hold fewer reads, as part a). Unspliced counts = the same samples with
the clock-plane state turned by the lead and multiplied by an amplitude factor (`simulate.clock_change`,
`simulate.thin_clock`), thinned to the unspliced depth. Covariates, permutations and bootstrap within cohort as
registered (cohort indicator and quadratic age; stage indicator in del Valle).

Settings: `a1` (lead 1.27 h, A337), `a0.75`, `a0.5`; `a1_lead0.61`, `a1_lead3.00` (A337 range); `steady` (no turn);
`shuffled` (a1, unspliced rows permuted within cohort); A338: `a338_24h` (turn 1.27 h, amplitude 1 / cos phi),
`a338_3h` (phi = arctan(2 pi tau / 3 h), tau from the 1.27 h lead at 24 h; amplitude 1 / cos phi); `greenwell` (the
A337 Greenwell AL + NR livers, real spliced and unspliced counts, n = min(n, 36), thinned to the compartment: the
model-free check).
Readings: A328 pass = gate p < 0.05 and the 90 % interval of delta above 0; decidable if pass rate >= 0.8 at a1.
A338 "daily" = gate passes, phi interval above 0, T_eff interval above 6 h, containing 24 h and below 96 h, with
T_eff = 24 h x tan(phi_ref) / tan(phi), phi_ref = 1.27 h (A337 fallback; r = 1 for stand-ins); decidable if the daily
rate is >= 0.8 under a338_24h and <= 0.1 under steady and a338_3h.

    python power_a328b.py [--n-sub 200] [--jobs 4] [--check]
Outputs (`../output/power_b/`): design.csv, detection.csv, pools.csv, subsamples.csv, passrates.csv, decisions.csv
(`_check` with --check). Controlled access: GTEx time of death in memory only; outputs are aggregates.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[2] / "A337_timed_lead" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from clocksync.simulate import clock_change, scale_clock, thin_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402
from common.mapping import load as mload  # noqa: E402
from score_a337 import boot, gate_p, interval, kappa, lead, within, wrap  # noqa: E402
import score_a337  # noqa: E402

OUT = HERE.parents[1] / "output" / "power_b"
SEED, N_PERM, N_BOOT = 3282, 500, 500
WF = cs.gtex()
GENES = list(WF.genes)
WAVEFORMS = {"all": WF, "core": WF.drop(["PER1", "PER2", "NR1D1"])}
HOURS = 24 / (2 * np.pi)
LEAD = 1.27                                   # h, A337 median timed lead
PHI_REF = LEAD / HOURS
TAU = HOURS * np.tan(PHI_REF)                 # h, mRNA lifetime implied by the lead at 24 h
PHI_3H = np.arctan(2 * np.pi * TAU / 3.0)
SETTINGS = {  # name: (clock scale, unspliced turn in h, unspliced amplitude factor)
    "a1": (1.0, LEAD, 1.0), "a0.75": (0.75, LEAD, 1.0), "a0.5": (0.5, LEAD, 1.0),
    "a1_lead0.61": (1.0, 0.61, 1.0), "a1_lead3.00": (1.0, 3.00, 1.0), "steady": (1.0, 0.0, 1.0),
    "a338_24h": (1.0, LEAD, 1 / np.cos(PHI_REF)), "a338_3h": (1.0, PHI_3H * HOURS, 1 / np.cos(PHI_3H))}
TISSUE = {"adrenal": "Adrenal Gland", "liver": "Liver", "heart": "Heart - Left Ventricle", "kidney": "Kidney - Cortex"}
CS_WPC = {"CS13": 4, "CS14": 4, "CS16": 5, "CS17": 6, "CS18": 6, "CS19": 7, "CS20": 7, "CS21": 7, "CS22": 8,
          "CS23": 8}  # Carnegie stage -> wpc as in the authors' library table
ESTIMATORS = ("plain", "matched", "pergene")
LOW_UNSPLICED = ("kameneva:10X_20_016",)  # note 2 addendum: labelled sensitivity (cortex-only unspliced shortfall)
S = {}


# ---- fetal designs (blind-safe) -----------------------------------------------------------------------------------

def _tenx(ds, cohort, cortex):
    out = []
    for layer in ("spliced", "unspliced"):
        c, t = mload.pseudobulk(ds, layer, group="a244_class" if cortex else None)
        if cortex:
            c, t = c[t.compartment == "cortex"], t[t.compartment == "cortex"]
        out.append((c.reindex(columns=GENES, fill_value=0), c.sum(axis=1), t))
    (s, ds_, t), (u, du, tu) = out
    assert (t.index == tu.index).all()
    obs = pd.DataFrame(dict(donor=t.donor.to_numpy(), cohort=cohort, age=t.age.to_numpy(float),
                            depth_s=ds_.to_numpy(float), depth_u=du.to_numpy(float)))
    return obs, s.to_numpy(float), u.to_numpy(float)


def _bulk(ds, cohort, comp):
    (ce, se), (ci, si) = mload.bulk(ds, "exon"), mload.bulk(ds, "intron")
    m = (se.compartment == comp).to_numpy()
    se, ce, ci = se[m], ce[m], ci.loc[ce.index[m]]
    if "stage" in se and se.stage.notna().all():
        age = [CS_WPC[v] if v in CS_WPC else float(v.rstrip("w")) for v in se.stage.astype(str)]
    else:
        age = se.age.to_numpy(float)
    obs = pd.DataFrame(dict(donor=se.donor.astype(str).to_numpy(), cohort=cohort, age=np.asarray(age, float),
                            depth_s=se.depth_exon.to_numpy(float), depth_u=si.loc[ce.index].depth_intron.to_numpy(float)))
    return obs, ce.reindex(columns=GENES, fill_value=0).to_numpy(float), ci.reindex(columns=GENES, fill_value=0).to_numpy(float)


def prep():
    """(designs, detection): one obs table per compartment and the cohort medians of each clock gene."""
    parts = {
        "cortex": ("adrenal", [_tenx("kameneva2021_fetal_adrenal", "kameneva", True),
                               _tenx("wang2022_fetal_adrenal", "wang", True),
                               _tenx("cheng2022_fetal_adrenal", "cheng", True),
                               _tenx("dong2020_fetal_adrenal", "dong", True)]),
        "delvalle": ("adrenal", [_bulk("delvalle2023_fetal_adrenal", "delvalle", "adrenal gland")]),
        "liver": ("liver", [_bulk("cardosomoreira2019_prenatal", "cardoso", "liver")]),
        "heart": ("heart", [_bulk("cardosomoreira2019_prenatal", "cardoso", "heart"),
                            _tenx("knightschrijver2023_fetal_heart", "knightschrijver", False)]),
        "kidney": ("kidney", [_bulk("cardosomoreira2019_prenatal", "cardoso", "kidney"),
                              _tenx("hochane2019_fetal_kidney", "hochane", False)]),
    }
    designs, det = [], []
    for comp, (tissue, members) in parts.items():
        for obs, s, u in members:
            for kind, x in (("spliced", s), ("unspliced", u)):
                det.append(dict(compartment=comp, cohort=obs.cohort.iloc[0], kind=kind,
                                **dict(zip(GENES, np.median(x, axis=0)))))
            designs.append(obs.assign(compartment=comp, tissue=tissue))
    d = pd.concat(designs, ignore_index=True)
    d["stage"] = d.age.map(lambda a: f"{a:g}")
    low = d[(d.compartment == "cortex") & d.donor.isin(LOW_UNSPLICED)]
    assert len(low) == len(LOW_UNSPLICED)
    d = pd.concat([d, d[(d.compartment == "cortex") & ~d.donor.isin(LOW_UNSPLICED)].assign(compartment="cortex_minus1")])
    det = pd.DataFrame(det)
    det = pd.concat([det, det[det.compartment == "cortex"].assign(compartment="cortex_minus1")])
    return d.reset_index(drop=True), det.reset_index(drop=True)


def covariates(obs) -> np.ndarray:
    if obs.cohort.iloc[0] == "delvalle":
        return dev.design(obs, age=None, factors=("stage",))
    return dev.design(obs, age="age", df=2, factors=("cohort",))


def targets(obs, det, kind) -> np.ndarray:
    """Cohort median count of each gene x fetal depth / cohort median depth (fetus x gene)."""
    col = "depth_s" if kind == "spliced" else "depth_u"
    m = det[det.kind == kind].set_index("cohort")[GENES]
    med = obs.groupby("cohort")[col].transform("median").to_numpy(float)
    return m.loc[obs.cohort].to_numpy(float) * (obs[col].to_numpy(float) / med)[:, None]


# ---- stand-ins -----------------------------------------------------------------------------------------------------

def pools(rng):
    out, rows = {}, []
    for key, tissue in TISSUE.items():
        pb, d = gtex.load([tissue], hardy=(0, 1, 2, 3, 4), genes="clock", depth="library_size")
        pb = pb[pb.donor.isin(set(d.index[in_window(d.tod.to_numpy())]))].reset_index(drop=True)
        clock = pb[["donor", "compartment", "depth", *GENES]].copy()
        depth = clock.depth.to_numpy(float)
        xs, xu = {}, {}
        for name, (a, turn, amp) in SETTINGS.items():
            s = scale_clock(clock, WF, a, depth=clock.depth, rng=rng)
            beta, _ = states(s, projector(s, WF), WF, depth)
            delta = clock_change(beta - beta.mean(0), amp, turn)
            xs[name] = s[GENES].to_numpy(np.int64)
            xu[name] = thin_clock(np.round(s[GENES].to_numpy(float)), WF.X, delta, rng).astype(np.int64)
        out[key] = dict(depth=depth, xs=xs, xu=xu)
        rows.append(dict(tissue=tissue, pool_samples=len(pb)))
    return out, pd.DataFrame(rows)


def thin(x, target, depth, rng):
    """Binomial thinning of x (n x G) to target (n x G); one common factor per gene where x holds fewer reads."""
    fmed = np.median(x / depth[:, None], axis=0)
    p = target / (np.maximum(fmed, 1e-12) * depth[:, None])
    capf = np.minimum(1.0, 1.0 / np.maximum(p.max(0), 1e-12))
    return rng.binomial(x, np.minimum(p * capf, 1.0)), float(np.mean(capf < 1))


def depth_match(s, u, ds, du, groups, rng):
    """Spliced counts thinned gene by gene to the unspliced count level (expected counts from cohort median fractions
    x depth), so that both count types carry the same counting noise and pseudocount compression; returns the thinned
    counts and their per-gene depths (n x G)."""
    es, eu = np.empty(s.shape), np.empty(u.shape)
    for g in np.unique(groups):
        m = groups == g
        es[m] = np.median(s[m] / ds[m, None], axis=0)[None] * ds[m, None]
        eu[m] = np.median(u[m] / du[m, None], axis=0)[None] * du[m, None]
    p = np.clip(eu / np.maximum(es, 1e-12), 0, 1)
    return rng.binomial(s.astype(np.int64), p), ds[:, None] * np.maximum(p, 1e-12)


def pergene(Ys, Yu, s, idx, wf, groups, rng) -> tuple:
    """Per-gene fit on the spliced state: for each gene g, the spliced state without g (leave-gene-out projector)
    is the regressor (its two coordinates: a Fourier fit on the inferred phase that stays defined on a partial arc);
    g's spliced and unspliced log fractions are regressed on it, b = beta_cos + i beta_sin, and the lead is
    arg sum_g b_s,g conj(b_u,g) (amplitude-weighted circular mean of the gene lags). Returns (lead, kappa, gate p,
    bootstrap leads)."""
    genes = list(wf.genes)
    P = within(groups, rng, N_PERM)
    B = boot(groups, rng, N_BOOT)
    c = cp = cb = 0
    ss = su = 0.0
    sup = np.zeros(N_PERM)
    for j, g in enumerate(genes):
        o = [k for k in range(len(genes)) if k != j]
        w = wf.drop([g])
        Lj = projector(pd.DataFrame(s[:, [idx[k] for k in o]], columns=[genes[k] for k in o]), w)
        X = Ys[:, o] @ Lj.T
        H = np.linalg.pinv(X)
        bs, bu = H @ Ys[:, j], H @ Yu[:, j]
        bs, bu = bs[0] + 1j * bs[1], bu[0] + 1j * bu[1]
        c = c + bs * np.conj(bu)
        ss, su = ss + abs(bs) ** 2, su + abs(bu) ** 2
        bp = np.einsum("in,pn->pi", H, Yu[P, j])
        bp = bp[:, 0] + 1j * bp[:, 1]
        cp = cp + bs * np.conj(bp)
        sup += np.abs(bp) ** 2
        Xb = X[B]
        XtX = np.einsum("bni,bnj->bij", Xb, Xb)
        bbs = np.linalg.solve(XtX, np.einsum("bni,bn->bi", Xb, Ys[B, j])[..., None])[..., 0]
        bbu = np.linalg.solve(XtX, np.einsum("bni,bn->bi", Xb, Yu[B, j])[..., None])[..., 0]
        cb = cb + (bbs[:, 0] + 1j * bbs[:, 1]) * np.conj(bbu[:, 0] + 1j * bbu[:, 1])
    k = abs(c) / np.sqrt(ss * su)
    kp = np.abs(cp) / np.sqrt(ss * sup)
    return (float(np.angle(c) * HOURS), float(k), float((1 + np.sum(kp >= k)) / (1 + N_PERM)),
            np.angle(cb) * HOURS)


def stats(s, u, obs, Z, groups, wf, rng, estimator="plain", depth_groups=None) -> dict:
    """estimator: 'plain' (registered), 'matched' (spliced thinned gene by gene to the unspliced level within the
    depth cohorts `depth_groups`, default `groups`), 'pergene' (per-gene fits on the leave-gene-out spliced state)."""
    match = estimator == "matched"
    idx = [GENES.index(g) for g in wf.genes]
    L = projector(pd.DataFrame(s[:, idx], columns=list(wf.genes)), wf)
    ds, du = obs.depth_s.to_numpy(float)[:, None], obs.depth_u.to_numpy(float)[:, None]
    if match:
        s, ds = depth_match(s, u, ds[:, 0], du[:, 0], groups if depth_groups is None else depth_groups, rng)

    def z(x, depth):
        depth = np.broadcast_to(depth, x.shape)
        y = dev.residual(np.log((x[:, idx] + 0.5) / depth[:, idx]), Z)
        b = (y - y.mean(0)) @ L.T
        return b[:, 0] + 1j * b[:, 1]

    if estimator == "pergene":
        def Y(x, depth):
            y = dev.residual(np.log((x[:, idx] + 0.5) / np.broadcast_to(depth, x.shape)[:, idx]), Z)
            return y - y.mean(0)
        d, kap, gp, draws = pergene(Y(s, ds), Y(u, du), s, idx, wf, groups, rng)
    else:
        zs, zu = z(s, ds), z(u, du)
        d, kap, gp = lead(zu, zs), kappa(zu, zs), gate_p(zu, zs, groups, rng, N_PERM)
        B = boot(groups, rng, N_BOOT)
        draws = np.angle((zu[B] * np.conj(zs[B])).sum(1)) * HOURS
    lo, hi = interval(d, draws)
    phi = (wrap(draws - d) + d) / HOURS
    T = np.where(phi <= 0, np.inf, np.where(phi >= np.pi / 2, 0.0,
                                            24 * np.tan(PHI_REF) / np.tan(np.clip(phi, 1e-9, np.pi / 2 - 1e-9))))
    return dict(kappa=kap, gate_p=gp, delta=d, lo=lo, hi=hi,
                T_lo=float(np.quantile(T, .05)), T_hi=float(np.quantile(T, .95)))


def one(task) -> list:
    comp, i = task
    rng = np.random.default_rng([SEED, i, sum(map(ord, comp))])
    obs, ts, tu, pool = S["obs"][comp], S["ts"][comp], S["tu"][comp], S["pools"][S["tissue"][comp]]
    n = min(len(obs), len(pool["depth"]))
    if n < len(obs):  # tissue pool smaller than the compartment: the pool is used whole (stated)
        keep = np.sort(rng.choice(len(obs), n, replace=False))
        obs, ts, tu = obs.iloc[keep].reset_index(drop=True), ts[keep], tu[keep]
    pick = rng.choice(len(pool["depth"]), n, replace=False)
    depth = pool["depth"][pick]
    Z, groups = covariates(obs), obs.cohort.to_numpy()
    out = []
    for name in [*SETTINGS, "shuffled"]:
        src = "a1" if name == "shuffled" else name
        s, cap_s = thin(pool["xs"][src][pick], ts, depth, rng)
        u, cap_u = thin(pool["xu"][src][pick], tu, depth, rng)
        if name == "shuffled":
            u = u[within(groups, rng, 1)[0]]
        for (w, wf), est in ((x, e) for x in WAVEFORMS.items() for e in ESTIMATORS):
            out.append(dict(compartment=comp, setting=name, waveform=w, estimator=est,
                            sub=i, n=n, capped_s=cap_s, capped_u=cap_u, **stats(s, u, obs, Z, groups, wf, rng, est)))
        if S["check"] and name in ("a1", "a338_3h"):  # direction check: no thinning
            x_s, x_u = pool["xs"][src][pick], pool["xu"][src][pick]
            deep = obs.assign(depth_s=depth, depth_u=depth)
            for est in ("plain", "pergene"):
                out.append(dict(compartment=comp, setting=f"{name}_deep", waveform="all", estimator=est, sub=i, n=n,
                                **stats(x_s, x_u, deep, Z, groups, WF, rng, est)))
    # model-free check: Greenwell AL + NR, real counts thinned to the compartment
    g = S["greenwell"]
    m = min(len(obs), len(g["obs"]))
    gp, fp = rng.choice(len(g["obs"]), m, replace=False), rng.choice(len(obs), m, replace=False)
    fo = obs.iloc[fp].reset_index(drop=True)
    s, _ = thin(g["s"][gp], ts[fp], g["ds"][gp], rng)
    u, _ = thin(g["u"][gp], tu[fp], g["du"][gp], rng)
    # covariates: regimen and the fetal depth cohort (its gene-wise depth pattern is removed as in the fetal design);
    # permutations and bootstrap within regimen x cohort
    go = fo.assign(regimen=g["obs"].group.to_numpy()[gp])
    Zg = dev.design(go, age=None, factors=("regimen", "cohort"))
    gg = (go.regimen + "|" + go.cohort).to_numpy()
    for (w, wf), est in ((x, e) for x in WAVEFORMS.items() for e in ESTIMATORS):
        out.append(dict(compartment=comp, setting="greenwell", waveform=w, estimator=est,
                        sub=i, n=m, **stats(s, u, go, Zg, gg, wf, rng, est, fo.cohort.to_numpy())))
    return out


# ---- readings ------------------------------------------------------------------------------------------------------

def passrates(res) -> pd.DataFrame:
    res = res.assign(a328=(res.gate_p < .05) & (res.lo > 0),
                     daily=(res.gate_p < .05) & (res.lo > 0) & (res.T_lo > 6) & (res.T_lo <= 24) & (res.T_hi >= 24)
                     & (res.T_hi < 96))
    rows = []
    for (c, s, w, e), r in res.groupby(["compartment", "setting", "waveform", "estimator"], sort=False):
        k, m = int(r.a328.sum()), len(r)
        ci = binomtest(k, m).proportion_ci(0.95, method="exact")
        rows.append(dict(compartment=c, setting=s, waveform=w, estimator=e, n=int(r.n.iloc[0]), n_subsamples=m, gate=r.gate_p.lt(.05).mean(),
                         a328_pass=k / m, ci_lo=ci.low, ci_hi=ci.high, daily=r.daily.mean(),
                         median_delta=r.delta.median(), median_width=(r.hi - r.lo).median()))
    return pd.DataFrame(rows)


def decisions(rates) -> pd.DataFrame:
    rows = []
    for (c, w, e), r in rates.groupby(["compartment", "waveform", "estimator"], sort=False):
        v = r.set_index("setting")
        rows.append(dict(compartment=c, waveform=w, estimator=e, n=int(v.loc["a1", "n"]),
                         a328_power=v.loc["a1", "a328_pass"], a328_steady=v.loc["steady", "a328_pass"],
                         a328_power_a075=v.loc["a0.75", "a328_pass"], a328_power_lead061=v.loc["a1_lead0.61", "a328_pass"],
                         a328_shuffled=v.loc["shuffled", "a328_pass"], a328_greenwell=v.loc["greenwell", "a328_pass"],
                         a328_decidable=v.loc["a1", "a328_pass"] >= 0.8,
                         a338_sensitivity=v.loc["a338_24h", "daily"],
                         a338_false=max(v.loc["steady", "daily"], v.loc["a338_3h", "daily"]),
                         a338_decidable=(v.loc["a338_24h", "daily"] >= 0.8)
                         & (max(v.loc["steady", "daily"], v.loc["a338_3h", "daily"]) <= 0.1)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0, sfx = time.time(), "_check" if a.check else ""
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 220)
    design, det = prep()
    design.to_csv(OUT / f"design{sfx}.csv", index=False)
    det.to_csv(OUT / f"detection{sfx}.csv", index=False)
    print(design.groupby(["compartment", "cohort"]).agg(n=("donor", "size"), depth_s=("depth_s", "median"),
                                                       depth_u=("depth_u", "median")).to_string(), flush=True)
    print(det.groupby(["compartment", "kind"])[GENES].median().median(axis=1).round(1).to_string(),
          f"prep {time.time() - t0:.0f}s", sep="\n", flush=True)
    S["pools"], pinfo = pools(np.random.default_rng([SEED, 0]))
    pinfo.to_csv(OUT / f"pools{sfx}.csv", index=False)
    print(pinfo.to_string(index=False), flush=True)
    gw = score_a337.load("greenwell")[0]
    S["greenwell"] = dict(obs=gw["obs"], s=gw["spliced"].reindex(columns=GENES, fill_value=0).to_numpy(np.int64),
                          u=gw["unspliced"].reindex(columns=GENES, fill_value=0).to_numpy(np.int64),
                          ds=gw["spliced"].sum(1).to_numpy(float), du=gw["unspliced"].sum(1).to_numpy(float))
    comps = list(design.compartment.unique())
    S.update(check=a.check,
             obs={c: design[design.compartment == c].reset_index(drop=True) for c in comps},
             tissue={c: design[design.compartment == c].tissue.iloc[0] for c in comps})
    S["ts"] = {c: targets(S["obs"][c], det[det.compartment == c], "spliced") for c in comps}
    S["tu"] = {c: targets(S["obs"][c], det[det.compartment == c], "unspliced") for c in comps}
    n = 3 if a.check else a.n_sub
    tasks, res = [(c, i) for i in range(n) for c in comps], []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, r in enumerate(pool.imap_unordered(one, tasks)):
            res.extend(r)
            if (k + 1) % 50 == 0 or k + 1 == len(tasks):
                print(f"{k + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res).sort_values(["compartment", "setting", "waveform", "estimator", "sub"])
    res.to_csv(OUT / f"subsamples{sfx}.csv", index=False)
    rates = passrates(res)
    rates.to_csv(OUT / f"passrates{sfx}.csv", index=False)
    dec = decisions(rates[~rates.setting.str.endswith("_deep")])
    dec.to_csv(OUT / f"decisions{sfx}.csv", index=False)
    print(rates.round(3).to_string(index=False), dec.round(3).to_string(index=False),
          f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
