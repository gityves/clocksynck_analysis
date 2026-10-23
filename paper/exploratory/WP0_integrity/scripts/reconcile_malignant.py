"""WP0 integrity, item 2 part 2: why do A86/A87 read malignant cells as less clock-shaped while the per-cell split (A18 part 2,
contrast 3) reads them as no less synchronous? All readouts side by side on the same 79 patients (the A87 set: malignant and
adjacent-normal epithelium and a stromal host), and what each readout returns when the clock is present at the cell level
but the genes have lost their common phase structure.

Exploratory, post hoc (integrity check; no registration). Readouts, per compartment (malignant, adjacent-normal epithelium,
tumour fibroblasts as reference):
  share      A86 pooled clock-plane share (Poisson diagonal subtracted), unit-bootstrap interval;
  coherence  A87 gene coherence R (leave-gene-out loadings against the waveform phases), unit bootstrap;
  alignment  A87 |complex correlation| of the compartment's state with the patient's stromal host (fibroblast: host =
             endothelial + mural, to keep the compartment out of its own host);
  coupling   anchored coupling c, learned anchors, per cohort refitted on the 79 patients (cohorts with >= 5 of them),
             median over cohorts, delete-one-cohort jackknife SE;
  D, lam, R  per-cell split of `compute_cancer_split` (D minus its acrophase-permuted null; rho^2 = s / (1 - R_W^2) with
             s = 2 psi c / (1 - c) of the refit above), delete-one-cohort jackknife SE.
Simulation (`incoherent_sim.csv`): patients, cell counts, mean depths and gene fractions of the real malignant and normal
groups; every patient has a host phase mu (uniform 08-18 h) shared by its stromal host rows, normal epithelium and malignant
cells (cell phases von Mises(mu, kappa=2), lam=3). Scenarios: `coherent` (malignant cells with the GTEx waveform) and
`incoherent` (malignant gene phases drawn uniformly on the circle, amplitudes unchanged, redrawn per replicate); each also with the malignant depth per cell x10 (`_umi10`: separates blindness from lack of power). Normal epithelium and
host are coherent in both. Same seeds across scenarios, so the replicates are paired.

    python reconcile_malignant.py [--subsample N] [--reps R] [--jobs 4]
Outputs: paper/exploratory/WP0_integrity/output/{malignant_readouts,incoherence_sim}.csv (.../subsample/ for checks).
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_5" / "scripts" / "cancer"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import special  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.cells import vonmises as VM  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.anchors import anchors  # noqa: E402

import _shared as sh  # noqa: E402
import compute_decomposition as cdec  # noqa: E402

spec = importlib.util.spec_from_file_location("ccs", PAPER / "results" / "figure_2" / "scripts" / "backbone" / "compute_cancer_split.py")
ccs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ccs)  # estimate, jk, cell-split constants (read only)

OUT = Path(__file__).resolve().parents[1] / "output"
GROUPS = PAPER / "results" / "figure_2" / "output" / "backbone" / "a18_cancer_groups.csv"
MAL, NORM, FIB = sh.MAL, sh.NORM, "tumour:fibroblast"
COMPS = {"malignant": MAL, "normal_epithelium": NORM, "fibroblast": FIB}
GENES = sh.GENES
SEED = 2026
LAM, KAPPA = 3.0, 2.0
MIN_CELLS, N_NULL, N_SPLITS = ccs.MIN_CELLS, ccs.N_NULL, ccs.N_SPLITS
WF = sh.WF
CW = WF.amp * np.exp(-1j * WF.psi)


def planes():
    """Split projectors of the waveform and of its 10 permuted planes, one stream as `compute_cancer_split`."""
    rng = np.random.default_rng(0)
    LA = VM._split_projectors(CW, N_SPLITS, rng)
    perms = [WF.permuted(rng) for _ in range(N_NULL)]
    return LA, [VM._split_projectors(p.amp * np.exp(-1j * p.psi), N_SPLITS, rng) for p in perms]


# ---- shared estimators -------------------------------------------------------------------------------------------------

def cohort_fit(tab: pd.DataFrame, study: str, anc_names) -> pd.DataFrame | None:
    """Anchored fit of one cohort's pseudobulks on the learned anchors; coupling c and s = 2 psi c / (1 - c) per compartment."""
    anc = [a for a in anc_names if a in set(tab.compartment)]
    if len(anc) < 3:
        return None
    try:
        t = cs.Dataset(tab.reset_index(drop=True), anchors=anc, name=study).fit(n_boot=0, gate=False).table
    except (ValueError, np.linalg.LinAlgError):
        return None
    c, psi = t.coupling.to_numpy(float), t.psi.to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(c < 1, 2 * psi * c / (1 - c), np.nan)
    return pd.DataFrame({"c": c, "s": s}, index=t.index)


def learned_anchors(study: str) -> list:
    """Learned anchors of the frozen A17 table (empty for cohorts that are not in it: gate-failed)."""
    try:
        return list(anchors(f"cancer_{study}", "learned"))
    except KeyError:
        return []


def fits_of(tables: dict, learned: dict) -> dict:
    """{study: DataFrame(c, s)} for the cohorts that fit."""
    out = {}
    for st, tab in tables.items():
        f = cohort_fit(tab, st, learned[st])
        if f is not None:
            out[st] = f
    return out


def med_over(fits: dict, comp: str, col: str, drop=None) -> float:
    v = [f.loc[comp, col] for st, f in fits.items() if st != drop and comp in f.index]
    return float(np.nanmedian(v)) if len(v) and np.isfinite(v).any() else np.nan


def jack(vals_fn, studies) -> float:
    v = np.array([vals_fn(s) for s in studies], float)
    v = v[np.isfinite(v)]
    n = len(v)
    return float(np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2))) if n > 2 else np.nan


def share_of(d: pd.DataFrame, P0) -> float:
    return sh.pooled_share(*sh.centred_units(d, "stratum"), P0) if d.stratum.nunique() and len(d) > 2 else np.nan


def split_estimates(G: pd.DataFrame, fits: dict, comp: str, LA) -> dict:
    """Pooled D (minus null), rho, lam, R of one compartment from its group table, and jackknife SEs (cohorts)."""
    s = {st: (f.loc[comp, "s"] if comp in f.index else np.nan) for st, f in fits.items()}
    g = G[G.compartment == comp]
    if not len(g):
        return {}
    e = ccs.estimate(g, s, CW, LA)
    sts = sorted(g.study.unique())
    se = dict(D=ccs.jk(g, s, lambda gg, ss: ccs.estimate(gg, ss, CW, LA)["D"]),
              R=ccs.jk(g, s, lambda gg, ss: ccs.estimate(gg, ss, CW, LA)["R"]),
              lam=ccs.jk(g, s, lambda gg, ss: np.log(ccs.estimate(gg, ss, CW, LA)["lam"])))
    return dict(e, groups=len(g), cohorts=len(sts), se=se)


# ---- real data ----------------------------------------------------------------------------------------------------------

def real_readouts(subsample, seed, n_boot, n_perm) -> pd.DataFrame:
    limit_threads(4)
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(SEED + seed)
    pb, _ = data.load("cancer_sc")
    pts = cdec.sh.sc_builder().paired_epithelium(pb)
    if subsample:
        pts = sorted(rng.permutation(pts)[:subsample])
    pb = pb[pb.donor.isin(pts)].reset_index(drop=True)
    P0 = sh.plane()
    rows = []

    def add(comp, readout, est, se=np.nan, lo=np.nan, hi=np.nan, n=np.nan, units=np.nan, src=""):
        rows.append(dict(compartment=comp, readout=readout, estimate=est, se=se, lo=lo, hi=hi, n_patients=n, n_units=units, source=src))

    # A86 share
    for name, comp in COMPS.items():
        d = pb[pb.compartment == comp]
        obs = share_of(d, P0)
        null = np.mean([sh.pooled_share(*sh.centred_units(d, "stratum"), sh.plane(WF.permuted(rng))) for _ in range(n_perm)])
        bs = [share_of(d.iloc[rng.integers(0, len(d), len(d))], P0) for _ in range(n_boot)]
        lo, hi = sh.mean_ci(bs)
        add(name, "share", obs, float(np.nanstd(bs)), lo, hi, d.donor.nunique(), len(d), "A86 pooled_share (same patients)")
        add(name, "share_null", null, np.nan, np.nan, np.nan, d.donor.nunique(), len(d), "mean of permuted-acrophase planes")

    # A87 coherence, alignment
    samples, host = cdec.single_cell_samples(pb)
    r = cs.decompose(samples, unit="donor", condition="condition", levels=("normal", "malignant"), strata="stratum", host=host,
                     n_boot=n_boot, seed=SEED + seed)
    n = samples.donor.nunique()
    for nm, key in (("malignant", "case"), ("normal_epithelium", "ref")):
        for out, k in (("coherence", f"coherence_{key}"), ("alignment", f"align_{key}")):
            add(nm, out, r.loc[k, "estimate"], np.nan, r.loc[k, "lo"], r.loc[k, "hi"], n, 2 * n, "A87 decompose")
    for k in ("d_coherence", "d_align", "d_log_plane", "d_log_off"):
        add("malignant-normal_epithelium", k, r.loc[k, "estimate"], np.nan, r.loc[k, "lo"], r.loc[k, "hi"], n, 2 * n, "A87 decompose")
    cls = pb.compartment.str.split(":").str[1]
    fs = pb[pb.compartment.isin([NORM, FIB])]
    fs = fs.assign(condition=fs.compartment.map({NORM: "normal", FIB: "fibroblast"}))
    both = set(fs.donor[fs.condition == "normal"]) & set(fs.donor[fs.condition == "fibroblast"])
    fs = fs[fs.donor.isin(both)]
    hostf = pb[pb.donor.isin(both) & cls.isin(("endothelial", "mural"))]
    fs = fs[fs.donor.isin(set(hostf.donor))]
    if fs.donor.nunique() >= 10:
        rf = cs.decompose(fs, unit="donor", condition="condition", levels=("normal", "fibroblast"), strata="stratum", host=hostf,
                          n_boot=n_boot, seed=SEED + seed)
        nf = fs.donor.nunique()
        for out, k in (("coherence", "coherence_case"), ("alignment", "align_case")):
            add("fibroblast", out, rf.loc[k, "estimate"], np.nan, rf.loc[k, "lo"], rf.loc[k, "hi"], nf, nf, "A87 decompose, host = endothelial + mural")

    # anchored coupling and the per-cell split on the refit
    tables, learned = {}, {}
    for st in sorted(pb.study.unique()):
        t, _ = data.load("cancer_sc", study=st, cnv=True)
        keep = {p.split("|", 1)[1] for p in pts if p.startswith(st + "|")}
        tables[st] = t[t.donor.isin(keep)]
        learned[st] = learned_anchors(st)
        t2, _ = data.load("cancer_sc", study=st, cnv=False)
        tables[st] = pd.concat([tables[st], t2[t2.donor.isin(keep) & ~t2.compartment.isin(set(tables[st].compartment))]])
    fits = fits_of(tables, learned)
    npat = {nm: pb[pb.compartment == comp].donor.nunique() for nm, comp in COMPS.items()}
    for nm, comp in COMPS.items():
        c = med_over(fits, comp, "c")
        sts = [st for st in fits if comp in fits[st].index]
        add(nm, "coupling", c, jack(lambda x: med_over(fits, comp, "c", drop=x), sts), np.nan, np.nan, npat[nm], len(sts), "anchored fit, cohorts with >= 5 patients")
    G = pd.read_csv(GROUPS)
    G = G[[f"{a}|{b}" in set(pts) for a, b in zip(G.study, G.resection)]]
    LA, _ = planes()
    est = {}
    for nm, comp in COMPS.items():
        e = split_estimates(G, fits, comp, LA)
        est[nm] = e
        if not e:
            continue
        n_pat = G[G.compartment == comp][["study", "resection"]].drop_duplicates().shape[0]
        for k, name, src in (("D", "D", "D minus permuted null"), ("rho", "rho", "rho"), ("lam", "lambda", "lambda"), ("R", "R", "R")):
            se = e["se"].get(k, np.nan)
            add(nm, name, e[k], se, np.nan, np.nan, n_pat, e["groups"],
                f"compute_cancer_split estimate ({src}; SE of lambda on the log scale)")
        g = G[G.compartment == comp]
        dn = float(np.sum(g.cells * g.D_null) / g.cells.sum())
        add(nm, "D_raw", e["D"] + dn, np.nan, np.nan, np.nan, np.nan, e["groups"], "pooled D before the null")
        add(nm, "D_null", dn, np.nan, np.nan, np.nan, np.nan, e["groups"], "cell-weighted mean permuted-acrophase D")
    # paired contrast on the common patients
    a, b = G[G.compartment == MAL], G[G.compartment == NORM]
    common = set(zip(a.study, a.resection)) & set(zip(b.study, b.resection))
    key = lambda g: [x in common for x in zip(g.study, g.resection)]
    ga, gb = a[key(a)], b[key(b)]
    both = pd.concat([ga, gb])
    sa = {st: (f.loc[MAL, "s"] if MAL in f.index else np.nan) for st, f in fits.items()}
    sb = {st: (f.loc[NORM, "s"] if NORM in f.index else np.nan) for st, f in fits.items()}
    ea, eb = ccs.estimate(ga, sa, CW, LA), ccs.estimate(gb, sb, CW, LA)

    def diffR(gg, ss):
        return ccs.estimate(gg[gg.compartment == MAL], {v: sa[v] for v in ss}, CW, LA)["R"] - ccs.estimate(gg[gg.compartment == NORM], {v: sb[v] for v in ss}, CW, LA)["R"]
    add("malignant-normal_epithelium", "d_R", ea["R"] - eb["R"], ccs.jk(both, sa, diffR), np.nan, np.nan, len(common), len(both), "paired patients, jackknife over cohorts")
    add("malignant-normal_epithelium", "d_log_lambda", np.log(ea["lam"]) - np.log(eb["lam"]), np.nan, np.nan, np.nan, len(common), len(both), "paired patients")
    cm, cn = med_over(fits, MAL, "c"), med_over(fits, NORM, "c")
    add("malignant-normal_epithelium", "d_coupling", cm - cn, np.nan, np.nan, np.nan, len(common), np.nan, "difference of the cohort medians")
    return pd.DataFrame(rows)


# ---- simulation ---------------------------------------------------------------------------------------------------------

def sim_context(pb: pd.DataFrame, pts) -> dict:
    """Real structure of the 79 patients: rows of every compartment that is a host (stromal) or an epithelium of interest."""
    cls = pb.compartment.str.split(":").str[1]
    use = pb[pb.donor.isin(pts) & (cls.isin(sh.STROMA) | pb.compartment.isin([MAL, NORM]))].reset_index(drop=True)
    return dict(rows=use)


def cells_of(row, mu, c, rng, umi_mult=1.0):
    n = int(row.n_cells)
    N = float(row.depth) / n * umi_mult
    p = row[GENES].to_numpy(float) / float(row.depth)
    (C, d), = VM.simulate(1, n, c, p, LAM, KAPPA, N, [mu], rng)
    return C, d


def sim_replicate(rep: int, scenario_psi: bool, umi_mult: float, rows: pd.DataFrame, learned: dict, G_real_depth, subsample) -> dict:
    """One replicate of one scenario; returns the readouts of malignant and normal epithelium."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    ss = np.random.SeedSequence([SEED, rep])
    r_mu, r_host, r_norm, r_mal, r_psi, r_stat = (np.random.default_rng(x) for x in ss.spawn(6))
    pts = sorted(rows.donor.unique())
    mu = dict(zip(pts, r_mu.uniform(8, 18, len(pts)) * 2 * np.pi / 24))
    if scenario_psi:
        psi = r_psi.uniform(0, 2 * np.pi, len(GENES))
        c_mal = WF.amp * np.exp(-1j * psi)
    else:
        c_mal = CW
        r_psi.uniform(0, 2 * np.pi, len(GENES))
    LA, LAp = planes()
    out_rows, groups = [], []
    for _, row in rows.iterrows():
        m = mu[row.donor]
        if row.compartment == MAL:
            C, N = cells_of(row, m, c_mal, r_mal, umi_mult)
        elif row.compartment == NORM:
            C, N = cells_of(row, m, CW, r_norm)
        else:  # host pseudobulk: Poisson counts with the coherent waveform at the patient's phase, amplitude lam
            p = row[GENES].to_numpy(float) / float(row.depth)
            e = np.exp(LAM * np.real(CW * np.exp(1j * m)))
            C, N = r_host.poisson(float(row.depth) * p * e)[None, :], np.array([float(row.depth)])
        new = row.copy()
        new[GENES] = C.sum(0)
        new["depth"] = float(N.sum())
        out_rows.append(new)
        if row.compartment in (MAL, NORM) and len(N) >= MIN_CELLS:
            st = VM.group_stats(C, N, LA, rng=r_stat)
            nul = float(np.mean([VM.group_stats(C, N, L, rng=r_stat)["D"] for L in LAp]))
            groups.append(dict(compartment=row.compartment, study=row.study, resection=row.donor, cells=st["n"], D=st["D"], D_null=nul,
                               clock_umi=float(C.sum() / len(N))))
    pb = pd.DataFrame(out_rows).reset_index(drop=True)
    G = pd.DataFrame(groups)
    P0 = sh.plane()
    res = {}
    samples, host = cdec.single_cell_samples(pb)
    r = cs.decompose(samples, unit="donor", condition="condition", levels=("normal", "malignant"), strata="stratum", host=host,
                     n_boot=3, seed=rep)
    tables = {}
    for st, d in pb.groupby("study"):
        tables[st] = d.assign(donor=d.donor.str.split("|").str[1])
    fits = fits_of(tables, learned)
    for nm, comp, key in (("malignant", MAL, "case"), ("normal_epithelium", NORM, "ref")):
        s = {st: (f.loc[comp, "s"] if comp in f.index else np.nan) for st, f in fits.items()}
        g = G[G.compartment == comp]
        est = ccs.estimate(g, s, CW, LA) if len(g) else {}
        res[nm] = dict(share=share_of(pb[pb.compartment == comp], P0), coherence=r.loc[f"coherence_{key}", "estimate"],
                       alignment=r.loc[f"align_{key}", "estimate"], coupling=med_over(fits, comp, "c"),
                       D=est.get("D", np.nan), lam=est.get("lam", np.nan), R=est.get("R", np.nan), rho=est.get("rho", np.nan),
                       D_null=float(np.sum(g.cells * g.D_null) / g.cells.sum()) if len(g) else np.nan)
    return res


def simulate_all(reps, jobs, subsample, seed) -> pd.DataFrame:
    pb, _ = data.load("cancer_sc")
    pts = sh.sc_builder().paired_epithelium(pb)
    if subsample:
        pts = sorted(np.random.default_rng(SEED + seed).permutation(pts)[:subsample])
    rows = sim_context(pb, pts)["rows"]
    learned = {st: learned_anchors(st) for st in rows.study.unique()}
    jobs_ = [(rep, scen, um) for rep in range(reps) for um in (1.0, 10.0) for scen in (False, True)]
    res = Parallel(n_jobs=jobs)(delayed(sim_replicate)(rep, scen, um, rows, learned, None, subsample) for rep, scen, um in jobs_)
    recs = []
    for (rep, scen, um), r in zip(jobs_, res):
        for nm, d in r.items():
            for k, v in d.items():
                recs.append(dict(scenario=("incoherent" if scen else "coherent") + ("_umi10" if um > 1 else ""), compartment=nm, readout=k, rep=rep, value=v))
    return pd.DataFrame(recs)


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--reps", type=int, default=50)
    args = ap.parse_args()
    run = cli.setup(args)
    quick = run["quick"]
    out = OUT / ("subsample" if quick else "")
    out.mkdir(parents=True, exist_ok=True)
    n_boot, n_perm = (30, 30) if quick else (500, 500)
    real = real_readouts(args.subsample, args.seed, n_boot, n_perm)
    real.to_csv(out / "malignant_readouts.csv", index=False)
    print(real.round(3).to_string(index=False), flush=True)
    reps = min(args.reps, 4) if quick else args.reps
    sim = simulate_all(reps, args.jobs, args.subsample, args.seed)
    sim.to_csv(out / "incoherence_runs.csv", index=False)
    s = (sim.groupby(["scenario", "compartment", "readout"]).value.agg(mean="mean", sd="std", n_finite=lambda x: int(np.isfinite(x).sum()),
                                                                          n_reps="size").reset_index())
    s.to_csv(out / "incoherence_sim.csv", index=False)
    print(s.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
