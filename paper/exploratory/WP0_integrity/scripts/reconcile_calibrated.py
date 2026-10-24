"""WP0 integrity, item 2 part 3: is the A87 gap (malignant vs normal epithelium) a depth effect, and a simulation calibrated to the
real normal epithelium. Post hoc, exploratory (no registration). Builds on `reconcile_malignant.py` (helpers imported, not modified).

    python reconcile_calibrated.py downsample [--subsample N] [--seeds 20] [--jobs 4]
    python reconcile_calibrated.py calibrate  [--subsample N] [--reps 4]   (grid search of the three free parameters)
    python reconcile_calibrated.py sim        [--subsample N] [--reps 50]  (coherent / incoherent / 3 h jitter malignant at real depth)

downsample: per patient (the 79 of A87), the cells of the normal epithelium (and of the tumour fibroblasts) are subsampled to the
number of cells of the patient's malignant group and thinned binomially cell by cell (clock-gene counts only; the cell depth is kept) so that the sorted clock-UMI-per-cell
distribution equals the malignant one (no upsampling: cells below target keep their counts); 20 seeds. Recomputed on the
thinned compartment: A86 share, A87 coherence and alignment (host as in `compute_decomposition`; fibroblast: host = endothelial + mural,
reference = full-depth normal epithelium), per-cell D (cell-weighted mean of D minus its 10-plane permuted-acrophase null,
groups of >= 100 cells) on the patients with >= 100 normal and malignant cells.
calibrate/sim: the simulation of `reconcile_malignant` plus three parameters absent there: B (between-patient modulation of the
epithelial gene fractions, B Re(c_g e^{i theta_p}) in the log fraction; per-cell dispersion untouched), H (amplitude of the host
pseudobulks) and sd_off (SD of the epithelial state phase around the host phase, h). Host phases mu_p uniform 08-18 h;
cells von Mises(theta_p = mu_p + offset, kappa 2), lam 3. Calibrated on the normal epithelium (share, coherence, alignment, coupling
within ~1 SE of the real values, D reported); alignment is the Poisson-corrected |ccorr|, which can exceed 1 by sampling noise of its corrected denominator; values above 3 (denominator near zero) are set to NA.
Outputs: output/{malignant_downsample,sim_calibration,sim_calibration_grid,incoherence_sim_calibrated}.csv (.../subsample/ for checks).
"""

from __future__ import annotations

import argparse
import itertools
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import reconcile_malignant as rm  # noqa: E402
from reconcile_malignant import (CW, GENES, MAL, NORM, FIB, WF, ccs, cdec, cs, data, sh, VM, limit_threads)  # noqa: E402

OUT = HERE.parents[0] / "output"
SEED = 3030
LAM, KAPPA = 3.0, 2.0
DEFAULT = dict(B=3.0, H=3.0, sd_off=2.0)


# ---- real cells ---------------------------------------------------------------------------------------------------------

def real_patients(subsample, seed):
    pb, _ = data.load("cancer_sc")
    pts = sh.sc_builder().paired_epithelium(pb)
    if subsample:
        pts = sorted(np.random.default_rng(SEED + seed).permutation(pts)[:subsample])
    return pb[pb.donor.isin(pts)].reset_index(drop=True), set(pts)


def load_cells(pts: set) -> dict:
    """{(donor, compartment): (counts, depth)} of the malignant, normal epithelium and tumour fibroblast cells of the patients."""
    cc, scb = data._builder("cells_cancer"), data._builder("cancer_sc")
    out = {}
    for s in sorted({p.split("|")[0] for p in pts}):
        cells = cc.load(s, "dispersion")
        if "barcode" not in cells:
            continue
        calls = scb.S.read_calls(s).rename("cnv_call").reset_index()
        for comp, r, g in ccs.cell_groups(cells, calls):
            if f"{s}|{r}" in pts and comp in (MAL, NORM, FIB):
                out[(f"{s}|{r}", comp)] = (g[GENES].to_numpy(float), g.depth.to_numpy(float))
    return out


def thin_to(Cn, Nn, Cm, rng):
    """Subsample the cells of a group to the malignant cell number and thin cell by cell to the malignant clock-UMI distribution."""
    n = min(len(Cm), len(Cn))
    idx = rng.choice(len(Cn), n, replace=False)
    C, N = Cn[idx], Nn[idx]
    u = C.sum(1)
    target = np.empty(n)
    target[np.argsort(u, kind="stable")] = np.quantile(Cm.sum(1), (np.arange(n) + 0.5) / n)
    f = np.where(u > 0, np.minimum(1.0, target / np.maximum(u, 1)), 1.0)
    return rng.binomial(C.astype(int), f[:, None]).astype(float), N  # depth kept: only the clock-gene counts are thinned (lower clock fraction)


def d_pooled(groups, LA, LAp, rng):
    """Cell-weighted mean of D minus its permuted-acrophase null over (C, N) groups of >= 100 cells."""
    num = den = 0.0
    for C, N in groups:
        if len(N) < ccs.MIN_CELLS:
            continue
        D = VM.group_stats(C, N, LA, rng=rng)["D"]
        Dn = float(np.mean([VM.group_stats(C, N, L, rng=rng)["D"] for L in LAp]))
        num += len(N) * (D - Dn)
        den += len(N)
    return num / den if den else np.nan, den


def row_of(base: pd.Series, C, N):
    r = base.copy()
    r[GENES] = C.sum(0)
    r["depth"] = float(N.sum())
    r["n_cells"] = len(N)
    return r


def down_seed(seed, pb, cells, target_comp):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(SEED + seed)
    LA, LAp = rm.planes()
    P0 = sh.plane()
    rows, gD, nclock, ncell, thinned = [], [], [], [], 0
    pbi = pb.set_index(["donor", "compartment"], drop=False)
    new = []
    for (donor, comp), row in pbi.iterrows():
        if comp != target_comp or (donor, comp) not in cells or (donor, MAL) not in cells:
            new.append(row)
            continue
        Cn, Nn = cells[(donor, comp)]
        Cm, _ = cells[(donor, MAL)]
        C, N = thin_to(Cn, Nn, Cm, rng)
        new.append(row_of(row, C, N))
        thinned += 1
        nclock.append(C.sum() / len(N))
        ncell.append(len(N))
        if len(N) >= ccs.MIN_CELLS and len(Cm) >= ccs.MIN_CELLS:
            gD.append((donor, C, N))
    pbs = pd.DataFrame(new).reset_index(drop=True)
    d = pbs[pbs.compartment == target_comp]
    out = dict(seed=seed, compartment=target_comp, share=rm.share_of(d, P0), clock_umi_per_cell=float(np.mean(nclock)),
               cells_per_patient=float(np.mean(ncell)), n_patients=len(d), n_thinned=thinned)
    if target_comp == NORM:
        samples, host = cdec.single_cell_samples(pbs)
        r = cs.decompose(samples, unit="donor", condition="condition", levels=("normal", "malignant"), strata="stratum", host=host, n_boot=3, seed=seed)
        out.update(coherence=r.loc["coherence_ref", "estimate"], alignment=r.loc["align_ref", "estimate"])
    else:
        cls = pbs.compartment.str.split(":").str[1]
        fs = pbs[pbs.compartment.isin([NORM, FIB])]
        fs = fs.assign(condition=fs.compartment.map({NORM: "normal", FIB: "fibroblast"}))
        both = set(fs.donor[fs.condition == "normal"]) & set(fs.donor[fs.condition == "fibroblast"])
        hostf = pbs[pbs.donor.isin(both) & cls.isin(("endothelial", "mural"))]
        fs = fs[fs.donor.isin(both & set(hostf.donor))]
        r = cs.decompose(fs, unit="donor", condition="condition", levels=("normal", "fibroblast"), strata="stratum", host=hostf, n_boot=3, seed=seed)
        out.update(coherence=r.loc["coherence_case", "estimate"], alignment=r.loc["align_case", "estimate"])
    out["D"], out["D_cells"] = d_pooled([(C, N) for _, C, N in gD], LA, LAp, rng)
    return out


def run_downsample(args):
    out = OUT / ("subsample" if args.subsample else "")
    out.mkdir(parents=True, exist_ok=True)
    pb, pts = real_patients(args.subsample, args.seed)
    cells = load_cells(pts)
    LA, LAp = rm.planes()
    rng = np.random.default_rng(0)
    # real references on the patients with >= 100 normal and malignant cells
    ok = [d for d in pts if (d, NORM) in cells and (d, MAL) in cells and len(cells[(d, NORM)][1]) >= 100 and len(cells[(d, MAL)][1]) >= 100]
    real = []
    P0 = sh.plane()
    for nm, comp in (("real_malignant", MAL), ("real_normal_full_depth", NORM)):
        grp = [cells[(d, comp)] for d in ok]
        D, w = d_pooled(grp, LA, LAp, rng)
        rows_ = pb[pb.compartment == comp]
        real.append(dict(condition=nm, compartment=comp, D=D, D_cells=w, share=rm.share_of(rows_, P0),
                         clock_umi_per_cell=float(np.mean([cells[(d, comp)][0].sum() / len(cells[(d, comp)][1]) for d in pts if (d, comp) in cells])),
                         cells_per_patient=float(np.mean([len(cells[(d, comp)][1]) for d in pts if (d, comp) in cells])), n_patients=len(rows_)))
    samples, host = cdec.single_cell_samples(pb)
    r = cs.decompose(samples, unit="donor", condition="condition", levels=("normal", "malignant"), strata="stratum", host=host, n_boot=3, seed=1)
    real[0].update(coherence=r.loc["coherence_case", "estimate"], alignment=r.loc["align_case", "estimate"])
    real[1].update(coherence=r.loc["coherence_ref", "estimate"], alignment=r.loc["align_ref", "estimate"])
    tabs = [pd.DataFrame(real)]
    for tc in (NORM, FIB):
        res = Parallel(n_jobs=args.jobs)(delayed(down_seed)(s, pb, cells, tc) for s in range(args.seeds))
        r = pd.DataFrame(res)
        num = [k for k in r.columns if k not in ("seed", "compartment")]
        tag = tc.split(":")[1]
        tabs.append(pd.DataFrame({**r[num].mean().to_dict(), "condition": f"downsampled_{tag}_mean", "compartment": tc}, index=[0]))
        tabs.append(pd.DataFrame({**r[num].std().to_dict(), "condition": f"downsampled_{tag}_sd", "compartment": tc}, index=[0]))
        tabs.append(r.assign(condition=f"seed_{tag}"))
    t = pd.concat(tabs, ignore_index=True)
    t.to_csv(out / "malignant_downsample.csv", index=False)
    pd.set_option("display.width", 220)
    print(t[~t.condition.str.startswith("seed")].drop(columns=["seed"], errors="ignore").round(3).to_string(index=False), flush=True)


# ---- calibrated simulation ----------------------------------------------------------------------------------------------

def sim_rep(rep: int, mode: str, par: dict, rows: pd.DataFrame, learned: dict) -> dict:
    """One replicate; mode of the malignant gene phases: coherent, incoherent (uniform) or jitter3h (SD 3 h)."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    ss = np.random.SeedSequence([SEED, rep])
    r_mu, r_off, r_host, r_norm, r_mal, r_psi, r_stat = (np.random.default_rng(x) for x in ss.spawn(7))
    pts = sorted(rows.donor.unique())
    mu = dict(zip(pts, r_mu.uniform(8, 18, len(pts)) * 2 * np.pi / 24))
    off = dict(zip(pts, r_off.normal(0, par["sd_off"] * 2 * np.pi / 24, len(pts))))
    jit = r_psi.normal(0, 3 * 2 * np.pi / 24, len(GENES))
    unif = r_psi.uniform(0, 2 * np.pi, len(GENES))
    psi = {"coherent": WF.psi, "incoherent": unif, "jitter3h": WF.psi + jit}[mode]
    c_mal = WF.amp * np.exp(-1j * psi)
    LA, LAp = rm.planes()
    new, groups = [], []
    for _, row in rows.iterrows():
        th = mu[row.donor] + off[row.donor]
        p = row[GENES].to_numpy(float) / float(row.depth)
        if row.compartment in (MAL, NORM):
            c = c_mal if row.compartment == MAL else CW
            p = p * np.exp(par["B"] * np.real(c * np.exp(1j * th)))
            n = int(row.n_cells)
            (C, N), = VM.simulate(1, n, c, p, LAM, KAPPA, float(row.depth) / n, [th], r_mal if row.compartment == MAL else r_norm)
        else:
            e = np.exp(par["H"] * np.real(CW * np.exp(1j * mu[row.donor])))
            C, N = r_host.poisson(float(row.depth) * p * e)[None, :], np.array([float(row.depth)])
        r_ = row.copy()
        r_[GENES] = C.sum(0)
        r_["depth"] = float(N.sum())
        new.append(r_)
        if row.compartment in (MAL, NORM) and len(N) >= ccs.MIN_CELLS:
            st = VM.group_stats(C, N, LA, rng=r_stat)
            nul = float(np.mean([VM.group_stats(C, N, L, rng=r_stat)["D"] for L in LAp]))
            groups.append(dict(compartment=row.compartment, study=row.study, resection=row.donor, cells=st["n"], D=st["D"], D_null=nul,
                               clock_umi=float(C.sum() / len(N))))
    pb = pd.DataFrame(new).reset_index(drop=True)
    G = pd.DataFrame(groups)
    P0 = sh.plane()
    samples, host = cdec.single_cell_samples(pb)
    r = cs.decompose(samples, unit="donor", condition="condition", levels=("normal", "malignant"), strata="stratum", host=host, n_boot=3, seed=rep)
    fits = rm.fits_of({st: d.assign(donor=d.donor.str.split("|").str[1]) for st, d in pb.groupby("study")}, learned)
    res = {}
    for nm, comp, key in (("malignant", MAL, "case"), ("normal_epithelium", NORM, "ref")):
        s = {st: (f.loc[comp, "s"] if comp in f.index else np.nan) for st, f in fits.items()}
        g = G[G.compartment == comp]
        est = ccs.estimate(g, s, CW, LA) if len(g) else {}
        al = r.loc[f"align_{key}", "estimate"]
        res[nm] = dict(share=rm.share_of(pb[pb.compartment == comp], P0), coherence=r.loc[f"coherence_{key}", "estimate"],
                       alignment=al if abs(al) <= 3 else np.nan, coupling=rm.med_over(fits, comp, "c"), D=est.get("D", np.nan),
                       lam=est.get("lam", np.nan), R=est.get("R", np.nan), rho=est.get("rho", np.nan))
    return res


def sim_rows(subsample, seed):
    pb, pts = real_patients(subsample, seed)
    rows = rm.sim_context(pb, pts)["rows"]
    f = rows[GENES].div(rows.depth, axis=0)
    rows = rows.copy()
    rows[GENES] = f.groupby([rows.study, rows.compartment]).transform("mean").mul(rows.depth, axis=0).to_numpy()
    return rows, {st: rm.learned_anchors(st) for st in rows.study.unique()}


def real_targets() -> pd.DataFrame:
    r = pd.read_csv(OUT / "malignant_readouts.csv")
    r = r[r.compartment == "normal_epithelium"].set_index("readout")
    se = {"share": r.loc["share", "se"], "coherence": (r.loc["coherence", "hi"] - r.loc["coherence", "lo"]) / 3.92,
          "alignment": (r.loc["alignment", "hi"] - r.loc["alignment", "lo"]) / 3.92, "coupling": r.loc["coupling", "se"],
          "D": r.loc["D", "se"], "lam": r.loc["lambda", "se"]}
    return pd.DataFrame({"real": {k: r.loc["lambda" if k == "lam" else k, "estimate"] for k in se}, "real_se": se})


def run_calibrate(args):
    out = OUT / ("subsample" if args.subsample else "")
    out.mkdir(parents=True, exist_ok=True)
    rows, learned = sim_rows(args.subsample, args.seed)
    axes = ((0.6, 0.8), (3.0, 5.0), (3.0, 5.0, 7.0)) if args.refine else ((0.4, 0.7, 1.0), (1.0, 2.0, 3.0), (2.0, 4.0, 6.0))
    grid = [dict(B=B, H=H, sd_off=o) for B, H, o in itertools.product(*axes)]
    jobs = [(i, rep) for i in range(len(grid)) for rep in range(args.reps)]
    res = Parallel(n_jobs=args.jobs)(delayed(sim_rep)(rep, "coherent", grid[i], rows, learned) for i, rep in jobs)
    recs = [dict(**grid[i], rep=rep, readout=k, value=v) for (i, rep), r in zip(jobs, res) for k, v in r["normal_epithelium"].items()]
    g = pd.DataFrame(recs)
    s = g.groupby(["B", "H", "sd_off", "readout"]).value.agg(mean="mean", sd="std").reset_index()
    tg = real_targets()
    fit = ["share", "coherence", "alignment", "coupling"]
    s = s.join(tg, on="readout")
    s["z"] = (s["mean"] - s.real) / s.real_se
    loss = s[s.readout.isin(fit)].assign(z2=lambda x: x.z ** 2).groupby(["B", "H", "sd_off"]).z2.sum().rename("loss")
    s = s.join(loss, on=["B", "H", "sd_off"])
    gf = out / "sim_calibration_grid.csv"
    if args.refine and gf.exists():  # second, finer grid: appended; the best of both grids is kept
        s = pd.concat([pd.read_csv(gf), s], ignore_index=True).drop_duplicates(["B", "H", "sd_off", "readout"], keep="last")
    s.to_csv(gf, index=False)
    best = s.loc[s.loss.idxmin(), ["B", "H", "sd_off"]].to_dict()
    cal = s[(s.B == best["B"]) & (s.H == best["H"]) & (s.sd_off == best["sd_off"])].assign(within_1se=lambda x: x.z.abs() <= 1)
    cal.to_csv(out / "sim_calibration.csv", index=False)
    pd.set_option("display.width", 200)
    print("best", best)
    print(cal.round(3).to_string(index=False))


def run_sim(args):
    out = OUT / ("subsample" if args.subsample else "")
    cal = pd.read_csv(out / "sim_calibration.csv")
    par = cal[["B", "H", "sd_off"]].iloc[0].to_dict()
    rows, learned = sim_rows(args.subsample, args.seed)
    reps = min(args.reps, 4) if args.subsample else args.reps
    jobs = [(rep, m) for rep in range(reps) for m in ("coherent", "incoherent", "jitter3h")]
    res = Parallel(n_jobs=args.jobs)(delayed(sim_rep)(rep, m, par, rows, learned) for rep, m in jobs)
    recs = [dict(scenario=m, compartment=nm, readout=k, rep=rep, value=v) for (rep, m), r in zip(jobs, res) for nm, d in r.items() for k, v in d.items()]
    sim = pd.DataFrame(recs)
    s = sim.groupby(["scenario", "compartment", "readout"]).value.agg(mean="mean", sd="std", median="median", n_finite=lambda x: int(np.isfinite(x).sum()), n_reps="size").reset_index()
    for k, v in par.items():
        s[k] = v
    s.to_csv(out / "incoherence_sim_calibrated.csv", index=False)
    sim.to_csv(out / "incoherence_runs_calibrated.csv", index=False)
    print(par)
    print(s[s.compartment == "malignant"].pivot(index="readout", columns="scenario", values="mean").round(3).to_string())
    print(s[s.compartment == "normal_epithelium"].pivot(index="readout", columns="scenario", values="mean").round(3).to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("task", choices=["downsample", "calibrate", "sim"])
    ap.add_argument("--subsample", type=int, default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--refine", action="store_true", help="calibrate: finer second grid around the first optimum")
    args = ap.parse_args()
    limit_threads(4)
    {"downsample": run_downsample, "calibrate": run_calibrate, "sim": run_sim}[args.task](args)


if __name__ == "__main__":
    main()
