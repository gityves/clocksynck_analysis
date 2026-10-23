"""A18 part 2: per-cell clock amplitude lambda and within-sample synchrony R of the compartments (exploratory).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A18, part 2; labelled exploratory there). For a target of the A17
re-analysis (its table, labels, frozen learned anchors) and every compartment k and group g:
- D: within-sample in-plane dispersion from the fraction-scale cross-gene moments of each (donor, k) group of >= 100 cells
  (`clocksync.cells.vonmises.group_stats`), pooled with cells as weights, minus the mean of its value on 10 planes of
  permuted acrophases;
- rho: the coherent per-sample amplitude, rho^2 = s / (1 - R_W^2), with s = z2 V the host-locked between-donor variance of
  the compartment's pseudobulk (= 2 psi c / (1 - c)) and R_W the window synchrony of the group's donors (0.738 for
  hospital-hours donors, 08-18 h; 0 for any-time donors; averaged over the group);
- lambda, R from lambda R = rho and the von Mises model D(lambda, kappa) = D (`split_lr`, local).
Donor bootstrap within groups (default 100): intervals of log lambda and R and, per contrast, of the differences case -
control of log lambda, R and D. Eligible: >= 20 cell groups, median clock UMI per cell >= 2, D above its null, rho defined;
a contrast is eligible when both groups are. Targets: hlca_fibrotic, lung_hlca, ild_natri, islets_t2d, pd, age_mssm, AD_mssm,
AD_radc, AD_seaad_dlpfc. Cells: `$SCCC_DATA/cells/<dataset>/cells.parquet` (per-cell counts of the 14 clock genes), mapped to
compartments with the class rules of the pseudobulk builders (`paper/common/datasets/_atlas`, `islets`, `parkinson`).

    python compute_per_cell.py [--subsample N] [--jobs 4] [--targets NAME ...] [--n-boot B]
Outputs: paper/results/figure_2/output/backbone/a18_cell_split.csv, a18_cell_split_contrasts.csv (legacy columns).
Legacy: legacy/desync/common/cell_split.py, scripts/x5_vm_human.py (split_lr).
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS, group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("atlas_fits")))
sys.path.insert(0, str(DATASETS / "_atlas"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import optimize, special  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.cells import vonmises as VM  # noqa: E402
from _shared import FOLDER_4, Target, frozen  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

MIN_CELLS, N_NULL, N_SPLITS = 100, 10, 30
R_HOSPITAL = float(np.sin(5 * 2 * np.pi / 24) / (5 * 2 * np.pi / 24))  # 08-18 h
HOSPITAL = {"surgical_resection", "biopsy", "brush", "biopsy_cryopreserved"}
DATASET = {"hlca_fibrotic": "hlca_full", "lung_hlca": "hlca_full", "ild_natri": "natri2024_ild", "islets_t2d": "islets_t2d",
           "pd": "parkinson_regions", "age_mssm": "psychad_mssm", "AD_mssm": "psychad_mssm", "AD_radc": "psychad_radc",
           "AD_seaad_dlpfc": "seaad_dlpfc"}
GENES = list(cs.gtex().genes)


# ---- cell -> compartment (the rules of the pseudobulk builders) -----------------------------------------------------

def _builder(name):
    f = DATASETS / name / "build.py"
    spec = importlib.util.spec_from_file_location(f"build_{name}", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _major(df: pd.DataFrame) -> pd.DataFrame:
    m = df.groupby(["donor", "assay"]).size().reset_index(name="n").sort_values("n").drop_duplicates(
        "donor", keep="last").set_index("donor").assay
    return df[df.assay.to_numpy() == m.reindex(df.donor).to_numpy()]


def map_cells(dataset: str, cols=None) -> pd.DataFrame:
    """Cells with `donor` (dataset-prefixed, as the pseudobulk tables), `compartment`, `depth`, `window`, genes."""
    import census as C
    import cortex as CX
    import lung as LU
    c = pd.read_parquet(DATA / "cells" / dataset / "cells.parquet")
    if dataset == "hlca_full":
        cls = c.ann_level_3.astype(str).map(LU.L3)
        cls[c.ann_level_2.astype(str) == "Smooth muscle"] = "smooth_muscle"
        keep = c.tissue_sampling_method.astype(str).isin(LU.HLCA_TISSUE) & (c.tissue.astype(str) != "nose")
        c["compartment"] = "lung:" + cls.astype(str)
        c = c[keep.to_numpy() & cls.notna().to_numpy()]
        c["window"] = np.where(c.tissue_sampling_method.astype(str).isin(HOSPITAL), "hospital", "any")
    elif dataset == "natri2024_ild":
        cls = c.annotation_level_1.astype(str).map(LU.L1)
        cls[c.lineage.astype(str).isin(LU.FIB_LINEAGES)] = "fibroblast"
        c["compartment"] = "lung:" + cls.astype(str)
        c = c[cls.notna().to_numpy()]
        c["window"] = "any"
    elif dataset == "islets_t2d":
        M = _builder("islets")
        cls = C.classify(c.cell_type.astype(str), M.CLASSES)
        c["compartment"] = "islet:" + cls.astype(str)
        c = c[cls.notna().to_numpy()]
        c["window"] = "any"
    elif dataset == "parkinson_regions":
        M = _builder("parkinson")
        cls = C.classify(c.cell_type.astype(str), M.CLASSES)
        reg = c.tissue.astype(str).map(M.REGIONS)
        ok = cls.notna().to_numpy() & reg.notna().to_numpy()
        c["compartment"] = "brain:" + reg.astype(str) + "_" + cls.astype(str)
        c = c[ok]
        c["window"] = "any"
    elif dataset in ("psychad_mssm", "psychad_radc", "seaad_dlpfc"):
        cls = C.classify(c.cell_type.astype(str), CX.BRAIN, ignore_case=False)
        c = c[c.assay.astype(str).str.startswith("10x").to_numpy() & cls.notna().to_numpy()]
        c["compartment"] = "brain:" + cls[c.index].astype(str)
        c["window"] = "any"
    else:
        raise KeyError(dataset)
    c["donor"] = f"{dataset}:" + c.donor_id.astype(str)
    c["assay"] = c.assay.astype(str)
    c = _major(c)
    return c.rename(columns={"raw_sum": "depth"})[["donor", "compartment", "depth", "window"] + GENES].reset_index(drop=True)


# ---- local: not in clocksync ---------------------------------------------------------------------------------------

def split_lr(D, rho, c, LA):
    """lambda, R with lambda R = rho and model D(lambda, kappa) = D; R -> 1 if D <= 0, nan if no rho (legacy
    `x5_vm_human.split_lr`; `clocksync.cells.vonmises` has `solve` for D and |Gamma|^2 but no split given rho)."""
    if not np.isfinite(rho) or rho <= 0:
        return np.nan, np.nan
    if D <= 0:
        return rho, 1.0

    def f(logk):
        k = np.exp(logk)
        R = special.ive(1, k) / special.ive(0, k)
        v = VM.model_stats(rho / R, k, c, LA, mu=np.linspace(0, 2 * np.pi, 12, endpoint=False))[0] - D
        return v if np.isfinite(v) else 1e6
    lo, hi = np.log(0.05), np.log(200.0)  # R from 0.025 to 0.9975
    if f(lo) < 0:  # even full desynchrony gives less dispersion than observed: report the boundary
        k = np.exp(lo)
    elif f(hi) > 0:
        k = np.exp(hi)
    else:
        k = np.exp(optimize.brentq(f, lo, hi, xtol=1e-4))
    R = float(special.ive(1, k) / special.ive(0, k))
    return float(rho / R), R


def coherent_s(pb: pd.DataFrame, anchors, k: str, donors) -> float:
    """s = z2 V of compartment k over `donors` (duplicates allowed, relabelled `<donor>#<j>`), the host refitted on these
    donors (legacy: the pooled host held fixed; see README, deviations)."""
    parts = []
    for j, r in enumerate(donors):
        parts.append(pb[pb.donor == r].assign(donor=f"{r}#{j}"))
    t = pd.concat(parts, ignore_index=True)
    anc = [a for a in anchors if a in set(t.compartment)]
    if len(anc) < 3 or k not in set(t.compartment) or t[t.compartment == k].donor.nunique() < 5:
        return np.nan
    try:
        tab = cs.Dataset(t, anchors=anc, name="cells").fit(n_boot=0, gate=False, targets=[k]).table
    except (ValueError, np.linalg.LinAlgError):
        return np.nan
    r = tab.loc[k]
    return float(r.z2 * r.V) if np.isfinite(r.z2) else np.nan


# ---- per target ------------------------------------------------------------------------------------------------------

def run_target(name: str, n_boot: int, seed: int, args) -> tuple[list, list]:
    limit_threads(1)
    warnings.filterwarnings("ignore")
    T = Target(name, args.subsample, args.seed)
    a = frozen()
    a = a[a.target == name]
    anc = a.compartment[a.anchor.astype(bool)].tolist()
    if len(anc) < 3 or str(a.status.iloc[0]).startswith("no anchor"):
        return [], []
    lab = T.labels if T.labels is not None else pd.Series("all", index=T.pb.donor.unique())
    groups = ("control", "case") if T.labels is not None else ("all",)
    pb = T.pb[T.pb.donor.isin(lab.index)]
    comps = sorted(pb.compartment.unique())
    cells = map_cells(DATASET[name])
    pairs = set(zip(pb.donor, pb.compartment))
    cells = cells[[p in pairs for p in zip(cells.donor, cells.compartment)]]
    wf = cs.gtex()
    cw = wf.amp * np.exp(-1j * wf.psi)
    rng = np.random.default_rng(seed)
    LA = VM._split_projectors(cw, N_SPLITS, rng)  # local: private constructor of the split projectors
    LAp = [VM._split_projectors((p := wf.permuted(rng)).amp * np.exp(-1j * p.psi), N_SPLITS, rng) for _ in range(N_NULL)]
    win = cells.groupby("donor").window.agg(lambda s: s.mode().iloc[0])
    stats = []
    for (r, k), g in cells.groupby(["donor", "compartment"]):
        if len(g) < MIN_CELLS:
            continue
        C_, N_ = g[GENES].to_numpy(float), g.depth.to_numpy(float)
        s0 = VM.group_stats(C_, N_, LA, rng=rng)
        stats.append(dict(donor=r, compartment=k, n=s0["n"], D=s0["D"],
                          D_null=float(np.mean([VM.group_stats(C_, N_, L, rng=rng)["D"] for L in LAp])),
                          clock_umi=float(C_.sum() / len(g))))
    S = pd.DataFrame(stats)
    if S.empty:
        return [], []

    def est(k, don):
        x = S[S.compartment == k].set_index("donor")
        x = x.loc[[r for r in don if r in x.index]]
        if len(x) < 5:
            return None
        D = float(np.sum(x.n * (x.D - x.D_null)) / x.n.sum())
        rw = np.mean([R_HOSPITAL if win.get(r) == "hospital" else 0.0 for r in don])
        s = coherent_s(pb, anc, k, list(don))
        rho = np.sqrt(s / (1 - rw ** 2)) if np.isfinite(s) and s > 0 else np.nan
        lam, R = split_lr(D, rho, cw, LA)
        umi = float(x.clock_umi.median())
        return dict(D=D, rho=rho, lam=lam, R=R, rw=rw, groups_cells=len(x), clock_umi=umi,
                    eligible=bool(len(x) >= 20 and umi >= 2 and D > 0 and np.isfinite(rho)))

    rows = []
    for k in comps:
        for g in groups:
            don = lab.index[lab == g].to_numpy()
            e = est(k, don)
            if e is None:
                continue
            boots = []
            for _ in range(n_boot):
                b = est(k, rng.choice(don, len(don)))
                if b is not None:
                    boots.append((np.log(b["lam"]) if b["lam"] > 0 else np.nan, b["R"], b["D"]))
            boots = np.array(boots, float) if boots else np.full((0, 3), np.nan)
            q = lambda j, p: np.nanquantile(boots[:, j], p) if len(boots) else np.nan
            rows.append(dict(target=name, compartment=k, group=g, **e, log_lam_lo=q(0, 0.025), log_lam_hi=q(0, 0.975),
                             R_lo=q(1, 0.025), R_hi=q(1, 0.975), boot=boots))
    return rows, contrasts(rows)


def contrasts(rows) -> list:
    df = pd.DataFrame(rows)
    out = []
    if df.empty:
        return out
    for (t, k), g in df.groupby(["target", "compartment"]):
        g = g.set_index("group")
        if not {"control", "case"} <= set(g.index):
            continue
        a, b = g.loc["control"], g.loc["case"]
        n = min(len(a.boot), len(b.boot))
        dl, dR, dD = (b.boot[:n, j] - a.boot[:n, j] for j in range(3))
        q = lambda x, p: np.nanquantile(x, p) if n else np.nan
        out.append(dict(target=t, compartment=k, d_log_lam=np.log(b.lam) - np.log(a.lam), d_R=b.R - a.R,
                        d_log_lam_lo=q(dl, 0.025), d_log_lam_hi=q(dl, 0.975), d_R_lo=q(dR, 0.025), d_R_hi=q(dR, 0.975),
                        R_control=a.R, R_case=b.R, lam_control=a.lam, lam_case=b.lam, rw_control=a.rw, rw_case=b.rw,
                        eligible=bool(a.eligible and b.eligible), D_control=a.D, D_case=b.D, d_D=b.D - a.D,
                        d_D_lo=q(dD, 0.025), d_D_hi=q(dD, 0.975)))
    return out


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--targets", nargs="*", default=list(DATASET))
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_boot = args.n_boot if args.n_boot is not None else (5 if run["quick"] else 100)
    out = out_dir(FOLDER_4, run["quick"])
    todo = [t for t in args.targets if (DATA / "cells" / DATASET[t] / "cells.parquet").exists()]
    res = Parallel(n_jobs=args.jobs)(delayed(run_target)(t, n_boot, args.seed, args) for t in todo)
    rows = [r for x, _ in res for r in x]
    new = pd.DataFrame(rows).drop(columns="boot")
    new.to_csv(out / "a18_cell_split.csv", index=False)
    c = pd.DataFrame([r for _, y in res for r in y])
    c.to_csv(out / "a18_cell_split_contrasts.csv", index=False)
    print(new.groupby("target").eligible.agg(["sum", "size"]), flush=True)


if __name__ == "__main__":
    main()
