"""A18 part 2, contrast 3: per-cell amplitude lambda and synchrony R of tumour and adjacent-normal compartments of the
single-cell cancer cohorts, on the anchored host (exploratory, as labelled in A18).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A18 part 2, contrast 3 "copy-number malignant vs adjacent-normal epithelium";
run note 2026-10-12/13). Estimand: for every compartment k pooled over cohorts, D_k = within-sample in-plane dispersion of the
cells of each (patient, k) group of >= 100 cells (`cells.vonmises.group_stats`), cell-weighted over groups, minus its mean
on 10 acrophase-permuted planes; rho_k^2 = s_k / (1 - R_W^2) with s_k = 2 psi c / (1 - c) of the cohort's anchored fit (learned
anchors; copy-number table for `tumour:malignant` / `tumour:epithelial_nonmalignant`, main table for the others),
R_W = sinc(5 omega) (surgery, 08-18 h), rho = median over cohorts; lambda, R from lambda R = rho and the von Mises model
D(lambda, kappa) = D. Paired contrasts (same patients, both groups present, >= 5 patients): differences of R and of log lambda.
Interval: delete-one-cohort jackknife SE (D and rho re-estimated). Eligible: >= 20 groups, median clock UMI per cell >= 2,
D > 0, rho finite. Cells: `data/cells_cancer` (`dispersion` tables; copy-number split by the stored calls).

    python compute_cancer_split.py [--subsample N] [--jobs 4]
Outputs: paper/results/figure_2/output/backbone/a18_cancer_groups.csv, a18_cancer_split.csv, a18_cancer_split_contrasts.csv
Legacy: legacy/desync/common/cancer_split.py with the groups of legacy/scripts/x5_vm_human.py --cnv.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import optimize, special  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.cells import vonmises as VM  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.anchors import TABLE, anchors  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "backbone"
MIN_CELLS, N_NULL, N_SPLITS = 100, 10, 30
XW = 5 * 2 * np.pi / 24
RW = np.sin(XW) / XW  # 08-18 h
GENES = list(cs.gtex().genes)
MAL, NMAL = "tumour:malignant", "tumour:epithelial_nonmalignant"
CNV_COMP = {"malignant": MAL, "normal": NMAL}
CONTRASTS = [("tumour:epithelial", "adjacent_normal:epithelial"), (MAL, "adjacent_normal:epithelial"),
             (NMAL, "adjacent_normal:epithelial"), (MAL, NMAL), ("tumour:fibroblast", "adjacent_normal:fibroblast"),
             (MAL, "tumour:fibroblast")]


# ---- local: not in clocksync ---------------------------------------------------------------------------------------

def split_lr(D, rho, c, LA):
    """lambda, R with lambda R = rho and model D(lambda, kappa) = D (legacy `x5_vm_human.split_lr`; same as
    backbone/compute_per_cell.py)."""
    if not np.isfinite(rho) or rho <= 0:
        return np.nan, np.nan
    if D <= 0:
        return rho, 1.0

    def f(logk):
        k = np.exp(logk)
        R = special.ive(1, k) / special.ive(0, k)
        v = VM.model_stats(rho / R, k, c, LA, mu=np.linspace(0, 2 * np.pi, 12, endpoint=False))[0] - D
        return v if np.isfinite(v) else 1e6
    lo, hi = np.log(0.05), np.log(200.0)
    k = np.exp(lo) if f(lo) < 0 else np.exp(hi) if f(hi) > 0 else np.exp(optimize.brentq(f, lo, hi, xtol=1e-4))
    R = float(special.ive(1, k) / special.ive(0, k))
    return float(rho / R), R


def study_s(study: str, subsample=None, seed=0) -> dict:
    """s = 2 psi c / (1 - c) per compartment on the study's learned anchors: copy-number table first, the main table
    for the compartments it lacks (`tumour:epithelial`)."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = {}
    for cnv in (True, False):
        pb, _ = data.load("cancer_sc", study=study, cnv=cnv)
        if subsample is not None:
            keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:max(subsample, 12)]
            pb = pb[pb.donor.isin(keep)]
        try:
            anc = [a for a in anchors(f"cancer_{study}", "learned") if a in set(pb.compartment)]
        except KeyError:
            return out
        if len(anc) < 3:  # legacy: a study without 3 anchors in its copy-number table contributes nothing
            return out
        try:
            t = cs.Dataset(pb.reset_index(drop=True), anchors=anc, name=study).fit(n_boot=0, gate=False).table
        except (ValueError, np.linalg.LinAlgError):
            continue
        c, psi = t.coupling.to_numpy(float), t.psi.to_numpy(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            s = np.where(c < 1, 2 * psi * c / (1 - c), np.nan)
        for k, v in zip(t.index, s):
            out.setdefault(k, float(v))
    return out


def estimate(g, s, c, LA):
    w = g.cells.to_numpy(float)
    D = float(np.sum(w * (g.D - g.D_null)) / w.sum())
    rho = np.sqrt(np.nanmedian([s.get(v, np.nan) for v in g.study.unique()]) / (1 - RW ** 2))
    lam, R = split_lr(D, rho, c, LA)
    return dict(D=D, rho=float(rho), lam=lam, R=R)


def jk(g, s, fn):
    st = g.study.unique()
    if len(st) < 3:
        return np.nan
    v = np.array([fn(g[g.study != x], {k: s[k] for k in s if k != x}) for x in st], float)
    v = v[np.isfinite(v)]
    n = len(v)
    return float(np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2))) if n > 2 else np.nan


# ---- groups ---------------------------------------------------------------------------------------------------------

def cell_groups(cells: pd.DataFrame, calls):
    """(compartment, patient, cells) per group; tumour epithelium also split by copy-number call."""
    for (comp, r), g in cells.groupby(["compartment", "resection"]):
        yield comp, r, g
    if calls is None:
        return
    e = cells[cells.compartment == "tumour:epithelial"].merge(calls, on=["sample", "barcode"], how="left")
    e = e[e.cnv_call.isin(list(CNV_COMP))]
    for (call, r), g in e.groupby(["cnv_call", "resection"]):
        yield CNV_COMP[call], r, g


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    out = out_dir(FOLDER, quick)
    cc = data._builder("cells_cancer")
    sc = data._builder("cancer_sc")
    reg = pd.read_csv(TABLE)
    studies = sorted(t.removeprefix("cancer_") for t in reg.target.unique() if t.startswith("cancer_"))
    studies = [s for s in studies if s in set(sc.CNV_STUDIES)]
    wf = cs.gtex()
    c = wf.amp * np.exp(-1j * wf.psi)
    rng = np.random.default_rng(0)  # one stream in the legacy order (D uses random half splits of the cells)
    LA = VM._split_projectors(c, N_SPLITS, rng)  # local: private constructor of the split projectors
    perms = [wf.permuted(rng) for _ in range(N_NULL)]
    LAp = [VM._split_projectors(p.amp * np.exp(-1j * p.psi), N_SPLITS, rng) for p in perms]
    groups = []
    used = []
    for s in studies:
        cells = cc.load(s, "dispersion")
        if "barcode" not in cells:  # tables of the first extraction have no barcodes: copy-number split impossible
            print(s, "no barcodes in the legacy table: rebuild with data/cells_cancer/build.py", flush=True)
            continue
        if quick:
            keep = np.random.default_rng(args.seed).permutation(sorted(cells.resection.unique()))[:max(args.subsample, 12)]
            cells = cells[cells.resection.isin(keep)]
            if len(used) >= 3:
                break
        calls = sc.S.read_calls(s).rename("cnv_call").reset_index()
        for comp, r, g in cell_groups(cells, calls):
            if len(g) < MIN_CELLS:
                continue
            C_, N_ = g[GENES].to_numpy(float), g.depth.to_numpy(float)
            st = VM.group_stats(C_, N_, LA, rng=rng)
            groups.append(dict(compartment=comp, study=s, resection=r, cells=st["n"], D=st["D"],
                               D_null=float(np.mean([VM.group_stats(C_, N_, L, rng=rng)["D"] for L in LAp])),
                               clock_umi=float(C_.sum() / len(g))))
        used.append(s)
        print(s, flush=True)
    G = pd.DataFrame(groups)
    G.to_csv(out / "a18_cancer_groups.csv", index=False)
    S = dict(zip(used, Parallel(n_jobs=args.jobs)(delayed(study_s)(v, args.subsample, args.seed) for v in used)))

    rows = []
    for k, g in G.groupby("compartment"):
        s = {v: S[v].get(k, np.nan) for v in S}
        e = estimate(g, s, c, LA)
        rows.append(dict(compartment=k, groups=len(g), studies=g.study.nunique(), clock_umi=float(g.clock_umi.median()), **e,
                         R_se=jk(g, s, lambda gg, ss: estimate(gg, ss, c, LA)["R"]),
                         log_lam_se=jk(g, s, lambda gg, ss: np.log(estimate(gg, ss, c, LA)["lam"])),
                         eligible=bool(len(g) >= 20 and g.clock_umi.median() >= 2 and e["D"] > 0 and np.isfinite(e["rho"]))))
    d = pd.DataFrame(rows)
    d.to_csv(out / "a18_cancer_split.csv", index=False)
    print(d.round(3).to_string(index=False))
    res = []
    for ka, kb in CONTRASTS:
        ga, gb = G[G.compartment == ka], G[G.compartment == kb]
        common = set(zip(ga.study, ga.resection)) & set(zip(gb.study, gb.resection))
        if len(common) < 5:
            continue
        key = lambda g: [x in common for x in zip(g.study, g.resection)]
        ga, gb = ga[key(ga)], gb[key(gb)]
        sa, sb = ({v: S[v].get(k, np.nan) for v in S} for k in (ka, kb))
        both = pd.concat([ga, gb])

        def diff(gg, ss, what):
            ea = estimate(gg[gg.compartment == ka], {v: sa[v] for v in ss}, c, LA)
            eb = estimate(gg[gg.compartment == kb], {v: sb[v] for v in ss}, c, LA)
            return ea[what] - eb[what] if what == "R" else np.log(ea["lam"]) - np.log(eb["lam"])

        ea, eb = estimate(ga, sa, c, LA), estimate(gb, sb, c, LA)
        res.append(dict(a=ka, b=kb, resections=len(common), studies=both.study.nunique(), R_a=ea["R"], R_b=eb["R"],
                        d_R=ea["R"] - eb["R"], d_R_se=jk(both, sa, lambda gg, ss: diff(gg, ss, "R")),
                        lam_a=ea["lam"], lam_b=eb["lam"], d_log_lam=np.log(ea["lam"]) - np.log(eb["lam"]),
                        d_log_lam_se=jk(both, sa, lambda gg, ss: diff(gg, ss, "lam")), rho_a=ea["rho"], rho_b=eb["rho"]))
    o = pd.DataFrame(res)
    o.to_csv(out / "a18_cancer_split_contrasts.csv", index=False)
    print(o.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
