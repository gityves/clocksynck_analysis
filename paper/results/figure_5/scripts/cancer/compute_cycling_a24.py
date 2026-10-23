"""A24, cancer arm: is the clock of cycling epithelial cells delayed?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A24 and its note of 2026-10-13). Cell-cycle phase per epithelial cell as
Seurat's CellCycleScoring (S and G2M module scores, Tirosh 2016 lists; S if S > G2M and S > 0, G2M if G2M >= S and G2M > 0, else
G1; cycling = S or G2M). The epithelial cells of tumour and adjacent-normal samples are split by phase (parts G1, S, G2M,
cycling; >= 20 cells per patient and part) and added to the study's non-epithelial pseudobulks; host = the study's learned
A17 anchors minus the epithelial ones (>= 3 are needed by `Dataset.fit`). Estimand: y5 = phase of a part minus phase of G1 of the
same compartment (h; < 0: the cycling clock reads earlier, "delayed"), per study for tumour and adjacent-normal epithelium
(patients with both parts >= 6), pooled over studies by DerSimonian-Laird (rows with SE <= 4 h), and tumour minus normal.
Interval: patient bootstrap (500, host refitted), shared replicates for the difference. Sensitivity of the note (b):
`--definition quartile` (cycling = top quartile of max(S, G2M) within study x tissue kind and > 0, S vs G2M by the larger score; G1 =
at or below the median and Seurat G1; other cells dropped). The colon arm with known time and the author-label definition (a)
are mouse analyses, not here.

    python compute_cycling_a24.py [--subsample N] [--jobs 4] [--definition seurat|quartile] [--build-missing]
Outputs (paper/results/figure_5/output/cancer/): a24_cancer[_quartile].csv, a24_cancer_pooled[_quartile].csv
Legacy: legacy/desync/common/a24_cycling.py (`cancer`, `_cancer_study`).
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

import clocksync as cs  # noqa: E402
from clocksync.anchors import family  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.anchors import TABLE, anchors  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer"
PARTS = ("S", "G2M", "cycling")
SE_MAX = 4.0
GENES = list(cs.gtex().genes)


def wrap(a):
    return (np.asarray(a, float) + 12.0) % 24.0 - 12.0


def call_phase(S, G2M) -> np.ndarray:
    S, G = np.asarray(S, float), np.asarray(G2M, float)
    return np.where((S > G) & (S > 0), "S", np.where((G >= S) & (G > 0), "G2M", "G1"))


def relabel(cells: pd.DataFrame, how: str) -> pd.DataFrame:
    """'seurat' (as A24) or 'quartile' (note b): top quartile of max(S, G2M) within tissue kind = cycling, at or below the
    median and Seurat G1 = G1, others dropped."""
    c = cells.copy()
    if how == "quartile":
        sg = np.where(c.S > c.G2M, "S", "G2M")
        m = np.maximum(c.S, c.G2M)
        q75 = m.groupby(c.kind).transform(lambda x: x.quantile(0.75))
        q50 = m.groupby(c.kind).transform("median")
        c["phase"] = np.where((m >= q75) & (m > 0), sg, np.where((m <= q50) & (c.phase == "G1"), "G1", "drop"))
    return c


def split_parts(cells: pd.DataFrame, key: str) -> pd.DataFrame:
    """Each cell once per part it belongs to: its phase (G1, S, G2M) and, if S or G2M, 'cycling'."""
    cells = cells[cells.phase != "drop"]
    out = pd.concat([cells, cells[cells.phase != "G1"].assign(phase="cycling")], ignore_index=True)
    out[key] = out[key] + "|" + out.phase
    return out


def statistics(ds: cs.Dataset, bases, comps):
    """y5 (h) and log amplitude ratio of every (base, part) against G1 of the same base on one fit of the host."""
    t = ds.fit(n_boot=0, gate=False, targets=comps).table
    out = {}
    for b in bases:
        for p in PARTS:
            a, g = f"{b}|{p}", f"{b}|G1"
            out[(b, p)] = ((float(wrap(t.offset_h[a] - t.offset_h[g])), float(t.log_amp_rel[a] - t.log_amp_rel[g]))
                           if a in t.index and g in t.index and t.z2[a] > 0 and t.z2[g] > 0  # undefined when z2 <= 0 (as legacy)
                           else (np.nan, np.nan))
    return out


def study_arm(study: str, seed: int, how: str, n_boot: int, subsample, sseed, build_missing: bool) -> pd.DataFrame:
    limit_threads(1)
    warnings.filterwarnings("ignore")
    name = f"cancer_{study}"
    cc = data._builder("cells_cancer")
    try:
        cells = cc.load(study, "epithelial", legacy=True) if build_missing or _has(study) else None
    except FileNotFoundError:
        cells = None
    if cells is None:
        return pd.DataFrame([dict(study=name, status="no state cells")])
    host = cc.load(study, "host")
    if subsample is not None:
        keep = np.random.default_rng(sseed).permutation(sorted(host.resection.unique()))[:max(subsample, 12)]
        cells, host = cells[cells.resection.isin(keep)], host[host.resection.isin(keep)]
    cells["phase"] = call_phase(cells.S, cells.G2M)
    cells = relabel(cells, how)
    cells["compartment"] = cells.kind + ":epithelial"
    sp_ = split_parts(cells, "compartment")
    epi = sp_.groupby(["resection", "compartment", "stratum"], as_index=False).agg(
        **{g: (g, "sum") for g in GENES}, depth=("depth", "sum"), n_cells=("depth", "size"))
    epi = epi[epi.n_cells >= 20]
    tab = pd.concat([host, epi], ignore_index=True).rename(columns={"resection": "donor"})
    tab["donor"] = tab.donor.astype(str)
    anc = [a for a in anchors(name, "learned") if family(a) != "epithelial" and a in set(host.compartment)]
    if len(anc) < 3:
        return pd.DataFrame([dict(study=name, status=f"{len(anc)} anchors (Dataset.fit needs 3)")])
    ds = cs.Dataset(tab, anchors=anc, name=study)
    bases = [b for b in ("tumour:epithelial", "adjacent_normal:epithelial") if f"{b}|G1" in set(ds.compartments)]
    comps = [c for c in ds.compartments if "|" in c]
    pt = statistics(ds, bases, comps)
    rng = np.random.default_rng(seed)
    reps = []
    for _ in range(n_boot):
        try:
            reps.append(statistics(ds.resample(rng), bases, comps))
        except (ValueError, np.linalg.LinAlgError):
            continue

    def interval(point, draws):
        x = np.asarray(draws, float)
        x = x[np.isfinite(x)]
        x = point + wrap(x - point) if np.isfinite(point) else np.array([])
        ok = len(x) > 20
        return (np.std(x, ddof=1) if ok else np.nan, np.quantile(x, 0.025) if ok else np.nan, np.quantile(x, 0.975) if ok else np.nan)

    rows = []
    for key, (u, lu) in pt.items():
        se, lo, hi = interval(u, [r[key][0] for r in reps])
        a, g = f"{key[0]}|{key[1]}", f"{key[0]}|G1"
        n = len(set(ds.table.loc[ds.table.compartment == a, "resection"]) & set(ds.table.loc[ds.table.compartment == g, "resection"])) \
            if a in set(ds.compartments) else 0
        rows.append(dict(study=name, status="ok" if n >= 6 else "< 6 patients", compartment=key[0], part=key[1], y5=u if n >= 6 else np.nan,
                         log_amp_ratio=lu, se=se, lo=lo, hi=hi, n_patients=n))
    for p in PARTS:
        kt, kn = ("tumour:epithelial", p), ("adjacent_normal:epithelial", p)
        if kt in pt and kn in pt and np.isfinite(pt[kt][0]) and np.isfinite(pt[kn][0]):
            v = float(wrap(pt[kt][0] - pt[kn][0]))
            x = [wrap(r[kt][0] - r[kn][0]) for r in reps]
            se, lo, hi = interval(v, x)
            rows.append(dict(study=name, status="ok", compartment="tumour_minus_normal", part=p, y5=v, se=se, lo=lo, hi=hi))
    return pd.DataFrame(rows)


def _has(study: str) -> bool:
    cc = data._builder("cells_cancer")
    return (cc.OUT / f"{study}_epithelial.csv.gz").exists() or (cc.LEGACY / "states" / "cells" / f"{study}.csv.gz").exists()


def pool(t: pd.DataFrame) -> dict:
    t = t[(t.se <= SE_MAX) & t.y5.notna()]
    if len(t) < 2:
        return dict(k=len(t), y5=np.nan, se=np.nan, lo=np.nan, hi=np.nan, tau=np.nan)
    m = cs.stats.meta(t.y5, t.se)
    return dict(k=m["k"], y5=m["random"], se=m["se_random"], lo=m["lo"], hi=m["hi"], tau=float(np.sqrt(m["tau2"])))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--definition", default="seurat", choices=("seurat", "quartile"))
    ap.add_argument("--build-missing", action="store_true", help="build the per-cell tables of studies without legacy ones")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = 10 if run["quick"] else 500
    reg = pd.read_csv(TABLE)
    studies = sorted(t.removeprefix("cancer_") for t in reg.target.unique() if t.startswith("cancer_"))
    if run["quick"]:
        studies = [s for s in studies if s in ("SRP273482", "GSE132465", "SRP327758")]
    res = Parallel(n_jobs=args.jobs)(delayed(study_arm)(s, 9500 + i, args.definition, n_boot, args.subsample, args.seed,
                                                       args.build_missing) for i, s in enumerate(studies))
    t = pd.concat(res, ignore_index=True)
    tag = "" if args.definition == "seurat" else f"_{args.definition}"
    t.to_csv(out / f"a24_cancer{tag}.csv", index=False)
    ok = t[t.status == "ok"]
    pooled = pd.DataFrame([dict(compartment=c, part=p, **pool(g)) for (c, p), g in ok.groupby(["compartment", "part"])])
    pooled.to_csv(out / f"a24_cancer_pooled{tag}.csv", index=False)
    print(t[t.status != "ok"][["study", "status"]].drop_duplicates().to_string(index=False))
    print(pooled.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
