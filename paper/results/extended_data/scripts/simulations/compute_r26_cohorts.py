"""R26: the A25 cohorts with the non-temporal host direction removed (A26).

Pre-registered as A26 (docs/PLAN_internal_desync.md (archive), 2026-10-14; run note: V26.1 fails, V26.3 fails with g_GTEx, V26.2 passes, A26 rejected,
so per-cohort centres stay uninterpreted). R26: the A25 cohorts and classes (`windows/compute_internal_time.py`: living,
deceased, other, cancer cohorts with an anchored host), corrected with g_GTEx (primary; g_Chen sensitivity), then the A25 estimator
(window axis / curvature centre, r, flip support, donor bootstrap; interpretable when the upper bound of r < 0.6).
  g learned on GTEx fast deaths (42 tissues) as in `compute_nontemporal.py` (`learn_g`); donor score s_d = mean over a donor's compartments of
  g_perp . (I - P) y / |g_perp|^2; the correction removes the regression of the clock coefficients on s_d per intercept group, here by
  binomial thinning of the counts (`compute_nontemporal.correct_counts`; A26 deviation of figure_S02). No recorded time is used on the data
  being corrected.
  T1-T3 as A25 (`compute_internal_time.tests`); T4 (descriptive): concentration of the axis centres (mod 12 h) of the deceased cohorts, A25
  (uncorrected, same run) against corrected, expected to fall.
Reading (A26): if V26.1 and V26.2 pass, per-cohort centres are interpreted for interpretable cohorts with flip support >= 0.8; V26.1 fails (and
V26.1 is not ported: `Dataset` takes counts, not Chen's arrays), so A25's cohort-level reading stands.

    python compute_r26_cohorts.py [--subsample N] [--jobs 4] [--g gtex chen]
`--subsample N`: eight cohorts, N donors each, 20 bootstraps.
Outputs: paper/results/extended_data/output/simulations/r26_{cohorts,tests_<g>,t4}.csv
Legacy: legacy/desync/common/a26_nontemporal.py (apply step; common_output/a26_cohorts.csv, a26_tests_{gtex,chen}.csv, a26_t4.csv)
"""

from __future__ import annotations

import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("windows")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import log_expr  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _internal_time as T  # noqa: E402  (windows)
import _shared as sh  # noqa: E402  (windows)
import compute_internal_time as CI  # noqa: E402
import compute_nontemporal as NT  # noqa: E402

FOLDER = "simulations"
GENES = NT.GENES


def learn(source: str, subsample, seed) -> pd.Series:
    """The gene vector g (14 clock genes) of the non-temporal state, learned on timed data."""
    if source == "gtex":
        pb, donors = sh.gtex_fast(30, subsample, seed)
        y = NT.gene_centred(log_expr(pb[GENES], pb.depth.to_numpy(float)).to_numpy())
        return NT.learn_g(y, GENES, pb.donor.to_numpy(), pb.compartment.to_numpy(), donors.tod.reindex(pb.donor).to_numpy())[0]
    s, y, genes, _ = NT.chen_inputs()
    return NT.learn_g(y, genes, s.donor.to_numpy(), s.region.astype(str).to_numpy(), (s.tod_raw.astype(float) + NT.CHEN_SUNRISE).to_numpy())[0]


def read(ds, e, n_boot, seed):
    reps = T.boot(ds, ds.anchors, e, n_boot, seed)
    return T.summarise(e, reps)


def one(entry, klass, g, src, n_boot, subsample, seed):
    limit_threads(1)
    try:
        pb, donors, anc, rule = CI.cohort(entry, subsample, seed)
        if anc is None:
            return [dict(entry=entry, klass=klass, g="none", status="no anchors"), dict(entry=entry, klass=klass, g=src, status="no anchors")]
        ds = cs.Dataset(pb, anchors=anc, name=entry)
        e0, _ = T.host_estimate(ds, anc)
        rows = [dict(entry=entry, klass=klass, g="none", rule=rule, n=len(ds.donor_ids), n_anchors=len(anc), status="ok") | read(ds, e0, n_boot, seed)]
        y = NT.gene_centred(log_expr(ds.table[GENES], ds.table.depth.to_numpy(float)).to_numpy())
        s = NT.donor_score(y, ds.data.group, ds.data.resection, ds.waveform, g)
        dsc = cs.Dataset(NT.correct_counts(ds, s, seed), anchors=anc, name=entry)
        e1, _ = T.host_estimate(dsc, anc)
        rows.append(dict(entry=entry, klass=klass, g=src, rule=rule, n=len(dsc.donor_ids), n_anchors=len(anc), status="ok") | read(dsc, e1, n_boot, seed))
        return rows
    except Exception as ex:
        return [dict(entry=entry, klass=klass, g=src, status=f"error: {type(ex).__name__}: {ex}")]


def concentration(c12) -> float:
    """Mean resultant length of axis centres (mod 12 h)."""
    c = np.asarray(c12, float)
    return float(abs(np.exp(1j * 2 * np.pi * c / 12).mean())) if len(c) else np.nan


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--g", nargs="+", default=["gtex", "chen"], choices=["gtex", "chen"])
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    cancer = sorted(t for t in pd.read_csv(sh.AN.TABLE).target.unique() if str(t).startswith("cancer_"))
    todo = [(e, "living") for e in CI.LIVING + cancer] + [(e, "deceased") for e in CI.DECEASED] + [(e, "other") for e in CI.OTHER]
    if run["quick"]:
        todo = [t for t in todo if t[0] in ("skin_ad", "kidney_sc", "heart", "lung_natri_control", "islets", "ckd_sc", "psychad_mssm", "lung_hlca|hospital")]
    sub = args.subsample
    gs = {src: learn(src, sub, args.seed) for src in args.g}
    allrows = []
    for src, g in gs.items():
        res = Parallel(n_jobs=args.jobs)(delayed(one)(e, k, g, src, run["n_boot"], sub, args.seed) for e, k in todo)
        rows = pd.DataFrame([r for rr in res for r in rr])
        allrows.append(rows)
    D = pd.concat(allrows, ignore_index=True)
    D = D[~((D.g == "none") & D.duplicated(["entry", "g"]))].reset_index(drop=True)  # the uncorrected rows are the same for every g
    D.to_csv(out / "r26_cohorts.csv", index=False)
    t4 = []
    for src in args.g:
        o = pd.concat([D[D.g == "none"], D[D.g == src]], ignore_index=True)
        TT = CI.tests(o[o.g == src].reset_index(drop=True))
        TT.to_csv(out / f"r26_tests_{src}.csv")
        print(f"--- g = {src}\n{TT.to_string()}")
        dec = lambda g_: D[(D.g == g_) & (D.status == "ok") & (D.klass == "deceased")]  # noqa: E731
        t4.append(dict(g=src, T4_R_deceased_before=concentration(dec("none").c_axis12), T4_R_deceased_after=concentration(dec(src).c_axis12)))
    t4 = pd.DataFrame(t4)
    t4.to_csv(out / "r26_t4.csv", index=False)
    print(t4.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
