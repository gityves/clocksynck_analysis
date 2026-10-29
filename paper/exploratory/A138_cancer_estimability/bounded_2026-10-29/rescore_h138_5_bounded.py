"""H138.5 re-scored with the bounded statistic (clocksync `stats.ccorr`, `nan_on_flag`). DESCRIPTIVE / POST HOC (2026-10-29,
roadmap 5.5): H138.5 was descriptive and out of the Holm family; nothing registered changes.

Same cohorts, units (host-referenced residual states, learned anchors as currently stored), null streams as `score.py`
(`one_study_obs`). A unit whose corrected variance (of either residual state) is <= 2 SE is NaN and drops out of the
pooled mean. The old statistic is recomputed on the same live units and compared with the stored `per_unit.csv` (the
anchors were relearned in Phase 3, so a mismatch is reported, not hidden). Primary (malignant vs normal epithelium) and
fibroblast companion.
    python rescore_h138_5_bounded.py [--jobs J] [--n-perm 10000]
Outputs: output/{per_unit,tests,breakdown}_bounded.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "scripts"))
sys.path.insert(0, str(HERE.parents[1] / "A149_paired_tumour_normal/bounded_2026-10-29"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import bounded_lib as B  # noqa: E402
import q3_lib as q  # noqa: E402
from common import cli  # noqa: E402
from power import SEED  # noqa: E402

OUT = HERE / "output"
OLD = HERE.parents[0] / "output" / "score"


def one(study, n_perm, i, seed):
    from clocksync.config import limit_threads
    limit_threads(1)
    pb, mal = q.study_table(study)
    out = {}
    for part, x in q.sc_units(study, pb, mal).items():
        if x is None:
            continue
        j = 0 if part == "primary" else 1
        rng = lambda: np.random.default_rng([SEED, 5, seed, i, j])  # noqa: E731
        raw = q.perm_null(x["za"], x["vA"], x["zb"], x["vB"], n_perm, rng())
        s = B.stat(x["za"], x["vA"], x["zb"], x["vB"])
        nl = B.perm_null(x, n_perm, rng())
        out[part] = dict(n=x["n"], null=nl, old_est=q.estimable(x["za"], x["vA"], x["zb"], x["vB"], raw), est=B.estimable(x, nl),
                         old=s["old"], new=s["new"], flag=s["flag"], malignant=mal, old_null_mean=float(np.nanmean(raw)))
    return study, out


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=10_000)
    args = ap.parse_args(argv)
    cli.setup(args)
    studies = list(q.sh.sc_builder().STUDIES)
    out = Parallel(n_jobs=args.jobs)(delayed(one)(s, args.n_perm, i, args.seed) for i, s in enumerate(studies))
    stored = pd.read_csv(OLD / "per_unit.csv")
    rows, keep = [], {"primary": [], "fibroblast": []}
    for study, res in out:
        for part, r in res.items():
            keep[part].append((study, r))
            h = "H138.5" if part == "primary" else "H138.5 fibroblast"
            s = stored[(stored.hypothesis == h) & (stored.unit == study)]
            rows.append(dict(part=part, unit=study, n=r["n"], stored_statistic=s.statistic.iloc[0] if len(s) else np.nan,
                             live_old_statistic=r["old"], old_gt1=r["old"] > 1, flagged=r["flag"], bounded_statistic=r["new"],
                             estimable_registered_rule=r["old_est"], stored_estimable=bool(s.estimable.iloc[0]) if len(s) else np.nan,
                             estimable_bounded=r["est"], p_bounded=B.perm_p(r["null"], r["new"]) if np.isfinite(r["new"]) else np.nan))
    tests = []
    st_old = pd.read_csv(OLD / "tests.csv")
    for part, name in (("primary", "H138.5"), ("fibroblast", "H138.5 fibroblast companion")):
        lst = [x for x in keep[part] if x[1]["est"]]
        if not lst:
            continue
        w = [r["n"] for _, r in lst]
        obs = float(q.pool(np.array([[r["new"]] for _, r in lst]), w)[0])
        nul = q.pool(np.array([r["null"] for _, r in lst]), w)
        o = st_old[st_old.hypothesis == name].iloc[0]
        tests.append(dict(hypothesis=name, n_units=sum(w), n_cohorts=len(lst), statistic_bounded=obs, p_bounded=B.perm_p(nul, obs),
                          null_mean=float(nul.mean()), null_sd=float(nul.std(ddof=1)), null_q95=float(np.quantile(nul, .95)),
                          mde=float(np.quantile(nul, .95) + .84 * nul.std(ddof=1)), cohorts=";".join(s for s, _ in lst),
                          excluded=";".join(s for s, r in keep[part] if not r["est"]),
                          old_n_units=o.n_units, old_n_cohorts=o.n_cohorts, old_statistic=o.statistic, old_p=o.p, old_null_mean=o.null_mean,
                          old_mde=o.detection_limit_mde, old_excluded=o.excluded_not_estimable))
    PER, T = pd.DataFrame(rows), pd.DataFrame(tests)
    PER.to_csv(OUT / "per_unit_bounded.csv", index=False)
    T.to_csv(OUT / "tests_bounded.csv", index=False)
    summ = PER.groupby("part").agg(units=("unit", "size"), flagged=("flagged", "sum"), old_gt1=("old_gt1", "sum"),
                                   estimable_registered=("estimable_registered_rule", "sum"), estimable_bounded=("estimable_bounded", "sum"))
    summ.to_csv(OUT / "breakdown_bounded.csv")
    print(summ, flush=True); print(PER.round(4).to_string(), flush=True); print(T.round(4).T.to_string(), flush=True)


if __name__ == "__main__":
    main()
