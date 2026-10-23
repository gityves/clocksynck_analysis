"""A20 and A23: hierarchical pooling of the disease contrasts (disease class x cell family), readouts R20.1-R20.4, R23.1,
and their sensitivities; the multiplicity of the per-contrast tests.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A20, A23). Reads the inputs of `compute_spread_inputs.py` (same folder:
a20_inputs*.csv, a23_inputs.csv). Model (`_hierarchy.fit`, NumPyro NUTS): y_ck ~ N(theta_ck, s_ck^2), theta_ck =
mu_{class(c), family(k)} + a_cohort + b_contrast + e_ck, half-normal priors on every tau (scale 0.2 for y1, 0.5 for
y2 and y3, 2 h for y4; S = 2 scale), 4 chains x 2,000 draws (1,000 warm-up).
  R20.1 fibrosis x lymphocyte (y1, y2, y3); R20.2 fibrosis x fibroblast_stromal / endothelial / mural (y1);
  R20.3 neurodegeneration, every family: is the 95 % interval of mu (y1) within +-0.10 (resilient)?;
  R20.4 every class x family with a 95 % interval excluding 0 (y1), with the y2 and y3 cells beside it;
  compartment calls: lfsr < 0.05, beside the BH-adjusted per-contrast permutation p* (multiplicity over all tests);
  R23.1 every class x family, pooled y4 (rows with SE <= 4 h).
Sensitivities (y1, the R20 cells): prior scale x0.5 and x2, Student-t likelihood, leave one contrast out, fixed-rule
anchors, SEs from the V15 null simulations (not ported).

    python compute_hierarchy.py [--subsample N] [--jobs 4] [--main-only]
Outputs (paper/results/figure_4/output/spread_multiplicity/): a20_{cells,rows,diag,readouts,sensitivity}.csv,
a23_{cells,rows,diag}.csv, multiplicity.csv
Legacy: legacy/desync/common/a20_pool.py, a23_phase.py (pool step). V16 / V18 simulations are not ported.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "spread_multiplicity"
SCALE = {"y1": 0.2, "y2": 0.5, "y3": 0.5}
SCALE_Y4 = 2.0
R20 = ["fibrosis|lymphocyte", "fibrosis|fibroblast_stromal", "fibrosis|endothelial", "fibrosis|mural"]
EQUIV = 0.10


def _fit(df, est, se, scale, tag, **kw):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _hierarchy as H
    return tag, est, H.fit(df, est, se, scale, **kw)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--main-only", action="store_true", help="main fits only, no sensitivities")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    df = pd.read_csv(out / "a20_inputs.csv")
    quick = run["quick"]
    main_kw = dict(chains=2, warmup=200, samples=300) if quick else dict(chains=4, warmup=1000, samples=2000)
    sens = dict(chains=2, warmup=100, samples=200) if quick else dict(chains=2, warmup=500, samples=1000)
    jobs = [(df, e, "se" + e[1], SCALE[e], "main", main_kw) for e in ("y1", "y2", "y3")]
    if not args.main_only:
        jobs += [(df, "y1", "se1", 0.1, "scale_x0.5", sens), (df, "y1", "se1", 0.4, "scale_x2", sens),
                 (df, "y1", "se1", 0.2, "student_t", dict(sens, robust=True))]
        jobs += [(df[df.contrast != c], "y1", "se1", 0.2, f"drop_{c}", sens) for c in sorted(df.contrast.unique())]
        fx = out / "a20_inputs_fixed.csv"
        if fx.exists():
            jobs.append((pd.read_csv(fx), "y1", "se1", 0.2, "fixed_rule", sens))
    res = Parallel(n_jobs=args.jobs)(delayed(_fit)(d, e, s, sc, t, **kw) for d, e, s, sc, t, kw in jobs)
    cells, rows, diag = [], [], []
    for tag, est, r in res:
        cells.append(r["cells"].assign(fit=tag, estimand=est))
        diag.append(dict(fit=tag, estimand=est, **r["diag"], **r["tau"]))
        if tag == "main":
            rows.append(r["rows"].assign(estimand=est))
    cells, diag, rows = pd.concat(cells, ignore_index=True), pd.DataFrame(diag), pd.concat(rows, ignore_index=True)
    # multiplicity: BH over every (contrast, compartment) permutation p* of the pooled contrasts, beside the lfsr calls
    key = ["contrast", "compartment"]
    p = df[key + ["p_star"]].drop_duplicates(key)
    p["q_bh_perm"] = cs.stats.bh(p.p_star.fillna(1.0).to_numpy())
    rows = rows.merge(p[key + ["q_bh_perm"]], on=key, how="left")
    cells.to_csv(out / "a20_cells.csv", index=False)
    rows.to_csv(out / "a20_rows.csv", index=False)
    diag.to_csv(out / "a20_diag.csv", index=False)
    main = cells[cells.fit == "main"]
    y1 = main[main.estimand == "y1"].set_index("cell")
    ro = []
    for c in R20:
        if c in y1.index:
            row = dict(readout="R20.1" if c.endswith("lymphocyte") else "R20.2", cell=c)
            for e in ("y1", "y2", "y3"):
                m = main[(main.estimand == e) & (main.cell == c)]
                if len(m):
                    row.update({f"{e}_median": m["median"].iloc[0], f"{e}_lo": m.lo.iloc[0], f"{e}_hi": m.hi.iloc[0],
                                f"{e}_p_pos": m.p_pos.iloc[0]})
            ro.append(row)
    for c, r in y1[y1.disease_class == "neurodegeneration"].iterrows():
        ro.append(dict(readout="R20.3", cell=c, y1_median=r["median"], y1_lo=r.lo, y1_hi=r.hi, y1_p_pos=r.p_pos,
                       resilient=bool(r.lo > -EQUIV and r.hi < EQUIV)))
    for c, r in y1[(y1.lo > 0) | (y1.hi < 0)].iterrows():
        row = dict(readout="R20.4", cell=c, y1_median=r["median"], y1_lo=r.lo, y1_hi=r.hi)
        for e in ("y2", "y3"):
            m = main[(main.estimand == e) & (main.cell == c)]
            if len(m):
                row.update({f"{e}_median": m["median"].iloc[0], f"{e}_lo": m.lo.iloc[0], f"{e}_hi": m.hi.iloc[0]})
        ro.append(row)
    pd.DataFrame(ro).to_csv(out / "a20_readouts.csv", index=False)
    cells[(cells.estimand == "y1") & cells.cell.isin(R20)].pivot_table(index="fit", columns="cell", values="median") \
        .to_csv(out / "a20_sensitivity.csv")
    m1 = rows[rows.estimand == "y1"]
    mult = dict(rows_y1=len(m1), lfsr_lt_005=int((m1.lfsr < 0.05).sum()), lfsr_specific_lt_005=int((m1.lfsr_specific < 0.05).sum()),
                p_star_lt_005=int((p.p_star < 0.05).sum()), tests=int(p.p_star.notna().sum()),
                expected_p_star_lt_005=0.05 * int(p.p_star.notna().sum()), bh_q_lt_005=int((p.q_bh_perm < 0.05).sum()))

    # ---- A23: relative phase -------------------------------------------------------------------------------------
    a = pd.read_csv(out / "a23_inputs.csv")
    a = a[a.y4.notna() & (a.se4 <= sh.SE_MAX_H)].reset_index(drop=True)
    kw4 = dict(chains=2, warmup=200, samples=300) if quick else dict(chains=4, warmup=1000, samples=2000)
    _, _, r4 = _fit(a, "y4", "se4", SCALE_Y4, "main", **kw4)
    r4["cells"].to_csv(out / "a23_cells.csv", index=False)
    r4["rows"].to_csv(out / "a23_rows.csv", index=False)
    pd.DataFrame([dict(**r4["diag"], **r4["tau"], rows=len(a))]).to_csv(out / "a23_diag.csv", index=False)
    mult.update(a23_rows=len(a), a23_lfsr_lt_005=int((r4["rows"].lfsr < 0.05).sum()),
                a23_cells_excluding_0=int(((r4["cells"].lo > 0) | (r4["cells"].hi < 0)).sum()), a23_cells=len(r4["cells"]))
    pd.Series(mult).to_csv(out / "multiplicity.csv", header=False)
    pd.set_option("display.width", 220)
    print(diag.round(3).to_string(index=False))
    print(pd.DataFrame(ro).round(3).to_string(index=False))
    print(pd.Series(mult).to_string())


if __name__ == "__main__":
    main()
