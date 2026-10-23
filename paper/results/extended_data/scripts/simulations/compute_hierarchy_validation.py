"""V16 and V18: validation of the A20 / A23 hierarchical model at the real structure, before the real fit.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A20: V16; A23: V18); ports of existing pre-registrations, not new tests.
The structure of the pooled tables is kept (rows = contrast x compartment, disease classes, families, cohorts, standard errors
s_i); true effects theta_i are simulated and y_i = theta_i + N(0, s_i^2).
  V16 (A20 inputs, estimands y1, y2, y3; scales 0.2, 0.5, 0.5):
    null: theta = 0, `--n-null` data sets per estimand. Criteria: share of rows called (lfsr < 0.05) <= 0.05; share of
      (class, family) cells whose 95 % interval covers 0 >= 0.93.
    planted (y1; `--n-planted`): mu = -0.20 for fibrosis x lymphocyte, +0.15 for fibrosis x endothelial, 0 elsewhere, plus
      N(0, 0.05^2) per row. Reported: posterior medians of the planted cells, coverage, share of planted rows called.
  V18 (A23 inputs, y4 in hours, rows with SE <= 4 h, scale 2 h): the same, null (100) and planted +3 h for fibrosis x epithelial
    and inflammation x epithelial (50). Criteria as V16 null.
Fits: 2 chains x 1000 draws (500 warm-up), `spread_multiplicity/_hierarchy.fit` (local: not in clocksync, imported, not copied).
Estimand: calling rate, interval coverage and recovery of the model; no null hypothesis test.

    python compute_hierarchy_validation.py [--subsample N] [--jobs 4] [--n-null 100] [--n-planted 50]
`--subsample N`: N data sets per scenario (null and planted), 1 chain pair as in the full run.
Inputs: the a20_inputs.csv and a23_inputs.csv of `spread_multiplicity` (`--inputs-dir`); if that folder has no full
run, the frozen tables of legacy/desync/common_output (the tables the legacy V16 / V18 used; identical structure).
Outputs: paper/results/extended_data/output/simulations/v16_{runs,summary}.csv, v18_{runs,summary}.csv
Legacy: legacy/scripts/v16_hierarchical.py (outputs/v16_hierarchical/), legacy/desync/common/a23_phase.py (v18; desync/common_output/a23_v18.csv)
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from common import cli  # noqa: E402
from common.paths import LEGACY, out_dir  # noqa: E402

FOLDER = "simulations"
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir, group_out  # noqa: E402
S06 = group_out("spread_multiplicity")
FROZEN = LEGACY / "desync" / "common_output"
SCALE = {"y1": 0.2, "y2": 0.5, "y3": 0.5, "y4": 2.0}
PLANT16 = {"fibrosis|lymphocyte": -0.20, "fibrosis|endothelial": 0.15}
PLANT18 = {"fibrosis|epithelial": 3.0, "inflammation|epithelial": 3.0}
SE_MAX = 4.0


def _hier():
    sys.path.insert(0, str(group_dir("spread_multiplicity")))
    import _hierarchy
    return _hierarchy


def read_input(name: str, inputs_dir) -> pd.DataFrame:
    for d in ([Path(inputs_dir)] if inputs_dir else [S06, FROZEN]):
        if (d / name).exists():
            return pd.read_csv(d / name)
    raise FileNotFoundError(name)


def one16(df, est, scenario, rep):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    se = "se" + est[1]
    df = df[np.isfinite(df[est]) & np.isfinite(df[se])].reset_index(drop=True)
    rng = np.random.default_rng([rep, 1 if scenario == "planted" else 0, int(est[1])])
    key = df.disease_class + "|" + df.family
    theta = np.zeros(len(df))
    if scenario == "planted":
        theta = key.map(PLANT16).fillna(0.0).to_numpy() + rng.normal(0, 0.05, len(df))
    sim = df.assign(**{est: theta + rng.normal(0, df[se].to_numpy())})
    r = _hier().fit(sim, est, se, SCALE[est], chains=2, warmup=500, samples=1000, seed=rep)
    cells, rows = r["cells"], r["rows"]
    out = dict(estimand=est, scenario=scenario, rep=rep, divergences=r["diag"]["divergences"],
               called=float((rows.lfsr < 0.05).mean()), cells_cover0=float(((cells.lo <= 0) & (cells.hi >= 0)).mean()))
    if scenario == "planted":
        for k, v in PLANT16.items():
            c = cells.set_index("cell").loc[k]
            m = (key == k).to_numpy()
            out |= {f"{k}_median": c["median"], f"{k}_covers": bool(c.lo <= v <= c.hi), f"{k}_called": float((rows.lfsr[m] < 0.05).mean())}
        out["called_unplanted"] = float((rows.lfsr[~key.isin(PLANT16).to_numpy()] < 0.05).mean())
    return out


def one18(df, scenario, rep):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    rng = np.random.default_rng([rep, scenario == "planted"])
    key = df.disease_class + "|" + df.family
    theta = key.map(PLANT18).fillna(0.0).to_numpy() if scenario == "planted" else np.zeros(len(df))
    sim = df.assign(y4=theta + rng.normal(0, df.se4.to_numpy()))
    r = _hier().fit(sim, "y4", "se4", SCALE["y4"], chains=2, warmup=500, samples=1000, seed=rep)
    c, rows = r["cells"].set_index("cell"), r["rows"]
    out = dict(scenario=scenario, rep=rep, called=float((rows.lfsr < 0.05).mean()),
               cells_cover0=float(((c.lo <= 0) & (c.hi >= 0)).mean()), divergences=r["diag"]["divergences"])
    for k, v in PLANT18.items():
        if k in c.index:
            out[f"{k}_median"], out[f"{k}_covers"] = c.loc[k, "median"], bool(c.loc[k, "lo"] <= v <= c.loc[k, "hi"])
    return out


def summarise(runs, by):
    s = runs.groupby(by).mean(numeric_only=True).drop(columns="rep")
    s["n_runs"] = runs.groupby(by).size()
    return s


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-null", type=int, default=100)
    ap.add_argument("--n-planted", type=int, default=50)
    ap.add_argument("--inputs-dir", default=None)
    ap.add_argument("--only", nargs="*", default=["v16", "v18"])
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_null, n_pl = (args.subsample, args.subsample) if run["quick"] else (args.n_null, args.n_planted)
    out = out_dir(FOLDER, run["quick"])
    if "v16" in args.only:
        df = read_input("a20_inputs.csv", args.inputs_dir)
        jobs = [(e, "null", i) for e in ("y1", "y2", "y3") for i in range(n_null)] + [("y1", "planted", i) for i in range(n_pl)]
        runs = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(one16)(df, *j) for j in jobs))
        runs.to_csv(out / "v16_runs.csv", index=False)
        s = summarise(runs, ["estimand", "scenario"])
        s.to_csv(out / "v16_summary.csv")
        print(s.round(3).T.to_string(), flush=True)
    if "v18" in args.only:
        t = read_input("a23_inputs.csv", args.inputs_dir)
        df = t[(t.se4 <= SE_MAX) & t.y4.notna()].reset_index(drop=True)
        jobs = [("null", i) for i in range(n_null)] + [("planted", i) for i in range(n_pl)]
        runs = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(one18)(df, *j) for j in jobs))
        runs.to_csv(out / "v18_runs.csv", index=False)
        s = summarise(runs, ["scenario"])
        s.to_csv(out / "v18_summary.csv")
        print(s.round(3).T.to_string(), flush=True)


if __name__ == "__main__":
    main()
