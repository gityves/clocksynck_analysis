"""CYCLOPS 2.0 (Hammarlund & Anafi; commit 56dbacb, Julia 1.6.7) run as a subprocess with its own depot and project
under `$SCCC_DATA/tools`. `CYCLOPS.Fit` fixes its RNG seed (1234), so `seed` shuffles the sample order instead.
Hyperparameters are those of `CYCLOPS_2_0_Template.jl` unless overridden.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pandas as pd

from common.paths import DATA

TOOLS = DATA / "tools"
JULIA = TOOLS / "julia-1.6.7" / "bin" / "julia"
JULIA_DEPOT = TOOLS / "julia_depot"
JULIA_PROJECT = TOOLS / "cyclops_env"
CYCLOPS_JL = TOOLS / "CYCLOPS-2.0" / "CYCLOPS.jl"
DRIVER = Path(__file__).with_name("cyclops_run.jl")


def _toml(d: dict) -> str:
    def fmt(v):
        if isinstance(v, bool):
            return str(v).lower()
        if isinstance(v, (list, tuple)):
            return "[" + ", ".join(fmt(x) for x in v) + "]"
        return json.dumps(v)
    return "".join(f"{k} = {fmt(v)}\n" for k, v in d.items())


def run_cyclops(expr: pd.DataFrame, seed_genes: list[str], out_dir: Path,
                times: pd.Series | None = None, seed: int = 0, threads: int = 8,
                covariates: pd.DataFrame | None = None, train_times: pd.Series | None = None, **hyper) -> pd.DataFrame:
    """Run CYCLOPS 2.0 on `expr` (genes x samples, non-log). Returns per-sample phases in radians.

    `times`: collection times in hours (index = sample IDs, any subset); used by CYCLOPS.Align to
    rotate/reflect the phases to clock time. Without it phases are aligned to mouse clock genes if present
    (otherwise arbitrary origin/direction). `covariates`: samples x discrete covariates (e.g. dataset), written as
    `<name>_D` rows; the first one is also the eigengene-regression covariate (template default). `train_times`:
    hours of time-stamped samples used during training (semi-supervised, `:train_sample_id`/`:train_sample_phase`).
    `hyper`: CYCLOPS hyperparameter overrides (template names, e.g. train_n_models=20, seed_min_CV=0.1).
    Returns columns phase_raw, phase, aligned_by (index = sample).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # CYCLOPS needs plain sample names (letters, digits, underscore): use S0..Sn and map back.
    ids = [f"S{i}" for i in range(expr.shape[1])]
    back = dict(zip(ids, expr.columns))
    x = expr.copy()
    x.columns = ids
    x.index.name = "Gene_Symbol"
    x = x.map(lambda v: f"{v:.6g}")
    if covariates is not None:  # covariate rows come first; values must start with a letter
        cov = covariates.loc[expr.columns].astype(str).apply(lambda c: "G" + c).T
        cov.index = [f"{c}_D" for c in covariates.columns]
        cov.columns = ids
        x = pd.concat([cov, x])
    x.rename_axis("Gene_Symbol").to_csv(out_dir / "expr.csv")
    (out_dir / "seeds.txt").write_text("\n".join(seed_genes) + "\n")
    cmd = [str(JULIA), "--project=" + str(JULIA_PROJECT), str(DRIVER), f"expr={out_dir / 'expr.csv'}",
           f"seeds={out_dir / 'seeds.txt'}", f"out={out_dir}", f"seed={seed}", f"workers={threads}",
           f"cyclops={CYCLOPS_JL}"]
    if hyper:
        (out_dir / "hyper.toml").write_text(_toml(hyper))
        cmd.append(f"hyper={out_dir / 'hyper.toml'}")
    if train_times is not None:
        tt = pd.Series(train_times).rename(index={v: k for k, v in back.items()})
        tt.rename_axis("sample").rename("time_hours").reset_index().to_csv(out_dir / "train_times.csv", index=False)
        cmd.append(f"train_times={out_dir / 'train_times.csv'}")
    if times is not None:
        t = pd.Series(times).rename(index={v: k for k, v in back.items()})
        t = t[t.index.isin(ids)]
        t.rename_axis("sample").rename("time_hours").reset_index().to_csv(out_dir / "times.csv", index=False)
        cmd.append(f"times={out_dir / 'times.csv'}")
    env = {**os.environ, "JULIA_DEPOT_PATH": str(JULIA_DEPOT), "JULIA_PROJECT": str(JULIA_PROJECT),
           "JULIA_NUM_THREADS": str(threads), "MPLBACKEND": "Agg", "OPENBLAS_NUM_THREADS": str(threads)}
    with open(out_dir / "cyclops.log", "w") as log:
        subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    res = pd.read_csv(out_dir / "phases.csv", index_col="sample")
    res.index = res.index.map(back)
    return res.loc[expr.columns]



def seeds() -> list[str]:
    """CYCLOPS's human seed genes (`cyclops_seed_genes.txt`, source on its first line)."""
    lines = Path(__file__).with_name("cyclops_seed_genes.txt").read_text().splitlines()
    return [g.strip() for g in lines if g.strip() and not g.startswith("#")]
