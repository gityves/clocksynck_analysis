"""A144 part 2 (D144.2, descriptive): does coupling follow within-sample synchrony R or per-cell amplitude lambda?

Registered in paper/docs/registered_experiments.md (A144). Rows: eligible rows of A18's `a18_cell_split.csv` in control groups
(and `lung_hlca` "all"), `age_mssm` dropped (MSSM controls counted once, via AD_mssm). c is the coupling of the same group's
anchored fit (host refitted on the group's donors, frozen learned A17 anchors; the fit of `compute_per_cell.coherent_s`, here
also keeping `coupling`). Consistency: rho = sqrt(s / (1 - R_W^2)), s = z2 V, against the table's rho (stop above 1e-6).
Models (OLS, target fixed effects): logit(c) ~ log lambda + log R, and logit(c) ~ log rho; 1,000 bootstraps of compartments
within target. Rows with c <= 0, c >= 1 or non-finite are dropped. No pass / fail. Compartment-level aggregates only.

    python compute_cell_decomposition.py
Outputs: paper/exploratory/A144_strong_clocks/output/{d144_2_rows,d144_2_coefficients}.csv
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

BB = PAPER / "results/figure_2/scripts/backbone"
_spec = importlib.util.spec_from_file_location("a18_per_cell", BB / "compute_per_cell.py")
PC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(PC)

AMENDMENT, SEED, N_BOOT = "A144_strong_clocks", 144, 1000


def fit_group(pb, anchors, k, donors):
    """(s, c) of compartment k over `donors`: as `PC.coherent_s` (host refitted, duplicates relabelled), plus coupling."""
    t = pd.concat([pb[pb.donor == r].assign(donor=f"{r}#{j}") for j, r in enumerate(donors)], ignore_index=True)
    anc = [a for a in anchors if a in set(t.compartment)]
    if len(anc) < 3 or k not in set(t.compartment) or t[t.compartment == k].donor.nunique() < 5:
        return np.nan, np.nan
    r = cs.Dataset(t, anchors=anc, name="cells").fit(n_boot=0, gate=False, targets=[k]).table.loc[k]
    return float(r.z2 * r.V), float(r.coupling)


def build_rows() -> tuple[pd.DataFrame, float]:
    a18 = pd.read_csv(BB.parents[1] / "output/backbone/a18_cell_split.csv")
    sel = a18[a18.eligible & ((a18.group == "control") | ((a18.target == "lung_hlca") & (a18.group == "all")))
              & (a18.target != "age_mssm")]
    out, anchors, tg = [], PC.frozen(), {}
    for _, r in sel.iterrows():
        if r.target not in tg:
            T = PC.Target(r.target)
            a = anchors[anchors.target == r.target]
            lab = T.labels if T.labels is not None else pd.Series("all", index=T.pb.donor.unique())
            tg[r.target] = (T.pb[T.pb.donor.isin(lab.index)], a.compartment[a.anchor.astype(bool)].tolist(), lab)
        pb, anc, lab = tg[r.target]
        s, c = fit_group(pb, anc, r.compartment, lab.index[lab == r.group].to_numpy())
        rho = np.sqrt(s / (1 - r.rw ** 2)) if np.isfinite(s) and s > 0 else np.nan
        out.append(dict(target=r.target, compartment=r.compartment, group=r.group, lam=r.lam, R=r.R, rho=rho, c=c, D=r.D,
                        rho_a18=r.rho))
    df = pd.DataFrame(out)
    return df, float(np.nanmax(np.abs(df.rho - df.rho_a18) / df.rho_a18))


def design(df, cols):
    T = pd.get_dummies(df.target).astype(float)
    return np.column_stack([df[cols].to_numpy(float), T.to_numpy()]), T.shape[1]


def ols(df, cols):
    X, _ = design(df, cols)
    return np.linalg.lstsq(X, df.y.to_numpy(float), rcond=None)[0][:len(cols)]


def coefs(df):
    b = ols(df, ["llam", "lR"])
    return np.array([b[0], b[1], b[1] - b[0], ols(df, ["lrho"])[0]])


def main():
    limit_threads(2)
    warnings.filterwarnings("ignore")
    rows, maxdiff = build_rows()
    print(f"rho reproduction: max relative difference {maxdiff:.3e}", flush=True)
    if not maxdiff < 1e-6:
        sys.exit("rho does not reproduce the A18 table: stopping")
    ok = np.isfinite(rows.c) & (rows.c > 0) & (rows.c < 1) & np.isfinite(rows.lam) & (rows.lam > 0) & np.isfinite(rows.R) & (rows.R > 0)
    print("rows per target (kept / dropped):\n" + rows.assign(kept=ok).groupby("target").kept.agg(["sum", "size"]).to_string())
    print("dropped: c <= 0", int((rows.c <= 0).sum()), "c >= 1", int((rows.c >= 1).sum()), "non-finite", int((~np.isfinite(rows.c)).sum()))
    d = rows[ok].copy().reset_index(drop=True)
    d["y"] = np.log(d.c / (1 - d.c))
    d["llam"], d["lR"], d["lrho"] = np.log(d.lam), np.log(d.R), np.log(d.rho)
    exploratory = exploratory_dir(AMENDMENT)
    rows.drop(columns=["rho_a18", "group"]).assign(kept=ok).to_csv(exploratory / "d144_2_rows.csv", index=False)

    est = coefs(d)
    rng = np.random.default_rng(SEED)
    idx = {t: g.index.to_numpy() for t, g in d.groupby("target")}
    B = np.array([coefs(d.loc[np.concatenate([rng.choice(v, len(v)) for v in idx.values()])].reset_index(drop=True))
                  for _ in range(N_BOOT)])
    names = ["beta_lambda", "beta_R", "beta_R_minus_beta_lambda", "beta_rho"]
    res = pd.DataFrame(dict(coefficient=names, estimate=est, lo=np.quantile(B, .025, axis=0), hi=np.quantile(B, .975, axis=0),
                            n_rows=len(d), n_targets=len(idx), n_boot=N_BOOT))
    res.to_csv(exploratory / "d144_2_coefficients.csv", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
