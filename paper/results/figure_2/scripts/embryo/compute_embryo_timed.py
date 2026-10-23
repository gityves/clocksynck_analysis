"""A108 / A108b: does the late mouse embryo's clock state follow the time of day? (Qiu 2024, E15.0-E18.75)

Pre-registered in docs/PLAN_internal_desync.md (archive) (A108: H108.1 and its sensitivity and controls; A108b: S1-S3, D1).
Estimand: R = |mean exp(i (phi_s - w t_s))|, the alignment of the untimed host phase phi_s of an embryo with its clock
time t_s = 24 x frac(E-day) up to one common offset (host on all cell types with >= 50 cells present in >= 75 % of the
window's 16 embryos, batch as stratum). Null: clock times permuted over embryos (10,000), one-sided (R >= 0); the
Fisher-Lee r beside it (the Jammalamadaka-SenGupta r is ill defined for four evenly spaced times). Supported if R > 0.5
and p < 0.05 (H108.1).
  H108.1               host phase vs clock time
  sensitivity          clock-plane states residualised on E-day (per cell type) before the host fit
  controls             3 pseudo-clock gene sets (expression-matched non-clock genes of the 67-gene panel)
  A108b S1             leave one E-day (batch of 4 embryos) out
  A108b S2 / S3        host without blood and liver lineages; CNS only; non-CNS non-blood only
  A108b D1             glucocorticoid-response score (FKBP5, TSC22D3, ZBTB16, DDIT4, KLF9): cosinor with an E-day
                       term, |a| permuted (descriptive)

    python compute_embryo_timed.py [--subsample N] [--seed S]
Outputs: paper/results/figure_2/output/embryo/a108_timed.csv (one row per test; add --subsample for checks)
Legacy:  legacy/desync/validation/natural_experiments/scripts/a108_embryo.py (H108.1 part) and a108b_embryo.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

FOLDER = "embryo"
GC = ["FKBP5", "TSC22D3", "ZBTB16", "DDIT4", "KLF9"]
KEEP = ["n", "R", "p_R", "r_fl", "p_fl", "offset_h", "median_abs_err_h"]


def timed(pb, donors, seed, n_perm, comps=None, residualise=False, waveform=None, name=""):
    """Host of the cell types `comps` (default all) and the agreement of its phase with clock time."""
    t = pb if comps is None else pb[pb.compartment.isin(comps)]
    kw = {} if waveform is None else dict(waveform=waveform)
    ds = cs.Dataset(t, anchors=sorted(t.compartment.unique()), donors=donors, name=name, **kw)
    if residualise:
        ds = S.residualised(ds, donors.day)
    r = ds.timed("t_clock", n_perm=n_perm, seed=seed)
    return {k: r[k] for k in KEEP} | dict(compartments=t.compartment.nunique())


def glucocorticoid(panel, don, donors, n_perm, seed):
    """D1: per-embryo score (mean of within-window z of log CPM, all cell types summed) and its cosinor on clock time
    with an E-day term. local score; the cosinor is `clocksync.timed.cosinor`."""
    w = panel[panel.donor.isin(don)].groupby("donor")[GC + ["depth"]].sum()
    z = np.log(w[GC].div(w.depth, axis=0) * 1e6 + 1)
    z = ((z - z.mean()) / z.std()).mean(axis=1)
    t, day = donors.t_clock.reindex(z.index).to_numpy(), donors.day.reindex(z.index).to_numpy()
    amp = cs.timed.cosinor(z.to_numpy(), t, covariates=day)
    rng = np.random.default_rng(seed)
    null = [cs.timed.cosinor(z.to_numpy(), rng.permutation(t), covariates=day)["amplitude"] for _ in range(n_perm)]
    return dict(n=len(z), amp_sd=float(amp["amplitude"]), peak_h=float(amp["acrophase_h"]),
                p_R=cs.stats.perm_p(null, amp["amplitude"], "greater", tie=0.0))


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = 500 if run["quick"] else 10_000
    pb, donors = data.load("mouse_qiu", subsample=args.subsample, seed=args.seed)
    win, don = S.window(pb, donors, *S.LATE)
    comps = sorted(win.compartment.unique())
    print(f"late window: {len(don)} embryos, {len(comps)} cell types", flush=True)
    rows = [dict(test="H108.1", **timed(win, donors, 1, n_perm)),
            dict(test="H108.1 sens: residualised on E-day", **timed(win, donors, 2, n_perm, residualise=True))]
    # pseudo-clock controls on the 67-gene panel
    panel = S.qiu_panel("celltype")
    panel = panel[panel.donor.isin(set(pb.donor))]
    pw, _ = S.window(panel, donors, *S.LATE)
    genes = [c for c in pw.columns if c not in ("donor", "compartment", "stratum", "n_cells", "depth")]
    lm = np.log(pw[genes].sum() / pw.depth.sum() + 1e-12)
    wf = cs.gtex()
    for j, s in enumerate(S.matched_sets(lm, wf.genes, 3, np.random.default_rng(7))):
        pseudo = cs.Waveform(tuple(s), wf.amp, wf.psi, None, None, f"pseudo{j}")
        rows.append(dict(test=f"H108.1 control: pseudo-clock set {j}", genes=";".join(s),
                         **timed(pw, donors, 10 + j, n_perm, waveform=pseudo)))
    # A108b
    batch = donors.batch.reindex(don)
    for j, b in enumerate(sorted(batch.unique())):
        keep = batch.index[batch != b]
        rows.append(dict(test=f"S1 without {b}", **timed(win[win.donor.isin(keep)], donors, 100 + j, n_perm,
                                                          comps=comps)))
    traj = S.trajectory_of().reindex(comps)
    sets = {"S2 without blood and liver": [c for c in comps if traj[c] not in S.BLOOD_LIVER],
            "S3 CNS only": [c for c in comps if traj[c] in S.CNS],
            "S3 periphery (non-CNS, non-blood)": [c for c in comps if traj[c] not in S.CNS | S.BLOOD_LIVER]}
    for j, (name, cc) in enumerate(sets.items()):
        try:
            rows.append(dict(test=name, **timed(win, donors, 200 + j, n_perm, comps=cc)))
        except ValueError as e:  # too few cell types in a subsample
            rows.append(dict(test=name, compartments=len(cc), status=str(e)[:80]))
    rows.append(dict(test="D1 glucocorticoid score", **glucocorticoid(panel, don, donors, n_perm, 5)))
    res = pd.DataFrame(rows)
    res.to_csv(out / "a108_timed.csv", index=False)
    print(res.drop(columns=["genes"], errors="ignore").round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
