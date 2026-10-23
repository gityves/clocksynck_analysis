"""V25.5: per-patient internal time against recorded resection time in a living hospital cohort (A25; post hoc).

Documented as the note V25.5 of A25 in docs/PLAN_internal_desync.md (archive) (2026-10-15, post hoc: first run exploratory, written unchanged; no
acceptance criterion set beforehand; hypothesis id V25.5). Data: GSE233242 (Li et al. 2024; 43 mastectomies, tumour and adjacent normal
bulk RNA-seq; `breast_timed_gse233242`: counts approximated from TPM x transcript length, rescaled to 30 M reads, 10 M and 60 M as
sensitivities; one patient with a one-digit id mismatch excluded, joined as a sensitivity). Host: the two compartments (normal, tumour) form
the host (a 2-compartment host: `Dataset.fit` returns the host states and Sigma_f but no estimates, so the frame shift is read from the
cross-covariance of the two compartments, local: not in clocksync). A25 estimator (`_internal_time.estimate`), flip resolved by the hospital
prior (c_day) and by curvature (c_curv). Truth: the recorded resection hour; the common offset (internal frame vs wall clock) is the circular mean
error and is removed before the per-patient error. Inference: the median absolute error against 2,000 shuffles of the times; patient
bootstrap (200, host refitted). Per tissue alone: tangent of each compartment's clock-plane state against the hour.
Estimands: Spearman of the tangent coordinate with the hour, circular correlation, median absolute error after the offset, flip agreement.
Reading (A25 note): under the hospital prior ordering and per-patient error pass; curvature can resolve the flip wrongly with a narrow window and
two compartments.

    python compute_breast_timed.py [--subsample N] [--seed S] [--n-boot B] [--n-perm P]
Outputs (paper/results/figure_1/output/windows/): v25_5_summary.csv, v25_5_patients.csv, v25_5_boot.csv, v25_5_gate.csv
Legacy: legacy/desync/validation/scripts/v25_5_breast_timed.py -> validation/output/v25_5_{summary,patients,boot}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as sps  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _internal_time as IT  # noqa: E402

FOLDER = "windows"
W, HOURS = IT.W, IT.HOURS
COMPS = ["normal", "tumour"]


def two_host_estimate(ds: cs.Dataset, comps=COMPS) -> tuple[dict, pd.DataFrame]:
    """local: not in clocksync. A25 window of a two-compartment host. Sigma_f and the donors' host states come from
    `Dataset.fit` (2 anchors); the frame shift (circular mean of arg z over the two compartments, in the frame of
    the first, `comps[0]`, the anchor with most donors) is half the phase of the second minus the first, read from their complex correlation."""
    ft = ds.fit(n_boot=0, gate=False)
    d = ds.data
    b = {k: pd.DataFrame(d.beta[d.compartment == k], index=d.resection[d.compartment == k]) for k in comps}
    common = b[comps[0]].index.intersection(b[comps[1]].index)
    o = np.angle(cs.stats.ccorr(b[comps[1]].loc[common].to_numpy(), b[comps[0]].loc[common].to_numpy()))
    e = IT.estimate(ft.host[["x", "y"]].to_numpy(), ft.sigma_f, o / 2, None, ft.host.index)
    return e, ft.host


def ccorr_h(a, b) -> float:
    return cs.stats.js_corr(np.asarray(a) * W, np.asarray(b) * W)


def offset_error(est, truth) -> tuple[float, float]:
    """Common offset (circular mean of est - truth, h) and the median |error| after removing it."""
    off = np.angle(np.exp(1j * W * IT.cdiff(est, truth)).mean()) * HOURS
    return float(off), float(np.median(np.abs(IT.cdiff(np.asarray(est) - off, truth))))


def score(e: dict, truth: np.ndarray) -> dict:
    tw = IT.true_window(truth)
    od, md = offset_error(e["donor_day_h"].to_numpy(), truth)
    oc, mc = offset_error(e["donor_h"].to_numpy(), truth)
    return dict(n=len(truth), spearman_t=sps.spearmanr(e["tangent_day"], truth)[0], ccorr_day=ccorr_h(e["donor_day_h"], truth),
                mae_day=md, offset_day=od, mae_curv=mc, offset_curv=oc, centre_day=e["c_day"], centre_curv=e["c_curv"],
                centre_true=tw["c_true"], flip_curv_agrees=bool(abs(IT.cdiff(e["c_day"], e["c_curv"])) < 6),
                w_est=e["w"], w_true=tw["w_true"], r=e["r"])


def run(pb, donors, waveform="human", n_boot=0, n_perm=2000, seed=0):
    ds = cs.Dataset(pb, anchors=COMPS, name="breast_timed", waveform=waveform)
    e, host = two_host_estimate(ds)
    hour = donors.resection_h
    truth = hour.reindex(host.index).to_numpy()
    row = score(e, truth)
    rng = np.random.default_rng(seed)
    row["p_spearman"] = sps.spearmanr(e["tangent_day"], truth)[1]
    null = np.array([offset_error(e["donor_day_h"].to_numpy(), rng.permutation(truth))[1] for _ in range(n_perm)])
    row["mae_day_null"] = float(np.median(null))
    row["p_mae_day"] = (1 + np.sum(null <= row["mae_day"])) / (n_perm + 1)
    # each tissue alone: its clock-plane state projected on the host's tangent (> 0 later); does it order the patients?
    ev, evec = np.linalg.eigh((e["sigma_f"] + e["sigma_f"].T) / 2)
    cand = np.arctan2(evec[1, 1], evec[0, 1]) + np.pi / 2
    use = cand if abs(IT.cdiff(e["c_day"], (cand + e["shift_h"] / HOURS) * HOURS)) < 6 else cand + np.pi
    e_t = np.array([-np.sin(use), np.cos(use)])
    d = ds.data
    for k in COMPS:
        i = np.flatnonzero(d.compartment == k)
        tk = (d.beta[i] - d.beta[i].mean(0)) @ e_t
        row[f"spearman_{k}_alone"] = sps.spearmanr(tk, hour.reindex(d.resection[i]).to_numpy())[0]
    boot = []
    t = ds.table.rename(columns={"resection": "donor"})
    by = {k: x for k, x in t.groupby("donor")}
    ids = np.array(sorted(by))
    for b in range(n_boot):
        pick = np.random.default_rng(1000 + b).choice(ids, len(ids))
        tb = pd.concat([by[p].assign(donor=f"{p}#{j}") for j, p in enumerate(pick)], ignore_index=True)
        try:
            eb, hb = two_host_estimate(cs.Dataset(tb, anchors=COMPS, name="breast_boot", waveform=waveform))
        except Exception:
            continue
        orig = hb.index.str.split("#").str[0]
        boot.append(score(eb, hour.reindex(orig).to_numpy()))
    return row, e, host, pd.DataFrame(boot)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--n-perm", type=int, default=None)
    args = ap.parse_args(argv)
    run_ = cli.setup(args)
    limit_threads(1)
    n_boot = args.n_boot if args.n_boot is not None else (20 if run_["quick"] else 200)
    n_perm = args.n_perm or (200 if run_["quick"] else 2000)
    out = out_dir(FOLDER, run_["quick"])
    b = data._builder("breast_timed_gse233242")
    pb, donors = b.build()
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(donors.index))[: args.subsample]
        pb, donors = pb[pb.donor.isin(keep)].reset_index(drop=True), donors.loc[keep]
    row, e, host, boot = run(pb, donors, n_boot=n_boot, n_perm=n_perm, seed=args.seed)
    per = pd.DataFrame(dict(patient=host.index, t=e["tangent_day"].to_numpy(), est_day=e["donor_day_h"].to_numpy(),
                            est_curv=e["donor_h"].to_numpy(), true_h=donors.resection_h.reindex(host.index).to_numpy()))
    per.to_csv(out / "v25_5_patients.csv", index=False)  # recorded hospital times of a living cohort: not controlled access
    boot.to_csv(out / "v25_5_boot.csv", index=False)
    for q in ("spearman_t", "mae_day", "mae_curv", "ccorr_day"):
        row[f"{q}_lo"], row[f"{q}_hi"] = np.nanpercentile(boot[q], [2.5, 97.5]) if len(boot) else (np.nan, np.nan)
    row["flip_curv_agrees_boot"] = boot.flip_curv_agrees.mean() if len(boot) else np.nan
    rows = [dict(run="primary (30 M reads, pooled waveform)") | row]
    if not args.subsample:
        for n in (10e6, 60e6):
            p2, d2 = b.build(n_reads=n)
            rows.append(dict(run=f"{n / 1e6:.0f} M reads") | run(p2, d2, n_perm=n_perm, seed=args.seed)[0])
        rows.append(dict(run="breast waveform") | run(pb, donors, waveform="Breast - Mammary Tissue", n_perm=n_perm, seed=args.seed)[0])
        p3, d3 = b.build(join_typo=True)
        rows.append(dict(run="W006063 joined as W006053") | run(p3, d3, n_perm=n_perm, seed=args.seed)[0])
    summ = pd.DataFrame(rows)
    summ.to_csv(out / "v25_5_summary.csv", index=False)
    g = cs.Dataset(pb, anchors=COMPS).gate(n_perm=200 if run_["quick"] else 1000, seed=args.seed)
    pd.DataFrame([dict(p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share, donors=g.n_donors)]).to_csv(
        out / "v25_5_gate.csv", index=False)
    pd.set_option("display.width", 220)
    print(summ.set_index("run").T.to_string())


if __name__ == "__main__":
    main()
