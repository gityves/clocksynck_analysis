"""T3 (HANDOFF_2026-10-31): one table of the true-time validations, each with its readout, error and prior-only baselines.

Prior-only baselines use no expression: "window centre" predicts the centre of the sampling window for every sample;
"best constant" is the single time minimising the median absolute circular error against the truth (an oracle, it uses
the truth). For an offset readout the prior-only prediction is no offset (0 h). Rows are read from the current outputs.
    python truetime_table.py   ->   output/truetime_table.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
from common.paths import exploratory_dir  # noqa: E402

F1 = PAPER / "results/figure_1/output"
EX = PAPER / "exploratory"


def circ_err(a, b):
    return np.abs((np.asarray(a) - np.asarray(b) + 12) % 24 - 12)


def best_constant(t):
    grid = np.arange(0, 24, 0.05)
    return min(np.median(circ_err(t, c)) for c in grid)


def main():
    rows = []
    a36 = pd.read_csv(F1 / "species/a36_readouts.csv").set_index("readout")
    host, off = a36.iloc[0], a36.loc["2 untimed offset vs timed -arg g (all tissues)"]
    rows += [dict(dataset="baboon (A36)", truth="ZT of sacrifice, full day", n="12 animals",
                  readout="host phase vs ZT (R)", value=host.R, error_h=host.median_abs_err_h,
                  window_centre_h=np.nan, best_constant_h=np.nan, note="full-day sampling: order and offsets are the readout"),
             dict(dataset="baboon (A36)", truth="timed organ offsets", n="63 tissues", readout="organ offsets (r_FL)",
                  value=off.r_fl, error_h=off.median_abs_err_h, window_centre_h=np.nan, best_constant_h=np.nan, note="")]
    v25 = pd.read_csv(EX / "A25_prior_baseline/output/prior_baseline.csv").set_index("estimator").median_abs_error_h
    w25 = pd.read_csv(F1 / "windows/v25_5_summary.csv").iloc[0]
    pick = lambda key: v25[[k for k in v25.index if k.startswith(key)][0]]
    rows.append(dict(dataset="V25.5 breast resections", truth="recorded resection time", n="42 patients",
                     readout="patient order (Spearman)", value=w25.spearman_t, error_h=pick("method, offset"),
                     window_centre_h=pick("always 13:00"), best_constant_h=pick("best constant"),
                     note=f"error with no use of truth (circular mean set to 13:00) {pick('method, circular'):.2f} h"))
    a113 = pd.read_csv(EX / "A113_living_time/output/a113_summary.csv")
    s = a113[(a113.dataset == "adipose GSE104674") & a113.arm.str.startswith("a")].iloc[0]
    ph = pd.read_csv(EX / "A113_living_time/output/a113_phases.csv")
    t = ph[ph.dataset == "adipose GSE104674"].time_h.to_numpy()
    centre = (t.min() + t.max()) / 2
    rows.append(dict(dataset="A113 adipose", truth="timed biopsies", n=f"{int(s.n_subjects)} men, {int(s.n_samples)} biopsies",
                     readout="internal time vs clock time (rho)", value=s.rho, error_h=s.med_abs_err_h,
                     window_centre_h=np.median(circ_err(t, centre)), best_constant_h=best_constant(t),
                     note=f"sampling {t.min():g}-{t.max():g} h; error after the circular offset fit"))
    sk = pd.read_csv(F1 / "method/skin_truth_summary.csv").iloc[0]
    rows.append(dict(dataset="A37 living skin", truth="timed epidermis-dermis offset", n=f"{int(sk.n_subjects)} subjects",
                     readout="untimed offset (h)", value=sk.offset_untimed, error_h=abs(sk.offset_untimed - sk.offset_timed),
                     window_centre_h=abs(sk.offset_timed), best_constant_h=np.nan,
                     note=f"timed offset {sk.offset_timed:.2f} h; prior-only = no offset"))
    v3 = pd.read_csv(EX / "A174_v3_heart/output/v3_scores.csv").set_index("analysis").loc["primary"]
    rows.append(dict(dataset="V3 heart surgery (A174)", truth="AM / PM window (4 h each)", n="16 patients",
                     readout="PM - AM host displacement angle (h)", value=v3.angle_h, error_h=abs(v3.error_h),
                     window_centre_h=np.nan, best_constant_h=np.nan,
                     note=f"separation p {v3.p_separation:.4f}; no prior-only predictor of a direction (chance error 6 h)"))
    out = exploratory_dir("T3_truetime")
    pd.DataFrame(rows).to_csv(out / "truetime_table.csv", index=False)
    print(pd.DataFrame(rows).drop(columns="note").round(3).to_string())


if __name__ == "__main__":
    main()
