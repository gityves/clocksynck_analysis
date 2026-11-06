"""A204b scoring: TimeTeller arm with matched tissue means (gate, H204b.4, power per window and c_T, MDE).

    $(conda info --base)/envs/sccc/bin/python paper/exploratory/A204_orderer_uncoupling/scripts/score_a204b.py
Reads output/a204b/runs (compute_a204.py --only tt --mesor-sd 0); writes output/a204b/{gate,tests,tt_per_rep,tt_cells}.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pandas as pd  # noqa: E402

import compute_a204 as A  # noqa: E402
import score_a204 as SC  # noqa: E402

out = A.OUT / "a204b"
tt = SC.tt_rows(out / "runs", 30)
tt.to_csv(out / "tt_per_rep.csv", index=False)
tc = tt[~tt.failed].groupby(["window", "c_T"]).agg(
    power=("slope_p", lambda x: float((x < 0.05).mean())), n_detect=("slope_p", lambda x: int((x < 0.05).sum())),
    median_slope=("slope", "median"), median_r=("r", "median"), median_circ_lin_r=("circ_lin_r", "median"),
    median_abs_err_h=("median_abs_err_h", "median"), share_7_20=("share_7_20", "median"), n_pred=("n", "median"),
    reps=("rep", "size")).reset_index()
tc.to_csv(out / "tt_cells.csv", index=False)
err = float(tt.query("window == 'F' and c_T == 0.9").median_abs_err_h.median())
gate = pd.DataFrame([dict(check="timeteller_F_c0.9_median_abs_err_h", value=err, threshold="<= 3", passes=err <= 3)])
gate.to_csv(out / "gate.csv", index=False)
T9 = tt.query("window == 'W' and c_T == 0.9 and not failed")
k = int((T9.slope_p < 0.05).sum())
verdict = ("PASS" if k >= 24 else "FAIL") if gate.passes.all() else "not scored (gate failed)"
tests = pd.DataFrame([dict(h="H204b.4", readout="W c0.9 n 108: TimeTeller slope p < 0.05", k=k, n=len(T9), rule=">= 24",
                           p=SC.binom_ge(k, len(T9)), verdict=verdict)])
tests.to_csv(out / "tests.csv", index=False)
print(gate.to_string()); print(tests.to_string()); print(tc.round(3).to_string())
print("MDE:", tc.query("power >= 0.8").groupby("window").c_T.min().to_dict())
