"""A258 evaluation: V258.1 (bias of the pooled estimator by design cell), V258.2 (acceptance of the arm prior), V258.3
(expected bias of every row of the equivalence families at its observed host-spread ratio).

    python evaluate_a258.py      (after simulate_a258.py --part bias|size|power)

Choices not fixed by the registration
- V258.2 uses the observed log V ratio of the pooled estimator (what a real table reports) to select the cells with log V
  ratio >= -1.0, and the target and the anchor ("fibroblast") alike; the bias criterion is applied to the cell means.
- V258.3: the expected bias of the pooled estimator at a row's observed log V ratio is interpolated linearly in |log V ratio|
  from the simulation cells of one anchor noise (0.05 and 0.15 reported; donors per arm 80), sign-reversed for a wider case
  arm (the narrower arm reads low). Rows of A164 / A202 / A208 carry no spread columns and are not mapped.
"""
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "output"
EXP = Path(__file__).resolve().parents[2]
TABLES = {"A231": EXP / "A231_age_coupling/output/score_tissues.csv", "A233": EXP / "A233_chronic_conditions/output/score_tissues.csv",
          "A235": EXP / "A235_metabolic_offset/output/score_tissues.csv"}


def main():
    b = pd.read_csv(OUT / "summary_bias.csv")
    key = ["psi_anchor", "case_window", "n", "comp"]
    wide = b.pivot_table(index=key, columns="prior", values=["d_star", "d_star_sem", "log_V", "rmse"]).reset_index()
    wide.columns = ["_".join(c).strip("_") for c in wide.columns]
    wide.to_csv(OUT / "v258_1_bias.csv", index=False)
    sel = wide[wide.log_V_pooled >= -1.0]
    crit_bias = bool((sel.d_star_arm.abs() <= 0.02).all())
    eq = wide[wide.case_window == "6-18"]
    rmse_ratio = float((eq.rmse_arm / eq.rmse_pooled).max())
    pw = pd.read_csv(OUT / "summary_power.csv")
    pw = pw[pw.comp == "T cell"].pivot_table(index=["psi_anchor", "case_window"], columns="prior", values="reject")
    power_loss = float((pw["pooled"] - pw["arm"]).max())
    sz = pd.read_csv(OUT / "summary_size.csv")
    size = sz.pivot_table(index=["psi_anchor", "case_window", "comp"], columns="prior", values="reject").reset_index()
    v2 = dict(cells_checked=len(sel), worst_abs_bias_arm=float(sel.d_star_arm.abs().max()),
              worst_cell=str(sel.loc[sel.d_star_arm.abs().idxmax(), key].to_dict()), cells_over_0_02=int((sel.d_star_arm.abs() > 0.02).sum()),
              cells_over_0_02_n_ge_80=int((sel[sel.n >= 80].d_star_arm.abs() > 0.02).sum()), bias_criterion=crit_bias,
              rmse_ratio_equal_spread_max=rmse_ratio, rmse_criterion=bool(rmse_ratio <= 1.10), power_loss_max=power_loss,
              power_criterion=bool(power_loss <= 0.05))
    v2["accepted"] = bool(v2["bias_criterion"] and v2["rmse_criterion"] and v2["power_criterion"])
    pd.DataFrame([v2]).to_csv(OUT / "v258_2_acceptance.csv", index=False)
    size.to_csv(OUT / "v258_size.csv", index=False)
    pw.reset_index().to_csv(OUT / "v258_power.csv", index=False)
    # V258.3
    rows = []
    for psi in (0.05, 0.15):
        cur = wide[(wide.psi_anchor == psi) & (wide.n == 80) & (wide.comp == "T cell")].sort_values("log_V_pooled")
        x, y = np.abs(cur.log_V_pooled.to_numpy()), cur.d_star_pooled.to_numpy()
        o = np.argsort(x)
        x, y = x[o], y[o]
        for fam, path in TABLES.items():
            t = pd.read_csv(path)
            lv = np.log(t.V_case / t.V_control)
            e = np.sign(-lv) * np.interp(np.abs(lv), x, -y) * -1
            for (_, r), l, ex in zip(t.iterrows(), lv, e):
                rows.append(dict(family=fam, psi_anchor=psi, design=r.get("design", ""), cohort=r.cohort, tissue=r.tissue, d_star=r.d_star,
                                 log_V_ratio=l, expected_bias=ex, flag=bool(abs(ex) > 0.05)))
    F = pd.DataFrame(rows)
    F.to_csv(OUT / "v258_3_rows.csv", index=False)
    S = F.groupby(["family", "psi_anchor"]).agg(rows=("flag", "size"), flagged=("flag", "sum"), median_abs_logV=("log_V_ratio", lambda v: float(np.median(np.abs(v)))),
                                               max_abs_expected=("expected_bias", lambda v: float(np.max(np.abs(v))))).reset_index()
    S.to_csv(OUT / "v258_3_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.Series(v2).to_string(), size.round(3).to_string(), pw.round(3).to_string(), S.round(3).to_string(), sep="\n\n")


if __name__ == "__main__":
    main()
