"""Multiplicity of the per-contrast tests: the distribution of the permutation p values of every A17 contrast x
compartment, Benjamini-Hochberg over all of them, and the spread-free readouts (d*, log ratios of V, z2, psi).

Descriptive (A17 reading rules; A20 motivation: "BH over 307 tests"). Reads the table of resilience
(`a17_contrast.csv`, learned rule, status ok, estimable p; run `resilience/compute_contrasts_a17.py`
first, same sample size). Per test set (d, d*): number of tests, p < 0.05 observed against expected, BH q < 0.05,
Holm-adjusted primaries, histogram of p in deciles; per target the number of tests and of p < 0.05; the spread
decomposition: median log ratio of V, z2 and psi, and the share of tests whose d and d* differ in sign.

    python compute_multiplicity.py [--subsample N]
Outputs (paper/results/figure_4/output/spread_multiplicity/): multiplicity_{summary,p_histogram,per_target,spread}.csv
Legacy: legacy/desync/common/unified.py (BH / Holm columns), a20_pool.py (q_bh_perm).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "spread_multiplicity"
SOURCE = "resilience"


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    t = pd.read_csv(out_dir(SOURCE, run["quick"]) / "a17_contrast.csv")
    out_rows, hist = [], []
    for rule in ("learned", "fixed"):
        u = t[(t.rule == rule) & (t.status == "ok") & t.p.notna()].reset_index(drop=True)
        for col, label in (("p", "d (coupling difference)"), ("p_star", "d* (common host spread)")):
            p = u[col].to_numpy(float)
            ok = np.isfinite(p)
            q = np.full(len(p), np.nan)
            q[ok] = cs.stats.bh(p[ok])
            out_rows.append(dict(rule=rule, test=label, tests=int(ok.sum()), targets=u.target.nunique(),
                                 p_lt_005=int((p[ok] < 0.05).sum()), expected_under_null=0.05 * ok.sum(),
                                 share_lt_005=float((p[ok] < 0.05).mean()), bh_q_lt_005=int((q[ok] < 0.05).sum()),
                                 min_p=float(np.nanmin(p)) if ok.any() else np.nan,
                                 primaries=int(u.primary.fillna(False).astype(bool).sum()),
                                 primaries_holm_lt_005=int((u.p_holm_primary < 0.05).sum())))
            h, _ = np.histogram(p[ok], bins=np.linspace(0, 1, 11))
            hist += [dict(rule=rule, test=label, bin_lo=i / 10, bin_hi=(i + 1) / 10, tests=int(n)) for i, n in enumerate(h)]
    pd.DataFrame(out_rows).to_csv(out / "multiplicity_summary.csv", index=False)
    pd.DataFrame(hist).to_csv(out / "multiplicity_p_histogram.csv", index=False)
    u = t[(t.rule == "learned") & (t.status == "ok") & t.p.notna()]
    per = u.groupby("target").agg(tests=("p", "size"), p_lt_005=("p", lambda x: int((x < 0.05).sum())),
                                  p_star_lt_005=("p_star", lambda x: int((x < 0.05).sum())),
                                  median_abs_d=("d", lambda x: float(np.nanmedian(np.abs(x)))),
                                  median_abs_d_star=("d_star", lambda x: float(np.nanmedian(np.abs(x)))))
    per.to_csv(out / "multiplicity_per_target.csv")
    sp = dict(tests=len(u), median_log_V_ratio=u.log_V_ratio.median(), median_log_z2_ratio=u.log_z2_ratio.median(),
              median_log_psi_ratio=u.log_psi_ratio.median(),
              share_d_dstar_opposite_sign=float((np.sign(u.d) != np.sign(u.d_star))[u.d.notna() & u.d_star.notna()].mean()),
              share_V_lower_in_cases=float((u.log_V_ratio < 0).mean()))
    pd.Series(sp).to_csv(out / "multiplicity_spread.csv", header=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(out_rows).round(3).to_string(index=False))
    print(pd.Series(sp).round(3).to_string())


if __name__ == "__main__":
    main()
