"""A129 gap: conditional distribution of (clocksync change - oracle change) among replicates whose oracle change is within
+-0.04 of the observed -0.186, pooled over the levels of each mechanism; share of those with gap <= observed -0.107 (-0.293 + 0.186).
Reads output/mech_raw.csv; writes output/gap_conditional.csv. Descriptive. Also copies the real-data d vs d* (A111, Liver, AL arm)."""
from pathlib import Path
import numpy as np, pandas as pd
O = Path(__file__).resolve().parent / "output"
R = pd.read_csv(O / "mech_raw.csv")
OBS_OR, OBS_CS = -0.1859, -0.2934
rows = []
for m, d in R.groupby("mech"):
    s = d[(d.chg_or - OBS_OR).abs() <= 0.04]
    rows.append(dict(mech=m, n_reps_in_band=len(s), n_levels=s.level.nunique(), mean_oracle_change=s.chg_or.mean() if len(s) else np.nan,
                     mean_cs_change=s.chg_cs.mean() if len(s) else np.nan, mean_gap=s.gap.mean() if len(s) else np.nan,
                     sd_gap=s.gap.std() if len(s) > 1 else np.nan, observed_gap=OBS_CS - OBS_OR,
                     share_gap_le_observed=float((s.gap <= OBS_CS - OBS_OR).mean()) if len(s) else np.nan,
                     ratio_of_means=s.chg_cs.mean() / s.chg_or.mean() if len(s) else np.nan))
G = pd.DataFrame(rows); G.to_csv(O / "gap_conditional.csv", index=False)
pd.set_option("display.width", 250); print(G.round(3).to_string())
