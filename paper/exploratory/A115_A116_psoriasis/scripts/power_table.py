"""One-sided (alpha 0.05) normal-approximation power for assumed true effects at the null SDs of a115_null_pooled.csv and a116_null.csv.
Output: output/power_table.csv"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm
O = Path(__file__).resolve().parents[1] / "output"
a = pd.read_csv(O / "a115_null_pooled.csv").query("analysis == 'psoriasis_pooled'")
b = pd.read_csv(O / "a116_null.csv")
rows = []
for _, r in a.iterrows(): rows.append(dict(design="A115 " + r.set + " " + r.scheme, sd_h=r.null_sd_h))
for _, r in b.iterrows(): rows.append(dict(design="A116 " + r.analysis, sd_h=r.null_sd_h))
t = pd.DataFrame(rows)
for e in (0.5, 1.0, 1.5): t[f"power_effect_{e}h"] = norm.cdf(e / t.sd_h - 1.645)
t["mde_80pct_h"] = (1.645 + 0.842) * t.sd_h
t.to_csv(O / "power_table.csv", index=False); print(t.round(3).to_string())
