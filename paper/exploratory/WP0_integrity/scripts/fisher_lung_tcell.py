"""WP0 item 1a: recompute the Fisher combined p of the fibrotic-lung T-cell coupling drop (A17 contrast tables)."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

ARCH = Path("/home/maxine/Documents/paychere/scCircadianCancer")
F = ARCH / "legacy/desync/common_output/a17_unified_contrast.csv"
OUT = Path(__file__).resolve().parents[1] / "output"
d = pd.read_csv(F)
rows = []
for target, cohort in [("hlca_fibrotic", "HLCA"), ("ild_natri", "Natri")]:
    i = d.index[(d.target == target) & (d.compartment == "lung:T")][0]
    r = d.loc[i]
    rows.append(dict(cohort=cohort, file=F.name, csv_row=int(i) + 2, target=target, column="p_perm_learned",
                     p=r.p_perm_learned, c_control=r.c_control_learned, c_case=r.c_case_learned,
                     diff=r.diff_learned, sidedness="two-sided (abs(null) >= abs(diff); clocksync _core/host.py:349)"))
p = np.array([r["p"] for r in rows])
chi = -2 * np.log(p).sum()
pf = stats.chi2.sf(chi, 2 * len(p))
rows.append(dict(cohort="Fisher combined", file=F.name, column="p_perm_learned", p=pf, diff=np.nan,
                 sidedness=f"chi2={chi:.3f}, df={2*len(p)}; with rounded 0.027, 0.009: {stats.chi2.sf(-2*np.log([.027,.009]).sum(), 4):.5f}; archive quotes 0.002"))
pd.DataFrame(rows).to_csv(OUT / "fisher_lung_tcell.csv", index=False)
print(pd.DataFrame(rows)[["cohort", "csv_row", "p"]])
