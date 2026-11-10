"""POST HOC feasibility look (review 2026-11-10), not a registered analysis.
Do d_star / log_z2_ratio stay at 0 when two groups have the same true loading and noise but different host spreads,
and the host is measured with realistic noise? Host fitted on the pooled anchors (Dataset.contrast)."""
import sys, warnings
import numpy as np, pandas as pd
import clocksync as cs
from clocksync.config import limit_threads
from clocksync.waveform import HOURS
limit_threads(2); warnings.simplefilter("ignore")
WF = cs.gtex(); ANCH = ["fibroblast", "endothelial", "pericyte"]
Z = {"fibroblast": 1.0, "endothelial": 0.8 * np.exp(1j * 2 / HOURS), "pericyte": 1.25 * np.exp(-1j * 1.5 / HOURS),
     "T cell": 0.6 * np.exp(-1j * 3 / HOURS)}

def cohort(n, seed, window, psi_a, prefix, psi_t=0.05):
    comps = {k: cs.Compartment(Z[k], psi=psi_a, mesor_sd=0.0) for k in ANCH}
    comps["T cell"] = cs.Compartment(Z["T cell"], psi=psi_t, mesor_sd=0.0)
    df, _ = cs.simulate(WF, comps, n=n, seed=seed, window_h=window, jitter_h=1.0, tau=0.0)
    df["resection"] = prefix + df.resection.str[1:]
    return df.rename(columns={"resection": "donor"})

rows = []
for psi_a in (0.002, 0.05, 0.15):
    for wb in ((6, 18), (8, 16), (10, 14)):
        for s in range(int(sys.argv[1])):
            a = cohort(80, 1000 + s, (6, 18), psi_a, "a"); b = cohort(80, 5000 + s, wb, psi_a, "b")
            g = pd.Series(["control"] * 80 + ["case"] * 80, index=sorted(a.donor.unique()) + sorted(b.donor.unique()))
            t = cs.Dataset(pd.concat([a, b], ignore_index=True), anchors=ANCH).contrast(
                g, n_perm=0, n_boot=0, blocks="none").table
            for k in ("T cell", "fibroblast"):
                r = t.loc[k]
                rows.append(dict(psi_anchor=psi_a, case_window=f"{wb[0]}-{wb[1]}", comp=k, d=r.d, d_star=r.d_star,
                                 log_V=r.log_V_ratio, log_z2=r.log_z2_ratio, log_psi=r.log_psi_ratio,
                                 c_ref=r.c_control, c_case=r.c_case))
o = pd.DataFrame(rows)
print(o.groupby(["psi_anchor", "case_window", "comp"]).agg(["mean", "sem"]).round(3).to_string())
o.to_csv(sys.argv[2], index=False)
