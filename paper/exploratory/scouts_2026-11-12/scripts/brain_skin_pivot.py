import pandas as pd
O="/home/maxine/Documents/paychere/scCircadianOverview/paper/exploratory/scouts_2026-11-12/output/"
G="ARNTL CLOCK NPAS2 PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 DBP RORA TEF CIART".split()
d=pd.read_csv(O+"brain_skin_clock_probe_detection.csv")
out=[]
for ds,g in d.groupby("dataset"):
    g=g[g.probe!="absent"]
    t=g.groupby(["gene","compartment"]).percentile.agg(["min","max"]).round(0).astype(int)
    t=t.assign(v=t["min"].astype(str)+"-"+t["max"].astype(str)).v.where(t["min"]!=t["max"],t["max"].astype(str)).unstack().reindex(G)
    t.insert(0,"n_probes",g.groupby("gene").probe.nunique().reindex(G).fillna(0).astype(int))
    t.to_csv(O+f"clock_percentile_{ds}.csv"); print("##",ds); print(t.to_string())
    print({k:list(v) for k,v in g.groupby("gene").probe.unique().items()} if ds.startswith("GSE7") else "")
