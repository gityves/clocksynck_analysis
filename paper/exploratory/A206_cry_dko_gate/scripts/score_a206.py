"""A206: single-set clock gate (`paired.reference_gate`, A52/A83 setting) in Cry1/Cry2 DKO restricted-feeding liver.

Registered in paper/docs/registered_experiments.md (A206). Data: $SCCC_DATA/vollmers2009_GSE13093 (parse_gse13093.py).
Outputs (../output/): gate.csv, cosinor.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2].parent))
sys.path.insert(0, str(HERE.parents[2]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import reference_gate  # noqa: E402
from common.paths import DATA  # noqa: E402

D = DATA / "vollmers2009_GSE13093"
OUT = HERE.parent / "output"
WT, DKO_RF, DKO_AL = ("Temporally restricted feeding in the liver of WT mice",
                      "Temporally restricted feeding in Liver of Cry1,Cry2 KO mice", "Ad libitum feeding in Liver of Cry1,Cry2 KO mice")

wf = cs.gtex()
E = pd.read_csv(D / "expression_clock.csv", index_col="gsm")[list(wf.genes)]
s = pd.read_csv(D / "samples.csv").set_index("gsm")
sets = {"WT_RF_even": (s.group == WT) & (s.ct % 2 == 0), "WT_RF_odd": (s.group == WT) & (s.ct % 2 == 1),
        "DKO_RF": s.group == DKO_RF, "DKO_AL": s.group == DKO_AL}
rows, cos = [], []
for k, m in sets.items():
    y = E.loc[s.index[m]].to_numpy(float)
    g = reference_gate(y, wf, n_perm=200, rng=206)
    rows.append(dict(set=k, n=int(m.sum()), **g, passes=g["p"] < 0.05))
    t = s.ct[m].to_numpy(float) * 2 * np.pi / 24
    X = np.column_stack([np.ones(len(t)), np.cos(t), np.sin(t)])
    B, *_ = np.linalg.lstsq(X, y, rcond=None)
    for j, gene in enumerate(wf.genes):
        cos.append(dict(set=k, gene=gene, amp_log2=float(np.hypot(B[1, j], B[2, j])),
                        acrophase_ct=float(np.arctan2(B[2, j], B[1, j]) * 24 / (2 * np.pi) % 24)))
OUT.mkdir(exist_ok=True)
gate = pd.DataFrame(rows)
gate.to_csv(OUT / "gate.csv", index=False)
pd.DataFrame(cos).to_csv(OUT / "cosinor.csv", index=False)
print(gate.round(4).to_string())
print(pd.DataFrame(cos).pivot(index="gene", columns="set", values="amp_log2").round(2).to_string())
