"""A192 post hoc (chosen after seeing that RORC keeps an association with G2M in the model): is RORC's signal clock or immune?

RORC's correlation with an immune score (PTPRC, CD3E, CD2, CD247 mean z), with the other clock genes (S without RORC) and with
G2M; and its pooled HR adjusted for G2M, ESR1, ERBB2, then + the immune score. Output `output/a192_posthoc_rorc.csv`.
"""
import numpy as np
import pandas as pd
import score_a192 as A

E, Z, c, _ = A.load()
imm = [g for g in ("PTPRC", "CD3E", "CD2", "CD247") if g in Z.columns]
c["immune"] = Z[imm].mean(axis=1)
c["S6"] = Z[[g for g in A.PROTECTIVE if g != "RORC"]].mean(axis=1)
rows = [dict(test=f"r(RORC, {k})", r=float(np.corrcoef(Z.RORC, c[k])[0, 1])) for k in ("immune", "S6", "g2m", "ESR1")]
cov = ["g2m", "ESR1", "ERBB2"]
for lab, cv in (("RORC + G2M, ESR1, ERBB2", cov), ("RORC + G2M, ESR1, ERBB2 + immune", cov + ["immune"]),
                ("immune + G2M, ESR1, ERBB2", None), ("S without RORC + G2M, ESR1, ERBB2", None)):
    x = c.immune if lab.startswith("immune") else c.S6 if lab.startswith("S without") else Z.RORC
    rows.append(dict(test=lab, **A.fit(c, x, cv if cv else cov)))
out = pd.DataFrame(rows)
out.to_csv(A.OUT / "a192_posthoc_rorc.csv", index=False)
print(out.round(4).to_string(index=False))
