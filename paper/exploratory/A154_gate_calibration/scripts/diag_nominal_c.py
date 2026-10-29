"""A154 post hoc diagnostic (written after V154.3 was scored): is the anchored estimate of c biased at large n when the anchors are
themselves weakly coupled? Target of nominal coupling 0.3 (z = 1), K anchors of the same size, 600 donors, no gene noise, no common
shift, depth 1e10 (no Poisson noise); anchor psi 0.002 (precise, as the A17 acceptance), 0.05, or equal to the target's.
Output: ../output/diag_nominal_c.csv (mean of 4 cohorts)."""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2])); sys.path.insert(0, str(HERE.parents[1] / "A155_anchor_self_inclusion/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402
from compute_a155 import psi_for  # noqa: E402

limit_threads(1)
warnings.simplefilter("ignore")
pt, rows = psi_for(0.3, 1.0), []
for apsi in (0.002, 0.05, pt):
    for na in (2, 3):
        est = []
        for s in range(4):
            comps = {f"a{k}": Compartment(z=1.0, psi=apsi) for k in range(na)} | {"k": Compartment(z=1.0, psi=pt)}
            tab, _ = cs.simulate(cs.gtex(), comps, n=600, depth=1e10, seed=s, tau=0.0, shift_sd=0.0)
            ds = cs.Dataset(tab.rename(columns={"resection": "donor"}), anchors=[f"a{k}" for k in range(na)])
            est.append(ds.fit(gate=False, targets=["k"], n_boot=0).table.loc["k", "coupling"])
        rows.append(dict(anchor_psi=apsi, n_anchors=na, nominal_c=0.3, c_hat_mean=float(np.mean(est)), cohorts=4))
pd.DataFrame(rows).to_csv(HERE.parent / "output" / "diag_nominal_c.csv", index=False)
print(pd.DataFrame(rows))
