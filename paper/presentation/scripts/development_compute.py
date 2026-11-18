"""Aggregates for the development slide (R4): pair-level view of the A307 / A309 W statistic and a post hoc check of the
mouse kidney's negative coherence.

pairs.csv: for every (signal-induced, core) clock-gene pair (PER1, PER2, NR1D1 x the 11 core genes) of the readouts
  where W was scored (human adrenal cortex, lung, liver: A307 rows and adjustment; mouse hepatocytes: A309), the
  correlation across individuals (r) and the relation expected from the adult clock, cos(adult phase difference).
  W = mean(r x expected). Reproduces A307 / A309 W (asserted).
kidney.csv (post hoc, descriptive): mouse nephron progenitors (A309 rows): signed coherence without each clock gene,
  each gene's log SD across embryos and its time-of-day fit (cosinor on the sampling hour, 00 / 06 / 12 / 18 h).

    python scripts/development_compute.py
Output: presentation/output/development/{pairs,kidney}.csv.
"""
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(EXP / "A307_fetal_gc_pattern" / "scripts"))
sys.path.insert(0, str(EXP / "A309_mouse_embryo_organs" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import lib_a309 as M  # noqa: E402
import score_a307 as H  # noqa: E402
from clocksync.decompose import signed_coherence  # noqa: E402
from common import development as dev  # noqa: E402

OUT = PAPER / "presentation" / "output" / "development"


def pairs(Y, genes, phi, si, core):
    R = pd.DataFrame(np.corrcoef(Y, rowvar=False), index=genes, columns=genes)
    return [dict(a=a, b=b, r=R.loc[a, b], expected=np.cos(phi[a] - phi[b])) for a in si for b in core]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    data = H.load()
    for org in H.ORGANS:
        counts, obs, cov = data[org]
        y, _ = H.adjusted(counts, obs, cov, H.GENES)
        rows += [dict(organ=org, species="human", **p) for p in pairs(y, H.GENES, H.PHI, H.SI, H.CORE)]
    m = M.late_rows()
    o = m[m.readout == "liver"]
    counts = M.genome(m).loc[o.index]
    y, _ = H.adjusted(counts, o, M.cov, M.GENES)
    rows += [dict(organ="liver", species="mouse", **p) for p in pairs(y, M.GENES, dict(zip(M.GENES, M.WF.psi)), M.SI, M.CORE)]
    d = pd.DataFrame(rows)
    W = d.assign(w=d.r * d.expected).groupby(["species", "organ"]).w.mean()
    ref = pd.read_csv(EXP / "A307_fetal_gc_pattern/output/tests.csv").query("statistic == 'W'").set_index("organ").value
    for org, v in ref.items():
        assert np.isclose(W["human", org], v), (org, W["human", org], v)
    ref = pd.read_csv(EXP / "A309_mouse_embryo_organs/output/gc_pattern.csv").query("statistic == 'W'").value.iloc[0]
    assert np.isclose(W["mouse", "liver"], ref), (W["mouse", "liver"], ref)
    d.to_csv(OUT / "pairs.csv", index=False)
    print(W.round(4).to_string())

    k = m[m.readout == "kidney"]
    tab, _ = dev.adjust(k, M.GENES, M.cov)
    y = dev._y(tab, M.WF)
    t = k.t_clock.to_numpy(float) * 2 * np.pi / 24
    X = np.c_[np.cos(t), np.sin(t)]
    out = []
    for j, g in enumerate(M.GENES):
        keep = [x for x in M.GENES if x != g]
        b, *_ = np.linalg.lstsq(X, y[:, j], rcond=None)
        out.append(dict(gene=g, C_without=signed_coherence(np.delete(y, j, 1), M.WF.subset(keep)),
                        log_sd=y[:, j].std(ddof=1), mean_count=k[g].mean(), tod_amplitude=np.hypot(*b),
                        tod_peak_h=np.arctan2(b[1], b[0]) * 24 / (2 * np.pi) % 24,
                        adult_peak_h=dict(zip(M.GENES, M.WF.psi))[g] * 24 / (2 * np.pi) % 24))
    out = pd.DataFrame(out).assign(C_all=signed_coherence(y, M.WF), embryos=len(k))
    out.to_csv(OUT / "kidney.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
