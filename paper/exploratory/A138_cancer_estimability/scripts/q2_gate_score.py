"""A138 Q2 simulator gate verdict (gate only; no H138.2-H138.7 readout).

Registered rule: in the full arm without heterogeneity (n = 300), CoFE and CYCLOPS recover the truth within their
published GTEx-like error; CoFE clock-gene median acrophase error <= 2 h. Statistic per replicate and method: the median
over the 15 clock genes of |acrophase error| after oracle alignment (`q2_gate.py`). Pass: median over replicates <= 2 h
(CYCLOPS: same bound; GTEx n = 80 full day 1.73 h). Failed fits count as failures (not dropped); the median is over
replicates with a fit, the failure count is reported beside it. Verdict: CoFE as published (genome-wide) and CYCLOPS
must both pass; clock-only CoFE is reported.

    python q2_gate_score.py [--dir output/q2_gate]      writes gate_summary.csv, gate_verdict.csv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BOUND = 2.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=str(HERE.parent / "output" / "q2_gate"))
    d = Path(ap.parse_args().dir)
    t = pd.concat([pd.read_csv(f) for f in sorted((d / "reps").glob("*.csv"))], ignore_index=True)
    rows = []
    for m, g in t.groupby("method"):
        ok = g[g.status == "ok"]
        e = ok.clock_err_med_h
        rng = np.random.default_rng(0)
        boot = [np.median(rng.choice(e.to_numpy(), len(e))) for _ in range(2000)] if len(e) else [np.nan]
        rows.append(dict(method=m, n_rep=len(g), n_fail=int((g.status != "ok").sum()), median_clock_err_h=e.median(),
                         ci_lo=np.percentile(boot, 2.5), ci_hi=np.percentile(boot, 97.5), q25=e.quantile(0.25),
                         q75=e.quantile(0.75), share_le_2h=float((e <= BOUND).sum() / len(g)),
                         median_err_vs_waveform_h=ok.clock_err_wf_med_h.median(), median_mae_h=ok.mae_h.median(),
                         median_circ_corr=ok.circ_corr.median(), median_fit_s=g.fit_s.median(),
                         share_interior_min=ok.interior_min.mean() if "interior_min" in ok else np.nan,
                         gate_pass=bool(len(e) and e.median() <= BOUND and (g.status != "ok").mean() < 0.5)))
    s = pd.DataFrame(rows)
    s.to_csv(d / "gate_summary.csv", index=False)
    p = s.set_index("method").gate_pass
    v = pd.DataFrame([dict(cofe_genome=p.get("cofe_genome"), cyclops=p.get("cyclops"), cofe_clock_only=p.get("cofe_clock"),
                           gate_pass=bool(p.get("cofe_genome") and p.get("cyclops")),
                           complete=bool((s.n_rep == s.n_rep.max()).all() and len(s) == 3))])
    v.to_csv(d / "gate_verdict.csv", index=False)
    print(s.round(3).to_string())
    print(v.to_string())


if __name__ == "__main__":
    main()
