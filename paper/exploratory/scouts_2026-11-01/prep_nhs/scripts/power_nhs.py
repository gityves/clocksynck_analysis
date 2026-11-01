"""NHS breast (GSE93601, GPL22920): gates and power for the replication of A175 H175.1 (G2M tertiles) and the ER contrast (H175.2),
labels only. Mirrors A177 `power_a177.py` (imported: `pair_gate`; A175 `Contrast`, `Cancer`; `clocksync.paired`). Microarray adaptation as
A177 (`counts=False`). NO observed pairing statistic, stratified or not, is evaluated (`Contrast.observed` is never called).
  H175.1  D = pairing(low G2M tertile) - pairing(high tertile)           (arms 0 / 2 of the tertile labels)
  H175.2  D_ER = pairing(ER+) - pairing(ER-), pairs with missing ER excluded (arms 1 / 0)
MDE = q95 + 0.84 SD of the label-permutation null (`q3_lib.summarise`; = 2.49 x SD for a normal null); TCGA references (A175 tests.csv):
D = +0.158 (H175.1, null SD 0.0772, MDE 0.1938), as the stop rule of A175 (ER contrast not estimable if MDE > 0.3; TCGA BRCA ER: null SD 0.169, MDE 0.408, not run; a replication
is only informative if MDE < the TCGA effect).
    python power_nhs.py [--n-perm 10000] [--n-gate 2000]
Outputs (output/): gates.csv, power.csv, power_arms.csv"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
EXP = HERE.parents[1]
sys.path.insert(0, str(EXP / "A177_c1_replication" / "scripts"))
import power_a177 as P  # noqa: E402  (pair_gate, log; puts A175 / A138 scripts and paper on sys.path)
import run_a175 as A  # noqa: E402
import q3_lib as q  # noqa: E402
from clocksync.paired import cohort_states, reference_gate  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE / "output"
KEY = "nhs_gpl22920"
GATE_SEED = 1775
TCGA_D_H1751 = 0.158


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-gate", type=int, default=2000)
    args = ap.parse_args(argv)
    cli.setup(args)
    t, n = pd.read_parquet(OUT / f"{KEY}_tumour.parquet"), pd.read_parquet(OUT / f"{KEY}_normal.parquet")
    p = pd.read_csv(OUT / f"{KEY}_pairs.csv")
    wf = q.WF.drop([g for g in q.WF.genes if g not in t.columns])
    g = list(wf.genes)
    bt, vt, bn, vn = cohort_states(t, n, False, wf)
    assert not vt.any() and not vn.any()
    C = A.Cancer(KEY, dict(n=len(t), raw=(bt, vt, bn, vn)))
    yn, yt = n[g].to_numpy(float) * np.log(2), t[g].to_numpy(float) * np.log(2)
    r = reference_gate(yn, wf, args.n_gate, np.random.default_rng([GATE_SEED, 0]))
    pg = P.pair_gate(yt, yn, wf, args.n_gate, np.random.default_rng([GATE_SEED, 10]))
    gates = [dict(cohort=KEY, n_pairs=len(t), n_clock_genes=len(g), missing_clock_genes=";".join(x for x in q.WF.genes if x not in g),
                  ref_share=r["share"], ref_share_null_mean=r["share_null_mean"], gate_ref_p=r["p"], gate_ref_pass=bool(r["p"] < 0.05),
                  **pg, gate_pair_pass=bool(pg["p_pairing"] < 0.05 and pg["p_acrophase"] < 0.05), n_perm_gate=args.n_gate)]
    P.log(KEY, {k: gates[0][k] for k in ("gate_ref_p", "p_pairing", "p_acrophase")})
    rows, arms = [], []
    for name, lab, a, b, tag in (("H175.1", p.tertile.to_numpy(), 0, 2, 1), ("H175.2_ER", p.er.to_numpy(), 1, 0, 2)):
        con = A.Contrast({KEY: C}, {KEY: lab}, a, b, 1775 + tag, {KEY: 0})
        res, _ = A.power_row(name, con, args.n_perm)
        res.update(cohort=KEY, hypothesis=name, null_sd_x2_49=2.49 * res["null_sd"], gate_ref_pass=gates[0]["gate_ref_pass"],
                   gate_pair_pass=gates[0]["gate_pair_pass"], n_arm_a=con.n[KEY][0], n_arm_b=con.n[KEY][1],
                   tcga_reference_D=TCGA_D_H1751 if name == "H175.1" else np.nan,
                   estimable=bool(res["mde"] < TCGA_D_H1751) if name == "H175.1" else bool(res["mde"] <= A.ER_STOP))
        rows.append(res)
        arms.append(dict(hypothesis=name, n_arm_a=con.n[KEY][0], n_arm_b=con.n[KEY][1], n_pairs=len(lab), n_missing_label=int((lab < 0).sum())))
        P.log(name, {k: res[k] for k in ("null_mean", "null_sd", "null_q95", "mde", "n_arm_a", "n_arm_b", "estimable")})
    pd.DataFrame(rows).to_csv(OUT / "power.csv", index=False)
    pd.DataFrame(arms).to_csv(OUT / "power_arms.csv", index=False)
    pd.DataFrame(gates).to_csv(OUT / "gates.csv", index=False)
    P.log("done")


if __name__ == "__main__":
    main()
