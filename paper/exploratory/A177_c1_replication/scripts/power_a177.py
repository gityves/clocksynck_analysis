"""A177 power step (replication of A175 H175.1; labels only). Per cohort and pooled over both cohorts (pairs-weighted):
pairs, tertile arm sizes, clock gate (estimability precondition) and the permutation null of D = pairing(low G2M tertile) -
pairing(high tertile), MDE = q95 + 0.84 SD (10,000 permutations of the tertile labels). NO observed pairing statistic, stratified
or not, is evaluated in this script (the `Contrast.observed` method of A175 is never called).

Reuses A175 (`run_a175.Contrast`, `Cancer`) and A138 (`q3_lib`) by import; states from `clocksync.paired.cohort_states`.
Microarray adaptation (the only one): `counts=False`, so the projector is unweighted, states are log2 x ln 2 and the Poisson
variances are 0 (`paired.states`); the statistic of `Cancer.stat` then reduces to the plain |complex correlation|. Clock genes
absent from a platform are dropped from the waveform (GSE14520: CIART). Gates, as A138 `tcga_states` on arrays:
  gate_ref   `paired.reference_gate` on the normals (A52 rule, p < 0.05, 2,000 permutations of acrophases);
  gate_pair  model-free tumour-normal cross-covariance gate of `Dataset.gate` (clock_cov of the pair, p_pairing against shuffled
             pairing of the normals, p_acrophase for the clock share against permuted acrophases), written here for arrays.
    python power_a177.py [--n-perm 10000] [--n-gate 2000]
Outputs (output/): power.csv, power_arms.csv, gates.csv"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parent / "A175_c1_pairing_strata" / "scripts"))
import run_a175 as A  # noqa: E402  (also puts A138 scripts and paper on sys.path)
import q3_lib as q  # noqa: E402
from clocksync.paired import cohort_states, reference_gate  # noqa: E402
from clocksync._core.study import plane  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE / "output"
KEYS = ("gse14520", "gse66229")
GATE_SEED = 177
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:6.1f} min]", *a, flush=True)


def load(key):
    t, n = pd.read_parquet(OUT / f"{key}_tumour.parquet"), pd.read_parquet(OUT / f"{key}_normal.parquet")
    p = pd.read_csv(OUT / f"{key}_pairs.csv")
    wf = q.WF.drop([g for g in q.WF.genes if g not in t.columns])
    return t, n, p, wf


def pair_gate(yt, yn, wf, n_perm, rng):
    """Tumour-normal cross-covariance gate (as `Dataset.gate` for two compartments): clock_cov = ||P C P||^2, C = cross-covariance
    of the gene-centred log expression of tumour and own normal (noise independent between samples). p_pairing: normals shuffled
    across patients; p_acrophase: clock share of C against planes of permuted acrophases."""
    yt, yn = yt - yt.mean(1, keepdims=True), yn - yn.mean(1, keepdims=True)
    cov = lambda b: (yt - yt.mean(0)).T @ (b - b.mean(0)) / (len(yt) - 1)
    P = plane(wf)
    cc = lambda C, Pm: np.sum((Pm @ C @ Pm) ** 2)
    C = cov(yn)
    obs, tot = cc(C, P), np.sum(C ** 2)
    null_cl = np.array([cc(cov(yn[rng.permutation(len(yn))]), P) for _ in range(n_perm)])
    null_sh = np.array([cc(C, plane(wf.permuted(rng))) / tot for _ in range(n_perm)])
    p = lambda null, o: float((1 + np.sum(null >= o)) / (1 + len(null)))
    return dict(pair_clock_share=float(obs / tot), pair_share_null_mean=float(null_sh.mean()), p_pairing=p(null_cl, obs), p_acrophase=p(null_sh, obs / tot))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-gate", type=int, default=2000)
    args = ap.parse_args(argv)
    cli.setup(args)
    cs, labels, gates, arms = {}, {}, [], []
    for i, key in enumerate(KEYS):
        t, n, p, wf = load(key)
        g = list(wf.genes)
        bt, vt, bn, vn = cohort_states(t, n, False, wf)
        assert not vt.any() and not vn.any()  # arrays: no Poisson term
        cs[key] = A.Cancer(key, dict(n=len(t), raw=(bt, vt, bn, vn)))
        labels[key] = p.tertile.to_numpy()
        yn, yt = n[g].to_numpy(float) * np.log(2), t[g].to_numpy(float) * np.log(2)
        r = reference_gate(yn, wf, args.n_gate, np.random.default_rng([GATE_SEED, i]))
        pg = pair_gate(yt, yn, wf, args.n_gate, np.random.default_rng([GATE_SEED, 10, i]))
        gates.append(dict(cohort=key, n_pairs=len(t), n_clock_genes=len(g), missing_clock_genes=";".join(x for x in q.WF.genes if x not in g),
                          ref_share=r["share"], ref_share_null_mean=r["share_null_mean"], gate_ref_p=r["p"], gate_ref_pass=bool(r["p"] < 0.05),
                          **pg, gate_pair_pass=bool(pg["p_pairing"] < 0.05 and pg["p_acrophase"] < 0.05), n_perm_gate=args.n_gate))
        log(key, {k: gates[-1][k] for k in ("gate_ref_p", "p_pairing", "p_acrophase")})
    rows = []
    cidx = {k: i for i, k in enumerate(KEYS)}
    for name, ks in (("gse14520", ["gse14520"]), ("gse66229", ["gse66229"]), ("pooled", list(KEYS))):
        con = A.Contrast({k: cs[k] for k in ks}, {k: labels[k] for k in ks}, 0, 2, 177, cidx)
        r, _ = A.power_row(name, con, args.n_perm)
        r["cohort"] = r.pop("hypothesis")
        r["gate_ref_pass"] = all(g["gate_ref_pass"] for g in gates if g["cohort"] in ks)
        rows.append(r)
        arms += [dict(cohort=name, member=k, n_pairs=cs[k].n, n_low=con.n[k][0], n_high=con.n[k][1], in_pool=k in con.use) for k in ks]
        log(name, {k: r[k] for k in ("null_mean", "null_sd", "null_q95", "mde")})
    pd.DataFrame(rows).to_csv(OUT / "power.csv", index=False)
    pd.DataFrame(arms).to_csv(OUT / "power_arms.csv", index=False)
    pd.DataFrame(gates).to_csv(OUT / "gates.csv", index=False)
    log("done")


if __name__ == "__main__":
    main()
