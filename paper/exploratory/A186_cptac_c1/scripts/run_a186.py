"""A186: CPTAC-3 replication of the A175 H175.1 proliferation pairing drop (TCGA D = 0.158).

Registered in paper/docs/registered_experiments.md, section "A186". Data: output/cache/<TYPE>.npz (prep_a186.py).
Stages (run with `python -I`):
  gates  : pairs, gate of the normals, G2M tertiles, arm sizes -> gates.csv, power_arms.csv       (labels only)
  power  : MDE of the label-permutation null (10,000, within type) -> power.csv; sets `stop` if MDE > 0.158
  tests  : refuses to run if power.csv says stop. Observed D, p, companions (200 pseudo-clock + 200 responsive sets,
           same labels), per-type D, A185 D_match -> tests.csv, per_type.csv, descriptive.csv, companion_D.csv
    python -I run_a186.py --stage gates|power|tests [--jobs 4] [--n-perm 10000] [--n-sets 200]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a186_lib as L  # noqa: E402
from a186_lib import r75, q, Cancer  # noqa: E402
sys.path.insert(0, str(L.HERE.parent / "A185_c1_spread_matched" / "scripts"))
import run_a185 as r85  # noqa: E402
from common import cli  # noqa: E402

TCGA_D = 0.1577395599995544  # A175 H175.1 (output/tests.csv)
STOP = 0.158
log = r75.log


def build(args):
    st, lab, g2m = {}, {}, {}
    for i, tp in enumerate(L.TYPES):
        st[tp] = L.states(tp, i)
        g2m[tp], lab[tp] = L.labels(tp)
    gated = [t for t in L.TYPES if st[t]["gated"]]
    cs = {t: Cancer(t, st[t]) for t in gated}
    cidx = {t: i for i, t in enumerate(L.TYPES)}
    con = r75.Contrast(cs, {t: lab[t] for t in gated}, 0, 2, 1, cidx)
    return st, lab, g2m, gated, cs, cidx, con


def stage_gates(args, out):
    st, lab, g2m, gated, cs, cidx, con = build(args)
    rows = []
    for t in L.TYPES:
        n = [int((lab[t] == k).sum()) for k in range(3)]
        rows.append(dict(type=t, n_pairs=st[t]["n"], gate_p=st[t]["gate_p"], gated=st[t]["gated"], n_low=n[0], n_mid=n[1], n_high=n[2],
                         in_pool=t in con.use, zero_clock_counts_normal=st[t]["zero_clock_counts_normal"],
                         zero_clock_counts_tumour=st[t]["zero_clock_counts_tumour"]))
    d = pd.DataFrame(rows)
    d.to_csv(out / "gates.csv", index=False)
    print(d.to_string(), flush=True)


def stage_power(args, out):
    stage_gates(args, out)
    st, lab, g2m, gated, cs, cidx, con = build(args)
    r, null = r75.power_row("H186.1", con, args.n_perm)
    r["stop_rule_mde_gt_0.158"] = bool(r["mde"] > STOP)
    r["total_pairs_in_pool"] = sum(st[c]["n"] for c in con.use)
    r["tcga_D"] = TCGA_D
    pd.DataFrame([r]).to_csv(out / "power.csv", index=False)
    log("power", {k: r[k] for k in ("null_mean", "null_sd", "null_q95", "mde", "stop_rule_mde_gt_0.158", "n_pairs_per_arm_per_cancer")})


def stage_tests(args, out):
    pw = pd.read_csv(out / "power.csv").iloc[0]
    if pw.mde > STOP:
        sys.exit(f"STOP: MDE {pw.mde:.4f} > {STOP}; not estimable in CPTAC. Test not run.")
    st, lab, g2m, gated, cs, cidx, con = build(args)
    null = con.null(args.n_perm)
    s = q.summarise(null)
    assert np.isclose(s["mde"], pw.mde, rtol=1e-9), "null differs from power.csv"
    x = null[np.isfinite(null)]
    per, pa, pb, D = con.observed()
    p = float((1 + np.sum(x >= D)) / (len(x) + 1))
    log(f"D {D:.4f} p {p:.4f} MDE {s['mde']:.4f}")
    res = Parallel(n_jobs=args.jobs)(delayed(L.companion_type)(t, cidx[t], args.n_sets, lab[t]) for t in con.use)
    res = dict(zip(con.use, res))
    log("companion states done")
    comp = {}
    for kind in ("matched", "responsive"):
        pa_, pb_ = [], []
        for t in con.use:
            C = [subset_stats(c, lab[t]) for c in res[t][kind]]
            pa_.append([a for a, b in C]); pb_.append([b for a, b in C])
        comp[kind] = (q.pool(np.array(pa_), [con.n[t][0] for t in con.use]) - q.pool(np.array(pb_), [con.n[t][1] for t in con.use]))
    row = dict(hypothesis="H186.1", D=D, pooled_low=pa, pooled_high=pb, p=p, null_mean=s["null_mean"], null_sd=s["null_sd"],
               null_q95=s["null_q95"], mde=s["mde"], n_perm=len(x), types=";".join(con.use),
               n_pairs_low=sum(con.n[t][0] for t in con.use), n_pairs_high=sum(con.n[t][1] for t in con.use), tcga_D=TCGA_D)
    for kind, Dk in comp.items():
        ps = r75.pseudo_summary(Dk, D)
        row |= {f"{kind}_q95": ps["q95"], f"{kind}_median": ps["median"], f"{kind}_percentile_of_D": ps["pct"],
                f"{kind}_D_exceeds_q95": ps["exceeds"], f"{kind}_n_sets": ps["n"]}
    row["clock_specific"] = bool(row["matched_D_exceeds_q95"] and row["responsive_D_exceeds_q95"])
    row["passes_p05"] = bool(p < 0.05 and D > 0)
    pd.DataFrame([row]).to_csv(out / "tests.csv", index=False)
    pd.DataFrame({f"{k}_D": v for k, v in comp.items()}).to_csv(out / "companion_D.csv", index=False)
    # per type -----------------------------------------------------------------------------------------------------
    pt = []
    for t in gated:
        t3 = [float(cs[t].stat((lab[t] == k)[None, :])[0]) for k in range(3)]
        pt.append(dict(type=t, n_pairs=cs[t].n, n_low=int((lab[t] == 0).sum()), n_high=int((lab[t] == 2).sum()), stat_low=t3[0],
                       stat_mid=t3[1], stat_high=t3[2], D_low_minus_high=t3[0] - t3[2], in_pool=t in con.use))
    pd.DataFrame(pt).to_csv(out / "per_type.csv", index=False)
    # D_match (A185 method) ----------------------------------------------------------------------------------------
    d = {t: dict(tert=lab[t]) for t in gated}
    amp = {t: np.abs(st[t]["za"]) for t in gated}
    m = r85.Match(cs, d, cidx, con.use, amp, 1)
    nm = m.null(args.n_perm)
    (Dm, nb, npairs), _ = m.observed()
    xm = nm[np.isfinite(nm)]
    sm = q.summarise(nm)
    desc = [dict(descriptive="D_match (A185 method, tumour amplitude tertiles within type)", value=Dm,
                 p=float((1 + np.sum(xm >= Dm)) / (len(xm) + 1)), null_sd=sm["null_sd"], mde=sm["mde"], n_bins=nb, n_pairs=npairs,
                 note=f"bins dropped {len(m.dropped)}: " + ";".join(f"{c}:{k}({a}/{b})" for c, k, _, a, b in m.dropped)),
            dict(descriptive="share of D removed (D - D_match)/D", value=(D - Dm) / D)]
    pd.DataFrame(desc).to_csv(out / "descriptive.csv", index=False)
    print(pd.DataFrame([row]).T.to_string(), flush=True)
    print(pd.DataFrame(pt).round(4).to_string(), flush=True)
    print(pd.DataFrame(desc).round(4).to_string(), flush=True)
    log("tests stage done")


def subset_stats(C: Cancer, lab: np.ndarray):
    return float(C.stat((lab == 0)[None, :])[0]), float(C.stat((lab == 2)[None, :])[0])


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["gates", "power", "tests"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-sets", type=int, default=200)
    args = ap.parse_args(argv)
    cli.setup(args)
    out = L.HERE / "output" / ("" if args.n_perm == 10_000 and args.n_sets == 200 else "debug")
    out.mkdir(parents=True, exist_ok=True)
    {"gates": stage_gates, "power": stage_power, "tests": stage_tests}[args.stage](args, out)


if __name__ == "__main__":
    main()
