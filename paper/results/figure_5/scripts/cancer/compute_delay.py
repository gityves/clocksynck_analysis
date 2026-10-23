"""A97 H97.1: phase of malignant cells vs normal epithelium relative to the patient's stromal host (single cell).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A97). Cohort, states and host as the single-cell part of
compute_decomposition.py (patients with malignant and adjacent-normal epithelium and a stromal host; one malignant and
one normal pseudobulk per patient; host = mean stromal state per patient; one projector; states centred within
study x stratum). c_m, c_n = Poisson-corrected complex correlation of the malignant / normal epithelium state with
the host, across patients; phase = arg c in hours (negative = delayed). Estimand delta = wrap(arg c_m - arg c_n).
H97.1: pooled delta < 0. Null: within-patient swaps of the malignant and normal state (2,000), one-sided p; interval:
patient bootstrap (2,000). Uninformative if |c_m| is not above the host-permuted-across-patients null (1,000): lower
bootstrap bound of |c_m| <= null 95th percentile. Per-study phases are descriptive (>= 5 patients).

    python compute_delay.py [--subsample N] [--seed S]
Outputs (paper/results/figure_5/output/cancer/): a97_delay.csv, a97_delay_per_study.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a97_delay.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync.paired import ccorr, projector, states  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer"
SEED = 97


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_swap, n_boot, n_host = (50, 50, 50) if run["quick"] else (2000, 2000, 1000)
    rng = np.random.default_rng(SEED + args.seed)
    pb, _ = sh.load_sc(args.subsample, args.seed)
    pts = sh.sc_builder().paired_epithelium(pb)
    cls = pb.compartment.str.split(":").str[1]
    ok = pb.donor.isin(pts)
    one = lambda x: x[x.donor.isin(pts)].groupby("donor").head(1).set_index("donor").loc[pts]
    mal, nrm = one(pb[pb.compartment == sh.MAL]), one(pb[pb.compartment == sh.NORM])
    st = pb[ok & cls.isin(sh.STROMA)]
    L = projector(pd.concat([mal, nrm, st]), sh.WF, counts=True)

    def state(x):
        b, v = states(x, L, sh.WF, depth=x.depth.to_numpy(float))
        return b - pd.DataFrame(b).groupby(x.stratum.to_numpy()).transform("mean").to_numpy(), v

    bs, vs = state(st)
    host = pd.DataFrame(bs).groupby(st.donor.to_numpy()).mean().loc[pts].to_numpy()
    k = st.groupby("donor").size().loc[pts].to_numpy()
    hv = pd.Series(vs).groupby(st.donor.to_numpy()).sum().loc[pts].to_numpy() / k ** 2
    bm, vm = state(mal)
    bn, vn = state(nrm)
    n, study = len(pts), mal.study.to_numpy()

    def cc(idx, swap=None):
        if swap is None:
            a, va, b, vb = bm[idx], vm[idx], bn[idx], vn[idx]
        else:
            s = swap[idx]
            a, va = np.where(s[:, None], bn[idx], bm[idx]), np.where(s, vn[idx], vm[idx])
            b, vb = np.where(s[:, None], bm[idx], bn[idx]), np.where(s, vm[idx], vn[idx])
        return ccorr(a, va, host[idx], hv[idx]), ccorr(b, vb, host[idx], hv[idx])

    def stat(idx, swap=None):
        cm, cn = cc(idx, swap)
        return dict(phase_m=np.angle(cm) * sh.HOURS, phase_n=np.angle(cn) * sh.HOURS,
                    delta=float(sh.wrap_h((np.angle(cm) - np.angle(cn)) * sh.HOURS)), abs_m=abs(cm), abs_n=abs(cn))

    all_idx = np.arange(n)
    obs = stat(all_idx)
    null = np.array([stat(all_idx, rng.random(n) < 0.5)["delta"] for _ in range(n_swap)])
    boot = pd.DataFrame([stat(rng.integers(0, n, n)) for _ in range(n_boot)])
    hq95 = float(np.quantile([abs(ccorr(bm, vm, host[p], hv[p])) for p in (rng.permutation(n) for _ in range(n_host))], 0.95))
    res = pd.DataFrame({k_: dict(estimate=v, lo=boot[k_].quantile(0.025), hi=boot[k_].quantile(0.975))
                        for k_, v in obs.items()}).T.assign(patients=n, studies=len(set(study)))
    res.loc["delta", "p_one_lower"] = float((1 + np.sum(null <= obs["delta"])) / (1 + n_swap))
    res.loc["delta", "null_mean"] = float(null.mean())
    res.loc["abs_m", "host_null_q95"] = hq95
    res.loc["abs_m", "informative"] = bool(res.loc["abs_m", "lo"] > hq95)
    res.to_csv(out / "a97_delay.csv")
    per = []
    for s in sorted(set(study)):
        i = np.flatnonzero(study == s)
        if len(i) >= 5:
            cm, cn = cc(i)
            per.append(dict(study=s, patients=len(i), phase_m=np.angle(cm) * sh.HOURS, phase_n=np.angle(cn) * sh.HOURS,
                            delta=float(sh.wrap_h((np.angle(cm) - np.angle(cn)) * sh.HOURS)), abs_m=abs(cm), abs_n=abs(cn)))
    pd.DataFrame(per).to_csv(out / "a97_delay_per_study.csv", index=False)
    print(res.round(3).to_string())
    print(pd.DataFrame(per).round(3).to_string())


if __name__ == "__main__":
    main()
