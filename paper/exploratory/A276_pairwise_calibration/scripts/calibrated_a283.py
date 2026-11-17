"""A283: the A278 test with a null of 5000 placebo draws (user's answer to Q36, option b).

Registered in paper/docs/registered_experiments.md (A283). pi = share of the 200 companion values below the clock's;
p_cal = (1 + #{null pi >= observed pi}) / (1 + draws). Null: 5000 placebo draws, seeds [283, d] (A277 `draw`).

  V283     with one half of the draws as the null, rejection rate (p_cal < 0.05) over the other half <= 0.0685 for each
           of the ten combinations, in both directions (null only)
  H283.1-5 only if V283 passes: H272.1, H272.2, H274.1, H274.2 (matched decides), H274.3 on the real pairing, null = all
           5000 draws; descriptive: the same per atlas and pair class

    python calibrated_a283.py [--draws 5000] [--jobs 6]
Outputs (../output/a283/): null_draws.npz, calibration.csv; if V283 passes: tests.csv, per_atlas.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import calibrate_a276 as C  # noqa: E402
import corrected_a277 as A277  # noqa: E402
from calibrated_a278 import FAMILIES, p_cal, pi  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

SEED, CHUNK, BAND = 283, 500, 0.0685
TESTS = {"H283.1": ("H272.1 T_SV", "both"), "H283.2": ("H272.2 D", "both"), "H283.3": ("H274.1 D", "both"),
         "H283.4": ("H274.2 T_SV", "matched"), "H283.5": ("H274.3 T_CX", "both")}


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--draws", type=int, default=5000)
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(C.AMENDMENT) / "a283"
    out.mkdir(parents=True, exist_ok=True)
    C.A150.ATLASES.update(C.A274.ALL)
    data, groups = A277.load_all(200)

    st, cm = [], []
    for a in range(0, args.draws, CHUNK):
        s, c = A277.run([[SEED, d] for d in range(a, min(a + CHUNK, args.draws))], data, groups, args.jobs)
        st.append(s), cm.append(c)
        C.log("draws", min(a + CHUNK, args.draws), "of", args.draws)
    null = {t: np.concatenate([s[t] for s in st]) for t in st[0]}
    cmn = {k: np.concatenate([c[k] for c in cm]) for k in cm[0]}
    np.savez_compressed(out / "null_draws.npz", **null)

    half = args.draws // 2
    cal = []
    for t, v in null.items():
        for kind, sl in FAMILIES.items():
            p = pi(v, sl)
            a, b = p[:half], p[half:]
            cal.append(dict(test=t, family=kind, draws=len(p), rate_null_first=float(np.mean(p_cal(b, a) < 0.05)),
                            rate_null_second=float(np.mean(p_cal(a, b) < 0.05)), null_rate_pi_1=float(np.mean(p == 1.0)),
                            null_pi_mean=float(p.mean())))
    cal = pd.DataFrame(cal)
    cal["V283_pass"] = (cal.rate_null_first <= BAND) & (cal.rate_null_second <= BAND)
    cal.to_csv(out / "calibration.csv", index=False)
    print(cal.round(4).to_string())
    ok = bool(cal.V283_pass.all())
    C.log("V283", "PASS" if ok else "FAIL")
    if not ok or args.draws < 5000:  # a short run is a code check: the real pairing is read only with the registered null
        C.log("restated tests not read")
        return

    obs, cmo = A277.draw(None, data, groups)
    tests = []
    for h, (t, rule) in TESTS.items():
        row = dict(test=h, restates=t, clock=float(obs[t][0]))
        for kind, sl in FAMILIES.items():
            o = pi(obs[t], sl)[0]
            row.update({f"{kind}_median": float(np.nanmedian(obs[t][sl])), f"{kind}_pi": float(o),
                        f"{kind}_p_cal": float(p_cal(o, pi(null[t], sl))[0])})
        row["passes"] = row["matched_p_cal"] < 0.05 and (rule == "matched" or row["responsive_p_cal"] < 0.05)
        tests.append(row)
    tests = pd.DataFrame(tests)
    tests.to_csv(out / "tests.csv", index=False)
    print(tests.round(4).to_string())
    pa = []
    for k, v in cmo.items():
        row = dict(study=k[0], atlas=k[1], klass=k[2], clock=float(v[0]))
        for kind, sl in FAMILIES.items():
            o = pi(v, sl)[0]
            row.update({f"{kind}_median": float(np.nanmedian(v[sl])), f"{kind}_pi": float(o),
                        f"{kind}_p_cal": float(p_cal(o, pi(cmn[k], sl))[0])})
        pa.append(row)
    pd.DataFrame(pa).to_csv(out / "per_atlas.csv", index=False)
    C.log("done")


if __name__ == "__main__":
    main()
