"""A278: the pairwise test calibrated by construction (after A276, A277).

Registered in paper/docs/registered_experiments.md (A278). For each statistic and family, pi = share of the 200 companion
values below the clock's; the null of pi comes from the 1000 placebo draws of A277 (`output/a277/null_draws.npz`).
p_cal = (1 + #{null pi >= observed pi}) / (1 + draws).

  V278     with N1 as the null, rejection rate (p_cal < 0.05) over the N2 draws <= 0.075 per combination (null only)
  --read   the restated tests H278.1-5 on the real pairing (registered: only after the user has seen the A276/A277
           report), and the same per atlas and pair class

    python calibrated_a278.py [--read] [--jobs 6]
Outputs (../output/a278/): calibration.csv; with --read: tests.csv, per_atlas.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import calibrate_a276 as C  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

FAMILIES = {"matched": slice(1, 201), "responsive": slice(201, 401)}
TESTS = {"H278.1": ("H272.1 T_SV", "both"), "H278.2": ("H272.2 D", "both"), "H278.3": ("H274.1 D", "both"),
         "H278.4": ("H274.2 T_SV", "matched"), "H278.5": ("H274.3 T_CX", "both")}


def pi(v, sl):
    """Share of companion values below the clock's (rows = draws or a single vector)."""
    v = np.atleast_2d(v)
    comp = v[:, sl]
    return np.nansum(comp < v[:, [0]], 1) / np.isfinite(comp).sum(1)


def p_cal(obs_pi, null_pi):
    return (1 + np.sum(null_pi[None, :] >= np.atleast_1d(obs_pi)[:, None], 1)) / (1 + len(null_pi))


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--read", action="store_true")
    args = ap.parse_args()
    cli.setup(args)
    src = exploratory_dir(C.AMENDMENT) / "a277"
    out = exploratory_dir(C.AMENDMENT) / "a278"
    out.mkdir(parents=True, exist_ok=True)
    z = np.load(src / "null_draws.npz")
    stats = sorted({k[3:] for k in z.files})

    cal = []
    for t in stats:
        for kind, sl in FAMILIES.items():
            n1, n2 = pi(z["N1 " + t], sl), pi(z["N2 " + t], sl)
            rate = float(np.mean(p_cal(n2, n1) < 0.05))
            cal.append(dict(test=t, family=kind, draws=len(n2), rejection_rate=rate,
                            null_rate_pi_1=float(np.mean(np.concatenate([n1, n2]) == 1.0))))
    cal = pd.DataFrame(cal)
    cal["V278_pass"] = cal.rejection_rate <= 0.075
    cal.to_csv(out / "calibration.csv", index=False)
    print(cal.round(3).to_string())
    print("V278", "PASS" if cal.V278_pass.all() else "FAIL")
    if not args.read:
        return

    import corrected_a277 as A277  # real pairing: loads every atlas
    C.A150.ATLASES.update(C.A274.ALL)
    data, groups = A277.load_all(200)
    obs, cmo = A277.draw(None, data, groups)
    tests = []
    for h, (t, rule) in TESTS.items():
        row = dict(test=h, restates=t, clock=float(obs[t][0]))
        for kind, sl in FAMILIES.items():
            null = np.concatenate([pi(z["N1 " + t], sl), pi(z["N2 " + t], sl)])
            o = pi(obs[t], sl)[0]
            row.update({f"{kind}_pi": float(o), f"{kind}_p_cal": float(p_cal(o, null)[0])})
        row["passes"] = row["matched_p_cal"] < 0.05 and (rule == "matched" or row["responsive_p_cal"] < 0.05)
        tests.append(row)
    pd.DataFrame(tests).to_csv(out / "tests.csv", index=False)
    print(pd.DataFrame(tests).round(4).to_string())
    # per atlas and pair class: null from fresh placebo draws of the class medians
    _, cm1 = A277.run([[A277.SEED, 1, d] for d in range(500)], data, groups, args.jobs)
    _, cm2 = A277.run([[A277.SEED, 2, d] for d in range(500)], data, groups, args.jobs)
    pa = []
    for k, v in cmo.items():
        row = dict(study=k[0], atlas=k[1], klass=k[2], clock=float(v[0]))
        for kind, sl in FAMILIES.items():
            null = np.concatenate([pi(cm1[k], sl), pi(cm2[k], sl)])
            o = pi(v, sl)[0]
            row.update({f"{kind}_median": float(np.nanmedian(v[sl])), f"{kind}_pi": float(o),
                        f"{kind}_p_cal": float(p_cal(o, null)[0])})
        pa.append(row)
    pd.DataFrame(pa).to_csv(out / "per_atlas.csv", index=False)


if __name__ == "__main__":
    main()
