"""A280 with the calibrated percentile test of A278 (one of the Q36 options): the clock's percentile among its 200
companions, read against the distribution of that percentile over placebo draws (donors permuted within compartment x
stratum; states unchanged).

Registered in paper/docs/registered_experiments.md (A280). `--placebo-only` computes the null draws alone (no statistic
on the real pairing: a code check). The real pairing is read only with `--q36` (after the user's answer to Q36).

    python calibrated_a280.py [--placebo-only | --q36] [--draws 1000] [--jobs 6]
Outputs (../output/calibrated/): null.npz, tests.csv (with --q36), per_atlas.csv (with --q36).
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


A280 = _load("a280c_a280", HERE.parent / "score_a280.py")  # A274 with the A280 atlases and seed
C276 = _load("a280c_a276", PAPER / "exploratory/A276_pairwise_calibration/scripts/calibrate_a276.py")
A274, A272, A150 = A280.A274, A280.A274.A272, A280.A274.A150
import cloudpickle  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

for _m in (A280, C276, A274, A272, A150):
    cloudpickle.register_pickle_by_value(_m)

FAMILIES = {"matched": slice(1, 201), "responsive": slice(201, 401)}


def load_all():
    A150.ATLASES.update(A274.ALL)
    data = {}
    for n in A274.ALL:
        try:
            A = A150.load_atlas(n, None, 0)
        except Exception as e:  # noqa: BLE001  (not built: named, left out)
            C276.log(n, "left out:", type(e).__name__)
            continue
        if A["rep"]["rows_with_count_diff"] != 0:
            C276.log(n, "left out: twin rule")
            continue
        pb = A["pb"].reset_index(drop=True)
        A["pb"] = pb
        wfs = [cs.gtex()]
        for kind in A272.KINDS:
            wfs += A274.pseudo_sets(A, list(A274.ALL).index(n), 200, kind)
        S = np.array([A272.atlas_states(pb[A272.META + list(w.genes)], w)[0] for w in wfs])
        key = pb[["donor", "compartment", "stratum"]].copy()
        pairs = A272.pair_index(key)
        meta = pd.DataFrame([dict(atlas=n, a=p["a"], b=p["b"], family_a=A150.FAM.family(p["a"]),
                                  family_b=A150.FAM.family(p["b"])) for p in pairs])
        cortex = n in A274.CORTEX
        meta["klass"] = [("CX_diff" if fa != fb else "CX_same") if cortex else A272.klass(A272.group(fa), A272.group(fb))
                         for fa, fb in zip(meta.family_a, meta.family_b)]
        data[n] = (key, S, meta)
        C276.log(n, "rows", len(key), "pairs", len(pairs))
    cnt = pd.concat([m for _, _, m in data.values()]).groupby(["atlas", "klass"]).size().unstack(fill_value=0)
    col = lambda c: cnt[c] if c in cnt else pd.Series(0, index=cnt.index)  # noqa: E731
    body = [n for n in data if n in A274.BODY and col("SV_X")[n] >= A274.MIN_PAIRS and col("LY_X")[n] >= A274.MIN_PAIRS]
    cx = [n for n in data if n in A274.CORTEX and col("CX_diff")[n] >= A274.MIN_PAIRS]
    return data, dict(body=body, cx=cx)


def statistics(cm, groups):
    out = {}
    if groups["body"]:
        out["H280.1 D"] = C276.mean_over(cm, groups["body"], "SV_X") - C276.mean_over(cm, groups["body"], "LY_X")
    if groups["cx"]:
        out["H280.2 T_CX"] = C276.mean_over(cm, groups["cx"], "CX_diff")
    return out


def draw(seed, data, groups):
    warnings.filterwarnings("ignore")
    rng = None if seed is None else np.random.default_rng(seed)
    rows, metas = [], []
    for n, (key, S, meta) in data.items():
        rows.append(C276.pair_abs(S, A272.pair_index(key if rng is None else A272.placebo(key, rng))))
        metas.append(meta)
    cm = C276.class_medians(np.concatenate(rows, 1), pd.concat(metas, ignore_index=True))
    return statistics(cm, groups), cm


def pi(v, sl):
    v = np.atleast_2d(v)
    comp = v[:, sl]
    return np.nansum(comp < v[:, [0]], 1) / np.isfinite(comp).sum(1)


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--draws", type=int, default=1000)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--placebo-only", action="store_true")
    g.add_argument("--q36", action="store_true", help="read the real pairing (only after the user's answer to Q36)")
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(A274.AMENDMENT) / "calibrated"
    out.mkdir(parents=True, exist_ok=True)
    data, groups = load_all()
    C276.log("groups", groups)
    res = Parallel(n_jobs=args.jobs)(delayed(draw)([A274.SEED, 9, d], data, groups) for d in range(args.draws))
    null = {t: np.array([r[0][t] for r in res]) for t in res[0][0]}
    np.savez_compressed(out / "null.npz", **null)
    for t, v in null.items():
        for kind, sl in FAMILIES.items():
            p = pi(v, sl)
            C276.log(t, kind, "null pi: mean", round(float(p.mean()), 3), "share pi = 1", round(float(np.mean(p == 1)), 4))
    if not args.q36:
        return
    obs, cmo = draw(None, data, groups)
    tests = []
    for t, v in obs.items():
        row = dict(test=t, clock=float(v[0]))
        for kind, sl in FAMILIES.items():
            o, nl = pi(v, sl)[0], pi(null[t], sl)
            row.update({f"{kind}_median": float(np.nanmedian(v[sl])), f"{kind}_pi": float(o),
                        f"{kind}_p_cal": float((1 + np.sum(nl >= o)) / (1 + len(nl)))})
        row["passes"] = row["matched_p_cal"] < 0.05 and row["responsive_p_cal"] < 0.05
        tests.append(row)
    pd.DataFrame(tests).to_csv(out / "tests.csv", index=False)
    print(pd.DataFrame(tests).round(4).to_string())


if __name__ == "__main__":
    main()
