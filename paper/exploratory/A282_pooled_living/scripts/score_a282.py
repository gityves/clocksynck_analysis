"""A282: the cell-type ranking of A272 and A274 pooled over the thirteen atlases, living against post-mortem donors.

Registered in paper/docs/registered_experiments.md (A282; descriptive, read with the A283 form). Per arm (donors of
that source class only; unknown and recall-based classes left out): pairs re-indexed on the arm's rows with the A272
pair rule, statistics as A272 / A274 (T_SV, T_LY, D per atlas; pooled means over atlases with >= 3 (SV, X) and >= 3
(LY, X) pairs of the arm; T_CX in cortex), each with the clock's pi among its 200 companions per family and p_cal against
5000 placebo draws of the same arm (seeds [282, arm, d]).

`--design` writes arms.csv only (no statistic). Scoring needs A283's calibration (V283) to have passed.

    python score_a282.py [--design] [--draws 5000] [--jobs 6]
Outputs (../output/): arms.csv (donors and pairs per atlas and arm), tests.csv, per_atlas.csv, family_matrix.csv.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER / "exploratory/A276_pairwise_calibration/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import calibrate_a276 as C  # noqa: E402
import corrected_a277 as A277  # noqa: E402
from calibrated_a278 import FAMILIES, p_cal, pi  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

AMENDMENT = "A282_pooled_living"
SEED = 282
ARMS = {"living": 1, "post_mortem": 2}
MIN_PAIRS = 3
SOURCE = HERE.parents[1] / "meta/living_postmortem.csv"
V283 = PAPER / "exploratory/A276_pairwise_calibration/output/a283/calibration.csv"
ALIAS = {"heart": "heart_global"}  # atlas key in the A150 loader -> name in the source table


def source_classes():
    s = pd.read_csv(SOURCE)
    s.loc[s.evidence.str.contains("not re-verified", na=False), "source_class"] = "unknown"  # recall only: left out
    return {(a, d): c for a, d, c in zip(s.atlas, s.donor.astype(str), s.source_class)}


def arm_data(data, cls, arm):
    """Rows of the arm only; pairs and class labels re-derived on them."""
    out, rows = {}, []
    for study, atl in data.items():
        out[study] = {}
        for n, (key, S, _) in atl.items():
            lab = np.array([cls.get((ALIAS.get(n, n), str(d)), "missing") for d in key.donor])
            m = lab == arm
            k = key[m].reset_index(drop=True)
            pairs = C.A272.pair_index(k)
            meta = pd.DataFrame([dict(atlas=n, a=p["a"], b=p["b"], family_a=C.A150.FAM.family(p["a"]),
                                      family_b=C.A150.FAM.family(p["b"])) for p in pairs])
            if len(meta):
                meta["klass"] = [C.klass(study, n, fa, fb) for fa, fb in zip(meta.family_a, meta.family_b)]
                out[study][n] = (k, S[:, m], meta)
            cnt = meta.klass.value_counts() if len(meta) else pd.Series(dtype=int)
            rows.append(dict(arm=arm, study=study, atlas=n, donors=k.donor.nunique(), rows=len(k),
                             missing_class=int((lab == "missing").sum()), pairs=len(meta),
                             **{c: int(cnt.get(c, 0)) for c in ("SV_X", "LY_X", "CX_diff")}))
    return out, pd.DataFrame(rows)


def draw(seed, data):
    """{(atlas, class): vector over sets} on the real pairing (seed None) or a placebo draw."""
    warnings.filterwarnings("ignore")
    rng = None if seed is None else np.random.default_rng(seed)
    cm = {}
    for atl in data.values():
        for n, (key, S, meta) in atl.items():
            pairs = C.A272.pair_index(key if rng is None else C.A272.placebo(key, rng))
            assert [(p["a"], p["b"]) for p in pairs] == list(zip(meta.a, meta.b)), (n, "pairs changed")
            cm.update({(n, k): v for (_, k), v in C.class_medians(C.pair_abs(S, pairs), meta).items()})
    return cm


def statistics(cm, sv, both, cx):
    out = {}
    if sv:
        out["T_SV"] = np.mean([cm[(n, "SV_X")] for n in sv], 0)
    if both:
        out["D"] = np.mean([cm[(n, "SV_X")] - cm[(n, "LY_X")] for n in both], 0)
        out["T_LY"] = np.mean([cm[(n, "LY_X")] for n in both], 0)
    if cx:
        out["T_CX"] = np.mean([cm[(n, "CX_diff")] for n in cx], 0)
    return out


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--draws", type=int, default=5000)
    ap.add_argument("--design", action="store_true")
    args = ap.parse_args()
    cli.setup(args)
    if not args.design and not (V283.exists() and pd.read_csv(V283).V283_pass.all()):
        sys.exit("A282: scoring waits for V283 (A283)")
    out = exploratory_dir(AMENDMENT)
    out.mkdir(parents=True, exist_ok=True)
    C.A150.ATLASES.update(C.A274.ALL)
    data, _ = A277.load_all(200)
    cls = source_classes()
    tests, per_atlas, fam, arms = [], [], [], []
    for arm, code in ARMS.items():
        d, a = arm_data(data, cls, arm)
        arms.append(a)
        C.log(arm, a[a.pairs > 0][["atlas", "donors", "pairs", "SV_X", "LY_X", "CX_diff"]].to_dict("records"))
        if args.design:
            continue
        body = a[~a.atlas.isin(list(C.A274.CORTEX))]
        sv = list(body.atlas[body.SV_X >= MIN_PAIRS])
        both = list(body.atlas[(body.SV_X >= MIN_PAIRS) & (body.LY_X >= MIN_PAIRS)])
        cx = list(a.atlas[a.atlas.isin(list(C.A274.CORTEX)) & (a.CX_diff >= MIN_PAIRS)])
        obs = draw(None, d)
        res = Parallel(n_jobs=args.jobs)(delayed(draw)([SEED, code, k], d) for k in range(args.draws))
        null = {k: np.array([r[k] for r in res]) for k in obs}
        C.log(arm, "draws done")
        st, nst = statistics(obs, sv, both, cx), {}
        for k in st:
            nst[k] = np.array([statistics(r, sv, both, cx)[k] for r in res])
        for k, v in st.items():
            row = dict(arm=arm, statistic=k, atlases=";".join(both if k in ("D", "T_LY") else sv if k == "T_SV" else cx),
                       clock=float(v[0]))
            for kind, sl in FAMILIES.items():
                o = pi(v, sl)[0]
                row.update({f"{kind}_median": float(np.nanmedian(v[sl])), f"{kind}_pi": float(o),
                            f"{kind}_p_cal": float(p_cal(o, pi(nst[k], sl))[0])})
            tests.append(row)
        for (n, k), v in obs.items():
            row = dict(arm=arm, atlas=n, klass=k, clock=float(v[0]))
            for kind, sl in FAMILIES.items():
                o = pi(v, sl)[0]
                row.update({f"{kind}_median": float(np.nanmedian(v[sl])), f"{kind}_pi": float(o),
                            f"{kind}_p_cal": float(p_cal(o, pi(null[(n, k)], sl))[0])})
            per_atlas.append(row)
        for atl in d.values():  # family matrix: median over pairs of each family pair
            for n, (key, S, meta) in atl.items():
                r = C.pair_abs(S, C.A272.pair_index(key))
                for i, m in meta.iterrows():
                    fam.append(dict(arm=arm, part="cortex" if n in C.A274.CORTEX else "body",
                                    f1=min(m.family_a, m.family_b), f2=max(m.family_a, m.family_b), atlas=n, clock=r[0, i],
                                    matched=np.nanmedian(r[FAMILIES["matched"], i]),
                                    responsive=np.nanmedian(r[FAMILIES["responsive"], i])))
    pd.concat(arms).to_csv(out / "arms.csv", index=False)
    if args.design:
        print(pd.concat(arms).to_string())
        return
    tests = pd.DataFrame(tests)
    tests.to_csv(out / "tests.csv", index=False)
    pd.DataFrame(per_atlas).to_csv(out / "per_atlas.csv", index=False)
    fm = pd.DataFrame(fam).groupby(["arm", "part", "f1", "f2"]).agg(
        pairs=("clock", "size"), atlases=("atlas", "nunique"), clock=("clock", "median"), matched=("matched", "median"),
        responsive=("responsive", "median")).reset_index()
    fm.to_csv(out / "family_matrix.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.concat(arms).to_string())
    print(tests.round(4).to_string())
    C.log("done")


if __name__ == "__main__":
    main()
