"""A276: calibration of the anchor-free pairwise test (A272, A274) over many placebo draws.

Registered in paper/docs/registered_experiments.md (A276). Atlases, rows, pairs, families, statistics and companion sets
exactly as scored by A272 and A274. States are computed once per row and gene set (row-wise, unchanged by relabelling
donors); each placebo draw permutes the donor labels within compartment x stratum (`score_a272.placebo`), rebuilds the
pairs and recomputes every statistic for the clock and the 2 x 200 companion sets.

  V276.1  rejection rate (clock above the companion q95) <= 0.075 for each of the ten statistic x family combinations
  V276.2  mean rejection rate over the ten <= 0.065
Draw 0 (no permutation) reproduces the scored tests.csv of A272 and A274 (code check, `repro.csv`).

    python calibrate_a276.py [--draws 500] [--jobs 6]
Outputs (../output/): repro.csv, draws.csv (one row per draw and combination), tests.csv, per_atlas.csv, pair_null.csv.
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

spec = importlib.util.spec_from_file_location("a276_a274", PAPER / "exploratory/A274_celltype_map/scripts/score_a274.py")
A274 = importlib.util.module_from_spec(spec)
sys.modules["a276_a274"] = A274
spec.loader.exec_module(A274)
A272, A150 = A274.A272, A274.A150
import cloudpickle  # noqa: E402

for _m in (A274, A272, A150):
    cloudpickle.register_pickle_by_value(_m)

AMENDMENT = "A276_pairwise_calibration"
SEED = 276
KINDS = A272.KINDS
A272_ATLASES = list(A272.ATLASES)
A274_ATLASES = ["organ_optic_nerve", "organ_eye_anterior", "kidney_sc", "kong_ileum", "kong_colon",
                "psychad_hbcc", "seaad_dlpfc", "seaad_mtg"]  # eligible by the A274 design addendum
log = A272.log


def load(study, n, n_pseudo):
    """Atlas rows and the state matrix S (1 + 2 n_pseudo sets x rows): clock, then matched, then responsive."""
    A = A150.load_atlas(n, None, 0)
    assert study == "A272" or A["rep"]["rows_with_count_diff"] == 0, n  # twin rule of A274 (A272 had none)
    A["pb"] = A["pb"].reset_index(drop=True)
    wfs = [cs.gtex()]
    for kind in KINDS:
        wfs += (A272.pseudo_sets(A, A272_ATLASES.index(n), n_pseudo, kind) if study == "A272"
                else A274.pseudo_sets(A, list(A274.ALL).index(n), n_pseudo, kind))
    S = np.array([A272.atlas_states(A["pb"][A272.META + list(w.genes)], w)[0] for w in wfs])
    return A["pb"][["donor", "compartment", "stratum"]].copy(), S


def klass(study, n, fa, fb):
    if study == "A274" and n in A274.CORTEX:
        return "CX_diff" if fa != fb else "CX_same"
    return A272.klass(A272.group(fa), A272.group(fb))


def pair_abs(S, pairs):
    """|r| (sets x pairs), the A272 statistic without Poisson correction, vectorised over sets."""
    out = np.full((len(S), len(pairs)), np.nan)
    for k, pr in enumerate(pairs):
        O = np.eye(pr["s"].max() + 1)[pr["s"]]
        C = np.eye(len(O)) - O @ (O / O.sum(0)).T
        a, b = S[:, pr["ia"]] @ C, S[:, pr["ib"]] @ C
        den = np.sqrt(np.sum(np.abs(a) ** 2, 1) * np.sum(np.abs(b) ** 2, 1))
        with np.errstate(invalid="ignore", divide="ignore"):
            out[:, k] = np.abs(np.sum(a * np.conj(b), 1)) / den
    return out


def class_medians(absr, meta):
    """{(atlas, class): vector over sets of the median |r| over the class's pairs}."""
    return {key: np.nanmedian(absr[:, idx], 1) for key, idx in meta.groupby(["atlas", "klass"]).indices.items()}


def mean_over(cm, atlases, cls):
    v = [cm[(n, cls)] for n in atlases if (n, cls) in cm]
    return np.mean(v, 0)


def statistics(cm, groups):
    """The ten registered statistics (vectors over sets)."""
    g2, g4 = groups["A272"], groups["A274"]
    return {"H272.1 T_SV": mean_over(cm["A272"], g2["sv"], "SV_X"),
            "H272.2 D": mean_over(cm["A272"], g2["ly"], "SV_X") - mean_over(cm["A272"], g2["ly"], "LY_X"),
            "H274.1 D": mean_over(cm["A274"], g4["both"], "SV_X") - mean_over(cm["A274"], g4["both"], "LY_X"),
            "H274.2 T_SV": mean_over(cm["A274"], g4["sv"], "SV_X"),
            "H274.3 T_CX": mean_over(cm["A274"], g4["cx"], "CX_diff")}


def one_draw(d, data, groups, n_pseudo):
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng([SEED, d])
    cm, absr_all = {}, {}
    for study, atl in data.items():
        rows, metas = [], []
        for n, (key, S, meta) in atl.items():
            pb = key if d == 0 else A272.placebo(key, rng)
            pairs = A272.pair_index(pb)
            assert [(p["a"], p["b"]) for p in pairs] == list(zip(meta.a, meta.b)), (n, "pairs changed")
            rows.append(pair_abs(S, pairs))
            metas.append(meta)
        absr = np.concatenate(rows, 1)
        absr_all[study] = absr
        cm[study] = class_medians(absr, pd.concat(metas, ignore_index=True))
    st = statistics(cm, groups)
    sl = {"matched": slice(1, 1 + n_pseudo), "responsive": slice(1 + n_pseudo, 1 + 2 * n_pseudo)}
    out = []
    for h, v in st.items():
        for kind, s in sl.items():
            c, comp = v[0], v[s]
            out.append(dict(draw=d, test=h, family=kind, clock=c, q95=float(np.nanquantile(comp, 0.95)),
                            median=float(np.nanmedian(comp)), percentile=float(np.mean(comp[np.isfinite(comp)] < c)),
                            reject=bool(c > np.nanquantile(comp, 0.95))))
    pa = []
    for study in cm:
        for (n, cls), v in cm[study].items():
            for kind, s in sl.items():
                pa.append(dict(draw=d, study=study, atlas=n, klass=cls, family=kind, clock=v[0],
                               reject=bool(v[0] > np.nanquantile(v[s], 0.95))))
    pn = {study: np.stack([a[0], np.nanmean(a[sl["matched"]], 0), np.nanmean(a[sl["responsive"]], 0)])
          for study, a in absr_all.items()}
    return out, pa, pn


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--draws", type=int, default=500)
    ap.add_argument("--n-pseudo", type=int, default=200)
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(AMENDMENT)
    out.mkdir(parents=True, exist_ok=True)
    A150.ATLASES.update(A274.ALL)

    data, groups = {}, {}
    for study, names in (("A272", A272_ATLASES), ("A274", A274_ATLASES)):
        data[study] = {}
        for n in names:
            key, S = load(study, n, args.n_pseudo)
            pairs = A272.pair_index(key)
            meta = pd.DataFrame([dict(atlas=n, a=p["a"], b=p["b"], family_a=A150.FAM.family(p["a"]),
                                      family_b=A150.FAM.family(p["b"])) for p in pairs])
            meta["klass"] = [klass(study, n, fa, fb) for fa, fb in zip(meta.family_a, meta.family_b)]
            data[study][n] = (key, S, meta)
            log(study, n, "rows", len(key), "pairs", len(pairs), "sets", len(S))
        cnt = pd.concat([m for _, _, m in data[study].values()]).groupby(["atlas", "klass"]).size().unstack(fill_value=0)
        col = lambda c: cnt[c] if c in cnt else pd.Series(0, index=cnt.index)  # noqa: E731
        if study == "A272":
            groups[study] = dict(sv=names, ly=[n for n in names if col("LY_X")[n] > 0])
        else:
            body = [n for n in names if n in A274.BODY]
            sv = [n for n in body if col("SV_X")[n] >= A274.MIN_PAIRS]
            groups[study] = dict(sv=sv, both=[n for n in sv if col("LY_X")[n] >= A274.MIN_PAIRS],
                                 cx=[n for n in names if n in A274.CORTEX and col("CX_diff")[n] >= A274.MIN_PAIRS])
        log(study, groups[study])

    # draw 0: reproduction of the scored tests
    t0, _, _ = one_draw(0, data, groups, args.n_pseudo)
    t0 = pd.DataFrame(t0)
    scored = pd.concat([pd.read_csv(exploratory_dir("A272_backbone_pairwise") / "tests.csv"),
                        pd.read_csv(exploratory_dir("A274_celltype_map") / "tests.csv")])
    scored = scored[scored.family.isin(list(KINDS))][["test", "family", "clock", "pseudo_q95"]]
    rep = t0.merge(scored, on=["test", "family"], suffixes=("", "_scored"))
    rep["clock_diff"], rep["q95_diff"] = (rep.clock - rep.clock_scored).abs(), (rep.q95 - rep.pseudo_q95).abs()
    rep.to_csv(out / "repro.csv", index=False)
    print(rep[["test", "family", "clock", "clock_scored", "q95", "pseudo_q95"]].to_string())
    if not (len(rep) == 10 and rep.clock_diff.max() < 1e-9 and rep.q95_diff.max() < 1e-9):
        log("draw 0 does not reproduce the scored tests: stop")
        sys.exit(1)

    res = Parallel(n_jobs=args.jobs, verbose=0)(delayed(one_draw)(d, data, groups, args.n_pseudo)
                                                for d in range(1, args.draws + 1))
    D = pd.DataFrame([r for o, _, _ in res for r in o])
    D.to_csv(out / "draws.csv", index=False)
    T = D.groupby(["test", "family"], sort=False).agg(draws=("reject", "size"), rejection_rate=("reject", "mean"),
                                                      mean_percentile=("percentile", "mean"),
                                                      clock_null_mean=("clock", "mean"), q95_mean=("q95", "mean")).reset_index()
    T["V276.1_pass"] = T.rejection_rate <= 0.075
    T.to_csv(out / "tests.csv", index=False)
    PA = pd.DataFrame([r for _, p, _ in res for r in p])
    PA.groupby(["study", "atlas", "klass", "family"]).agg(draws=("reject", "size"), rejection_rate=("reject", "mean"),
                                                          clock_null_mean=("clock", "mean")).reset_index().to_csv(
        out / "per_atlas.csv", index=False)
    pn = []
    for study, atl in data.items():
        meta = pd.concat([m for _, _, m in atl.values()], ignore_index=True)
        m = np.mean([p[study] for _, _, p in res], 0)
        pn.append(meta.assign(study=study, clock_null=m[0], matched_null=m[1], responsive_null=m[2]))
    pd.concat(pn).to_csv(out / "pair_null.csv", index=False)
    print(T.round(3).to_string())
    print("V276.1", "PASS" if T["V276.1_pass"].all() else "FAIL", "| V276.2 mean rate", round(T.rejection_rate.mean(), 4),
          "PASS" if T.rejection_rate.mean() <= 0.065 else "FAIL")
    log("done")


if __name__ == "__main__":
    main()
