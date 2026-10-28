"""A152 fits: A17 reproduction, observed state couplings, gates, and the within-donor label permutation null (H152.1).

Estimator: `clocksync.Dataset(pb, anchors=A17 learned set).fit(n_boot=0, gate=False)`; the host is fitted on the anchors,
which the lymphocyte rows do not touch, so a label swap only changes the two target rows.
Tables (prepare_lung.py, prepare_part2.py): lung b_/m_ T50, T20, state, cd8_state; heart, breast (+ _blockLV, breast_cd8).
Permutation: within each donor with both states, the two labels swapped with probability 1/2, refit, 1,000 times, rng
[152, table index] (donors with one state keep their label). Statistic of H152.1: mean over heart and breast of c_res - c_circ.

    python fit_a152.py [--jobs 4] [--n-perm 1000]
Outputs (output/): repro_a17.csv, fits_observed.csv, gates.csv, perm_null.csv, perm_summary.csv, h1521.csv, states_summary.csv.
"""
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a152_lib as L  # noqa: E402
from _shared import Target  # noqa: E402
from common import cli  # noqa: E402

# name -> (A17 target, table file, (resident label, circulating label) or None)
S = {"lung_state": ("lung_hlca", "lung_b_state", ("lung:T_resident_like", "lung:T_circulating_like")),
     "lung_m_state": ("lung_hlca", "lung_m_state", ("lung:T_resident_like", "lung:T_circulating_like")),
     "lung_cd8": ("lung_hlca", "lung_b_cd8_state", ("lung:T_CD8_resident_like", "lung:T_CD8_circulating_like")),
     "heart_state": ("heart", "part2_heart", ("heart:lymph_resident_like", "heart:lymph_circulating_like")),
     "heart_blockLV": ("heart", "part2_heart_blockLV", ("heart:lymph_resident_like", "heart:lymph_circulating_like")),
     "breast_state": ("organ_breast", "part2_breast", ("breast:lymph_resident_like", "breast:lymph_circulating_like")),
     "breast_blockLV": ("organ_breast", "part2_breast_blockLV", ("breast:lymph_resident_like", "breast:lymph_circulating_like")),
     "breast_cd8": ("organ_breast", "part2_breast_cd8", ("breast:CD8_resident_like", "breast:CD8_circulating_like"))}
ONE = {"lung_T50_b": ("lung_hlca", "lung_b_T50"), "lung_T20_b": ("lung_hlca", "lung_b_T20"),
       "lung_T50_m": ("lung_hlca", "lung_m_T50"), "lung_T20_m": ("lung_hlca", "lung_m_T20")}
PERM = ["lung_state", "heart_state", "breast_state", "lung_cd8", "breast_cd8"]  # tables with the permutation null


def table(name: str, subset_both: bool = False):
    tgt, f, pair = (S[name] if name in S else (*ONE[name], None))
    T = Target(tgt)
    pb = pd.read_parquet(L.OUT / f"{f}_pb.parquet")
    keep = set(L.compartments_of(tgt)) | set(pair or []) | ({"lung:T"} if pair is None else set())
    pb = pb[pb.donor.isin(T.donors.index) & pb.compartment.isin(keep)]
    if subset_both and pair:
        both = pb[pb.compartment.isin(pair)].groupby("donor").compartment.nunique()
        pb = pb[pb.donor.isin(both.index[both == 2])]
    return T, pb.reset_index(drop=True), pair


def _perm_chunk(name, ks, jobs_seed):
    L_ = L
    from clocksync.config import limit_threads
    limit_threads(1)
    warnings.filterwarnings("ignore")
    T, pb, pair = table(name)
    anc = L.anchors(T.name)
    out = []
    for k in ks:
        q = L.permute_pairs(pb, *pair, np.random.default_rng([L.SEED, jobs_seed, k]))
        t = L.fit_table(q, T.donors, anc, T.name, targets=list(pair))
        out.append((k, float(t.loc[pair[0], "coupling"]), float(t.loc[pair[1], "coupling"])))
    return out


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=1000)
    args = ap.parse_args()
    from clocksync.config import limit_threads
    limit_threads(1)
    warnings.filterwarnings("ignore")
    t00 = time.time()
    # 1. A17 reproduction (point estimates, n_boot = 0) of the three targets, learned rule
    csv = pd.read_csv(L.PAPER / "results/figure_2/output/atlas_fits/a17_healthy_fits.csv")
    rows = []
    for tgt in ("lung_hlca", "heart", "organ_breast"):
        T = Target(tgt)
        t = L.fit_table(T.pb, T.donors, L.anchors(tgt), tgt)
        ref = csv[(csv.target == tgt) & (csv.rule == "learned")].set_index("compartment")
        for c in t.index:
            rows.append(dict(target=tgt, compartment=c, coupling=t.loc[c, "coupling"], a17_csv=ref.coupling.get(c, np.nan),
                             n=t.loc[c, "n"], a17_n=ref.n.get(c, np.nan)))
    rep = pd.DataFrame(rows).assign(diff=lambda d: d.coupling - d.a17_csv)
    rep.to_csv(L.OUT / "repro_a17.csv", index=False)
    print("repro max |diff|", rep["diff"].abs().max(), f"{time.time() - t00:.0f}s", flush=True)
    # 2. observed fits and gates
    obs, gates = [], []
    for i, name in enumerate(list(S) + list(ONE)):
        for both in ((False, True) if name in S and name in ("lung_state", "heart_state", "breast_state") else (False,)):
            T, pb, pair = table(name, both)
            anc = L.anchors(T.name)
            t = L.fit_table(pb, T.donors, anc, T.name).reset_index()
            t.insert(0, "table", name + ("_bothonly" if both else ""))
            obs.append(t)
            if not both:
                g = L.cs.Dataset(pb, anchors=None, donors=T.donors, name=name).gate(n_perm=args.n_perm, seed=L.SEED + i)
                gates.append(dict(table=name, clock_cov=g.clock_cov, clock_share=g.clock_share, p_pairing=g.p_pairing,
                                  p_acrophase=g.p_acrophase, n_donors=g.n_donors, n_compartments=g.n_compartments,
                                  passed=bool(g.passes)))
        print(name, f"{time.time() - t00:.0f}s", flush=True)
    pd.concat(obs).to_csv(L.OUT / "fits_observed.csv", index=False)
    pd.DataFrame(gates).to_csv(L.OUT / "gates.csv", index=False)
    # 3. permutation nulls
    nulls = []
    for j, name in enumerate(PERM):
        chunks = np.array_split(np.arange(args.n_perm), args.jobs)
        res = Parallel(n_jobs=args.jobs)(delayed(_perm_chunk)(name, c.tolist(), list(S).index(name)) for c in chunks)
        d = pd.DataFrame([r for c in res for r in c], columns=["perm", "c_res", "c_circ"]).sort_values("perm")
        d.insert(0, "table", name)
        nulls.append(d)
        print("perm", name, f"{time.time() - t00:.0f}s", flush=True)
    N = pd.concat(nulls).assign(diff=lambda d: d.c_res - d.c_circ)
    N.to_csv(L.OUT / "perm_null.csv", index=False)
    O = pd.concat(obs).set_index(["table", "compartment"])
    summ, dobs = [], {}
    for name in PERM:
        pair = S[name][2]
        d = O.loc[(name, pair[0]), "coupling"] - O.loc[(name, pair[1]), "coupling"]
        nl = N[N.table == name]["diff"].to_numpy()
        dobs[name] = d
        summ.append(dict(table=name, c_res=O.loc[(name, pair[0]), "coupling"], c_circ=O.loc[(name, pair[1]), "coupling"],
                         diff=d, null_mean=nl.mean(), null_sd=nl.std(ddof=1), null_q95=np.quantile(nl, 0.95),
                         mde=np.quantile(nl, 0.95) + 0.84 * nl.std(ddof=1), p_one_sided=(1 + (nl >= d).sum()) / (1 + len(nl)),
                         min_p=1 / (1 + len(nl)), n_perm=len(nl)))
    pd.DataFrame(summ).to_csv(L.OUT / "perm_summary.csv", index=False)
    # H152.1: mean over heart and breast
    nm = (N[N.table == "heart_state"].set_index("perm")["diff"] + N[N.table == "breast_state"].set_index("perm")["diff"]).to_numpy() / 2
    stat = (dobs["heart_state"] + dobs["breast_state"]) / 2
    pd.DataFrame([dict(statistic=stat, null_mean=nm.mean(), null_sd=nm.std(ddof=1), null_q95=np.quantile(nm, 0.95),
                       mde=np.quantile(nm, 0.95) + 0.84 * nm.std(ddof=1), p_one_sided=(1 + (nm >= stat).sum()) / (1 + len(nm)),
                       min_p=1 / (1 + len(nm)), n_perm=len(nm))]).to_csv(L.OUT / "h1521.csv", index=False)
    # 4. states summary: donors, cells, depth, resident share
    rows = []
    for name, (tgt, f, pair) in S.items():
        T, pb, pair = table(name)
        s = pb[pb.compartment.isin(pair)]
        w = s.pivot(index="donor", columns="compartment", values="n_cells")
        share = (w[pair[0]] / (w[pair[0]] + w[pair[1]])).dropna()
        for lab, c in zip(("resident_like", "circulating_like"), pair):
            x = s[s.compartment == c]
            rows.append(dict(table=name, state=lab, donors=x.donor.nunique(), donors_both=int(w.notna().all(axis=1).sum()),
                             median_cells=x.n_cells.median(), q25_cells=x.n_cells.quantile(.25), q75_cells=x.n_cells.quantile(.75),
                             median_depth=x.depth.median(), median_share_res=share.median(), q25_share=share.quantile(.25),
                             q75_share=share.quantile(.75)))
    pd.DataFrame(rows).to_csv(L.OUT / "states_summary.csv", index=False)
    print("DONE", f"{time.time() - t00:.0f}s", flush=True)


if __name__ == "__main__":
    main()
