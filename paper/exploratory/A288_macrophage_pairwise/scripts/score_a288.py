"""A288: do macrophages keep the tissue's shared clock state as stromal-vascular cells do?

Registered in paper/docs/registered_experiments.md (A288). Anchor-free pairwise readout (`clocksync.pairwise`), A283
form. Per atlas: T_MAC = median clock |r| over (macrophage, X) pairs, T_LY over (lymphocyte, X) pairs, X = non-immune
classes; D_MAC = T_MAC - T_LY; pooled = mean over the atlases with >= 1 pair. 200 matched + 200 responsive companion
sets per atlas (seeds [288, atlas, family]); null = 5000 placebo draws (seeds [288, d, atlas]); p_cal as A283.

  H288.1  pooled T_MAC: p_cal < 0.05 against both families
  H288.2  pooled D_MAC > 0 with p_cal < 0.05 against both families          (Holm over H288.1-2)
  H288.3  HLCA alveolar against monocyte-derived macrophages, A284 design (thinning, swap null, 2000 draws);
          IBDverse CD163hi against CD163lo descriptive

    python score_a288.py [--design] [--draws 5000] [--swap-draws 2000] [--jobs 6]
Outputs (../output/): design.csv, design_h3.csv, tests.csv, per_atlas.csv, h288_3.csv, null_draws.npz
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(PAPER / "exploratory/A284_residency/scripts"))
import cloudpickle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import build_a288 as B  # noqa: E402
import score_a284 as S  # noqa: E402
from clocksync import pairwise as pw  # noqa: E402
from common import cli  # noqa: E402

cloudpickle.register_pickle_by_value(S)
SEED, N_SETS, CHUNK = 288, 200, 500
FAMILIES = {"matched": slice(1, 1 + N_SETS), "responsive": slice(1 + N_SETS, 1 + 2 * N_SETS)}
IMMUNE = ("lymphocyte", "myeloid")
META = ["donor", "compartment", "stratum", "n_cells", "depth"]
MAC = {"lung_hlca": ["lung:macrophage"], "organ_breast": ["breast:macrophage_DC"],
       "organ_colon": ["colon:macrophage_DC"], "organ_muscle": ["muscle:macrophage_DC"],
       "organ_optic_nerve": ["optic_nerve:head_macrophage_DC", "optic_nerve:nerve_macrophage_DC"],
       "organ_eye_anterior": ["eye_anterior:ciliary_macrophage_DC", "eye_anterior:trabecular_macrophage_DC"],
       "kong_ileum": ["ileum:macrophage_DC"], "kong_colon": ["colon:macrophage_DC"],
       "ibd_ileum": ["ileum:macrophage"], "ibd_rectum": ["rectum:macrophage"]}
GENOME = {**S.ATLAS_GENOME, "organ_muscle": "organ_muscle_genome", "organ_optic_nerve": "organ_optic_nerve_genome",
          "organ_eye_anterior": "organ_eye_anterior_genome"}
H3 = {"lung_hlca": ("lung:macrophage_alveolar", "lung:macrophage_monocyte_derived", True),
      "ibd_ileum": ("ileum:mac_CD163hi", "ileum:mac_CD163lo", False),
      "ibd_rectum": ("rectum:mac_CD163hi", "rectum:mac_CD163lo", False)}


def table(n, A150):
    """Healthy rows of the atlas (A17 targets, genome-wide); Kong ileum epithelial without IELs (A284); IBDverse with
    the macrophage rows (`build_a288.ibd_macrophage`) on the same donors and strata."""
    A150.ATLASES.update(GENOME)
    if n == "kong_ileum":
        pb, _ = S.atlas_table(n, A150)
        pb = pb[~pb.compartment.str.endswith(("T_resident", "T_circulating"))]
    else:
        pb = A150.load_atlas(n, None, 0)["pb"]
    pb = pb.reset_index(drop=True)
    sub = None
    if n.startswith("ibd_"):
        st = pb.drop_duplicates("donor").set_index("donor").stratum
        genes = [c for c in pb.columns if c not in META]
        rows = B.ibd_macrophage(n.removeprefix("ibd_"), set(st.index), genes)
        rows["stratum"] = rows.donor.map(st)
        rows = rows.reindex(columns=pb.columns).fillna(0.0)
        sub = rows[rows.compartment.str.contains(":mac_")]
        pb = pd.concat([pb, rows[~rows.compartment.str.contains(":mac_")]], ignore_index=True)
    return pb, sub


def classes(pb, n, family):
    fam = {c: family(c) for c in pb.compartment.unique()}
    X = sorted(c for c, f in fam.items() if f not in IMMUNE and c not in MAC[n])
    LY = sorted(c for c, f in fam.items() if f == "lymphocyte")
    return X, LY


def prepare(n, k, pb, family):
    """States of the clock and the 400 companion sets; pair masks."""
    from clocksync.config import limit_threads
    limit_threads(1)
    warnings.filterwarnings("ignore")
    import clocksync as cs
    X, LY = classes(pb, n, family)
    keep = MAC[n] + X + LY
    t = pb[pb.compartment.isin(keep)].reset_index(drop=True)
    genes = [c for c in t.columns if c not in META]
    wfs = [cs.gtex()]
    for j, kind in enumerate(FAMILIES):
        wfs += pw.companions(t[META + genes], N_SETS, kind, np.random.default_rng([SEED, k, j]))
    Z = np.array([pw.states(t, w) for w in wfs])
    return dict(atlas=n, key=t[["donor", "compartment", "stratum"]], Z=Z, X=set(X), LY=set(LY), MAC=set(MAC[n]),
                donors=t.donor.nunique())


def stats(p, key):
    """T_MAC and T_LY (sets) of one atlas on the pairing of `key`, and their pair counts."""
    pairs = pw.pair_index(key)
    a = np.array([q["a"] for q in pairs], dtype=object)
    b = np.array([q["b"] for q in pairs], dtype=object)
    inn = lambda arr, s: np.array([x in s for x in arr], bool)  # noqa: E731
    mac = (inn(a, p["MAC"]) & inn(b, p["X"])) | (inn(b, p["MAC"]) & inn(a, p["X"]))
    ly = (inn(a, p["LY"]) & inn(b, p["X"])) | (inn(b, p["LY"]) & inn(a, p["X"]))
    sel = mac | ly
    A = np.full((len(p["Z"]), len(pairs)), np.nan)
    if sel.any():
        A[:, sel] = pw.abs_r(p["Z"], [q for q, s in zip(pairs, sel) if s])
    med = lambda m: np.nanmedian(A[:, m], 1) if m.any() else np.full(len(p["Z"]), np.nan)  # noqa: E731
    return med(mac), med(ly), int(mac.sum()), int(ly.sum())


def pooled(per):
    """Pooled T_MAC and D_MAC (sets) from per-atlas (T_MAC, T_LY) pairs: mean over atlases with a value."""
    tm = np.array([x[0] for x in per])
    d = np.array([x[0] - x[1] for x in per])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmean(tm, 0), np.nanmean(d, 0), tm, d


def placebo_chunk(ps, ds, use_mac, use_d):
    from clocksync.config import limit_threads
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out_t, out_d, at, ad_ = [], [], [], []
    for d in ds:
        per = []
        for i, p in enumerate(ps):
            tm, tl, _, _ = stats(p, pw.placebo(p["key"], np.random.default_rng([SEED, d, i])))
            per.append((np.where(use_mac[i], tm, np.nan), np.where(use_d[i], tl, np.nan)))
        t, dd, tm_a, d_a = pooled(per)
        out_t.append(t), out_d.append(dd), at.append(tm_a), ad_.append(d_a)
    return np.array(out_t), np.array(out_d), np.array(at), np.array(ad_)


def read(obs, null):
    row = {}
    for f, sl in FAMILIES.items():
        o = pw.percentile(obs, sl)[0]
        row.update({f"{f}_median": float(np.nanmedian(obs[sl])), f"{f}_pi": float(o),
                    f"{f}_p_cal": float(pw.calibrated_p(o, pw.percentile(null, sl))[0])})
    return row


def h3(tabs, A150, family, swap_draws, jobs):
    """H288.3: A284 units (thinning, swap null) for the macrophage splits."""
    S.SEED = SEED
    rows, design = [], []
    for k, (n, (res, circ, decides)) in enumerate(H3.items()):
        if n == "lung_hlca":
            sp_, _ = B.hlca_split()
            base = tabs[n]
            st = base.drop_duplicates("donor").set_index("donor").stratum
            add = sp_[sp_.compartment.isin([res, circ]) & sp_.donor.isin(st.index)].copy()
            add["stratum"] = add.donor.map(st)
            pb = pd.concat([base[~base.compartment.isin(MAC[n])], add.reindex(columns=base.columns).fillna(0.0)],
                           ignore_index=True)
        else:
            base, sub = tabs[n], tabs[n + "_sub"]
            pb = pd.concat([base[~base.compartment.isin(MAC[n])], sub], ignore_index=True)
        both = set(pb.donor[pb.compartment == res]) & set(pb.donor[pb.compartment == circ])
        design.append(dict(atlas=n, resident=res, monocyte_derived=circ, donors_both=len(both), decides=decides,
                           X=";".join(sorted(c for c in pb.compartment.unique()
                                             if family(c) not in IMMUNE and c not in (res, circ)))))
        if len(both) < 15:
            continue
        u = S.unit(n, res, circ, pb, family, 100 + k)
        d, nr, nc = S.delta(u)
        if nr == 0 or nc == 0:
            design[-1]["eligible"] = False
            continue
        design[-1]["eligible"] = True
        chunks = np.array_split(np.arange(swap_draws), jobs)
        null = np.concatenate(Parallel(n_jobs=jobs)(delayed(S.draws)([u], c) for c in chunks))
        row = dict(atlas=n, decides=decides, donors=len(both), pairs_res=nr, pairs_md=nc, delta=float(d[0]),
                   p_swap=float((1 + np.sum(null[:, 0] >= d[0])) / (1 + len(null))))
        row.update(read(d, null))
        row["passes"] = bool(decides and row["delta"] > 0 and row["p_swap"] < 0.05) if decides else None
        rows.append(row)
    return pd.DataFrame(design), pd.DataFrame(rows)


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--design", action="store_true")
    ap.add_argument("--draws", type=int, default=5000)
    ap.add_argument("--swap-draws", type=int, default=2000)
    args = ap.parse_args()
    cli.setup(args)
    out = HERE.parents[1] / "output"
    out.mkdir(parents=True, exist_ok=True)
    A150 = S._a280()
    family = A150.FAM.family
    tabs = {}
    for n in MAC:
        tabs[n], sub = table(n, A150)
        if sub is not None:
            tabs[n + "_sub"] = sub
    ps = Parallel(n_jobs=args.jobs)(delayed(prepare)(n, k, tabs[n], family) for k, n in enumerate(MAC))
    design, obs_per = [], []
    for p in ps:
        tm, tl, nm, nl = stats(p, p["key"])
        obs_per.append((tm, tl))
        design.append(dict(atlas=p["atlas"], donors=p["donors"], MAC=";".join(sorted(p["MAC"])), n_X=len(p["X"]),
                           n_LY=len(p["LY"]), pairs_MAC_X=nm, pairs_LY_X=nl,
                           X=";".join(sorted(p["X"])), LY=";".join(sorted(p["LY"]))))
    D = pd.DataFrame(design)
    D["in_T_MAC"] = D.pairs_MAC_X > 0
    D["in_D_MAC"] = (D.pairs_MAC_X > 0) & (D.pairs_LY_X > 0)
    D.to_csv(out / "design.csv", index=False)
    pd.set_option("display.width", 250)
    print(D.drop(columns=["X", "LY"]).to_string(), flush=True)
    if args.design:
        d3 = []
        for n, (res, circ, decides) in H3.items():
            if n == "lung_hlca":
                sp_, _ = B.hlca_split()
                pb = sp_[sp_.donor.isin(set(tabs[n].donor))]
            else:
                pb = tabs[n + "_sub"]
            both = set(pb.donor[pb.compartment == res]) & set(pb.donor[pb.compartment == circ])
            d3.append(dict(atlas=n, resident=res, monocyte_derived=circ, donors_both=len(both), decides=decides,
                           res_donors=pb.donor[pb.compartment == res].nunique(),
                           md_donors=pb.donor[pb.compartment == circ].nunique()))
        d3 = pd.DataFrame(d3)
        d3.to_csv(out / "design_h3.csv", index=False)
        print(d3.to_string())
        return

    use_mac = [bool(x) for x in D.in_T_MAC]
    use_d = [bool(x) for x in D.in_D_MAC]
    per = [(np.where(use_mac[i], tm, np.nan), np.where(use_d[i], tl, np.nan)) for i, (tm, tl) in enumerate(obs_per)]
    t_obs, d_obs, tm_a, d_a = pooled(per)
    nt, nd, nta, nda = [], [], [], []
    for a in range(0, args.draws, CHUNK):
        ds = np.arange(a, min(a + CHUNK, args.draws))
        res = Parallel(n_jobs=args.jobs)(delayed(placebo_chunk)(ps, c, use_mac, use_d)
                                         for c in np.array_split(ds, args.jobs))
        for r in res:
            nt.append(r[0]), nd.append(r[1]), nta.append(r[2]), nda.append(r[3])
        print("draws", ds[-1] + 1, "of", args.draws, flush=True)
    nt, nd, nta, nda = (np.concatenate(x) for x in (nt, nd, nta, nda))
    np.savez_compressed(out / "null_draws.npz", T_MAC=nt, D_MAC=nd)
    tests = []
    for h, obs, null in (("H288.1 T_MAC", t_obs, nt), ("H288.2 D_MAC", d_obs, nd)):
        row = dict(test=h, atlases=int(np.isfinite([x[0] for x in (tm_a if h.endswith("T_MAC") else d_a)]).sum()),
                   clock=float(obs[0]), **read(obs, null))
        row["p_max"] = max(row["matched_p_cal"], row["responsive_p_cal"])
        tests.append(row)
    T = pd.DataFrame(tests)
    order = np.argsort(T.p_max.to_numpy())
    holm, run = np.empty(len(T)), 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(T) - r) * T.p_max[i]))
        holm[i] = run
    T["p_holm"] = holm
    T["passes"] = (T.clock > 0) & (T.p_holm < 0.05)
    T.to_csv(out / "tests.csv", index=False)
    pa = []
    for i, p in enumerate(ps):
        for stat, obs, null in (("T_MAC", tm_a[i], nta[:, i]), ("D_MAC", d_a[i], nda[:, i])):
            if np.isfinite(obs[0]):
                pa.append(dict(atlas=p["atlas"], statistic=stat, clock=float(obs[0]), **read(obs, null)))
    pd.DataFrame(pa).to_csv(out / "per_atlas.csv", index=False)
    print(T.round(4).to_string(), pd.DataFrame(pa).round(3).to_string(), sep="\n", flush=True)
    d3, r3 = h3(tabs, A150, family, args.swap_draws, args.jobs)
    d3.to_csv(out / "design_h3.csv", index=False)
    r3.to_csv(out / "h288_3.csv", index=False)
    print(d3.drop(columns="X").to_string(), r3.round(4).to_string(), sep="\n")


if __name__ == "__main__":
    main()
