"""A336 part 1 wrapper: A288 H288.1 (pooled T_MAC), H288.2 (pooled D_MAC = T_MAC - T_LY) and H288.3 (HLCA alveolar vs
monocyte-derived macrophages, A284 design) on the 14 genes, the core loop, the trio contribution, the gene-drop null and
the trio alone.

Original: A288_macrophage_pairwise/scripts/score_a288.py (imported, not edited; `table`, `stats`, `pooled`,
`placebo_chunk`, `read`) and A284_residency/scripts/score_a284.py (`delta`, `draws`; `unit` copied with extra waveforms
appended after the 14-gene block so that the thinning draws are identical). Companions drawn as originally (`pw.companions`,
seeds [288, atlas, family], 14-gene waveform) and restricted to the core positions. Stacked-state design of ct_lib.py.
The 5000-draw placebo null (seeds [288, d, atlas]) and the 2000-draw swap null (H288.3) are redrawn for the core
projector; the 14-gene p_cal of H288.1-2 come from the stored null_draws.npz after 4 recomputed draws are shown identical
(H288.3's swap null was not saved: both projectors redrawn). Core p_cal informative only if V283 passed on the core
projector (output/part1/a283.csv). Drops: 100 random + trio (clock only). Trio contribution: 200 pair-level donor
bootstraps within stratum, rng [336, 1] (H288.3: the same on the thinned replicates). Gate: trio gate (1000 permutations)
on every atlas entering the statistic.

    python p1_a288.py [--check] [--draws 5000] [--swap-draws 2000] [--jobs 4]
Output: ../output/part1/a288.csv
"""
import argparse
import importlib.util
import sys
import warnings

import ct_lib as L
from ct_lib import *  # noqa: F401,F403
import a336_common as C

sys.path.insert(0, str(L.PAPER / "exploratory/A288_macrophage_pairwise/scripts"))
sys.path.insert(0, str(L.PAPER / "exploratory/A284_residency/scripts"))
import cloudpickle  # noqa: E402
import score_a288 as A  # noqa: E402
import score_a284 as S  # noqa: E402
from clocksync import pairwise as pw  # noqa: E402

SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a288.py"
OUT = L.PAPER / "exploratory/A288_macrophage_pairwise/output"
V283_ROW = L.EXP / "A336_two_readouts/output/part1/a283.csv"
META = A.META
S.SEED = A.SEED
ALL = {"H288.1 T_MAC": ("A288_H288.1", "pooled T_MAC (macrophage vs other non-immune classes, 10 atlases)"),
       "H288.2 D_MAC": ("A288_H288.2", "pooled D_MAC = T_MAC - T_LY (8 atlases)")}


def prepare(n, k, pb, family, with_extra=None):
    """A288.prepare with the waveform list kept: states of the 14-gene block (and of `with_extra` waveforms after it)."""
    X, LY = A.classes(pb, n, family)
    keep = A.MAC[n] + X + LY
    t = pb[pb.compartment.isin(keep)].reset_index(drop=True)
    genes = [c for c in t.columns if c not in META]
    wfs = [cs.gtex()]
    for j, kind in enumerate(A.FAMILIES):
        wfs += pw.companions(t[META + genes], A.N_SETS, kind, np.random.default_rng([A.SEED, k, j]))
    return t, wfs, X, LY


def boot_atlas(p, Z2, rng):
    """Pair-level donor bootstrap of T_MAC and T_LY (2 waveforms): ((2,), (2,))."""
    pairs = pw.pair_index(p["key"])
    a = np.array([q["a"] for q in pairs], dtype=object)
    b = np.array([q["b"] for q in pairs], dtype=object)
    inn = lambda arr, s: np.array([x in s for x in arr], bool)  # noqa: E731
    mac = (inn(a, p["MAC"]) & inn(b, p["X"])) | (inn(b, p["MAC"]) & inn(a, p["X"]))
    ly = (inn(a, p["LY"]) & inn(b, p["X"])) | (inn(b, p["LY"]) & inn(a, p["X"]))
    sel = mac | ly
    Aa = np.full((len(Z2), len(pairs)), np.nan)
    if sel.any():
        Aa[:, sel] = L.boot_abs(Z2, [q for q, s in zip(pairs, sel) if s], rng)
    med = lambda m: np.nanmedian(Aa[:, m], 1) if m.any() else np.full(len(Z2), np.nan)  # noqa: E731
    return med(mac), med(ly)


def null_task(ds):
    return A.placebo_chunk(G["ps"], ds, G["use_mac"], G["use_d"])[:2]


def unit_w(n, res, circ, pb, family, k, wfs_extra):
    """score_a284.unit with extra waveforms after the 14-gene block (same thinning draws: `use` from the 14-gene sets)."""
    from clocksync.config import limit_threads as lt
    lt(1)
    warnings.filterwarnings("ignore")
    fam = {c: family(c) for c in pb.compartment.unique()}
    X = [c for c in fam if fam[c] not in S.IMMUNE and c not in (res, circ)]
    both = sorted(set(pb.donor[pb.compartment == res]) & set(pb.donor[pb.compartment == circ]))
    t = pb[pb.donor.isin(both) & pb.compartment.isin([res, circ] + X)].reset_index(drop=True)
    meta = ["donor", "compartment", "stratum", "n_cells", "depth"]
    genes = [c for c in t.columns if c not in meta]
    clock = list(cs.gtex().genes)
    wfs = [cs.gtex()]
    for j, kind in enumerate(("matched", "responsive")):
        wfs += pw.companions(t[meta + genes], S.N_SETS, kind, np.random.default_rng([S.SEED, k, j]))
    use = sorted({g for w in wfs for g in w.genes})
    all_w = wfs + (extra_waveforms(wfs, wfs_extra) if wfs_extra is not None else [])
    t = t[meta + use]
    ir = t.index[t.compartment == res].to_numpy()
    ic = t.index[t.compartment == circ].to_numpy()
    r_d, c_d = t.donor[ir].to_numpy(), t.donor[ic].to_numpy()
    ic = ic[pd.Index(c_d).get_indexer(r_d)]
    reps = []
    for rep in range(S.N_REP):
        rng = np.random.default_rng([S.SEED, k, 100 + rep])
        tt = t.copy()
        cnt = np.array(tt[use].to_numpy(), dtype=np.int64)
        s = tt[clock].sum(1).to_numpy(float)
        for a_, b_ in zip(ir, ic):
            hi, lo = (a_, b_) if s[a_] > s[b_] else (b_, a_)
            p = s[lo] / max(s[hi], 1.0)
            cnt[hi] = rng.binomial(cnt[hi], p)
            tt.loc[hi, "depth"] = tt.depth[hi] * p
        tt[use] = cnt
        reps.append(np.array([pw.states(tt, w) for w in all_w]))
    return dict(atlas=n, res=res, circ=circ, X=X, key=t[["donor", "compartment", "stratum"]], ir=ir, ic=ic, Z=reps,
                donors=len(both), table=t)


def delta_boot(u, rng):
    """Pair-level donor bootstrap of Delta for the 14-gene clock and the core clock (rows 0 and IC.start of Z)."""
    pairs = pw.pair_index(u["key"])
    sel_r = [i for i, p in enumerate(pairs) if u["res"] in (p["a"], p["b"]) and u["circ"] not in (p["a"], p["b"])]
    sel_c = [i for i, p in enumerate(pairs) if u["circ"] in (p["a"], p["b"]) and u["res"] not in (p["a"], p["b"])]
    r = np.full((S.N_REP, 2, len(pairs)), np.nan)
    for i in sel_r + sel_c:
        idx = L.resample(pairs[i]["s"], rng)
        pr = dict(ia=pairs[i]["ia"][idx], ib=pairs[i]["ib"][idx], s=pairs[i]["s"][idx])
        for q, Z in enumerate(u["Z"]):
            r[q, :, i] = L.fast_abs(Z[[0, IC.start]], pr)
    return np.mean([np.nanmedian(x[:, sel_r], 1) - np.nanmedian(x[:, sel_c], 1) for x in r], 0)


def swap_task(ds):
    return S.draws([G["u"]], ds)


def holm(tab):
    p = np.array([max(t["matched_p_cal"], t["responsive_p_cal"]) for t in tab])
    order, out, run = np.argsort(p), np.empty(len(p)), 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - r) * p[i]))
        out[i] = run
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--draws", type=int, default=5000)
    ap.add_argument("--swap-draws", type=int, default=2000)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    A150 = S._a280()
    family = A150.FAM.family
    names = list(A.MAC)[:3] if chk else list(A.MAC)
    tabs = {}
    for n in names:
        tabs[n], sub = A.table(n, A150)
        if sub is not None:
            tabs[n + "_sub"] = sub
    ref = pd.read_csv(OUT / "tests.csv").set_index("test")
    ref3 = pd.read_csv(OUT / "h288_3.csv").set_index("atlas")

    # ---- 14 genes (stage 1; no core waveform exists yet) ---------------------------------------------------------
    ps, aux = [], []
    for k, n in enumerate(names):
        t, wfs, X, LY = prepare(n, k, tabs[n], family)
        Z = np.array([pw.states(t, w) for w in wfs])
        ps.append(dict(atlas=n, key=t[["donor", "compartment", "stratum"]], Z=Z, X=set(X), LY=set(LY), MAC=set(A.MAC[n]), donors=t.donor.nunique()))
        aux.append((t[META + sorted({g for w in wfs for g in w.genes} | set(cs.gtex().genes))].copy(), wfs))
    obs_per, use_mac, use_d = [], [], []
    for p in ps:
        tm, tl, nm, nl = A.stats(p, p["key"])
        obs_per.append((tm, tl))
        use_mac.append(nm > 0)
        use_d.append(nm > 0 and nl > 0)
    per = [(np.where(use_mac[i], tm, np.nan), np.where(use_d[i], tl, np.nan)) for i, (tm, tl) in enumerate(obs_per)]
    t_obs, d_obs, tm_a, d_a = A.pooled(per)
    G.update(ps=ps, use_mac=use_mac, use_d=use_d)
    ndraw = 8 if chk else a.draws
    z = np.load(OUT / "null_draws.npz")
    n14d = 8 if chk else 4  # draws recomputed for the 14-gene projector (the rest is read from the stored null)
    t0 = time.time()
    res = L.run_pool(null_task, [np.arange(n14d)], 1)
    n14 = (res[0][0], res[0][1])
    log(f"14-gene null check: {n14d} draws recomputed in {(time.time() - t0) / 60:.1f} min")
    if not chk:
        assert np.allclose(n14[0], z["T_MAC"][:n14d], equal_nan=True, atol=1e-9) and np.allclose(n14[1], z["D_MAC"][:n14d], equal_nan=True, atol=1e-9), "stored null differs"
    chunks = np.array_split(np.arange(ndraw), max(ndraw // 25, 1))
    nulls14 = {"H288.1 T_MAC": z["T_MAC"] if not chk else n14[0], "H288.2 D_MAC": z["D_MAC"] if not chk else n14[1]}
    obs14 = {"H288.1 T_MAC": t_obs, "H288.2 D_MAC": d_obs}
    bad, tab14 = [], {}
    for t, o in obs14.items():
        tab14[t] = A.read(o, nulls14[t])
    if not chk:
        for t, o in obs14.items():
            row = ref.loc[t]
            chk_ = [("clock", o[0], row["clock"])] + [(k, tab14[t][k], row[k]) for k in tab14[t]]
            bad += [(t, nm, g, w) for nm, g, w in chk_ if not np.isclose(g, w, atol=1e-6, rtol=0)]
        # H288.3 14-gene
        bp, _ = A.B.hlca_split()
        st = tabs["lung_hlca"].drop_duplicates("donor").set_index("donor").stratum
        add = bp[bp.compartment.isin(["lung:macrophage_alveolar", "lung:macrophage_monocyte_derived"]) & bp.donor.isin(st.index)].copy()
        add["stratum"] = add.donor.map(st)
        pb3 = pd.concat([tabs["lung_hlca"][~tabs["lung_hlca"].compartment.isin(A.MAC["lung_hlca"])],
                         add.reindex(columns=tabs["lung_hlca"].columns).fillna(0.0)], ignore_index=True)
        res_, circ_, _ = A.H3["lung_hlca"]
        u14 = unit_w("lung_hlca", res_, circ_, pb3, family, 100, None)
        d14, nr, nc = S.delta(u14)
        G["u"] = u14
        sw14 = np.concatenate(L.run_pool(swap_task, np.array_split(np.arange(a.swap_draws), 40), a.jobs))
        r14 = A.read(d14, sw14)
        p_swap14 = float((1 + np.sum(sw14[:, 0] >= d14[0])) / (1 + len(sw14)))
        for nm, g, w in (("delta", d14[0], ref3.loc["lung_hlca", "delta"]), ("p_swap", p_swap14, ref3.loc["lung_hlca", "p_swap"])) + tuple(
                (k, r14[k], ref3.loc["lung_hlca", k]) for k in r14):
            if not np.isclose(g, w, atol=1e-6, rtol=0):
                bad.append(("H288.3", nm, g, w))
        log("14-gene reproduction:", "OK" if not bad else bad)
    if bad:
        C.write_rows("a288", [dict(id=i, figure_panel="Fig 2", statistic=s, reproduced=False,
                                   notes=f"14-gene values not reproduced (tolerance 1e-6): {bad}", script=SCRIPT)
                              for i, s in [(v[0], v[1]) for v in ALL.values()] + [("A288_H288.3", "HLCA alveolar vs monocyte-derived macrophages")]])
        return

    # ---- core, drops, trio states ------------------------------------------------------------------------------------
    dl, tl = C.drops(n=100)
    if chk:
        dl = dl[:3] + [x for x in dl if x[0] == tl]
    lay = L.layout(len(dl))
    for p, (tt, wfs) in zip(ps, aux):
        ex = np.array([pw.states(tt, w) for w in extra_waveforms(wfs, dl)])
        p["Z14"], p["Z"] = p["Z"], np.vstack([p["Z"], ex])
    obs_per = [A.stats(p, p["key"])[:2] for p in ps]
    per = [(np.where(use_mac[i], tm, np.nan), np.where(use_d[i], tl_, np.nan)) for i, (tm, tl_) in enumerate(obs_per)]
    t_all, d_all, _, _ = A.pooled(per)
    G["ps"] = [dict(p, Z=p["Z"][IC]) for p in ps]
    t0 = time.time()
    res = L.run_pool(null_task, chunks, a.jobs)
    nc_ = (np.concatenate([r[0] for r in res]), np.concatenate([r[1] for r in res]))
    log(f"core null: {ndraw} draws in {(time.time() - t0) / 60:.1f} min")
    nullc = {"H288.1 T_MAC": nc_[0], "H288.2 D_MAC": nc_[1]}
    obs_all = {"H288.1 T_MAC": t_all, "H288.2 D_MAC": d_all}
    v283 = pd.read_csv(V283_ROW).set_index("id").loc["A283_V283"] if V283_ROW.exists() else None
    v283_ok = (v283 is not None and v283.core_status == "holds") or chk
    unin = "uninformative (V283 fails on the core projector or is not available)"

    # ---- trio contribution and gates ---------------------------------------------------------------------------------
    brng = np.random.default_rng([336, 1])
    diffs = {t: [] for t in ALL}
    for b in range(3 if chk else 200):
        per = []
        for i, p in enumerate(ps):
            tm, tl_ = boot_atlas(p, p["Z"][[0, IC.start]], brng)
            per.append((np.where(use_mac[i], tm, np.nan), np.where(use_d[i], tl_, np.nan)))
        tb_, db_, _, _ = A.pooled(per)
        diffs["H288.1 T_MAC"].append(tb_[0] - tb_[1])
        diffs["H288.2 D_MAC"].append(db_[0] - db_[1])
    gates = {}
    for i, p in enumerate(ps):
        if use_mac[i]:
            tt = aux[i][0]
            gates[p["atlas"]] = L.gate_atlas(tt, C.trio(), 50 if chk else 1000, p["atlas"])
    not chk and log("gates", {n: (g[0], round(g[1], 3), round(g[2], 3)) for n, g in gates.items()})
    dwf = dict(dl)

    tabc = {t: A.read(obs_all[t][IC], nullc[t]) for t in ALL}
    t14 = {t: A.read(obs_all[t][I14], nulls14[t]) for t in ALL}
    hc = holm([tabc[t] for t in ALL])
    h14 = holm([t14[t] for t in ALL])
    rows = []
    for j, (t, (rid, nm)) in enumerate(ALL.items()):
        v = obs_all[t]
        dv = {lab: float(v[lay["drops"]][q]) for q, (lab, _) in enumerate(dl)}
        assert np.isclose(dv[tl], v[IC][0], atol=1e-9)
        pr = C.drop_percentiles(dv, tl, dwf, +1)
        used = [p["atlas"] for i, p in enumerate(ps) if (use_mac[i] if t.startswith("H288.1") else use_d[i])]
        gp = [gates[n] for n in used]
        gate_all = bool(all(g[0] for g in gp))
        lo, hi = L.ci90(diffs[t])
        c, r14 = tabc[t], t14[t]
        if t.startswith("H288.1"):
            met = [c["matched_p_cal"] < 0.05, c["responsive_p_cal"] < 0.05]
            crit = "p_cal < 0.05 against both families (core, 5000 redrawn placebo draws; Holm p reported)"
        else:
            met = [v[IC][0] > 0, hc[j] < 0.05]
            crit = "D_MAC > 0 with p_cal < 0.05 against both families, Holm over H288.1-2 (core, 5000 redrawn placebo draws)"
        st = L.status(met) if v283_ok else unin
        rows.append(dict(
            id=rid, figure_panel="Fig 2 (macrophages)", statistic=nm + ", calibrated p", direction=1,
            value_14_stored=float(ref.loc[t, "clock"]), value_14=float(v[I14][0]), reproduced=True, value_core=float(v[IC][0]),
            core_criterion=crit + ", given V283 passes on the core projector", core_status=st,
            diff_14_minus_core=float(v[I14][0] - v[IC][0]), diff_lo90=lo, diff_hi90=hi,
            diff_draws="200 pair-level donor bootstraps within stratum, rng [336,1], shared (clock statistic)",
            **{k: pr[k] for k in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")},
            trio_gate_passes=gate_all, trio_alone_value=float(v[lay["trio"]]) if gate_all else np.nan,
            reading=C.reading(st, pr["trio_drop_below_q10"]) if v283_ok else st, blind="core blind",
            notes=(f"value = clock statistic (the p_cal are here): matched p_cal 14 {r14['matched_p_cal']:.4f} (stored {ref.loc[t, 'matched_p_cal']:.4f}) / core "
                   f"{c['matched_p_cal']:.4f} (pi {c['matched_pi']:.3f}); responsive p_cal 14 {r14['responsive_p_cal']:.4f} (stored {ref.loc[t, 'responsive_p_cal']:.4f}) "
                   f"/ core {c['responsive_p_cal']:.4f} (pi {c['responsive_pi']:.3f}); Holm p 14 {h14[j]:.4f} / core {hc[j]:.4f}; V283 core "
                   f"{'passes' if v283_ok else 'fails or missing'}; trio gate passes in {sum(g[0] for g in gp)} of {len(gp)} atlases"), script=SCRIPT))

    # ---- H288.3 (HLCA), core, drops, trio ---------------------------------------------------------------------------
    bp, _ = A.B.hlca_split()
    base = tabs["lung_hlca"]
    st_ = base.drop_duplicates("donor").set_index("donor").stratum
    res_, circ_, _ = A.H3["lung_hlca"]
    add = bp[bp.compartment.isin([res_, circ_]) & bp.donor.isin(st_.index)].copy()
    add["stratum"] = add.donor.map(st_)
    pb3 = pd.concat([base[~base.compartment.isin(A.MAC["lung_hlca"])], add.reindex(columns=base.columns).fillna(0.0)], ignore_index=True)
    u = unit_w("lung_hlca", res_, circ_, pb3, family, 100, dl)
    dall, nr, nc = S.delta(u)
    G["u"] = dict(u, Z=[Z[:IC.stop] for Z in u["Z"]])
    sw = np.concatenate(L.run_pool(swap_task, np.array_split(np.arange(4 if chk else a.swap_draws), 2 if chk else 40), a.jobs))
    sw14c, swc = sw[:, I14], sw[:, IC]
    rc, r14 = A.read(dall[IC], swc), A.read(dall[I14], sw14c)
    p_sw = lambda d_, nl: float((1 + np.sum(nl[:, 0] >= d_[0])) / (1 + len(nl)))  # noqa: E731
    dd = []
    for b in range(3 if chk else 200):
        dd.append(delta_boot(u, brng))
    dd = np.array(dd)
    dv = {lab: float(dall[lay["drops"]][q]) for q, (lab, _) in enumerate(dl)}
    assert np.isclose(dv[tl], dall[IC][0], atol=1e-9)
    pr = C.drop_percentiles(dv, tl, dwf, +1)
    g3 = L.gate_atlas(u["table"], C.trio(), 50 if chk else 1000, "hlca_h288_3")
    lo, hi = L.ci90(dd[:, 0] - dd[:, 1])
    met = [rc["matched_p_cal"] < 0.05, rc["responsive_p_cal"] < 0.05]
    st = L.status(met) if v283_ok else unin
    rows.append(dict(
        id="A288_H288.3", figure_panel="Fig 2 (macrophages)", statistic="HLCA alveolar vs monocyte-derived macrophages, Delta of clock |r| (A284 design), calibrated p",
        direction=1, value_14_stored=float(ref3.loc["lung_hlca", "delta"]), value_14=float(dall[I14][0]), reproduced=True,
        value_core=float(dall[IC][0]),
        core_criterion="resident above monocyte-derived: p_cal < 0.05 against both families (core, 2000 redrawn swap draws), given V283 passes on the core projector",
        core_status=st, diff_14_minus_core=float(dall[I14][0] - dall[IC][0]), diff_lo90=lo, diff_hi90=hi,
        diff_draws="200 pair-level donor bootstraps within stratum on the 10 thinned replicates, rng [336,1], shared",
        **{k: pr[k] for k in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")},
        trio_gate_passes=g3[0], trio_alone_value=float(dall[lay["trio"]]) if g3[0] else np.nan,
        reading=C.reading(st, pr["trio_drop_below_q10"]) if v283_ok else st, blind="core blind",
        notes=(f"value = clock Delta; p_swap 14 {p_sw(dall[I14], sw14c):.4f} / core {p_sw(dall[IC], swc):.4f}; matched p_cal 14 {r14['matched_p_cal']:.4f} "
               f"(stored {ref3.loc['lung_hlca', 'matched_p_cal']:.4f}) / core {rc['matched_p_cal']:.4f}; responsive p_cal 14 {r14['responsive_p_cal']:.4f} (stored "
               f"{ref3.loc['lung_hlca', 'responsive_p_cal']:.4f}) / core {rc['responsive_p_cal']:.4f}; the 14-gene criterion itself fails (stored passes=False); "
               f"V283 core {'passes' if v283_ok else 'fails or missing'}; trio gate p_pairing {g3[1]:.4f}, p_acrophase {g3[2]:.4f}; swap null redrawn for both projectors"),
        script=SCRIPT))
    if chk:
        assert len(rows) == 3
        print("check ok")
        return
    df = C.write_rows("a288", rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes"]).round(4).T.to_string(), flush=True)
    log("done")


if __name__ == "__main__":
    main()
