"""A122: controls for the corrected A85 immune cue-following result (locked-set control, pseudo-cue sets).

Registered in paper/docs/registered_experiments.md (A122). Reuses the A85 functions (cancer/compute_cue_following.py).
  H122.1: (immune - structural, cue index) - (immune - structural, locked index); atlas bootstrap (2000), one-sided.
  H122.2: 200 expression-matched random gene sets (each cue gene replaced by a random non-clock/cue/locked gene among
          the 25 nearest in mean log CPM in that atlas, carrying the cue gene's amplitude and acrophase); the observed
          immune - structural cue difference must exceed their 95th percentile.
Weights: inverse-variance with unit-bootstrap SEs (500 draws) per set; pseudo sets use the cue-set weights.
    python a122_cue_controls.py [--jobs 4] [--quick]
Outputs (../output): reproduction.csv, h122.csv, per_atlas.csv, depth.csv, composition.csv (+ classes_*.csv, pseudo.csv)
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER.parent))
sys.path.insert(0, str(PAPER))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA, group_out  # noqa: E402

spec = importlib.util.spec_from_file_location("a85", PAPER / "results/figure_5/scripts/cancer/compute_cue_following.py")
a85 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a85)

OUT = HERE.parent / "output"
SEED, N_BOOT, N_PSEUDO, N_ATLAS_BOOT, POOL, NEAR = 122, 500, 200, 2000, 3000, 25
W = a85.W
IMMUNE, STRUCTURAL = a85.IMMUNE, a85.STRUCTURAL


def wave(ref, name):
    x = ref[ref.set == name]
    return x.index.tolist(), x.amp.to_numpy(float), x.acro_h.to_numpy(float) * W


def mkwf(genes, amp, psi):
    return cs.Waveform(tuple(genes), np.asarray(amp), np.asarray(psi), None, None, "cue")


def group_of(k):
    return "immune" if k in IMMUNE else "structural" if k in STRUCTURAL else "other"


def work(args):
    name, draws85, draws122, quick = args
    ref = pd.read_csv(group_out("gtex_cues") / "a81b_cue_reference.csv", index_col=0)
    core = cs.gtex()
    core_g = list(core.genes)
    sets = {s: wave(ref, s) for s in ("cue", "locked")}
    sets = {s: ([canonical(g) or g for g in v[0]], v[1], v[2]) for s, v in sets.items()}
    allg = set(core_g) | set(sets["cue"][0]) | set(sets["locked"][0])
    cu = data._builder("cue_units")
    a = ad.read_h5ad(DATA / "atlas_wg" / name / "pseudobulk.h5ad", backed="r")
    names = a.var.feature_name.astype(str) if "feature_name" in a.var else a.var_names.astype(str)
    rng = np.random.default_rng(SEED + sum(map(ord, name)))
    cand = [g for g in pd.unique(names) if g not in allg and not g.startswith(("ENSG", "MT-"))]
    cand = list(rng.choice(cand, min(POOL, len(cand)), replace=False))
    t = cu.table(name, list(allg) + cand)
    if quick:
        t = t[t.unit.isin(np.random.default_rng(0).permutation(t.unit.unique())[:150])]
    cand = [g for g in cand if g in t.columns]
    lcpm = lambda g: float(np.mean(np.log2(1 + 1e6 * t[g].to_numpy(float) / t.depth.to_numpy(float))))
    expr = pd.Series({g: lcpm(g) for g in cand})
    # core state once
    t0 = t.assign(c=a85.state(t, core))
    t0["c"] = t0.groupby("cls", group_keys=False).apply(
        lambda x: pd.Series(a85.centre(x.c.to_numpy(), x.stratum), index=x.index))

    def design(genes, amp, psi):
        pres = [i for i, g in enumerate(genes) if g in t.columns and g not in core_g]
        gp = [genes[i] for i in pres]
        wf = mkwf(gp, amp[pres], psi[pres])
        out = {}
        for k in t0.cls.unique():
            own = t0[t0.cls == k].set_index("unit")
            rest = t0[(t0.cls != k) & t0.unit.isin(own.index)]
            cnt = rest.groupby("unit").size()
            units = cnt[cnt >= 2].index
            if len(units) < a85.MIN_UNITS:
                continue
            H = rest.groupby("unit").c.mean().loc[units]
            summed = rest.groupby("unit")[gp + ["depth"]].sum().loc[units]
            Q = a85.centre(a85.state(summed, wf), own.loc[units, "stratum"])
            out[k] = (own.loc[units, "c"].to_numpy(), a85.centre(H.to_numpy(), own.loc[units, "stratum"]), Q)
        return out, wf

    res = dict(name=name, rows=[], pseudo=[], comp=[])
    n_boot = 50 if quick else N_BOOT
    for s, draws in (("cue", draws85), ("locked", draws122)):
        des, wf = design(*sets[s])
        for k, (c, H, Q) in des.items():
            n = len(c)
            bs = np.array([a85.fit(c[i], H[i], Q[i]) for i in
                           (np.random.default_rng(d).integers(0, n, n) for d in draws[:n_boot])])
            res["rows"].append(dict(set=s, atlas=name, cls=k, units=n, index=a85.fit(c, H, Q), se=float(bs.std()),
                                    corr_host=a85.corr(c, H), corr_cue=a85.corr(c, Q), corr_host_cue=a85.corr(H, Q),
                                    group=group_of(k), boot=bs))
        # own-class state variance (raw state centred within stratum, not scaled)
        for k in des:
            x = t0[t0.cls == k]
            z = pd.Series(a85.state(x, wf)).reset_index(drop=True)
            z = z - z.groupby(x.stratum.to_numpy()).transform("mean")
            res["comp"].append(dict(set=s, atlas=name, cls=k, own_state_var=float(np.mean(np.abs(z.to_numpy()) ** 2)),
                                    log_depth=float(np.median(np.log10(x.depth))), n_cells=float(x.n_cells.median())))
    # pseudo sets
    cg, camp, cpsi = sets["cue"]
    pres = [i for i, g in enumerate(cg) if g in t.columns and g not in core_g]
    ce = {cg[i]: lcpm(cg[i]) for i in pres}
    prng = np.random.default_rng(SEED * 1000 + sum(map(ord, name)))
    for j in range(N_PSEUDO if not quick else 10):
        used, rep = set(), []
        for i in pres:
            order = (expr - ce[cg[i]]).abs().sort_values().index
            pick = [g for g in order if g not in used][:NEAR]
            g = pick[prng.integers(len(pick))]
            used.add(g)
            rep.append(g)
        des, _ = design(rep, camp[pres], cpsi[pres])
        for k, (c, H, Q) in des.items():
            res["pseudo"].append((j, k, a85.fit(c, H, Q)))
    return res


def pool(d, idx, w):
    out = {}
    for g in ("immune", "structural"):
        m = (d.group == g).to_numpy()
        out[g] = np.sum(w[m] * idx[m]) / np.sum(w[m])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    from clocksync.config import limit_threads
    limit_threads(1)
    cu = data._builder("cue_units")
    atl = ["lung_hlca", "heart_global"] if args.quick else cu.ATLASES
    nb = 50 if args.quick else N_BOOT
    r85, r122 = np.random.default_rng(85), np.random.default_rng(SEED)
    jobs = [(n, [r85.integers(0, 10 ** 9) for _ in range(nb)], [r122.integers(0, 10 ** 9) for _ in range(nb)], args.quick)
            for n in atl]
    with Pool(args.jobs) as p:
        results = p.map(work, jobs, chunksize=1)
    out = OUT / "quick" if args.quick else OUT
    out.mkdir(exist_ok=True, parents=True)
    rows = pd.DataFrame([r for x in results for r in x["rows"]])
    boots = {s: np.array([b for b in rows[rows.set == s].boot]) for s in ("cue", "locked")}
    rows.drop(columns="boot").to_csv(out / "classes_cue_locked.csv", index=False)
    D = {s: rows[rows.set == s].reset_index(drop=True) for s in ("cue", "locked")}
    assert (D["cue"][["atlas", "cls"]].values == D["locked"][["atlas", "cls"]].values).all()
    d = D["cue"]
    wc = 1 / np.maximum(D["cue"].se.to_numpy(), 1e-3) ** 2
    wl = 1 / np.maximum(D["locked"].se.to_numpy(), 1e-3) ** 2
    # reproduction
    obs = pool(d, d["index"].to_numpy(), wc)
    diffs = np.array([(lambda p: p["immune"] - p["structural"])(pool(d, boots["cue"][:, r], wc)) for r in range(nb)])
    rep = dict(immune=obs["immune"], structural=obs["structural"], diff=obs["immune"] - obs["structural"],
               p_upper=(1 + np.sum(diffs <= 0)) / (1 + nb), classes=len(d), atlases=d.atlas.nunique())
    ref = pd.read_csv(PAPER / "results/figure_5/output/cancer/a85_verdict.csv", index_col=0, header=None)[1]
    pd.DataFrame([dict(quantity=k, reproduced=v, a85_output=float(ref[k])) for k, v in rep.items()]).to_csv(
        out / "reproduction.csv", index=False)
    print(pd.read_csv(out / "reproduction.csv").to_string(), flush=True)
    # H122.1
    ol = pool(D["locked"], D["locked"]["index"].to_numpy(), wl)
    dc, dl = obs["immune"] - obs["structural"], ol["immune"] - ol["structural"]
    atlases = d.atlas.unique()
    ai = pd.Series(range(len(atlases)), index=atlases)[d.atlas].to_numpy()
    rng = np.random.default_rng(SEED)
    stats = []
    while len(stats) < N_ATLAS_BOOT:
        m = np.bincount(rng.integers(0, len(atlases), len(atlases)), minlength=len(atlases))[ai].astype(float)
        try:
            pc = pool(d, d["index"].to_numpy(), wc * m)
            pl = pool(D["locked"], D["locked"]["index"].to_numpy(), wl * m)
        except Exception:
            continue
        v = (pc["immune"] - pc["structural"]) - (pl["immune"] - pl["structural"])
        if np.isfinite(v):
            stats.append(v)
    stats = np.array(stats)
    st = dc - dl
    p1 = (1 + np.sum(stats <= 0)) / (1 + N_ATLAS_BOOT)
    # H122.2
    ps = pd.DataFrame([(j, k, v, x["name"]) for x in results for (j, k, v) in x["pseudo"]],
                      columns=["j", "cls", "index", "atlas"])
    ps = ps.merge(d[["atlas", "cls", "group"]].assign(w=wc), on=["atlas", "cls"])
    pd_ = ps.groupby("j").apply(lambda x: pool(x.assign(group=x.group), x["index"].to_numpy(), x.w.to_numpy()))
    pdiff = np.array([v["immune"] - v["structural"] for v in pd_])
    pdiff_locked_ref = pdiff
    pd.DataFrame(dict(set=range(len(pdiff)), immune=[v["immune"] for v in pd_], structural=[v["structural"] for v in pd_],
                      diff=pdiff)).to_csv(out / "pseudo.csv", index=False)
    q95 = float(np.quantile(pdiff, 0.95))
    p2 = (1 + np.sum(pdiff >= dc)) / (1 + len(pdiff))
    # Holm
    ps_ = {"H122.1": p1, "H122.2": p2}
    order = sorted(ps_, key=ps_.get)
    holm, run = {}, 0.0
    for i, h in enumerate(order):
        run = max(run, min(1, (len(order) - i) * ps_[h]))
        holm[h] = run
    h = pd.DataFrame([
        dict(hypothesis="H122.1", statistic="(imm-str cue) - (imm-str locked)", value=st, cue_diff=dc, locked_diff=dl,
             immune_locked=ol["immune"], structural_locked=ol["structural"], lo=np.quantile(stats, .025),
             hi=np.quantile(stats, .975), p=p1, holm=holm["H122.1"], verdict="PASS" if st > 0 and p1 < 0.05 else "FAIL",
             n_draws=N_ATLAS_BOOT),
        dict(hypothesis="H122.2", statistic="observed imm-str cue diff vs 95th pct of pseudo sets", value=dc, cue_diff=dc,
             locked_diff=np.nan, pseudo_q95=q95, pseudo_median=float(np.median(pdiff)), pseudo_max=float(pdiff.max()),
             pseudo_mean=float(pdiff.mean()), p=p2, holm=holm["H122.2"], verdict="PASS" if dc > q95 else "FAIL",
             n_draws=len(pdiff))])
    h.to_csv(out / "h122.csv", index=False)
    print(h.T.to_string(), flush=True)
    # per atlas
    pa = []
    for a_, x in d.groupby("atlas"):
        xl = D["locked"].loc[x.index]
        r = dict(atlas=a_, classes=len(x), immune_classes=int((x.group == "immune").sum()),
                 structural_classes=int((x.group == "structural").sum()))
        if r["immune_classes"] and r["structural_classes"]:
            pc = pool(x, x["index"].to_numpy(), wc[x.index]); pl = pool(xl, xl["index"].to_numpy(), wl[x.index])
            r.update(immune_cue=pc["immune"], structural_cue=pc["structural"], diff_cue=pc["immune"] - pc["structural"],
                     immune_locked=pl["immune"], structural_locked=pl["structural"],
                     diff_locked=pl["immune"] - pl["structural"])
            r["cue_minus_locked"] = r["diff_cue"] - r["diff_locked"]
            r["sign_cue"], r["sign_cue_minus_locked"] = np.sign(r["diff_cue"]), np.sign(r["cue_minus_locked"])
        pa.append(r)
    pd.DataFrame(pa).to_csv(out / "per_atlas.csv", index=False)
    # composition + depth
    comp = pd.DataFrame([r for x in results for r in x["comp"]])
    cc = comp[comp.set == "cue"].drop(columns="set").rename(columns={"own_state_var": "own_cue_var"})
    cl = comp[comp.set == "locked"][["atlas", "cls", "own_state_var"]].rename(columns={"own_state_var": "own_locked_var"})
    cc = cc.merge(cl, on=["atlas", "cls"]).merge(d[["atlas", "cls", "group", "units"]], on=["atlas", "cls"])
    summ = cc.groupby("group")[["own_cue_var", "own_locked_var", "log_depth", "n_cells"]].median().reset_index()
    summ.insert(0, "cls", "MEDIAN"); summ.insert(0, "atlas", "ALL")
    pd.concat([cc, summ]).to_csv(out / "composition.csv", index=False)
    dp = d[["atlas", "cls", "group", "units", "index"]].rename(columns={"index": "index_cue"}).assign(
        index_locked=D["locked"]["index"].to_numpy()).merge(cc[["atlas", "cls", "log_depth"]], on=["atlas", "cls"])
    from scipy.stats import spearmanr
    rowsd = []
    for a_, x in list(dp.groupby("atlas")) + [("ALL_demeaned", None)]:
        if x is None:
            x = dp.copy()
            for c in ("index_cue", "index_locked", "log_depth"):
                x[c] = x[c] - x.groupby("atlas")[c].transform("mean")
        for s in ("index_cue", "index_locked"):
            rr = spearmanr(x.log_depth, x[s]) if len(x) > 2 else (np.nan, np.nan)
            rowsd.append(dict(scope=a_, kind="spearman_index_vs_logdepth", set=s, n=len(x), value=rr[0], p=rr[1]))
    # difference adjusted for depth: index ~ atlas FE + logdepth + immune (immune vs structural classes)
    x = dp[dp.group.isin(["immune", "structural"])].copy()
    X0 = pd.get_dummies(x.atlas).astype(float)
    for s in ("index_cue", "index_locked"):
        for lab, cols in (("unadjusted", [x.group.eq("immune").astype(float)]),
                          ("adjusted_logdepth", [x.group.eq("immune").astype(float), x.log_depth])):
            M = np.column_stack([X0.to_numpy()] + [c.to_numpy() for c in cols])
            beta = np.linalg.lstsq(M, x[s].to_numpy(), rcond=None)[0]
            res_ = x[s].to_numpy() - M @ beta
            s2 = res_ @ res_ / (len(x) - M.shape[1])
            se = np.sqrt(np.diag(s2 * np.linalg.inv(M.T @ M)))
            rowsd.append(dict(scope="ALL_atlasFE", kind=f"immune_minus_structural_{lab}", set=s, n=len(x),
                              value=beta[X0.shape[1]], se=se[X0.shape[1]]))
            if lab != "unadjusted":
                rowsd.append(dict(scope="ALL_atlasFE", kind="logdepth_slope", set=s, n=len(x), value=beta[-1], se=se[-1]))
    pd.DataFrame(rowsd).to_csv(out / "depth.csv", index=False)
    print("done", flush=True)


if __name__ == "__main__":
    main()
