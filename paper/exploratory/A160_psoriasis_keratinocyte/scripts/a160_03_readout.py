"""A160 step 3: readout. Run only after a160_02_blind_power.py (the variants that qualify are read from
output/a160_blind_power.csv).

H160.1 (kc) and H160.2 (bins): PP minus PN within-group core-minus-cue offset (h), one-sided exact within-pair label swaps
(all 2^n pair labelings, the identity included: p = share of labelings with a statistic <= the observed), Holm over the
variants that qualify (MDE <= 1.7 h). Non-qualifying variants: estimate and swap p reported as descriptive.
Descriptives: PN minus healthy (trait, as A115; bootstrap interval, libraries resampled within group), whole-skin
pseudobulk of the same libraries, core offset alone (paired arg ccorr of the core states, PP relative to PN), NR1D1 leave-out
(core waveform without NR1D1; cue set unchanged).

    python a160_03_readout.py
Outputs: output/a160_tests.csv, a160_descriptives.csv, a160_per_bin.csv
"""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import _lib as L

X, obs = L.load()
plan = json.load(open(L.OUT / "a160_plan.json"))
power = pd.read_csv(L.OUT / "a160_blind_power.csv").drop_duplicates("variant").set_index("variant")
status = lambda k: obs.loc[k, "status"]
BINU = [u for u in plan if u != "kc"]
SEED = 1601


def unit_states(u, pl, drop=None):
    """states of the PP rows then PN rows of the pairs of unit u."""
    n = len(pl["pairs"]["pair"])
    return L.states(X, pl["pairs"]["case"] + pl["pairs"]["ref"], pl["genes"], drop), n


def rows_flip(n, f):
    i = np.arange(n); return np.where(f, n + i, i), np.where(f, i, n + i)


def stat(prep, flipmap, kind="diff"):
    """Pooled statistic for a variant: prep = {unit: (states, n, pairs)}; flipmap: pair -> bool."""
    d, se = [], []
    for u, (st, n, pairs) in prep.items():
        f = np.array([flipmap.get(p, False) for p in pairs])
        ia, ib = rows_flip(n, f)
        if kind == "diff":
            d.append(L.diff(st, ia, ib))
            se.append(L.boot_se(st, ia, ib, np.random.default_rng(SEED)) if len(prep) > 1 else np.nan)
        else:  # core offset alone, paired (PP relative to PN)
            d.append(L.paired_core(st["core"], ia, ib)); se.append(np.nan)
    if len(prep) == 1:
        return d[0], d, se
    d, se = np.array(d), np.array(se); ok = np.isfinite(se) if kind == "diff" else np.ones(len(d), bool)
    return (L.wmean(d[ok], se[ok]) if kind == "diff" else L.wmean(d, np.ones(len(d)))), d, se


def swap_test(prep, kind="diff"):
    pairs = sorted({p for _, (_, _, pr) in prep.items() for p in pr})
    obs_v, d, se = stat(prep, {}, kind)
    null = np.array([stat(prep, dict(zip(pairs, f)), kind)[0] for f in L.all_flips(len(pairs))])
    p = float(np.mean(null <= obs_v + 1e-9))
    return obs_v, p, len(null), len(pairs), d, se


def prep_of(units, drop=None, plan_=None):
    plan_ = plan_ or plan
    out = {}
    for u in units:
        st, n = unit_states(u, plan_[u], drop)
        out[u] = (st, n, plan_[u]["pairs"]["pair"])
    return out


def plan_group(group):
    """The kc pairs / pool re-pointed to another group of the same libraries (whole skin)."""
    pl = plan["kc"]; sw = lambda k: k.split("|")[0] + "|" + group
    pool = [sw(k) for k in pl["pool"]]
    pr = dict(pair=pl["pairs"]["pair"], case=[sw(k) for k in pl["pairs"]["case"]], ref=[sw(k) for k in pl["pairs"]["ref"]])
    return dict(pairs=pr, pool=pool, genes=L.cue_genes(X, pool))


def trait(units, plan_, label, n_boot=1000):
    """PN minus healthy within-group core-cue offset (h) with a bootstrap interval; bins: weighted over the eligible bins."""
    rng = np.random.default_rng(SEED)
    per = {}
    for u in units:
        pl = plan_[u]
        pn = [k for k in pl["pool"] if status(k) == "nonlesional"]; hl = [k for k in pl["pool"] if status(k) == "healthy"]
        st = L.states(X, pn + hl, pl["genes"]); ia, ib = np.arange(len(pn)), len(pn) + np.arange(len(hl))
        per[u] = (st, ia, ib)
    val = lambda ps: (L.diff(ps[0], ps[1], ps[2]) if len(units) == 1 else None)
    if len(units) == 1:
        u = units[0]; st, ia, ib = per[u]; est = L.diff(st, ia, ib)
        b = [L.diff(st, rng.choice(ia, len(ia)), rng.choice(ib, len(ib))) for _ in range(n_boot)]
    else:
        def comb(rs):
            d, se = [], []
            for u in units:
                st, ia, ib = per[u]
                a, bb = (ia, ib) if rs is None else (rs[0][u], rs[1][u])
                d.append(L.diff(st, a, bb)); se.append(L.boot_se(st, ia, ib, np.random.default_rng(SEED), 50))
            return L.wmean(d, se)
        est = comb(None)
        b = [comb(({u: rng.choice(per[u][1], len(per[u][1])) for u in units}, {u: rng.choice(per[u][2], len(per[u][2])) for u in units})) for _ in range(n_boot // 5)]
    b = np.asarray(L.S.wrap(np.array(b) - est))
    n_pn = len([k for k in plan_[units[0]]["pool"] if status(k) == "nonlesional"]); n_h = len(plan_[units[0]]["pool"]) - n_pn
    return dict(item=label, estimate_h=est, lo_h=est + np.quantile(b, .025), hi_h=est + np.quantile(b, .975), n_boot=len(b), n_PN=n_pn, n_healthy=n_h)


tests, desc, per_bin = [], [], []
variants = {"kc": ["kc"], "bins": BINU}
for v, units in variants.items():
    prep = prep_of(units)
    est, p, nperm, npair, d, se = swap_test(prep)
    tests.append(dict(variant=v, hypothesis={"kc": "H160.1", "bins": "H160.2"}[v], estimate_h=est, n_pairs=npair, n_labelings=nperm, p_one_sided=p,
                      mde_gate_h=power.loc[v, "mde_gate_h"], confirmatory=bool(power.loc[v, "confirmatory"])))
    if v == "bins":
        for u, dd, ss in zip(units, d, se):
            per_bin.append(dict(unit=u, pairs=len(plan[u]["pairs"]["pair"]), diff_h=dd, boot_se_h=ss))
    # descriptives on the same variant: core offset alone, NR1D1 leave-out
    est_c, p_c, *_ = swap_test(prep, "core")
    desc.append(dict(variant=v, item="core_offset_alone_PP_vs_PN_h", estimate_h=est_c, p_one_sided=p_c))
    est_n, p_n, *_ = swap_test(prep_of(units, drop=["NR1D1"]))
    desc.append(dict(variant=v, item="core_cue_without_NR1D1_h", estimate_h=est_n, p_one_sided=p_n))
T = pd.DataFrame(tests)
q = T.confirmatory
T["p_holm"] = np.nan
if q.any():
    ps = T.loc[q, "p_one_sided"].to_numpy(); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0
    for r, i in enumerate(o):
        run = max(run, (m - r) * ps[i]); adj[i] = min(1.0, run)
    T.loc[q, "p_holm"] = adj
T["H160_pass"] = np.where(T.confirmatory, (T.estimate_h < 0) & (T.p_holm < 0.05), np.nan)
T.to_csv(L.OUT / "a160_tests.csv", index=False)
pd.DataFrame(per_bin).to_csv(L.OUT / "a160_per_bin.csv", index=False)

# whole-skin pseudobulk of the same libraries (kc pairs), and of all complete pairs
pw = plan_group("whole")
w_prep = {"whole": (L.states(X, pw["pairs"]["case"] + pw["pairs"]["ref"], pw["genes"]), len(pw["pairs"]["pair"]), pw["pairs"]["pair"])}
est, p, nperm, npair, *_ = swap_test(w_prep)
desc.append(dict(variant="whole_skin_same_libraries", item="core_cue_PP_minus_PN_h", estimate_h=est, p_one_sided=p, n_pairs=npair))
est, p, *_ = swap_test(w_prep, "core")
desc.append(dict(variant="whole_skin_same_libraries", item="core_offset_alone_PP_vs_PN_h", estimate_h=est, p_one_sided=p, n_pairs=npair))
prw, poolw, _ = L.pairs_for(X, obs, "whole", 0)
pl_all = dict(pairs=prw.to_dict("list"), pool=poolw, genes=L.cue_genes(X, poolw))
a_prep = prep_of(["whole"], plan_={"whole": pl_all})
est, p, nperm, npair, *_ = swap_test(a_prep)
desc.append(dict(variant="whole_skin_all_complete_pairs", item="core_cue_PP_minus_PN_h", estimate_h=est, p_one_sided=p, n_pairs=npair))
est, p, *_ = swap_test(a_prep, "core")
desc.append(dict(variant="whole_skin_all_complete_pairs", item="core_offset_alone_PP_vs_PN_h", estimate_h=est, p_one_sided=p, n_pairs=npair))
# trait PN - healthy
tr = [dict(variant="kc", **trait(["kc"], plan, "PN_minus_healthy_h")),
      dict(variant="bins", **trait(BINU, plan, "PN_minus_healthy_h")),
      dict(variant="whole_skin_same_libraries", **trait(["kc"], {"kc": pw}, "PN_minus_healthy_h"))]
D = pd.concat([pd.DataFrame(desc), pd.DataFrame(tr)], ignore_index=True)
D.to_csv(L.OUT / "a160_descriptives.csv", index=False)
print(T.round(4).to_string(index=False)); print(pd.DataFrame(per_bin).round(3).to_string(index=False)); print(D.round(3).to_string(index=False))
