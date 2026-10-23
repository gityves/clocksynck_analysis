"""A92 and A93: what keeps the organs of a ventilated donor coupled, and does the clock follow a gradient of
pre-death disturbance? Inter-organ coherence and alignment of the organ clocks with the wall time in GTEx.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A92, A93; run notes there). Per tissue, the clock state is the
projection of the 14 clock genes on the fixed GTEx waveform (one projector per tissue, counts, Poisson variance:
`clocksync.paired.projector / states`); the inter-organ coherence C(S) of a donor set S is the mean over tissue pairs
with >= 15 (A92) or >= 20 (A93) shared donors of the Poisson-corrected |complex correlation| of the two states.
  A92  H92.1 phase memory: C(shortest) - C(longest tertile of ventilation hours) > 0, hours permuted within cause of
       death (2,000); H92.2 brain-lesion proxy: C(no primary brain lesion) - C(cerebrovascular, trauma) > 0, labels
       permuted (2,000). Descriptive: glucocorticoid adjustment (GC score regressed out of every state), wall-time
       alignment R = |mean exp(i (phase - w t_death))| in ventilated vs fast donors.
  A93  groups fast (Hardy 1-2), H3, H4, ventilated (Hardy 0) at equal n (larger group subsampled): H93.1 R(H4) - R(vent)
       > 0, labels permuted between H4 and ventilated (1,000); H93.2 the same for C (two-sided); position index;
       organ-specific and brain-body coherence; H93.3 season: host phase minus clock time, summer vs winter deaths of
       fast donors, two-sided permutation (2,000).
Controlled access: aggregates only.

    python compute_coupling_sources.py [--subsample N]
Outputs: paper/results/figure_3/output/gtex_cues/a92_{summary,walltime}.csv, a93_{groups,summary}.csv
Legacy: legacy/desync/validation/gtex_covariates/scripts/{a92_coupling_sources,a93_disturbance}.py
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _covariates as K  # noqa: E402
import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
import importlib  # noqa: E402
P = importlib.import_module("clocksync.paired")  # the module (the package namespace holds the function `paired`)
from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

W = S.W
BODY = [t for t in S.TISSUES if t != "Pituitary" and not t.startswith("Brain")]
EXTRA = ["Liver", "Adrenal Gland"]
WF = cs.gtex()


def sub_donors(pb, don, n, seed):
    if not n:
        return pb, don
    keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:n]
    return pb[pb.donor.isin(keep)].reset_index(drop=True), don.loc[don.index.isin(keep)]


def tissue_state(rows: pd.DataFrame):
    """(b, v): clock-plane state and Poisson variance of one tissue's rows, one projector for the tissue."""
    x = rows[list(WF.genes)]
    L = P.projector(x, WF, counts=True)
    return P.states(x, L, WF, depth=rows.depth.to_numpy(float))


# ---- A92 ---------------------------------------------------------------------------------------------------------------

def states92(pb, don):
    """tissue -> cohort -> DataFrame(index donor; b0, b1, v, gc, tod): states centred within cohort, GC score = mean
    within-cohort z of log CPM of the six GC genes."""
    out = {}
    for t in BODY:
        r = pb[pb.compartment == t].reset_index(drop=True)
        b, v = tissue_state(r)
        lg = np.log(r[list(S.GC)].to_numpy(float) / r.depth.to_numpy(float)[:, None] * 1e6 + 0.5)
        coh = np.where(don.hardy.reindex(r.donor).to_numpy() == 0, "vent", "fast")
        out[t] = {}
        for c in ("vent", "fast"):
            m = coh == c
            z = (lg[m] - lg[m].mean(0)) / lg[m].std(0)
            bb = b[m] - b[m].mean(0)
            x = pd.DataFrame(dict(b0=bb[:, 0], b1=bb[:, 1], v=v[m], gc=z.mean(1),
                                  tod=don.tod.reindex(r.donor[m]).to_numpy(float)), index=r.donor[m].to_numpy())
            out[t][c] = x[~x.index.duplicated()]
    return out


def gc_adjust(x: pd.DataFrame) -> pd.DataFrame:
    c = x.b0.to_numpy() + 1j * x.b1.to_numpy()
    D = np.column_stack([np.ones(len(x)), x.gc.to_numpy()])
    r = c - D @ np.linalg.lstsq(D, c, rcond=None)[0]
    return x.assign(b0=r.real, b1=r.imag)


def coherence92(st: dict, donors, min_pair=15) -> float:
    donors = set(donors)
    vals = []
    for t1, t2 in itertools.combinations(st, 2):
        a, b = st[t1], st[t2]
        common = sorted(set(a.index) & set(b.index) & donors)
        if len(common) < min_pair:
            continue
        a, b = a.loc[common], b.loc[common]
        vals.append(abs(P.ccorr(a[["b0", "b1"]].to_numpy(), a.v.to_numpy(), b[["b0", "b1"]].to_numpy(), b.v.to_numpy())))
    return float(np.mean(vals)) if vals else np.nan


def a92(pb, don, ph, run, rng) -> tuple[pd.Series, pd.DataFrame]:
    n_perm, n_boot = (50, 20) if run["quick"] else (2000, 500)
    st = states92(pb, don)
    vent = {t: st[t]["vent"] for t in BODY}
    fast = {t: st[t]["fast"] for t in BODY}
    donors = sorted(set().union(*[set(x.index) for x in vent.values()]) & set(ph.index[ph.hours.notna()]))
    hrs = ph.hours.reindex(donors).to_numpy(float)
    cause = ph.cause.reindex(donors).to_numpy()
    dn = np.asarray(donors)
    tert = pd.qcut(hrs, 3, labels=False)
    out = {"C_fast": coherence92(fast, set().union(*[set(x.index) for x in fast.values()])),
           "C_vent": coherence92(vent, dn)}
    for k in range(3):
        out[f"C_tertile{k + 1}"] = coherence92(vent, dn[tert == k])
        out[f"days_median_tertile{k + 1}"] = float(np.median(hrs[tert == k]) / 24)
    d1 = lambda lab: coherence92(vent, dn[lab == 0]) - coherence92(vent, dn[lab == 2])
    obs1 = d1(tert)
    null1 = np.array([d1(_perm_within(tert, cause, rng)) for _ in range(n_perm)])
    out.update(H92_1_diff=obs1, H92_1_p=float((1 + np.sum(null1 >= obs1)) / (1 + n_perm)))
    lesion = np.isin(cause, ["cerebrovascular", "trauma"]).astype(int)
    d2 = lambda lab: coherence92(vent, dn[lab == 0]) - coherence92(vent, dn[lab == 1])
    obs2 = d2(lesion)
    null2 = np.array([d2(rng.permutation(lesion)) for _ in range(n_perm)])
    out.update(C_no_brain_lesion=coherence92(vent, dn[lesion == 0]), C_brain_lesion=coherence92(vent, dn[lesion == 1]),
               n_no_brain_lesion=int((lesion == 0).sum()), n_brain_lesion=int(lesion.sum()), H92_2_diff=obs2,
               H92_2_p=float((1 + np.sum(null2 >= obs2)) / (1 + n_perm)))
    for c in sorted(set(cause)):
        if (cause == c).sum() >= 30:
            out[f"C_cause_{c}"] = coherence92(vent, dn[cause == c])
    vadj = {t: gc_adjust(x) for t, x in vent.items()}
    out["C_vent_gc_adjusted"] = coherence92(vadj, dn)
    bs = []
    for _ in range(n_boot):
        s = set(rng.choice(dn, len(dn)))  # a set, as the legacy run: duplicates collapse (a 63 % subsample)
        bs.append(coherence92(vent, s) - coherence92(vadj, s))
    out.update(gc_drop=out["C_vent"] - out["C_vent_gc_adjusted"], gc_drop_lo=float(np.quantile(bs, 0.025)),
               gc_drop_hi=float(np.quantile(bs, 0.975)))
    wt = []
    for t in BODY:
        for c, x in (("vent", vent[t]), ("fast", fast[t])):
            x = x.dropna(subset=["tod"])
            ang = np.arctan2(x.b1, x.b0) - W * x.tod
            wt.append(dict(tissue=t, cohort=c, n=len(x), R=float(abs(np.exp(1j * ang).mean()))))
    wt = pd.DataFrame(wt)
    out["R_walltime_vent_median"] = float(wt[wt.cohort == "vent"].R.median())
    out["R_walltime_fast_median"] = float(wt[wt.cohort == "fast"].R.median())
    out["H92_1"] = "PASS" if (obs1 > 0 and out["H92_1_p"] < 0.05) else "FAIL"
    out["H92_2"] = "PASS" if (obs2 > 0 and out["H92_2_p"] < 0.05) else "FAIL"
    return pd.Series(out), wt


def _perm_within(labels, groups, rng):
    out = np.asarray(labels).copy()
    for g in np.unique(groups):
        i = np.flatnonzero(groups == g)
        out[i] = np.asarray(labels)[rng.permutation(i)]
    return out


# ---- A93 ---------------------------------------------------------------------------------------------------------------
GROUPS = {"fast": [1, 2], "H3": [3], "H4": [4], "vent": [0]}
MIN_N = 20


def states93(pb, don, tissues):
    """group -> tissue -> (donors, centred complex state, Poisson variance, tod)."""
    out = {g: {} for g in GROUPS}
    for t in tissues:
        r = pb[pb.compartment == t].reset_index(drop=True)
        if len(r) < 20:
            continue
        b, v = tissue_state(r)
        h = don.hardy.reindex(r.donor).to_numpy()
        for g, hd in GROUPS.items():
            m = np.isin(h, hd)
            subj = r.donor.to_numpy()[m]
            _, first = np.unique(subj, return_index=True)
            c = (b[m, 0] + 1j * b[m, 1])[first]
            out[g][t] = (subj[first], c - c.mean(), v[m][first], don.tod.reindex(subj[first]).to_numpy(float))
    return out


def align_R(st, donors, tissues):
    r = []
    for t in tissues:
        if t not in st:
            continue
        s, c, _, tod = st[t]
        m = np.isin(s, donors)
        if m.sum() >= MIN_N:
            r.append(abs(np.exp(1j * (np.angle(c[m]) - W * tod[m])).mean()))
    return float(np.median(r)) if r else np.nan


def pair_corr(st, t1, t2, donors):
    """Poisson-corrected |complex correlation| (local: the legacy A93 form, normalised by the corrected variances)."""
    s1, c1, v1, _ = st[t1]
    s2, c2, v2, _ = st[t2]
    common = np.intersect1d(np.intersect1d(s1, s2), donors)
    if len(common) < MIN_N:
        return np.nan
    i1, i2 = np.searchsorted(s1, common), np.searchsorted(s2, common)
    a, b = c1[i1] - c1[i1].mean(), c2[i2] - c2[i2].mean()
    k = len(common) - 1
    ea, eb = np.sum(abs(a) ** 2) / k - v1[i1].mean(), np.sum(abs(b) ** 2) / k - v2[i2].mean()
    return float(abs(np.sum(a * np.conj(b)) / k) / np.sqrt(max(ea, 1e-12) * max(eb, 1e-12)))


def coherence93(st, donors, set1, set2=None):
    set1 = [t for t in set1 if t in st]
    set2 = None if set2 is None else [t for t in set2 if t in st]
    pairs = itertools.combinations(set1, 2) if set2 is None else itertools.product(set1, set2)
    v = [pair_corr(st, a, b, donors) for a, b in pairs]
    v = [x for x in v if np.isfinite(x)]
    return float(np.mean(v)) if v else np.nan


def donors_of(st):
    return np.unique(np.concatenate([s for s, *_ in st.values()]))


def at_n(fn, st, donors, n, rng, n_sub):
    if len(donors) <= n:
        return fn(st, donors)
    return float(np.nanmedian([fn(st, rng.choice(donors, n, replace=False)) for _ in range(n_sub)]))


def season_test(pbf, donf, ph, rng, n_perm) -> dict:
    f = S.host_states(pbf, S.TISSUES, name="fast")
    z = pd.Series(np.exp(1j * (np.arctan2(f.y, f.x) - W * donf.tod.reindex(f.index))), index=f.index)
    d = pd.DataFrame({"z": z, "season": ph.season.reindex(f.index)}).dropna()
    mean_h = lambda zz: float(np.angle(np.mean(zz)) / W)
    lab = d.season.to_numpy()
    zz = d.z.to_numpy()
    stat = lambda lb: (mean_h(zz[lb == "Summer"]) - mean_h(zz[lb == "Winter"]) + 12) % 24 - 12
    obs = stat(lab)
    null = np.array([stat(rng.permutation(lab)) for _ in range(n_perm)])
    out = dict(summer_minus_winter_h=obs, p_two=float((1 + np.sum(np.abs(null) >= abs(obs))) / (1 + len(null))))
    for s in ("Spring", "Summer", "Fall", "Winter"):
        q = zz[lab == s]
        out[f"n_{s}"], out[f"offset_{s}_h"], out[f"R_{s}"] = int(len(q)), mean_h(q), float(abs(np.mean(q)))
    return out


def a93(pb, don, ph, pbf, donf, run, rng):
    n_perm, n_sub, n_sub_perm = (30, 5, 2) if run["quick"] else (1000, 20, 5)
    brain = sorted(t for t in pb.compartment.unique() if t.startswith("Brain"))
    st = states93(pb, don, BODY + EXTRA + brain)
    dn = {g: donors_of(st[g]) for g in GROUPS}
    n4 = len(dn["H4"])
    rows = []
    for g in GROUPS:
        n = min(len(dn[g]), n4)
        row = dict(group=g, donors=len(dn[g]), n_used=n,
                   R=at_n(lambda s, d: align_R(s, d, BODY), st[g], dn[g], n, rng, n_sub),
                   C=at_n(lambda s, d: coherence93(s, d, BODY), st[g], dn[g], n, rng, n_sub))
        for t in ["Liver", "Adrenal Gland", "Muscle - Skeletal"]:
            row[f"C_{t}"] = at_n(lambda s, d: coherence93(s, d, [t], [x for x in BODY if x != t]), st[g], dn[g], n, rng,
                                 n_sub)
        if g != "vent":
            row["C_brain_body"] = at_n(lambda s, d: coherence93(s, d, brain, BODY), st[g], dn[g], n, rng, n_sub)
            row["C_brain_brain"] = at_n(lambda s, d: coherence93(s, d, brain), st[g], dn[g], n, rng, n_sub)
        rows.append(row)
        print(row, flush=True)
    G = pd.DataFrame(rows).set_index("group")
    pool = {}
    for t in BODY:
        if t in st["H4"] and t in st["vent"]:
            cat = tuple(np.concatenate([st["H4"][t][i], st["vent"][t][i]]) for i in range(4))
            o = np.argsort(cat[0])
            pool[t] = tuple(x[o] for x in cat)
    allp = np.union1d(dn["H4"], dn["vent"])
    lab4 = np.isin(allp, dn["H4"])

    def diffs(l4):
        d4, dv = allp[l4], allp[~l4]
        sub = [rng.choice(dv, len(d4), replace=False) for _ in range(n_sub_perm)]
        r = align_R(pool, d4, BODY) - np.median([align_R(pool, s, BODY) for s in sub])
        c = coherence93(pool, d4, BODY) - np.median([coherence93(pool, s, BODY) for s in sub])
        return r, c

    obs_r = G.loc["H4", "R"] - G.loc["vent", "R"]
    obs_c = G.loc["H4", "C"] - G.loc["vent", "C"]
    null = np.array([diffs(rng.permutation(lab4)) for _ in range(n_perm)])
    v = dict(R_fast=G.loc["fast", "R"], R_H3=G.loc["H3", "R"], R_H4=G.loc["H4", "R"], R_vent=G.loc["vent", "R"],
             H93_1_diff=obs_r, H93_1_p=float((1 + np.sum(null[:, 0] >= obs_r)) / (1 + n_perm)),
             position_index=(G.loc["H4", "R"] - G.loc["vent", "R"]) / (G.loc["fast", "R"] - G.loc["vent", "R"]),
             H93_2_diff=obs_c, H93_2_p=float((1 + np.sum(np.abs(null[:, 1]) >= abs(obs_c))) / (1 + n_perm)))
    v["H93_1"] = "PASS" if (obs_r > 0 and v["H93_1_p"] < 0.05) else "FAIL"
    v.update({f"season_{k}": x for k, x in season_test(pbf, donf, ph, rng, 200 if run["quick"] else 2000).items()})
    return G, pd.Series(v)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(K.FOLDER, run["quick"])
    ph = K.protected()
    genes = list(WF.genes) + list(S.GC)
    pb92, don92 = gtex.load(tissues=BODY, hardy=(0, 1, 2), genes=genes, depth="library_size", require_tod=False)
    pb92, don92 = sub_donors(pb92, don92, args.subsample, args.seed)
    s92, wt = a92(pb92, don92, ph, run, np.random.default_rng(92))
    s92.to_csv(out / "a92_summary.csv", header=False)
    wt.to_csv(out / "a92_walltime.csv", index=False)
    print(s92.to_string(), "\n", wt.round(2).to_string(index=False), flush=True)

    pb93, don93 = gtex.load(hardy=(0, 1, 2, 3, 4), genes="clock", depth="library_size", require_tod=True)
    pb93, don93 = sub_donors(pb93, don93, args.subsample, args.seed)
    pbf, donf = K.cohort((1, 2), tissues=S.TISSUES, subsample=args.subsample, seed=args.seed, require_tod=False)
    G, v = a93(pb93, don93, ph, pbf, donf, run, np.random.default_rng(93))
    G.to_csv(out / "a93_groups.csv")
    v.to_csv(out / "a93_summary.csv", header=False)
    print(G.round(3).to_string(), "\n", v.to_string(), flush=True)


if __name__ == "__main__":
    main()
