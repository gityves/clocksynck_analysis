"""A48: ageing with known time in the mouse (Wolff 2023, GSE201207): do old tissues stay synchronised?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A48). 6 tissues x 3 ages (6, 18, 27 months) x 12 circadian times; a unit
is a pool of 2 mice giving all six tissues. Human GTEx waveform, host anchored on the five peripheral tissues (adrenal,
heart, kidney, lung, muscle). Order as legacy: timed truth, gates, untimed estimates, hypotheses.
  Timed truth: per tissue and age, adjusted R^2 of the clock-plane state on [1, cos, sin] of the true circadian time
    (CT folded to 24 h) and the timed phase (cosinor of the clock genes on the waveform, A36).
  H48.1 (validation): at each age, circular (Jammalamadaka-SenGupta) correlation of the untimed tissue offsets with the
    timed tissue phases > 0, p < 0.05, in >= 2 of 3 ages.
  H48.2 (ageing): peripheral-module change 27 - 6 months, timed (adjusted R^2) and untimed (d* at the common host
    spread, one pooled young + old host, age is the stratum); one-sided permutation p (age labels permuted within CT).
  H48.3 (hypothalamus): coupling and timed R^2 of the hypothalamus against the peripheral median at each age.
MDE = 2.486 x SD of the permutation null.

    python compute_mouse_ageing.py [--subsample N]      (--subsample N: number of permutations of H48.2)
Outputs (paper/results/figure_4/output/resilience/): a48_{tissue_qc,gates,timed_age,untimed_age,h481,h482,h482_tissues,
h483}.csv
Legacy: legacy/desync/validation/mouse_ageing/scripts/a48_wolff.py (+ common.py).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "resilience"
W = 2 * np.pi / 24
PERIPH = ["adrenal", "heart", "kidney", "lung", "muscle"]
AGES = ["young", "aged", "old"]
Z_MDE = 1.645 + 0.842
WF = cs.gtex()
GENES = list(WF.genes)


def load_wolff():
    """Standard table; built from the public counts when missing (paper/common/datasets/mouse_ageing_wolff/build.py)."""
    name = "mouse_ageing_wolff"
    try:
        return data.load(name), data._builder(name)
    except FileNotFoundError:
        mod = data._builder(name)
        pb, donors = mod.build()
        data.save(name, pb, donors)
        return data.load(name), mod


# ---- local estimators of the timed truth (local: not in clocksync) -----------------------------------------------

def design(hours):
    t = np.asarray(hours, float) * W
    return np.c_[np.ones(len(t)), np.cos(t), np.sin(t)]


def adj_r2(beta, hours) -> float:
    """Adjusted R^2 of the 2-d clock-plane state on [1, cos, sin] of the true time (both coordinates pooled)."""
    D = design(hours)
    n = len(D)
    if n < 6:
        return np.nan
    B = np.linalg.lstsq(D, beta, rcond=None)[0]
    rss = ((beta - D @ B) ** 2).sum()
    tss = ((beta - beta.mean(0)) ** 2).sum()
    return float(1 - (rss / tss) * (n - 1) / (n - 3))


def timed_gain(sub: pd.DataFrame, hours):
    """Cosinor of log(CPM + 0.5) of the clock genes with median CPM >= 5, complex gain on the waveform (A36)."""
    cpm = sub[GENES].to_numpy(float) / sub.depth.to_numpy(float)[:, None] * 1e6
    ex = np.median(cpm, 0) >= 5
    if ex.sum() < 3:
        return np.nan + 0j
    B = np.linalg.lstsq(design(hours), np.log(cpm[:, ex] + 0.5), rcond=None)[0]
    x = (WF.amp * np.exp(1j * WF.psi))[ex]
    return (np.conj(x) @ (B[1] + 1j * B[2])) / (np.conj(x) @ x)


def states(ds: cs.Dataset) -> pd.DataFrame:
    """Clock-plane state of every (donor, compartment) row of the Dataset (the public `Dataset.data` arrays)."""
    d = ds.data
    return pd.DataFrame(d.beta, columns=["x", "y"], index=pd.MultiIndex.from_arrays([d.resection, d.compartment],
                                                                                    names=["donor", "compartment"]))


def timed_truth(ds, pb, ct, groups, comps) -> pd.DataFrame:
    b = states(ds)
    p = pb.set_index(["donor", "compartment"])
    rows = []
    for k in comps:
        for g in sorted(groups.unique()):
            donors = groups.index[groups == g]
            idx = [(u, k) for u in donors if (u, k) in b.index]
            if not idx:
                continue
            h = ct.reindex([u for u, _ in idx]).to_numpy(float)
            gain = timed_gain(p.loc[idx], h)
            rows.append(dict(tissue=k, group=g, n=len(idx), r2=adj_r2(b.loc[idx].to_numpy(), h),
                             phase=-np.angle(gain), g_abs=abs(gain)))
    return pd.DataFrame(rows)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_gate = 100 if run["quick"] else 2000
    n_perm = args.subsample if args.subsample else 1000
    (pb, donors), mod = load_wolff()
    q = mod.tissue_qc(pb)
    q.to_csv(out / "a48_tissue_qc.csv")
    pb = pb[pb.compartment.isin(q.index[q.kept])]
    comps = sorted(pb.compartment.unique())
    pb, bad = mod.drop_units(pb, PERIPH)
    age = donors.age.reindex(sorted(pb.donor.unique()))
    ct = (donors.ct % 24).reindex(age.index)
    print("tissues", comps, "units dropped", bad, flush=True)

    ds_age = {a: cs.Dataset(pb[pb.donor.isin(age.index[age == a])], anchors=PERIPH, name=a) for a in AGES}
    timed = pd.concat([timed_truth(ds_age[a], pb, ct, pd.Series(a, index=age.index[age == a]), comps).assign(age=a)
                       for a in AGES])
    timed.to_csv(out / "a48_timed_age.csv", index=False)

    gates = []
    for a in AGES:
        for label, cc in ((f"{a} (5 peripheral)", PERIPH), (f"{a} (all 6, descriptive)", comps)):
            sub = pb[pb.donor.isin(age.index[age == a]) & pb.compartment.isin(cc)]
            g = cs.Dataset(sub, anchors=PERIPH).gate(n_perm=n_gate, seed=0)
            gates.append(dict(host=label, clock_share=g.clock_share, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase,
                              passes=bool(g.passes)))
    yo = pb[pb.donor.isin(age.index[age.isin(["young", "old"])])]
    g = cs.Dataset(yo[yo.compartment.isin(PERIPH)], anchors=PERIPH, strata=age.to_dict()).gate(n_perm=n_gate, seed=0)
    gates.append(dict(host="young+old (5 peripheral)", clock_share=g.clock_share, p_pairing=g.p_pairing,
                      p_acrophase=g.p_acrophase, passes=bool(g.passes)))
    gates = pd.DataFrame(gates)
    gates.to_csv(out / "a48_gates.csv", index=False)
    print(gates.round(4).to_string(index=False), flush=True)

    rows, h481, h483 = [], [], []
    for a in AGES:
        f = ds_age[a].fit(n_boot=0, gate=False).table
        T = timed[timed.age == a].set_index("tissue")
        tab = pd.DataFrame(dict(age=a, tissue=f.index, anchor=f.anchor.to_numpy(), coupling=f.coupling.to_numpy(),
                                offset_h=f.offset_h.to_numpy(), boundary=f.boundary.to_numpy(), n=f.n.to_numpy(),
                                timed_r2=T.r2.reindex(f.index).to_numpy(), timed_phase=T.phase.reindex(f.index).to_numpy()))
        rows.append(tab)
        gp = bool(gates[gates.host == f"{a} (5 peripheral)"].passes.iloc[0])
        for tag, m in ((f"H48.1 {a} all 6 tissues", np.ones(len(tab), bool)),
                       (f"H48.1 {a} peripheral 5 (sensitivity)", tab.anchor.to_numpy())):
            r = cs.timed.agreement(tab.offset_h[m] * W, tab.timed_phase[m], n_perm=10000)
            h481.append(dict(test=tag, age=a, gate=gp, n=int(r.n), r=r.r_js, p=r.p_js, r_fl=r.r_fl,
                             offset_h=r.offset_h, median_abs_err_h=r.median_abs_err_h))
        hyp, pr = tab[tab.tissue == "hypothalamus"].iloc[0], tab[tab.anchor]
        h483.append(dict(age=a, gate=gp, c_hypothalamus=hyp.coupling, c_periph_median=pr.coupling.median(),
                         diff_untimed=hyp.coupling - pr.coupling.median(), r2_hypothalamus=hyp.timed_r2,
                         r2_periph_median=pr.timed_r2.median(), diff_timed=hyp.timed_r2 - pr.timed_r2.median(),
                         c_rank_of_6=int((tab.coupling >= hyp.coupling).sum()),
                         r2_rank_of_6=int((tab.timed_r2 >= hyp.timed_r2).sum())))
    pd.concat(rows).to_csv(out / "a48_untimed_age.csv", index=False)
    h481 = pd.DataFrame(h481)
    h481.to_csv(out / "a48_h481.csv", index=False)
    pd.DataFrame(h483).to_csv(out / "a48_h483.csv", index=False)
    prim = h481[h481.test.str.contains("all 6")]
    print(h481.round(4).to_string(index=False))
    print("H48.1 criterion (r > 0 and p < 0.05 in >= 2 of 3 ages):",
          "PASS" if int(((prim.r > 0) & (prim.p < 0.05)).sum()) >= 2 else "FAIL", flush=True)

    # ---- H48.2: young vs old on one pooled host (age = stratum), labels permuted within CT -------------------
    ds = cs.Dataset(yo, anchors=PERIPH, strata=age.to_dict(), name="young+old")
    lab0 = age[age.isin(["young", "old"])]
    ctu = ct.reindex(lab0.index)
    b = states(ds)
    pyo = pb.set_index(["donor", "compartment"])
    sets = {k: [u for u in lab0.index if (u, k) in b.index] for k in comps}

    def stat(lab):
        c = ds.contrast(lab, levels=("young", "old"), n_perm=0, n_boot=0).table
        out_ = {}
        for k in comps:
            r2 = {}
            for g in ("young", "old"):
                idx = [(u, k) for u in sets[k] if lab[u] == g]
                r2[g] = adj_r2(b.loc[idx].to_numpy(), ctu.reindex([u for u, _ in idx]).to_numpy(float))
            out_[k] = dict(d=c.d.get(k, np.nan), dstar=c.d_star.get(k, np.nan), dr2=r2["old"] - r2["young"])
        return pd.DataFrame(out_).T

    obs = stat(lab0)
    obs.rename_axis("tissue").to_csv(out / "a48_h482_tissues.csv")
    mo = obs.loc[PERIPH].mean()
    rng = np.random.default_rng(0)
    null = []
    for _ in range(n_perm):
        new = lab0.copy()
        for _, ix in lab0.groupby(ctu).groups.items():
            new[ix] = rng.permutation(lab0[ix].to_numpy())
        null.append(stat(new).loc[PERIPH].mean())
    null = pd.DataFrame(null)
    summ = []
    for m, nm in (("dr2", "timed: adjusted R^2 (27 - 6 mo)"), ("dstar", "untimed: d* (27 - 6 mo)"),
                  ("d", "untimed: d (27 - 6 mo, descriptive)")):
        x = null[m].to_numpy()
        x = x[np.isfinite(x)]
        summ.append(dict(statistic=nm, observed=mo[m], null_mean=x.mean(), null_sd=x.std(ddof=1),
                         p_lower=(1 + np.sum(x <= mo[m])) / (1 + len(x)), n_perm=len(x), MDE=Z_MDE * x.std(ddof=1)))
    summ = pd.DataFrame(summ)
    summ.to_csv(out / "a48_h482.csv", index=False)
    print(summ.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
