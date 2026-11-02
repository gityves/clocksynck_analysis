"""A195 scoring (registered_experiments.md, A195, commit afc2515): is the proliferation clock-signal loss confined to cycling
cells? L1 clock-plane variance (`clocksync.decompose.components` on Poisson-corrected log fractions) of malignant pseudobulks
from A187 (bins NC / S / G2M per patient; cycling = S + G2M, ALL = NC + S + G2M).

Order (stop rules binding): signal check on NC -> POWER (MDE = 2.49 x null SD; a test with MDE > 0.40 is not read) ->
H195.1 (cycling - NC, within patient), H195.2 (NC, high - low cycling-fraction tertile; ALL as reference), Holm over two ->
specificity (200 matched pseudo sets, observed labels) -> artefact checks (c1, c2, d).
Run: python -I score_a195.py [--n_perm 10000] [--n_pseudo 200] [--power-only]
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd

from clocksync import paired as P
from clocksync.config import limit_threads
from clocksync.decompose import components, logfrac

HERE = Path(__file__).resolve()
C4 = HERE.parents[2]
A187 = C4 / "A187" / "output"
OUT = HERE.parents[1] / "output"
sys.path.insert(0, str(C4 / "A187" / "scripts"))
from score_a187 import STUDIES, cycling_fraction  # noqa: E402

MIN_CELLS, MIN_ARM, LIMIT, SEED = 50, 3, 0.40, 20261102
WF = P.gtex()
CLOCK = list(WF.genes)
META = ["study", "patient", "bin", "n_cells", "depth"]


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def load():
    """One row per study x patient x bin (NC, CYC, ALL): n_cells, depth, all genes (clock genes from the clock tables)."""
    parts = []
    for s in STUDIES:
        g = pd.read_csv(A187 / f"genome_{s}.csv.gz", index_col=[0, 1])
        g.index = pd.MultiIndex.from_arrays([g.index.get_level_values(0).astype(str), g.index.get_level_values(1)])
        pb = pd.read_csv(A187 / f"pseudobulk_clock_{s}.csv.gz")
        pb["patient"] = pb.patient.astype(str)
        pb = pb.set_index(["patient", "bin"])
        assert np.allclose(g.loc[pb.index, CLOCK].to_numpy(float), pb[CLOCK].to_numpy(float))
        x = g.astype(float).join(pb[["n_cells", "depth"]])
        b = x.index.get_level_values(1)
        pat = x.index.get_level_values(0)
        nc = x[b == "NC"].droplevel(1)
        cyc = x[b.isin(["S", "G2M"])].groupby(pat[b.isin(["S", "G2M"])]).sum(min_count=1)
        al = x.groupby(pat).sum(min_count=1)
        for lab, t in (("NC", nc), ("CYC", cyc), ("ALL", al)):
            t = t.copy()
            t.index = t.index.rename("patient")
            parts.append(pd.concat([pd.DataFrame({"study": s, "bin": lab}, index=t.index), t], axis=1).reset_index())
    d = pd.concat(parts, ignore_index=True)
    d["patient"] = d.patient.astype(str)
    return d


def lp(rows, wf=WF):
    """Pooled components over studies (log fractions centred within study)."""
    y, pv, k = logfrac(rows, list(wf.genes), rows.study.to_numpy())
    return components(y, pv, k, wf)


def thin_to(rows, target, rng, genes):
    """Binomially thin each row's `genes` to depth `target` (array) where the row is deeper; depth set accordingly."""
    y = rows.copy()
    p = np.minimum(target / y.depth.to_numpy(float), 1.0)
    y[genes] = rng.binomial(y[genes].to_numpy(float).round().astype(np.int64), p[:, None])
    y["depth"] = np.minimum(y.depth.to_numpy(float), target)
    return y


class Paired:
    """H195.1: patients with >= MIN_CELLS in NC and in CYC; label swap within patient."""
    def __init__(self, d, case="CYC"):
        w = d[d.bin.isin(["NC", case])].pivot_table(index=["study", "patient"], columns="bin", values="n_cells")
        ok = w[(w.NC >= MIN_CELLS) & (w[case] >= MIN_CELLS)].index
        key = pd.MultiIndex.from_frame(d[["study", "patient"]])
        self.a = d[(d.bin == "NC") & key.isin(ok)].set_index(["study", "patient"]).loc[ok].reset_index()
        self.b = d[(d.bin == case) & key.isin(ok)].set_index(["study", "patient"]).loc[ok].reset_index()

    def stat(self, swap=None, wf=WF, a=None, b=None):
        a, b = (self.a if a is None else a), (self.b if b is None else b)
        if swap is not None:
            a, b = a.copy(), b.copy()
            cols = ["depth"] + list(wf.genes)
            av, bv = a[cols].to_numpy(), b[cols].to_numpy()
            a.loc[swap, cols], b.loc[swap, cols] = bv[swap], av[swap]
        return lp(b, wf)["log_plane"] - lp(a, wf)["log_plane"]


class Tertile:
    """H195.2: patients with >= MIN_CELLS NC cells; tertiles of cycling fraction within study; studies with >= MIN_ARM per arm."""
    def __init__(self, d, fr):
        nc = d[(d.bin == "NC") & (d.n_cells >= MIN_CELLS)].merge(fr, on=["study", "patient"])
        nc["tert"] = nc.groupby("study").frac_cycling.transform(lambda f: pd.qcut(f.rank(method="first"), 3, labels=False))
        keep = nc.groupby("study").tert.apply(lambda t: min((t == 0).sum(), (t == 2).sum()) >= MIN_ARM)
        nc = nc[nc.study.isin(keep[keep].index)].reset_index(drop=True)
        self.keys = nc[["study", "patient"]]
        self.lab = nc.tert.to_numpy()
        self.rows = {"NC": nc}
        al = d[d.bin == "ALL"].set_index(["study", "patient"])
        self.rows["ALL"] = al.loc[pd.MultiIndex.from_frame(self.keys)].reset_index()
        self.strata = nc.study.to_numpy()

    def perm(self, rng):
        lab = self.lab.copy()
        for s in np.unique(self.strata):
            i = np.flatnonzero(self.strata == s)
            lab[i] = rng.permutation(lab[i])
        return lab

    def stat(self, which="NC", lab=None, wf=WF, rows=None):
        r = self.rows[which] if rows is None else rows
        l = self.lab if lab is None else lab
        return lp(r[l == 2], wf)["log_plane"] - lp(r[l == 0], wf)["log_plane"]


def pseudo_wfs(d, n, seed):
    """Expression-matched pseudo-clock waveforms from the NC pseudobulks (genes present in every study)."""
    nc = d[d.bin == "NC"]
    g = nc.drop(columns=META).dropna(axis=1, how="any")
    g = g.loc[:, g.sum() > 0]
    pm = np.log(g.div(nc.depth.to_numpy(float), axis=0).mean() + 1e-12)[g.mean() >= 1]
    return P.pseudo_waveforms(pm, n, np.random.default_rng([seed, 7]), WF, "matched")


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p)
    adj = np.maximum.accumulate((m - np.arange(m)) * p[o]).clip(max=1)
    out = np.empty(m); out[o] = adj
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n_perm", type=int, default=10000)
    ap.add_argument("--n_pseudo", type=int, default=200)
    ap.add_argument("--power-only", action="store_true")
    a = ap.parse_args()
    limit_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    full = load()
    fr = cycling_fraction()
    fr["patient"] = fr.patient.astype(str)
    pw = pseudo_wfs(full, a.n_pseudo, SEED)
    d = full[META + CLOCK].copy()  # narrow frames: permutations copy rows
    dp = full[META + CLOCK + sorted({g for w in pw for g in w.genes} - set(CLOCK))].copy()
    del full
    H1, H2 = Paired(d), Tertile(d, fr)
    H1p, H2p = Paired(dp), Tertile(dp, fr)
    prep = pd.DataFrame([dict(test="H195.1", studies=H1.a.study.nunique(), patients=len(H1.a), cells_NC=int(H1.a.n_cells.sum()),
                              cells_CYC=int(H1.b.n_cells.sum())),
                         dict(test="H195.2", studies=len(np.unique(H2.strata)), patients=len(H2.lab), n_low=int((H2.lab == 0).sum()),
                              n_high=int((H2.lab == 2).sum()), studies_in=";".join(np.unique(H2.strata)))])
    prep.to_csv(OUT / "prep.csv", index=False)
    say("\n" + prep.to_string(index=False))

    # ---- signal check: NC coherence vs matched pseudo sets ----
    nc_all = dp[(dp.bin == "NC") & (dp.n_cells >= MIN_CELLS)]
    R = lp(nc_all)["coherence"]
    Rp = np.array([lp(nc_all, w)["coherence"] for w in pw])
    sig = dict(patients=len(nc_all), coherence_R=R, pseudo_q95=float(np.quantile(Rp, .95)), pctile=float((Rp < R).mean() * 100),
               passes=bool(R > np.quantile(Rp, .95)))
    pd.DataFrame([sig]).to_csv(OUT / "signal_check.csv", index=False)
    say("signal check", sig)

    # ---- power (null distributions) ----
    rng = np.random.default_rng(SEED)
    n1 = np.array([H1.stat(rng.random(len(H1.a)) < .5) for _ in range(a.n_perm)])
    labs = [H2.perm(rng) for _ in range(a.n_perm)]
    n2 = np.array([H2.stat("NC", l) for l in labs])
    nall = np.array([H2.stat("ALL", l) for l in labs])
    power = pd.DataFrame([dict(test=t, null_sd=float(n.std(ddof=1)), mde=2.49 * float(n.std(ddof=1)), limit=LIMIT,
                               estimable=2.49 * float(n.std(ddof=1)) <= LIMIT, n_perm=a.n_perm)
                          for t, n in (("H195.1", n1), ("H195.2", n2), ("dALL (reference)", nall))])
    power.to_csv(OUT / "power.csv", index=False)
    say("\n" + power.round(4).to_string(index=False))
    if a.power_only:
        return

    # ---- tests ----
    p2 = lambda n, o: float((1 + np.sum(np.abs(n - n.mean()) >= abs(o - n.mean()))) / (1 + len(n)))
    obs = {"H195.1": H1.stat(), "H195.2": H2.stat("NC"), "dALL (reference)": H2.stat("ALL")}
    nulls = {"H195.1": n1, "H195.2": n2, "dALL (reference)": nall}
    tests = pd.DataFrame([dict(test=t, d_log_plane=o, p_two_sided=p2(nulls[t], o), null_mean=float(nulls[t].mean()))
                          for t, o in obs.items()]).merge(power[["test", "mde", "estimable"]], on="test")
    est = tests.test.isin(["H195.1", "H195.2"]) & tests.estimable
    tests["p_holm"] = np.nan
    if est.any():
        tests.loc[est, "p_holm"] = holm(tests.loc[est, "p_two_sided"])
    # specificity: pseudo sets on the observed labels
    ps = {"H195.1": [H1p.stat(wf=w) for w in pw], "H195.2": [H2p.stat("NC", wf=w) for w in pw],
          "dALL (reference)": [H2p.stat("ALL", wf=w) for w in pw]}
    tests["pseudo_q05"] = [float(np.quantile(ps[t], .05)) for t in tests.test]
    tests["pseudo_median"] = [float(np.median(ps[t])) for t in tests.test]
    tests["pseudo_pctile"] = [float((np.array(ps[t]) < o).mean() * 100) for t, o in zip(tests.test, tests.d_log_plane)]
    tests.to_csv(OUT / "tests.csv", index=False)

    # ---- artefact checks and descriptives ----
    rng = np.random.default_rng(SEED + 1)
    art = []
    a_th = thin_to(H1.a, np.minimum(H1.a.depth.to_numpy(float), H1.b.depth.to_numpy(float)), rng, CLOCK)
    b_th = thin_to(H1.b, a_th.depth.to_numpy(float), rng, CLOCK)
    art.append(dict(check="c1 H195.1, NC thinned to cycling depth", d_log_plane=H1.stat(a=a_th, b=b_th)))
    for case in ("S", "G2M"):
        Hc = Paired(d, case)
        art.append(dict(check=f"d H195.1 {case} - NC", d_log_plane=Hc.stat(), patients=len(Hc.a)))
    nc = H2.rows["NC"]
    tgt = nc.groupby("study").depth.transform(lambda x: x.quantile(.10)).to_numpy(float)
    art.append(dict(check="c2 H195.2, NC thinned to study 10th-percentile depth", d_log_plane=H2.stat("NC", rows=thin_to(nc, tgt, rng, CLOCK))))
    for s in np.unique(H2.strata):
        m = H2.strata == s
        r, l = nc[m], H2.lab[m]
        cmp = lambda rr: components(*logfrac(rr, CLOCK), WF)["log_plane"]
        art.append(dict(check=f"d per study H195.2 {s}", d_log_plane=cmp(r[l == 2]) - cmp(r[l == 0]),
                        d_all=cmp(H2.rows["ALL"][m][l == 2]) - cmp(H2.rows["ALL"][m][l == 0]), n_low=int((l == 0).sum()), n_high=int((l == 2).sum())))
    for k, (pa, pb) in {"H195.1 NC": (H1.a, None), "H195.1 CYC": (H1.b, None)}.items():
        c = lp(pa)
        art.append(dict(check=f"components {k}", log_plane=c["log_plane"], log_off=c["log_off"], coherence=c["coherence"]))
    for which in ("NC", "ALL"):
        for t, nm in ((0, "low"), (2, "high")):
            c = lp(H2.rows[which][H2.lab == t])
            art.append(dict(check=f"components H195.2 {which} {nm}", log_plane=c["log_plane"], log_off=c["log_off"], coherence=c["coherence"]))
    art = pd.DataFrame(art)
    art.to_csv(OUT / "artefacts.csv", index=False)
    pd.set_option("display.width", 250)
    say("\n" + tests.round(4).to_string(index=False))
    say("\n" + art.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
