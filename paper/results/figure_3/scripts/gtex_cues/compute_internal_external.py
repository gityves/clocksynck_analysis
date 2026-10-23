"""A91 and A27 D5: are the rhythmic genes of a tissue driven by internal (clock) time or by external (wall) time?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A91; A27 part 2 D5 is post hoc, exploratory).
A91 (GTEx fast deaths, 26 tissues, genes with median CPM >= 5): per tissue and gene, log(CPM + 0.5) ~ D2 covariates (sex,
age bracket, Hardy, ischaemic time and its square) + f (host clock state of the other `v2_gtex` tissues, 2 columns) +
cos/sin(w t_death) (2 columns). Partial F: internal given external (p_int), external given internal (p_ext). Calls:
clock-following p_int < 0.01 and p_ext > 0.05, time-following p_ext < 0.01 and p_int > 0.05, both p_int and p_ext < 0.01.
Noise floor per tissue: the median count of non-clock genes called with f permuted over donors (clock calls), with t
permuted (time calls), larger of the two ("both"), 20 draws each. Positive control: share of the clock-plane genes called
clock-following or both. H91: in tissues with >= 3 rhythmic clock-plane genes and >= 50 rhythmic calls above chance, the
ventilated lock (A81b/A84/A89 table) is higher for clock-following than time-following genes; statistic = median over
tissues of the difference of median lock; class labels permuted within tissue (2,000), one-sided. Descriptive: timed
amplitude on t_death in donors >= 50 y vs < 50 y at equal n (20 draws), ratio per gene.
A27 D5 (ventilated donors, 14 non-pituitary tissues, host from the other `v2_gtex` tissues of the ventilated donors): the same
model; number of non-clock genes with p < 0.001 against f- and t-permuted nulls (50), internal-only / external-only / both
counts, the phase of internal-only genes against the fast-death timed cosinor, and the host-vs-time alignment of the
cohorts. (The CHIRAL-DIP arm is `compute_dip_arm.py`.)
Hosts: `clocksync` `Dataset.fit(...).host`. Aggregates only (GTEx controlled access).

    python compute_internal_external.py [--subsample N] [--jobs 2]
Outputs: a91_tissues.csv, a91_genes.csv.gz (gene calls; used by A100-A102b), a91_lock_by_class.csv, a91_age.csv,
a91_verdict.csv, diag5_summary.csv, diag5_genes.csv.gz, diag5_host_vs_time.csv.
Legacy: legacy/desync/validation/paired_damping/scripts/a91_internal_external.py; glucocorticoid/scripts/diag5_int_ext.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _genes as G  # noqa: E402
import _shared as S  # noqa: E402
from _genes import FOLDER, W  # noqa: E402

SEED, N_NULL, N_SUB = 91, 20, 20


def calls(Y, Zc, f, E):
    D = np.c_[np.ones(len(Y)), Zc, f, E]
    k = D.shape[1]
    _, p_int = S.partial_f(Y, D, [k - 4, k - 3])
    _, p_ext = S.partial_f(Y, D, [k - 2, k - 1])
    cls = np.full(Y.shape[1], "none", dtype=object)
    cls[(p_int < 0.01) & (p_ext > 0.05)] = "clock"
    cls[(p_ext < 0.01) & (p_int > 0.05)] = "time"
    cls[(p_int < 0.01) & (p_ext < 0.01)] = "both"
    return cls, p_int, p_ext


def amp_t(Y, t, Zc):
    B, _ = S.lstsq_fit(Y, np.c_[np.cos(W * t), np.sin(W * t)], Zc)
    return np.abs(B[-2] + 1j * B[-1])


def a91_tissue(tis, fv, fast_sub, seed, cov=None):
    """`cov(Y, ob, f)` -> covariate matrix (default: the D2 covariates; A102 passes the extended set)."""
    X, genes, ob = G.load_rows(tis, (1, 2), fast_sub)
    expr = np.median(X, 0) >= 5
    Y, genes = np.log(X[:, expr] + 0.5), genes[expr]
    rng = np.random.default_rng(seed)
    subj = ob.subject.to_numpy()
    m = np.isin(subj, fv.index)
    Y, ob = Y[m], ob[m]
    f = fv[["x", "y"]].loc[subj[m]].to_numpy()
    t = ob.tod.to_numpy(float) % 24
    E = np.c_[np.cos(W * t), np.sin(W * t)]
    Zc = S.d2_covariates(ob) if cov is None else cov(Y, ob, f)
    f = f - f.mean(0)
    Zc = Zc[:, Zc.std(0) > 0]
    cls, p_int, p_ext = calls(Y, Zc, f, E)
    nonclock = ~np.isin(genes, list(G.CLOCK))
    count = lambda c: {k: int(((c == k) & nonclock).sum()) for k in ("clock", "time", "both")}
    obs = count(cls)
    nf = pd.DataFrame([count(calls(Y, Zc, f[rng.permutation(len(f))], E)[0]) for _ in range(N_NULL)]).median()
    nt = pd.DataFrame([count(calls(Y, Zc, f, E[rng.permutation(len(E))])[0]) for _ in range(N_NULL)]).median()
    nm = {"clock": nf["clock"], "time": nt["time"], "both": max(nf["both"], nt["both"])}
    fr = np.isin(genes, list(G.CLOCK))
    old = ob.age_bracket.astype(str).str.extract(r"(\d+)")[0].astype(float).to_numpy() >= 50
    io, iy = np.flatnonzero(old), np.flatnonzero(~old)
    n = min(len(io), len(iy))
    ratio = np.full(len(genes), np.nan)
    if n >= 20:
        aa = []
        for _ in range(N_SUB):
            so, sy = rng.choice(io, n, replace=False), rng.choice(iy, n, replace=False)
            zo, zy = Zc[so][:, Zc[so].std(0) > 0], Zc[sy][:, Zc[sy].std(0) > 0]
            aa.append(amp_t(Y[so], t[so], zo) / amp_t(Y[sy], t[sy], zy))
        ratio = np.median(aa, axis=0)
    summ = dict(tissue=tis, donors=len(Y), genes=int(nonclock.sum()), **{f"n_{k}": v for k, v in obs.items()},
                **{f"null_{k}": float(nm[k]) for k in obs},
                clock_genes_called=float(np.isin(cls[fr], ["clock", "both"]).mean()), age_n_per_group=int(n))
    return summ, pd.DataFrame(dict(tissue=tis, gene=genes, cls=cls, p_int=p_int, p_ext=p_ext, age_ratio=ratio))


def h91(per: pd.DataFrame, summ: pd.DataFrame, g26: pd.DataFrame, n_perm: int, rng) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    frame = g26[g26.gene.isin(G.CLOCK) & (g26.p_fast < 0.01)].groupby("tissue").size()
    excess = g26.groupby("tissue").apply(lambda d: (d.p_fast < 0.01).sum() - 0.01 * len(d), include_groups=False)
    both = [t for t in summ.tissue if frame.get(t, 0) >= 3 and excess.get(t, 0) >= 50]
    x = per[per.tissue.isin(both) & per.cls.isin(["clock", "time"])].merge(
        g26[["tissue", "gene", "lock"]], on=["tissue", "gene"]).dropna(subset=["lock"])
    stat = lambda d: np.nanmedian([g[g.cls == "clock"].lock.median() - g[g.cls == "time"].lock.median()
                                   for _, g in d.groupby("tissue") if (g.cls == "clock").any() and (g.cls == "time").any()])
    obs = stat(x)
    null = np.array([stat(x.assign(cls=x.groupby("tissue").cls.transform(lambda c: rng.permutation(c.to_numpy()))))
                     for _ in range(n_perm)])
    v = dict(tissues=len(both), median_diff=float(obs), p_upper=cs.stats.perm_p(null, obs, "greater", tie=0.0),
             n_perm=n_perm)
    v["H91"] = "PASS" if (obs > 0 and v["p_upper"] < 0.05) else "FAIL"
    age = (per[per.cls.isin(["clock", "time"])].groupby(["tissue", "cls"]).age_ratio.median().unstack()
           if "age_ratio" in per else None)
    return v, x.groupby(["tissue", "cls"]).lock.agg(["size", "median"]).unstack(), age


# ---- A27 D5: ventilated donors --------------------------------------------------------------------------------

def d5_tissue(tis, fv, seed, n_perm=50):
    Xf, gf, obf = G.load_rows(tis, (1, 2))
    Xv, gv, obv = G.load_rows(tis, (0,))
    expr = (np.median(Xf, 0) >= 5) & (np.median(Xv, 0) >= 5)
    Yf, Yv, genes = np.log(Xf[:, expr] + 0.5), np.log(Xv[:, expr] + 0.5), gf[expr]
    rng = np.random.default_rng(seed)
    subj = obv.subject.to_numpy()
    m = np.isin(subj, fv.index)
    Yv, obv, subj = Yv[m], obv[m], subj[m]
    f = fv[["x", "y"]].loc[subj].to_numpy()
    f = f - f.mean(0)
    t = obv.tod.to_numpy(float) % 24
    E = np.c_[np.cos(W * t), np.sin(W * t)]
    Z = S.d2_covariates(obv)
    Z = Z[:, Z.std(0) > 0]

    def run(fm, Em):
        D = np.c_[np.ones(len(Yv)), Z, fm, Em]
        k = D.shape[1]
        B, p_int = S.partial_f(Yv, D, [k - 4, k - 3])
        _, p_ext = S.partial_f(Yv, D, [k - 2, k - 1])
        return B, p_int, p_ext

    B, p_int, p_ext = run(f, E)
    c_int, c_ext = B[-4] + 1j * B[-3], B[-2] + 1j * B[-1]
    c_fast, p_fast = S.timed_cosinor(Yf, obf.tod.to_numpy(float) % 24, S.d2_covariates(obf))
    fr = np.isin(genes, list(G.CLOCK))
    nonf = ~fr
    n_int0, n_ext0 = [], []
    for _ in range(n_perm):
        _, pi0, _ = run(f[rng.permutation(len(f))], E)
        _, _, pe0 = run(f, E[rng.permutation(len(E))])
        n_int0.append(int(((pi0 < 0.001) & nonf).sum()))
        n_ext0.append(int(((pe0 < 0.001) & nonf).sum()))
    rot = np.angle(np.exp(1j * (np.angle(c_int[fr]) - np.angle(c_fast[fr]))).mean())
    io = (p_int < 0.001) & (p_ext > 0.1) & nonf
    sel = io & (p_fast < 0.05)
    Ri = np.exp(1j * (np.angle(c_int[sel]) - rot - np.angle(c_fast[sel]))).mean() if sel.sum() else np.nan
    phi = np.angle(f[:, 0] + 1j * f[:, 1])
    dec = np.exp(1j * (phi - W * t)).mean()
    n_i, n_e = int(((p_int < 0.001) & nonf).sum()), int(((p_ext < 0.001) & nonf).sum())
    summ = dict(tissue=tis, n=len(Yv), genes=int(nonf.sum()), f_vs_t_R=abs(dec),
                int_sig=n_i, int_null95=float(np.quantile(n_int0, 0.95)), int_null_med=float(np.median(n_int0)),
                p_int_count=(1 + sum(x >= n_i for x in n_int0)) / (n_perm + 1),
                ext_sig=n_e, ext_null95=float(np.quantile(n_ext0, 0.95)), ext_null_med=float(np.median(n_ext0)),
                p_ext_count=(1 + sum(x >= n_e for x in n_ext0)) / (n_perm + 1),
                int_only=int(io.sum()), ext_only=int(((p_ext < 0.001) & (p_int > 0.1) & nonf).sum()),
                both=int(((p_ext < 0.001) & (p_int < 0.001) & nonf).sum()),
                int_only_fast_rhythmic=int(sel.sum()), int_only_phase_R=abs(Ri) if sel.sum() else np.nan,
                int_only_phase_rot_h=float(np.angle(Ri) / W % 24) if sel.sum() else np.nan,
                frame_int_sig=int(((p_int < 0.001) & fr).sum()), frame_ext_sig=int(((p_ext < 0.001) & fr).sum()),
                frame_n=int(fr.sum()))
    per = pd.DataFrame(dict(tissue=tis, gene=genes, frame=fr, p_int=p_int, p_ext=p_ext, p_fast=p_fast,
                            amp_int=np.abs(c_int), amp_ext=np.abs(c_ext), amp_fast=np.abs(c_fast),
                            int_h=(np.angle(c_int) - rot) / W % 24, ext_h=np.angle(c_ext) / W % 24,
                            fast_h=np.angle(c_fast) / W % 24))
    return summ, per


def host_vs_time(o: pd.DataFrame, sub) -> pd.DataFrame:
    """Alignment of the 15-tissue host phase with clock time in fast deaths and ventilated donors (R, circular correlation)."""
    rows = []
    for lab, hardy in (("fast deaths", (1, 2)), ("ventilated", (0,))):
        pb, _ = S.clock_table(hardy=hardy)
        h = S.host_states(pb, list(S.TISSUES), donors=sub)
        d = o[o.DTHHRDY.isin(hardy) & o.tod.notna()].drop_duplicates("subject").set_index("subject")
        d = d.loc[d.index.intersection(h.index)]
        phi = np.angle(h.x.loc[d.index] + 1j * h.y.loc[d.index])
        tt = d.tod.to_numpy(float) % 24
        R = np.exp(1j * (phi - W * tt)).mean()
        rows.append(dict(cohort=lab, donors=len(d), R=float(abs(R)),
                         circ_corr=float(np.corrcoef(np.cos(phi - np.angle(R)), np.cos(W * tt))[0, 1])))
    return pd.DataFrame(rows)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    o = G.obs_table()
    sub = G.pick_donors(o, args.subsample, args.seed)
    H = G.halves(o, sub)
    n_perm = 100 if run["quick"] else 2000
    g26 = G.selection_tables(out, args.subsample, args.seed, args.jobs)

    # A91: host from the fast deaths (all donors) of the other 15 tissues
    fast = H["fast"]["sel"] | H["fast"]["eval"]
    hosts = G.hosts_for(G.ALL26, (1, 2), donors=sub and fast)
    res = Parallel(n_jobs=args.jobs)(delayed(a91_tissue)(t, hosts[t], fast if sub else None, SEED + i)
                                     for i, t in enumerate(G.ALL26))
    summ = pd.DataFrame([r[0] for r in res])
    per = pd.concat([r[1] for r in res], ignore_index=True)
    summ.to_csv(out / "a91_tissues.csv", index=False)
    per.to_csv(out / "a91_genes.csv.gz", index=False)
    v, lock_by, age = h91(per, summ, g26, n_perm, np.random.default_rng(SEED))
    lock_by.to_csv(out / "a91_lock_by_class.csv")
    age.to_csv(out / "a91_age.csv")
    pd.Series(v).to_csv(out / "a91_verdict.csv", header=False)
    pd.set_option("display.width", 250)
    print(summ.round(2).to_string(index=False), "\n", pd.Series(v).to_string(), flush=True)

    # A27 D5: ventilated donors, all ventilated donors form the host
    targets = [t for t in S.TISSUES if t != "Pituitary"]
    vent = None if sub is None else sub
    hv = G.hosts_for(targets, (0,), donors=vent)
    res = Parallel(n_jobs=args.jobs)(delayed(d5_tissue)(t, hv[t], 50 + i, 10 if run["quick"] else 50)
                                     for i, t in enumerate(targets))
    pd.DataFrame([r[0] for r in res]).to_csv(out / "diag5_summary.csv", index=False)
    pd.concat([r[1] for r in res]).to_csv(out / "diag5_genes.csv.gz", index=False)
    host_vs_time(o, sub).to_csv(out / "diag5_host_vs_time.csv", index=False)
    print(pd.DataFrame([r[0] for r in res]).round(2).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
