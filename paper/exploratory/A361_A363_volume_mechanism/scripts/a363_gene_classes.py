"""A363 part 2: when a person's clock volume is low, do the E-box / PAR-bZip targets that cytokines suppress (Cavadini
2007) lose more amplitude than the positive limb, or do all clock genes scale together?

GTEx fast deaths, A360 rows (25 tissues). Per halving of the tissues (20, seed [363, 1]) and role (instrument half X,
measured half Y): instrument l_p = log host radius of half X, residualised on the first two harmonics of its phase and
centred; for each tissue t of Y, the tissue's complex gain Z_t on the half-X direction (state at the mean radius) gives
each gene's predicted clock deviation P = X_g . (Z_t m u_p); per gene class C, pooled least squares of the gene's
centred log expression y on P and P l: e_C = b2 / b1 (elasticity of the class's amplitude to the person volume).
T = mean over halvings and roles of e_supp - e_comp. Null (uniform scaling) and power by simulation at the GTEx design.

    python a363_gene_classes.py --stage power [--jobs 4]   # simulations only (no GTEx expression)
    python a363_gene_classes.py --stage score               # GTEx, registered test
Outputs (`../output/`): a363_sims.csv, a363_power.csv | a363_tests.csv. Aggregates only.
"""
import argparse
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

import vm_common as V
import clocksync as cs
from clocksync._core.factor import conformal
from clocksync.config import limit_threads
from clocksync.elasticity import _design

CLASSES = {"core": (V.CORE, ["PER3", "DBP", "TEF"], ["ARNTL", "NPAS2", "CLOCK"]),
           "14 genes": (V.WF, ["PER1", "PER2", "PER3", "DBP", "TEF"], ["ARNTL", "NPAS2", "CLOCK", "NR1D1"])}
N_NULL, N_ALT, N_SPLITS = 40, 20, 20
A312 = V.EXP / "A312_elasticity" / "output" / "tissues.csv"


def residuals(pb, d, genes) -> np.ndarray:
    """Per Dataset row: log((count + 0.5) / depth) of `genes`, centred within tissue x intercept group."""
    key = pd.Series(np.arange(len(pb)), index=pd.MultiIndex.from_arrays([pb.donor.astype(str), pb.compartment]))
    rows = key.loc[list(zip(d.resection.astype(str), d.compartment))].to_numpy()
    y = np.log((pb[genes].to_numpy(float)[rows] + 0.5) / pb.depth.to_numpy(float)[rows, None])
    grp = (pd.Series(d.compartment.astype(str)) + "|" + pd.Series(d.group.astype(str))).to_numpy()
    return y - pd.DataFrame(y).groupby(grp).transform("mean").to_numpy()


def class_T(pb, ds, w, supp, comp, n_splits, seed):
    """T and the two class elasticities (means over halvings and roles); per-gene elasticities (descriptive)."""
    d = ds.data
    comps = sorted(set(d.compartment))
    genes = list(w.genes)
    Y = residuals(pb, d, genes)
    X = np.asarray(w.X)
    beta = d.beta.copy()
    grp = (pd.Series(d.compartment.astype(str)) + "|" + pd.Series(d.group.astype(str))).to_numpy()
    beta -= pd.DataFrame(beta).groupby(grp).transform("mean").to_numpy()
    res = np.asarray(d.resection).astype(str)
    comp_of = np.asarray(d.compartment)
    gi = {g: genes.index(g) for g in genes}
    rng = np.random.default_rng(seed)
    T, es, ec, per_gene = [], [], [], {g: [] for g in genes}
    for _ in range(n_splits):
        p = list(rng.permutation(comps))
        halves = (sorted(p[: len(p) // 2]), sorted(p[len(p) // 2:]))
        hs = [ds.host(h) for h in halves]
        for xi, yi in ((0, 1), (1, 0)):
            h = hs[xi]
            th = np.angle(h.x.to_numpy() + 1j * h.y.to_numpy())
            H = np.column_stack([np.ones(len(h))] + [f(kk * th) for kk in (1, 2) for f in (np.cos, np.sin)])
            lr = np.log(h.radius.to_numpy())
            ell = pd.Series(lr - H @ np.linalg.lstsq(H, lr, rcond=None)[0], index=h.index.astype(str))
            ell -= ell.mean()
            m = float(h.radius.mean())
            u = pd.DataFrame(np.column_stack([h.x, h.y]) / h.radius.to_numpy()[:, None], index=h.index.astype(str))
            acc = {g: [] for g in genes}
            for t in halves[yi]:
                i = np.flatnonzero((comp_of == t) & np.isin(res, ell.index))
                if len(i) < 10:
                    continue
                f = m * u.loc[res[i]].to_numpy()
                zc = np.linalg.lstsq(_design(f), beta[i].reshape(-1), rcond=None)[0]
                pred = f @ conformal(zc[0] + 1j * zc[1]).T
                P = pred @ X.T                                   # rows x genes: predicted clock deviation
                L = ell.loc[res[i]].to_numpy()
                for g in genes:
                    acc[g].append((Y[i, gi[g]], P[:, gi[g]], L))

            def fit(gs):
                y = np.concatenate([a[0] for g in gs for a in acc[g]])
                P_ = np.concatenate([a[1] for g in gs for a in acc[g]])
                L_ = np.concatenate([a[2] for g in gs for a in acc[g]])
                b = np.linalg.lstsq(np.column_stack([P_, P_ * L_]), y, rcond=None)[0]
                return b[1] / b[0]
            e_s, e_c = fit(supp), fit(comp)
            es.append(e_s)
            ec.append(e_c)
            T.append(e_s - e_c)
            for g in genes:
                per_gene[g].append(fit([g]))
    return float(np.mean(T)), float(np.mean(es)), float(np.mean(ec)), {g: float(np.median(v)) for g, v in per_gene.items()}


def sim_one(args):
    kind, draw = args
    V.quiet()
    limit_threads(1)
    w, supp, comp = CLASSES["core"]
    t = pd.read_csv(A312, index_col="tissue").drop(index=["Whole Blood", "Testis"])
    table = V.tissue_table(t.index)
    presence = (t.n / 244).clip(upper=1.0)
    extra = {"null": 0.0, "x0.5": 0.5, "x0.25": 0.25}[kind]
    rng = np.random.default_rng([363, {"null": 0, "x0.5": 1, "x0.25": 2}[kind], draw])
    pb = V.simulate(table, presence, w, rng, gene_extra={g: extra for g in supp})
    comps = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb, anchors=comps, waveform=w, name=f"a363{kind}{draw}")
    T, e_s, e_c, _ = class_T(pb, ds, w, supp, comp, N_SPLITS, [363, 1])
    return dict(kind=kind, draw=draw, T=T, e_supp=e_s, e_comp=e_c)


def main():
    V.quiet()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    V.OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "power":
        jobs = [("null", d) for d in range(N_NULL)] + [(k, d) for k in ("x0.5", "x0.25") for d in range(N_ALT)]
        rows = []
        with Pool(a.jobs) as pool:
            for r in pool.imap_unordered(sim_one, jobs):
                rows.append(r)
                print(f"{len(rows)}/{len(jobs)} simulated cohorts ({time.time() - t0:.0f}s)", flush=True)
        s = pd.DataFrame(rows).sort_values(["kind", "draw"])
        s.to_csv(V.OUT / "a363_sims.csv", index=False)
        nl = s[s.kind == "null"].T
        res = dict(T_null_median=float(nl.median()), T_null_q05=float(nl.quantile(0.05)), T_null_q95=float(nl.quantile(0.95)))
        for k in ("x0.5", "x0.25"):
            v = s[s.kind == k].T
            res[f"T_{k}_median"] = float(v.median())
            res[f"power_{k}"] = float(np.mean(v > res["T_null_q95"]))
        res["gate"] = bool(res["power_x0.5"] >= 0.8)
        pd.Series(res).to_csv(V.OUT / "a363_power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(V.OUT / "a363_power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("A363 part 2 not powered: not decidable at this design; not scored")
    pb, _, _, _ = V.a359.load()
    comps = sorted(pb.compartment.unique())
    rows = []
    for gs in ("core", "14 genes"):
        w, supp, comp = CLASSES[gs]
        ds = cs.Dataset(pb[V.META + list(w.genes)], anchors=comps, waveform=w, name=f"a363_{gs}")
        T, e_s, e_c, pg = class_T(pb, ds, w, supp, comp, N_SPLITS, [363, 1])
        r = dict(test="H363.2" if gs == "core" else "H363.2 (14 genes, reported)", T=T, e_supp=e_s, e_comp=e_c)
        if gs == "core":
            r.update(T_null_q05=float(pw["T_null_q05"]), T_null_q95=float(pw["T_null_q95"]))
            r["reading"] = ("E-box / PAR-bZip targets follow the volume more (cytokine-like)" if T > r["T_null_q95"] else
                            "positive limb follows the volume more" if T < r["T_null_q05"] else
                            "uniform scaling of the clock genes")
        r.update({f"gene_{g}": v for g, v in pg.items()})
        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(V.OUT / "a363_tests.csv", index=False)
    print(out.round(3).T.to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
