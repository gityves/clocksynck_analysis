"""A334 (registered, revised before data): is a loss of clock loudness carried by cells dephasing (R) or by quieter
cells (lambda)?

Per group (sample x cell class): D, the within-group in-plane dispersion (`clocksync.cells.vonmises.group_stats`), and
rho^2, the squared coherent amplitude of the group's pseudobulk (product of two cell halves' waveform-plane projections
of the mesor-centred log pseudobulk). Pooled D and rho^2 are solved for lambda and kappa (R = I1/I0(kappa)) with their
exact von Mises counterparts. share_R = dlog R / dlog rho (1 pure dephasing, 0 pure damping), rho = lambda R.

    python score_a334.py --stage calibrate   # skin (GSE223109): V334.0-4 and the gate
    python score_a334.py --stage power       # aorta aligned + Bmal1 KO only: clock check, placebo MDEs
    python score_a334.py --stage score       # aorta misaligned vs aligned: H334.1-2 (only if the gates passed)
Outputs (`../output/`): calibration.csv, calib_mix.csv, calib_sim.csv, calib_zt.csv | aorta_power.csv,
aorta_placebo.csv | aorta_tests.csv, aorta_boot.csv. Public data, parsed only.
"""
import argparse
import glob
import gzip
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from joblib import Parallel, delayed
from scipy import optimize, special, stats

import clocksync as cs
from clocksync.cells import vonmises as VM
from clocksync.cells.dispersion import _projector
from clocksync.config import limit_threads
from clocksync.waveform import fit_waveform

HERE = Path(__file__).resolve()
OUT = HERE.parents[1] / "output"
DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
SKIN = DATA / "mouse_skin_timed/GSE223109/raw"
SKIN_CLASSES = Path("/home/maxine/Documents/paychere/scCircadianCancer/legacy/desync/validation/paired_damping/output/"
                    "mouse_skin_cells.csv.gz")
AORTA = DATA / "aorta_jetlag_auerbach2026/adata_qc_filtered_compressed.h5ad"
MIN_CELLS, MIN_UMI, N_SPLITS, N_NULL = 100, 1000, 30, 20
MU = VM.MU_GRID[::6]
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.0f}s] {msg}", flush=True)


class Model:
    """Waveform c_g = A_g exp(-i psi_g), gene-split projectors and the plane projector; exact model D and rho^2."""

    def __init__(self, c, seed):
        self.c = c
        self.LA = VM._split_projectors(c, N_SPLITS, np.random.default_rng(seed))
        self.L = _projector(np.column_stack([c.real, -c.imag]), np.ones(len(c)))

    def D(self, lam, kap):
        return VM.model_stats(lam, kap, self.c, self.LA, mu=MU)[0]

    def rho2(self, lam, kap):
        q = kap * np.exp(-1j * MU)[:, None]
        a = VM.log_i0(q + lam * self.c[None, :]) - VM.log_i0(kap)
        z = (a - a.mean(0)) @ self.L.T
        return float(np.mean((z ** 2).sum(1)))

    def solve(self, D, r2):
        if not (D > 0 and r2 > 0):
            return np.nan, np.nan, np.nan

        def loss(p):
            lam, kap = np.exp(p)
            return np.log(self.D(lam, kap) / D) ** 2 + np.log(max(self.rho2(lam, kap), 1e-12) / r2) ** 2
        best = min((optimize.minimize(loss, np.log([l0, k0]), method="Nelder-Mead",
                                      options=dict(xatol=1e-4, fatol=1e-10, maxiter=400))
                    for l0 in (0.5, 1.0, 2.0) for k0 in (0.5, 2.0, 6.0)), key=lambda r: r.fun)
        lam, kap = np.exp(best.x)
        return float(lam), float(kap), float(special.ive(1, kap) / special.ive(0, kap))


def record(C, N, M, rng, **meta):
    """Per-group statistics: D (all cells), two halves' log pseudobulks, cells."""
    st = VM.group_stats(C, N, M.LA, rng=rng)
    h = rng.permutation(len(N)) % 2 == 0
    y = [np.log(C[s].sum(0) / N[s].sum() + 1e-9) for s in (h, ~h)]
    return dict(D=st["D"], G2=st["G2"], n=len(N), y1=y[0], y2=y[1], **meta)


def mesor(recs):
    """ZT-balanced mean of the halves' log pseudobulks (plain mean if no ZT)."""
    zt = np.array([r.get("zt", 0) for r in recs])
    out = []
    for k in ("y1", "y2"):
        Y = np.array([r[k] for r in recs])
        out.append(np.mean([Y[zt == z].mean(0) for z in np.unique(zt)], 0))
    return out


def estimate(recs, M, mes=None):
    m1, m2 = mesor(recs) if mes is None else mes
    n = np.array([r["n"] for r in recs], float)
    D = float(np.sum(n * [r["D"] for r in recs]) / n.sum())
    r2s = [float((M.L @ (r["y1"] - m1)) @ (M.L @ (r["y2"] - m2))) for r in recs]
    r2 = float(np.sum(n * r2s) / n.sum())
    lam, kap, R = M.solve(D, r2)
    return dict(D=D, rho2=r2, lam=lam, kappa=kap, R=R, rho=lam * R, groups=len(recs), cells=int(n.sum()))


def share(ref, new):
    dl = np.log(new["rho"] / ref["rho"])
    return float(np.log(new["R"] / ref["R"]) / dl) if dl != 0 else np.nan


def d_null(groups, wf, seed):
    """Pooled D (cell-weighted) on N_NULL planes of permuted acrophases. groups: [(C, N)]."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(N_NULL):
        w = wf.permuted(rng)
        LAp = VM._split_projectors(w.amp * np.exp(-1j * w.psi), N_SPLITS, rng)
        st = [VM.group_stats(C, N, LAp, rng=rng) for C, N in groups]
        n = np.array([s["n"] for s in st], float)
        out.append(float(np.sum(n * [s["D"] for s in st]) / n.sum()))
    return np.array(out)


def waveform():
    wf = cs.gtex()
    return wf, wf.amp * np.exp(-1j * wf.psi)


# ---------------------------------------------------------------- skin
def load_skin(genes):
    lab = pd.read_csv(SKIN_CLASSES, index_col=0)
    groups = []
    for m in sorted(glob.glob(str(SKIN / "*matrix.mtx.gz"))):
        p = m.replace("matrix.mtx.gz", "")
        gsm = Path(m).name.split("_")[0]
        feat = [l.split("\t")[1].upper() for l in gzip.open(p + "features.tsv.gz", "rt").read().splitlines()]
        bc = [l.strip() for l in gzip.open(p + "barcodes.tsv.gz", "rt").read().splitlines()]
        X = sp.csr_matrix(sio.mmread(m)).T.tocsr()
        tot = np.asarray(X.sum(1)).ravel()
        Xc = X[:, [feat.index(g) for g in genes]].toarray().astype(float)
        cls = lab.reindex([f"{gsm}_{b}" for b in bc])
        for k in ("fibroblast", "myeloid"):
            sel = (cls.final.to_numpy() == k) & (tot >= MIN_UMI)
            if sel.sum() >= MIN_CELLS:
                s = cls[sel].iloc[0]
                groups.append(dict(C=Xc[sel], N=tot[sel].astype(float), sample=s["sample"], zt=int(s.zt),
                                   rep=int(s.replicate), cls=k, clock_umi=float(Xc[sel].sum() / sel.sum())))
    return groups


def simulate(ng, nc, c, p, lam, kap, mus, depths, rng):
    """Von Mises cells, Poisson counts, depths drawn from the observed ones; no within-group renormalisation."""
    out = []
    for s in range(ng):
        th = rng.vonmises(mus[s], kap, nc) if kap > 0 else rng.uniform(0, 2 * np.pi, nc)
        e = np.exp(lam * np.real(c[None, :] * np.exp(1j * th)[:, None]))
        N = rng.choice(depths, nc)
        out.append((rng.poisson(N[:, None] * p[None, :] * e).astype(float), N))
    return out


def calibrate(jobs):
    wf, c = waveform()
    M = Model(c, [334, 1])
    G = load_skin(list(wf.genes))
    log(f"skin groups: {pd.Series([g['cls'] for g in G]).value_counts().to_dict()}")
    rng = np.random.default_rng([334, 2])
    recs = [record(g["C"], g["N"], M, rng, zt=g["zt"], rep=g["rep"], cls=g["cls"]) for g in G]
    fib = [i for i, g in enumerate(G) if g["cls"] == "fibroblast"]
    rows = []
    for cls in ("fibroblast", "myeloid"):
        idx = [i for i, g in enumerate(G) if g["cls"] == cls]
        e = estimate([recs[i] for i in idx], M)
        nul = d_null([(G[i]["C"], G[i]["N"]) for i in idx], wf, [334, 3])
        eg = VM.solve([dict(D=recs[i]["D"], G2=recs[i]["G2"], n=recs[i]["n"]) for i in idx], c, M.LA)
        rows.append(dict(check="V334.0" if cls == "fibroblast" else "V334.0 (myeloid, secondary)", cls=cls, **e,
                         clock_umi_median=float(np.median([G[i]["clock_umi"] for i in idx])),
                         D_null_max=float(nul.max()), D_null_mean=float(nul.mean()),
                         gamma_route_lam=eg["lam"], gamma_route_R=eg["R"], passes=bool(e["D"] > nul.max())))
        log(f"{cls}: {rows[-1]}")
    base = rows[0]

    # V334.1 empirical dephasing by mixing (fibroblasts)
    key = {(G[i]["zt"], G[i]["rep"]): i for i in fib}
    mix_rows = []
    for delta in (4, 8, 12):
        r = np.random.default_rng([334, 10, delta])
        mrec, rrec, tags = [], [], []
        for i in fib:
            j = key.get(((G[i]["zt"] + delta) % 24, G[i]["rep"]))
            if j is None:
                continue
            n = min(len(G[i]["N"]), len(G[j]["N"]))
            a = r.choice(len(G[i]["N"]), n // 2, replace=False)
            b = r.choice(len(G[j]["N"]), n - n // 2, replace=False)
            mrec.append(record(np.vstack([G[i]["C"][a], G[j]["C"][b]]), np.concatenate([G[i]["N"][a], G[j]["N"][b]]),
                               M, r, zt=G[i]["zt"]))
            k = r.choice(len(G[i]["N"]), n, replace=False)
            rrec.append(record(G[i]["C"][k], G[i]["N"][k], M, r, zt=G[i]["zt"]))
            tags.append((G[i]["zt"], G[j]["zt"]))
        em, er = estimate(mrec, M), estimate(rrec, M)
        sh = share(er, em)
        jk = []
        for z in sorted({t[0] for t in tags}):
            km = [k for k, t in enumerate(tags) if z not in t]
            kr = [k for k, t in enumerate(tags) if t[0] != z]
            jk.append(share(estimate([rrec[k] for k in kr], M), estimate([mrec[k] for k in km], M)))
        jk = np.array(jk)
        nb = len(jk)
        se = float(np.sqrt((nb - 1) / nb * np.nansum((jk - np.nanmean(jk)) ** 2)))
        mix_rows.append(dict(delta_h=delta, groups=len(tags), D_ref=er["D"], D_mix=em["D"], lam_ref=er["lam"],
                             R_ref=er["R"], rho_ref=er["rho"], lam_mix=em["lam"], R_mix=em["R"], rho_mix=em["rho"],
                             rho_ratio=em["rho"] / er["rho"], rho_ratio_expected=float(abs(np.cos(np.pi * delta / 24))),
                             dlog_lam=float(np.log(em["lam"] / er["lam"])), dlog_R=float(np.log(em["R"] / er["R"])),
                             share_R=sh, share_R_jk_se=se, share_R_lo90=sh - 1.645 * se, share_R_hi90=sh + 1.645 * se))
        log(f"mix {delta} h: {mix_rows[-1]}")
    mx = pd.DataFrame(mix_rows)
    mx.to_csv(OUT / "calib_mix.csv", index=False)
    m8 = mx[mx.delta_h == 8].iloc[0]
    rows.append(dict(check="V334.1", cls="fibroblast", share_R=m8.share_R, share_R_lo90=m8.share_R_lo90,
                     passes=bool(m8.share_R >= 0.75 and m8.share_R_lo90 >= 0.5)))

    # V334.2 simulated damping and dephasing at the skin's depth
    p = np.sum([G[i]["C"].sum(0) for i in fib], 0) / np.sum([G[i]["N"].sum() for i in fib])
    depths = np.concatenate([G[i]["N"] for i in fib])
    ncell = int(np.median([len(G[i]["N"]) for i in fib]))
    lam0, kap0, R0 = base["lam"], base["kappa"], base["R"]
    kap1 = optimize.brentq(lambda k: special.ive(1, k) / special.ive(0, k) - 0.5 * R0, 1e-6, 500)

    def sim(rep):
        r = np.random.default_rng([334, 20, rep])
        mus = (np.arange(len(fib)) % 6) * 2 * np.pi / 6 + r.normal(0, 0.1, len(fib))
        zt = list(np.arange(len(fib)) % 6)
        out = {}
        for arm, (lam, kap) in dict(ref=(lam0, kap0), damp=(0.5 * lam0, kap0), dephase=(lam0, kap1)).items():
            g = simulate(len(fib), ncell, c, p, lam, kap, mus, depths, r)
            out[arm] = estimate([record(a, b, M, r, zt=z) for (a, b), z in zip(g, zt)], M)
        return dict(rep=rep, share_damp=share(out["ref"], out["damp"]), share_dephase=share(out["ref"], out["dephase"]),
                    lam_ref=out["ref"]["lam"], R_ref=out["ref"]["R"])
    sims = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(sim)(k) for k in range(100)))
    sims.to_csv(OUT / "calib_sim.csv", index=False)
    v2 = dict(check="V334.2", cls="simulated", lam=lam0, kappa=kap0, R=R0,
              share_damp_median=float(sims.share_damp.median()), share_dephase_median=float(sims.share_dephase.median()),
              wrong_damp=float((sims.share_damp >= 0.5).mean()), wrong_dephase=float((sims.share_dephase <= 0.5).mean()),
              sim_lam_ref_median=float(sims.lam_ref.median()), sim_R_ref_median=float(sims.R_ref.median()))
    v2["passes"] = bool(v2["share_damp_median"] <= 0.25 and v2["share_dephase_median"] >= 0.75
                        and v2["wrong_damp"] < 0.1 and v2["wrong_dephase"] < 0.1)
    rows.append(v2)
    log(f"V334.2: {v2}")

    # V334.3 time invariance (cell bootstrap within groups; mesor of all 18 samples)
    mes = mesor([recs[i] for i in fib])

    def zt_boot(z, b):
        r = np.random.default_rng([334, 30, z, b + 1])
        rr = []
        for i in fib:
            if G[i]["zt"] != z:
                continue
            k = np.arange(len(G[i]["N"])) if b < 0 else r.integers(0, len(G[i]["N"]), len(G[i]["N"]))
            rr.append(record(G[i]["C"][k], G[i]["N"][k], M, r, zt=z))
        e = estimate(rr, M, mes)
        return dict(zt=z, b=b, log_lam=np.log(e["lam"]), R=e["R"], D=e["D"], rho2=e["rho2"])
    zts = sorted({G[i]["zt"] for i in fib})
    zb = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(zt_boot)(z, b) for z in zts for b in range(-1, 200)))
    zt_tab = zb[zb.b < 0].set_index("zt").drop(columns="b").join(
        zb[zb.b >= 0].groupby("zt")[["log_lam", "R"]].std(), rsuffix="_se")
    zt_tab.to_csv(OUT / "calib_zt.csv")
    v3 = dict(check="V334.3", cls="fibroblast")
    for k in ("log_lam", "R"):
        w = 1 / zt_tab[f"{k}_se"] ** 2
        mu = np.sum(w * zt_tab[k]) / w.sum()
        Q = float(np.sum(w * (zt_tab[k] - mu) ** 2))
        v3[f"Q_{k}"], v3[f"p_{k}"] = Q, float(stats.chi2.sf(Q, len(zt_tab) - 1))
    v3["passes"] = bool(v3["p_log_lam"] > 0.01 and v3["p_R"] > 0.01)
    rows.append(v3)
    log(f"V334.3: {v3}")

    # V334.4 depth invariance
    r = np.random.default_rng([334, 40])
    th = [record(r.binomial(G[i]["C"].astype(int), 0.5).astype(float), G[i]["N"] * 0.5, M, r, zt=G[i]["zt"]) for i in fib]
    et = estimate(th, M)
    v4 = dict(check="V334.4", cls="fibroblast", lam_thin=et["lam"], R_thin=et["R"],
              dlog_lam=float(np.log(et["lam"] / base["lam"])), dlog_R=float(np.log(et["R"] / base["R"])))
    v4["passes"] = bool(abs(v4["dlog_lam"]) <= 0.1 and abs(v4["dlog_R"]) <= 0.1)
    rows.append(v4)

    # sensitivity: waveform fitted on the skin fibroblast pseudobulks
    y = np.log(np.array([G[i]["C"].sum(0) / G[i]["N"].sum() for i in fib]) + 1e-7)
    wfit = fit_waveform(y, np.array([G[i]["zt"] for i in fib], float), list(wf.genes), name="skin_fib")
    Mf = Model(wfit.amp * np.exp(-1j * wfit.psi), [334, 50])
    rf = np.random.default_rng([334, 51])
    ef = estimate([record(G[i]["C"], G[i]["N"], Mf, rf, zt=G[i]["zt"]) for i in fib], Mf)
    rows.append(dict(check="sensitivity: skin-fitted waveform", cls="fibroblast", **ef))
    gate = all(x["passes"] for x in rows if x["check"] in ("V334.0", "V334.1", "V334.2"))
    rows.append(dict(check="gate (V334.0-2)", passes=gate))
    pd.DataFrame(rows).to_csv(OUT / "calibration.csv", index=False)
    log("gate " + str(gate))


# ---------------------------------------------------------------- aorta
def load_aorta(genes, arms, M, seed):
    import anndata as ad
    a = ad.read_h5ad(AORTA, backed="r")
    up = {v.upper(): v for v in a.var_names}
    o = a.obs
    arm = np.where(o.bmal1_ko.astype(bool), "ko", np.where(o.misaligned.astype(bool), "misaligned", "aligned"))
    keep = np.isin(arm, arms) & o.cell_type.isin(["SMC", "Fibroblast"]).to_numpy() & (o.total_counts.to_numpy() >= MIN_UMI)
    sub = a[np.flatnonzero(keep), [up[g] for g in genes]].to_memory()
    X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    ob = sub.obs.assign(arm=arm[keep])
    rng = np.random.default_rng(seed)
    recs, raw = [], []
    for (lib, ct), d in ob.groupby(["sample_id", "cell_type"], observed=True):
        if len(d) < MIN_CELLS:
            continue
        k = ob.index.get_indexer(d.index)
        C, N = X[k].astype(float), d.total_counts.to_numpy(float)
        recs.append(record(C, N, M, rng, library=lib, cls=ct, arm=d.arm.iloc[0], zt=int(d.zt.iloc[0]),
                           sex=d.sex.iloc[0], clock_umi=float(C.sum() / len(N))))
        raw.append((C, N))
    return recs, raw


def aorta_power(jobs):
    wf, c = waveform()
    M = Model(c, [334, 60])
    recs, raw = load_aorta(list(wf.genes), ["aligned", "ko"], M, [334, 61])
    log(f"aorta groups: {pd.Series([(r['arm'], r['cls']) for r in recs]).value_counts().to_dict()}")
    rows, plc = [], []
    for cls in ("SMC", "Fibroblast"):
        al = [i for i, r in enumerate(recs) if r["cls"] == cls and r["arm"] == "aligned"]
        ko = [i for i, r in enumerate(recs) if r["cls"] == cls and r["arm"] == "ko"]
        ra = [recs[i] for i in al]
        ea = estimate(ra, M)
        mes = mesor(ra)
        nul = d_null([raw[i] for i in al], wf, [334, 62])
        nko = np.array([recs[i]["n"] for i in ko], float)
        D_ko = float(np.sum(nko * [recs[i]["D"] for i in ko]) / nko.sum()) if ko else np.nan
        clock_ok = bool(ea["D"] > nul.max() and D_ko <= 0.5 * ea["D"])
        by = {}
        for i in al:
            by.setdefault(recs[i]["zt"], []).append(i)

        def split(s):
            r = np.random.default_rng([334, 63, s])
            a8 = [int(i) for z, v in sorted(by.items()) for i in r.choice(v, min(2, len(v)), replace=False)]
            b6 = [i for i in al if i not in a8]
            e8, e6 = estimate([recs[i] for i in a8], M, mes), estimate([recs[i] for i in b6], M, mes)
            return dict(cls=cls, split=s, dlog_rho=np.log(e8["rho"] / e6["rho"]), dlog_lam=np.log(e8["lam"] / e6["lam"]),
                        dlog_R=np.log(e8["R"] / e6["R"]))
        pl = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(split)(s) for s in range(200)))
        plc.append(pl)
        mde = {k: float(2.49 * pl[k].std()) for k in ("dlog_rho", "dlog_lam", "dlog_R")}
        powered = bool(all(v <= np.log(1 / 0.8) for v in mde.values()))
        rows.append(dict(cls=cls, libraries_aligned=len(al), libraries_ko=len(ko), cells_aligned=ea["cells"],
                         clock_umi_median=float(np.median([recs[i]["clock_umi"] for i in al])),
                         D_aligned=ea["D"], D_null_max=float(nul.max()), D_ko=D_ko, ko_over_aligned=D_ko / ea["D"],
                         rho2_aligned=ea["rho2"], lam_aligned=ea["lam"], R_aligned=ea["R"], rho_aligned=ea["rho"],
                         clock_check=clock_ok, **{f"MDE_{k}": v for k, v in mde.items()}, powered=powered,
                         gate=bool(clock_ok and powered)))
        log(str(rows[-1]))
    pd.DataFrame(rows).to_csv(OUT / "aorta_power.csv", index=False)
    pd.concat(plc).to_csv(OUT / "aorta_placebo.csv", index=False)


def aorta_score(jobs):
    cal = pd.read_csv(OUT / "calibration.csv").set_index("check")
    pw = pd.read_csv(OUT / "aorta_power.csv").set_index("cls")
    if str(cal.loc["gate (V334.0-2)", "passes"]) != "True":
        raise SystemExit("calibration gate failed: decomposition not run")
    wf, c = waveform()
    M = Model(c, [334, 60])
    recs, _ = load_aorta(list(wf.genes), ["aligned", "misaligned"], M, [334, 61])
    rows, boots = [], []
    for cls in ("SMC", "Fibroblast"):
        if str(pw.loc[cls, "clock_check"]) != "True":
            rows.append(dict(cls=cls, verdict="D not clock-specific in the aorta"))
            continue
        idx = [i for i, r in enumerate(recs) if r["cls"] == cls]
        al = [i for i in idx if recs[i]["arm"] == "aligned"]
        mi = [i for i in idx if recs[i]["arm"] == "misaligned"]
        zt = {i: recs[i]["zt"] for i in idx}
        zts = sorted(set(zt.values()))

        def contrast(a, b):
            ea, eb = estimate([recs[i] for i in a], M), estimate([recs[i] for i in b], M)
            return dict(dlog_rho=np.log(eb["rho"] / ea["rho"]), dlog_lam=np.log(eb["lam"] / ea["lam"]),
                        dlog_R=np.log(eb["R"] / ea["R"]), dD=eb["D"] - ea["D"], ea=ea, eb=eb)
        obs = contrast(al, mi)

        def perm(s):
            r = np.random.default_rng([334, 70, s])
            a, b = [], []
            for z in zts:
                v = [i for i in idx if zt[i] == z]
                pick = set(int(x) for x in r.choice(v, sum(zt[i] == z for i in mi), replace=False))
                b += [i for i in v if i in pick]
                a += [i for i in v if i not in pick]
            return contrast(a, b)["dlog_rho"]

        def boot(s):
            r = np.random.default_rng([334, 71, s])
            a = [int(i) for z in zts for i in r.choice([j for j in al if zt[j] == z], sum(zt[j] == z for j in al))]
            b = [int(i) for z in zts for i in r.choice([j for j in mi if zt[j] == z], sum(zt[j] == z for j in mi))]
            d = contrast(a, b)
            return dict(cls=cls, boot=s, **{k: d[k] for k in ("dlog_rho", "dlog_lam", "dlog_R", "dD")})
        null = np.array(Parallel(n_jobs=jobs)(delayed(perm)(s) for s in range(2000)))
        bt = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(boot)(s) for s in range(2000)))
        bt["share_R"] = bt.dlog_R / bt.dlog_rho
        boots.append(bt)
        p = float((1 + np.sum(null <= obs["dlog_rho"])) / 2001)
        sh = float(obs["dlog_R"] / obs["dlog_rho"])
        lo, hi = np.nanquantile(bt.share_R, [0.05, 0.95])
        t = dict(cls=cls, primary=cls == "SMC", powered=str(pw.loc[cls, "powered"]) == "True",
                 libraries_aligned=len(al), libraries_misaligned=len(mi),
                 **{f"{k}_aligned": obs["ea"][k] for k in ("D", "rho2", "lam", "R", "rho")},
                 **{f"{k}_misaligned": obs["eb"][k] for k in ("D", "rho2", "lam", "R", "rho")},
                 dlog_rho=obs["dlog_rho"], p_one_sided=p, dlog_lam=obs["dlog_lam"], dlog_R=obs["dlog_R"], dD=obs["dD"],
                 **{f"{k}_{q}": float(bt[k].quantile(v)) for k in ("dlog_rho", "dlog_lam", "dlog_R", "dD")
                    for q, v in (("lo90", 0.05), ("hi90", 0.95))},
                 share_R=sh, share_R_lo90=float(lo), share_R_hi90=float(hi))
        h1 = obs["dlog_rho"] < 0 and p < 0.05
        t["H334.1"] = "PASS" if h1 else "FAIL"
        t["H334.2"] = ("not scored (H334.1 failed)" if not h1 else "dephasing" if lo > 0.5 else
                       "damping" if hi < 0.5 else "not separated")
        if not t["powered"]:
            t["note"] = "not powered: descriptive"
        rows.append(t)
        log(str(t))
    pd.DataFrame(rows).to_csv(OUT / "aorta_tests.csv", index=False)
    if boots:
        pd.concat(boots).to_csv(OUT / "aorta_boot.csv", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["calibrate", "power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    {"calibrate": calibrate, "power": aorta_power, "score": aorta_score}[a.stage](a.jobs)
    log("done")


if __name__ == "__main__":
    main()
