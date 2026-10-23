"""A32 and A35: blood composition rhythms (positive control) and does blood's own clock read clock time in living blood?

A32 (pre-registered): marker-set scores (neutrophil, T, B, NK, monocyte, eosinophil, lymphocyte, L-N; fixed in the legacy
`blood_comp.SETS`) in GTEx whole blood. Timed: score ~ covariates + cos / sin(time of death). Untimed: score ~ covariates + the
host state of the 14 other tissues' clock genes (`Dataset.fit(...).host`), peak in the clock frame. Primary: L-N is timed-rhythmic
(p < 0.05) and its untimed peak lies within 3 h of the timed peak. External truth: GSE39445 (living blood every 3 h, 26 subjects,
438 samples), score ~ block (subject x protocol) + cos / sin(clock time) [+ hours awake, or the melatonin circadian phase],
peak CIs by subject bootstrap.
A35 (pre-registered): phase of blood's 14-gene clock (projection on the global GTEx waveform) vs clock time in GSE39445; arms
(a) single sample, (b) after regressing four marker scores (neutrophil, lymphocyte, monocyte, NK), (c) block-centred, (c+b), (c+b_w).
Statistics: circular correlation rho with clock time, mean resultant R, offset, median |error|, per-block rho (>= 6 samples),
subject bootstrap (B = 1000, seed 35). Primary: rho > 0.3 with bootstrap p < 0.05 in (b), (c) or (c+b). Secondary: (b) vs (a).
Per gene: cosinor with block intercepts, BH, raw and composition-adjusted.

    python compute_blood_clock.py [--subsample N]   (N: GTEx whole-blood donors; fewer bootstrap replicates)
Outputs: paper/results/figure_2/output/gene_rhythms/a32_{gtex,external,external_vs_gtex}.csv, a35_{arms,per_block,genes,phases,verdict}.csv/txt
Legacy: blood_composition/scripts/{blood_comp,external}.py (GTEx part without bootstrap, oracle frame, cross-fit, ventilated donors),
blood_clock/scripts/blood_clock.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as sps  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import projector  # noqa: E402
from common import cli, data, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "gene_rhythms"
W = 2 * np.pi / 24
SETS = {
    "neutrophil": "FCGR3B CXCR2 CXCR1 CSF3R FPR1 FPR2 S100A12 CMTM2 MME AQP9 ALPL KCNJ15 PROK2 MMP25 FFAR2 TNFRSF10C",
    "T": "CD3E CD3D CD3G CD2 CD5 CD6 CD28 LCK CD40LG ITK TRAT1 LEF1 TCF7 BCL11B CD8A CD8B",
    "B": "CD19 MS4A1 CD79A CD79B CD22 BANK1 BLK FCRL1 FCRL2 PAX5 TCL1A VPREB3 CR2 FCRLA TNFRSF13C",
    "NK": "KLRF1 NCR1 SH2D1B KIR2DL1 KIR2DL3 KIR2DL4 KIR3DL1 KIR3DL2 KLRC1 NCAM1 GNLY PRF1 NKG7",
    "monocyte": "CD14 CSF1R VCAN CD300E MS4A6A MS4A7 CD163 CCR2 MAFB CPVL LGALS2 CD86 KCTD12",
    "eosinophil": "CLC IL5RA CCR3 SIGLEC8 PRG2 EPX RNASE2 RNASE3 ALOX15 HRH4 GATA1 ADORA3 ADGRE1 SMPD3",
}
SETS = {k: v.split() for k, v in SETS.items()}
SETS["lymphocyte"] = SETS["T"] + SETS["B"] + SETS["NK"]
ORDER = ["L-N", "lymphocyte", "neutrophil", "T", "B", "NK", "monocyte", "eosinophil"]
COMP = ["neutrophil", "lymphocyte", "monocyte", "NK"]
ARMS = ["a", "b", "c", "c+b", "c+b_w"]
BLOOD = "Whole Blood"
HOST_T = ["Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)",
          "Esophagus - Mucosa", "Nerve - Tibial", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
          "Adipose - Subcutaneous", "Lung", "Pituitary", "Brain - Cerebellum", "Brain - Cortex"]


def wrap(h):
    return (np.asarray(h) + 12) % 24 - 12


def cmean(a):
    return np.angle(np.exp(1j * a).mean())


def circcorr(a, b):
    return float(cs.stats.js_corr(a, b))


def regress(y, cols, Z):
    """OLS of y on [1, Z, cols]: complex loading of the last two columns and their partial-F p. (local: not in clocksync)"""
    Z = Z[:, Z.std(0) > 0] if Z.shape[1] else Z
    D = np.c_[np.ones(len(y)), Z, cols]
    B, *_ = np.linalg.lstsq(D, y, rcond=None)
    rss1 = ((y - D @ B) ** 2).sum(0)
    D0 = D[:, :-2]
    rss0 = ((y - D0 @ np.linalg.lstsq(D0, y, rcond=None)[0]) ** 2).sum(0)
    dfr = len(y) - D.shape[1]
    return B[-2] + 1j * B[-1], sps.f.sf(((rss0 - rss1) / 2) / (rss1 / dfr), 2, dfr)


def block_cosinor(Y, blk, t, extra=None):
    """Cosinor of every column of Y with block intercepts (and `extra` columns): amplitude, peak (h), partial-F p."""
    Bd = pd.get_dummies(pd.Series(blk).astype(str)).to_numpy(float)
    D0 = Bd if extra is None else np.hstack([Bd, extra])
    D1 = np.hstack([D0, np.c_[np.cos(W * t), np.sin(W * t)]])
    B1 = np.linalg.lstsq(D1, Y, rcond=None)[0]
    rss1 = ((Y - D1 @ B1) ** 2).sum(0)
    rss0 = ((Y - D0 @ np.linalg.lstsq(D0, Y, rcond=None)[0]) ** 2).sum(0)
    dfr = len(Y) - np.linalg.matrix_rank(D1)
    c = B1[-2] + 1j * B1[-1]
    return np.abs(c), np.angle(c) / W % 24, sps.f.sf(((rss0 - rss1) / 2) / (rss1 / dfr), 2, dfr)


# ---- GSE39445 -----------------------------------------------------------------------------------------------------

def load_external():
    mod = data._builder("blood_timed_gse39445")
    G, ob = mod.gene_matrix(alias=False)   # A32: platform symbols
    Ga, _ = mod.gene_matrix(alias=True)    # A35: panel symbols
    return G, Ga, ob.reset_index(drop=True)


def marker_scores(G, kept=None):
    thr = G.mean(1).median()
    kept = kept or {k: [g for g in v if g in G.index and G.loc[g].mean() >= thr] for k, v in SETS.items()}
    Z = G.sub(G.mean(1), axis=0).div(G.std(1) + 1e-12, axis=0)
    S = pd.DataFrame({k: Z.loc[v].mean(axis=0) for k, v in kept.items() if len(v) >= 5})
    S["L-N"] = S["lymphocyte"] - S["neutrophil"]
    return S[[c for c in ORDER if c in S]], kept


def external(G, ob, n_boot, rng_seed=3):
    S, kept = marker_scores(G)
    blk, subs = ob.block.to_numpy(), ob.subject.to_numpy()
    out = []

    def fit(idx, tc, extra):
        o = ob.iloc[idx]
        ex = o[["awake_h"]].to_numpy(float) if extra else None
        return block_cosinor(S.iloc[idx].to_numpy(), o.block.to_numpy(), o[tc].to_numpy(float), ex)

    rng = np.random.default_rng(rng_seed)
    for lab, tc, ex in (("clock time", "clock_h", False), ("clock time + hours awake", "clock_h", True),
                        ("circadian phase (melatonin)", "cphase", False)):
        amp, pk, p = fit(np.arange(len(ob)), tc, ex)
        bs = np.array([fit(np.concatenate([np.flatnonzero(subs == s) for s in rng.choice(np.unique(subs), len(np.unique(subs)))]),
                           tc, ex)[1] for _ in range(n_boot)])
        for j, k in enumerate(S.columns):
            dd = wrap(bs[:, j] - pk[j])
            out.append(dict(timescale=lab, score=k, n_genes=len(kept.get(k, [])) if k != "L-N" else np.nan, amp=amp[j], peak_h=pk[j],
                            p=p[j], peak_lo=(pk[j] + np.quantile(dd, .025)) % 24, peak_hi=(pk[j] + np.quantile(dd, .975)) % 24))
    return pd.DataFrame(out), S, kept


# ---- A35 ----------------------------------------------------------------------------------------------------------

def scores_a35(Mk):
    cols = []
    for k in COMP:
        m = Mk[k]
        cols.append(((m - m.mean(0)) / (m.std(0, ddof=1) + 1e-12)).mean(1))
    return np.column_stack(cols)


def block_centre(A, blk):
    A = A.copy()
    for b in np.unique(blk):
        i = blk == b
        A[i] -= A[i].mean(0)
    return A


def ols_resid(Y, D):
    return Y - D @ np.linalg.lstsq(D, Y, rcond=None)[0]


def arm_y(arm, Y, S, blk):
    one = np.ones((len(Y), 1))
    return {"a": lambda: Y - Y.mean(0), "b": lambda: ols_resid(Y, np.hstack([one, S])), "c": lambda: block_centre(Y, blk),
            "c+b": lambda: block_centre(ols_resid(Y, np.hstack([one, S])), blk),
            "c+b_w": lambda: ols_resid(block_centre(Y, blk), block_centre(S, blk))}[arm]()


def score(ph, tr):
    e = ph - tr
    m = np.exp(1j * e).mean()
    err = np.abs(np.angle(np.exp(1j * (e - np.angle(m))))) / W
    return dict(rho=circcorr(ph, tr), R=float(abs(m)), offset_h=float(np.angle(m) / W % 24), med_abs_err_h=float(np.median(err)))


def bh(p):
    return cs.stats.bh(np.asarray(p))


def a35(Ga, G, ob, n_boot, out):
    wf = cs.gtex()
    genes = list(wf.genes)
    Y = (Ga.loc[genes].T * np.log(2)).to_numpy()  # natural log, samples x genes
    S_, kept = marker_scores(G)
    Mk = {k: G.loc[kept[k]].to_numpy().T for k in COMP}
    L = projector(pd.DataFrame(Y, columns=genes), wf, counts=False)  # unweighted projector (arrays)
    blk, t, sub = ob.block.to_numpy(), ob.clock_h.to_numpy(), ob.subject.to_numpy()
    tr, cr = W * t, W * ob.cphase.to_numpy()
    S = scores_a35(Mk)

    def phases(arm, Yy, Ss, bb):
        beta = arm_y(arm, Yy, Ss, bb) @ L.T
        return np.angle(beta[:, 0] + 1j * beta[:, 1])

    rows, pbs, ph_all = [], [], {}
    for arm in ARMS:
        ph = phases(arm, Y, S, blk)
        ph_all[arm] = ph
        pb = pd.Series({b: circcorr(ph[blk == b], tr[blk == b]) for b in np.unique(blk) if (blk == b).sum() >= 6})
        pbs.append(pb.rename(arm))
        rows.append(dict(arm=arm, truth="clock time") | score(ph, tr) | dict(
            block_rho_median=float(pb.median()), block_rho_frac_pos=float((pb > 0).mean()), n_blocks=len(pb)))
        rows.append(dict(arm=arm, truth="circadian phase (melatonin)") | score(ph, cr))
    res = pd.DataFrame(rows)
    rng = np.random.default_rng(35)
    subs = np.unique(sub)
    bs = {arm: [] for arm in ARMS}
    for _ in range(n_boot):
        pick = rng.choice(subs, len(subs), replace=True)
        idx, nb = [], []
        for j, s in enumerate(pick):
            i = np.flatnonzero(sub == s)
            idx.append(i)
            nb.append(np.char.add(blk[i].astype(str), f"#{j}"))
        idx, nb = np.concatenate(idx), np.concatenate(nb)
        Sb = scores_a35({k: v[idx] for k, v in Mk.items()})
        for arm in ARMS:
            bs[arm].append(score(phases(arm, Y[idx], Sb, nb), tr[idx]))
    bsd = {arm: pd.DataFrame(v) for arm, v in bs.items()}
    for arm in ARMS:
        i = (res.arm == arm) & (res.truth == "clock time")
        d = bsd[arm]
        for k in ("rho", "R", "med_abs_err_h"):
            res.loc[i, f"{k}_lo"], res.loc[i, f"{k}_hi"] = d[k].quantile(.025), d[k].quantile(.975)
        res.loc[i, "p_boot_rho_le0"] = float((d.rho <= 0).mean())
    res.to_csv(out / "a35_arms.csv", index=False)
    pd.concat(pbs, axis=1).to_csv(out / "a35_per_block.csv")
    pd.DataFrame({f"phase_h_{a}": (v / W) % 24 for a, v in ph_all.items()}).assign(
        block=blk, clock_h=t).to_csv(out / "a35_phases.csv", index=False)
    sec = {k: (float(res.loc[(res.arm == "b") & (res.truth == "clock time"), k].iloc[0]
                     - res.loc[(res.arm == "a") & (res.truth == "clock time"), k].iloc[0]),
               (bsd["b"][k] - bsd["a"][k]).quantile(.025), (bsd["b"][k] - bsd["a"][k]).quantile(.975)) for k in ("R", "rho")}
    pr = res[(res.truth == "clock time") & res.arm.isin(["b", "c", "c+b"])]
    passed = pr[(pr.rho > 0.3) & (pr.p_boot_rho_le0 < 0.05)].arm.tolist()
    verdict = (f"PRIMARY {'PASS' if passed else 'FAIL'}: arms passing (rho > 0.3, boot p < 0.05): {passed or 'none'}\n"
               + "\n".join(f"  {r.arm}: rho {r.rho:.3f} [{r.rho_lo:.3f}, {r.rho_hi:.3f}], p {r.p_boot_rho_le0:.3f}" for r in pr.itertuples())
               + "\nSECONDARY (b)-(a): " + "; ".join(f"{k} {v[0]:+.3f} [{v[1]:+.3f}, {v[2]:+.3f}]" for k, v in sec.items()))
    (out / "a35_verdict.txt").write_text(verdict + "\n")
    # per gene
    g_rows, wb = [], cs.gtex(BLOOD)
    for lab, ex in (("raw", None), ("composition-adjusted", S)):
        amp, pk, p = block_cosinor(Y, blk, t, ex)
        q = bh(p)
        for j, g in enumerate(genes):
            gp, bp = wf.psi[j] / W % 24, wb.psi[list(wb.genes).index(g)] / W % 24
            g_rows.append(dict(fit=lab, gene=g, amp_ln=amp[j], peak_h=pk[j], p=p[j], q=q[j], gtex_global_psi_h=gp,
                               gtex_global_amp=wf.amp[j], d_vs_global_h=wrap(pk[j] - gp), gtex_blood_psi_h=bp, d_vs_blood_h=wrap(pk[j] - bp)))
    gt = pd.DataFrame(g_rows)
    gt.to_csv(out / "a35_genes.csv", index=False)
    return res, verdict, gt


# ---- A32 on GTEx whole blood --------------------------------------------------------------------------------------

def d2_cov(d):
    Z = pd.get_dummies(pd.DataFrame(dict(sex=d.male.astype(str), age=d.age.astype(str), hardy=d.hardy.astype(str))),
                       drop_first=True).to_numpy(float)
    isch = d.ischaemic_h.astype(float).fillna(d.ischaemic_h.median()).to_numpy()
    isch = (isch - isch.mean()) / (isch.std() + 1e-9)
    return np.c_[Z, isch, isch ** 2]


def a32_gtex(args, run):
    markers = sorted({g for s in SETS.values() for g in s})
    pb, d = gtex.load([BLOOD], genes=markers, depth="library_size")
    pb = pb.set_index("donor")
    d = d.loc[pb.index]
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(len(pb))[:args.subsample]
        pb, d = pb.iloc[keep], d.iloc[keep]
    cpm = pb[[g for g in markers if g in pb.columns]].to_numpy(float) / pb.depth.to_numpy(float)[:, None] * 1e6
    cols = [g for g in markers if g in pb.columns]
    y = np.log(cpm + 0.5)
    z = (y - y.mean(0)) / (y.std(0) + 1e-12)
    med = pd.Series(np.median(cpm, 0), index=cols)
    kept = {k: [g for g in v if g in med.index and med[g] >= 1] for k, v in SETS.items()}
    S = pd.DataFrame({k: z[:, [cols.index(g) for g in v]].mean(1) for k, v in kept.items()}, index=pb.index)
    S["L-N"] = S["lymphocyte"] - S["neutrophil"]
    S = S[ORDER]
    # host of the 14 other tissues (clock genes only), same donors
    hp, _ = gtex.load(HOST_T, depth="library_size")
    hp = hp[hp.donor.isin(pb.index)]
    ds = cs.Dataset(hp, anchors=None, name="gtex_host")
    f = ds.fit(n_boot=0, gate=False).host
    f = f[["x", "y"]].reindex(pb.index)
    ok = f.notna().all(axis=1).to_numpy()
    # frame: the shift from the host frame to clock time, read from one gene's `gene_rhythms` (acrophase vs raw loading)
    g1 = hp[["donor", "compartment", "depth", "PER1"]].assign(stratum="", n_cells=1)
    r = ds.gene_rhythms(g1.rename(columns={"donor": "donor"}), genes=["PER1"], min_mean=0.0)
    rot = float(np.angle(np.exp(1j * (r.acrophase_h * W - np.angle(r.c_re + 1j * r.c_im))).mean()))
    ff = f.to_numpy()[ok]
    ff = ff - ff.mean(0)
    dd, Sm = d[ok], S[ok]
    t = dd.tod.to_numpy(float) % 24
    Z = d2_cov(dd)
    rows = []
    for k in Sm.columns:
        yy = Sm[k].to_numpy()
        ct, pt = regress(yy, np.c_[np.cos(W * t), np.sin(W * t)], Z)
        cu, pu = regress(yy, ff, Z)
        tp, up = np.angle(ct) / W % 24, (np.angle(cu) + rot) / W % 24
        rows.append(dict(group="fast deaths", score=k, n=len(yy), timed_amp=abs(ct), timed_peak_h=tp, timed_p=pt, untimed_peak_h=up,
                         untimed_p=pu, d_h=float(wrap(up - tp)), rot_h=rot / W % 24))
    return pd.DataFrame(rows)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot35 = 20 if run["quick"] else 1000
    G, Ga, ob = load_external()
    ext, S, kept = external(G, ob, 20 if run["quick"] else 500)
    ext.to_csv(out / "a32_external.csv", index=False)
    gt = a32_gtex(args, run)
    gt.to_csv(out / "a32_gtex.csv", index=False)
    e = ext[ext.timescale == "clock time"].set_index("score")
    g = gt.set_index("score")
    cmp = pd.DataFrame(dict(external_peak_h=e.peak_h, external_p=e.p, gtex_timed_peak_h=g.timed_peak_h, gtex_timed_p=g.timed_p,
                            gtex_untimed_peak_h=g.untimed_peak_h))
    cmp["d_gtex_timed_h"] = wrap(cmp.gtex_timed_peak_h - cmp.external_peak_h)
    cmp["d_gtex_untimed_h"] = wrap(cmp.gtex_untimed_peak_h - cmp.external_peak_h)
    cmp.to_csv(out / "a32_external_vs_gtex.csv")
    pd.set_option("display.width", 250)
    print(gt.round(3).to_string(index=False))
    print(ext[ext.timescale == "clock time"].round(3).to_string(index=False))
    res, verdict, genes = a35(Ga, G, ob, n_boot35, out)
    print(res.round(3).to_string(index=False))
    print(verdict)


if __name__ == "__main__":
    main()
