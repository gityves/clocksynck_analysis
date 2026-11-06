"""A225 POST HOC descriptives (not registered): why does A224's dosage adjustment cut the real LumA - others D (0.276 -> 0.101) while adjusted ~ raw in the A225 simulation?
Reuses compute_a225.py up to the stage branch (scratch output dir, no registered output touched); reads output/sim_calibration.csv (lam). Writes output/posthoc_*.csv.
  posthoc_arm_coherence_real.csv / _sim.csv : raw and adjusted coherence per cell and arm (sim: scenario N, mean over replicates)
  posthoc_offset_plane.csv                  : CN offset (beta_g CN_g) vs clock plane, per arm, real vs sim
  posthoc_cn_vs_clock_residual.csv          : per gene and arm, corr(CN_g, residual of gene g after the clock-plane fit from the other genes), real vs sim
  posthoc_variants_real.csv                 : adjustment variants (cohort vs pooled beta, with / without phase covariates), D and arm coherences
    python -I posthoc_a225.py [--nrep 100]"""
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUTD = HERE.parent / "output"
nrep = int(sys.argv[sys.argv.index("--nrep") + 1]) if "--nrep" in sys.argv else 100
src = (HERE / "compute_a225.py").read_text(); src = src[: src.index('if args.stage == "sim":')]
sys.argv = ["compute_a225.py", "--stage", "power", "--out", tempfile.mkdtemp(prefix="a225_ph_"), "--jobs", "1"]
P = str(HERE / "compute_a225.py"); ns_ = {"__file__": P, "__name__": "a225"}
exec(compile(src, P, "exec"), ns_); globals().update({k: v for k, v in ns_.items() if not k.startswith("__")})
import numpy as np, pandas as pd  # noqa: E402

ONE = lambda Y: np.ones_like(Y)
acoh = lambda Y, idx: A.coh_masked(Y, ONE(Y), idx, WF)
XP = WF.X; PL = XP @ np.linalg.inv(XP.T @ XP) @ XP.T
cal = pd.read_csv(OUTD / "sim_calibration.csv").set_index("cohort")
beta_real = {n: np.nan_to_num(np.array([dos[(dos.cohort == ("METABRIC" if n == "METABRIC" else "TCGA")) & (dos.gene == g)].beta_ln.iloc[0] for g in GENES], float)) for n in COH}


def arm_rows(Yr, Ya, name, tag):
    rows = []
    for cl in COH[name]["cells"]:
        for lab, arm in ((1, "LumA"), (0, "others")):
            ii = cl["idx"][cl["lab"] == lab]; r, a = acoh(Yr, ii), acoh(Ya, ii)
            rows.append(dict(source=tag, cohort=name, cell=cl["key"], arm=arm, n=len(ii), coh_raw=r, coh_adjusted=a, change=a - r))
    return rows


def plane_stats(Y, o, CN, idx):
    Yc, oc = Y[idx] - Y[idx].mean(0), o[idx] - o[idx].mean(0); z, s = Yc @ CW, oc @ CW
    Co, Cy = np.cov(oc.T), np.cov(Yc.T)
    return dict(n=len(idx), offset_sd_total=float(np.sqrt(np.trace(Co))), offset_plane_share=float(np.trace(PL @ Co @ PL) / np.trace(Co)),
                offset_plane_var_over_expr_plane_var=float(np.trace(PL @ Co @ PL) / np.trace(PL @ Cy @ PL)),
                expr_plane_share=float(np.trace(PL @ Cy @ PL) / np.trace(Cy)),
                complex_corr=float(abs(np.sum(np.conj(z) * s)) / np.sqrt(np.sum(abs(z) ** 2) * np.sum(abs(s) ** 2))),
                corr_re=float(np.corrcoef(z.real, s.real)[0, 1]), corr_im=float(np.corrcoef(z.imag, s.imag)[0, 1]))


def resid_corr(Y, CN, idx):
    """per gene: corr(CN_g, residual of y_g after regressing on the 2-d clock state of the other genes), and corr(CN_g, that fitted part); arm-centred."""
    y = Y[idx] - Y[idx].mean(0); out = []
    for g in range(y.shape[1]):
        keep = np.arange(y.shape[1]) != g; Xk = XP[keep]; st = y[:, keep] @ Xk @ np.linalg.inv(Xk.T @ Xk)
        co = np.linalg.lstsq(st, y[:, g], rcond=None)[0]; fit = st @ co; res = y[:, g] - fit; cn = CN[idx, g]; ok = np.isfinite(cn)
        if ok.sum() < 10 or cn[ok].std() == 0:
            out.append((np.nan, np.nan)); continue
        out.append((np.corrcoef(cn[ok], res[ok])[0, 1], np.corrcoef(cn[ok], fit[ok])[0, 1]))
    return out


def per_arm(Y, o, CN, name, tag, off_rows, rc_rows):
    for cl in COH[name]["cells"]:
        for lab, arm in ((1, "LumA"), (0, "others")):
            ii = cl["idx"][cl["lab"] == lab]
            off_rows.append(dict(source=tag, cohort=name, cell=cl["key"], arm=arm, **plane_stats(Y, o, CN, ii)))
            for g, (r, f) in zip(GENES, resid_corr(Y, CN, ii)):
                rc_rows.append(dict(source=tag, cohort=name, cell=cl["key"], arm=arm, gene=g, corr_cn_residual=r, corr_cn_clockfit=f))


# ---- real ---------------------------------------------------------------------------------------------------------------------------------------------
rr, off, rc = [], [], []
for name, c in COH.items():
    rr += arm_rows(c["Y"], c["Yadj"], name, "real")
    per_arm(c["Y"], np.nan_to_num(c["CN"]) * beta_real[name], c["CN"], name, "real", off, rc)
say("real done")
# ---- simulation, scenario N ---------------------------------------------------------------------------------------------------------------------------
sr, sim_off, sim_rc = [], [], []
for k, (name, c) in enumerate(COH.items()):
    par = make_par(name)
    for r in range(nrep):
        Y, _ = A.gen_clock(par, cal.loc[name, "lam"], 1.0, np.random.default_rng([SEED, 500, r, k]))
        Ya, _ = A.adjust(Y, par["cn"], A.leave_out_phases(Y, par["cw"]))
        sr += [dict(x, rep=r) for x in arm_rows(Y, Ya, name, "sim_N")]
        o2, r2 = [], []; per_arm(Y, par["cn0"] * par["beta"], par["cn"], name, "sim_N", o2, r2); sim_off += [dict(x, rep=r) for x in o2]; sim_rc += [dict(x, rep=r) for x in r2]
    say("sim", name, "done")
key = ["source", "cohort", "cell", "arm"]
simm = pd.DataFrame(sr).groupby(key, sort=False).agg(n=("n", "first"), coh_raw=("coh_raw", "mean"), coh_adjusted=("coh_adjusted", "mean"), change=("change", "mean"), change_sd=("change", "std")).reset_index()
real = pd.DataFrame(rr); real.to_csv(OUTD / "posthoc_arm_coherence_real.csv", index=False); simm.to_csv(OUTD / "posthoc_arm_coherence_sim.csv", index=False)
agg = lambda d, cols: d.groupby(["source", "cohort", "arm"], sort=False)[cols].mean().reset_index()
oc = ["n", "offset_sd_total", "offset_plane_share", "offset_plane_var_over_expr_plane_var", "expr_plane_share", "complex_corr", "corr_re", "corr_im"]
pd.concat([agg(pd.DataFrame(off), oc), agg(pd.DataFrame(sim_off), oc)]).to_csv(OUTD / "posthoc_offset_plane.csv", index=False)
rcc = ["corr_cn_residual", "corr_cn_clockfit"]
pd.concat([pd.DataFrame(rc).groupby(["source", "cohort", "arm", "gene"], sort=False)[rcc].mean().reset_index(),
           pd.DataFrame(sim_rc).groupby(["source", "cohort", "arm", "gene"], sort=False)[rcc].mean().reset_index()]).to_csv(OUTD / "posthoc_cn_vs_clock_residual.csv", index=False)
print(real.groupby(["cohort", "arm"])[["coh_raw", "coh_adjusted", "change"]].mean().round(4).to_string()); print(simm.groupby(["cohort", "arm"])[["coh_raw", "coh_adjusted", "change"]].mean().round(4).to_string())

# ---- 5. variants on real data -------------------------------------------------------------------------------------------------------------------------
def slope_nophase(Ylist, CNlist):
    """per-gene beta (ln) of y ~ project + beta CN, no phase covariates (rows with finite CN)."""
    b = []
    for j in range(14):
        num = den = 0.0
        for Y, CN in zip(Ylist, CNlist):
            ok = np.isfinite(CN[:, j]); y, x = Y[ok, j], CN[ok, j]
            if ok.sum() > 1:
                num += ((x - x.mean()) * (y - y.mean())).sum(); den += ((x - x.mean()) ** 2).sum()
        b.append(num / den if den > 0 else 0.0)
    return np.array(b)


adjb = lambda Y, CN, b: Y - b * np.nan_to_num(CN)
vrows, dcell = [], {}
for name, c in COH.items():
    Y, CN = c["Y"], c["CN"]; ph = A.leave_out_phases(Y, CW)
    V = {"raw": Y, "A224 adjusted (BRCA pooled-23 beta, METABRIC cohort beta; with phase)": c["Yadj"], "cohort beta with phase": A.adjust(Y, CN, ph)[0],
         "cohort beta no phase": adjb(Y, CN, slope_nophase([Y], [CN]))}
    V["pooled-23 beta no phase (METABRIC: cohort)"] = adjb(Y, CN, slope_nophase([Y0[p] for p in PROJECTS], [CNM(CNt, tc[p]["sample"].to_numpy()) for p in PROJECTS])) if name == "BRCA" else V["cohort beta no phase"]
    for v, Yv in V.items():
        arms = {a: np.mean([acoh(Yv, cl["idx"][cl["lab"] == l]) for cl in c["cells"]]) for l, a in ((1, "LumA"), (0, "others"))}
        Dc = A.cell_stats(Yv, ONE(Yv), c["cells"], WF, 0, None)[0]; dcell.setdefault(v, []).append(Dc)
        vrows.append(dict(cohort=name, variant=v, D=Dc.mean(), coh_LumA_mean_cells=arms["LumA"], coh_others_mean_cells=arms["others"], n_cells=len(Dc), **{f"D_cell{i}": d for i, d in enumerate(Dc)}))
vrows += [dict(cohort="pooled", variant=v, D=np.concatenate(d).mean(), n_cells=sum(len(x) for x in d)) for v, d in dcell.items()]
vr = pd.DataFrame(vrows)
vr.to_csv(OUTD / "posthoc_variants_real.csv", index=False); print(vr[["cohort", "variant", "D", "coh_LumA_mean_cells", "coh_others_mean_cells"]].round(4).to_string())
say("done posthoc")
