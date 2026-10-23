"""Remove the non-temporal host direction before reading the internal window (A26).

Pre-registered as A26 (2026-10-14, written before running; run note: V26.1 fails, V26.3 fails with g_GTEx, V26.2 passes
unchanged from A25, A26 rejected). The A25 diagnostic found a donor-level state shared across compartments, unrelated to
recorded time, whose clock-plane projection points at about 4 h (mod 12) in GTEx and Chen 2016 and biases the A25 axis.
Handle: (1) learn on timed data the gene vector g: per compartment the gene-centred log expression is regressed on (1, cos t,
sin t), the residuals are averaged over a donor's compartments, g is the leading principal component over donors (GTEx
fast deaths, 42 tissues; Chen, 2 regions, separately); reported: variance share, out-of-plane fraction, in-plane angle
(expected about 4 h mod 12), cosine between g_GTEx and g_Chen. (2) Donor score s_d: mean over the donor's compartments of
g_perp . (I - P) y_i / |g_perp|^2 (y_i gene-centred, centred within its intercept group; g_perp the part of g outside the
clock plane); no recorded time is used on the data being corrected. (3) Correction: per intercept group the clock
coefficients are regressed on s_d and the fitted part removed. Here the correction is applied to the counts, by binomial
thinning with the clock-plane change -s_d a (`simulate.thin_clock`), because `Dataset` takes counts, not states.
Validation (cross-learned): score free of time (R^2 of s_d on cos, sin of recorded time <= 0.10 in GTEx, Chen, colon);
V26.2 GTEx corrected with g_Chen (median axis error <= 2 h and flip correct >= 75 % over the A25 windows with r < 0.6);
V26.3 mouse colon corrected with g_GTEx (the four 3-ZT subsets with r < 0.6 within 2 h of their common offset). Not ported:
V26.1 (Chen windows corrected with g_GTEx: arrays cannot enter `Dataset`) and R26 (the A25 cohorts corrected).

    python compute_nontemporal.py [--subsample N] [--seed S] [--jobs 4]
Outputs (paper/results/figure_1/output/windows/): nontemporal_{g,g_summary,score_r2,gtex,colon}.csv
Legacy: legacy/desync/common/a26_nontemporal.py -> common_output/a26_{g,g_summary,...}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import log_expr  # noqa: E402
from clocksync.simulate import thin_clock  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _internal_time as T  # noqa: E402
import _shared as sh  # noqa: E402
import compute_internal_time as CI  # noqa: E402

FOLDER = "windows"
W, HOURS = T.W, T.HOURS
WF = cs.gtex()
GENES = list(WF.genes)
CHEN_SUNRISE = 6.5


def plane_projector(wf) -> np.ndarray:
    X = wf.X - wf.X.mean(0)
    return X @ np.linalg.solve(X.T @ X, X.T)


def gene_centred(y):
    return y - y.mean(1, keepdims=True)


def centre_groups(y, groups):
    y = y.copy()
    for g in np.unique(groups):
        i = groups == g
        y[i] -= y[i].mean(0)
    return y


def learn_g(y, genes, donor, comp, hours):
    """local: not in clocksync. Leading PC of the donor-averaged, time-removed residuals; its variance share."""
    t = np.asarray(hours, float) * W
    res = np.zeros_like(y)
    for k in np.unique(comp):
        i = comp == k
        X = np.column_stack([np.ones(i.sum()), np.cos(t[i]), np.sin(t[i])])
        res[i] = y[i] - X @ np.linalg.lstsq(X, y[i], rcond=None)[0]
    D = pd.DataFrame(res).groupby(donor).mean().to_numpy().copy()
    D -= D.mean(0)
    _, sv, vt = np.linalg.svd(D, full_matrices=False)
    return pd.Series(vt[0] / np.linalg.norm(vt[0]), index=genes), float(sv[0] ** 2 / (sv ** 2).sum())


def g_perp(g, wf):
    v = g.reindex(list(wf.genes)).fillna(0.0).to_numpy()
    return v - plane_projector(wf) @ v


def in_plane_angle_h(g, wf):
    v = g.reindex(list(wf.genes)).fillna(0.0).to_numpy()
    Xc = wf.X - wf.X.mean(0)
    b = np.linalg.lstsq(Xc, v - v.mean(), rcond=None)[0]
    return float(np.mod(np.arctan2(b[1], b[0]) * HOURS, 12.0))


def donor_score(y, groups, donor, wf, g):
    gp = g_perp(g, wf)
    if not np.any(gp):
        raise ValueError("g has no out-of-plane part on these genes")
    yc = centre_groups(y, groups)
    r = yc - yc @ plane_projector(wf).T
    return pd.Series(r @ gp / (gp @ gp)).groupby(np.asarray(donor)).mean()


def time_r2(s, hours):
    h = hours.reindex(s.index).to_numpy(float) * W
    X = np.column_stack([np.ones(len(h)), np.cos(h), np.sin(h)])
    e = s.to_numpy() - X @ np.linalg.lstsq(X, s.to_numpy(), rcond=None)[0]
    return float(1 - e.var() / s.to_numpy().var())


def correct_counts(ds: cs.Dataset, s: pd.Series, seed: int) -> pd.DataFrame:
    """Pseudobulk table with the clock-plane regression on the donor score removed per intercept group, by thinning."""
    d = ds.data
    x = s.reindex(d.resection).to_numpy(float)
    delta = np.zeros((len(x), 2))
    for g in np.unique(d.group):
        i = d.group == g
        xi = x[i] - np.nanmean(x[i])
        if i.sum() < 4 or not np.any(xi) or np.isnan(xi).any():
            continue
        bi = d.beta[i] - d.beta[i].mean(0)
        delta[i] = -np.outer(xi, xi @ bi / (xi @ xi))
    t = ds.table.copy()
    t[GENES] = thin_clock(t[GENES].to_numpy(float), ds.waveform.X, delta, np.random.default_rng(seed))
    return t.rename(columns={"resection": "donor"})


def gtex_inputs(args):
    pb, donors = sh.gtex_fast(30, args.subsample, args.seed)
    y = gene_centred(log_expr(pb[GENES], pb.depth.to_numpy(float)).to_numpy())
    return pb, donors, y


def chen_inputs():
    s, y, wf = sh.chen_load()
    yc = (y - y.groupby(s.region.to_numpy()).transform("mean")).to_numpy() * np.log(2.0)  # log2 -> natural log
    return s, gene_centred(yc), list(wf.genes), wf


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    pd.set_option("display.width", 220)
    pb, donors, y = gtex_inputs(args)
    tod = donors.tod
    gg, sh_g = learn_g(y, GENES, pb.donor.to_numpy(), pb.compartment.to_numpy(), tod.reindex(pb.donor).to_numpy())
    s, yc, cgenes, cwf = chen_inputs()
    gc, sh_c = learn_g(yc, cgenes, s.donor.to_numpy(), s.region.astype(str).to_numpy(),
                       (s.tod_raw.astype(float) + CHEN_SUNRISE).to_numpy())
    G = pd.DataFrame({"gtex": gg, "chen": gc})
    common = G.dropna().index
    cos = float(G.loc[common, "gtex"] @ G.loc[common, "chen"] / (np.linalg.norm(G.loc[common, "gtex"]) * np.linalg.norm(G.loc[common, "chen"])))
    G.to_csv(out / "nontemporal_g.csv")
    summ = pd.DataFrame([dict(source="gtex", share=sh_g, perp_frac=float(np.linalg.norm(g_perp(gg, WF)) ** 2), angle_h=in_plane_angle_h(gg, WF)),
                         dict(source="chen", share=sh_c, perp_frac=float(np.linalg.norm(g_perp(gc, cwf)) ** 2), angle_h=in_plane_angle_h(gc, cwf))]).assign(cosine_gtex_chen=cos)
    summ.to_csv(out / "nontemporal_g_summary.csv", index=False)
    print(summ.round(3).to_string(index=False))
    # score free of time
    r2 = []
    sg = donor_score(y, pb.compartment.to_numpy(), pb.donor.to_numpy(), WF, gc)
    r2.append(dict(data="gtex_fast", g="chen", r2_time=time_r2(sg, tod)))
    sc = donor_score(yc, s.region.astype(str).to_numpy(), s.donor.to_numpy(), cwf, gg)
    r2.append(dict(data="chen", g="gtex", r2_time=time_r2(sc, (s.groupby("donor").tod_raw.first().astype(float) + CHEN_SUNRISE))))
    cpb, cdon = data.load("colon_chrono", subsample=args.subsample, seed=args.seed)
    cy = gene_centred(log_expr(cpb[GENES], cpb.depth.to_numpy(float)).to_numpy())
    cgroup = (cpb.compartment + "|" + cpb.donor.map(cdon.arm)).to_numpy()
    scol = donor_score(cy, cgroup, cpb.donor.to_numpy(), WF, gg)
    r2.append(dict(data="colon", g="gtex", r2_time=time_r2(scol, cdon.zt.astype(float))))
    R = pd.DataFrame(r2).assign(accept=lambda d: d.r2_time <= 0.10)
    R.to_csv(out / "nontemporal_score_r2.csv", index=False)
    print(R.round(3).to_string(index=False))
    # V26.2: GTEx corrected with g_Chen, the A25 windows
    ds = cs.Dataset(pb, anchors=None, donors=donors, name="gtex_fast")
    tab = correct_counts(ds, sg, args.seed)
    dsc = cs.Dataset(tab, anchors=None, donors=donors, name="gtex_fast_corrected")
    t = tod.dropna() % 24
    anchors = list(dsc.compartments)
    centres, widths = (range(0, 24, 6), (9, 18)) if run["quick"] else (range(0, 24, 3), (6, 9, 12, 18))
    sets = [("all", np.nan, np.nan, t.index.to_numpy())] + [
        (f"{c:g}h+-{w / 2:g}", c, w, t.index[np.abs(T.cdiff(t.to_numpy(), c)) <= w / 2].to_numpy()) for c in centres for w in widths]
    rows = []
    for tag, d_ in (("uncorrected", ds), ("corrected_g_chen", dsc)):
        rr = Parallel(n_jobs=args.jobs)(delayed(CI.window_row)(d_, keep, t, dict(data="gtex_fast", version=tag, window=nm, c_set=c, w_set=w), anchors)
                                        for nm, c, w, keep in sets)
        rows += rr
    D = pd.DataFrame(rows)
    D.to_csv(out / "nontemporal_gtex.csv", index=False)
    for tag, x in D[D.status == "ok"].groupby("version"):
        g = x[(x.r < T.R_MAX) & x.window.ne("all")]
        if len(g):
            conc = float(abs(np.exp(1j * 2 * np.pi * x[x.window.ne("all")].c_axis12 / 12).mean()))
            print(f"V26.2 {tag}: windows with r < 0.6: {len(g)}; median |axis error| {g.err_axis12.abs().median():.2f} h (<= 2); "
                  f"flip correct {g.flip_ok.mean():.2f} (>= 0.75); concentration of axis centres {conc:.2f}")
    # V26.3: mouse colon corrected with g_GTEx
    cds = cs.Dataset(cpb, anchors=CI.COLON_ANCHORS, donors=cdon, strata="arm", name="colon_chrono")
    # the score is computed on all compartments of the table; the host uses the anchors
    ctab = correct_counts(cds, scol, args.seed)
    cdc = cs.Dataset(ctab, anchors=CI.COLON_ANCHORS, donors=cdon, strata="arm", name="colon_corrected")
    zt = cdon.zt.astype(float)
    crow, zts = [], (1.0, 7.0, 13.0, 19.0)
    for tag, d_ in (("uncorrected", cds), ("corrected_g_gtex", cdc)):
        for i in range(4):
            pick = [zts[(i + j) % 4] for j in range(3)]
            keep = zt.index[zt.isin(pick)].to_numpy()
            crow.append(CI.window_row(d_, keep, zt, dict(version=tag, zts="/".join(f"{p:g}" for p in pick)), CI.COLON_ANCHORS,
                                      min_donors=6, min_comp_donors=4))
    C = pd.DataFrame(crow)
    for tag, x in C[C.status == "ok"].groupby("version"):
        off = np.angle(np.exp(1j * W * x.err_curv).mean()) * HOURS
        C.loc[x.index, "err_relative"] = T.cdiff(x.err_curv, off)
    C.to_csv(out / "nontemporal_colon.csv", index=False)
    print(C[["version", "zts", "status", "r", "err_curv", "err_relative"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
