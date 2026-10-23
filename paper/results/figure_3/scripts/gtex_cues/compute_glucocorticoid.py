"""A27 part 2 (post hoc, disclosed): ventilated donors keep an internal clock that is decoupled from the recorded time;
which genes keep following it, and is a glucocorticoid (GC) axis the remaining zeitgeber?

Amendment A27 (docs/PLAN_internal_desync.md (archive); exploratory, not pre-registered, written after running). Parts here:
  D6   a priori gene classes among the genes rhythmic in fast deaths (timed cosinor p < 0.01, not clock plane), read in
       ventilated donors on the host of the other tissues: retention (still follows internal time, p_int < 1e-3),
       phase kept (|internal - fast phase| < 3 h after the frame rotation), amplitude ratio relative to the clock-plane
       genes of the tissue, follows time of death (p_ext < 0.01), each against matched random sets (same tissue, same
       p_fast decile; 2,000 draws; two-sided empirical p). Classes: clock output, GC response, heat shock, SREBP,
       mTORC1, immediate-early, hypoxia, UPR (host arm only).
  D6b  the GC score (mean within-tissue z of FKBP5, TSC22D3, ZBTB16, KLF9, KLF15, DDIT4, log CPM) by Hardy class; R^2 of
       the donor GC score on the host state and on the time of death; host-time-of-death alignment by GC tertile.
  D6c  GC effect on the 14 clock genes (log CPM ~ tissue + GC score, pooled): fraction of the effect inside the GTEx
       waveform plane (span{1, u, v}) against random donor scores, and the phase it pulls towards; without PER1.
Aggregates only (time of death is controlled access; only per-class, per-tissue and pooled quantities are written).

    python compute_glucocorticoid.py [--subsample N] [--jobs 1]
Outputs: paper/results/figure_3/output/gtex_cues/diag6_classes.csv, diag6_tissues.csv, diag6b_gc.csv,
diag6c_gc_summary.csv, diag6c_gc_clock_effects.csv
Legacy: legacy/desync/validation/glucocorticoid/scripts/{diag5_int_ext (gene table, host arm),diag6_classes,diag6b_gc,
diag6c_gc_zeitgeber}.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _covariates as K  # noqa: E402
import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

W = S.W
WF = cs.gtex()
CK = list(WF.genes)
FRAME = set(CK) | {"BMAL1"}
GMT = [Path("/home/maxine/Documents/paychere/CoPhaser/test/hazel/data/h.all.v2026.1.Hs.symbols.gmt"),
       DATA / "reference" / "h.all.v2024.1.Hs.symbols.gmt"]
TARGETS = [t for t in S.TISSUES if t != "Pituitary" and not t.startswith("Brain")]


def hallmark() -> dict:
    f = next(p for p in GMT if p.exists())
    return {l.split("\t")[0]: set(l.strip().split("\t")[2:]) for l in open(f)}


def classes() -> dict:
    H = hallmark()
    return {
        "clock output (non-frame)": {"HLF", "NFIL3", "BHLHE40", "BHLHE41", "RORC", "NPAS2", "ARNTL2", "CRY2", "TSPAN4",
                                     "USP2", "LONRF3", "TEF", "NAMPT", "KLF10", "DDIT4L", "WEE1", "PPARA", "CIPC"},
        "glucocorticoid response": {"FKBP5", "TSC22D3", "ZBTB16", "KLF9", "KLF15", "DDIT4", "SGK1", "DUSP1", "CDKN1C",
                                    "MT2A", "MT1X", "PDK4", "ERRFI1", "IRS2", "CEBPD", "SESN1", "FOXO3", "GLUL", "TXNIP",
                                    "PER1", "ALOX15B", "MAOA", "IL1R2", "CD163", "SLC19A2", "RHOB", "FAM107A", "SULT1A1"},
        "temperature / heat shock": {"HSPA1A", "HSPA1B", "HSPA8", "HSPH1", "HSP90AA1", "HSP90AB1", "DNAJB1", "DNAJA1",
                                     "HSPB1", "HSPD1", "HSPE1", "BAG3", "CIRBP", "RBM3", "SERPINH1", "CHORDC1", "HSPA4L",
                                     "HSPA6", "DNAJA4", "CRYAB", "FKBP4", "STIP1"},
        "feeding: SREBP / cholesterol": H["HALLMARK_CHOLESTEROL_HOMEOSTASIS"],
        "nutrient: mTORC1": H["HALLMARK_MTORC1_SIGNALING"],
        "immediate-early / activity": {"NR4A1", "NR4A2", "NR4A3", "FOS", "FOSB", "JUN", "JUNB", "EGR1", "EGR2", "EGR3",
                                       "ATF3", "PPARGC1A", "MYC", "ANKRD1", "XIRP1", "HES1", "IER2", "IER3", "IER5", "MAFF",
                                       "KLF2", "KLF4", "CYR61", "CCN1", "CTGF", "CCN2", "GADD45B", "PPP1R15A", "BTG2",
                                       "ZFP36"},
        "hypoxia": H["HALLMARK_HYPOXIA"],
        "unfolded protein response": H["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"],
    }


# ---- D6: gene table of the ventilated arm (host), as the legacy diag5 ------------------------------------------------

def gene_table(tissue: str, pbv: pd.DataFrame, subs=None) -> pd.DataFrame:
    """Per gene of one tissue (local: not in clocksync): p_int / p_ext (partial F of the host state given the time of
    death, and the reverse) and loadings in ventilated donors; the timed cosinor of fast deaths (p_fast, amplitude,
    phase). Genes with median CPM >= 5 in both cohorts; log(CPM + 0.5); D2 covariates. The host is the anchored host of
    the other 14 tissues in ventilated donors (`Dataset.fit(...).host`)."""
    Yf, gf, obf = S.expression(tissue, hardy=(1, 2), median_cpm=None)
    Yv, gv, obv = S.expression(tissue, hardy=(0,), median_cpm=None)
    if subs is not None:
        mf, mv = obf.subject.isin(subs).to_numpy(), obv.subject.isin(subs).to_numpy()
        Yf, obf, Yv, obv = Yf[mf], obf[mf].reset_index(drop=True), Yv[mv], obv[mv].reset_index(drop=True)
    cpm_f, cpm_v = np.exp(Yf) - 0.5, np.exp(Yv) - 0.5
    ok = (np.median(cpm_f, 0) >= 5) & (np.median(cpm_v, 0) >= 5)
    Yf, Yv, genes = Yf[:, ok], Yv[:, ok], gf[ok]
    host = S.host_states(pbv, [x for x in S.TISSUES if x != tissue], name=f"vent_wo_{tissue}")
    subj = obv.subject.to_numpy()
    m = np.isin(subj, host.index)
    Yv, obv, subj = Yv[m], obv[m].reset_index(drop=True), subj[m]
    f = host.loc[subj, ["x", "y"]].to_numpy()
    f = f - f.mean(0)
    t = obv.tod.to_numpy(float)
    E = np.c_[np.cos(W * t), np.sin(W * t)]
    Z = S.d2_covariates(obv)
    Z = Z[:, Z.std(0) > 0]
    D = np.c_[np.ones(len(Yv)), Z, f, E]
    k = D.shape[1]
    B, p_int = S.partial_f(Yv, D, [k - 4, k - 3])
    _, p_ext = S.partial_f(Yv, D, [k - 2, k - 1])
    c_int = B[-4] + 1j * B[-3]
    c_fast, p_fast = S.timed_cosinor(Yf, obf.tod.to_numpy(float), S.d2_covariates(obf))
    fr = np.isin(genes, list(FRAME))
    rot = np.angle(np.exp(1j * (np.angle(c_int[fr]) - np.angle(c_fast[fr]))).mean())
    return pd.DataFrame(dict(tissue=tissue, gene=genes, frame=fr, p_int=p_int, p_ext=p_ext, p_fast=p_fast,
                             amp_int=np.abs(c_int), amp_fast=np.abs(c_fast), int_h=(np.angle(c_int) - rot) / W % 24,
                             fast_h=np.angle(c_fast) / W % 24))


def metrics(d: pd.DataFrame) -> pd.Series:
    dphi = np.abs(np.angle(np.exp(1j * W * (d.int_h - d.fast_h)))) / W
    return pd.Series(dict(n=len(d), retained=(d.p_int < 1e-3).mean(), phase_kept=((d.p_int < 1e-3) & (dphi < 3)).mean(),
                          amp_rel=np.nanmedian(d.amp_rel), ext=(d.p_ext < 0.01).mean()))


def class_tests(g: pd.DataFrame, n_draw: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    g = g.copy()
    g["amp_ratio"] = g.amp_int / g.amp_fast
    fr = g[g.frame & (g.p_fast < 0.01)].groupby("tissue").amp_ratio.median()
    g["amp_rel"] = g.amp_ratio / g.tissue.map(fr)
    r = g[(g.p_fast < 0.01) & ~g.frame].copy().reset_index(drop=True)
    r["bin"] = r.groupby("tissue").p_fast.transform(lambda x: pd.qcut(x.rank(method="first"), 10, labels=False))
    rows = []
    clockset = set(CK)
    for cname, genes in list(classes().items()) + [("all rhythmic", None)]:
        if genes is not None:
            genes = set(genes) - (clockset - {"NPAS2", "CRY2", "TEF"} if cname == "clock output (non-frame)" else clockset)
        inc = r if genes is None else r[r.gene.isin(genes)]
        obs = metrics(inc)
        out = dict(arm="host", cls=cname, **obs.to_dict(), tissues=inc.tissue.nunique())
        if genes is not None and len(inc):
            idx = {k: np.asarray(v) for k, v in r[~r.gene.isin(genes)].groupby(["tissue", "bin"]).groups.items()}
            keys = [k for k in zip(inc.tissue, inc.bin) if k in idx]
            null = pd.DataFrame([metrics(r.loc[[rng.choice(idx[k]) for k in keys]]) for _ in range(n_draw)])
            for m in ("retained", "phase_kept", "amp_rel", "ext"):
                out[f"{m}_null"] = null[m].median()
                lo, hi = (1 + (null[m] <= obs[m]).sum()) / (n_draw + 1), (1 + (null[m] >= obs[m]).sum()) / (n_draw + 1)
                out[f"{m}_p"] = min(1.0, 2 * min(lo, hi))
        rows.append(out)
    return pd.DataFrame(rows)


# ---- D6b / D6c ---------------------------------------------------------------------------------------------------------

def gc_scores(pb: pd.DataFrame) -> pd.Series:
    """Donor GC score: mean over tissues of the within-tissue z of mean log CPM of the six GC genes (z over all samples
    of the tissue, Hardy 0-2)."""
    lg = pd.DataFrame(np.log(pb[list(S.GC)].to_numpy(float) / pb.depth.to_numpy(float)[:, None] * 1e6 + 0.5),
                      columns=list(S.GC))
    z = lg.groupby(pb.compartment.to_numpy()).transform(lambda x: (x - x.mean()) / x.std()).mean(axis=1)
    return z.groupby(pb.donor.to_numpy()).mean()


def d6b(pb, don, gc) -> pd.DataFrame:
    rows = []
    for h in (0, 1, 2):
        x = gc[gc.index.isin(don.index[don.hardy == h])]
        rows.append(dict(quantity="GC score by Hardy", group=f"Hardy {h}", n=len(x), mean=x.mean(), sd=x.std()))
    for lab, hd in (("fast deaths", (1, 2)), ("ventilator", (0,))):
        sub = pb[pb.donor.isin(don.index[don.hardy.isin(hd)])]
        host = S.host_states(sub, S.TISSUES, name=lab)
        s = host.index.intersection(gc.index)
        f = host.loc[s, ["x", "y"]].to_numpy()
        g = gc.loc[s].to_numpy()
        t = W * don.tod.reindex(s).to_numpy(float)
        r2 = lambda D, y: 1 - ((y - D @ np.linalg.lstsq(D, y, rcond=None)[0]) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        ok = np.isfinite(t)
        phi = np.arctan2(f[:, 1], f[:, 0])
        terc = pd.qcut(g, 3, labels=False)
        row = dict(quantity="GC score vs host and time of death", group=lab, n=len(s), r2_host=r2(np.c_[np.ones(len(g)), f], g),
                   r2_time_of_death=r2(np.c_[np.ones(ok.sum()), np.cos(t[ok]), np.sin(t[ok])], g[ok]),
                   R_host_tod_overall=float(abs(np.exp(1j * (phi[ok] - t[ok])).mean())))
        for k in range(3):
            m = (terc == k) & ok
            row[f"R_host_tod_gc_tertile{k + 1}"] = float(abs(np.exp(1j * (phi[m] - t[m])).mean()))
        rows.append(row)
    return pd.DataFrame(rows)


def d6c(pb, don, gc, n_null, seed):
    """GC effect on the 14 clock genes, pooled within tissue; non-brain tissues (as the legacy run)."""
    rng = np.random.default_rng(seed)
    U = np.asarray(WF.X)
    Pb = np.c_[np.ones(len(CK)), U]
    Q, _ = np.linalg.qr(Pb)
    inplane = lambda b: float(((Q @ (Q.T @ b)) ** 2).sum() / (b ** 2).sum())
    pull = lambda P_, b: (lambda c: (float(np.arctan2(c[2], c[1]) / W % 24), float(np.hypot(c[1], c[2]))))(
        np.linalg.lstsq(P_, b, rcond=None)[0])
    lg = pd.DataFrame(np.log(pb[CK].to_numpy(float) / pb.depth.to_numpy(float)[:, None] * 1e6 + 0.5), columns=CK)
    lg["tissue"], lg["donor"] = pb.compartment.to_numpy(), pb.donor.to_numpy()
    lg["vent"] = (don.hardy.reindex(lg.donor).to_numpy() == 0)
    lg["gc"] = gc.reindex(lg.donor).to_numpy()
    summ, genes = [], []
    for lab, v in (("fast deaths", False), ("ventilator", True)):
        d = lg[lg.vent == v].dropna(subset=["gc"])
        Y = d[CK].to_numpy() - d.groupby("tissue")[CK].transform("mean").to_numpy()
        g = d.gc.to_numpy() - d.groupby("tissue").gc.transform("mean").to_numpy()
        b = (g @ Y) / (g @ g)
        us = d.donor.unique()
        sub = d.donor.to_numpy()
        nul = []
        for _ in range(n_null):
            rs = pd.Series(rng.standard_normal(len(us)), index=us)
            gr = rs.loc[sub].to_numpy()
            gr = gr - pd.Series(gr).groupby(d.tissue.to_numpy()).transform("mean").to_numpy()
            nul.append(inplane((gr @ Y) / (gr @ gr)))
        ph_h, mag = pull(Pb, b)
        resid = Y - np.outer(g, b)
        se = np.sqrt((resid ** 2).sum(0) / (len(g) - 1 - d.tissue.nunique()) / (g @ g))
        genes += [dict(group=lab, gene=c, estimate=b[i], se=se[i], n_samples=len(g), n_donors=len(us))
                  for i, c in enumerate(CK)]
        m = np.array([c != "PER1" for c in CK])
        Qm, _ = np.linalg.qr(Pb[m])
        bm = b[m]
        summ.append(dict(group=lab, n_donors=len(us), inplane=inplane(b), inplane_null_median=float(np.median(nul)),
                         inplane_null_q95=float(np.quantile(nul, 0.95)), pull_h=ph_h, pull_magnitude=mag,
                         inplane_without_PER1=float(((Qm @ (Qm.T @ bm)) ** 2).sum() / (bm ** 2).sum()),
                         pull_h_without_PER1=float(pull(Pb[m], bm)[0])))
    return pd.DataFrame(summ), pd.DataFrame(genes)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(K.FOLDER, run["quick"])
    n_draw, n_null = (100, 20) if run["quick"] else (2000, 200)

    # donors: all Hardy 0-2 donors (or a random subsample), clock and GC genes
    pb, don = gtex.load(tissues=S.TISSUES, hardy=(0, 1, 2), genes=CK + list(S.GC), depth="library_size",
                        require_tod=False)
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[:args.subsample]
        pb, don = pb[pb.donor.isin(keep)].reset_index(drop=True), don.loc[don.index.isin(keep)]
    subs = set(pb.donor) if args.subsample else None
    cols = ["donor", "compartment", "depth"] + CK
    pbv = pb[pb.donor.isin(don.index[don.hardy == 0])][cols].reset_index(drop=True)

    # D6: the ventilated arm gene table over the 12 body tissues
    tabs = []
    for t in TARGETS:
        g = gene_table(t, pbv, subs=subs)
        tabs.append(g)
        print(t, len(g), int(((g.p_fast < 0.01) & ~g.frame).sum()), "rhythmic non-frame genes", flush=True)
    G = pd.concat(tabs, ignore_index=True)
    res = class_tests(G, n_draw, 0)
    res.to_csv(out / "diag6_classes.csv", index=False)
    pd.set_option("display.width", 250)
    print(res[["cls", "n", "tissues", "retained", "retained_null", "retained_p", "phase_kept", "phase_kept_null",
               "phase_kept_p", "amp_rel", "amp_rel_p", "ext", "ext_p"]].round(3).to_string(index=False), flush=True)
    per_t = G[(G.p_fast < 0.01) & ~G.frame].groupby("tissue").agg(
        rhythmic_fast=("gene", "size"), retained=("p_int", lambda p: (p < 1e-3).mean()), ext=("p_ext", lambda p: (p < 0.01).mean()))
    per_t.to_csv(out / "diag6_tissues.csv")

    # D6b, D6c: GC score on the Hardy 0-2 samples of the host tissues
    gc = gc_scores(pb)
    d6b(pb[cols + list(S.GC)], don, gc).to_csv(out / "diag6b_gc.csv", index=False)
    nb = pb[~pb.compartment.str.startswith("Brain")]
    summ, genes = d6c(nb, don, gc, n_null, 0)
    summ.to_csv(out / "diag6c_gc_summary.csv", index=False)
    genes.to_csv(out / "diag6c_gc_clock_effects.csv", index=False)
    print(pd.read_csv(out / "diag6b_gc.csv").round(3).to_string(index=False), "\n", summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
