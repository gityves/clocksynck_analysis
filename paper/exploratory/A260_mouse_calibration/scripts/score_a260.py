"""A260: does the untimed procedure put mouse liver genes in the right class, with and without a local clock and with
feeding moved? (Manella 2021, GSE159135)

Registered in paper/docs/registered_experiments.md (A260). Truth: timed cosinor per arm and tissue (the scout's code,
`scouts_2026-11-15/T0c_calibration/scout_t0c.py`). Untimed: `Dataset.gene_rhythms` of the liver genome on the host of
WAT + lung (V260.1, V260.2) or lung alone (V260.3), one arm at a time, 50 host-permuted nulls. D260.4: A259's kappa of
the liver and per-gene loadings on the liver's deviation from its WAT + lung prediction.

    python score_a260.py [--jobs 4]
Outputs (output/): untimed_genes.csv.gz, truth_genes.csv.gz, tests.csv, descriptive.csv.

Choices not fixed by the registration
- Host of two anchors (WAT + lung) passed to `rhythms.gene_rhythms` (the liver's host; the dataset's own leave-one-out
  hosts of the anchors would have one compartment each). V260.3, "host from the lung alone": a host needs >= 2
  compartments, so every gene is regressed (OLS, intercept) on the lung's own clock-plane state (`deviation.states`),
  with the same Wald p and 50 permutations of the lung states as null.
- Clock genes failing the A111 QC (median count >= 5 in every tissue over all mice) or absent from a file (CLOCK in the
  knockout lung, as A111) are dropped from the waveform.
- Genome depth = sum of the retained genes of the file (as the builder's clock table and the scout's truth).
- AUC by the Mann-Whitney statistic; 90 % interval from 2,000 bootstraps of genes within class (seed 260).
- Circular correlation: Fisher-Lee (`clocksync.stats.fisher_lee`).
- D260.4 per-gene deviation loadings: liver log fraction ~ stratum-free OLS on the host state (WAT + lung composite, both
  rotated as in `clocksync.deviation`) and the liver's residual state (liver state minus its least-squares prediction
  from the host); Wald p of the residual's two columns.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "common/datasets/manella2021"))
sys.path.insert(0, str(PAPER / "exploratory/scouts_2026-11-15/T0c_calibration"))
import clocksync as cs  # noqa: E402
from clocksync import deviation as DV  # noqa: E402
from clocksync import rhythms  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync.rhythms import frame_shift  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import data  # noqa: E402
import build as M  # noqa: E402
from scout_t0c import CLOCK, cosinor  # noqa: E402

OUT = HERE.parent / "output"
SEED = 260
W = 2 * np.pi / 24


def genome(genotype: str, feeding: str) -> pd.DataFrame:
    x = M.read_tissue("Liver", genotype)
    x = x.groupby(x.index.map(canonical)).sum()
    x = x[[c for c in x.columns if f"_{feeding}_" in c]]
    g = x.T.reset_index(names="donor")
    g.insert(1, "compartment", "Liver")
    g.insert(2, "stratum", f"{genotype}_{feeding}")
    g.insert(3, "depth", x.sum().to_numpy(float))
    return g


def truth() -> dict:
    res = {}
    for t in M.TISSUES:
        for g in M.GENOTYPES:
            x = M.read_tissue(t, g)
            x = x.groupby(x.index.map(canonical)).sum()
            for f in ("AL", "DF"):
                cols = [c for c in x.columns if f"_{f}_" in c]
                zt = np.array([int(c.split("_ZT")[1].split("_")[0]) for c in cols], float)
                res[(t, g, f)] = cosinor(x[cols], zt)
    return res


def classes(res) -> pd.Series:
    j = res[("Liver", "AlbCre", "AL")].join(res[("Liver", "BLKO", "AL")], lsuffix="_wt", rsuffix="_ko", how="inner")
    j = j[~j.index.isin(CLOCK)]
    ratio = j.amp_ko / j.amp_wt
    c = pd.Series("other", index=j.index)
    c[(j.q_wt < 0.05) & (j.p_ko > 0.2) & (ratio < 0.5)] = "local"
    c[(j.q_wt < 0.05) & (j.q_ko < 0.05) & (ratio >= 0.5)] = "systemic"
    c[(j.q_ko < 0.05) & (j.p_wt > 0.2)] = "ko_only"
    return c


def auc(score_pos, score_neg) -> float:
    u = stats.mannwhitneyu(score_pos, score_neg, alternative="two-sided").statistic
    return float(u / (len(score_pos) * len(score_neg)))


def auc_ci(pos, neg, rng, n=2000):
    pos, neg = np.asarray(pos), np.asarray(neg)
    b = [auc(rng.choice(pos, len(pos)), rng.choice(neg, len(neg))) for _ in range(n)]
    return auc(pos, neg), *np.quantile(b, [0.05, 0.95])


def untimed(pb, wf, arm, anchors, n_null) -> pd.DataFrame:
    """Per-gene untimed loadings of the liver. Two anchors: `Dataset.gene_rhythms` with the anchors' host passed in
    (its own leave-one-out hosts of the anchors would have one anchor each). One anchor (lung): OLS of every gene on
    the anchor's own clock-plane state (a host needs >= 2 compartments)."""
    g, f = arm.split("_")
    sub = pb[pb.stratum == arm].reset_index(drop=True)
    gen = genome(g, f)
    if len(anchors) >= 2:
        ds = cs.Dataset(sub, anchors=anchors, waveform=wf, name=arm)
        full = A.host(ds.data, anchors)
        H = ({k: full for k in ds.compartments}, frame_shift(full))
        r = rhythms.gene_rhythms(ds, gen, n_null=n_null, seed=SEED, hosts=H)
        r = r[r.compartment == "Liver"]
    else:
        S = DV.states(sub, wf).unstack("compartment")[anchors[0]]
        q = gen.set_index("donor").loc[S.index]
        cols = [c for c in q.columns if c not in ("compartment", "stratum", "depth")]
        ok = (q[cols].mean() >= 5).to_numpy()
        Y = np.log((q[cols].loc[:, ok].to_numpy(float) + 0.5) / q.depth.to_numpy(float)[:, None])
        Sv = S.to_numpy()
        fr = np.c_[Sv.real - Sv.real.mean(), Sv.imag - Sv.imag.mean()]
        rng = np.random.default_rng(SEED)
        out = []
        for b in range(-1, n_null):
            fb = fr if b < 0 else fr[rng.permutation(len(fr))]
            D = np.c_[np.ones(len(fb)), fb]
            B, *_ = np.linalg.lstsq(D, Y, rcond=None)
            s2 = ((Y - D @ B) ** 2).sum(0) / (len(Y) - 3)
            Xi = np.linalg.pinv(D.T @ D)
            se = np.sqrt(s2 * (Xi[1, 1] + Xi[2, 2]) / 2)
            c = B[1] + 1j * B[2]
            out.append(pd.DataFrame(dict(compartment="Liver", gene=np.asarray(cols)[ok], amplitude=np.abs(c),
                                         acrophase_h=np.mod(np.angle(c), 2 * np.pi) / W, se=se,
                                         p=np.exp(-np.abs(c) ** 2 / (2 * se ** 2)), null=b)))
        r = pd.concat(out, ignore_index=True)
    obs = r[r.null == -1].set_index("gene")
    nul = r[r.null >= 0]
    obs["fdr"] = [min(1.0, (nul.p <= p).sum() / max(n_null, 1) / max((obs.p <= p).sum(), 1)) for p in obs.p]
    return obs.assign(arm=arm, anchors="+".join(anchors))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=50)
    a = ap.parse_args()
    limit_threads(4)
    warnings.simplefilter("ignore", UserWarning)
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    pb, _ = data.load("manella2021", allow_stale=True)
    genes = [c for c in pb.columns if c not in ("donor", "compartment", "stratum", "depth")]
    med = pb.groupby("compartment")[genes].median()
    keep = [g for g in genes if (med[g] >= 5).all() and pb[g].notna().all()]
    wf = cs.gtex().subset([g for g in cs.gtex().genes if g in keep])
    print("clock genes kept:", len(wf.genes), "dropped:", sorted(set(genes) - set(keep)), flush=True)

    res = truth()
    cls = classes(res)
    tt = pd.concat({f"{t}|{g}|{f}": r for (t, g, f), r in res.items()}).rename_axis(["arm", "gene"]).reset_index()
    tt[tt.arm.str.startswith("Liver") | tt.gene.isin(CLOCK)].to_csv(OUT / "truth_genes.csv.gz", index=False)

    U = {}
    for arm in ("AlbCre_AL", "BLKO_AL"):
        U[(arm, "WAT+Lung")] = untimed(pb, wf, arm, ["WAT", "Lung"], a.n_null)
        print(arm, "WAT+Lung done", flush=True)
    for arm in ("AlbCre_AL", "AlbCre_DF", "BLKO_AL", "BLKO_DF"):
        U[(arm, "Lung")] = untimed(pb, wf, arm, ["Lung"], a.n_null)
        print(arm, "Lung done", flush=True)
    pd.concat(U.values()).reset_index().to_csv(OUT / "untimed_genes.csv.gz", index=False)

    rows = []
    # V260.1: knockout AL, -log10 p separates systemic from local
    ko = U[("BLKO_AL", "WAT+Lung")]
    s = -np.log10(ko.p.clip(lower=1e-300))
    pos, neg = s.reindex(cls.index[cls == "systemic"]).dropna(), s.reindex(cls.index[cls == "local"]).dropna()
    v, lo, hi = auc_ci(pos, neg, rng)
    rows.append(dict(test="V260.1", readout="AUC systemic vs local, -log10 p untimed (knockout AL)", value=v, lo90=lo,
                     hi90=hi, n_pos=len(pos), n_neg=len(neg), passed=bool(v >= 0.80 and lo > 0.65)))
    # V260.2: untimed log amplitude ratio knockout / control
    wt = U[("AlbCre_AL", "WAT+Lung")]
    lr_u = np.log(ko.amplitude / wt.amplitude.reindex(ko.index))
    j = res[("Liver", "AlbCre", "AL")].join(res[("Liver", "BLKO", "AL")], lsuffix="_wt", rsuffix="_ko", how="inner")
    j = j[~j.index.isin(CLOCK)]
    lr_t = np.log(j.amp_ko / j.amp_wt)
    sel = j.index[(j.q_wt < 0.05) | (j.q_ko < 0.05)]
    d = pd.DataFrame(dict(u=lr_u.reindex(sel), t=lr_t.reindex(sel))).replace([np.inf, -np.inf], np.nan).dropna()
    rho = stats.spearmanr(d.u, d.t).statistic
    bs = [stats.spearmanr(*d.sample(len(d), replace=True, random_state=int(rng.integers(2**31))).T.to_numpy()).statistic
          for _ in range(2000)]
    rows.append(dict(test="V260.2a", readout="Spearman untimed vs timed log amplitude ratio knockout / control",
                     value=rho, lo90=np.quantile(bs, 0.05), hi90=np.quantile(bs, 0.95), n_pos=len(d),
                     passed=bool(rho >= 0.60)))
    pos, neg = lr_u.reindex(cls.index[cls == "systemic"]).dropna(), lr_u.reindex(cls.index[cls == "local"]).dropna()
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    v, lo, hi = auc_ci(pos, neg, rng)
    rows.append(dict(test="V260.2b", readout="AUC systemic vs local, untimed log amplitude ratio", value=v, lo90=lo,
                     hi90=hi, n_pos=len(pos), n_neg=len(neg), passed=bool(v >= 0.80 and lo > 0.65)))
    # V260.3: shift reader on the lung host
    for g in ("AlbCre", "BLKO"):
        ual, udf = U[(f"{g}_AL", "Lung")], U[(f"{g}_DF", "Lung")]
        tal, tdf = res[("Liver", g, "AL")], res[("Liver", g, "DF")]
        la, ld = res[("Lung", g, "AL")], res[("Lung", g, "DF")]
        lk = la.join(ld, lsuffix="_al", rsuffix="_df", how="inner")
        lk = lk[(lk.q_al < 0.05) & (lk.q_df < 0.05) & lk.index.isin(CLOCK)]
        lung_shift = float(np.median((lk.phase_h_df - lk.phase_h_al + 12) % 24 - 12)) if len(lk) else np.nan
        k = tal.join(tdf, lsuffix="_al", rsuffix="_df", how="inner")
        k = k[(k.q_al < 0.05) & (k.q_df < 0.05) & ~k.index.isin(CLOCK)]
        k = k[k.index.isin(ual.index) & k.index.isin(udf.index)]
        t_shift = (k.phase_h_df - k.phase_h_al - lung_shift + 12) % 24 - 12
        u_shift = (udf.acrophase_h.reindex(k.index) - ual.acrophase_h.reindex(k.index) + 12) % 24 - 12
        ok = t_shift.notna() & u_shift.notna()
        fl = cs.stats.fisher_lee(t_shift[ok].to_numpy() * W, u_shift[ok].to_numpy() * W) if ok.sum() >= 5 else np.nan
        idx = np.flatnonzero(ok.to_numpy())
        bs = [cs.stats.fisher_lee(t_shift[ok].to_numpy()[b] * W, u_shift[ok].to_numpy()[b] * W)
              for b in (rng.choice(idx.size, idx.size) for _ in range(2000))] if ok.sum() >= 5 else [np.nan]
        rows.append(dict(test="V260.3" if g == "AlbCre" else "V260.3 (knockout, reported)",
                         readout=f"circular correlation untimed vs timed DF-AL shift relative to lung ({g})", value=fl,
                         lo90=np.nanquantile(bs, 0.05), hi90=np.nanquantile(bs, 0.95), n_pos=int(ok.sum()),
                         lung_clock_shift_h=lung_shift,
                         median_abs_diff_h=float(np.median(np.abs(((u_shift - t_shift)[ok] + 12) % 24 - 12))) if ok.sum() else np.nan,
                         passed=bool(fl >= 0.70) if g == "AlbCre" else None))
    tests = pd.DataFrame(rows)
    tests.to_csv(OUT / "tests.csv", index=False)

    # D260.4: kappa of the liver and per-gene loadings on its deviation
    desc = []
    for arm in ("AlbCre_AL", "BLKO_AL"):
        sub = pb[pb.stratum == arm].reset_index(drop=True)
        r = DV.incoherence(sub, "Liver", ["WAT", "Lung"], waveform=wf, n_boot=1000, seed=SEED)
        desc.append(dict(item=f"kappa liver vs WAT+Lung, {arm}", value=r.kappa, lo90=r.lo, hi90=r.hi, hours=r.hours,
                         hours_hi=r.hours_hi, n=r.donors))
        S = DV.states(sub, wf).unstack("compartment").dropna()
        rot = lambda c: (lambda v: v / abs(v))(np.mean((S["WAT"] - S["WAT"].mean()) * np.conj(S[c] - S[c].mean())))
        h = (S["WAT"] + S["Lung"] * rot("Lung")) / 2
        X = np.c_[h.to_numpy().real, h.to_numpy().imag]
        L = np.c_[S["Liver"].to_numpy().real, S["Liver"].to_numpy().imag]
        Bh, *_ = np.linalg.lstsq(np.c_[np.ones(len(X)), X], L, rcond=None)
        dev = L - np.c_[np.ones(len(X)), X] @ Bh
        g, f = arm.split("_")
        gt = genome(g, f).set_index("donor").loc[S.index]
        cols = [c for c in gt.columns if c not in ("compartment", "stratum", "depth")]
        ok = gt[cols].mean() >= 5
        Y = np.log((gt[cols].loc[:, ok].to_numpy(float) + 0.5) / gt.depth.to_numpy(float)[:, None])
        D = np.c_[np.ones(len(X)), X, dev]
        B, *_ = np.linalg.lstsq(D, Y, rcond=None)
        rss1 = ((Y - D @ B) ** 2).sum(0)
        B0, *_ = np.linalg.lstsq(D[:, :3], Y, rcond=None)
        rss0 = ((Y - D[:, :3] @ B0) ** 2).sum(0)
        dfr = len(Y) - D.shape[1]
        p = stats.f.sf(((rss0 - rss1) / 2) / (rss1 / dfr), 2, dfr)
        pg = pd.Series(p, index=np.asarray(cols)[ok.to_numpy()])
        for c in ("local", "systemic", "other"):
            gg = pg.reindex(cls.index[cls == c]).dropna()
            desc.append(dict(item=f"share p_deviation < 0.01, {c} genes, {arm}", value=float((gg < 0.01).mean()),
                             n=len(gg)))
    pd.DataFrame(desc).to_csv(OUT / "descriptive.csv", index=False)
    pd.set_option("display.width", 250)
    print(tests.round(3).to_string(index=False), "\n", pd.DataFrame(desc).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
