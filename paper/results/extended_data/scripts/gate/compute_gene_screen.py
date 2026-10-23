"""A28 gene screen: clock-target genes that share the clock's phase across tissues (selection half of GTEx fast deaths).

Pre-registered as A28 (2026-10-15; docs/PLAN_internal_desync.md (archive); the screen was run on the selection half before the amendment was written and
fixes the gene sets W18 / Wcons / Wgen). Hypothesis id: none (screen). Split: GTEx fast deaths (Hardy 1-2) by subject, 50 / 50, seed 0. Per host
tissue (the 15 tissues of A28): cosinor of log(CPM + 0.5) on the time of death with the D2 covariates (sex, age bracket, Hardy, ischaemic time and its
square) -> complex coefficient c_gt; the tissue's complex clock gain g_t from the 14 clock genes against the global waveform (`gtex()`), c'_gt =
c_gt / g_t. Per gene, over the tissues where it is expressed (median CPM >= 5): shared c' (inverse-variance mean), z, heterogeneity I^2 and Cochran Q p,
acrophase SD over the tissues where it is significant (|c'| / se > 2), amplitude-ratio CV. Ventilator support: tissues where the gene follows the host
(p_int < 0.001, p_ext > 0.1, phase within 3 h of the shared phase) in the host frame and in the CHIRAL-DIP frame (tables diag5_genes and diag5_dip_genes of
`gtex_cues`; run it first). Thresholds from the 14 clock genes: tissues >= 10, z >= 6, amp >= 0.1, I^2 <= 0.5, phase SD <= 2 h (strict set), support in >= 2
tissues for both frames (W18; without the support: Wcons); the generous set (tissues >= 10, z >= 4, amp >= 0.05, I^2 <= 0.5, phase SD <= 3 h with support) is reported as a
count. Null: none (descriptive screen). Aggregates over genes only (GTEx controlled access): no per-donor rows.

    python compute_gene_screen.py [--subsample N] [--seed S]
Outputs (paper/results/extended_data/output/gate/): gene_screen_genes.csv (per gene), gene_screen_gains.csv (per tissue), gene_screen_sets.csv (the sets)
Legacy: legacy/desync/validation/extended_clock/scripts/screen.py -> output/screen_genes.csv, selection_subjects.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.gtex import FAST, LABELLED  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "gate"
S07 = "gtex_cues"
W = 2 * np.pi / 24
TISSUES = ("Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)",
           "Esophagus - Mucosa", "Nerve - Tibial", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
           "Adipose - Subcutaneous", "Lung", "Pituitary", "Brain - Cerebellum", "Brain - Cortex", "Whole Blood")


def covariates(ob) -> np.ndarray:
    """local: not in clocksync. D2 covariates: sex, age bracket, Hardy dummies, standardised ischaemic time and its square."""
    Z = pd.get_dummies(pd.DataFrame(dict(sex=ob.sex.astype(str), age=ob.age_bracket.astype(str), hardy=ob.DTHHRDY.astype(str))),
                       drop_first=True).to_numpy(float)
    isch = ob.ischemic_time.astype(float).fillna(ob.ischemic_time.astype(float).median()).to_numpy()
    isch = (isch - isch.mean()) / (isch.std() + 1e-9)
    return np.c_[Z, isch, isch ** 2]


def summarise(x: pd.DataFrame) -> pd.Series:
    w = 1 / x.se.to_numpy() ** 2
    cc = x.c.to_numpy()
    m = (w * cc).sum() / w.sum()
    se_m = np.sqrt(1 / w.sum())
    Q = (w * np.abs(cc - m) ** 2).sum() / 2
    k, dfq = len(x), len(x) - 1
    I2 = max(0.0, (Q - dfq) / Q) if Q > 0 else 0.0
    sig = np.abs(cc) / x.se.to_numpy() > 2
    ph = np.angle(cc[sig])
    spread = np.sqrt(-2 * np.log(abs(np.exp(1j * ph).mean()))) / W if sig.sum() >= 2 else np.nan
    ratio = np.abs(cc[sig]) / abs(m) if sig.sum() else np.array([np.nan])
    return pd.Series(dict(tissues=k, amp=abs(m), psi_h=np.angle(m) / W % 24, z=abs(m) / se_m, I2=I2,
                          p_het=stats.chi2.sf(2 * Q, 2 * dfq), n_sig=int(sig.sum()), phase_sd_h=spread,
                          amp_ratio_cv=np.std(ratio) / np.mean(ratio) if sig.sum() >= 2 else np.nan))


def support(path, S):
    v = pd.read_csv(path)
    v = v[(v.p_int < 1e-3) & (v.p_ext > 0.1)].join(S[["psi_h"]], on="gene")
    v = v[np.abs(np.angle(np.exp(1j * W * (v.int_h - v.psi_h)))) / W <= 3]
    return v.groupby("gene").tissue.nunique()


sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_out  # noqa: E402
def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    sub = "subsample" if run["quick"] else ""
    a = ad.read_h5ad(LABELLED, backed="r")
    o = a.obs
    subj = np.array(sorted(o[o.DTHHRDY.isin(list(FAST))].subject.astype(str).unique()))
    sel = set(np.random.default_rng(0).choice(subj, len(subj) // 2, replace=False))
    wf = cs.gtex()
    Xc = wf.amp * np.exp(1j * wf.psi)
    var = np.asarray(a.var_names)
    rows, gains = [], []
    for tis in TISSUES:
        keep = np.flatnonzero((o.DTHHRDY.isin(list(FAST)) & (o.tissue.astype(str) == tis) & o.tod.notna()
                               & o.subject.astype(str).isin(sel)).to_numpy())
        if args.subsample and len(keep) > args.subsample:
            keep = np.sort(np.random.default_rng(args.seed).permutation(keep)[: args.subsample])
        ob = o.iloc[keep]
        X = a.X[keep]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        cpm = X / ob.library_size.to_numpy(float)[:, None] * 1e6
        ex = np.median(cpm, 0) >= 5
        Y = np.log(cpm[:, ex] + 0.5)
        t = ob.tod.to_numpy(float) % 24
        D = np.c_[np.ones(len(Y)), covariates(ob), np.cos(W * t), np.sin(W * t)]
        B = np.linalg.lstsq(D, Y, rcond=None)[0]
        res = Y - D @ B
        s2 = (res ** 2).sum(0) / (len(Y) - D.shape[1])
        XtXi = np.linalg.pinv(D.T @ D)
        se = np.sqrt(s2 * (XtXi[-2, -2] + XtXi[-1, -1]) / 2)
        c = B[-2] + 1j * B[-1]
        gnames = var[ex]
        names = list(gnames)
        ci = [names.index(g) for g in wf.genes if g in names]
        xc = np.array([Xc[list(wf.genes).index(names[i])] for i in ci])
        g_t = (np.conj(xc) @ c[ci]) / (np.conj(xc) @ xc)
        gains.append(dict(tissue=tis, samples=len(Y), gain_abs=abs(g_t), shift_h=np.angle(g_t) / W))
        rows.append(pd.DataFrame(dict(tissue=tis, gene=gnames, c=c / g_t, se=se / abs(g_t), n=len(Y), mean=Y.mean(0))))
    a.file.close()
    df = pd.concat(rows)
    S = df.groupby("gene").apply(summarise, include_groups=False)
    f1, f2 = group_out(S07) / sub / "diag5_genes.csv.gz", group_out(S07) / sub / "diag5_dip_genes.csv.gz"
    S["vent_host"], S["vent_dip"] = support(f1, S), support(f2, S)
    S = S.fillna({"vent_host": 0, "vent_dip": 0})
    S["clock"] = S.index.isin(list(wf.genes))
    S.to_csv(out / "gene_screen_genes.csv")
    pd.DataFrame(gains).to_csv(out / "gene_screen_gains.csv", index=False)
    base = (S.tissues >= 10) & (S.I2 <= 0.5) & ~S.clock
    sup = (S.vent_host >= 2) & (S.vent_dip >= 2)
    strict = base & (S.z >= 6) & (S.amp >= 0.1) & (S.phase_sd_h <= 2)
    sets = {"strict_with_support (W18)": strict & sup, "strict_without_support (Wcons)": strict,
            "generous (before the D6 exclusion)": base & (S.z >= 4) & (S.amp >= 0.05) & sup & (S.phase_sd_h <= 3)}
    rep = pd.DataFrame([dict(set=k, genes=int(v.sum()), members=";".join(S[v].sort_values("z", ascending=False).index[:40])) for k, v in sets.items()])
    rep.to_csv(out / "gene_screen_sets.csv", index=False)
    print(f"selection subjects {len(sel)}")
    print(S[S.clock].round(2).to_string())
    print(rep.to_string(index=False))


if __name__ == "__main__":
    main()
