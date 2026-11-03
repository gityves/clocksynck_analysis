"""A201 (registered in `paper/docs/registered_experiments.md`, commit 745b167): after a 6 h light advance, do aortic cell types keep
their clock phase relations while losing alignment to the light cycle? (Auerbach et al., bioRxiv 2026, Zenodo 22847203.)

Per library x cell type, pseudobulk counts of the 14 clock genes (mouse orthologs of the GTEx set), log(CPM + 1). Per condition
and cell type, per gene cosinor on ZT with a sex term; z_g = a_g - i b_g. Delta_c = arg sum_g z_g^mis conj(z_g^ali) (h).
H201.0: mean of Delta_SMC and Delta_fib differs from 0 (90 % interval excludes 0). H201.1: D = Delta_SMC - Delta_fib within
+/- 1.5 h (90 % interval inside). SD = larger of a Poisson-weight cell bootstrap (1,000) and a library jackknife. Power: 1.645 SD
<= 0.75 h. Secondary: amplitude ratios, EC and macrophage, per sex, Bmal1 KO amplitude collapse, clocksync gate per condition.

    python score_a201.py [--nboot 1000] [--power-only]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

H5 = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/aorta_jetlag_auerbach2026/adata_qc_filtered_compressed.h5ad")
OUT = HERE.parents[1] / "output"
GENES = [g.capitalize() for g in cs.gtex().genes]  # ARNTL -> Arntl, ...
MAIN, EXTRA = ["SMC", "Fibroblast"], ["EC", "Macrophage"]
BOUND, LIMIT, MIN_CELLS, SEED = 1.5, 0.75, 30, 20261103
W = 2 * np.pi / 24


def wrap(h):
    return (np.asarray(h) + 12) % 24 - 12


def load():
    a = ad.read_h5ad(H5, backed="r")
    genes = [g for g in GENES if g in a.var_names]
    o = a.obs[["sample_id", "sex", "misaligned", "zt", "bmal1_ko", "cell_type", "total_counts"]].copy()
    o["cond"] = np.where(o.bmal1_ko, "ko", np.where(o.misaligned, "mis", "ali"))
    X = a[:, genes].to_memory().X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return o.reset_index(drop=True), X.astype(float), genes


def pseudobulk(o, X, w=None):
    """Per (library, cell type): clock-gene counts, total counts, cells; `w` = per-cell weights (bootstrap)."""
    w = np.ones(len(o)) if w is None else w
    key = o.sample_id.astype(str) + "|" + o.cell_type.astype(str)
    codes, uniq = pd.factorize(key)
    k = len(uniq)
    C = np.column_stack([np.bincount(codes, X[:, j] * w, k) for j in range(X.shape[1])])
    T, N = np.bincount(codes, o.total_counts.to_numpy() * w, k), np.bincount(codes, w, k)
    lib = o.groupby(key.to_numpy(), observed=True)[["sample_id", "cell_type", "sex", "zt", "cond"]].first().loc[uniq]
    return lib.assign(total=T, cells=N), C


def fit(lib, C, cond, ct, keep=None):
    """Complex amplitudes z_g (genes) of a condition x cell type: cosinor on ZT with a sex term."""
    m = (lib.cond == cond).to_numpy() & (lib.cell_type == ct).to_numpy()
    if keep is not None:
        m &= keep
    Y = np.log(C[m] / lib.total.to_numpy()[m, None] * 1e6 + 1)
    w = W * lib.zt.to_numpy(float)[m]
    D = np.column_stack([np.ones(m.sum()), (lib.sex.to_numpy()[m] == "female").astype(float), np.cos(w), np.sin(w)])
    b = np.linalg.lstsq(D, Y, rcond=None)[0]
    return b[2] - 1j * b[3]


def shift(lib, C, ct, keep=None):
    zm, za = fit(lib, C, "mis", ct, keep), fit(lib, C, "ali", ct, keep)
    return float(np.angle(np.sum(zm * np.conj(za))) / W), float(np.abs(zm).sum() / np.abs(za).sum())


def stats(lib, C, keep=None):
    s = {ct: shift(lib, C, ct, keep) for ct in MAIN}
    d_smc, d_fib = s["SMC"][0], s["Fibroblast"][0]
    mean = float(np.angle(np.exp(1j * W * d_smc) + np.exp(1j * W * d_fib)) / W)
    return dict(delta_smc=d_smc, delta_fib=d_fib, mean_shift=mean, D=float(wrap(d_smc - d_fib)),
                amp_ratio_smc=s["SMC"][1], amp_ratio_fib=s["Fibroblast"][1])


def uncertainty(o, X, lib, C, est, nboot):
    rng = np.random.default_rng(SEED)
    boot = []
    for _ in range(nboot):
        lb, Cb = pseudobulk(o, X, rng.poisson(1.0, len(o)).astype(float))
        boot.append(stats(lb, Cb))
    boot = pd.DataFrame(boot)
    jack = []
    libs = lib.loc[lib.cond.isin(["ali", "mis"]), "sample_id"].unique()
    for s in libs:
        jack.append(stats(lib, C, (lib.sample_id != s).to_numpy()))
    jack = pd.DataFrame(jack)
    n = len(jack)
    sd = {}
    for k in ("D", "mean_shift"):
        bs = np.std(wrap(boot[k] - est[k]), ddof=1)
        js = np.sqrt((n - 1) / n * np.sum(wrap(jack[k] - jack[k].mean()) ** 2))
        sd[k] = dict(sd_boot=float(bs), sd_jack=float(js), sd=float(max(bs, js)))
    return sd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nboot", type=int, default=1000); ap.add_argument("--power-only", action="store_true")
    a = ap.parse_args()
    limit_threads(8)
    OUT.mkdir(parents=True, exist_ok=True)
    o, X, genes = load()
    keep = o.cell_type.isin(MAIN + EXTRA).to_numpy()
    o, X = o[keep].reset_index(drop=True), X[keep]
    lib, C = pseudobulk(o, X)
    print(f"{len(genes)} clock genes: {genes}", flush=True)
    est = stats(lib, C)
    sd = uncertainty(o[o.cell_type.isin(MAIN)].reset_index(drop=True), X[o.cell_type.isin(MAIN).to_numpy()],
                     lib, C, est, a.nboot)
    estimable = 1.645 * sd["D"]["sd"] <= LIMIT
    pd.DataFrame([dict(genes=len(genes), **{f"{k}_{q}": v for k, d in sd.items() for q, v in d.items()},
                       limit=LIMIT, estimable=estimable)]).to_csv(OUT / "power.csv", index=False)
    print(pd.read_csv(OUT / "power.csv").round(3).to_string(index=False), flush=True)
    if a.power_only:
        return

    lo0, hi0 = est["mean_shift"] - 1.645 * sd["mean_shift"]["sd"], est["mean_shift"] + 1.645 * sd["mean_shift"]["sd"]
    lo1, hi1 = est["D"] - 1.645 * sd["D"]["sd"], est["D"] + 1.645 * sd["D"]["sd"]
    g0 = not (lo0 <= 0 <= hi0)
    tests = pd.DataFrame([
        dict(test="H201.0 misalignment acted (mean shift != 0)", estimate=est["mean_shift"], lo90=lo0, hi90=hi0,
             verdict="pass" if g0 else "fail"),
        dict(test="H201.1 SMC - fibroblast shift within +/-1.5 h", estimate=est["D"], lo90=lo1, hi90=hi1,
             verdict="not estimable" if not estimable else ("uninformative (H201.0 failed)" if not g0 else
                                                            ("pass" if (lo1 > -BOUND and hi1 < BOUND) else "fail")))])
    tests.to_csv(OUT / "tests.csv", index=False)

    # ---- secondary ----
    sec = [dict(analysis="estimates", **est)]
    for ct in EXTRA:
        ok = lib[(lib.cell_type == ct) & lib.cond.isin(["ali", "mis"])].cells.min() >= MIN_CELLS
        d, r = shift(lib, C, ct)
        sec.append(dict(analysis=f"{ct} shift (min cells ok: {ok})", delta=d, amp_ratio=r))
    for sx in ("female", "male"):
        s = stats(lib, C, (lib.sex == sx).to_numpy())
        sec.append(dict(analysis=f"{sx} only", **s))
    for ct in MAIN:  # Bmal1 KO: amplitude collapse relative to aligned controls
        sec.append(dict(analysis=f"Bmal1 KO / aligned amplitude, {ct}",
                        amp_ratio=float(np.abs(fit(lib, C, "ko", ct)).sum() / np.abs(fit(lib, C, "ali", ct)).sum())))
    # clocksync untimed: libraries as donors, cell types as compartments, per condition
    g = []
    for cond in ("ali", "mis", "ko"):
        m = (lib.cond == cond).to_numpy() & lib.cell_type.isin(MAIN).to_numpy()
        cnt = pd.DataFrame(C[m], columns=[x.upper() if x.upper() in cs.gtex().genes else x for x in genes])
        key = lib[m][["sample_id", "cell_type"]].astype(str).reset_index(drop=True)
        pb = cs.from_bulk(cnt.round(), key, donor="sample_id", compartment="cell_type",
                          depth=pd.Series(lib.total.to_numpy()[m]))
        try:
            gt = cs.Dataset(pb, anchors=MAIN).gate(n_perm=1000, seed=SEED)
            g.append(dict(condition=cond, donors=gt.n_donors, p_pairing=gt.p_pairing, p_acrophase=gt.p_acrophase,
                          clock_share=gt.clock_share))
        except Exception as e:
            g.append(dict(condition=cond, error=f"{type(e).__name__}: {e}"[:150]))
    pd.DataFrame(sec).to_csv(OUT / "secondary.csv", index=False)
    pd.DataFrame(g).to_csv(OUT / "gate.csv", index=False)
    pd.set_option("display.width", 250)
    print(tests.round(3).to_string(index=False)); print(pd.DataFrame(sec).round(3).to_string(index=False))
    print(pd.DataFrame(g).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
