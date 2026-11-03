"""A199 (registered in `paper/docs/registered_experiments.md`, commit 356c89c): does CYCLOPS magnitude measure clock strength in
samples taken in working hours?

GTEx v10 fast deaths with recorded time of death. Tissues: Breast - Mammary Tissue plus the 11 other tissues with most window
donors (A123b counts), Whole Blood excluded. Designs per tissue: W = deaths 08-18 h; Fs = as many fast deaths drawn at random
from the whole day. CYCLOPS 2.0 as in A123b (template hyperparameters, CYCLOPS seed list, no covariate), sample-order seeds
0-2; magnitude = mean over seeds of the within-fit z-scored `Magnitude`. Truth: timed clock strength, the projection of a
donor's centred clock genes on the 5-fold cross-fitted cosinor profile at its recorded time.
H199.1: median over tissues of rho_Fs - rho_W > 0 (rho = Spearman(magnitude, timed strength)); H199.2: median rho_Fs > 0.
One-sided Wilcoxon signed-rank over tissues, Holm over two.

    python a199.py fit [--jobs 3] [--threads 4]       # 72 CYCLOPS fits, cached in output/work/ (a rerun resumes)
    python a199.py score [--power-only]

Controlled access: recorded times and per-donor values stay in memory or under the git-ignored output/work/; output/*.csv hold
per-tissue aggregates only.
"""
import argparse
import concurrent.futures as cf
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path[:0] = [str(PAPER), str(PAPER / "results/figure_1/scripts/method")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr, wilcoxon  # noqa: E402

from _shared import load  # noqa: E402
from clocksync import paired as P  # noqa: E402
from common.benchmark import cyclops as C  # noqa: E402
from common.benchmark.data import in_window, log2_normalised  # noqa: E402

OUT = HERE.parents[1] / "output"
WORK = OUT / "work"
A123B = PAPER / "exploratory/A123_working_hours/output/a123b/runtime.csv"
BREAST, N_OTHER, SEED, SEEDS, K, LIMIT = "Breast - Mammary Tissue", 11, 20261103, (0, 1, 2), 2.49, 0.20
CLOCK = list(P.gtex().genes)
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:6.1f} min]", *a, flush=True)


def tissues():
    r = pd.read_csv(A123B)
    w = r[(r.design == "W") & (r.method == "cyclops")].set_index("tissue").n_samples.drop("Whole Blood", errors="ignore")
    return [BREAST] + [t for t in w.sort_values(ascending=False).index if t != BREAST][:N_OTHER]


def designs(tod):
    """{design: donor index}: W = 08-18 h deaths, Fs = as many random fast deaths from the whole day."""
    w = tod.index[in_window(tod)]
    rng = np.random.default_rng([SEED, len(tod), len(w)])
    return {"W": w, "Fs": pd.Index(rng.choice(tod.index.to_numpy(), len(w), replace=False))}


def slug(t):
    return t.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "-")


def fit_task(t, design, seed, threads):
    from clocksync.config import limit_threads

    limit_threads(threads)
    d = WORK / f"{design}_{slug(t)}_s{seed}"
    if (d / "fit_output.csv").exists():
        return t, design, seed, "cached", 0.0
    ts = time.time()
    counts, tod = load([t])[t]
    tod = tod.dropna()
    idx = designs(tod)[design]
    expr = 2 ** log2_normalised(counts.loc[tod.index]).loc[idx] - 1
    try:
        C.run_cyclops(expr.T, C.seeds(), d, seed=seed, threads=threads)
        status = "ok"
    except Exception as e:  # a failed fit drops the tissue (registered)
        status = f"failed: {type(e).__name__}: {e}"[:200]
    return t, design, seed, status, time.time() - ts


def fit(a):
    WORK.mkdir(parents=True, exist_ok=True)
    tasks = [(t, d, s) for t in tissues() for d in ("W", "Fs") for s in SEEDS]
    log(f"{len(tasks)} fits, {a.jobs} lanes x {a.threads} threads")
    rows = []
    with cf.ProcessPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(fit_task, *x, a.threads) for x in tasks]
        for f in cf.as_completed(futs):
            t, d, s, st, rt = f.result()
            rows.append(dict(tissue=t, design=d, seed=s, status=st, runtime_s=rt))
            log(f"{len(rows)}/{len(tasks)} {d} {t} s{s}: {st} ({rt / 60:.1f} min)")
            pd.DataFrame(rows).to_csv(OUT / "fits.csv", index=False)
    log("done")


# ---- scoring ------------------------------------------------------------------------------------------------------------

def timed_strength(L, tod):
    """5-fold cross-fitted cosinor of the clock genes; projection of each donor on the expected profile at its recorded time."""
    Y = L[[g for g in CLOCK if g in L.columns]].to_numpy(float)
    w = 2 * np.pi * tod.loc[L.index].to_numpy(float) / 24
    X = np.column_stack([np.ones(len(w)), np.cos(w), np.sin(w)])
    fold = np.random.default_rng(SEED).permutation(len(w)) % 5
    s = np.empty(len(w))
    for f in range(5):
        tr, te = fold != f, fold == f
        b = np.linalg.lstsq(X[tr], Y[tr], rcond=None)[0]
        p = X[te, 1:] @ b[1:]
        s[te] = ((Y[te] - b[0]) * p).sum(1) / (p ** 2).sum(1)
    return pd.Series(s, index=L.index)


def magnitude(t, design, idx):
    z = []
    for s in SEEDS:
        f = pd.read_csv(WORK / f"{design}_{slug(t)}_s{s}" / "fit_output.csv")
        m = pd.Series(f.Magnitude.astype(float).to_numpy(), index=f.ID.map({f"S{i}": x for i, x in enumerate(idx)}))
        ph = pd.Series(f.Phase.astype(float).to_numpy(), index=m.index)
        z.append(((m - m.mean()) / m.std(), ph))
    mags = pd.concat([x for x, _ in z], axis=1)
    return mags.mean(1), mags, z[0][1]


def circ_r(a, b):
    """|Jammalamadaka-SenGupta circular correlation|, best of the two directions (rotation-invariant by construction)."""
    def js(x, y):
        sx, sy = np.sin(x - np.angle(np.exp(1j * x).mean())), np.sin(y - np.angle(np.exp(1j * y).mean()))
        return (sx * sy).sum() / np.sqrt((sx ** 2).sum() * (sy ** 2).sum())
    return max(abs(js(a, b)), abs(js(-a, b)))


def r2_hour(m, h):
    w = 2 * np.pi * h / 24
    X = np.column_stack([np.ones(len(w)), np.cos(w), np.sin(w)])
    res = m - X @ np.linalg.lstsq(X, m, rcond=None)[0]
    return 1 - (res ** 2).sum() / ((m - m.mean()) ** 2).sum()


def score(a):
    ts = tissues()
    data = load(ts)
    power = []
    for t in ts:
        counts, tod = data[t]
        tod = tod.dropna()
        power.append(dict(tissue=t, n_all=len(tod), n_W=len(designs(tod)["W"])))
    pw = pd.DataFrame(power)
    pw = pw[pw.n_W > 3]
    sd = np.sqrt((2 / (pw.n_W - 3)).mean())
    mde = K * sd * 1.05 / np.sqrt(len(pw))
    pd.DataFrame([dict(tissues=len(pw), median_n_W=pw.n_W.median(), min_n_W=pw.n_W.min(), mde=mde, limit=LIMIT,
                       estimable=mde <= LIMIT)]).to_csv(OUT / "power.csv", index=False)
    log(f"power: {len(pw)} tissues, median n_W {pw.n_W.median():.0f}, MDE {mde:.3f} (limit {LIMIT})")
    if a.power_only:
        return

    fits = pd.read_csv(OUT / "fits.csv")
    bad = set(fits.loc[~fits.status.isin(["ok", "cached"]), "tissue"])
    rows = []
    for t in ts:
        if t in bad:
            rows.append(dict(tissue=t, dropped="fit failed")); continue
        counts, tod = data[t]
        tod = tod.dropna()
        L = log2_normalised(counts.loc[tod.index])
        strength = timed_strength(L, tod)
        expressed = L.loc[:, (L.mean() >= 3) & ~L.columns.isin(CLOCK)]
        depth = counts.loc[tod.index].sum(1)
        r = dict(tissue=t, n_W=len(designs(tod)["W"]))
        for d, idx in designs(tod).items():
            m, mags, ph = magnitude(t, d, idx)
            Z = expressed.loc[idx]
            Z = (Z - Z.mean()) / Z.std().replace(0, 1)
            pc1 = pd.Series(np.linalg.svd(Z.to_numpy(), full_matrices=False)[0][:, 0], index=idx)
            h = tod.loc[idx].to_numpy(float)
            cc = mags.corr(method="spearman").to_numpy()
            r.update({f"rho_{d}": spearmanr(m, strength.loc[idx])[0], f"r2_hour_{d}": r2_hour(m.to_numpy(), h),
                      f"rho_depth_{d}": spearmanr(m, depth.loc[idx])[0], f"abs_rho_pc1_{d}": abs(spearmanr(m, pc1)[0]),
                      f"seed_concordance_min_{d}": cc[np.triu_indices(3, 1)].min(),
                      f"phase_circ_r_{d}": circ_r(ph.to_numpy(), 2 * np.pi * h / 24)})
        rows.append(r)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "per_tissue.csv", index=False)
    ok = res.dropna(subset=["rho_W", "rho_Fs"])
    diff = ok.rho_Fs - ok.rho_W
    p1 = wilcoxon(diff, alternative="greater").pvalue
    p2 = wilcoxon(ok.rho_Fs, alternative="greater").pvalue
    hol = np.minimum(1, np.maximum.accumulate(np.array([2, 1]) * np.sort([p1, p2])))
    holm = dict(zip(np.argsort([p1, p2]), hol))
    tests = pd.DataFrame([
        dict(test="H199.1 median rho_Fs - rho_W", tissues=len(ok), statistic=float(diff.median()), p=p1, p_holm=holm[0]),
        dict(test="H199.2 median rho_Fs", tissues=len(ok), statistic=float(ok.rho_Fs.median()), p=p2, p_holm=holm[1])])
    tests["verdict"] = np.where(mde > LIMIT, "not estimable", np.where(tests.p_holm < .05, "pass", "fail"))
    tests.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(tests.round(4).to_string(index=False)); print(res.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["fit", "score"])
    ap.add_argument("--jobs", type=int, default=3); ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--power-only", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    fit(a) if a.step == "fit" else score(a)


if __name__ == "__main__":
    main()
