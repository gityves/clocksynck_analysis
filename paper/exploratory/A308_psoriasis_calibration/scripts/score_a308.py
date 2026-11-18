"""A308 (registered in `paper/docs/registered_experiments.md`): is the psoriatic lesion delay (A183 H183.1) beyond the
lesion-wide rotation of random expression-matched genes, once the companions' own estimation noise is removed?

A183's H183.1 statistic (pooled `clocksync.paired` offset over the used batches, clock-fit weights) for the clock and
A183's 200 matched and 200 responsive sets (same seeds). Draws: per-person label swaps (`paired.swap_offset_null` form,
states centred within sample type), the same swaps for every set; draw 0 = no swap, reproduces A183.
  V308  share of draws with the clock below the matched q05 of the draw (A285 form; responsive reported).
  H308  q = Phi((o_c - m) / sqrt(tau^2 + s_c^2)), m and tau^2 = var(o_j) - mean(s_j^2) over the matched sets, s_j the
        SD of set j's placebo offsets; passes if q < 0.05. 90 % interval over resampled sets (reported).

    python score_a308.py [--draws 2000]
Outputs (`../output/`): sets.csv (per set: observed offset, noise SD), tests.csv.
"""
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
A183_DIR = PAPER / "exploratory/A183_skin_gse154200"
spec = importlib.util.spec_from_file_location("a308_a183", A183_DIR / "scripts/score_a183.py")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 308
N_SETS = 200


def complex_pairs(case, ref, wf):
    """Centred complex states (lesional a, non-lesional b) and their measurement variances, one projector."""
    ba, va, bb, vb = P.cohort_states(case, ref, True, wf)
    ba, bb = ba - ba.mean(0), bb - bb.mean(0)
    return ba[:, 0] + 1j * ba[:, 1], va, bb[:, 0] + 1j * bb[:, 1], vb


def ccorr_draws(a, va, b, vb, F):
    """`stats.ccorr` of (lesional, non-lesional) under each row of swaps F (draws x persons), vectorised."""
    A_ = np.where(F, b, a)
    B_ = np.where(F, a, b)
    VA, VB = np.where(F, vb, va), np.where(F, va, vb)
    A_ = A_ - A_.mean(1, keepdims=True)
    B_ = B_ - B_.mean(1, keepdims=True)
    k = a.size - 1
    ea = np.sum(abs(A_) ** 2, 1) / k - VA.mean(1)
    eb = np.sum(abs(B_) ** 2, 1) / k - VB.mean(1)
    return np.sum(A_ * np.conj(B_), 1) / k / np.sqrt(np.maximum(ea, 1e-12) * np.maximum(eb, 1e-12))


def offsets(pairs, w, F):
    """Pooled offset [h] per draw: sum over batches of w * ccorr."""
    z = sum(wk * ccorr_draws(*pairs[k], F[k]) for k, wk in w.items())
    return np.angle(z) * P.HOURS


def wsd(h):
    return float(np.sqrt(np.mean(A.wrap(h - np.angle(np.mean(np.exp(1j * h * np.pi / 12))) * 12 / np.pi) ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=2000)
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    s = A.load()
    means = A.gene_means(s)
    wf = P.gtex()
    D = A.Design(s)
    coh = pd.read_csv(A183_DIR / "output/h183_1_paired_cohorts.csv", dtype={"cohort": str})
    used = str(pd.read_csv(A183_DIR / "output/h183_1_paired_pooled.csv", index_col=0).loc["cohorts", "value"]).split(";")
    w = dict(zip(used, 1 / coh.set_index("cohort").loc[used, "se"].to_numpy() ** 2))
    ps = s[s.diagnosis == "psoriasis"]
    data = {}
    for k in used:
        case = ps[(ps.cohort == k) & (ps.condition == "lesional")].set_index("unit")
        ref = ps[(ps.cohort == k) & (ps.condition == "nonlesional")].set_index("unit")
        u = sorted(set(case.index) & set(ref.index))
        data[k] = (case.loc[u], ref.loc[u])
    sets = [("clock", 0, {k: wf for k in used})]
    for kind, sub in (("matched", 7), ("responsive", 8)):
        crng = np.random.default_rng([A.SEED, sub])
        per = {k: P.pseudo_waveforms(means[k], N_SETS, crng, wf, kind) for k in D.keys}
        sets += [(kind, j, {k: per[k][j] for k in used}) for j in range(N_SETS)]
    rng = np.random.default_rng(SEED)
    F = {k: np.vstack([np.zeros(len(data[k][0]), bool), rng.random((a.draws, len(data[k][0]))) < 0.5]) for k in used}
    O = np.empty((len(sets), a.draws + 1))
    for i, (_, _, wfs) in enumerate(sets):
        O[i] = offsets({k: complex_pairs(*data[k], wfs[k]) for k in used}, w, F)
    kinds = np.array([k for k, _, _ in sets])
    tab = pd.DataFrame(dict(kind=kinds, set=[j for _, j, _ in sets], offset_h=O[:, 0], noise_sd_h=[wsd(o[1:]) for o in O]))
    tab.to_csv(OUT / "sets.csv", index=False)
    ref = pd.read_csv(A183_DIR / "output/tests.csv").iloc[0]
    clock = O[0]
    rows = []
    for kind in ("matched", "responsive"):
        M = O[kinds == kind]
        q05 = np.quantile(M, 0.05, axis=0)
        repro = bool(np.isclose(clock[0], ref.statistic_h, atol=1e-3) and np.isclose(q05[0], ref[f"{kind}_q05"], atol=1e-3))
        o, sd = M[:, 0], tab.loc[tab.kind == kind, "noise_sd_h"].to_numpy()
        sc = tab.noise_sd_h.iloc[0]

        def q_of(idx):
            tau2 = max(0.0, np.var(o[idx], ddof=1) - np.mean(sd[idx] ** 2))
            return norm.cdf((clock[0] - o[idx].mean()) / np.sqrt(tau2 + sc ** 2)), tau2

        q, tau2 = q_of(np.arange(len(o)))
        brng = np.random.default_rng([SEED, 1])
        qb = [q_of(brng.integers(0, len(o), len(o)))[0] for _ in range(2000)]
        rows.append(dict(family=kind, repro=repro, clock_offset=clock[0], clock_noise_sd=sc, q05_obs=q05[0],
                         pct_obs=float(np.mean(o <= clock[0])), draws=a.draws,
                         v308_rejection=float(np.mean(clock[1:] < q05[1:])), sets_mean=float(o.mean()),
                         sets_var=float(np.var(o, ddof=1)), noise_var_mean=float(np.mean(sd ** 2)), tau2=tau2,
                         tau=float(np.sqrt(tau2)), noise_share=float(np.mean(sd ** 2) / np.var(o, ddof=1)), q=float(q),
                         q_lo90=float(np.quantile(qb, 0.05)), q_hi90=float(np.quantile(qb, 0.95)),
                         passes=bool(q < 0.05) if kind == "matched" else np.nan))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "tests.csv", index=False)
    print(T.T.to_string())


if __name__ == "__main__":
    main()
