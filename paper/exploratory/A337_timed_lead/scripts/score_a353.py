"""A353 scoring (registered in `paper/docs/registered_experiments.md`): is the pre-mRNA lead of human muscle recovered from
daytime biopsies alone? Perrin 2018 vastus lateralis (adapter-trimmed re-map), direct introns, A337 units.

Subsets of clock times (fixed in the registration): full (00-20 h, six times), D3 (08, 12, 16 h; primary), D4 (08, 12,
16, 20 h), D2 (08, 12 h), N3 (20, 00, 04 h; descriptive). In each subset only subjects with every time of the subset
are kept. Per subset and waveform (14 genes; core loop descriptive):
- untimed: projector from the subset's spliced counts, states residualised on the subject within the subset
  (`score_a337.centred_states`); lead delta_U = arg sum z_u conj(z_s), gate kappa with its within-subject permutation p;
- subset timed lead (descriptive): arg(beta_u conj(beta_s)), beta = sum z exp(-i 2 pi t / 24);
- per-sample phase (descriptive): angle of each sample's spliced state against its clock time after one common
  rotation; concentration R and median absolute error (h).
Truth: the full-cycle timed lead (all subjects, full projector). Paired 90 % intervals from 2,000 subject-cluster
bootstraps: subjects of the subset drawn with replacement; the subset lead and the truth are recomputed on the same
subjects (the truth on their full-cycle rows). The full row also carries A337's within-subject interval.
H353.1 (D3, 14 genes): gate p < 0.05 and the paired interval of delta_U - truth inside +-1.5 h (pass); interval beyond
the bound on both sides: undecidable; otherwise fail. H353.2: the same on D4 and D2.
`--check` replaces the unspliced counts with thinned spliced counts and shuffles times within subject (no real
statistic; expected lead near 0) to test the code.

    python score_a353.py [--check]
Outputs (`../output/a353[_check]/`): subsets.csv, verdicts.csv.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_a337 as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector  # noqa: E402

SEED, N_PERM, N_BOOT, BOUND = 353, 10000, 2000, 1.5
SUBSETS = {"full": (0, 4, 8, 12, 16, 20), "D3": (8, 12, 16), "D4": (8, 12, 16, 20), "D2": (8, 12), "N3": (20, 0, 4)}
H = 24 / (2 * np.pi)


def states(sp, un, obs, wf):
    L = projector(sp, wf)
    return A.centred_states(sp, obs, wf, L), A.centred_states(un, obs, wf, L)


def timed(zu, zs, t):
    e = np.exp(-1j * A.W * t)
    return float(np.angle(np.sum(zu * e) * np.conj(np.sum(zs * e))) * H)


def visit_phase(zs, t):
    d = np.angle(zs) - A.W * t
    rot = np.angle(np.mean(np.exp(1j * d)))
    return float(np.abs(np.mean(np.exp(1j * d)))), float(np.median(np.abs(A.wrap((d - rot) * H))))


def verdict(gate_p, lo, hi):
    if gate_p < .05 and -BOUND < lo and hi < BOUND:
        return "pass"
    return "undecidable" if lo <= -BOUND and hi >= BOUND else "fail"


def score(unit, wname, rng):
    wf = A.WAVEFORMS[wname]
    obs, sp, un = unit["obs"], unit["spliced"], unit["direct"]
    t, grp = obs.t.to_numpy(float), obs.group.to_numpy()
    zs_f, zu_f = states(sp, un, obs, wf)
    truth = timed(zu_f, zs_f, t)
    rows = []
    for name, times in SUBSETS.items():
        # full keeps every library (reproduces A337 Perrin); a subset keeps subjects with all of its times
        have = obs.groupby("group").t.apply(lambda x: name == "full" or set(times) <= set(x))
        subj = have.index[have].to_numpy()
        m = obs.group.isin(subj).to_numpy() & np.isin(t, times)
        o = obs[m]
        zs, zu = states(sp[m], un[m], o, wf)
        ts, gs = o.t.to_numpy(float), o.group.to_numpy()
        dU = A.lead(zu, zs)
        R, err = visit_phase(zs, ts)
        row = dict(waveform=wname, subset=name, times=" ".join(f"{x:02d}" for x in times), subjects=len(subj),
                   samples=int(m.sum()), kappa=A.kappa(zu, zs), gate_p=A.gate_p(zu, zs, gs, rng, N_PERM),
                   delta_U=dU, delta_T_subset=timed(zu, zs, ts), truth=truth, diff=float(A.wrap(dU - truth)),
                   phase_R=R, phase_median_err_h=err)
        idx_s = {g: np.flatnonzero(gs == g) for g in subj}
        idx_f = {g: np.flatnonzero(grp == g) for g in subj}
        bU, bT = np.empty(N_BOOT), np.empty(N_BOOT)
        for b in range(N_BOOT):
            draw = rng.choice(subj, len(subj))
            i = np.concatenate([idx_s[g] for g in draw])
            j = np.concatenate([idx_f[g] for g in draw])
            bU[b], bT[b] = A.lead(zu[i], zs[i]), timed(zu_f[j], zs_f[j], t[j])
        row["dU_lo"], row["dU_hi"] = A.interval(dU, bU)
        row["truth_lo"], row["truth_hi"] = A.interval(truth, bT)
        row["diff_lo"], row["diff_hi"] = A.interval(row["diff"], A.wrap(bU - bT))
        if name == "full":  # A337's within-subject resampling, for comparison of interval widths
            B = A.boot(gs, rng, N_BOOT)
            row["dU_lo_within"], row["dU_hi_within"] = A.interval(dU, [A.lead(zu[b], zs[b]) for b in B])
        rows.append(row)
        print(f"  {wname} {name}: {len(subj)} subjects, {int(m.sum())} samples", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    global N_PERM, N_BOOT
    limit_threads(4)
    t0 = time.time()
    unit = A.load("perrin")[0]
    if a.check:  # no real statistic: unspliced := thinned spliced, times shuffled within subject
        N_PERM, N_BOOT = 200, 100
        rng = np.random.default_rng(0)
        unit["direct"] = unit["spliced"].apply(lambda c: rng.binomial(c.astype(np.int64), 0.1))
        o = unit["obs"]
        unit["obs"] = o.assign(t=o.groupby("group").t.transform(lambda x: rng.permutation(x.to_numpy())))
    print(len(unit["obs"]), "samples;", unit["obs"].group.nunique(), "subjects;",
          dict(unit["obs"].t.value_counts().sort_index()), flush=True)
    rows = []
    for k, w in enumerate(A.WAVEFORMS):
        rows += score(unit, w, np.random.default_rng([SEED, k]))
    t = pd.DataFrame(rows)
    out = A.OUT / ("a353_check" if a.check else "a353")
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "subsets.csv", index=False)
    v = t[(t.waveform == "all") & t.subset.isin(["D3", "D4", "D2"])].assign(
        test=lambda d: np.where(d.subset == "D3", "H353.1", "H353.2"),
        verdict=lambda d: [verdict(g, lo, hi) for g, lo, hi in zip(d.gate_p, d.diff_lo, d.diff_hi)])
    v = v[["test", "subset", "times", "subjects", "samples", "gate_p", "delta_U", "truth", "diff", "diff_lo", "diff_hi",
           "verdict"]]
    v.to_csv(out / "verdicts.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.round(3).to_string(index=False), v.round(3).to_string(index=False), sep="\n")
    print(f"runtime {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
