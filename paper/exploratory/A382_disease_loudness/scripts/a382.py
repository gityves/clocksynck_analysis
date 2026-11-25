"""A382 (registered 2026-11-25): do chronic diseases quiet the clocks of the body in absolute terms?

A373 statistic per condition: per cohort and tissue, A164 / A233 matched sets (case + up to 3 controls; obesity 1:1),
debiased clock-plane variance W per arm (A315), L_t = log(W_case / W_control), Q = median over tissues of both
cohorts; companions (200 expression-matched, 200 responsive sets, A373 seeds per cohort); label null: the case swapped
with a random member of its matched set.

Stages:
  power  matching, label null of Q (1000), companion Q; no clock Q for any condition.
  score  (written after the power note)

    RUN_CAP_GB=40 paper/common/run.sh paper/exploratory/A382_disease_loudness/scripts/a382.py --stage power
Outputs (`../output/`): power.csv, null.csv, companions.csv, matching.csv. Controlled access: labels, ages and matched
sets stay in memory; aggregates only.
"""
import argparse
import sys
import time
import warnings
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A373_age_loudness" / "scripts"))
sys.path.insert(0, str(EXP / "A233_chronic_conditions" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a373_age_loudness as L3  # noqa: E402
import compute_a233 as C233  # noqa: E402
from a373_age_loudness import ML, WF, Cohort  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
CONDITIONS = ("t2d", "obese", "htn", "ihd", "copd", "depression", "smoking")
N_SETS, N_PERM, SEED, MARGIN = 200, 1000, 382, 0.21


def design(c):
    if c == "t2d":
        return ML.Design("t2d", lambda q: (q.t2d == 1, q.t2d == 0), ML.match_ratio, min_cases=10)
    if c == "obese":
        return ML.Design("obese", lambda q: (q.bmi >= 30, (q.bmi >= 18.5) & (q.bmi < 25)),
                         partial(ML.match_ratio, ratio=1), min_cases=10)
    return C233.design_of(c)


def W_arm(z, v, idx):
    if len(idx) < L3.MIN_N:
        return np.nan
    zz = z[idx]
    return float(np.mean(np.abs(zz - zz.mean()) ** 2) - np.mean(v[idx]))


def Q_of(st, arms):
    L = []
    for k, (a, b) in arms.items():
        wa, wb = W_arm(*st[k], a), W_arm(*st[k], b)
        L.append(np.log(wa / wb) if (wa > 0 and wb > 0) else np.nan)
    return float(np.nanmedian(L)) if np.isfinite(L).any() else np.nan


class Sets:
    """Matched sets of one tissue as a padded row matrix (n_sets x 4, -1 = empty) for fast relabelling."""

    def __init__(self, sets):
        k = max(len(s) for s in sets)
        self.S = np.full((len(sets), k), -1)
        for i, s in enumerate(sets):
            self.S[i, :len(s)] = s
        self.n = (self.S >= 0).sum(1)

    def arms(self):
        return self.S[:, 0], self.S[:, 1:][self.S[:, 1:] >= 0]

    def relabel(self, rng):
        pick = (rng.random(len(self.n)) * self.n).astype(int)
        case = self.S[np.arange(len(self.n)), pick]
        cols = np.arange(self.S.shape[1])[None, :]
        ctl = self.S[(self.S >= 0) & (cols != pick[:, None])]
        return case, ctl


def setup():
    warnings.simplefilter("ignore")
    cohorts = ML.load_cohorts()
    plans, mrows = {}, []
    for c in CONDITIONS:
        d = design(c)
        for coh, (pb, q, _) in cohorts.items():
            plan, m, _ = ML.make_plan(d, q, pb, coh)
            mrows += m
            plans[(c, coh)] = plan
    tissues = {coh: sorted({t for (c, co), p in plans.items() if co == coh for t in p}) for coh in cohorts}
    co, sets = {}, {}
    for coh, h in ML.COH.items():
        pb, _ = gtex.load(tissues=tissues[coh], hardy=h, genes="all", depth="library_size", require_tod=(coh == "fast"))
        genes = [g for g in pb.columns if g not in ("donor", "compartment", "depth")]
        x = pb[genes].to_numpy(np.float32)
        lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
        lm = lm[~lm.index.duplicated()]
        del x
        sets[coh] = {kind: pseudo_waveforms(lm, N_SETS, np.random.default_rng([373, tag]), WF, kind=kind)
                     for kind, tag in (("matched", 1), ("responsive", 2))}
        keep = sorted({g for ws in sets[coh].values() for w in ws for g in w.genes} | set(WF.genes))
        co[coh] = Cohort(pb[["donor", "compartment", "depth"] + keep].reset_index(drop=True))
        del pb
        L3.say(coh, "genome loaded")
    designs = {}
    for c in CONDITIONS:
        arms = {}
        for coh in cohorts:
            for t, ss in plans[(c, coh)].items():
                pos = co[coh].pos.get(t, {})
                # a set enters with its case and at least one control present in the tissue's rows
                rows = [[pos[str(d)] for d in s if str(d) in pos] for s in ss if str(s[0]) in pos]
                rows = [r for r in rows if len(r) >= 2]
                if len(rows) >= 10:
                    arms[(coh, t)] = Sets(rows)
        designs[c] = arms
    return co, sets, designs, pd.DataFrame(mrows)


def power():
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    co, sets, designs, M = setup()
    M.to_csv(OUT / "matching.csv", index=False)
    keys = sorted({k for a in designs.values() for k in a})
    clock = {k: co[k[0]].zv(WF, k[1]) for k in keys}
    L3.say("clock states", f"({time.time() - t0:.0f}s)")
    comp_states = {(kind, j): {k: co[k[0]].zv(sets[k[0]][kind][j], k[1]) for k in keys}
                   for kind in ("matched", "responsive") for j in range(N_SETS)}
    L3.say("companion states", f"({time.time() - t0:.0f}s)")
    rows, nrows, crows = [], [], []
    for i, (c, arms) in enumerate(designs.items()):
        rng = np.random.default_rng([SEED, i])
        null = []
        for _ in range(N_PERM):
            null.append(Q_of(clock, {k: s.relabel(rng) for k, s in arms.items()}))
        null = np.array(null)
        base = {k: s.arms() for k, s in arms.items()}
        cq = {"matched": [], "responsive": []}
        for (kind, j), st in comp_states.items():
            v = Q_of(st, base)
            cq[kind].append(v)
            crows.append(dict(condition=c, kind=kind, set=j, Q=v))
        nrows += [dict(condition=c, Q=v) for v in null]
        sd = float(np.nanstd(null))
        rows.append(dict(condition=c, tissues=len(arms), tissues_vent=sum(k[0] == "vent" for k in arms),
                         cases_median=float(np.median([len(s.n) for s in arms.values()])),
                         label_null_sd=sd, label_null_q05=float(np.nanquantile(null, 0.05)),
                         expected_half_width=1.645 * sd, gate_H382_2=bool(1.645 * sd < MARGIN),
                         matched_median=float(np.nanmedian(cq["matched"])),
                         matched_q05=float(np.nanquantile(cq["matched"], 0.05)),
                         responsive_median=float(np.nanmedian(cq["responsive"])),
                         responsive_q05=float(np.nanquantile(cq["responsive"], 0.05))))
        L3.say(c, "done", f"({time.time() - t0:.0f}s)")
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "power.csv", index=False)
    pd.DataFrame(nrows).to_csv(OUT / "null.csv", index=False)
    pd.DataFrame(crows).to_csv(OUT / "companions.csv", index=False)
    pd.set_option("display.width", 220)
    print(out.round(3).to_string(index=False), f"\ndone {time.time() - t0:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power"], required=True)
    ap.parse_args()
    power()


if __name__ == "__main__":
    main()
