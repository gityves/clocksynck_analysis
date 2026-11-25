"""A382 (registered 2026-11-25): do chronic diseases quiet the clocks of the body in absolute terms?

A373 statistic per condition: per cohort and tissue, A164 / A233 matched sets (case + up to 3 controls; obesity 1:1),
debiased clock-plane variance W per arm (A315), L_t = log(W_case / W_control), Q = median over tissues of both
cohorts; companions (200 expression-matched, 200 responsive sets, A373 seeds per cohort); label null: the case swapped
with a random member of its matched set.

Stages:
  power  matching, label null of Q (1000), companion Q; no clock Q for any condition.
  score  (written after the power note, before any case-against-control loudness) per condition: Q (14 genes, primary),
         core loop, pooled (unmatched) arms, per cohort; 90 % matched-set bootstrap (2000, sets resampled within each
         tissue cell); label p from the power-stage null; Q* against the power-stage companions; H382.1 (Holm over the
         seven) and H382.2; glucocorticoid-target check: per cell, case-minus-control mean of the per-sample GC score
         (mean z-scored log CPM of `controls.GC_TARGETS` present), median over cells.

    RUN_CAP_GB=40 paper/common/run.sh paper/exploratory/A382_disease_loudness/scripts/a382.py --stage power | score
Outputs (`../output/`): power.csv, null.csv, companions.csv, matching.csv | tests.csv. Controlled access: labels, ages and matched
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
from clocksync._core.controls import GC_TARGETS  # noqa: E402
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
        keep = sorted({g for ws in sets[coh].values() for w in ws for g in w.genes} | set(WF.genes)
                      | {g for g in GC_TARGETS if g in pb.columns})
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


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]))
    out = np.empty_like(adj)
    out[o] = adj
    return out


def gc_score(cohort, t):
    """Per-sample glucocorticoid-target score of a tissue: mean z-scored log CPM of the GC targets present."""
    x = cohort.pb.loc[cohort.rows[t]]
    g = [c for c in GC_TARGETS if c in x.columns]
    y = np.log(x[g].to_numpy(float) / x.depth.to_numpy(float)[:, None] * 1e6 + 1)
    return ((y - y.mean(0)) / np.where(y.std(0) > 0, y.std(0), 1)).mean(1)


def pooled_arms(cohorts, c, co):
    """Every donor of the case and control pools with data in the tissue (unmatched), per eligible cell."""
    d = design(c)
    out = {}
    for coh, (pb, q, _) in cohorts.items():
        for t, pos in co[coh].pos.items():
            ids = [x for x in pos if x in q.index]
            cas, con = d.pools(q.loc[ids])
            a = np.array([pos[x] for x in np.array(ids)[cas.to_numpy(bool)]])
            b = np.array([pos[x] for x in np.array(ids)[con.to_numpy(bool)]])
            out[(coh, t)] = (a, b)
    return out


def score():
    limit_threads(4)
    t0 = time.time()
    co, sets, designs, M = setup()
    cohorts = ML.load_cohorts()
    pw = pd.read_csv(OUT / "power.csv").set_index("condition")
    comp = pd.read_csv(OUT / "companions.csv")
    null = pd.read_csv(OUT / "null.csv")
    keys = sorted({k for a in designs.values() for k in a})
    st = {gs: {k: co[k[0]].zv(w, k[1]) for k in keys} for gs, w in (("14 genes", WF), ("core", L3.CORE))}
    gc = {k: gc_score(co[k[0]], k[1]) for k in keys}
    L3.say("states", f"({time.time() - t0:.0f}s)")
    rows, crows = [], []
    for i, (c, arms) in enumerate(designs.items()):
        base = {k: s.arms() for k, s in arms.items()}
        q14, qc = Q_of(st["14 genes"], base), Q_of(st["core"], base)
        rng = np.random.default_rng([SEED, i, 3])
        bt = []
        for _ in range(2000):
            ab = {}
            for k, s in arms.items():
                j = rng.integers(0, len(s.n), len(s.n))
                S = s.S[j]
                ab[k] = (S[:, 0], S[:, 1:][S[:, 1:] >= 0])
            bt.append(Q_of(st["14 genes"], ab))
        lo, hi = np.nanquantile(bt, [0.05, 0.95])
        nl = null[null.condition == c].Q.to_numpy()
        mt = comp[(comp.condition == c) & (comp.kind == "matched")].Q
        rs = comp[(comp.condition == c) & (comp.kind == "responsive")].Q
        pa = pooled_arms(cohorts, c, co)
        pa = {k: v for k, v in pa.items() if k in st["14 genes"]}
        dgc = [float(np.mean(gc[k][a]) - np.mean(gc[k][b])) for k, (a, b) in base.items()]
        r = dict(condition=c, cells=len(arms), Q=q14, amp_ratio=float(np.exp(q14 / 2)), lo90=float(lo), hi90=float(hi),
                 amp_lo90=float(np.exp(lo / 2)), amp_hi90=float(np.exp(hi / 2)),
                 Q_star=q14 - float(mt.median()), pct_matched=float((mt < q14).mean() * 100),
                 pct_responsive=float((rs < q14).mean() * 100),
                 p_label=float((1 + np.sum(nl <= q14)) / (1 + len(nl))),
                 matched_q05=float(mt.quantile(0.05)), label_q05=float(np.quantile(nl, 0.05)),
                 Q_core=qc, amp_ratio_core=float(np.exp(qc / 2)),
                 Q_pooled=Q_of(st["14 genes"], pa), Q_vent=Q_of(st["14 genes"], {k: v for k, v in base.items() if k[0] == "vent"}),
                 Q_fast=Q_of(st["14 genes"], {k: v for k, v in base.items() if k[0] == "fast"}),
                 gc_case_minus_control_median=float(np.nanmedian(dgc)),
                 decidable_H382_2=bool(pw.loc[c, "gate_H382_2"]))
        rows.append(r)
        L3.say(c, "done", f"({time.time() - t0:.0f}s)")
    t = pd.DataFrame(rows)
    t["p_holm"] = holm(t.p_label)
    t["H382_1"] = (t.Q < t.matched_q05) & (t.Q < t.label_q05) & (t.p_holm < 0.05) & (t.hi90 < 0)
    t["H382_2"] = t.decidable_H382_2 & (t.lo90 > -MARGIN) & (t.hi90 < MARGIN)
    t.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.round(3).T.to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    stage = ap.parse_args().stage
    OUT.mkdir(parents=True, exist_ok=True)
    power() if stage == "power" else score()


if __name__ == "__main__":
    main()
