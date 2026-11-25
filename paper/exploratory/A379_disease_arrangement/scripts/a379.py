"""A379 (registered 2026-11-25): does a chronic disease move the organ arrangement on the brain-death scale?

A317 machinery (A109 `Tissues` states, pair phases, Fisher-Lee rho against a fast-death reference) on donor-matched
case and control arms of seven conditions (A164 / A233 labels and matching rule), reference = fast deaths in neither
arm, equal-n control draws (one control per case, 20 draws).

Stages:
  power  matching counts, pairs passing, SD of D_rho under 2000 within-set relabellings, MDE = 2.49 SD; no observed
         arrangement value of any arm.
  score  (written after the power note, before any arm value) per condition: the same matched sets and the same first
         control draw (hence the same pairs) as the power stage, 20 control draws; D_rho, one-sided relabelling p
         (2000), Holm over the seven; the same D_rho for 200 matched and 200 responsive sets (same arms and draws);
         H379.1 and H379.2; descriptive: each arm's rho and the noise-only expectation (A374 rule within the reduced
         reference: n = min(cases, half the reference) fast deaths drawn against the rest of the reference, 200 draws).

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A379_disease_arrangement/scripts/a379.py --stage power | score
Outputs (`../output/`): power.csv | tests.csv, companions.csv. Controlled access: labels, ages and matched sets stay in memory; aggregates only.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A317_braindeath_confirm" / "scripts"))
sys.path.insert(0, str(EXP / "A231_age_coupling" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import matched_label as ML  # noqa: E402
import score_a317 as A  # noqa: E402
from score_a317 import B, H, raw, rho  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
N_DRAWS, N_PERM, SEED = 20, 2000, 379
CONDITIONS = ("t2d", "obese", "htn", "ihd", "copd", "depression", "smoking")
BRAIN_DEATH = 0.621  # A317 D_rho


def pools(q, c):
    """(case mask, control mask, controls per case): A164 definitions for t2d / obese, A233 for the others."""
    if c == "obese":
        return q.bmi >= 30, (q.bmi >= 18.5) & (q.bmi < 25), 1
    return q[c] == 1, q[c] == 0, 3


def covariates(T):
    """Donor table (hardy, male, ischaemic time, age, bmi, labels, cohort) for the donors of the A109 tissue set."""
    rows = []
    for coh, (_, q, _) in ML.load_cohorts().items():
        rows.append(q.assign(cohort=coh))
    q = pd.concat(rows)
    q = q[~q.index.duplicated()]
    return q.reindex(T.donors)


def matched_sets(q, T, c):
    pos = {d: i for i, d in enumerate(T.donors)}
    sets = []
    for coh in ("vent", "fast"):
        qc = q[q.cohort == coh]
        cas, con, ratio = pools(qc, c)
        sets += [[pos[d] for d in s] for s in ML.match_ratio(qc, cas, con, ratio=ratio)]
    return sets


def arms(sets, rng):
    """(cases, one control per case) for one draw."""
    return np.array([s[0] for s in sets]), np.array([rng.choice(s[1:]) for s in sets])


def relabel(sets, rng):
    out_c, out_k = [], []
    for s in sets:
        p = rng.permutation(s)
        out_c.append(p[0])
        out_k.append(p[1])
    return np.array(out_c), np.array(out_k)


def power():
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    Z, M = raw(T, "clock", 0, T.R)
    q = covariates(T)
    fast = T.group_donors("fast")
    print(f"states and covariates ({time.time() - t0:.0f}s)", flush=True)
    rows = []
    for c in CONDITIONS:
        rng = np.random.default_rng([SEED, CONDITIONS.index(c)])
        sets = matched_sets(q, T, c)
        used = {i for s in sets for i in s}
        ref = np.array([i for i in fast if i not in used])
        cases, ctl0 = arms(sets, rng)
        mask = B.pair_mask(M, ref, cases, ctl0)
        null = []
        for _ in range(N_PERM):
            a, b = relabel(sets, rng)
            null.append(rho(Z, M, ref, b, mask) - rho(Z, M, ref, a, mask))
        null = np.array(null)
        sd = float(np.nanstd(null, ddof=1))
        coh = q.cohort.to_numpy()
        rows.append(dict(condition=c, sets=len(sets), cases_vent=int(np.sum(coh[cases] == "vent")),
                         cases_fast=int(np.sum(coh[cases] == "fast")), controls=len(used) - len(sets),
                         reference=len(ref), pairs=int(mask.sum()), null_sd=sd, mde=2.49 * sd,
                         q95=float(np.nanquantile(null, 0.95)),
                         decidable_h2=bool(mask.sum() >= 30 and 2.49 * sd <= BRAIN_DEATH)))
        print(pd.Series(rows[-1]).to_dict(), f"({time.time() - t0:.0f}s)", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.round(3).to_string(index=False))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]))
    out = np.empty_like(adj)
    out[o] = adj
    return out


def d_of(Z, M, ref, cases, draws, mask):
    return float(np.median([rho(Z, M, ref, d, mask) for d in draws])) - rho(Z, M, ref, cases, mask)


def score():
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    A.RS.add_sets(T, 200, 0)
    Z, M = raw(T, "clock", 0, T.R)
    sets_kz = {(kind, k): raw(T, kind, k, T.R) for kind in ("pseudo", "resp") for k in range(200)}
    q = covariates(T)
    fast = T.group_donors("fast")
    print(f"states, sets and covariates ({time.time() - t0:.0f}s)", flush=True)
    rows, comp = [], []
    for c in CONDITIONS:
        i = CONDITIONS.index(c)
        rng = np.random.default_rng([SEED, i])  # as the power stage: cases, first draw, pairs
        sets = matched_sets(q, T, c)
        used = {j for s in sets for j in s}
        ref = np.array([j for j in fast if j not in used])
        cases, ctl0 = arms(sets, rng)
        mask = B.pair_mask(M, ref, cases, ctl0)
        draws = [ctl0] + [arms(sets, rng)[1] for _ in range(N_DRAWS - 1)]
        D = d_of(Z, M, ref, cases, draws, mask)
        rng2 = np.random.default_rng([SEED, i, 2])
        null = []
        for _ in range(N_PERM):
            a, b = relabel(sets, rng2)
            null.append(rho(Z, M, ref, b, mask) - rho(Z, M, ref, a, mask))
        null = np.array(null)
        sd = float(np.nanstd(null, ddof=1))
        cd = {}
        for (kind, k), (Zk, Mk) in sets_kz.items():
            v = d_of(Zk, Mk, ref, cases, draws, mask)
            comp.append(dict(condition=c, family=kind, set=k, D_rho=v))
            cd.setdefault(kind, []).append(v)
        rng3 = np.random.default_rng([SEED, i, 3])
        n = min(len(cases), len(ref) // 2)
        noise = []
        for _ in range(200):
            p = rng3.permutation(ref)
            noise.append(rho(Z, M, p[n:], p[:n], mask))
        rows.append(dict(condition=c, cases=len(cases), reference=len(ref), pairs=int(mask.sum()),
                         rho_cases=rho(Z, M, ref, cases, mask),
                         rho_controls=float(np.median([rho(Z, M, ref, d, mask) for d in draws])), D_rho=D,
                         p=float((1 + np.sum(null >= D)) / (1 + N_PERM)), null_sd=sd,
                         upper_90=D + 1.645 * sd,
                         matched_q95=float(np.nanquantile(cd["pseudo"], 0.95)),
                         responsive_q95=float(np.nanquantile(cd["resp"], 0.95)),
                         noise_n=n, noise_median=float(np.median(noise)), noise_q05=float(np.quantile(noise, 0.05)),
                         noise_q95=float(np.quantile(noise, 0.95))))
        print(c, f"done ({time.time() - t0:.0f}s)", flush=True)  # values printed only in the final table
    t = pd.DataFrame(rows)
    t["p_holm"] = holm(t.p)
    t["H379_1"] = (t.D_rho > 0) & (t.p_holm < 0.05) & (t.D_rho > t.matched_q95) & (t.D_rho > t.responsive_q95)
    t["H379_2"] = t.upper_90 < BRAIN_DEATH
    t.to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(comp).to_csv(OUT / "companions.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.round(3).to_string(index=False), f"\n({time.time() - t0:.0f}s)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    stage = ap.parse_args().stage
    OUT.mkdir(parents=True, exist_ok=True)
    power() if stage == "power" else score()


if __name__ == "__main__":
    main()
