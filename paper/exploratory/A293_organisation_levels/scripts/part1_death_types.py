"""A293 Part 1 (registered in `paper/docs/registered_experiments.md`): is the preserved organ arrangement of slow deaths and
of ventilated living-brain donors clock-specific?

A109 / A132 data (GTEx v10, 12 body tissues, donors with >= 2 tissues and a time of death). Per group, Fisher-Lee rho of
the pair phases (pairs with >= 20 shared donors in the group and in fast deaths) with those of fast deaths, for the
clock genes, 200 expression-matched pseudo-clock sets per tissue (A109 `Tissues(200, seed 0)`) and 200 stress-responsive
sets (`compute_h109_responsive.add_sets`), same donors and pairs. Hardy 0 groups at the A132 equal-n draw.
Aggregates only (GTEx time of death and phenotypes are controlled access).

    python part1_death_types.py [--subsample N]
Outputs (`../output/part1/`): reproduction.csv, groups.csv (one row per group), sets.csv (per set and group).
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A109_coupling_nulls" / "scripts"))
sys.path.insert(0, str(EXP / "A132_brain_death" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
import compute_brain_death as B  # noqa: E402
import compute_coupling_nulls as N  # noqa: E402
import compute_h109_responsive as R  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output" / "part1"
SEED = 0
N_FL = 10_000


def raw(T, kind, k):
    Z = np.zeros((len(T.donors), len(H.BODY)), complex)
    M = np.zeros((len(T.donors), len(H.BODY)))
    for j, t in enumerate(H.BODY):
        r = T.R[t]
        m = r >= 0
        Z[m, j] = T.Z[kind][t][k][r[m]]
        M[m, j] = 1
    return Z, M


def rho(Z, M, fast, g, mask):
    return float(ST.fisher_lee(B.phase(Z, M, fast, mask), B.phase(Z, M, g, mask)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--subsample", type=int, default=None)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    n_sets = 20 if a.subsample else 200
    T = H.Tissues(n_sets, SEED, a.subsample)
    T.build_states(np.random.default_rng(1), 0, real=True)
    R.add_sets(T, n_sets, SEED)
    print(f"states built {time.time() - t0:.0f}s", flush=True)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    ns = min(len(od), len(pm))
    rd = np.random.default_rng(132 + SEED)
    od_ = rd.choice(od, ns, replace=False) if len(od) > ns else od  # A132 score(): organ donors first
    pm_ = rd.choice(pm, ns, replace=False) if len(pm) > ns else pm
    groups = {"slow_hardy4": T.group_donors("H4"), "ventilated_postmortem": pm_, "hardy3": T.group_donors("H3"),
              "organ_donor": od_}
    Zc, Mc = raw(T, "clock", 0)
    masks = {g: B.pair_mask(Mc, fast, d) for g, d in groups.items()}

    rep = []
    for g, d in groups.items():
        r_mask = rho(Zc, Mc, fast, d, masks[g])
        r_all = float(ST.fisher_lee(H.pair_phase(N._H_of(Zc, Mc, fast)), H.pair_phase(N._H_of(Zc, Mc, d))))
        rep.append(dict(group=g, n=len(d), pairs=int(masks[g].sum()), rho_mask=r_mask, rho_all_pairs=r_all))
    rep = pd.DataFrame(rep)
    rep.to_csv(OUT / "reproduction.csv", index=False)
    print(rep.round(3).to_string(index=False), flush=True)
    want = {"slow_hardy4": 0.784, "ventilated_postmortem": 0.690, "organ_donor": 0.161}
    ok = all(round(float(rep.set_index("group").loc[g, c]), 3) == v
             for g, v in want.items() for c in (("rho_all_pairs",) if g == "slow_hardy4" else ("rho_mask",)))
    if not ok and not a.subsample:
        raise SystemExit("A109 / A132 clock rho not reproduced: stop (registered check)")

    rng = np.random.default_rng([293, 1])
    rows, sets = [], []
    for g, d in groups.items():
        m = masks[g]
        r0 = rho(Zc, Mc, fast, d, m)
        Hg = N._H_of(Zc, Mc, d)
        th_f = B.phase(Zc, Mc, fast, m)
        null = np.array([ST.fisher_lee(th_f, H.permute_labels(Hg, rng.permutation(len(H.BODY)))[m]) for _ in range(N_FL)])
        row = dict(group=g, n=len(d), pairs=int(m.sum()), rho_clock=r0, p_labels=N.pval(r0, null))
        for kind, lab in (("pseudo", "matched"), ("resp", "responsive")):
            rs = []
            for k in range(n_sets):
                Z, M = raw(T, kind, k)
                rs.append(rho(Z, M, fast, d, m))
                sets.append(dict(group=g, family=lab, set=k, rho=rs[-1]))
            rs = np.array(rs)
            row.update({f"{lab}_median": float(np.nanmedian(rs)), f"{lab}_q95": float(np.nanquantile(rs, 0.95)),
                        f"{lab}_pi": float(np.mean(rs < r0)), f"{lab}_n": int(np.isfinite(rs).sum())})
        row["passes"] = bool(g in ("slow_hardy4", "ventilated_postmortem") and r0 > row["matched_q95"] and row["p_labels"] < 0.05)
        row["test"] = {"slow_hardy4": "H293.1", "ventilated_postmortem": "H293.2"}.get(g, "descriptive")
        rows.append(row)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, f"{time.time() - t0:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "groups.csv", index=False)
    pd.DataFrame(sets).to_csv(OUT / "sets.csv", index=False)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
