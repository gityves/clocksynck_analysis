"""A335 part 1b: A308 (psoriasis, H308 q against the matched family) rerun with the core loop.

Reuses `score_a308.py` (cohorts, pooling weights, swaps F, offsets, q). A183's seeds draw the same 200 matched and 200
responsive sets per batch from the 14-gene waveform; each is restricted to its genes at the core positions (A329 rule).
    python a335_1b_psoriasis.py [--draws 2000] [--check]   # --check: 3 sets, 20 draws, assertions only
Output: ../output/a335_1b_tests.csv (14-gene recomputed + original + core), a335_1b_sets.csv.
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
sys.path.insert(0, str(PAPER))
spec = importlib.util.spec_from_file_location("a308", PAPER / "exploratory/A308_psoriasis_calibration/scripts/score_a308.py")
M8 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M8)
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

A = M8.A
OUT = HERE.parents[1] / "output"
DROP = ["PER1", "PER2", "NR1D1"]


def tests(tab, O, kinds, draws, label):
    rows = []
    clock = O[0]
    for kind in ("matched", "responsive"):
        o = O[kinds == kind][:, 0]
        sd = tab.loc[tab.kind == kind, "noise_sd_h"].to_numpy()
        sc = tab.noise_sd_h.iloc[0]

        def q_of(idx):
            tau2 = max(0.0, np.var(o[idx], ddof=1) - np.mean(sd[idx] ** 2))
            return norm.cdf((clock[0] - o[idx].mean()) / np.sqrt(tau2 + sc ** 2)), tau2

        q, tau2 = q_of(np.arange(len(o)))
        brng = np.random.default_rng([M8.SEED, 1])
        qb = [q_of(brng.integers(0, len(o), len(o)))[0] for _ in range(2000)]
        M = O[kinds == kind]
        rows.append(dict(waveform=label, family=kind, clock_offset=clock[0], clock_noise_sd=sc,
                         q05_obs=float(np.quantile(M, 0.05, axis=0)[0]), pct_obs=float(np.mean(o <= clock[0])),
                         v308_rejection=float(np.mean(clock[1:] < np.quantile(M, 0.05, axis=0)[1:])), sets_mean=float(o.mean()),
                         sets_var=float(np.var(o, ddof=1)), tau=float(np.sqrt(tau2)), q=float(q),
                         q_lo90=float(np.quantile(qb, 0.05)), q_hi90=float(np.quantile(qb, 0.95)),
                         shift_holds=bool(q < 0.05) if kind == "matched" else np.nan))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    n_sets, draws = (3, 20) if a.check else (M8.N_SETS, a.draws)
    limit_threads(1)
    s = A.load()
    means = A.gene_means(s)
    wf = P.gtex()
    core = wf.drop(DROP)
    core_idx = [k for k, g in enumerate(wf.genes) if g in core.genes]
    assert len(core_idx) == len(wf.genes) - 3
    D = A.Design(s)
    coh = pd.read_csv(M8.A183_DIR / "output/h183_1_paired_cohorts.csv", dtype={"cohort": str})
    used = str(pd.read_csv(M8.A183_DIR / "output/h183_1_paired_pooled.csv", index_col=0).loc["cohorts", "value"]).split(";")
    w = dict(zip(used, 1 / coh.set_index("cohort").loc[used, "se"].to_numpy() ** 2))
    ps = s[s.diagnosis == "psoriasis"]
    data = {}
    for k in used:
        case = ps[(ps.cohort == k) & (ps.condition == "lesional")].set_index("unit")
        ref = ps[(ps.cohort == k) & (ps.condition == "nonlesional")].set_index("unit")
        u = sorted(set(case.index) & set(ref.index))
        data[k] = (case.loc[u], ref.loc[u])
    restrict = lambda x, j: Waveform(tuple(x.genes[k] for k in core_idx), core.amp, core.psi, None, None, f"core{j}")  # noqa: E731
    full = [("clock", 0, {k: wf for k in used})]
    cor = [("clock", 0, {k: core for k in used})]
    for kind, sub in (("matched", 7), ("responsive", 8)):
        crng = np.random.default_rng([A.SEED, sub])
        per = {k: P.pseudo_waveforms(means[k], n_sets, crng, wf, kind) for k in D.keys}
        full += [(kind, j, {k: per[k][j] for k in used}) for j in range(n_sets)]
        cor += [(kind, j, {k: restrict(per[k][j], j) for k in used}) for j in range(n_sets)]
    rng = np.random.default_rng(M8.SEED)
    F = {k: np.vstack([np.zeros(len(data[k][0]), bool), rng.random((draws, len(data[k][0]))) < 0.5]) for k in used}
    kinds = np.array([k for k, _, _ in full])

    def run(sets):
        O = np.empty((len(sets), draws + 1))
        for i, (_, _, wfs) in enumerate(sets):
            O[i] = M8.offsets({k: M8.complex_pairs(*data[k], wfs[k]) for k in used}, w, F)
        tab = pd.DataFrame(dict(kind=kinds, set=[j for _, j, _ in sets], offset_h=O[:, 0], noise_sd_h=[M8.wsd(o[1:]) for o in O]))
        return tab, O

    tf, Of = run(full)
    if a.check:
        assert all(len(wfs[used[0]].genes) == 11 for _, _, wfs in cor)
        run(cor)  # code path only, result discarded
        print("check ok")
        return
    Tf = tests(tf, Of, kinds, draws, "14gene")
    ref = pd.read_csv(M8.HERE.parents[1] / "output/tests.csv").set_index("family")
    for _, r in Tf.iterrows():
        o = ref.loc[r.family]
        for c in ("clock_offset", "q", "q05_obs", "pct_obs", "tau", "q_lo90", "q_hi90"):
            if not np.isclose(r[c], o[c], atol=1e-6):
                raise SystemExit(f"A308 not reproduced: stop before reading the core loop: {r.family} {c} {r[c]} vs {o[c]}")
    print("14-gene reproduction OK", flush=True)
    tc, Oc = run(cor)
    Tc = tests(tc, Oc, kinds, draws, "core")
    T = pd.concat([Tf, Tc], ignore_index=True)
    T.to_csv(OUT / "a335_1b_tests.csv", index=False)
    pd.concat([tf.assign(waveform="14gene"), tc.assign(waveform="core")]).to_csv(OUT / "a335_1b_sets.csv", index=False)
    pd.set_option("display.width", 250)
    print(T.round(4).T.to_string())


if __name__ == "__main__":
    main()
