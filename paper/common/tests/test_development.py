"""`common.development` on a synthetic developmental series: an age trend shared by the compartments must not pass
as a clock, and a clock riding on it must be found and be clock-specific (about 1 min)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "clocksync"))
import clocksync as cs  # noqa: E402
from clocksync._core.controls import RESPONSIVE_GENES  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import development as dev  # noqa: E402

limit_threads(2)
WF = cs.gtex()
N, COMPS = 40, ("heart", "liver", "kidney")
COV = lambda g: dev.design(g, "age", 2)  # noqa: E731


def series(amplitude, seed=0):
    """Clock genes (amplitude 0 = no clock; `simulate` adds a shift common to the clock genes of each sample) and 300
    other genes, all with a gene-specific trend in age shared by the three compartments of a donor."""
    rng = np.random.default_rng(seed)
    comps = {k: cs.Compartment(1.0, psi=0.002, mesor_sd=0.0) for k in COMPS}
    tab, _ = cs.simulate(WF, comps, n=N, seed=seed, window_h=(0, 24), tau=0.0, amplitude=amplitude)
    tab = tab.rename(columns={"resection": "donor"})
    age = pd.Series(rng.uniform(4, 20, N), index=sorted(tab.donor.unique()))
    a = ((tab.donor.map(age) - 12) / 5).to_numpy()
    other = [f"G{j:03d}" for j in range(300 - len(RESPONSIVE_GENES))] + list(RESPONSIVE_GENES)
    lam = np.exp(rng.uniform(2, 7, len(other)) + np.outer(a + 0.3 * a ** 2, rng.normal(0, 0.5, len(other))))
    genes = list(WF.genes)
    clock = rng.poisson(tab[genes].to_numpy() * np.exp(np.outer(a + 0.3 * a ** 2, rng.normal(0, 0.5, len(genes)))))
    counts = pd.DataFrame(np.column_stack([clock, rng.poisson(lam)]), columns=genes + other)
    obs = tab[["donor", "compartment", "depth"]].assign(age=tab.donor.map(age).to_numpy(),
                                                        stratum=np.where(tab.donor.map(age) < 12, "early", "late"))
    return counts, obs


def gate(counts, obs, cov):
    tab, _ = dev.adjust(dev.table(counts, obs, WF.genes), WF.genes, cov)
    return dev.shared(tab, WF, 200, 1)


def test_age_trend_alone_is_removed():
    counts, obs = series(amplitude=0.0)
    raw = gate(counts, obs.assign(stratum="one"), lambda g: np.empty((len(g), 0)))
    adj = gate(counts, obs, COV)
    assert raw.p_pairing < 0.01 and adj.p_pairing > 0.05 and adj.pairing_z < 3
    r = dev.l1_test(counts, obs[obs.compartment == "liver"], COV, WF, 200, 20, [1])
    assert r["p"] > 0.05 and r["k"] == 3
    assert r["R_unsigned"] > 0.6 > r["C"]  # the common shift passes as unsigned coherence, not as signed


def test_clock_on_age_trend_is_found_and_specific():
    counts, obs = series(amplitude=1.0)
    r = dev.l1_test(counts, obs[obs.compartment == "liver"], COV, WF, 200, 20, [1])
    assert r["p"] < 0.05 and r["clock_specific"] and r["logo_min"] > 0.5
    tab, ws = dev.l2_prepare(counts, obs, COV, WF, 10, [2])
    g = dev.shared(tab, WF, 200, 1)
    pseudo = {k: pd.DataFrame([dev.l2_pseudo(counts, obs, COV, w, 200, 1) for w in v]) for k, v in ws.items()}
    row = dev.l2_row(g, pseudo)
    assert row["gate"] and row["clock_specific"] and row["n"] == N


def test_not_estimable():
    counts, obs = series(amplitude=1.0)
    few = obs[obs.donor.isin(sorted(obs.donor.unique())[:4])]
    assert gate(counts, few, COV) is None and dev.l2_row(None, {}) == dict(estimable=False)


def test_score_runs_l1_and_l2_with_companions():
    counts, obs = series(amplitude=1.0)
    l1 = {"liver": (obs.compartment == "liver", COV)}
    l2 = {"organs": (obs.compartment.notna(), COV), "few": (obs.donor.isin(sorted(obs.donor.unique())[:4]), COV)}
    t1, t2, pseudo = dev.score(counts, obs, l1, l2, WF, seed=3, n_perm=200, n_pseudo=4, n_perm_pseudo=200, jobs=2)
    t2 = t2.set_index("test")
    assert t1.p[0] < 0.05 and t2.gate["organs"] and not t2.estimable["few"]
    assert len(pseudo) == 16 and pseudo[pseudo.test == "few"].pairing_z.isna().all()
    v = dev.verdicts(t1, t2.reset_index(), {"H1": ["liver"], "H2": ["organs", "few"]}).set_index("test")
    assert list(v.verdict) == ["PASS", "PASS", "not estimable"] and list(v.level) == ["L1", "L2", "L2"]
    assert v.p_holm["organs"] == 2 * v.p["organs"] and v.p["organs"] == t2.loc["organs", ["p_pairing", "p_acrophase"]].max()


def test_holm():
    assert np.allclose(dev.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])


def test_verdicts_without_l2():
    l1 = pd.DataFrame(dict(test=["a", "b"], n=9, C=[0.8, 0.5], p=[0.01, 0.03], matched_pct=[100.0, 90.0],
                           responsive_pct=100.0, clock_specific=[True, False]))
    v = dev.verdicts(l1, pd.DataFrame(), {"H": ["a", "b"]})
    assert list(v.verdict) == ["PASS", "not clock-specific"] and np.allclose(v.p_holm, [0.02, 0.03])


def test_pooled_cohorts_are_permuted_within_cohort():
    counts, obs = series(amplitude=1.0)
    obs = obs[obs.compartment == "liver"].copy()
    rng = np.random.default_rng(5)
    obs["cohort"] = rng.permutation(np.repeat(list("abcd"), N // 4))
    platform = pd.DataFrame(rng.normal(0, 1, (4, counts.shape[1])), index=list("abcd"), columns=counts.columns)
    counts = np.round(counts.loc[obs.index] * np.exp(platform.loc[obs.cohort].to_numpy()))  # capture per platform
    cov = lambda g: dev.design(g, "age", 2, ("cohort",))  # noqa: E731
    r = dev.l1_test(counts, obs, cov, WF, 200, 20, [1], "cohort")
    assert r["p"] < 0.05 and r["clock_specific"] and r["k"] == 6
    assert dev.l1_test(counts, obs, COV, WF, 200, 20, [1])["C"] < r["C"]  # platform offsets left in hide the clock
    y, groups = rng.normal(size=(N, 3)), list(obs.reset_index().groupby("cohort").indices.values())
    yp = dev._permuted(y, groups, rng)
    assert not np.allclose(yp, y) and all(np.allclose(np.sort(yp[i], 0), np.sort(y[i], 0)) for i in groups)


def test_surrogates_recover_an_unrecorded_trend():
    counts, obs = series(amplitude=1.0)
    obs = obs[obs.compartment == "liver"]
    a = (obs.age - 12) / 5
    pc = dev.surrogates(counts.loc[obs.index], obs.depth, n=2, exclude=WF.genes)
    assert pc.shape == (N, 2) and np.allclose(pc.std(0), 1)
    assert abs(np.corrcoef(pc[:, 0], a + 0.3 * a ** 2)[0, 1]) > 0.9
