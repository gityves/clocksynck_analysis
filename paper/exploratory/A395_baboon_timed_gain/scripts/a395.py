"""A395 (baboon: does an organ's timed own clock noise grow with its gain?), registered 2026-11-26.

Baboon (Mure 2018, A263 table, 12 animals, ZT known), peripheral tissues, one `Dataset` (peripheral anchors);
`Dataset.timed_loadings("zt")`: timed gain |a| and ML own noise psi (beyond measurement noise) per organ. Precondition:
|a| above the q95 of 1,000 ZT permutations. Animal correction: organ k refitted with hours ZT_s + delta_{s,-k}, the mean
linear phase departure of animal s over its other passing organs. b_t = OLS slope of log psi on log |a|^2.

Stages:
  power  passing organs; organ-bootstrap SE of b_t (14 genes, core; corrected and uncorrected); no slope printed.
  score  b_t with 90 % organ-bootstrap interval (H395.1: corrected, 14 genes), uncorrected, core loop; phase scatter.

    paper/common/run.sh paper/exploratory/A395_baboon_timed_gain/scripts/a395.py --stage power
Outputs (`../output/`): organs.csv (per organ: |a|, psi, scatter; written in the power stage, read only at scoring),
power.csv; tests.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory" / "A263_companions" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import score_a263 as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.timed import W, departures  # noqa: E402
from common import data  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 395
KEYS = ["donor", "compartment", "depth"]
TRIO = ("PER1", "PER2", "NR1D1")
N_PERM, N_BOOT, MIN_ORGANS = 1000, 2000, 15


def load():
    pb = A.table("baboon", A.WF.genes)
    pb["donor"] = pb.donor.astype(str)
    _, per, _ = A.groups("baboon", sorted(pb.compartment.unique()))
    _, donors = data.load("baboon")
    zt = donors.zt.astype(float)
    zt.index = zt.index.astype(str)
    return pb[pb.compartment.isin(per)].reset_index(drop=True), per, zt


def timed(ds, hours, residuals=False):
    return ds.timed_loadings(hours, return_residuals=residuals)


def precondition(ds, zt, per):
    obs = timed(ds, zt).abs_amplitude.reindex(per)
    rng = np.random.default_rng([SEED, 1])
    null = np.array([timed(ds, pd.Series(rng.permutation(zt.to_numpy()), index=zt.index)).abs_amplitude.reindex(per)
                     .to_numpy(float) for _ in range(N_PERM)])
    return pd.Series(obs.to_numpy(float) > np.nanquantile(null, 0.95, axis=0), index=per)


def organ_table(ds, zt, organs):
    """Per passing organ: uncorrected and animal-corrected |a|, psi and phase scatter (h)."""
    tab, res = timed(ds, zt, residuals=True)
    load_ = {k: tab.abs_amplitude[k] * np.exp(1j * W * tab.clock_offset_h[k]) for k in organs}
    r = res.reset_index()
    r.columns = ["donor", "compartment", "res"]
    r = r[r.compartment.isin(organs)].copy()
    r["donor"] = r.donor.astype(str).str.split("#").str[0]
    r["dep"] = departures(r.res.to_numpy(), r.compartment.map(load_).to_numpy(), zt.reindex(r.donor).to_numpy())[0]
    rows = []
    for k in organs:
        oth = r[r.compartment != k].groupby("donor").dep.mean()
        tk, rk = timed(ds, zt + oth.reindex(zt.index).fillna(0.0), residuals=True)
        rr = rk.reset_index()
        rr.columns = ["donor", "compartment", "res"]
        rr = rr[rr.compartment == k]
        dk = departures(rr.res.to_numpy(), np.full(len(rr), tk.abs_amplitude[k] * np.exp(1j * W * tk.clock_offset_h[k])),
                        (zt + oth.reindex(zt.index).fillna(0.0)).reindex(rr.donor.astype(str).str.split("#").str[0])
                        .to_numpy())[0]
        rows.append(dict(organ=k, amp=tab.abs_amplitude[k], psi=tab.psi[k], amp_corr=tk.abs_amplitude[k],
                         psi_corr=tk.psi[k], scatter_h=float(np.std(r[r.compartment == k].dep, ddof=1)),
                         scatter_corr_h=float(np.std(dk, ddof=1))))
    return pd.DataFrame(rows).set_index("organ")


def slope(x, y):
    return float(np.polyfit(x, y, 1)[0])


def boot(x, y, rng, n=N_BOOT):
    out = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        if np.ptp(x[i]) > 0:
            out.append(slope(x[i], y[i]))
    return np.array(out)


def readouts(T):
    """(label, x, y) of the slopes: log psi on log |a|^2 (corrected, uncorrected); log scatter on log |a|."""
    ok = (T.psi_corr > 0) & (T.amp_corr > 0)
    ok_u = (T.psi > 0) & (T.amp > 0)
    return [("b_t corrected", np.log(T.amp_corr[ok] ** 2).to_numpy(), np.log(T.psi_corr[ok]).to_numpy()),
            ("b_t uncorrected", np.log(T.amp[ok_u] ** 2).to_numpy(), np.log(T.psi[ok_u]).to_numpy()),
            ("scatter slope corrected", np.log(T.amp_corr[ok]).to_numpy(), np.log(T.scatter_corr_h[ok]).to_numpy())]


def tables():
    pb, per, zt = load()
    out = {}
    for lab, wf in (("14 genes", A.WF), ("core", A.WF.drop(list(TRIO)))):
        ds = cs.Dataset(pb[KEYS + list(wf.genes)], anchors=per, waveform=wf, name=f"a395 {lab}", min_donors=5)
        passing = precondition(ds, zt, per)
        organs = list(passing.index[passing])
        out[lab] = (organ_table(ds, zt, organs), len(per), len(organs))
        print(lab, f"{len(organs)} of {len(per)} organs pass", flush=True)
    return out


def power():
    t0 = time.time()
    rows = []
    for lab, (T, n_per, n_pass) in tables().items():
        T.to_csv(OUT / f"organs_{lab.replace(' ', '_')}.csv")
        for j, (name, x, y) in enumerate(readouts(T)):
            se = float(np.std(boot(x, y, np.random.default_rng([SEED, 2, j])), ddof=1))
            rows.append(dict(readout=lab, slope=name, organs_tissues=n_per, organs_pass=n_pass, organs_used=len(x),
                             se=se, mde_1_minus_b=2.49 * se, decidable=bool(n_pass >= MIN_ORGANS and 2.49 * se <= 1)))
    P = pd.DataFrame(rows)
    P.to_csv(OUT / "power.csv", index=False)
    print(P.round(3).to_string(index=False), f"({time.time() - t0:.0f}s)")


def score():
    P = pd.read_csv(OUT / "power.csv")
    rows = []
    for lab in ("14 genes", "core"):
        T = pd.read_csv(OUT / f"organs_{lab.replace(' ', '_')}.csv", index_col="organ")
        for j, (name, x, y) in enumerate(readouts(T)):
            b = slope(x, y)
            bs = boot(x, y, np.random.default_rng([SEED, 2, j]))
            lo, hi = np.quantile(bs, [0.05, 0.95])
            p = P[(P.readout == lab) & (P.slope == name)].iloc[0]
            row = dict(readout=lab, slope=name, value=b, lo90=lo, hi90=hi, organs=len(x), decidable=p.decidable,
                       excludes_1=bool(hi < 1 or lo > 1), excludes_0=bool(lo > 0 or hi < 0))
            if lab == "14 genes" and name == "b_t corrected":
                row["test"] = "H395.1"
                row["passes"] = bool(p.decidable and row["excludes_1"])
            rows.append(row)
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "tests.csv", index=False)
    print(T.round(3).to_string(index=False))


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    power() if a.stage == "power" else score()


if __name__ == "__main__":
    main()
