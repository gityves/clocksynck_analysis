"""A157: is the Jammalamadaka-SenGupta (JS) correlation usable on an evenly spread 4-h grid? Simulation of the A113 designs.

Registered in paper/docs/registered_experiments.md (A157). Designs: muscle GSE108539 (10 subjects, 57 samples, clock times
12 16 20 0 4 8 h, array input) and adipose GSE104674 (12 subjects, 48 samples, ZT 15.5, 0.25, 5.75, 11.25, counts).
True internal time = recorded time + subject offset (SD 1 h) + sample noise (SD 1 h) [noise SD: choice, = clocksync.simulate jitter_h].
Samples: clock-plane state A exp(i W t_true) (A = observed clock-plane amplitude), gene log expression = mesor + X beta + N(0, tau^2)
(muscle: log2 intensities of the 13-gene panel, unweighted projector; adipose: Poisson counts at each sample's observed depth and
the observed per-gene median rate, plus clocksync.simulate's per-gene offset sd 0.5 and common shift sd 0.2). Estimator = A113's:
`clocksync.paired.projector` + `states`, cohort-centred, phase = arg of the state. The gene noise tau is calibrated BEFORE any
correlation is computed: bisection (pilot of 100 replicates, seed 1157, correlations never evaluated) so that the median
|error| of the estimated vs true internal time after the common offset equals the observed A113 median error (muscle 1.30 h,
adipose 1.11 h). Then 500 replicates (seed 157): JS and Fisher-Lee (clocksync.stats) of estimated vs true internal time (and vs
recorded time), shares of JS < 0, FL > 0.5.
Muscle waveform: the published GTEx waveform with the ARNTL row dropped (the A113 sensitivity arm), not the 13-gene refit
(which needs the controlled-access time of death).
Outputs: output/a157_calibration.csv, a157_replicates.csv, a157_summary.csv
"""
from __future__ import annotations

import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from common import data  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
W = 2 * np.pi / 24
N_REP, SEED, PILOT_SEED, N_PILOT = 500, 157, 1157, 100
TARGET = {"muscle": 1.30, "adipose": 1.11}  # observed A113 median |error| (a113_summary.csv: 1.2972, 1.1052)
cdiff = lambda a, b: (np.asarray(a) - np.asarray(b) + np.pi) % (2 * np.pi) - np.pi


def design(name):
    if name == "muscle":
        tab, _ = data._builder("muscle_timed_gse108539").build("exon")
        wf = cs.gtex().subset([g for g in cs.gtex().genes if g != "ARNTL"])
        return tab, "clock_h", wf, False
    tab, _ = data._builder("adipose_timed_gse104674").build()
    return tab, "zt_h", cs.gtex(), True


class Sim:
    def __init__(self, name):
        self.name = name
        tab, tcol, wf, self.counts = design(name)
        self.wf, self.genes = wf, list(wf.genes)
        self.blk = tab.block.astype(str).to_numpy()
        self.subj = np.unique(self.blk)
        self.hours = tab[tcol].to_numpy(float)
        self.depth = tab.depth.to_numpy(float) if self.counts else None
        x = tab[self.genes].reset_index(drop=True)
        self.x_obs = x
        L = projector(x, wf, counts=self.counts)
        beta, _ = states(x, L, wf, depth=self.depth)
        beta = beta - beta.mean(0)
        z = beta[:, 0] + 1j * beta[:, 1]
        self.amp = float(abs(np.mean(z * np.exp(-1j * W * self.hours))) * 1.0)  # observed clock-plane amplitude vs recorded time
        self.amp_obs_rms = float(np.sqrt(np.mean(abs(z) ** 2)))
        if self.counts:
            r = np.log(np.maximum(np.median(x.to_numpy(float) / self.depth[:, None], 0), 1e-12))
            self.mesor = r
        self.X = wf.X

    def draw(self, tau, rng):
        n = len(self.hours)
        off = dict(zip(self.subj, rng.normal(0, 1.0, len(self.subj))))
        t_true = self.hours + np.array([off[b] for b in self.blk]) + rng.normal(0, 1.0, n)
        phase = W * t_true
        beta = self.amp * np.column_stack([np.cos(phase), np.sin(phase)])
        G = len(self.genes)
        if self.counts:
            lr = (self.mesor + rng.normal(0, 0.5, G) + beta @ self.X.T + rng.normal(0, tau, (n, G))
                  + rng.normal(0, 0.2, (n, 1)))
            x = pd.DataFrame(rng.poisson(self.depth[:, None] * np.exp(lr)).astype(float), columns=self.genes)
        else:
            x = pd.DataFrame((beta @ self.X.T + rng.normal(0, tau, (n, G))) / np.log(2.0), columns=self.genes)
        L = projector(x, self.wf, counts=self.counts)
        b, _ = states(x, L, self.wf, depth=self.depth)
        b = b - b.mean(0)
        # the projector's sign / orientation is fixed by the waveform: phase = arg of the state
        return np.angle(b[:, 0] + 1j * b[:, 1]), phase, W * self.hours

    def med_err(self, ph, tr):
        e = ph - tr
        m = np.angle(np.exp(1j * e).mean())
        return float(np.median(np.abs(cdiff(e, m))) / W)

    def calibrate(self, hi=0.6):
        """Bisection on tau of the median error (pilot seed; no correlation evaluated). Error is monotone in tau."""
        def f(tau):
            rng = np.random.default_rng(PILOT_SEED)
            return float(np.mean([self.med_err(*self.draw(tau, rng)[:2]) for _ in range(N_PILOT)]))
        lo, rows = 1e-4, []
        while f(hi) < TARGET[self.name]:
            hi *= 2
        for _ in range(20):
            mid = (lo + hi) / 2
            v = f(mid)
            rows.append(dict(design=self.name, tau=mid, pilot_median_error_h=v, target_h=TARGET[self.name]))
            lo, hi = (mid, hi) if v < TARGET[self.name] else (lo, mid)
            if abs(v - TARGET[self.name]) < 0.01:
                break
        return mid, pd.DataFrame(rows)


def main():
    OUT.mkdir(exist_ok=True)
    cal, rep = [], []
    for name in ("muscle", "adipose"):
        s = Sim(name)
        tau, tab = s.calibrate()
        tab.to_csv(OUT / f"a157_calibration_{name}.csv", index=False)
        cal.append(dict(design=name, n_samples=len(s.hours), n_subjects=len(s.subj), amplitude_observed=s.amp,
                        rms_radius_observed=s.amp_obs_rms, tau=tau, target_median_error_h=TARGET[name],
                        pilot_median_error_h=tab.pilot_median_error_h.iloc[-1], calibration_steps=len(tab)))
        print(cal[-1], flush=True)
        rng = np.random.default_rng(SEED)
        for r in range(N_REP):
            ph, tr, rec = s.draw(tau, rng)
            rep.append(dict(design=name, rep=r, js=cs.stats.js_corr(ph, tr), fl=cs.stats.fisher_lee(ph, tr),
                            js_recorded=cs.stats.js_corr(ph, rec), fl_recorded=cs.stats.fisher_lee(ph, rec),
                            median_error_h=s.med_err(ph, tr)))
    pd.DataFrame(cal).to_csv(OUT / "a157_calibration.csv", index=False)
    R = pd.DataFrame(rep)
    R.to_csv(OUT / "a157_replicates.csv", index=False)
    rows = []
    for name, g in R.groupby("design"):
        rows.append(dict(design=name, n_rep=len(g), js_median=g.js.median(), js_q05=g.js.quantile(.05), js_q95=g.js.quantile(.95),
                         fl_median=g.fl.median(), fl_q05=g.fl.quantile(.05), fl_q95=g.fl.quantile(.95),
                         share_js_lt0=(g.js < 0).mean(), share_fl_gt05=(g.fl > 0.5).mean(),
                         share_js_lt0_and_fl_gt05=((g.js < 0) & (g.fl > 0.5)).mean(),
                         js_recorded_median=g.js_recorded.median(), share_js_recorded_lt0=(g.js_recorded < 0).mean(),
                         fl_recorded_median=g.fl_recorded.median(), median_error_h_mean=g.median_error_h.mean()))
    S = pd.DataFrame(rows).set_index("design")
    m, a = S.loc["muscle"], S.loc["adipose"]
    S["V157_1_muscle_part"] = [bool(m.share_js_lt0 >= .2 and m.share_fl_gt05 >= .8) if d == "muscle" else np.nan for d in S.index]
    S["V157_1_adipose_part"] = [bool(a.share_js_lt0 < .05) if d == "adipose" else np.nan for d in S.index]
    S["V157_1_passes"] = bool(m.share_js_lt0 >= .2 and m.share_fl_gt05 >= .8 and a.share_js_lt0 < .05)
    S.to_csv(OUT / "a157_summary.csv")
    pd.set_option("display.width", 250)
    print(S.round(3).T.to_string())


if __name__ == "__main__":
    main()
