"""A340 scoring (registered in `paper/docs/registered_experiments.md`, A340 and its power note): the improperness iota of
the fetal clock-plane states (spliced counts). One run.

iota = |sum_i (z_i^2 - nu_i)| / sum_i (|z_i|^2 - t_i) (`power_a340.iota`; residualisation of A328 part b, exact log-count
noise); references noise-matched to each unit (power note). 90 % interval from 2,000 bootstrap resamples within cohort. H340.1 (del Valle 14 genes and core, liver 14
genes): one-sided p = share of the power stage's iota_line draws at 0.75 of the adult clock that are <= the observed
iota (Holm over the two units on the 14 genes); position against the iota_daily 5-95 % range. Descriptive: 10x cortex,
`cortex_minus1`, liver core; the angle between the main axis (arg of the pseudo-variance / 2) and the glucocorticoid-
target direction (regression of the states on the A307 GC score: FKBP5, TSC22D3, ZBTB16, KLF9), as an axis (0-90 deg).
`--check`: every gene permuted within cohort before any statistic (no real iota).

    python score_a340.py [--jobs 4] [--check]
Outputs (`../output/score[_check]/`): iota.csv, verdicts.csv.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2] / "A328_fetal_running_loop" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a328b as B  # noqa: E402  (sets the import paths of the analysis repo)
import lognoise  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector  # noqa: E402
from common import development as dev  # noqa: E402
from power_a340 import iota  # noqa: E402
from score_a337 import boot  # noqa: E402

OUT = HERE.parents[1] / "output"
POWER = OUT / "power"
SEED, N_BOOT = 3401, 2000
UNITS = ("delvalle", "liver", "cortex", "cortex_minus1")
SCORED = {("delvalle", "all"), ("delvalle", "core"), ("liver", "all")}
GC = ("FKBP5", "TSC22D3", "ZBTB16", "KLF9")


def axis_angles(s, depth, Z, wf, groups, gc) -> tuple[float, float]:
    """Main-axis angle (deg, mod 180) and its angle to the GC direction (deg, 0-90)."""
    L = projector(pd.DataFrame(s, columns=list(wf.genes)), wf)
    y = dev.residual(np.log((s + 0.5) / depth[:, None]), Z)
    b = (y - y.mean(0)) @ L.T
    z = b[:, 0] + 1j * b[:, 1]
    dof = 1 - (1 + Z.shape[1]) / len(s)
    V = np.einsum("ag,ng,bg->nab", L, dof * lognoise.var_log(lognoise.expected(s, depth, groups)), L)
    main = np.angle(np.sum(z ** 2) - np.sum(V[:, 0, 0] - V[:, 1, 1] + 2j * V[:, 0, 1])) / 2
    g = dev.residual(gc[:, None], Z)[:, 0]
    g = g - g.mean()
    d = np.angle(np.sum(z * g))
    diff = np.degrees(np.abs((main - d + np.pi / 2) % np.pi - np.pi / 2))
    return float(np.degrees(main) % 180), float(diff)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    global N_BOOT
    limit_threads(4)
    t0 = time.time()
    out = OUT / ("score_check" if a.check else "score")
    out.mkdir(parents=True, exist_ok=True)
    units = [u for u in B.load_units() if u["name"] in UNITS]
    rng = np.random.default_rng(SEED)
    if a.check:
        N_BOOT = 200
        for u in units:
            g = u["obs"].cohort.to_numpy()
            s = u["s"].to_numpy(np.int64).copy()
            for c in np.unique(g):
                i = np.flatnonzero(g == c)
                for j in range(s.shape[1]):
                    s[i, j] = s[rng.permutation(i), j]
            u["s"] = pd.DataFrame(s, columns=u["s"].columns)
    print(f"load {time.time() - t0:.0f}s", flush=True)
    ref = pd.read_csv(POWER / "subsamples.csv")
    ref = ref[ref.scale == 0.75]
    q = pd.read_csv(POWER / "quantiles.csv")
    q = q[q.scale == 0.75].set_index(["unit", "waveform", "kind"])
    rows = []
    for u in units:
        obs = u["obs"]
        Z, groups, ds = B.covariates(obs), obs.cohort.to_numpy(), obs.depth_s.to_numpy(float)
        gc = np.log((u["s"][list(GC)].to_numpy(float) + 0.5) / ds[:, None]).mean(1)
        for w, wf in B.WAVEFORMS.items():
            s = u["s"][list(wf.genes)].to_numpy(np.int64)
            val = iota(s, ds, Z, wf, groups)
            bs = np.array([iota(s[b], ds[b], Z[b], wf, groups[b]) for b in boot(groups, rng, N_BOOT)])
            main_deg, gc_deg = axis_angles(s, ds, Z, wf, groups, gc)
            unit_ref = "cortex" if u["name"] == "cortex_minus1" else u["name"]
            line = ref[(ref.unit == unit_ref) & (ref.waveform == w) & (ref.kind == "line")].iota.to_numpy()
            rows.append(dict(unit=u["name"], waveform=w, n=len(obs), iota=val, iota_lo=float(np.quantile(bs, .05)),
                             iota_hi=float(np.quantile(bs, .95)),
                             p_line=float((1 + np.sum(line <= val)) / (1 + len(line))),
                             line_q05=q.loc[(unit_ref, w, "line"), "q05"],
                             daily_q05=q.loc[(unit_ref, w, "daily"), "q05"],
                             daily_q95=q.loc[(unit_ref, w, "daily"), "q95"],
                             main_axis_deg=main_deg, axis_to_gc_deg=gc_deg))
            print(u["name"], w, f"{time.time() - t0:.0f}s", flush=True)
    t = pd.DataFrame(rows)
    t["in_daily_range"] = (t.iota >= t.daily_q05) & (t.iota <= t.daily_q95)
    t.to_csv(out / "iota.csv", index=False)
    v = t[[(r.unit, r.waveform) in SCORED for r in t.itertuples()]].copy()
    v["holm"] = np.nan
    m = v.waveform == "all"
    v.loc[m, "holm"] = dev.holm(v.loc[m, "p_line"].to_numpy())
    v.loc[~m, "holm"] = v.loc[~m, "p_line"]  # del Valle core: a single secondary test
    v["H340_1"] = v.holm < .05
    v.to_csv(out / "verdicts.csv", index=False)
    pd.set_option("display.width", 220)
    print(t.round(3).to_string(index=False), v.round(4).to_string(index=False), f"runtime {time.time() - t0:.0f}s",
          sep="\n", flush=True)


if __name__ == "__main__":
    main()
