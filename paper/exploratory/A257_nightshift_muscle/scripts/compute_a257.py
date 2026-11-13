"""A257: after a 12 h shift of behaviour, does the human muscle clock stay on central time or follow behaviour? (GSE106800)

Registered in paper/docs/registered_experiments.md (A257; not blind, A124 scored the same data). Unit = subject x condition;
delta = z(19:00) - z(07:00), z = L y with the unweighted GTEx projector over all arrays (y = log2 x ln 2). H257.1: rho of A252
for the aligned units; H257.2: shift s = arg(sum delta_mis conj(delta_ctl)) / w; H257.3: rho_mis / rho_ctl.

    python compute_a257.py [--placebo]

Choices not fixed by the registration
- Genes: main probe sets of the A124 probe map, mean log2 per canonical gene (with the eight A124 aliases); genome-wide means
  for the companions from the same table.
- Exact nulls enumerate all 2^n label swaps (n = 12 aligned units for H257.1, 11 subjects for H257.3).
- Companion sets are drawn with `clocksync.paired.pseudo_waveforms` from the mean log expression over all arrays (seed 257).
- `--placebo`: the two time labels are swapped at random once within every unit before anything else (output/placebo/).
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

limit_threads(4)
A124 = PAPER / "exploratory/A124_night_shift"
GEO = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/geo/GSE106800")
OUT = HERE.parent / "output"
SEED, N_SETS, NB = 257, 200, 2000
WF = cs.gtex()
W = 2 * np.pi / 24
EXTRA = {"FAM46B": "TENT5B", "MTERFD3": "MTERF2", "HIST2H2BE": "H2BC21", "C12ORF23": "TMEM263", "SEPT4": "SEPTIN4",
         "C10ORF116": "ADIRF-AS1", "CENPC1": "CENPC", "FAM134B": "RETREG1"}


def canon(s):
    c = canonical(s)
    return EXTRA.get(c, c)


def load():
    smp = pd.read_csv(A124 / "output/samples.csv")
    smp["cond"] = np.where(smp["experimental condition"].str.contains("misaligned"), "mis", "ctl")
    smp["pm"] = (smp["time of sampling"] == "7 PM (19h)").astype(int)
    pm = pd.read_csv(A124 / "output/probe_map.csv.gz", keep_default_na=False)
    pm = pm[(pm.category == "main") & (pm.symbols_all != "")]
    pm = pm.assign(gene=pm.symbols_all.str.split(";").map(lambda s: canon(s[0])))
    expr = pd.read_csv(GEO / "GSE106800_series_matrix.txt.gz", sep="\t", comment="!", index_col=0)
    expr.index = expr.index.astype(str)
    pm = pm[pm.probe.astype(str).isin(expr.index)]
    x = expr.loc[pm.probe.astype(str), smp.gsm].groupby(pm.gene.to_numpy()).mean()
    return smp, (x.T * np.log(2.0)).astype(float)


def zstates(y, wf):
    g = [x for x in wf.genes if x in y.columns]
    w = wf.subset(g)
    L = P.projector(y[g], w, counts=False)
    b = y[g].to_numpy(float) @ L.T
    return b[:, 0] + 1j * b[:, 1]


def deltas(smp, z):
    """(subject, cond) -> delta (19:00 minus 07:00)."""
    d = {}
    for (s, c), g in smp.groupby(["subjectid", "cond"]):
        if set(g.pm) == {0, 1}:
            d[(s, c)] = z[g.index[g.pm == 1][0]] - z[g.index[g.pm == 0][0]]
    return d


CHORD = np.exp(1j * W * 19) - np.exp(1j * W * 7)


def rho(ds):
    ds = np.asarray(ds)
    return float(np.sum((ds * np.conj(CHORD)).real) / (len(ds) * abs(CHORD) ** 2))


def signs(n):
    return np.array(list(itertools.product([1.0, -1.0], repeat=n)))


def main():
    placebo = "--placebo" in sys.argv
    out = OUT / ("placebo" if placebo else "")
    out.mkdir(parents=True, exist_ok=True)
    smp, y = load()
    if placebo:
        rng = np.random.default_rng([SEED, 9])
        for _, g in smp.groupby(["subjectid", "cond"]):
            if len(g) == 2 and rng.random() < 0.5:
                smp.loc[g.index, "pm"] = 1 - smp.loc[g.index, "pm"].to_numpy()
    z = zstates(y, WF)
    d = deltas(smp, z)
    ctl = {s: v for (s, c), v in d.items() if c == "ctl"}
    both = sorted(s for s in ctl if (s, "mis") in d)
    dc, dm = np.array([ctl[s] for s in both]), np.array([d[(s, "mis")] for s in both])
    al = np.array(list(ctl.values()))
    # H257.1: swapping the two time labels of a unit flips the sign of its delta
    S1 = signs(len(al))
    nul1 = (S1 * (al * np.conj(CHORD)).real[None, :]).sum(1) / (len(al) * abs(CHORD) ** 2)
    obs1 = rho(al)
    lm = y.mean(0)
    lm = lm[~lm.index.duplicated()]
    comp = {}
    for kind, tag in (("matched", 3), ("responsive", 4)):
        g = [x for x in WF.genes if x in y.columns]
        sets = P.pseudo_waveforms(lm, N_SETS, np.random.default_rng([SEED, tag]), WF.subset(g), kind)
        comp[kind] = []
        for w in sets:
            zz = zstates(y, w)
            dd = deltas(smp, zz)
            a = np.array([v for (s, c), v in dd.items() if c == "ctl"])
            m = np.array([dd[(s, "mis")] for s in both])
            c2 = np.array([dd[(s, "ctl")] for s in both])
            comp[kind].append((rho(a), rho(c2) - rho(m)))
        comp[kind] = np.array(comp[kind])
    row1 = dict(test="H257.1 aligned rotation", units=len(al), rho=obs1, null_sd=float(nul1.std()), p=float(np.mean(nul1 >= obs1)))
    for kind in comp:
        row1.update({f"{kind}_median": float(np.median(comp[kind][:, 0])), f"{kind}_q95": float(np.quantile(comp[kind][:, 0], .95)),
                     f"p_{kind}": float((1 + np.sum(comp[kind][:, 0] >= obs1)) / (1 + N_SETS))})
    row1["passes"] = bool(row1["p"] < 0.05 and obs1 > row1["matched_q95"] and obs1 > row1["responsive_q95"])
    # H257.2
    shift = lambda a, b: float(np.angle(np.sum(a * np.conj(b))) / W)  # noqa: E731
    s_obs = shift(dm, dc)
    rng = np.random.default_rng([SEED, 5])
    bs = np.array([shift(dm[i], dc[i]) for i in (rng.integers(0, len(both), len(both)) for _ in range(NB))])
    bs = (bs - s_obs + 12) % 24 - 12 + s_obs
    lo, hi = float(np.quantile(bs, .05)), float(np.quantile(bs, .95))
    if -3 <= lo and hi <= 3:
        reading = "stays on central time"
    elif (lo <= 12 <= hi) or (lo <= -12 <= hi):
        reading = "follows behaviour"
    elif not (lo <= 0 <= hi):
        reading = f"partial shift of {s_obs:.1f} h"
    else:
        reading = "not resolved"
    row2 = dict(test="H257.2 shift of the muscle clock (h)", subjects=len(both), s=s_obs, lo90=lo, hi90=hi, reading=reading)
    # H257.3
    rc, rm = rho(dc), rho(dm)
    S3 = signs(len(both))
    pc, pmis = (dc * np.conj(CHORD)).real, (dm * np.conj(CHORD)).real
    nul3 = []
    for sg in S3:
        sw = sg < 0
        a, b = np.where(sw, pmis, pc), np.where(sw, pc, pmis)
        nul3.append((a.sum() - b.sum()) / (len(both) * abs(CHORD) ** 2))
    nul3 = np.array(nul3)
    br = np.array([rho(dm[i]) / rho(dc[i]) for i in (rng.integers(0, len(both), len(both)) for _ in range(NB))])
    row3 = dict(test="H257.3 damping rho_mis / rho_ctl", subjects=len(both), rho_ctl=rc, rho_mis=rm, R=rm / rc, lo90=float(np.quantile(br, .05)),
                hi90=float(np.quantile(br, .95)), p_ctl_minus_mis=float(np.mean(nul3 >= rc - rm)),
                matched_diff_q95=float(np.quantile(comp["matched"][:, 1], .95)), responsive_diff_q95=float(np.quantile(comp["responsive"][:, 1], .95)))
    z2 = zstates(y, WF.drop(["PER1", "NR1D1"]))
    d2 = deltas(smp, z2)
    desc = [dict(item="s without PER1 and NR1D1 (h)", value=shift(np.array([d2[(s, "mis")] for s in both]), np.array([d2[(s, "ctl")] for s in both])))]
    desc += [dict(item=f"per-subject angle of delta_mis against delta_ctl (h), {s}", value=float(np.angle(a * np.conj(b)) / W)) for s, a, b in zip(both, dm, dc)]
    pd.DataFrame([row1, row2, row3]).to_csv(out / "score_tests.csv", index=False)
    pd.DataFrame(desc).to_csv(out / "score_descriptive.csv", index=False)
    pd.set_option("display.width", 250)
    for r in (row1, row2, row3):
        print(pd.Series(r).to_string(), "\n")
    print(pd.DataFrame(desc).round(3).to_string())


if __name__ == "__main__":
    main()
