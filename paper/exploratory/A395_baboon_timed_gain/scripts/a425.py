"""A425 (mouse: the timed own-noise law of A397 read with mouse waveforms), registered 2026-11-28.

A397 exactly (Deota 2023, 48 mice, 16 peripheral tissues, arm as stratum; `a395.py` estimator: timed |a| and psi from
`Dataset.timed_loadings("zt")`, animal correction, precondition with 1,000 ZT permutations, organ bootstrap), with the
clock read through four waveforms: E (Wolff 2023, primary), E2 (Manella 2021 Alb-Cre), Mf (fitted on Deota itself,
descriptive) and H (GTEx, A397 reproduced). Waveforms E / E2 as A406 (`a406.external_waveforms`).

Stages:
  power  per waveform: passing tissues, per-organ tables (written, read only at scoring), organ-bootstrap SE of b_t;
         no slope printed.
  score  H425.1 (b_t corrected, 14 genes, E: interval excludes 1), H425.2 (paired b_t(H) - b_t(E) > 0), descriptive
         rows (E2, Mf, core, scatter slope, per-organ ellipse misfit under H and E).

    paper/common/run.sh paper/exploratory/A395_baboon_timed_gain/scripts/a425.py --stage power
Outputs (`../output/a425/`): organs_<wf>_<readout>.csv, power.csv, tests.csv. Aggregates only.
"""
import argparse
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(EXP / "A400_mouse_person" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a395 as B  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import fit_waveform  # noqa: E402

OUT = HERE.parents[1] / "output" / "a425"
WAVEFORMS = ("E", "E2", "Mf", "H")
SEED = 425


def waveforms(pb):
    import a406
    w = a406.external_waveforms()
    genes = list(B.A.WF.genes)
    lf = np.log((pb[genes].to_numpy(float) + 0.5) / pb.depth.to_numpy(float)[:, None])
    zt = B.S_ZT.reindex(pb.donor).to_numpy(float)
    w["Mf"] = fit_waveform(lf, zt, genes, group=pb.compartment.astype(str) + "|" + pb.stratum.astype(str), name="Mf")
    w["H"] = B.A.WF
    return w


def ellipse(ds, zt, organs):
    """Per organ: second-harmonic share of the squared residual modulus over ZT (ellipse misfit), descriptive."""
    _, res = B.timed(ds, zt, residuals=True)
    r = res.reset_index()
    r.columns = ["donor", "compartment", "res"]
    r["t"] = zt.reindex(r.donor.astype(str).str.split("#").str[0]).to_numpy(float) * B.W
    out = {}
    for k in organs:
        x = r[r.compartment == k]
        m = np.abs(x.res.to_numpy()) ** 2
        out[k] = float(np.abs(np.mean(m * np.exp(-2j * x.t.to_numpy()))) / np.mean(m))
    return pd.Series(out, name="ellipse2")


def tables():
    pb, per, zt = B.load_mouse()
    B.S_ZT = zt
    wfs = waveforms(pb)
    rows = []
    for name in WAVEFORMS:
        for lab, wf in (("14 genes", wfs[name]), ("core", wfs[name].drop(list(B.TRIO)))):
            ds = cs.Dataset(pb[B.KEYS + list(wf.genes)], anchors=per, waveform=wf, name=f"a425 {name} {lab}",
                            min_donors=5)
            passing = B.precondition(ds, zt, per)
            organs = list(passing.index[passing])
            T = B.organ_table(ds, zt, organs).join(ellipse(ds, zt, organs))
            T.to_csv(OUT / f"organs_{name}_{lab.replace(' ', '_')}.csv")
            x, y = B.readouts(T)[0][1:]
            se = float(np.std(B.boot(x, y, np.random.default_rng([SEED, 2])), ddof=1))
            rows.append(dict(waveform=name, readout=lab, tissues=len(per), passing=len(organs), se_b_t=se,
                             mde_1_minus_b=2.49 * se, decidable=bool(len(organs) >= 10 and 2.49 * se <= 1)))
            print(name, lab, f"{len(organs)} of {len(per)} pass", flush=True)
    P = pd.DataFrame(rows)
    P.to_csv(OUT / "power.csv", index=False)
    print(P.round(3).to_string(index=False))


def score():
    P = pd.read_csv(OUT / "power.csv").set_index(["waveform", "readout"])
    T = {(w, r): pd.read_csv(OUT / f"organs_{w}_{r.replace(' ', '_')}.csv", index_col="organ")
         for w in WAVEFORMS for r in ("14 genes", "core")}
    rows = []
    for (w, r), t in T.items():
        for j, (name, x, y) in enumerate(B.readouts(t)):
            bs = B.boot(x, y, np.random.default_rng([SEED, 3, j]))
            lo, hi = np.quantile(bs, [0.05, 0.95])
            row = dict(waveform=w, readout=r, slope=name, value=B.slope(x, y), lo90=lo, hi90=hi, organs=len(x),
                       ellipse2_median=float(t.ellipse2.median()), decidable=bool(P.loc[(w, r), "decidable"]))
            if (w, r, name) == ("E", "14 genes", "b_t corrected"):
                row.update(test="H425.1", passes=bool(row["decidable"] and hi < 1))
            rows.append(row)
    # H425.2: paired b_t(H) - b_t(E), organs passing under both, organ bootstrap
    h, e = T[("H", "14 genes")], T[("E", "14 genes")]
    common = h.index.intersection(e.index)
    h, e = h.loc[common], e.loc[common]
    xh, yh = np.log(h.amp_corr ** 2).to_numpy(), np.log(h.psi_corr).to_numpy()
    xe, ye = np.log(e.amp_corr ** 2).to_numpy(), np.log(e.psi_corr).to_numpy()
    ok = np.isfinite(xh + yh + xe + ye)
    xh, yh, xe, ye = xh[ok], yh[ok], xe[ok], ye[ok]
    rng = np.random.default_rng([SEED, 4])
    d = []
    for _ in range(B.N_BOOT):
        i = rng.integers(0, len(xh), len(xh))
        if np.ptp(xh[i]) > 0 and np.ptp(xe[i]) > 0:
            d.append(B.slope(xh[i], yh[i]) - B.slope(xe[i], ye[i]))
    lo, hi = np.quantile(d, [0.05, 0.95])
    rows.append(dict(waveform="H - E", readout="14 genes", slope="b_t corrected (paired)",
                     value=B.slope(xh, yh) - B.slope(xe, ye), lo90=lo, hi90=hi, organs=int(ok.sum()),
                     test="H425.2", passes=bool(lo > 0)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "tests.csv", index=False)
    print(out.round(3).to_string(index=False))


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    B.use_mouse()
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    tables() if a.stage == "power" else score()


if __name__ == "__main__":
    main()
