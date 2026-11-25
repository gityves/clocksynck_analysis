"""A389 note 2 (POST HOC): companion comparison on standardised pooled offsets (z = offset / SE), clock against the 200
matched sets per organ (power-stage caches), inverse-variance pooling over the organs of the scored table.
Output: output/score/note2_standardised.csv.

    python a389_note2.py
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import a389 as Q  # noqa: E402

H = Q.HOURS
OUT = Q.OUT / "score"


def se_h(R, n):
    return np.sqrt((1 - R ** 2) / (2 * n * R ** 2)) * H


def main():
    with open(Q.OUT / "cache/units.pkl", "rb") as f:
        cu = pickle.load(f)
    with open(Q.OUT / "cache/units_atlas.pkl", "rb") as f:
        au = pickle.load(f)
    units = cu + au
    org = pd.read_csv(OUT / "organs.csv")
    rows = []
    for (st, pn, rd), T in org.groupby(["sample", "pair", "readout"]):
        se_c = se_h(T.R.to_numpy(), T.n.to_numpy())
        w = 1 / se_c ** 2
        mu = float(np.sum(w * T.peak_h) / w.sum())
        z_clock = mu / float(1 / np.sqrt(w.sum()))
        offs, ses = [], []
        for organ, n in zip(T.organ, T.n):
            num = ea = eb = None
            for r in units:
                if "comp" not in r or (st, pn) not in r["comp"] or Q.ORGAN.get(r["cancer"], r["cancer"]) != organ:
                    continue
                c, x, y = r["comp"][(st, pn)][rd]
                ok = np.isfinite(x)
                c, x, y = np.where(ok, c, 0), np.where(ok, x, 0), np.where(ok, y, 0)
                num, ea, eb = (c, x, y) if num is None else (num + c, ea + x, eb + y)
            R = np.abs(num) / np.sqrt(ea * eb)
            offs.append(-np.angle(num) * H)
            ses.append(se_h(np.clip(R, 1e-3, 0.999), n))
        offs, ses = np.array(offs), np.array(ses)
        wc = 1 / ses ** 2
        mu_c = np.sum(wc * offs, 0) / wc.sum(0)
        z_c = mu_c * np.sqrt(wc.sum(0))
        rows.append(dict(sample=st, pair=pn, readout=rd, organs=len(T), clock_peak_h=mu, clock_se_h=1 / np.sqrt(w.sum()),
                         z_clock=z_clock, comp_median_h=float(np.median(mu_c)), comp_z_median=float(np.median(z_c)),
                         comp_z_q05=float(np.quantile(z_c, .05)), comp_z_q95=float(np.quantile(z_c, .95)),
                         share_abs_z_ge_clock=float(np.mean(np.abs(z_c) >= abs(z_clock))),
                         clock_minus_comp_median_h=mu - float(np.median(mu_c))))
    D = pd.DataFrame(rows)
    D.to_csv(OUT / "note2_standardised.csv", index=False)
    pd.set_option("display.width", 220)
    print(D.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
