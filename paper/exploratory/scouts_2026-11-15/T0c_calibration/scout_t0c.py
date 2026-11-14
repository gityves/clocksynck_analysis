"""Scout for the calibration set (T0c): design and timed class sizes in Manella 2021 (GSE159135). No untimed readout.

Four arms (Alb-Cre control / hepatocyte Bmal1 knockout x ad libitum / daytime feeding), liver, WAT and lung, 12 ZT x 2
mice. Per arm, tissue and gene (mean count >= 5): timed cosinor on ZT, log((count + 0.5) / depth), BH q. Liver classes
from the recorded times (the truth a calibration would be scored against):
  local    rhythmic in control AL liver (q < 0.05), not in knockout AL liver (p > 0.2) and amplitude ratio KO / control < 0.5
  systemic rhythmic in control AL (q < 0.05) and in knockout AL liver (q < 0.05), amplitude ratio >= 0.5
  ko_only  rhythmic in knockout AL liver only (q < 0.05; control p > 0.2)
Feeding: per gene rhythmic in both AL and DF of one genotype, the timed phase shift DF minus AL (h), relative to the lung
clock shift of that genotype (median shift of the lung's rhythmic clock genes).
Also the shift of each tissue's core clock (the 14 GTEx clock genes, median per-gene shift) DF minus AL.

    python scout_t0c.py
Outputs (output/): design.csv, classes.csv (counts), shifts.csv (aggregate shifts).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "common/datasets/manella2021"))
import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
import build as M  # noqa: E402

OUT = HERE / "output"
W = 2 * np.pi / 24
CLOCK = {canonical(g) for g in cs.gtex().genes}


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def cosinor(x: pd.DataFrame, zt: np.ndarray) -> pd.DataFrame:
    keep = x.mean(axis=1) >= 5
    x = x[keep]
    y = np.log((x.to_numpy(float) + 0.5) / x.sum().to_numpy(float)[None, :]).T  # depth = sum of retained genes
    D = np.c_[np.ones(len(zt)), np.cos(W * zt), np.sin(W * zt)]
    B, *_ = np.linalg.lstsq(D, y, rcond=None)
    rss1 = ((y - D @ B) ** 2).sum(0)
    rss0 = ((y - y.mean(0)) ** 2).sum(0)
    dfr = len(zt) - 3
    p = stats.f.sf(((rss0 - rss1) / 2) / (rss1 / dfr), 2, dfr)
    return pd.DataFrame(dict(amp=np.hypot(B[1], B[2]), phase_h=np.mod(np.arctan2(B[2], B[1]), 2 * np.pi) / W, p=p,
                             q=bh(p)), index=[canonical(g) for g in x.index])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res, design = {}, []
    for t in M.TISSUES:
        for g in M.GENOTYPES:
            x = M.read_tissue(t, g)
            x = x.groupby(x.index.map(canonical)).sum()
            for feed in ("AL", "DF"):
                cols = [c for c in x.columns if f"_{feed}_" in c]
                zt = np.array([int(c.split("_ZT")[1].split("_")[0]) for c in cols], float)
                r = cosinor(x[cols], zt)
                res[(t, g, feed)] = r
                design.append(dict(tissue=t, genotype=g, feeding=feed, mice=len(cols), zt=len(set(zt)),
                                   genes_tested=len(r), median_depth=float(x[cols].sum().median()),
                                   rhythmic_q05=int((r.q < 0.05).sum()),
                                   clock_rhythmic_q05=int(((r.q < 0.05) & r.index.isin(CLOCK)).sum())))
    pd.DataFrame(design).to_csv(OUT / "design.csv", index=False)
    wt, ko = res[("Liver", "AlbCre", "AL")], res[("Liver", "BLKO", "AL")]
    j = wt.join(ko, lsuffix="_wt", rsuffix="_ko", how="inner")
    j = j[~j.index.isin(CLOCK)]
    ratio = j.amp_ko / j.amp_wt
    cls = pd.Series("other", index=j.index)
    cls[(j.q_wt < 0.05) & (j.p_ko > 0.2) & (ratio < 0.5)] = "local"
    cls[(j.q_wt < 0.05) & (j.q_ko < 0.05) & (ratio >= 0.5)] = "systemic"
    cls[(j.q_ko < 0.05) & (j.p_wt > 0.2)] = "ko_only"
    cls.value_counts().rename_axis("class").rename("genes").to_csv(OUT / "classes.csv")
    rows = []
    for g in M.GENOTYPES:
        for t in M.TISSUES:
            a, d = res[(t, g, "AL")], res[(t, g, "DF")]
            k = a.join(d, lsuffix="_al", rsuffix="_df", how="inner")
            both = k[(k.q_al < 0.05) & (k.q_df < 0.05)]
            sh = (both.phase_h_df - both.phase_h_al + 12) % 24 - 12
            ck = sh[both.index.isin(CLOCK)]
            rows.append(dict(genotype=g, tissue=t, genes_rhythmic_both=len(both), clock_genes_both=len(ck),
                             clock_shift_median_h=float(np.median(ck)) if len(ck) else np.nan,
                             all_shift_median_h=float(np.median(sh)) if len(sh) else np.nan,
                             share_shift_gt_6h=float((sh.abs() > 6).mean()) if len(sh) else np.nan))
    pd.DataFrame(rows).to_csv(OUT / "shifts.csv", index=False)
    pd.set_option("display.width", 200)
    for f in ("design", "classes", "shifts"):
        print(pd.read_csv(OUT / f"{f}.csv").round(2).to_string(index=False), "\n")


if __name__ == "__main__":
    main()
