"""A97 H97.2: does gene coherence of the clock genes fall with proliferation in TCGA tumours? (+ pre-specified exploratory
survival and breast ER analyses.)

Pre-registered in docs/PLAN_internal_desync.md (archive) (A97). Tumours: all primary tumours of each project (one per patient,
any histology; the ten paired cancers plus OV). Proliferation score = mean z of Hallmark G2M checkpoint genes (clock
genes removed); purity proxy = mean z of COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRB, PTPRC, CD3E, CD14, LYZ, CD74
(high = low purity); z within cancer on log2 CPM; tertiles within cancer by rank.
H97.2: T = mean over cancers of gene coherence (`decompose.components`, A87) in the high minus the low proliferation
tertile; null: tumour order permuted within cancer (1,000), one-sided p (loss of coherence). Companion on the same
permutations: purity tertiles (low-purity minus high-purity tertile); check: proliferation tertiles within each half of
the purity score. Descriptive per tertile: log plane, log off-plane variance, coherence, clock-shaped share.
Exploratory (pre-specified): per tumour the Poisson-corrected off-plane fraction of its clock-gene deviation from the
cancer's normal mean; Cox per cancer (age, ordinal stage), inverse-variance pooling of the log HR per SD; breast
components by ER status (patient bootstrap, 1,000).

    python compute_proliferation.py [--subsample N] [--seed S]
Outputs (paper/results/figure_5/output/cancer/): a97_proliferation.csv, a97_proliferation_per_cancer.csv, a97_tertiles.csv,
a97_survival_per_cancer.csv, a97_survival_pooled.csv, a97_breast_er.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a97_tcga_load.py, a97_proliferation.py, a97_survival.py.
Local (not in clocksync): the off-plane survival score and the Cox fit (statsmodels).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.decompose import components  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer"
SEED = 97
G = len(sh.GENES)
PL = sh.plane()
ROMAN = {"I": 1, "II": 2, "III": 3, "IV": 4}


class Cohort:
    """Tumours of one cancer: log fractions of the clock genes (centred over the chosen subset) and the two scores."""

    def __init__(self, c: str, tum: pd.DataFrame):
        d = tum[np.isfinite(tum.prolif) & np.isfinite(tum.purity) & (tum.depth > 0)]
        self.c, self.d = c, d
        n = d[sh.GENES].to_numpy(float)
        self.y = np.log((n + 0.5) / d.depth.to_numpy(float)[:, None])
        self.iv = 1 / (n + 0.5)
        self.prolif, self.purity = d.prolif.to_numpy(), d.purity.to_numpy()

    def comp(self, idx) -> dict:
        y = self.y[idx] - self.y[idx].mean(0)
        pv = self.iv[idx].mean(0)
        r = components(y, pv, 1, sh.WF)
        r["share"] = sh.pooled_share(y, pv, 1, PL)
        return r

    def coh(self, idx) -> float:
        return self.comp(idx)["coherence"]


def tertile(x, k=3):
    return pd.qcut(pd.Series(x).rank(method="first"), k, labels=False).to_numpy()


def stat(co: Cohort, lab, perm=None) -> float:
    lab = lab if perm is None else lab[perm]
    return co.coh(np.flatnonzero(lab == 2)) - co.coh(np.flatnonzero(lab == 0))


def proliferation(cos, rng, n_perm, out):
    rows, per = [], []
    for co in cos:
        lp, lu = tertile(co.prolif), tertile(co.purity)
        for name, lab in (("prolif", lp), ("purity", lu)):
            for t in range(3):
                per.append(dict(cancer=co.c, score=name, tertile=["low", "mid", "high"][t], n=int((lab == t).sum()),
                                **co.comp(np.flatnonzero(lab == t))))
        rows.append(dict(cancer=co.c, n=len(co.d), spearman_prolif_purity=pd.Series(co.prolif).corr(pd.Series(co.purity), method="spearman"),
                         d_coh_prolif=stat(co, lp), d_coh_purity=stat(co, lu)))  # tertile 2 = low purity
    pc = pd.DataFrame(rows)
    pd.DataFrame(per).to_csv(out / "a97_tertiles.csv", index=False)
    obs = dict(prolif=pc.d_coh_prolif.mean(), purity=pc.d_coh_purity.mean())
    lab = [(tertile(co.prolif), tertile(co.purity)) for co in cos]
    nullp, nullu = [], []
    for _ in range(n_perm):
        ap, au = [], []
        for co, (lp, lu) in zip(cos, lab):
            perm = rng.permutation(len(lp))
            ap.append(stat(co, lp, perm))
            au.append(stat(co, lu, perm))
        nullp.append(np.mean(ap))
        nullu.append(np.mean(au))
    nullp, nullu = np.array(nullp), np.array(nullu)
    pl = lambda null, o: cs.stats.perm_p(null, o, "less", tie=0.0)
    res = [dict(test="proliferation (high - low)", estimate=obs["prolif"], null_mean=nullp.mean(),
                null_q05=np.quantile(nullp, .05), p_one_lower=pl(nullp, obs["prolif"])),
           dict(test="purity proxy (low-purity - high-purity)", estimate=obs["purity"], null_mean=nullu.mean(),
                null_q05=np.quantile(nullu, .05), p_one_lower=pl(nullu, obs["purity"]))]
    for h, hname in ((0, "high-purity half"), (1, "low-purity half")):
        sub = []
        for co in cos:
            idx = np.flatnonzero(tertile(co.purity, 2) == h)
            sub.append((co, idx, tertile(co.prolif[idx])))

        def f(perms=None):
            return np.mean([co.coh(idx[(l if p is None else l[p]) == 2]) - co.coh(idx[(l if p is None else l[p]) == 0])
                            for (co, idx, l), p in zip(sub, perms or [None] * len(sub))])

        o = f()
        nl = np.array([f([rng.permutation(len(idx)) for _, idx, _ in sub]) for _ in range(n_perm)])
        res.append(dict(test=f"proliferation within {hname}", estimate=o, null_mean=nl.mean(),
                        null_q05=np.quantile(nl, .05), p_one_lower=pl(nl, o)))
    res = pd.DataFrame(res).assign(n_cancers=len(cos), n_perm=n_perm)
    res.to_csv(out / "a97_proliferation.csv", index=False)
    pc.to_csv(out / "a97_proliferation_per_cancer.csv", index=False)
    print(res.round(4).to_string())
    print(pc.round(3).to_string(), flush=True)


# ---- exploratory: survival -----------------------------------------------------------------------------------------

def stage_num(s):
    if not isinstance(s, str) or s.endswith("X"):
        return np.nan
    m = re.match(r"Stage (IV|III|II|I)(?![IV])", s)
    return ROMAN[m.group(1)] if m else np.nan


def off_scores(tum: pd.DataFrame, nor: pd.DataFrame):
    """local: Poisson-corrected off-plane fraction (clipped to [0, 1]) and the raw fraction of each tumour's clock-gene
    deviation from the mean of the cancer's normals (A97 exploratory)."""
    M, Q = np.eye(G) - 1 / G, np.eye(G) - PL
    n = tum[sh.GENES].to_numpy(float)
    y = np.log((n + 0.5) / tum.depth.to_numpy(float)[:, None])
    nn = nor[sh.GENES].to_numpy(float)
    yn = np.log((nn + 0.5) / nor.depth.to_numpy(float)[:, None])
    d = (y - yn.mean(0)) @ M
    var = 1 / (n + 0.5) + (1 / (nn + 0.5)).mean(0) / len(nn)
    tot, off = np.sum(d ** 2, 1), np.sum((d @ Q) ** 2, 1)
    ntot, noff = var @ (M * M).sum(0), var @ ((Q @ M) ** 2).sum(0)
    return np.clip((off - noff) / np.maximum(tot - ntot, 1e-6), 0, 1), off / np.maximum(tot, 1e-12)


def cox(df: pd.DataFrame, cols):
    m = PHReg(df.time.to_numpy(float), df[cols].astype(float).to_numpy(), status=df.event.to_numpy(int),
              ties="efron").fit(disp=False)
    return m.params[0], m.bse[0]


def survival(raw: dict, out, quick: bool):
    rows = []
    for c, (tum, nor) in raw.items():
        if nor is None or not len(nor):
            print(c, "no normals: skipped")
            continue
        tum = tum[tum.depth > 0].copy()
        tum["off_c"], tum["off_raw"] = off_scores(tum, nor)
        tum["event"] = (tum.vital == "Dead").astype(int)
        tum["time"] = np.where(tum.event == 1, tum.t_death, tum.t_fu).astype(float)
        tum["age10"] = pd.to_numeric(tum.age, errors="coerce") / 10
        tum["stage_n"] = tum.stage.map(stage_num).fillna(tum.clinical_stage.map(stage_num))
        use_stage = tum.stage_n.notna().mean() > 0.5
        for sc in ("off_c", "off_raw"):
            d = tum.assign(z=(tum[sc] - tum[sc].mean()) / tum[sc].std())
            cols = ["z", "age10"] + (["stage_n"] if use_stage else [])
            d = d.dropna(subset=cols + ["time", "event"])
            d = d[d.time > 0]
            if d.event.sum() < (4 if quick else 10):
                continue
            b, se = cox(d, cols)
            rows.append(dict(cancer=c, score=sc, n=len(d), events=int(d.event.sum()), stage_used=use_stage, logHR=b, se=se,
                             HR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se),
                             p=float(2 * (1 - norm.cdf(abs(b / se)))), mean_score=d[sc].mean()))
    r = pd.DataFrame(rows)
    r.to_csv(out / "a97_survival_per_cancer.csv", index=False)
    pool = []
    for sc, g in (r.groupby("score") if len(r) else []):
        m = cs.stats.meta(g.logHR, g.se)
        pool.append(dict(score=sc, cancers=";".join(g.cancer), HR_fixed=np.exp(m["fixed"]),
                         lo_fixed=np.exp(m["fixed"] - 1.96 * m["se_fixed"]), hi_fixed=np.exp(m["fixed"] + 1.96 * m["se_fixed"]),
                         p_fixed=float(2 * (1 - norm.cdf(abs(m["fixed"] / m["se_fixed"])))), HR_re=np.exp(m["random"]),
                         lo_re=np.exp(m["lo"]), hi_re=np.exp(m["hi"]), tau2=m["tau2"], k=m["k"]))
    pool = pd.DataFrame(pool)
    pool.to_csv(out / "a97_survival_pooled.csv", index=False)
    print(pool.round(3).to_string(), flush=True)


def breast_er(tum: pd.DataFrame, rng, n_boot, out):
    b = tum[tum.depth > 0]
    y = np.log((b[sh.GENES].to_numpy(float) + 0.5) / b.depth.to_numpy(float)[:, None])
    iv = 1 / (b[sh.GENES].to_numpy(float) + 0.5)

    def comp(idx):
        return components(y[idx] - y[idx].mean(0), iv[idx].mean(0), 1, sh.WF)

    pos, neg = np.flatnonzero(b.er == "Positive"), np.flatnonzero(b.er == "Negative")
    op, on = comp(pos), comp(neg)
    boots = []
    for _ in range(n_boot):
        a, c = comp(rng.choice(pos, len(pos))), comp(rng.choice(neg, len(neg)))
        boots.append({**{f"{k}_pos": a[k] for k in a}, **{f"{k}_neg": c[k] for k in c}, **{f"d_{k}": a[k] - c[k] for k in a}})
    bb = pd.DataFrame(boots)
    est = {**{f"{k}_pos": op[k] for k in op}, **{f"{k}_neg": on[k] for k in on}, **{f"d_{k}": op[k] - on[k] for k in op}}
    res = pd.DataFrame({k: dict(estimate=v, lo=bb[k].quantile(.025), hi=bb[k].quantile(.975)) for k, v in est.items()}).T
    res.assign(n_pos=len(pos), n_neg=len(neg)).to_csv(out / "a97_breast_er.csv")
    print(res.round(3).to_string())


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm, n_boot = (30, 30) if run["quick"] else (1000, 1000)
    raw = {c: sh.load_tumours(c, args.subsample, args.seed) for c in sh.PROJECTS}
    proliferation([Cohort(c, t) for c, (t, _) in raw.items()], np.random.default_rng(972 + args.seed), n_perm, out)
    survival(raw, out, run["quick"])
    breast_er(raw["BRCA"][0], np.random.default_rng(973 + args.seed), n_boot, out)


if __name__ == "__main__":
    main()
