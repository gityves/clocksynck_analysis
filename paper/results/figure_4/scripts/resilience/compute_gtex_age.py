"""A94: ageing in GTEx, amplitude or phase? Individual amplitude, alignment between people and brain-body coherence by
age tertile, and the timed amplitude and phase slopes per decade.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A94-A99 common definitions, A94). Clock state per sample: the 14 clock
genes projected on the fixed GTEx waveform (counts, one projector per tissue, Poisson variance, library-size depth).
Radius r2 = |b - centre|^2 - v (phase-free individual amplitude; centre: timed-fit intercept in fast deaths, cohort mean in
ventilated donors); wall-time alignment R = |mean exp(i (arg(b - centre) - w tod))|; inter-organ coherence C
(Poisson-corrected |complex correlation| of centred states over shared donors). Age tertiles within each cohort (fast
deaths Hardy 1-2, ventilated Hardy 0); oldest minus youngest at equal n per tissue (median of 5 subsamples), median over
the 14 body tissues. Nulls: age permuted within Hardy class x sex (n_perm).
  H94.1 (amplitude): log mean radius, fast and ventilated, then their mean; two-sided.
  H94.2 (alignment between people): wall-time R, fast deaths; two-sided.
  H94.3 (brain-body): C(brain, body), fast deaths, one-sided (< 0).
  Secondary: timed amplitude and phase slopes per decade (complex regression b = Z g + (a + a1 x) exp(i w t); covariates
  sex, Hardy, ischaemic time), per tissue and inverse-variance pooled, fast deaths and Hardy 2 only.
Controlled access (dbGaP phs000424): exact age, Hardy class and time of death are held in memory; only aggregates are
written. Tertile cut points are age aggregates (as in the legacy summary).

    python compute_gtex_age.py [--subsample N] [--seed S]      (--subsample N: N random donors per cohort)
Outputs (paper/results/figure_4/output/resilience/): a94_summary.csv, a94_groups.csv, a94_timed_slopes.csv,
a94_age_by_hardy.csv
Legacy: legacy/desync/validation/gtex_covariates/scripts/{gtex_states,a94_age_sex_bmi}.py (A94 part; A99, A96 and the
cue-core secondary are not here). The per-sample state, radius, alignment and C estimators are local (clocksync has no
per-sample state API): local: not in clocksync.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import projector as _proj, states as _states  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "resilience"
W = 2 * np.pi / 24
K_SUB, MIN_N = 5, 20
BODY = ["Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)",
        "Esophagus - Mucosa", "Nerve - Tibial", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
        "Adipose - Subcutaneous", "Lung", "Whole Blood", "Liver", "Adrenal Gland"]
WF = cs.gtex()


# ---- per-sample states -------------------------------------------------------------------------------------------

def load_states(tissues, subsample, seed) -> pd.DataFrame:
    """One row per (tissue, donor): complex state b, Poisson variance v, tod (h), Hardy class, sex (female = 1),
    ischaemic time (h), exact age. Hardy 0-4 with a time of death. In memory only."""
    pb, don = gtex.load(tissues=tissues, hardy=(0, 1, 2, 3, 4), depth="library_size")
    don = don.join(gtex.phenotypes(["AGE"]).AGE.astype(float).rename("age_exact"))
    if subsample:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb = pb[pb.donor.isin(keep)]
    rows = []
    for t, x in pb.groupby("compartment"):
        if len(x) < 20:
            continue
        x = x.reset_index(drop=True)
        L = _proj(x, WF, counts=True)
        b, v = _states(x, L, WF, depth=x.depth.to_numpy(float))
        d = don.reindex(x.donor)
        rows.append(pd.DataFrame(dict(tissue=t, subject=x.donor.to_numpy(), b=b[:, 0] + 1j * b[:, 1], v=v,
                                      tod=d.tod.to_numpy(float), hardy=d.hardy.to_numpy(int),
                                      sex=1 - d.male.to_numpy(float), isch=d.ischaemic_h.to_numpy(float),
                                      age=d.age_exact.to_numpy(float))))
    s = pd.concat(rows, ignore_index=True)
    s["isch"] = s.isch.fillna(s.isch.median())
    return s


def timed_fit(b, t, Z, x=None):
    """Complex least squares b = Z g + (a + a1 x) exp(i w t): a, a1 with SEs and the fitted centre Z g."""
    e = np.exp(1j * W * t)
    cols = [Z.astype(complex), e[:, None]] + ([(x * e)[:, None]] if x is not None else [])
    D = np.hstack(cols)
    beta, *_ = np.linalg.lstsq(D, b, rcond=None)
    res = b - D @ beta
    s2 = np.sum(abs(res) ** 2) / max(len(b) - D.shape[1], 1)
    se = np.sqrt(np.real(np.diag(s2 * np.linalg.pinv(D.conj().T @ D))))
    k = Z.shape[1]
    return dict(a=beta[k], a1=beta[k + 1] if x is not None else np.nan, se_a1=se[k + 1] if x is not None else np.nan,
                centre=Z @ beta[:k].real + 1j * (Z @ beta[:k].imag))


def covariates(s: pd.DataFrame) -> np.ndarray:
    isch = (s.isch - s.isch.mean()) / (s.isch.std() or 1)
    Z = np.column_stack([np.ones(len(s)), s.sex, (s.hardy == 2).astype(float), isch])
    return Z[:, [0] + [j for j in range(1, Z.shape[1]) if Z[:, j].std() > 0]]


class Cohort:
    """tissue -> arrays sorted by donor: centred state z, variance v, radius r2, phase angle minus w tod."""

    def __init__(self, s: pd.DataFrame, fast: bool):
        self.t = {}
        for tis, x in s.groupby("tissue"):
            x = x.sort_values("subject")
            c = timed_fit(x.b.to_numpy(), x.tod.to_numpy(), covariates(x))["centre"] if fast else x.b.mean()
            z = x.b.to_numpy() - c
            self.t[tis] = dict(subj=x.subject.to_numpy(), z=z, v=x.v.to_numpy(), r2=abs(z) ** 2 - x.v.to_numpy(),
                               ang=np.angle(z) - W * x.tod.to_numpy())
        self.donors = s.drop_duplicates("subject").set_index("subject")

    def sub(self, tis, donors):
        d = self.t[tis]
        m = np.isin(d["subj"], donors)
        return {k: v[m] for k, v in d.items()}


def radius(c, tis, donors):
    m = c.sub(tis, donors)["r2"].mean()
    return np.log(m) if m > 0 else np.nan


def walltime_R(c, tis, donors):
    return float(abs(np.exp(1j * c.sub(tis, donors)["ang"]).mean()))


def pair_corr(c, t1, t2, donors):
    a, b = c.sub(t1, donors), c.sub(t2, donors)
    common = np.intersect1d(a["subj"], b["subj"])
    if len(common) < MIN_N:
        return np.nan
    i1, i2 = np.searchsorted(a["subj"], common), np.searchsorted(b["subj"], common)
    x, y = a["z"][i1] - a["z"][i1].mean(), b["z"][i2] - b["z"][i2].mean()
    k = len(common) - 1
    ex, ey = np.sum(abs(x) ** 2) / k - a["v"][i1].mean(), np.sum(abs(y) ** 2) / k - b["v"][i2].mean()
    return float(abs(np.sum(x * np.conj(y)) / k) / np.sqrt(max(ex, 1e-12) * max(ey, 1e-12)))


def coherence(c, donors, set1, set2=None):
    pairs = itertools.combinations(set1, 2) if set2 is None else itertools.product(set1, set2)
    v = [pair_corr(c, a, b, donors) for a, b in pairs if a in c.t and b in c.t]
    v = [x for x in v if np.isfinite(x)]
    return float(np.mean(v)) if v else np.nan


def tertiles(x: pd.Series, rng) -> pd.Series:
    r = (x + rng.uniform(0, 1e-6, len(x))).rank(method="first").to_numpy() - 1
    return pd.Series(np.floor(3 * r / len(x)).astype(int), index=x.index)


def contrast(c, lab, A, B, fn, tissues, seed, per_tissue=True):
    """Median over tissues of fn(A) - fn(B) at equal n (median of K_SUB subsamples); or one value over the set."""
    rng = np.random.default_rng(seed)
    dA, dB = lab.index[lab == A].to_numpy(), lab.index[lab == B].to_numpy()
    if not per_tissue:
        n = min(len(dA), len(dB))
        return float(np.nanmedian([fn(c, rng.choice(dA, n, replace=False)) - fn(c, rng.choice(dB, n, replace=False))
                                   for _ in range(K_SUB)]))
    out = []
    for t in tissues:
        if t not in c.t:
            continue
        s = c.t[t]["subj"]
        a, b = np.intersect1d(dA, s), np.intersect1d(dB, s)
        n = min(len(a), len(b))
        if n < MIN_N:
            continue
        out.append(np.nanmedian([fn(c, t, rng.choice(a, n, replace=False)) - fn(c, t, rng.choice(b, n, replace=False))
                                 for _ in range(K_SUB)]))
    return float(np.nanmedian(out)) if out else np.nan


def perm_null(stat, lab, strata, rng, n):
    null = []
    for _ in range(n):
        new = lab.copy()
        for _, idx in lab.groupby(strata.reindex(lab.index)).groups.items():
            new.loc[idx] = rng.permutation(lab.loc[idx].to_numpy())
        null.append(stat(new))
    null = np.array(null)
    return null[np.isfinite(null)]


def p_two(obs, null):
    m = np.median(null)
    return float((1 + np.sum(abs(null - m) >= abs(obs - m))) / (1 + len(null)))


def timed_slope(x, col):
    f = timed_fit(x.b.to_numpy(), x.tod.to_numpy(), covariates(x), x[col].to_numpy(float))
    q = f["a1"] / f["a"]
    se = f["se_a1"] / abs(f["a"])
    return dict(log_amp=q.real, phase_h=q.imag / W, se_log_amp=se, se_phase_h=se / W, n=len(x))


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = 50 if run["quick"] else 1000
    rng = np.random.default_rng(94 + args.seed)
    seed = 94
    s = load_states(BODY + ["Pancreas"] + _brain_tissues(), args.subsample, args.seed)
    brain = sorted(t for t in s.tissue.unique() if t.startswith("Brain"))
    s = s.dropna(subset=["tod"])
    fast_s, vent_s = s[s.hardy.isin([1, 2])], s[s.hardy == 0]
    C = {"fast": Cohort(fast_s, True), "vent": Cohort(vent_s, False)}
    res, rows, labs, nulls = {}, [], {}, {}
    comp = s.drop_duplicates("subject").assign(dec=lambda d: d.age // 10 * 10).groupby(["hardy", "dec"]).size()
    comp.unstack(fill_value=0).to_csv(out / "a94_age_by_hardy.csv")
    for k, c in C.items():
        dn = c.donors.dropna(subset=["age"])
        labs[k] = tertiles(dn.age, rng)
        strata = dn.hardy.astype(str) + "|" + dn.sex.astype(str)
        stat = lambda lab, c=c: contrast(c, lab, 2, 0, radius, BODY, seed)
        o = stat(labs[k])
        nulls[k] = perm_null(stat, labs[k], strata, rng, n_perm)
        res[f"H94_1_{k}"], res[f"H94_1_{k}_p"] = o, p_two(o, nulls[k])
        res[f"age_tertile_cuts_{k}"] = ";".join(str(int(dn.age[labs[k] == g].min())) for g in range(3))
    m = min(len(nulls["fast"]), len(nulls["vent"]))
    o = (res["H94_1_fast"] + res["H94_1_vent"]) / 2
    res["H94_1"], res["H94_1_p"] = o, p_two(o, (nulls["fast"][:m] + nulls["vent"][:m]) / 2)
    c, dn = C["fast"], C["fast"].donors.dropna(subset=["age"])
    strata = dn.hardy.astype(str) + "|" + dn.sex.astype(str)
    st2 = lambda lab: contrast(c, lab, 2, 0, walltime_R, BODY, seed)
    o = st2(labs["fast"])
    nl = perm_null(st2, labs["fast"], strata, rng, n_perm)
    res["H94_2"], res["H94_2_p"] = o, p_two(o, nl)
    st3 = lambda lab: contrast(c, lab, 2, 0, lambda cc, d: coherence(cc, d, brain, BODY), None, seed, per_tissue=False)
    o = st3(labs["fast"])
    nl = perm_null(st3, labs["fast"], strata, rng, n_perm)
    res["H94_3"], res["H94_3_p_lower"] = o, float((1 + np.sum(nl <= o)) / (1 + len(nl)))
    for k, cc in C.items():
        lab = labs[k]
        for g in range(3):
            d = lab.index[lab == g].to_numpy()
            row = dict(cohort=k, group=f"age_t{g}", donors=len(d),
                       log_radius=np.nanmedian([radius(cc, t, d) for t in BODY if t in cc.t]),
                       C_body=coherence(cc, d, BODY))
            if k == "fast":
                row.update(R_wall=np.nanmedian([walltime_R(cc, t, d) for t in BODY if t in cc.t]),
                           C_brain_body=coherence(cc, d, brain, BODY),
                           C_hypothalamus_body=coherence(cc, d, ["Brain - Hypothalamus"], BODY))
            rows.append(row)
    ts = []
    for sel, x0 in (("fast", fast_s), ("hardy2", fast_s[fast_s.hardy == 2])):
        x0 = x0.dropna(subset=["age"]).assign(dec=lambda d: (d.age - d.age.mean()) / 10)
        for t, x in x0.groupby("tissue"):
            if len(x) >= (40 if not run["quick"] else 20):
                ts.append(dict(selection=sel, tissue=t, **timed_slope(x, "dec")))
    ts = pd.DataFrame(ts)
    ts.to_csv(out / "a94_timed_slopes.csv", index=False)
    for sel, x in ts.groupby("selection"):
        body = x[x.tissue.isin(BODY)]
        for col, se in (("log_amp", "se_log_amp"), ("phase_h", "se_phase_h")):
            r = cs.stats.meta(body[col], body[se])
            res[f"timed_{col}_per_decade_{sel}"], res[f"timed_{col}_se_{sel}"] = r["fixed"], r["se_fixed"]
        br = x[x.tissue.str.startswith("Brain")]
        if len(br):
            res[f"timed_log_amp_per_decade_brain_{sel}"] = cs.stats.meta(br.log_amp, br.se_log_amp)["fixed"]
    v = pd.Series(res)
    v["n_perm"] = n_perm
    v["H94_1_verdict"] = "PASS" if v.H94_1_p < 0.05 else "FAIL"
    v["H94_2_verdict"] = "PASS" if v.H94_2_p < 0.05 else "FAIL"
    v["H94_3_verdict"] = "PASS" if (v.H94_3 < 0 and v.H94_3_p_lower < 0.05) else "FAIL"
    v.to_csv(out / "a94_summary.csv", header=False)
    pd.DataFrame(rows).to_csv(out / "a94_groups.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).round(3).to_string(index=False), "\n", v.to_string())


def _brain_tissues():
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        return sorted(t for t in a.obs.tissue.astype(str).unique() if t.startswith("Brain"))
    finally:
        a.file.close()


if __name__ == "__main__":
    main()
