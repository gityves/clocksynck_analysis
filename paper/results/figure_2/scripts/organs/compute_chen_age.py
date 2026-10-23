"""A3: timed test of the age effect on the prefrontal clock (Chen et al. 2016, GSE71620, BA11 and BA47, time of death known).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A3; `desync/05_brain_age_timed`). Per region a sample's clock-plane state
is the projection of its clock-gene log expression (natural log, equal gene weights, microarray: no counting-noise model)
on the plane of the GTEx waveform (`clocksync.paired.projector / states`). Then by complex least squares,
beta = mu + d x + (a + c x) e^{it} (t = time of death, ZT):
  - primary: log2 amplitude ratio log2|a + c| / |a| of donors > 60 y over donors < 40 y (donors in between excluded), regions
    pooled by inverse-variance weights from the bootstrap SEs, donor bootstrap (both regions of a donor together);
  - secondaries: age per decade (log2 amplitude change and phase advance, h; the GTEx covariate model), the SD of the phase
    departure from time of death (> 60 over < 40), the phase difference old minus young;
  - sensitivity (`--adjust ph rin pmi`): agonal-state and RNA-quality proxies as mean and amplitude modifiers;
  - per gene cosinor amplitude (log2) in young and old, per region.
The pairing of BA11 and BA47 samples to donors is rebuilt from (age, sex, race, PMI, pH, time of death), as in the legacy.

    python compute_chen_age.py [--subsample N] [--adjust ph rin pmi] [--n-boot 500]
Outputs: paper/results/figure_2/output/organs/chen_age[_adj-...]_{region,pooled}.csv, chen_age_genes.csv
Legacy: legacy/desync/05_brain_age_timed/scripts/chen_age.py (data: $SCCC_DATA/brain_tod/GSE71620/, legacy prepare_gse71620.py)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

FOLDER = "organs"
RAW = DATA / "brain_tod" / "GSE71620"
W = 2 * np.pi / 24
KEY = ["age", "sex", "race", "pmi", "ph", "tod_raw"]
YOUNG, OLD = 40, 60
COLS = ["log2_amp_old_vs_young", "phase_old_minus_young_h", "log2_dep_sd_old_vs_young", "log2_amp_per_decade",
        "phase_advance_per_decade_h"]


def load():
    """Samples with a time of death (donor id rebuilt), log2 expression of the waveform genes, the waveform."""
    s = pd.read_csv(RAW / "samples.csv")
    s = s[s.tod_raw.notna()].copy()
    s["donor"] = s[KEY].astype(str).agg("|".join, axis=1)
    n = s.groupby("donor").region.nunique()
    if (n != 2).any() or s.groupby(["donor", "region"]).size().max() > 1:
        raise ValueError("BA11/BA47 pairing is not one-to-one")
    s["donor"] = s.donor.map({d: f"d{i:03d}" for i, d in enumerate(sorted(s.donor.unique()))})
    e = pd.read_csv(RAW / "expression_clock.csv").set_index("gsm").loc[s.gsm]
    e.columns = [canonical(c.split("|")[0]) for c in e.columns]
    wf = cs.gtex()
    genes = [g for g in wf.genes if g in e.columns]
    y = e[genes].T.groupby(level=0).mean().T[genes]
    return s.reset_index(drop=True), y.reset_index(drop=True), wf.subset(genes)


def clock_states(s, y, wf) -> np.ndarray:
    """Complex clock-plane state per sample: one projector (equal gene weights) on the log expression centred within
    region."""
    yc = y - y.groupby(s.region.to_numpy()).transform("mean")
    L = projector(yc, wf, counts=False)
    b, _ = states(yc, L, wf)  # log2 -> natural log inside `states`
    return b[:, 0] + 1j * b[:, 1]


def timed(beta, t, x, z=None):
    """local: not in clocksync. beta = mu + d x + (a + c x) e^{it} (+ covariates z as mean and amplitude modifiers) by
    complex least squares; returns a, c and the residuals."""
    e = np.exp(1j * t * W)
    z = np.empty((len(beta), 0)) if z is None else z
    D = np.column_stack([np.ones(len(beta)), x, e, x * e, z, z * e[:, None]])
    coef = np.linalg.lstsq(D, beta, rcond=None)[0]
    return coef[2], coef[3], beta - D @ coef


def measures(s, beta, adjust=()) -> pd.DataFrame:
    rows = []
    for r in ("BA11", "BA47"):
        i = (s.region == r).to_numpy()
        age, t = s.age.to_numpy(float)[i], s.tod_raw.to_numpy(float)[i]
        g = np.where(age > OLD, 1.0, np.where(age < YOUNG, 0.0, np.nan))
        z = s.loc[i, list(adjust)].to_numpy(float) if adjust else None
        if z is not None:
            z = np.where(np.isnan(z), np.nanmean(z, 0), z)
            z = (z - z.mean(0)) / z.std(0)
        k = ~np.isnan(g)
        a, c, res = timed(beta[i][k], t[k], g[k], None if z is None else z[k] - z[k].mean(0))
        u = res / ((a + c * g[k]) * np.exp(1j * t[k] * W))
        dep = np.angle(1 + u) / W
        a2, c2, _ = timed(beta[i], t, (age - age.mean()) / 10, z)
        rows.append(dict(region=r, log2_amp_old_vs_young=np.log2(abs(a + c) / abs(a)),
                         phase_old_minus_young_h=np.angle((a + c) / a) / W,
                         log2_dep_sd_old_vs_young=np.log2(dep[g[k] == 1].std() / dep[g[k] == 0].std()),
                         log2_amp_per_decade=(c2 / a2).real / np.log(2), phase_advance_per_decade_h=(c2 / a2).imag / W,
                         amp_young=abs(a), amp_old=abs(a + c), n_young=int((g == 0).sum()), n_old=int((g == 1).sum()),
                         n_all=int(i.sum())))
    return pd.DataFrame(rows).set_index("region")


def gene_cosinor(s, y) -> pd.DataFrame:
    rows = []
    for r in ("BA11", "BA47"):
        for name, sel in (("young", s.age < YOUNG), ("old", s.age > OLD)):
            i = ((s.region == r) & sel).to_numpy()
            t = s.tod_raw.to_numpy(float)[i]
            for gname in y.columns:
                c = cs.timed.cosinor(y[gname].to_numpy()[i], t)
                rows.append(dict(region=r, group=name, gene=gname, amp_log2=float(c["amplitude"]),
                                 peak_zt=float(c["acrophase_h"]), n=int(i.sum())))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--adjust", nargs="*", default=[], help="sample covariates to adjust for (ph, rin, pmi)")
    ap.add_argument("--n-boot", type=int, default=None, help="donor bootstrap replicates (default 500; 40 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = args.n_boot if args.n_boot is not None else (40 if run["quick"] else 500)
    rng = np.random.default_rng(args.seed)
    s, y, wf = load()
    if args.subsample:
        keep = rng.permutation(sorted(s.donor.unique()))[:args.subsample]
        m = s.donor.isin(keep).to_numpy()
        s, y = s[m].reset_index(drop=True), y[m].reset_index(drop=True)
    beta = clock_states(s, y, wf)
    tag = "_adj-" + "-".join(args.adjust) if args.adjust else ""
    point = measures(s, beta, args.adjust)
    donors = s.donor.unique()
    rows = {d: np.flatnonzero(s.donor.to_numpy() == d) for d in donors}
    reps = []
    for _ in range(n_boot):
        idx = np.concatenate([rows[d] for d in rng.choice(donors, len(donors))])
        reps.append(measures(s.iloc[idx].reset_index(drop=True), beta[idx], args.adjust))
    reps = pd.concat(reps, keys=range(n_boot), names=["rep"])
    for c in COLS:
        q = reps[c].groupby(level="region").quantile([0.025, 0.975]).unstack()
        point[f"{c}_lo"], point[f"{c}_hi"] = q[0.025], q[0.975]
        point[f"{c}_se"] = reps[c].groupby(level="region").std()
    point.to_csv(out / f"chen_age{tag}_region.csv")
    pooled = []
    for c in COLS:
        w = 1 / point[f"{c}_se"] ** 2
        w = w / w.sum()
        b = (reps[c].unstack("region") * w).sum(axis=1)
        est = (point[c] * w).sum()
        pooled.append(dict(measure=c, pooled=est, lo=b.quantile(0.025), hi=b.quantile(0.975), se=b.std(),
                           p_two_sided=min(1, 2 * min((b <= 0).mean(), (b >= 0).mean()) + 1 / n_boot)))
    pooled = pd.DataFrame(pooled)
    pooled.to_csv(out / f"chen_age{tag}_pooled.csv", index=False)
    gene_cosinor(s, y).to_csv(out / "chen_age_genes.csv", index=False)
    print(point[COLS + ["amp_young", "amp_old", "n_young", "n_old"]].round(3).T.to_string())
    print(pooled.round(3).to_string())


if __name__ == "__main__":
    main()
