"""A103 (GTEx): is the brain-periphery split of the shared phase departure a collection or agonal artefact?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A103). Base analysis = `compute_organs_timed.py` (fast deaths, Hardy 1-2;
tissues with >= 40 donors; timed fit on recorded time; linear phase departures; brain vs periphery correlation of the mean
standardised departure disattenuated by split-half reliability, r*). Variants:
  - base: the 14-gene plane (sanity: r 0.57, r* 0.64);
  - H103.1 `adj`: per tissue, departures regressed on the sample's ischaemic time (h) and Hardy 2 vs 1; donor-level slope of
    (brain score - periphery score) on (brain - periphery mean ischaemic time, h) with bootstrap interval. PASS (split not
    explained) if r*_adj < 1 with bootstrap upper bound < 1;
  - H103.2 `noPER12`: clock plane without PER1 and PER2 (glucocorticoid / CREB-inducible within the hour), `cs.gtex().subset`;
  - D1: both, exact departures, amplitude departures; D2: Hardy 1 only, Hardy 2 only;
  - D3: immediate-early score (FOS, FOSB, JUNB, EGR1, NR4A1, ATF3; within-tissue z of log CPM) in brain and periphery vs the
    brain-specific departure (residual of the brain score on the periphery score);
  - D4 (Hardy 0) only when >= 4 brain tissues have >= 40 ventilated donors (not met in GTEx v10: not run).
Donor bootstrap of r* (200; reliabilities recomputed with 20 splits each). H103.3 (sheep) is in
`compute_brain_periphery_species.py`. Aggregates only (controlled access).

    python compute_brain_artefact.py [--subsample N] [--n-boot 200]
Outputs: paper/results/figure_2/output/organs/a103_gtex_{rstar,delay,ieg,isch_tissue}.csv
Legacy: legacy/desync/validation/brain_artefact/scripts/a103_gtex.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402

FOLDER = "organs"
IEG = ["FOS", "FOSB", "JUNB", "EGR1", "NR4A1", "ATF3"]
ACUTE = ("PER1", "PER2")
N_SPLIT = 50


def phase_departures(pb, t, genes, exact=False):
    """Phase and amplitude departures of every row from the timed fit of the clock plane built on `genes`."""
    wf = cs.gtex().subset(list(genes))
    ds = S.dataset(pb, waveform=wf)
    tissue = ds.data.compartment
    t_h = t.tod.to_numpy()
    lab, _, a, _, e = S.timed_fit(ds.data.beta, ds.data.V, t_h, tissue)
    return S.departures(e, a[np.searchsorted(lab, tissue)], t_h, exact)


def adjust(values, t: pd.DataFrame):
    """Per tissue OLS residuals on ischaemic time (h) and Hardy 2 vs 1 (when both classes are present)."""
    out = np.asarray(values, float).copy()
    for _, idx in t.groupby("compartment").indices.items():
        X = [np.ones(len(idx)), t.isch.to_numpy()[idx]]
        h2 = (t.hardy.to_numpy()[idx] == 2).astype(float)
        if 0 < h2.mean() < 1:
            X.append(h2)
        X = np.column_stack(X)
        ok = np.isfinite(X).all(1) & np.isfinite(out[idx])
        B, *_ = np.linalg.lstsq(X[ok], out[idx][ok], rcond=None)
        out[idx[ok]] = out[idx][ok] - X[ok] @ B + B[0]
        out[idx[~ok]] = np.nan
    return out


def block(Z, label, rng, n_boot):
    pt = S.block_stat(Z, rng, N_SPLIT)
    lo, hi = S.block_boot(Z, rng, n_boot)
    row = dict(analysis=label, donors=pt["donors"], r=pt["r"], rel_periph=pt["rel_periph"], rel_brain=pt["rel_brain"],
               r_star=pt["r_star"], r_star_lo=lo, r_star_hi=hi, n_boot=n_boot)
    row["pass"] = bool(np.isfinite(pt["r_star"]) and pt["r_star"] < 1 and hi < 1) if n_boot else np.nan
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    return row


def donor_scores(Z):
    isb = Z.columns.str.startswith("Brain")
    s = pd.concat([S._score(Z.loc[:, ~isb]).rename("periph"), S._score(Z.loc[:, isb]).rename("brain")], axis=1).dropna()
    B = np.c_[np.ones(len(s)), s.periph]
    coef, *_ = np.linalg.lstsq(B, s.brain, rcond=None)
    s["brain_specific"] = s.brain - B @ coef
    return s


def delay_slope(Z, t, label, rng, n_boot):
    """Slope of (brain - periphery score) on (brain - periphery mean ischaemic time, h) over donors."""
    s = donor_scores(Z)
    isch = t.assign(b=t.compartment.str.startswith("Brain")).groupby(["donor", "b"]).isch.mean().unstack()
    d = pd.DataFrame(dict(y=s.brain - s.periph, x=(isch[True] - isch[False]).reindex(s.index))).dropna()
    fit = lambda dd: np.polyfit(dd.x, dd.y, 1)[0]
    bs = [fit(d.iloc[rng.integers(0, len(d), len(d))]) for _ in range(n_boot)]
    row = dict(analysis=label, donors=len(d), slope_per_h=fit(d), slope_lo=np.quantile(bs, 0.025),
               slope_hi=np.quantile(bs, 0.975), r=d.corr().iloc[0, 1], delay_median_h=d.x.median(),
               delay_iqr_h=d.x.quantile(0.75) - d.x.quantile(0.25))
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    return row


def ieg_assoc(Z, t, ieg, rng, n_boot):
    """Correlation over donors of the brain-specific departure with brain, periphery and brain - periphery IEG score."""
    zi = ieg.groupby(t.compartment.to_numpy()).transform(lambda x: (x - x.mean()) / x.std()).mean(1)
    d = pd.DataFrame(dict(donor=t.donor.to_numpy(), b=t.compartment.str.startswith("Brain").to_numpy(), z=zi.to_numpy()))
    m = d.groupby(["donor", "b"]).z.mean().unstack()
    s = donor_scores(Z).join(m.rename(columns={True: "ieg_brain", False: "ieg_periph"}), how="inner").dropna()
    s["ieg_diff"] = s.ieg_brain - s.ieg_periph
    rows = []
    for y in ("brain_specific", "brain", "periph"):
        for x in ("ieg_brain", "ieg_periph", "ieg_diff"):
            r = s[[x, y]].corr().iloc[0, 1]
            bs = [s.iloc[i][[x, y]].corr().iloc[0, 1] for i in (rng.integers(0, len(s), len(s)) for _ in range(n_boot))]
            rows.append(dict(departure=y, ieg=x, donors=len(s), r=r, r_lo=np.quantile(bs, 0.025), r_hi=np.quantile(bs, 0.975)))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None, help="bootstrap replicates (default 200; 10 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = args.n_boot if args.n_boot is not None else (10 if run["quick"] else 200)
    rng = np.random.default_rng(0)
    pb, donors = S.cohort("fast", genes=list(S.CLOCK) + IEG, subsample=args.subsample, seed=args.seed)
    d = donors.reindex(pb.donor)
    t = pd.DataFrame({"donor": pb.donor.to_numpy(), "compartment": pb.compartment.to_numpy(), "tod": d.tod.to_numpy(),
                      "hardy": d.hardy.to_numpy(), "isch": S.sample_ischaemic_h(pb)})
    ieg = pd.DataFrame(np.log1p(pb[IEG].to_numpy(float) / pb.depth.to_numpy(float)[:, None] * 1e6), columns=IEG)
    genes14, genes12 = list(S.CLOCK), [g for g in S.CLOCK if g not in ACUTE]
    print(f"fast: {t.donor.nunique()} donors, {t.compartment.nunique()} tissues "
          f"({t.compartment.str.startswith('Brain').groupby(t.compartment).first().sum()} brain)", flush=True)
    zw = lambda v: S.zwide(v, t.donor.to_numpy(), t.compartment.to_numpy())

    rows, slopes = [], []
    ph14, am14 = phase_departures(pb, t, genes14)
    ph12, am12 = phase_departures(pb, t, genes12)
    Z = {"base": zw(ph14), "adj": zw(adjust(ph14, t)), "noPER12": zw(ph12), "adj+noPER12": zw(adjust(ph12, t))}
    rows.append(block(Z["base"], "base (01 sanity)", rng, n_boot))
    rows.append(block(Z["adj"], "PRIMARY H103.1: ischaemic time + Hardy adjusted", rng, n_boot))
    rows.append(block(Z["noPER12"], "PRIMARY H103.2: no PER1/PER2", rng, n_boot))
    rows.append(block(Z["adj+noPER12"], "D1: adjusted, no PER1/PER2", rng, n_boot))
    for lab, genes in (("exact", genes14), ("exact noPER12", genes12)):
        ph, _ = phase_departures(pb, t, genes, True)
        rows.append(block(zw(ph), f"D1: {lab}", rng, n_boot))
        rows.append(block(zw(adjust(ph, t)), f"D1: {lab}, adjusted", rng, n_boot))
    rows.append(block(zw(am14), "D1: amplitude departures", rng, n_boot))
    rows.append(block(zw(adjust(am14, t)), "D1: amplitude departures, adjusted", rng, n_boot))
    rows.append(block(zw(am12), "D1: amplitude departures, no PER1/PER2", rng, n_boot))
    for h in (1, 2):
        don = t.donor[t.hardy == h].unique()
        rows.append(block(Z["base"].loc[Z["base"].index.isin(don)], f"D2: Hardy {h} only", rng, n_boot))
        rows.append(block(Z["noPER12"].loc[Z["noPER12"].index.isin(don)], f"D2: Hardy {h} only, no PER1/PER2", rng, n_boot))
    slopes.append(delay_slope(Z["base"], t, "H103.1 donor slope: fast, 14 genes", rng, n_boot))
    slopes.append(delay_slope(Z["noPER12"], t, "fast, no PER1/PER2", rng, n_boot))
    ieg14 = ieg_assoc(Z["base"], t, ieg, rng, n_boot).assign(plane="14")
    ieg12 = ieg_assoc(Z["noPER12"], t, ieg, rng, n_boot).assign(plane="noPER12")
    pd.DataFrame(rows).to_csv(out / "a103_gtex_rstar.csv", index=False)
    pd.DataFrame(slopes).to_csv(out / "a103_gtex_delay.csv", index=False)
    pd.concat([ieg14, ieg12]).to_csv(out / "a103_gtex_ieg.csv", index=False)
    eff = []
    for tis, g in t.assign(ph=ph14).groupby("compartment"):
        ok = np.isfinite(g.ph) & np.isfinite(g.isch)
        eff.append(dict(tissue=tis, brain=tis.startswith("Brain"), n=int(ok.sum()),
                        slope_h_per_h=np.polyfit(g.isch[ok], g.ph[ok], 1)[0]))
    pd.DataFrame(eff).to_csv(out / "a103_gtex_isch_tissue.csv", index=False)


if __name__ == "__main__":
    main()
