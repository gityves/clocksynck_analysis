"""Does extending the 14-gene waveform with clock-target genes make the untimed phase estimate more precise? (A28)

Pre-registered as A28 (2026-10-15, written before the evaluation; run note: the primary fails, the 14-gene waveform
stands). GTEx fast deaths (Hardy 1-2) are split 50/50 by subject (seed 0): selection half (every waveform fit), evaluation
half. Waveforms fitted on the selection half with `waveform.fit_waveform` (15 tissues, tissue as group, y = log((count +
0.5) / library size)): W14 (the clock genes, refitted), W18 (+ BHLHE41, HLF, PPL, ADIRF-AS1: primary), Wcons (+ 14 genes),
Wgen (+ 19 genes); the original `gtex()` waveform for reference. The gene sets are those fixed in the amendment (the
screen that produced them is not rerun).
  E1 (primary) single sample, single tissue: per evaluation sample the phase of its clock-plane state (projector per
     tissue x half), error = phase - 2 pi t_death / 24 after removing a per-tissue offset estimated on the selection half;
     score = median |error| pooled over samples and per tissue; donor-bootstrap CI (1000) of the change against W14,
     tissues improved of 15. Pass: CI of the change below 0 and >= 9 of 15 tissues improve; useful: gain >= 10 %.
  E2 multi-tissue host (`Dataset.timed`): evaluation donors with all 15 tissues and with three tissues (Muscle, Heart AA,
     Lung): R of host phase vs time of death and median error after the common offset.
  E3 ventilated donors (Hardy 0): host (15 tissues) vs CHIRAL DIP phase (`data/chiral_DIP.csv`) and vs time of death
     (reference); no harm: R drops by <= 0.05. Not independent for W18 / Wgen.
E4 (Li 2024 breast, A25 estimator) is not ported: the host of two compartments needs `Dataset.fit` with 2 anchors.
Controlled access: times of death stay in memory; the output has pooled and per-tissue medians, correlations, no donor rows.

    python compute_extended_clock.py [--subsample N] [--seed S] [--jobs 4] [--n-boot N]
Outputs (paper/results/extended_data/output/gate/): extended_clock_{waveforms,e1_summary,e1_tissues,e2_e3}.csv
Legacy: legacy/desync/validation/extended_clock/scripts/evaluate.py -> output/{waveforms,e1_summary,e1_tissues,e2_e3}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from clocksync.waveform import fit_waveform  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import GTEX_PHENOTYPES, out_dir  # noqa: E402

import _gate_shared as sh  # noqa: E402

FOLDER = "gate"
W = sh.W
CORE = list(cs.gtex().genes)
STRICT = ["BHLHE41", "HLF", "PPL", "ADIRF-AS1"]
CONS = STRICT + ["CFAP70", "CYS1", "MTHFD1L", "MUC20", "PLXDC1", "POMZP3", "RHOBTB3", "TXNIP", "VWF", "WASH8P"]
GEN = ["ADIRF-AS1", "BHLHE40", "BHLHE41", "CAMKK1", "CHST7", "FAM76A", "HCP5", "HLF", "IL1R1", "LAMA2", "LFNG", "LONRF1",
       "MAMDC2", "NFIL3", "NIPSNAP1", "PIK3IP1", "PPL", "RAB3IL1", "SLC43A3"]
ARMS = {"W14": [], "W18": STRICT, "Wcons": CONS, "Wgen": GEN}
TIS = ("Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)",
       "Esophagus - Mucosa", "Nerve - Tibial", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
       "Adipose - Subcutaneous", "Lung", "Pituitary", "Brain - Cerebellum", "Brain - Cortex", "Whole Blood")
THREE = ["Muscle - Skeletal", "Heart - Atrial Appendage", "Lung"]
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import PAPER_DATA  # noqa: E402
DIP = PAPER_DATA / "chiral_DIP.csv"


def cd(a):
    return np.angle(np.exp(1j * a))


def load(args):
    """Rows of the 15 tissues, Hardy 0-2, library-size depth, clock genes + the extra genes; `half` = sel / eval / vent."""
    genes = CORE + sorted(set(CONS) | set(GEN))
    pb, donors = gtex.load(tissues=TIS, hardy=(0, 1, 2), genes=genes, depth="library_size", require_tod=False)
    allfast, _ = gtex.load(hardy=(1, 2), genes=CORE[:1], require_tod=False)  # all fast-death subjects, any tissue
    subj = np.array(sorted(allfast.donor.unique()))
    sel = set(np.random.default_rng(0).choice(subj, len(subj) // 2, replace=False))
    pb = pb.merge(donors[["hardy", "tod"]], left_on="donor", right_index=True)
    pb["half"] = np.where(pb.hardy == 0, "vent", np.where(pb.donor.isin(sel), "sel", "eval"))
    if args.subsample is not None:
        keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[:args.subsample]
        pb = pb[pb.donor.isin(keep)]
    return pb.reset_index(drop=True)


def waveforms(df):
    s = df[(df.half == "sel") & df.tod.notna()]
    out = {}
    for arm, extra in ARMS.items():
        g = CORE + extra
        y = np.log((s[g].to_numpy(float) + 0.5) / s.depth.to_numpy(float)[:, None])
        out[arm] = fit_waveform(y, s.tod.to_numpy(float) % 24, g, group=s.compartment.to_numpy(), name=arm)
    out["gtex_orig"] = cs.gtex()
    return out


def e1(df, wf):
    """local: not in clocksync (single-sample phase error). Per tissue x half one projector, per tissue the offset and
    rotation from the selection half (timed), error on the evaluation half."""
    d = df[df.half.isin(["sel", "eval"]) & df.tod.notna()].reset_index(drop=True)
    g = list(wf.genes)
    beta = np.zeros((len(d), 2))
    for (_, _), i in d.groupby(["compartment", "half"]).indices.items():
        sub = d.iloc[i]
        L = projector(sub[g], wf)
        beta[i], _ = states(sub[g], L, wf, depth=sub.depth.to_numpy(float))
    wt = W * (d.tod.to_numpy(float) % 24)
    err = np.full(len(d), np.nan)
    for t in TIS:
        s = np.flatnonzero(((d.compartment == t) & (d.half == "sel")).to_numpy())
        e = np.flatnonzero(((d.compartment == t) & (d.half == "eval")).to_numpy())
        if len(s) < 10 or not len(e):
            continue
        D = np.c_[np.ones(len(s)), np.cos(wt[s]), np.sin(wt[s])]
        B, *_ = np.linalg.lstsq(D, beta[s], rcond=None)
        b0 = B[0]
        ph_s = np.angle((beta[s] - b0) @ [1, 1j])
        rot = np.angle(np.exp(1j * (ph_s - wt[s])).mean())
        ph_e = np.angle((beta[e] - b0) @ [1, 1j])
        err[e] = cd(ph_e - wt[e] - rot)
    d["err_h"] = err / W
    return d.loc[(d.half == "eval") & d.err_h.notna(), ["donor", "compartment", "err_h"]]


def host_dataset(df, wf, tissues, half):
    d = df[(df.half == half) & df.compartment.isin(tissues)]
    g = ["donor", "compartment", "depth"] + list(wf.genes)
    return cs.Dataset(d[g].reset_index(drop=True), anchors=list(tissues), waveform=wf, name=f"{half}")


def e2e3(df, wf, arm):
    limit_threads(1)
    rows = []
    tod = df.drop_duplicates("donor").set_index("donor").tod
    for lab, tissues in (("15 tissues", list(TIS)), ("3 tissues", THREE)):
        ds = host_dataset(df, wf, tissues, "eval")
        t = tod.reindex(ds.donor_ids).dropna() % 24
        r = ds.timed(t, n_perm=200)
        rows.append(dict(arm=arm, readout="E2", host=lab, n=int(r.n), R=r.R, offset_h=r.offset_h,
                         med_abs_err_h=r.median_abs_err_h))
    ds = host_dataset(df, wf, list(TIS), "vent")
    dip = pd.read_csv(DIP, sep=";", decimal=",")
    dip.index = "GTEX-" + dip["SUBJ.ID"].astype(str)
    for lab, t in (("E3 host vs DIP", dip.hour.reindex(ds.donor_ids).dropna() % 24),
                   ("E3 host vs time of death (ref.)", tod.reindex(ds.donor_ids).dropna() % 24)):
        r = ds.timed(t, n_perm=200)
        rows.append(dict(arm=arm, readout=lab, host="15 tissues", n=int(r.n), R=r.R, offset_h=r.offset_h,
                         med_abs_err_h=r.median_abs_err_h))
    return rows


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None, help="E1 bootstrap draws (default 1000; 100 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    nb = args.n_boot if args.n_boot is not None else (100 if run["quick"] else 1000)
    df = load(args)
    print(f"{df.donor.nunique()} donors; halves: {df.drop_duplicates('donor').half.value_counts().to_dict()}", flush=True)
    wfs = waveforms(df)
    pd.concat([w.table().assign(arm=k).reset_index(names="gene") for k, w in wfs.items()]).to_csv(
        out / "extended_clock_waveforms.csv", index=False)
    errs = {k: e1(df, w) for k, w in wfs.items()}
    e = pd.concat({k: v.set_index(["donor", "compartment"]).err_h for k, v in errs.items()}, axis=1).dropna()
    ab = e.abs()
    rng = np.random.default_rng(args.seed)
    donors = ab.index.get_level_values(0).unique().to_numpy()
    by = {dn: ab.xs(dn, level=0).to_numpy() for dn in donors}
    boots = pd.DataFrame([np.median(np.vstack([by[dn] for dn in rng.choice(donors, len(donors))]), axis=0)
                          for _ in range(nb)], columns=ab.columns)
    pt = ab.groupby(level=1).median()
    rows = []
    for k in ab.columns:
        dlt = boots[k] - boots["W14"]
        rows.append(dict(arm=k, samples=len(ab), tissues=pt.shape[0], med_abs_err_h=ab[k].median(),
                         change_h=ab[k].median() - ab["W14"].median(),
                         change_pct=100 * (ab[k].median() / ab["W14"].median() - 1),
                         ci_lo=dlt.quantile(0.025), ci_hi=dlt.quantile(0.975),
                         tissues_better=int((pt[k] < pt["W14"]).sum()), within2h=(ab[k] <= 2).mean()))
    s = pd.DataFrame(rows)
    s["primary_pass"] = (s.ci_hi < 0) & (s.tissues_better >= 9) & s.arm.isin(["W18"])
    s.to_csv(out / "extended_clock_e1_summary.csv", index=False)
    pt.to_csv(out / "extended_clock_e1_tissues.csv")
    res = Parallel(n_jobs=args.jobs)(delayed(e2e3)(df, w, k) for k, w in wfs.items())
    e23 = pd.DataFrame([r for rr in res for r in rr])
    e23.to_csv(out / "extended_clock_e2_e3.csv", index=False)
    pd.set_option("display.width", 220)
    print(s.round(3).to_string(index=False))
    print(e23.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
