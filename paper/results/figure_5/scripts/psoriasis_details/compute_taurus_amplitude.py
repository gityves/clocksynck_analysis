"""A53 H53.1: tissue-wide clock amplitude of inflamed vs non-inflamed gut biopsies of the same endoscopy (TAURUS).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A53 H53.1). Builder `taurus` (Thomas et al. 2024, Zenodo 13768607): unit =
patient x timepoint with >= 1 A45-eligible inflamed and non-inflamed sample (T cells and >= 2 of fibroblast, endothelial,
mural; the sample with most anchor cells is kept per status). Every class x status is a compartment of the unit (`I:<class>`,
`N:<class>`), so the amplitudes of both biopsies are read on one free host (`Dataset(anchors=None)`). Estimand: the mean over
classes present in both statuses in >= 10 units of log|z_I:k| - log|z_N:k|; H53.1: < 0, one-sided p < 0.05 against statuses
swapped per unit (n_perm) and below the 5th percentile of 200 expression-matched pseudo-clock sets (all-gene sample x class
pseudobulks, `Paired` building block `pseudo_waveforms`); the gate must pass. Descriptive: per-class log amplitude ratio,
coupling difference (leave-one-out host) and offset difference.

    python compute_taurus_amplitude.py [--subsample N] [--n-perm 2000] [--n-pseudo 200] [--jobs 4]
Outputs (paper/results/figure_5/output/psoriasis_details/): a53_taurus_{counts,gate,classes,nulls,verdict}.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a53_taurus.py (+ ibd_taurus/scripts/a45_taurus.py loaders)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import _details as D  # noqa: E402
import _shared as H  # noqa: E402
import _singlecell as C  # noqa: E402
import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA  # noqa: E402

SEED, MIN_UNITS = 53, 10
TCELL = "gut:T"
ANCHORS = ["gut:fibroblast", "gut:endothelial", "gut:mural"]
ALL_GENES = DATA / "ibd" / "TAURUS" / "pseudobulk" / "pseudobulk_all.h5ad"


def units(pb, sm) -> pd.DataFrame:
    """Status-split table: donor = patient|timepoint, compartment = I:/N:<class>; one sample per status."""
    sets = pb.groupby("donor").compartment.apply(set)
    ok = sets[sets.apply(lambda x: TCELL in x and len(x & set(ANCHORS)) >= 2)].index
    s = sm[sm.index.isin(ok) & sm.arm.notna()]
    both = s.groupby("patient").arm.nunique()
    s = s[s.patient.isin(both.index[both == 2])].copy()
    s["unit"] = s.patient.astype(str) + "|" + s.timepoint.astype(str)
    anc = pb[pb.compartment.isin(ANCHORS)].groupby("donor").n_cells.sum()
    s["anchor_cells"] = anc.reindex(s.index).fillna(0)
    both = s.groupby("unit").arm.nunique()
    s = s[s.unit.isin(both.index[both == 2])]
    pick = s.sort_values("anchor_cells", ascending=False).groupby(["unit", "arm"]).head(1)
    t = pb[pb.donor.isin(pick.index)].copy()
    t["status"] = t.donor.map(pick.arm).map({"case": "I", "control": "N"})
    t["sample"] = t.donor
    t["donor"] = t.donor.map(pick.unit)
    t["compartment"] = t.status + ":" + t.compartment.str.split(":").str[1]
    return t.drop(columns=["stratum"]).reset_index(drop=True)


def classes(t) -> list[str]:
    n = t.groupby("compartment").donor.nunique()
    return [c for c in sorted({c[2:] for c in n.index}) if n.get(f"I:{c}", 0) >= MIN_UNITS and n.get(f"N:{c}", 0) >= MIN_UNITS]


def stat(t, ks, waveform="human", min_donors=MIN_UNITS) -> dict:
    t = t[t.compartment.str[2:].isin(ks)].reset_index(drop=True)
    tab = cs.Dataset(t, anchors=None, waveform=waveform, min_donors=min_donors).fit(n_boot=0, gate=False).table
    d = {k: float(tab.log_amp_rel[f"I:{k}"] - tab.log_amp_rel[f"N:{k}"]) for k in ks}
    c = {k: float(tab.coupling[f"I:{k}"] - tab.coupling[f"N:{k}"]) for k in ks}
    o = {k: float(C.wrap(tab.offset_h[f"I:{k}"] - tab.offset_h[f"N:{k}"])) for k in ks}
    return dict(stat=float(np.mean(list(d.values()))), log_amp=d, d_coupling=c, d_offset_h=o)


def swapped(t, rng):
    u = t.donor.unique()
    flip = set(u[rng.random(len(u)) < 0.5])
    s = t.copy()
    f = s.donor.isin(flip)
    s.loc[f, "compartment"] = np.where(s.loc[f, "compartment"].str.startswith("I:"), "N:" + s.loc[f, "compartment"].str[2:],
                                       "I:" + s.loc[f, "compartment"].str[2:])
    return s


def pseudo_tables(t, n_sets, rng):
    """Status-split tables over pseudo-clock genes (counts from the all-gene pseudobulks)."""
    A = ad.read_h5ad(ALL_GENES)
    key = A.obs.donor_id.astype(str) + "|" + A.obs.cell_type.astype(str)
    row = pd.Series(np.arange(A.n_obs), index=key.to_numpy())
    want = (t["sample"].str.replace("taurus:", "", regex=False) + "|" + t.compartment.str[2:]).to_numpy()
    allc = pd.DataFrame(np.asarray(A.X[row.loc[want].to_numpy()].todense() if hasattr(A.X, "todense")
                                   else A.X[row.loc[want].to_numpy()]), columns=A.var_names)
    ref = t.compartment.str.startswith("N:").to_numpy()
    lm = np.log(allc[ref].div(t.depth[ref].to_numpy(), axis=0).mean() + 1e-12)
    lm = lm[allc[ref].mean() >= 1]
    out = []
    for wf in H.PD.pseudo_waveforms(lm, n_sets, rng, H.WF):
        tp = t[["donor", "compartment", "depth", "n_cells"]].copy()
        tp[list(wf.genes)] = allc[list(wf.genes)].to_numpy()
        out.append((tp, wf))
    return out


def main(argv=None):
    ap = H.parser(__doc__)
    ap.add_argument("--n-pseudo", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    n_perm = args.n_perm or (30 if run["quick"] else 2000)
    n_pseudo = args.n_pseudo or (10 if run["quick"] else 200)
    pb, sm = C.load("taurus", None)  # builder: samples are the donors
    if args.subsample:
        pat = np.random.default_rng(args.seed).permutation(sorted(sm.patient.unique()))[:max(args.subsample // 2, 8)]
        sm = sm[sm.patient.isin(pat)]
        pb = pb[pb.donor.isin(sm.index)]
    t = units(pb, sm)
    ks = classes(t)
    t = t[t.compartment.str[2:].isin(ks)].reset_index(drop=True)
    t.groupby("compartment").donor.nunique().to_csv(out / "a53_taurus_counts.csv")
    gate = C.gate(t, SEED, min_donors=MIN_UNITS)
    pd.Series(gate).to_csv(out / "a53_taurus_gate.csv", header=False)
    obs = stat(t, ks)
    print("units", t.donor.nunique(), "classes", ks, "\n", gate, "\nstat", obs["stat"], flush=True)
    rng = np.random.default_rng(SEED)
    par = Parallel(n_jobs=args.jobs)
    null = np.array([x if x is not None else np.nan for x in par(
        delayed(C.safe)(lambda s: stat(swapped(t, np.random.default_rng(s)), ks)["stat"], s)
        for s in rng.integers(0, 2**31, n_perm))], float)
    ps = pseudo_tables(t, n_pseudo, rng)
    pnull = np.array([x if x is not None else np.nan for x in par(
        delayed(C.safe)(lambda tp, wf: stat(tp, ks, wf)["stat"], tp, wf) for tp, wf in ps)], float)
    pd.DataFrame(dict(swap=pd.Series(null), pseudo=pd.Series(pnull))).to_csv(out / "a53_taurus_nulls.csv", index=False)
    per = pd.DataFrame(dict(log_amp=obs["log_amp"], d_coupling=obs["d_coupling"], d_offset_h=obs["d_offset_h"]))
    per.index.name = "class"
    per.to_csv(out / "a53_taurus_classes.csv")
    p_swap, p_ps = C.p_lower(null, obs["stat"]), C.p_lower(pnull, obs["stat"])
    q05 = float(np.nanquantile(pnull, 0.05))
    v = dict(units=t.donor.nunique(), classes=";".join(ks), gate_pass=gate["gate_pass"], stat=obs["stat"], p_swap=p_swap,
             n_swap=int(np.isfinite(null).sum()), pseudo_q05=q05, p_pseudo=p_ps, n_pseudo=int(np.isfinite(pnull).sum()),
             stromal_vascular=float(np.mean([obs["log_amp"][k] for k in ks if k in ("fibroblast", "endothelial", "mural")])),
             immune=float(np.mean([obs["log_amp"][k] for k in ks if k in ("T", "B_plasma", "myeloid")])),
             H53_1=("UNINFORMATIVE (gate)" if not gate["gate_pass"] else
                    "PASS" if (obs["stat"] < 0 and p_swap < 0.05 and obs["stat"] < q05) else "FAIL"))
    pd.Series(v).to_csv(out / "a53_taurus_verdict.csv", header=False)
    print(per.round(3).to_string(), "\n", pd.Series(v).to_string())


if __name__ == "__main__":
    main()
