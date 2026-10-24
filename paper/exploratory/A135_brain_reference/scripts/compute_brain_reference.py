"""A135 part 1: which brain compartment family, if any, should carry the host? Control donors only.

Registered in paper/docs/registered_experiments.md (A135). No disease contrast is computed here.
Per dataset (control donors) and compartment:
  1. clock strength: clock-plane share of the compartment's own Poisson-corrected covariance (rows = donor x
     compartment, centred within stratum; as A86) against 1,000 permuted-acrophase planes; excess = share - null mean.
  2. gate: share of the other compartments with which the model-free pairwise cross gate (`Dataset.gate` on the two-compartment
     subset; p_pairing and p_acrophase < 0.05) passes.
  3. centrality: mean pairwise C (|complex correlation| of the Poisson-corrected clock states over shared donors,
     `Dataset.coherence`, the A109 C) with the other compartments.
Rule (A135): rank by the rank sum of 1-3 (rank 1 = best); a family is the reference if it ranks first in >= 2/3 of the
datasets and passes the gate with >= 80 % of the others in each of those datasets; ties -> larger median centrality.

Choices not fixed by the registration (stated, not tuned): datasets that share donors are collapsed to one
analysis on the largest control set (nested control sets: psychad_mssm = SCZ_mssm controls, 257 donors, a superset of
the 163 of AD/DLBD/PD/VaD_mssm; psychad_hbcc, psychad_radc = DM_radc, 26, a superset of AD_radc, 25); 5 independent
datasets, so a reference must rank first in >= 4. Compartments need >= 15 control donors; pairs >= 15 shared donors
(pd has 24 controls). Within a dataset, equal rank sums are broken by centrality. `pd`: regions pooled by cell type for the
family map, compartments stay region-specific. 200 gate permutations, 100 coherence permutations.

    python compute_brain_reference.py [--subsample N] [--seed S] [--jobs 4]
Outputs (output/, subsample: output/subsample/): a135_criteria.csv, a135_ranks.csv, a135_decision.csv,
a135_family_map.csv, a135_datasets.csv
"""

from __future__ import annotations

import itertools
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "clocksync" / "src"))
sys.path.insert(0, str(ROOT / "paper"))
import clocksync as cs  # noqa: E402
from clocksync.decompose import logfrac  # noqa: E402
from common import cli, data  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
WF = cs.gtex()
GENES = list(WF.genes)
MIN_DONORS, MIN_SHARED, ALPHA = 15, 15, 0.05
# dataset -> (control label columns of the contrasts it serves, the contrast whose controls (largest set) are used)
DATASETS = {
    "psychad_mssm": (["AD_mssm", "DLBD_mssm", "PD_mssm", "SCZ_mssm", "VaD_mssm"], "SCZ_mssm"),
    "psychad_hbcc": (["ASCVD_hbcc", "BD_hbcc", "SCZ_hbcc"], "ASCVD_hbcc"),
    "psychad_radc": (["AD_radc", "DM_radc"], "DM_radc"),
    "seaad_dlpfc": (["AD_seaad_dlpfc"], "AD_seaad_dlpfc"),
    "parkinson": (["pd"], "pd"),
}
REGION = re.compile(r"^(DMV|GPi|M1|PFC|V1)_")


def family(comp: str) -> str:
    """Compartment -> family (cell type; regions of `pd` pooled; GABAergic = inhibitory, glutamatergic = excitatory)."""
    c = REGION.sub("", comp.removeprefix("brain:"))
    if c.startswith("inhibitory") or c == "GABAergic":
        return "inhibitory"
    return {"excitatory": "excitatory", "glutamatergic": "excitatory", "neuron_other": "neuron_other"}.get(c, c)


def plane(w) -> np.ndarray:
    Xc = w.X - w.X.mean(0)
    return Xc @ np.linalg.solve(Xc.T @ Xc, Xc.T)


def share_cov(S, P) -> float:
    return float(np.sum((P @ S @ P) ** 2) / np.sum(S ** 2))


def clock_strength(x: pd.DataFrame, n_perm: int, rng) -> dict:
    """Poisson-corrected clock-plane share of the rows of one compartment against permuted-acrophase planes (A86)."""
    y, pv, k = logfrac(x, GENES, x.stratum.astype(str).to_numpy())
    S = y.T @ y / (len(y) - k) - np.diag(pv)
    obs = share_cov(S, plane(WF))
    null = np.array([share_cov(S, plane(WF.permuted(rng))) for _ in range(n_perm)])
    return dict(share=obs, share_null=float(null.mean()), excess=obs - float(null.mean()),
                p_share=float((1 + np.sum(null >= obs)) / (1 + n_perm)))


def gate_row(ds, k, partners, n_perm, seed, ms):
    warnings.simplefilter("ignore")
    out = []
    for l in partners:
        try:
            g = ds.subset(compartments=[k, l], anchors=[k, l]).gate(n_perm=n_perm, seed=seed, min_shared=ms)
            out.append((k, l, bool(g.passes), g.p_pairing, g.p_acrophase))
        except ValueError:  # fewer than MIN_SHARED shared donors: pair not testable
            out.append((k, l, np.nan, np.nan, np.nan))
    return out


def one_dataset(name, run, subsample, seed, jobs):
    warnings.simplefilter("ignore")
    contrasts, ref = DATASETS[name]
    pb, donors = data.load(name)
    ctrl = sorted(donors.index[donors[ref] == "control"])
    if subsample:
        ctrl = sorted(np.random.default_rng(seed).permutation(ctrl)[:subsample])
    pb = pb[pb.donor.isin(ctrl)].reset_index(drop=True)
    n_d = pb.groupby("compartment").donor.nunique()
    comps = sorted(n_d.index[n_d >= (5 if subsample else MIN_DONORS)])
    pb = pb[pb.compartment.isin(comps)].reset_index(drop=True)
    ds = cs.Dataset(pb, anchors=comps, name=name)
    comps = ds.compartments
    rng = np.random.default_rng(135 + seed)
    rows = []
    for k in comps:
        rows.append(dict(dataset=name, compartment=k, family=family(k), n_control_donors=int(n_d[k]),
                         **clock_strength(pb[pb.compartment == k], 50 if run["quick"] else 1000, rng)))
    crit = pd.DataFrame(rows).set_index("compartment")
    # gate: each unordered pair once, target = first
    pairs = list(itertools.combinations(comps, 2))
    byk = {}
    for k, l in pairs:
        byk.setdefault(k, []).append(l)
    ms = 5 if subsample else MIN_SHARED
    ng = 30 if run["quick"] else 200
    res = Parallel(n_jobs=jobs)(delayed(gate_row)(ds, k, ls, ng, seed, ms) for k, ls in byk.items())
    g = pd.DataFrame([r for rr in res for r in rr], columns=["k", "l", "passes", "p_pairing", "p_acrophase"])
    g["dataset"] = name
    full = pd.concat([g, g.rename(columns={"k": "l", "l": "k"})])
    gs = full.groupby("k").passes.agg(gate_pass_share="mean", n_gate_pairs="count")
    crit["gate_pass_share"], crit["n_gate_pairs"] = gs.gate_pass_share, gs.n_gate_pairs
    # centrality: A109 C
    pr, _ = ds.coherence(n_perm=10 if run["quick"] else 100, n_boot=0, seed=seed, min_shared=ms)
    both = pd.concat([pr, pr.rename(columns={"k": "l", "l": "k"})])
    cc = both.groupby("k").agg(centrality=("rho", "mean"), centrality_excess=("excess", "mean"),
                               n_partners=("rho", "count"))
    crit = crit.join(cc)
    crit["contrasts_served"] = ";".join(contrasts)
    crit["n_comp"] = len(comps)
    return crit.reset_index(), g


def rank_dataset(c: pd.DataFrame) -> pd.DataFrame:
    r = c[["dataset", "compartment", "family"]].copy()
    for col, new in (("excess", "rank_clock"), ("gate_pass_share", "rank_gate"), ("centrality", "rank_centrality")):
        r[new] = c[col].rank(ascending=False, method="average")
    r["rank_sum"] = r[["rank_clock", "rank_gate", "rank_centrality"]].sum(1)
    r["centrality"] = c.centrality
    r["rank"] = r.sort_values(["rank_sum", "centrality"], ascending=[True, False]).assign(i=range(len(r))).i.reindex(r.index) + 1
    r["tie_on_rank_sum"] = r.rank_sum.map(r.rank_sum.value_counts()) > 1
    return r.drop(columns="centrality").sort_values("rank")


def decide(crit: pd.DataFrame, ranks: pd.DataFrame, merge_neurons: bool = False) -> pd.DataFrame:
    fam = (lambda f: "neuron" if f in ("excitatory", "inhibitory", "neuron_other") else f) if merge_neurons else (lambda f: f)
    n_ds = ranks.dataset.nunique()
    top = ranks[ranks["rank"] == 1].merge(crit[["dataset", "compartment", "gate_pass_share", "centrality"]],
                                            on=["dataset", "compartment"])
    top["family"] = top.family.map(fam)
    rows = []
    for f, g in top.groupby("family"):
        rows.append(dict(family=f, n_first=len(g), n_datasets=n_ds, frac_first=len(g) / n_ds,
                         datasets_first=";".join(g.dataset), gate_min_where_first=g.gate_pass_share.min(),
                         median_centrality_top=g.centrality.median()))
    d = pd.DataFrame(rows).sort_values("n_first", ascending=False)
    d["qualifies"] = (d.frac_first >= 2 / 3) & (d.gate_min_where_first >= 0.8)
    q = d[d.qualifies].sort_values("median_centrality_top", ascending=False)
    d["decision"] = ""
    d["merged_neurons"] = merge_neurons
    d.loc[d.index == (q.index[0] if len(q) else None), "decision"] = "reference"
    return d


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = OUT / "subsample" if run["quick"] else OUT
    out.mkdir(parents=True, exist_ok=True)
    crits, gates = [], []
    for name in DATASETS:
        c, g = one_dataset(name, run, args.subsample, args.seed, args.jobs)
        print(name, len(c), "compartments", flush=True)
        crits.append(c)
        gates.append(g)
    crit = pd.concat(crits, ignore_index=True)
    crit.to_csv(out / "a135_criteria.csv", index=False)
    pd.concat(gates).to_csv(out / "a135_gate_pairs.csv", index=False)
    ranks = pd.concat([rank_dataset(c) for _, c in crit.groupby("dataset", sort=False)], ignore_index=True)
    ranks.to_csv(out / "a135_ranks.csv", index=False)
    dec = decide(crit, ranks)
    dec2 = decide(crit, ranks, merge_neurons=True)  # post hoc sensitivity, not the registered rule
    final = "reference: " + ", ".join(dec[dec.decision == "reference"].family) if (dec.decision == "reference").any() \
        else "no single reference"
    dec["final"] = final
    pd.concat([dec, dec2.assign(final="post hoc sensitivity (neurons merged)")], ignore_index=True).to_csv(
        out / "a135_decision.csv", index=False)
    pd.DataFrame(sorted({(c, family(c)) for c in crit.compartment}), columns=["compartment", "family"]).to_csv(
        out / "a135_family_map.csv", index=False)
    pd.DataFrame([dict(dataset=k, contrasts_served=";".join(v[0]), control_set_of=v[1],
                       n_control_donors=int(crit[crit.dataset == k].n_control_donors.max()))
                  for k, v in DATASETS.items()]).to_csv(out / "a135_datasets.csv", index=False)
    pd.set_option("display.width", 250)
    print(ranks[ranks["rank"] <= 3].to_string(index=False))
    print(dec.to_string(index=False), "\n", final)


if __name__ == "__main__":
    main()
