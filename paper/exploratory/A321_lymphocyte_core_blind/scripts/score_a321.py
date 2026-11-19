"""A321 (registered): blind core-loop test of tissue lymphocytes on atlases not used by A281 / A287 / A294.

Per atlas (genome-wide table): lymphocyte compartments with >= 15 donors shared with a partner (fibroblast-stromal, else
endothelial, else clock-carrying parenchyma). Per lymphocyte compartment: core-loop signed coherence C (11 genes,
A294 `CORE`) with per-gene permutations within strata, the 14-gene C, a power control (partner thinned to the
lymphocytes' clock-gene count, read on the core loop) and the percentile among 200 expression-matched 11-gene sets.
H321.1: Stouffer combination of the primary compartments' one-sided core-loop p over informative atlases.

    python score_a321.py --stage design | score [--jobs 8]
Outputs (`../output/`): design.csv | compartments.csv, tests.csv.
"""
from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

from clocksync.decompose import signed_coherence  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import cli  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


A294 = _load("a321_a294", PAPER / "exploratory/A294_lymphocyte_hormonal/scripts/score_a294.py")
A281 = A294.A281
import cloudpickle  # noqa: E402

for _m in (A294, A281, A281.A150):
    cloudpickle.register_pickle_by_value(_m)

OUT = HERE.parents[1] / "output"
SEED, MIN_SHARED, N_SETS = 321, 15, 200
ATLASES = ["skin_ganier_genome", "lung_natri_genome", "liver_xu_genome", "organ_adipose_genome"]
CORE = A294.CORE
CLOCK = A281.CLOCK
PARENCHYMA_OK = ("parenchyma", "epithelial")  # clock-carrying classes (A272-A283)


def load(name):
    pb, donors = A281.A150.GRS.genome_table(name)
    if name == "lung_natri_genome":  # control donors only (registered)
        ctrl = donors.index[_is_control(donors)]
        pb = pb[pb.donor.isin(ctrl)].reset_index(drop=True)
    return pb


def _is_control(donors):
    for c in ("disease", "condition", "status", "group", "is_control", "control"):
        if c in donors.columns:
            v = donors[c].astype(str).str.lower()
            if c in ("is_control", "control"):
                return v.isin(("true", "1"))
            return v.str.contains("control|normal|healthy")
    raise KeyError(f"no control column among {list(donors.columns)}")


def design(name, pb):
    nd = pb.groupby("compartment").donor.nunique()
    fam = {c: A281.A150.FAM.family(c) for c in nd.index}
    don = {c: set(pb.donor[pb.compartment == c]) for c in nd.index}
    rows = []
    for ly in [c for c in nd.index if fam[c] == "lymphocyte"]:
        partner = None
        for want in ("fibroblast_stromal", "endothelial") + PARENCHYMA_OK:
            cands = [c for c in nd.index if fam[c] == want and len(don[c] & don[ly]) >= MIN_SHARED]
            if cands:
                partner = max(cands, key=lambda c: len(don[c] & don[ly]))
                break
        rows.append(dict(atlas=name, lymphocyte=ly, partner=partner, partner_family=fam.get(partner),
                         shared_donors=len(don[partner] & don[ly]) if partner else 0, lymph_donors=int(nd[ly])))
    d = pd.DataFrame(rows)
    if len(d):
        d = d[d.shared_donors >= MIN_SHARED].copy()
        d["primary"] = False
        if len(d):
            d.loc[d.shared_donors.idxmax(), "primary"] = True
    return d


def task(name, ly, partner, tl, ts, matched, seed):
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(seed)
    d = sorted(set(tl.donor) & set(ts.donor))
    tl = tl[tl.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    ts = ts[ts.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    target = float(tl[CLOCK].sum(1).median())
    core = list(CORE.genes)
    c11, p11 = A281.perm_test(A281.y_of(tl, core), A281.groups_of(tl), A281.N_PERM, rng, wf=CORE)
    c14, p14 = A281.perm_test(A281.y_of(tl, CLOCK), A281.groups_of(tl), A281.N_PERM, rng)
    thin_p = []
    for _ in range(A281.N_THIN):
        tt = A281.thin(ts, target, rng)
        thin_p.append(A281.perm_test(A281.y_of(tt, core), A281.groups_of(tt), A281.N_PERM_THIN, rng, wf=CORE)[1])
    pm = np.array([signed_coherence(A281.y_of(tl, list(w.genes)), w) for w in matched])
    pm = pm[np.isfinite(pm)]
    return dict(atlas=name, lymphocyte=ly, partner=partner, shared_donors=len(d), lymph_clock_count_median=target,
                C11=c11, p11=p11, C14=c14, p14=p14, core_power=float(np.mean(np.array(thin_p) < 0.05)),
                matched_median=float(np.median(pm)), matched_q95=float(np.quantile(pm, 0.95)),
                matched_pct=float(np.mean(pm < c11)), matched_valid=len(pm))


def matched_sets(pb, k):
    genes = [c for c in pb.columns if c not in ("donor", "compartment", "stratum", "n_cells", "depth")]
    X = pb[genes].to_numpy(float)
    frac = (X / pb.depth.to_numpy(float)[:, None]).mean(0)
    keep = (X.mean(0) >= 1) | np.isin(genes, CLOCK)
    lm = pd.Series(np.log(frac + 1e-12), index=genes)[keep]
    return pseudo_waveforms(lm, N_SETS, np.random.default_rng([SEED, k]), CORE, "matched")


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["design", "score"], required=True)
    args = ap.parse_args()
    cli.setup(args)
    OUT.mkdir(parents=True, exist_ok=True)
    tabs, des = {}, []
    for name in ATLASES:
        pb = load(name)
        tabs[name] = pb
        d = design(name, pb)
        print(name, "donors", pb.donor.nunique(), "\n", d.to_string(), flush=True)
        des.append(d)
    des = pd.concat(des, ignore_index=True)
    des.to_csv(OUT / "design.csv", index=False)
    if args.stage == "design":
        return
    tasks = []
    for k, name in enumerate(ATLASES):
        d = des[des.atlas == name]
        if not len(d):
            continue
        pb = tabs[name]
        matched = matched_sets(pb, k)
        cols = sorted({g for w in matched for g in w.genes} | set(CLOCK))
        sub = pb[["donor", "compartment", "stratum", "n_cells", "depth"] + cols]
        for j, r in enumerate(d.itertuples()):
            tasks.append(delayed(task)(name, r.lymphocyte, r.partner, sub[sub.compartment == r.lymphocyte],
                                       sub[sub.compartment == r.partner], matched, [SEED, k, j]))
    comp = pd.DataFrame(Parallel(n_jobs=args.jobs)(tasks)).merge(des[["atlas", "lymphocyte", "primary"]],
                                                                   on=["atlas", "lymphocyte"])
    comp.to_csv(OUT / "compartments.csv", index=False)
    prim = comp[comp.primary & (comp.core_power >= 0.8)]
    res = dict(informative_primary=len(prim), atlases=";".join(prim.atlas))
    if len(prim) >= 2:
        z = norm.isf(prim.p11.to_numpy()).sum() / np.sqrt(len(prim))
        p = float(norm.sf(z))
        h2 = bool(prim.matched_pct.mean() > 0.5 and (prim.C11 > prim.matched_median).all())
        res.update(stouffer_z=float(z), p_combined=p, H321_1=p < 0.05, mean_matched_pct=float(prim.matched_pct.mean()),
                   H321_2=h2 if p < 0.05 else None)
    else:
        res.update(H321_1=None, note="fewer than 2 informative primary compartments: not decidable")
    pd.DataFrame([res]).to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(comp.round(4).to_string(index=False))
    print(pd.Series(res).to_string(), flush=True)


if __name__ == "__main__":
    main()
