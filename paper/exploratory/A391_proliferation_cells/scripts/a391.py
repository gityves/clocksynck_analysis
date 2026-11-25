"""A391: does proliferation quiet the tumour's stroma? (and single-cell accounting of the bulk gain loss, A304)

Registered in paper/docs/registered_experiments.md (A391). Host and units: the A387 Part 2 table (34 cohorts; host = the
patient's normal fibroblasts, endothelium, mural). Proliferation index: mean within-cohort z of log CPM of the Hallmark
G2M genes (clock genes removed; A127 score) in the patient's malignant / tumour-epithelium pseudobulk, from the A388
genome tables (`A388_malignant_tumour_host/output/cache/genome_<cohort>.parquet`). Arms: high vs low, median split within
cohorts with >= 4 labelled patients.

  power  (no arm contrast of a real clock quantity)
         counts per arm and compartment; power of the one-sided loudness test (log z2 ratio high / low < 0, labels
         permuted within cohort, 200 permutations) at damping s of the compartment's host-locked part in the "high" arm
         (`clocksync.power.damped`, labels themselves permuted first: exchangeable arms), s in {1, 0.8, 0.63, 0.5};
         matched-set comparators: the same contrast on the 200 A387 matched sets under one random label permutation
         each (q05).                                                   -> output/power/*.csv
  score  (after the power note)

    python a391.py --stage power [--jobs 6] [--genome-jobs 3] [--smoke]     (paper/common/run.sh, RUN_CAP_GB 30)
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(EXP / "A387_tumour_full_model/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a387 as R  # noqa: E402
from a387 import ANCHORS, BASE, C, NORMAL, TARGETS, log  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync import power as PW  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.stats import perm_p  # noqa: E402

OUT = HERE.parents[1] / "output"
A388_CACHE = EXP / "A388_malignant_tumour_host/output/cache"
STROMA = {"fibroblast": TARGETS["fibroblast"], "endothelial": TARGETS["endothelial"]}
SCALES, N_REP, N_PERM, MIN_COHORT, SEED = (1.0, 0.8, 0.63, 0.5), 40, 200, 4, 391


def g2m_genes():
    spec = importlib.util.spec_from_file_location("fig5_shared", PAPER / "results/figure_5/scripts/cancer/_shared.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.hallmark_g2m()


def proliferation(table):
    """Patient -> proliferation index (A387 donor ids `<cohort>|<patient>`), from the A388 malignant genome tables."""
    genes = g2m_genes()
    out = {}
    for f in sorted(A388_CACHE.glob("genome_*.parquet")):
        study = f.stem.replace("genome_", "")
        G = pd.read_parquet(f)
        mal = G[G.compartment.str.startswith("tumour:") & ~G.compartment.isin(
            ["tumour:fibroblast", "tumour:endothelial", "tumour:mural"])]
        if mal.empty:
            continue
        mal = mal.groupby("donor").first()                     # one malignant pseudobulk per patient
        g = [x for x in genes if x in mal.columns]
        lcpm = np.log1p(mal[g].div(mal.depth, axis=0) * 1e6)
        z = (lcpm - lcpm.mean()) / lcpm.std(ddof=0).replace(0, np.nan)
        for d, v in z.mean(axis=1).items():
            out[f"{study}|{d}"] = float(v)
    return pd.Series(out, name="prolif")


def arms(table, prolif):
    hosted = set(table.donor[table.compartment.isin(ANCHORS)])
    p = prolif[prolif.index.isin(hosted)]
    lab = {}
    for study, s in p.groupby(p.index.str.split("|").str[0]):
        if len(s) < MIN_COHORT:
            continue
        r = s.rank(method="first")
        for d, v in r.items():
            lab[d] = "high" if v > len(s) / 2 else "low"
    return pd.Series(lab, name="arm")


def power_one(table, lab, k, rd, s, rep):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng([SEED, int(1000 * s), rep, rd == "core", hash(k) % 1000])
    lr = lab.copy()
    for _, idx in lab.groupby(lab.index.str.split("|").str[0]).groups.items():
        lr[idx] = rng.permutation(lab[idx].to_numpy())
    ds = cs.Dataset(table[table.donor.isin(lab.index)], anchors=ANCHORS, waveform=R.READOUTS[rd], mixed_strata=True)
    h = A._hosts(A.Data.of(ds.study), ANCHORS, [k])[k]
    t = PW.damped(ds.study, lr, ANCHORS, k, s, rng, case="high", h=h)
    if t is None:
        return dict(target=k, readout=rd, s=s, rep=rep, p=np.nan)
    t = t.rename(columns={"resection": "donor"})
    dd = cs.Dataset(t, anchors=ANCHORS, waveform=R.READOUTS[rd], mixed_strata=True)
    blocks = pd.Series(lr.index.str.split("|").str[0], index=lr.index)
    c = dd.contrast(lr, levels=("low", "high"), targets=[k], n_perm=N_PERM, n_boot=0, seed=int(rng.integers(1e9)),
                    blocks=blocks, host_prior="arm")
    obs = float(c.observed["log_z2_ratio"][k])
    return dict(target=k, readout=rd, s=s, rep=rep, log_z2=obs, p=perm_p(c.null["log_z2_ratio"][k], obs, "less"))


def comp_one(G, lab, wfs, k, rd):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = []
    blocks = pd.Series(lab.index.str.split("|").str[0], index=lab.index)
    for i, w in wfs:
        rng = np.random.default_rng([SEED, 9, i, rd == "core"])
        lr = lab.copy()
        for _, idx in lab.groupby(blocks).groups.items():
            lr[idx] = rng.permutation(lab[idx].to_numpy())
        try:
            dd = cs.Dataset(G[G.donor.isin(lab.index)], anchors=ANCHORS, waveform=w, mixed_strata=True)
            c = dd.contrast(lr, levels=("low", "high"), targets=[k], n_perm=1, n_boot=0, seed=i, blocks=blocks,
                            host_prior="arm")
            out.append(dict(target=k, readout=rd, set=i, log_z2=float(c.observed["log_z2_ratio"][k])))
        except Exception as e:  # noqa: BLE001
            out.append(dict(target=k, readout=rd, set=i, error=str(e)[:150]))
    return out


def stage_power(args):
    import a387_score as S
    from clocksync.paired import pseudo_waveforms
    out = OUT / ("power_smoke" if args.smoke else "power")
    out.mkdir(parents=True, exist_ok=True)
    table = R.pooled_table(args.jobs, args.smoke)
    prolif = proliferation(table)
    lab = arms(table, prolif)
    prolif.to_csv(out / "proliferation_index.csv")
    lab.to_csv(out / "arms.csv")
    cnt = []
    for kind, v in {**STROMA, "epithelium": TARGETS["epithelium"]}.items():
        don = set(table.donor[table.compartment == v])
        for a in ("high", "low"):
            cnt.append(dict(target=v, arm=a, patients=int(sum(d in don for d in lab.index[lab == a]))))
    cnt = pd.DataFrame(cnt)
    cnt.to_csv(out / "counts.csv", index=False)
    log(f"labelled patients {len(lab)}; cohorts {lab.index.str.split('|').str[0].nunique()}\n" + cnt.to_string(index=False))
    nrep = 4 if args.smoke else N_REP
    jobs = [(k, rd, s, r) for k in STROMA.values() for rd in R.READOUTS for s in SCALES for r in range(nrep)]
    res = Parallel(n_jobs=args.jobs)(delayed(power_one)(table, lab, *j) for j in jobs)
    P = pd.DataFrame(res)
    if "log_z2" not in P:
        P["log_z2"] = np.nan
    P.to_csv(out / "power_reps.csv", index=False)
    # a replicate whose log z2 ratio is undefined (an arm's z2 truncated at 0) counts as not detected (conservative)
    curve = P.groupby(["target", "readout", "s"]).agg(power=("p", lambda x: float(np.mean(x.fillna(1.0) < 0.05))),
                                                       n=("p", "size"), undefined=("p", lambda x: int(x.isna().sum())),
                                                       log_z2_median=("log_z2", "median")).reset_index()
    curve["log_z2_expected"] = 2 * np.log(curve.s)
    curve.to_csv(out / "power_curve.csv", index=False)
    log("power curve\n" + curve.round(3).to_string(index=False))
    # matched-set comparators under permuted labels
    keep = set(zip(table.donor, table.compartment))
    studies = sorted({d.split("|", 1)[0] for d in lab.index})
    Gs = Parallel(n_jobs=args.genome_jobs)(delayed(S.genome_cohort)(s, {x for x in keep if x[0].startswith(s + "|")})
                                           for s in studies)
    G = pd.concat([g for g in Gs if g is not None and len(g)], ignore_index=True)
    G = G.drop(columns=["stratum"]).merge(table[["donor", "compartment", "stratum"]], on=["donor", "compartment"])
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth", "n_cells")]
    genes = [g for g in genes if G[g].notna().all()]
    ref = G[G.compartment.isin(ANCHORS)]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
    n_set = 10 if args.smoke else R.N_SET
    sets = pseudo_waveforms(lm, n_set, np.random.default_rng([387, 0]), BASE, "matched")
    used = sorted({g for w in sets for g in w.genes} | set(BASE.genes))
    Gu = G[["donor", "compartment", "stratum", "depth"] + used]
    tasks = [(k, rd, [(i, w if rd == "14" else C.restrict(w, BASE, C.core(BASE))) for i, w in enumerate(sets)])
             for k in STROMA.values() for rd in R.READOUTS]
    comp = Parallel(n_jobs=args.jobs)(delayed(comp_one)(Gu, lab, ch, k, rd) for k, rd, ws in tasks
                                      for ch in S.chunks(ws, args.jobs))
    comp = pd.DataFrame([r for x in comp for r in x])
    comp.to_csv(out / "companions_permuted.csv", index=False)
    summ = comp.groupby(["target", "readout"]).log_z2.agg(q05=lambda x: float(np.nanquantile(x, .05)),
                                                          median="median", n=lambda x: int(x.notna().sum())).reset_index()
    summ.to_csv(out / "companions_summary.csv", index=False)
    log("companions under permuted labels\n" + summ.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power",), required=True)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--genome-jobs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    stage_power(args)
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
