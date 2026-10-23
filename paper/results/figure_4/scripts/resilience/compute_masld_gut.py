"""A43-A45: fatty liver (pooled open MASLD cohorts) and inflamed vs non-inflamed gut within patients (TAURUS).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A43-A44 and the A44 addition; A45). Both are case-control contrasts
of one target class on a host of fixed stromal-vascular anchors, with a specificity companion on the same permutations.
  A43-A44 (`masld_pooled`, unit = patient, stratum = study, labels permuted within study): D1 = no MASLD / steatosis
    only vs MASH (primary), D2 = fibrosis F0-2 vs F3-4 (secondary); anchors endothelial, stellate, macrophage;
    primary hepatocyte d = c_high - c_low < 0, one-sided p < 0.05 (here p_less); companion hepatocyte minus
    cholangiocyte (A44 addition: the primary is "specific" only if this is < 0 too); leave-one-study-out, without the
    living-donor controls of GSE289173, per study with a DerSimonian-Laird pool; MDE of hepatocyte and cholangiocyte.
  A45 (`taurus`, unit = biopsy sample, patient = permutation block, stratum = gut site): inflamed vs non-inflamed;
    anchors fibroblast, endothelial, mural; primary T-cell d < 0 with companion d_T - d_fibroblast < 0; sensitivities
    before treatment, same-site pairs, Crohn's, UC; MDE of the T cells.
The A22 own-d* detector is applied (an anchor with p* < 0.05 and |d*| >= 0.2 is dropped and the contrast rerun). The
gate is run on the control, case and pooled donors of every analysis.

    python compute_masld_gut.py [--subsample N] [--seed S] [--jobs 4] [--mde]
Outputs (paper/results/figure_4/output/resilience/): a44_{analyses,gates,contrast,companion,meta,mde}.csv,
a45_{analyses,gates,contrast,companion,mde}.csv
Legacy: legacy/desync/validation/masld/scripts/{masld,a43_gse289173,a44_pooled}.py,
legacy/desync/validation/ibd_taurus/scripts/a45_taurus.py.
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
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "resilience"
MIN_GROUP = 15
HEP, CHOL = "liver:hepatocyte", "liver:cholangiocyte"
LIVER_ANCHORS = ["liver:endothelial", "liver:stellate", "liver:macrophage"]
STUDIES = ("GSE202379", "GSE289173", "GSE244832", "GSE256398")
TCELL, FIB = "gut:T", "gut:fibroblast"
GUT_ANCHORS = ["gut:fibroblast", "gut:endothelial", "gut:mural"]
TIE = 1e-4


def one_analysis(name, pb, donors, lab, anchors, companion, blocks, run, seed, targets=None):
    """Gates (control / case / pooled), contrast with companion and the A22 own-d* detector, for one analysis."""
    limit_threads(1)
    lab = lab[lab.isin(sh.GROUPS)]
    pb = pb[pb.donor.isin(lab.index)].reset_index(drop=True)
    lab = lab[lab.index.isin(set(pb.donor))]
    anc = [a for a in anchors if a in set(pb.compartment)]
    sets = pb.groupby("donor").compartment.apply(set)
    feas = {}
    for g in sh.GROUPS:
        s = sets.reindex(lab.index[lab == g]).dropna()
        feas[f"donors_{g}"] = int(len(s))
        feas[f"target_2anchors_{g}"] = int(s.apply(lambda x: (HEP in x or TCELL in x) and len(x & set(anc)) >= 2).sum())
    head = dict(analysis=name, anchors=";".join(anc), **feas)
    out = dict(head=head, gates=[], contrast=None, companion=None)
    if len(anc) < 2:
        out["head"]["status"] = "no anchor"
        return out
    ds = cs.Dataset(pb, anchors=anc, donors=donors.assign(lab=lab), name=name)
    for g, groups in (("control", ["control"]), ("case", ["case"]), ("pooled", ["control", "case"])):
        try:
            gt = ds.subset(donors=lab.index[lab.isin(groups)]).gate(n_perm=run["n_perm"], seed=seed)
            out["gates"].append(dict(analysis=name, group=g, p_pairing=gt.p_pairing, p_acrophase=gt.p_acrophase,
                                     clock_share=gt.clock_share, passed=bool(gt.passes)))
        except Exception as e:  # an arm with too few gateable compartments
            out["gates"].append(dict(analysis=name, group=g, passed=False, note=f"not gateable: {e}"))
    kw = dict(levels=sh.GROUPS, blocks=blocks, n_perm=run["n_perm"], n_boot=run["n_boot"], seed=seed)
    c = ds.contrast("lab", companion=companion, targets=targets, **kw)
    t = c.table
    flag = [a for a in anc if a in t.index and t.loc[a, "p_star"] < 0.05 and abs(t.loc[a, "d_star"]) >= 0.2]
    out["head"]["a22_dropped"] = ";".join(flag)
    if flag and len(anc) - len(flag) >= 2:
        anc = [a for a in anc if a not in flag]
        ds = cs.Dataset(pb, anchors=anc, donors=donors.assign(lab=lab), name=name)
        c = ds.contrast("lab", companion=companion, targets=targets, **kw)
    out["head"].update(status="ok", anchors_final=";".join(anc))
    out["contrast"] = c.table.reset_index().assign(analysis=name)
    out["companion"] = c.companion.assign(analysis=name) if c.companion is not None else None
    return out


def merge(results, out, prefix):
    tabs = {"analyses": [r["head"] for r in results], "gates": [g for r in results for g in r["gates"]]}
    pd.DataFrame(tabs["analyses"]).to_csv(out / f"{prefix}_analyses.csv", index=False)
    pd.DataFrame(tabs["gates"]).to_csv(out / f"{prefix}_gates.csv", index=False)
    for key in ("contrast", "companion"):
        parts = [r[key] for r in results if r[key] is not None]
        if parts:
            pd.concat(parts, ignore_index=True).to_csv(out / f"{prefix}_{key}.csv", index=False)


def subsample_lab(lab, n, seed, block=None):
    if n is None:
        return lab
    rng = np.random.default_rng(seed)
    if block is None:
        return sh.subsample_donors(lab, n, seed)
    keep = rng.permutation(block.reindex(lab.index).dropna().unique())[:n]  # whole blocks
    return lab[block.reindex(lab.index).isin(keep)]


def liver(args, run, out):
    pb, donors = data.load("masld_pooled")
    d1, d2 = donors.arm_d1, donors.arm_d2
    if args.subsample:  # per-study subsample so that every study keeps both arms
        keep = []
        for s in STUDIES:
            keep += list(subsample_lab(d1[(donors.study == s) & d1.isin(sh.GROUPS)], args.subsample // 4, args.seed).index)
        d1, d2 = d1.loc[keep], d2.reindex(keep)
    living = donors.living_donor.astype(bool)
    jobs = {"a44_d1": (d1, STUDIES), "a44_d2": (d2, STUDIES)}
    for s in STUDIES:
        jobs[f"a44_loso_{s}"] = (d1, tuple(x for x in STUDIES if x != s))
    jobs["a44_no_living_donors"] = (d1[~living.reindex(d1.index, fill_value=False)], STUDIES)
    tasks = []
    for i, (n, (lab, studies)) in enumerate(jobs.items()):
        lab = lab[donors.study.reindex(lab.index).isin(studies)]
        tasks.append((n, lab, f"pooled {i}"))
    for s in STUDIES:  # per study, the pooled primary's anchors
        tasks.append((f"a44_study_{s}", d1[donors.study.reindex(d1.index) == s], "study"))
    tasks.append(("a44_study_GSE289173_boramae",
                  d1[(donors.study.reindex(d1.index) == "GSE289173") & (donors.hospital.reindex(d1.index) == "Boramae")],
                  "study"))
    res = Parallel(n_jobs=args.jobs)(
        delayed(one_analysis)(n, pb, donors, lab, LIVER_ANCHORS, CHOL, None, run, args.seed + 100 * i)
        for i, (n, lab, _) in enumerate(tasks))
    merge(res, out, "a44")
    c = pd.read_csv(out / "a44_contrast.csv")
    ps = c[c.analysis.isin([f"a44_study_{s}" for s in STUDIES])]
    rows = []
    for k in (HEP, CHOL):
        t = ps[(ps.compartment == k) & ps.d.notna() & ps.lo.notna()]
        r = cs.stats.meta(t.d, (t.hi - t.lo) / 3.92)  # SE from the bootstrap 95 % width, as legacy
        rows.append(dict(compartment=k, per_study=";".join(f"{a.split('_')[-1]}:{b:+.2f}" for a, b in zip(t.analysis, t.d)),
                         **r))
    pd.DataFrame(rows).to_csv(out / "a44_meta.csv", index=False)
    if args.mde:
        lab = jobs["a44_d1"][0]
        pbl = pb[pb.donor.isin(lab.index)]
        ds = cs.Dataset(pbl, anchors=LIVER_ANCHORS, donors=donors.assign(lab=lab), name="a44_d1")
        rows = []
        for k in (HEP, CHOL):
            m = ds.contrast("lab", targets=[k], levels=sh.GROUPS, n_perm=0, n_boot=0, mde=True, seed=args.seed,
                            n_mde=10 if run["quick"] else 100).mde.iloc[0]
            rows.append(dict(compartment=k, mde=m["mde"], c_pooled=m["c_pooled"]))
        pd.DataFrame(rows).to_csv(out / "a44_mde.csv", index=False)


def eligible_samples(pb, sm, keep, block):
    """Samples with T cells + >= 2 anchors, labelled, in a block with >= 1 sample of each label."""
    sets = pb.groupby("donor").compartment.apply(set)
    ok = sets[sets.apply(lambda x: TCELL in x and len(x & set(GUT_ANCHORS)) >= 2)].index
    s = sm[sm.index.isin(ok) & sm.arm.isin(sh.GROUPS) & sm.index.isin(keep)]
    both = s.groupby(s[block]).arm.nunique()
    return s.index[s[block].isin(both.index[both == 2])]


def gut(args, run, out):
    pb, sm = data.load("taurus")
    if args.subsample:
        pat = np.random.default_rng(args.seed).permutation(sorted(sm.patient.unique()))[:max(args.subsample // 4, 6)]
        sm = sm[sm.patient.isin(pat)]
        pb = pb[pb.donor.isin(sm.index)]
    sens = {"a45_primary": (sm.index, "patient"), "a45_pretreatment": (sm.index[sm.pretreatment.astype(bool)], "patient"),
            "a45_same_site": (sm.index, "patient_site"), "a45_cd": (sm.index[sm.disease == "CD"], "patient"),
            "a45_uc": (sm.index[sm.disease == "UC"], "patient")}
    tasks = []
    for n, (keep, block) in sens.items():
        el = eligible_samples(pb, sm, keep, block)
        tasks.append((n, sm.arm.reindex(el), sm[block].reindex(el).astype(str)))
    res = []
    for i, (n, lab, blk) in enumerate(tasks):
        d = sm.assign(block=blk.reindex(sm.index))
        res.append(one_analysis(n, pb, d, lab, GUT_ANCHORS, FIB, "block", run, args.seed + 100 * i))
    merge(res, out, "a45")
    if args.mde:
        n, lab, blk = tasks[0]
        d = sm.assign(block=blk.reindex(sm.index), lab=lab)
        ds = cs.Dataset(pb[pb.donor.isin(lab.index)], anchors=GUT_ANCHORS, donors=d, name=n)
        m = ds.contrast("lab", targets=[TCELL], levels=sh.GROUPS, blocks="block", n_perm=0, n_boot=0, mde=True,
                        seed=args.seed, n_mde=10 if run["quick"] else 100).mde.iloc[0]
        pd.DataFrame([dict(compartment=TCELL, mde=m["mde"], c_pooled=m["c_pooled"])]).to_csv(out / "a45_mde.csv",
                                                                                           index=False)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--mde", action="store_true")
    ap.add_argument("--only", choices=("liver", "gut"), default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    if args.only in (None, "liver"):
        liver(args, run, out)
    if args.only in (None, "gut"):
        gut(args, run, out)
    for f in ("a44_analyses", "a44_contrast", "a44_companion", "a44_meta", "a45_analyses", "a45_companion"):
        p = out / f"{f}.csv"
        if p.exists():
            print(f, "\n", pd.read_csv(p).round(3).head(30).to_string(index=False))


if __name__ == "__main__":
    main()
