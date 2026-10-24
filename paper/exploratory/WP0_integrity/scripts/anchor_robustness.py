"""WP0 item 8: anchor robustness of every host-referenced disease contrast cited in Figures 4-5.

Per contrast and target compartment, two readouts (helpers: `paper/common/anchor_robustness.py`, public clocksync API):
  1. r*: agreement between the host of the anchor block and the host of the remaining compartments of the same
     donors (`Dataset.agreement`, as A103). With a block of < 4 compartments the split-half reliability is undefined:
     only the raw r (a lower bound of r*) is reported.
  2. the case-minus-control contrast on the anchored host (the snapshot value of the cited table, read from file,
     and recomputed here) and on the host of all compartments except the target (leave-one-out), one set of
     label permutations per host. A change of sign of d between the snapshot and the all-LOO host is flagged.
Contrast sets: A17-A22 forest (22 A20 contrasts, learned A17 anchors; Fig. 4, and Fig. 5 for the lung T cells),
A43-A45 (fatty liver, TAURUS gut), A105b (COVID-19 lung classes on two fixed anchors) and the tumour-minus-normal
epithelial coupling of the cancer cohorts (A17 cancer, learned anchors; Fig. 5). Permutations 200 (A17 registered
1,000), bootstrap 100, r*: 50 bootstrap replicates and 20 splits. Not recomputed: see the report.

    python anchor_robustness.py [--subsample N] [--seed S] [--jobs 4] [--only forest liver gut covid cancer]
Output: output/anchor_robustness/anchor_robustness.csv (subsample: output/anchor_robustness/subsample/)
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[4]
PAPER = ROOT / "paper"
sys.path.insert(0, str(ROOT / "clocksync" / "src"))
sys.path.insert(0, str(PAPER))
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import anchor_robustness as ar, cli, data  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output" / "anchor_robustness"
ARCH = Path("/home/maxine/Documents/paychere/scCircadianCancer/legacy/desync/common_output")
A17 = ARCH / "a17_unified_contrast.csv"
A17C = ARCH / "a17_cancer.csv"
A17CD = ARCH / "a17_cancer_tumour_minus_normal.csv"
RES = PAPER / "results/figure_4/output/resilience"
NAT = PAPER / "results/figure_4/output/natural_experiments"
GROUPS = ("control", "case")
FIGPANEL = {"forest": ("4", "Forest plot of all disease contrasts with MDEs (coupling d*)", "A17-A22")}
FIG5_LUNG = ("5", "Fibrotic-lung T cells lose relative amplitude, two cohorts", "A98, A5, A17, A19")
HEADLINE = {("ild_natri", "lung:T"), ("hlca_fibrotic", "lung:T")}


def _mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def _rel(p):
    try:
        return str(Path(p).relative_to(ROOT))
    except ValueError:
        return str(p)


def _rstar(ds, n_boot, n_split, seed):
    r = ar.block_agreement(ds, n_boot=n_boot, n_split=n_split, seed=seed)
    g = lambda k: r.get(k, np.nan)
    return dict(n_anchors=g("n_anchors"), n_rest=g("n_rest"), r=g("r"), r_lo=g("r_lo"), r_hi=g("r_hi"), r_star=g("r_star"),
                r_star_lo=g("r_star_lo"), r_star_hi=g("r_star_hi"), r_note=g("note"))


def _assemble(t, meta, snap, rs, n_perm, n_boot):
    """Rows of one contrast: contrast table of `host_contrast`, snapshot columns, r*."""
    t = t.copy()
    for k, v in meta.items():
        t[k] = v
    s = snap.reindex(t.target)
    t["contrast_anchored"] = s.d.to_numpy()
    t["se_anchored_snapshot"] = s.se.to_numpy()
    t["p_anchored_snapshot"] = s.p.to_numpy()
    t["d_star_anchored_snapshot"] = s.d_star.to_numpy()
    t["log_z2_anchored_snapshot"] = s.log_z2.to_numpy()
    t["snapshot_path"] = _rel(meta["snapshot_file"])
    t["snapshot_column"] = meta["snapshot_column"]
    t = t.rename(columns={"d_anchored": "contrast_anchored_recomputed", "se_d_anchored": "se_anchored_recomputed",
                          "p_anchored": "p_anchored_recomputed", "d_loo": "contrast_all_loo", "se_d_loo": "se_all_loo",
                          "p_loo": "p_all_loo", "sign_change": "sign_change_recomputed"})
    t["sign_change"] = ar.sign_change(t.contrast_anchored, t.contrast_all_loo)
    t["sign_change_z2_snapshot"] = ar.sign_change(t.log_z2_anchored_snapshot, t.log_z2_ratio_loo)
    t["n_control"], t["n_case"] = t.donors_control_loo, t.donors_case_loo
    for k, v in rs.items():
        t[k] = v
    t["n_perm"], t["n_boot"] = n_perm, n_boot
    return t


def _snap(df, d, se=None, p=None, d_star=None, log_z2=None, index="compartment"):
    out = pd.DataFrame(dict(d=df[d].to_numpy(), se=df[se].to_numpy() if se is not None else np.nan,
                            p=df[p].to_numpy() if p is not None else np.nan,
                            d_star=df[d_star].to_numpy() if d_star is not None else np.nan,
                            log_z2=df[log_z2].to_numpy() if log_z2 is not None else np.nan), index=df[index].to_numpy())
    return out[~out.index.duplicated()]


# ---- A17-A22 forest -----------------------------------------------------------------------------------------------

def forest_one(name, run, subsample, seed):
    limit_threads(1)
    warnings.simplefilter("ignore")
    res = _mod("res_shared", PAPER / "results/figure_4/scripts/resilience/_shared.py")
    reg = res.targets("contrast").set_index("name", drop=False)
    row = reg.loc[name]
    row = pd.Series(row, name=name)
    pb, donors, lab, blocks = res.load_contrast(row, subsample, seed)
    anc = [a for a in (res.anchor_set(name, "learned") or []) if a in set(pb.compartment)]
    if len(anc) < 3:
        return pd.DataFrame()
    fam = res.family_map()
    ds = cs.Dataset(pb, anchors=anc, donors=donors.assign(lab=lab), name=name)
    prim = {p.removesuffix("_coupling") for p in str(row.primary).split(";") if p and p != "nan"}
    counts = pb.assign(g=pb.donor.map(lab)).groupby(["compartment", "g"]).donor.nunique().unstack(fill_value=0)
    targets = [k for k in ds.compartments if fam.get(k, "other") != "other" and counts.get(GROUPS[0], pd.Series()).get(k, 0) >= 5
               and counts.get(GROUPS[1], pd.Series()).get(k, 0) >= 5]
    if not targets:
        return pd.DataFrame()
    t = ar.host_contrast(ds, "lab", targets, levels=GROUPS, blocks=blocks, n_perm=run["n_perm"], n_boot=run["n_boot"], seed=seed)
    rs = _rstar(ds, run["n_rstar_boot"], run["n_split"], seed)
    a = pd.read_csv(A17)
    a = a[a.target == name]
    snap = _snap(a.assign(se=(a.hi_learned - a.lo_learned) / 3.92), "diff_learned", "se", "p_perm_learned", "d_star",
                 "log_z2_ratio")
    fig, panel, am = FIGPANEL["forest"]
    meta = dict(figure=fig, panel=panel, amendment=am, dataset=row.dataset, contrast=name, case_control="case vs control",
                snapshot_file=A17, snapshot_column="diff_learned", target_is_primary=False)
    t = _assemble(t, meta, snap, rs, run["n_perm"], run["n_boot"])
    t["target_is_primary"] = t.target.map(lambda k: k.split(":")[-1] in prim or any(k.endswith(":" + p) for p in prim))
    t["target_is_anchor"] = t.target.isin(anc)
    extra = t[[(name, k) in HEADLINE for k in t.target]].copy()
    if len(extra):
        extra[["figure", "panel", "amendment"]] = FIG5_LUNG
        t = pd.concat([t, extra], ignore_index=True)
    return t


# ---- A43-A45 ------------------------------------------------------------------------------------------------------

def liver_one(run, subsample, seed):
    limit_threads(1)
    warnings.simplefilter("ignore")
    res = _mod("res_shared", PAPER / "results/figure_4/scripts/resilience/_shared.py")
    pb, donors = data.load("masld_pooled")
    lab = donors.arm_d1
    lab = lab[lab.isin(GROUPS)]
    lab = res.subsample_donors(lab, subsample, seed)
    pb = pb[pb.donor.isin(lab.index)].reset_index(drop=True)
    lab = lab[lab.index.isin(set(pb.donor))]
    anc = [a for a in ["liver:endothelial", "liver:stellate", "liver:macrophage"] if a in set(pb.compartment)]
    ds = cs.Dataset(pb, anchors=anc, donors=donors.assign(lab=lab), name="a44_d1")
    counts = pb.assign(g=pb.donor.map(lab)).groupby(["compartment", "g"]).donor.nunique().unstack(fill_value=0)
    targets = [k for k in ds.compartments if counts.loc[k, "control"] >= 5 and counts.loc[k, "case"] >= 5]
    t = ar.host_contrast(ds, "lab", targets, levels=GROUPS, blocks=None, n_perm=run["n_perm"], n_boot=run["n_boot"], seed=seed)
    rs = _rstar(ds, run["n_rstar_boot"], run["n_split"], seed)
    a = pd.read_csv(RES / "a44_contrast.csv")
    a = a[a.analysis == "a44_d1"]
    snap = _snap(a.assign(se=(a.hi - a.lo) / 3.92), "d", "se", "p", "d_star", "log_z2_ratio")
    meta = dict(figure="4", panel="Fatty liver: liver-wide loss carried by procurement (ED)", amendment="A43, A44, A96, A45",
                dataset="masld_pooled (MASH vs no MASLD/steatosis, D1)", contrast="a44_d1", case_control="MASH vs no MASH",
                snapshot_file=RES / "a44_contrast.csv", snapshot_column="d (analysis a44_d1)")
    t = _assemble(t, meta, snap, rs, run["n_perm"], run["n_boot"])
    t["target_is_primary"] = t.target == "liver:hepatocyte"
    t["target_is_anchor"] = t.target.isin(anc)
    return t


def gut_one(run, subsample, seed):
    limit_threads(1)
    warnings.simplefilter("ignore")
    res = _mod("res_shared", PAPER / "results/figure_4/scripts/resilience/_shared.py")
    pb, sm = data.load("taurus")
    if subsample:
        pat = np.random.default_rng(seed).permutation(sorted(sm.patient.unique()))[:max(subsample // 4, 6)]
        sm = sm[sm.patient.isin(pat)]
        pb = pb[pb.donor.isin(sm.index)]
    anc = ["gut:fibroblast", "gut:endothelial", "gut:mural"]
    sets = pb.groupby("donor").compartment.apply(set)
    ok = sets[sets.apply(lambda x: "gut:T" in x and len(x & set(anc)) >= 2)].index
    s = sm[sm.index.isin(ok) & sm.arm.isin(GROUPS)]
    both = s.groupby(s.patient).arm.nunique()
    el = s.index[s.patient.isin(both.index[both == 2])]
    lab = sm.arm.reindex(el)
    d = sm.assign(block=sm.patient.astype(str), lab=lab)
    pb = pb[pb.donor.isin(el)].reset_index(drop=True)
    ds = cs.Dataset(pb, anchors=[a for a in anc if a in set(pb.compartment)], donors=d, name="a45_primary")
    counts = pb.assign(g=pb.donor.map(lab)).groupby(["compartment", "g"]).donor.nunique().unstack(fill_value=0)
    targets = [k for k in ds.compartments if counts.loc[k, "control"] >= 5 and counts.loc[k, "case"] >= 5]
    t = ar.host_contrast(ds, "lab", targets, levels=GROUPS, blocks="block", n_perm=run["n_perm"], n_boot=run["n_boot"], seed=seed)
    rs = _rstar(ds, run["n_rstar_boot"], run["n_split"], seed)
    a = pd.read_csv(RES / "a45_contrast.csv")
    a = a[a.analysis == "a45_primary"]
    snap = _snap(a.assign(se=(a.hi - a.lo) / 3.92), "d", "se", "p", "d_star", "log_z2_ratio")
    meta = dict(figure="4", panel="Inflamed gut stroma intact (ED)", amendment="A43, A44, A96, A45",
                dataset="taurus (inflamed vs non-inflamed, patient blocks)", contrast="a45_primary",
                case_control="inflamed vs non-inflamed", snapshot_file=RES / "a45_contrast.csv",
                snapshot_column="d (analysis a45_primary)")
    t = _assemble(t, meta, snap, rs, run["n_perm"], run["n_boot"])
    t["target_is_primary"] = t.target == "gut:T"
    t["target_is_anchor"] = t.target.isin(anc)
    return t


# ---- A105b COVID-19 lung classes ----------------------------------------------------------------------------------

def covid_one(run, subsample, seed):
    limit_threads(1)
    warnings.simplefilter("ignore")
    S = _mod("embryo_shared", PAPER / "results/figure_2/scripts/embryo/_shared.py")
    melms, ms = data.load("covid_melms", subsample=subsample, seed=seed)
    lung = S.eligible(melms, 20)
    anc = ["Endothelial cells", "Fibroblasts"]
    ds = S.dataset(lung, ms, anchors=anc)
    targets = [k for k in ds.compartments]
    t = ar.host_contrast(ds, "group", targets, levels=("Control", "COVID-19"), blocks="none", n_perm=run["n_perm"],
                         n_boot=run["n_boot"], seed=seed)
    rs = _rstar(ds, run["n_rstar_boot"], run["n_split"], seed)
    a = pd.read_csv(NAT / "a105_arms.csv")
    a = a[a.anchor.notna()]
    snap = _snap(a.assign(se=(a.hi - a.lo) / 3.92), "d", "se", "p", "d_star", "log_z2_ratio")
    meta = dict(figure="4", panel="COVID-19 decedents keep organ clocks and endothelial sharing (lung arms)", amendment="A105, A105b",
                dataset="covid_melms (COVID-19 vs 7 living surgical controls)", contrast="melms_covid",
                case_control="COVID-19 vs surgical control", snapshot_file=NAT / "a105_arms.csv", snapshot_column="d")
    # t: contrast columns are named levels[0]/levels[1]; host_contrast keys on control/case, so map here
    t = _assemble(t.rename(columns=lambda c: c), meta, snap, rs, run["n_perm"], run["n_boot"])
    t["target_is_primary"] = False
    t["target_is_anchor"] = t.target.isin(anc)
    return t


# ---- cancer cohorts -----------------------------------------------------------------------------------------------

def cancer_one(study, run, subsample, seed):
    limit_threads(1)
    warnings.simplefilter("ignore")
    from common.anchors import anchors as get_anchors
    pb, _ = data.load("cancer_sc", study=study, cnv=False)
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:max(subsample, 12)]
        pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
    try:
        anc = [c for c in get_anchors(f"cancer_{study}", "learned") if c in set(pb.compartment)]
    except KeyError:
        return pd.DataFrame()
    pair = ("tumour:epithelial", "adjacent_normal:epithelial")
    if len(anc) < 3 or not set(pair) <= set(pb.compartment):
        return pd.DataFrame()
    ds = cs.Dataset(pb, anchors=anc, name=study)
    rows = {}
    for host, d in (("anchored", ds), ("loo", ar.loo_dataset(ds))):
        f = d.fit(n_boot=run["n_boot"], seed=seed, gate=False)
        x = (f.replicates[pair[0]] - f.replicates[pair[1]]).to_numpy()
        x = x[np.isfinite(x)]
        diff = f.table.coupling[pair[0]] - f.table.coupling[pair[1]]
        rows[host] = (diff, float(np.std(x, ddof=1)) if len(x) > 2 else np.nan, int(pb[pb.compartment == pair[0]].donor.nunique()))
    rs = _rstar(ds, run["n_rstar_boot"], run["n_split"], seed)
    snap = pd.read_csv(A17CD)
    snap = snap[(snap.target == f"cancer_{study}") & (snap.rule == "learned") & (snap.cls == "epithelial")]
    t = pd.DataFrame([dict(
        figure="5", panel="Malignant cells: clock-shaped share; decomposition amplitude/coherence (A87 on A17 anchors)",
        amendment="A86, A87, A97 (A17 cancer coupling)", dataset=f"cancer_sc:{study}", contrast=f"cancer_{study}",
        case_control="tumour vs adjacent normal (epithelial)", target="tumour:epithelial - adjacent_normal:epithelial",
        target_is_primary=True, target_is_anchor=False, anchors=";".join(anc),
        contrast_anchored=snap["diff"].iloc[0] if len(snap) else np.nan,
        se_anchored_snapshot=snap.se.iloc[0] if len(snap) else np.nan, p_anchored_snapshot=np.nan,
        snapshot_path=_rel(A17CD), snapshot_column="diff (rule learned, cls epithelial)",
        contrast_anchored_recomputed=rows["anchored"][0], se_anchored_recomputed=rows["anchored"][1],
        contrast_all_loo=rows["loo"][0], se_all_loo=rows["loo"][1], p_all_loo=np.nan,
        n_control=rows["loo"][2], n_case=rows["loo"][2], n_perm=0, n_boot=run["n_boot"], **rs)])
    t["sign_change"] = ar.sign_change(t.contrast_anchored, t.contrast_all_loo)
    t["sign_change_recomputed"] = ar.sign_change(t.contrast_anchored_recomputed, t.contrast_all_loo)
    return t


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="+", default=["forest", "liver", "gut", "covid", "cancer"])
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    run.update(n_perm=20 if quick else 200, n_boot=10 if quick else 100, n_rstar_boot=5 if quick else 50,
               n_split=5 if quick else 20)
    out = OUT / "subsample" if quick else OUT
    out.mkdir(parents=True, exist_ok=True)
    res = _mod("res_shared", PAPER / "results/figure_4/scripts/resilience/_shared.py")
    sp = _mod("spread_shared", PAPER / "results/figure_4/scripts/spread_multiplicity/_shared.py")
    jobs = []
    if "forest" in args.only:
        names = list(sp.CLASSES)
        if quick:
            names = [n for n in names if n in ("ild_natri", "dcm", "skin_ad", "islets_t2d")]
        jobs += [delayed(forest_one)(n, run, args.subsample, args.seed) for n in names]
    if "liver" in args.only:
        jobs.append(delayed(liver_one)(run, args.subsample, args.seed))
    if "gut" in args.only:
        jobs.append(delayed(gut_one)(run, args.subsample, args.seed))
    if "covid" in args.only:
        jobs.append(delayed(covid_one)(run, args.subsample, args.seed))
    if "cancer" in args.only:
        reg = pd.read_csv(PAPER / "common/anchors_a17.csv")
        studies = sorted(t.removeprefix("cancer_") for t in reg.target.unique() if t.startswith("cancer_"))
        passed = pd.read_csv(A17C)
        passed = set(passed[(passed.rule == "learned") & passed.gate_passed.fillna(False).astype(bool)].target.str.removeprefix("cancer_"))
        studies = [s for s in studies if s in passed]
        if quick:
            studies = [s for s in studies if s in ("SRP273482", "GSE132465", "GSE181919")]
        jobs += [delayed(cancer_one)(s, run, args.subsample, args.seed) for s in studies]
    parts = Parallel(n_jobs=args.jobs, verbose=5)(jobs)
    df = pd.concat([p for p in parts if len(p)], ignore_index=True)
    cols = ["figure", "panel", "amendment", "dataset", "contrast", "target", "target_is_primary", "target_is_anchor",
            "case_control", "anchors", "n_anchors", "n_rest", "r", "r_lo", "r_hi", "r_star", "r_star_lo", "r_star_hi", "r_note",
            "contrast_anchored", "se_anchored_snapshot", "p_anchored_snapshot", "d_star_anchored_snapshot",
            "snapshot_path", "snapshot_column", "contrast_anchored_recomputed", "se_anchored_recomputed",
            "p_anchored_recomputed", "contrast_all_loo", "se_all_loo", "p_all_loo", "d_star_loo", "p_star_loo",
            "sign_change", "sign_change_recomputed", "log_z2_anchored_snapshot", "log_z2_ratio_loo",
            "sign_change_z2_snapshot", "n_control", "n_case", "status_loo", "n_perm", "n_boot"]
    df = df.reindex(columns=cols)
    df.to_csv(out / "anchor_robustness.csv", index=False)
    s = df[df.target_is_primary == True]  # noqa: E712
    print(s[["contrast", "target", "r_star", "r", "contrast_anchored", "contrast_all_loo", "sign_change"]].round(3).to_string(), flush=True)
    print(f"{len(df)} rows; sign_change {int((df.sign_change == True).sum())}", flush=True)  # noqa: E712


if __name__ == "__main__":
    main()
