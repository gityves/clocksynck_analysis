"""A88: are cue rhythms decoupled by fasting or surgery? Within-sample core-cue coupling relative to core-locked.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A88; sets from A81b). Per unit: r_cue = |ccorr(core, cue)| and r_lock =
|ccorr(core, locked)| over the unit's samples (one projector per set, Poisson-corrected states; floor 0.02), rho = r_cue /
r_lock; Delta = log rho minus log rho of the matched GTEx fast-death tissue. Units: outpatient skin (psoriasis and atopic
dermatitis, gated A54 / A53 series, uninvolved and lesional samples), endoscopy (inflamed vs non-inflamed IBD, gated A55
series), surgery after fasting (TCGA adjacent normal and tumour, gated A52 cancers); reference: GTEx fast deaths (all donors
with a time of death, Hardy 1-2) of both skins, colon transverse and sigmoid, lung, liver, breast, prostate, thyroid, kidney
cortex, stomach, oesophagus mucosa, and (descriptive) ventilated evaluation donors (A81 half, seed 81). Matched tissue: skin
-> both skins, IBD and COAD -> colon, other cancers by organ; BLCA and UCEC unmatched (excluded). H88: mean Delta of surgical
and endoscopy units (uninvolved) minus mean Delta of outpatient skin units (uninvolved) < 0, unit bootstrap (10,000),
one-sided p < 0.05. Secondary: digestive (COAD, LIHC, IBD) vs other surgical units. GTEx values stay in memory (aggregates
only in the outputs).
Needs compute_other_diseases.py, compute_replication.py and cancer/compute_tcga_paired.py outputs with the same
--subsample.

    python compute_cue_context.py [--subsample N] [--n-boot 10000]
Outputs (paper/results/figure_5/output/psoriasis_details/): a88_gtex.csv, a88_units.csv, a88_verdict.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a88_cue_context.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("psoriasis")))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _cues as Q  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
from common import cli, gtex  # noqa: E402
from compute_cues_diseases import cancer_table  # noqa: E402

PD = H.PD
FLOOR, SEED, SEED_HALVES = 0.02, 88, 81
GTEX = ["Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)", "Colon - Transverse", "Colon - Sigmoid",
        "Lung", "Liver", "Breast - Mammary Tissue", "Prostate", "Thyroid", "Kidney - Cortex", "Stomach", "Esophagus - Mucosa"]
MATCH = {"skin": GTEX[:2], "colon": GTEX[2:4], "COAD": GTEX[2:4], "LIHC": ["Liver"], "BRCA": ["Breast - Mammary Tissue"],
         "PRAD": ["Prostate"], "KIRC": ["Kidney - Cortex"], "KIRP": ["Kidney - Cortex"], "LUAD": ["Lung"],
         "THCA": ["Thyroid"]}


def coupling(x: pd.DataFrame, counts: bool, wfs: dict) -> dict:
    """|ccorr| of the core state with the cue and the locked states over the rows of `x` (counts + depth, or arrays)."""
    st = {}
    for k, wf in wfs.items():
        g = [y for y in wf.genes if y in x.columns]
        w = wf.subset(g) if k == "core" else H.restrict(wf, g)
        L = PD.projector(x, w, counts=counts)
        st[k] = PD.states(x, L, w, depth=x.depth.to_numpy(float) if counts else None)
    r = {k: max(abs(PD.ccorr(*st["core"], *st[k])), FLOOR) for k in ("cue", "locked")}
    return dict(n=len(x), r_cue=r["cue"], r_lock=r["locked"], log_rho=float(np.log(r["cue"] / r["locked"])))


def halves(subsample, seed):
    """A81 donor halves (selection / evaluation per Hardy group, from one generator seeded 81); `subsample`: that many
    random donors per group first (check runs)."""
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs[["DTHHRDY", "tod", "subject"]].copy()
    finally:
        a.file.close()
    o["subject"] = o.subject.astype(str)
    rng, out = np.random.default_rng(SEED_HALVES), {}
    pick = np.random.default_rng(seed)
    for grp, hardy in (("fast", [1, 2]), ("vent", [0])):
        s = np.array(sorted(o.loc[o.DTHHRDY.isin(hardy) & o.tod.notna(), "subject"].unique()))
        s = s[rng.permutation(len(s))]
        sel, ev = set(s[: len(s) // 2]), set(s[len(s) // 2:])
        if subsample:
            keep = set(pick.permutation(s)[:subsample])
            sel, ev = sel & keep, ev & keep
        out[grp] = {"sel": sel, "eval": ev}
    return out


def gtex_reference(wfs, subsample, seed) -> pd.DataFrame:
    H_ = halves(subsample, seed)
    fast, vent = H_["fast"]["sel"] | H_["fast"]["eval"], H_["vent"]["eval"]
    genes = Q.all_set_genes(wfs)
    rows = []
    for tis in GTEX:
        for grp, hardy, sub in (("gtex_fast", (1, 2), fast), ("gtex_vent_eval", (0,), vent)):
            d = Q.gtex_rows(tis, hardy, sub, genes)
            if len(d) >= (8 if subsample else 20):
                rows.append(dict(context=grp, unit=tis, **coupling(d.drop(columns=["hardy", "tod", "sex", "age_bracket", "subject"]), True, wfs)))
    return pd.DataFrame(rows)


def unit_rows(context, name, tissue, samples, w, wfs):
    rows = []
    for k, x in samples.groupby("cohort", sort=False):
        if k not in w:
            continue
        case, ref = H.split_pairs(x, H.LEVELS if "lesional" in set(x.condition) else ("normal", "tumour"))
        counts = H.is_counts(x)
        for status, y in (("uninvolved", ref), ("lesional", case)):
            rows.append(dict(context=context, unit=f"{name}:{k}", tissue=tissue(k), status=status,
                             **coupling(y, counts, wfs)))
    return rows


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    out6 = H.out_dir(H.FOLDER, run["quick"])
    n_boot = H.sizes(args, run, 10000, 10000)["n_boot"]
    wfs = {k: v for k, v in Q.cue_waveforms().items() if k in ("core", "cue", "locked")}
    g = gtex_reference(wfs, args.subsample, args.seed)
    g.to_csv(out / "a88_gtex.csv", index=False)
    print(g.round(3).to_string(), flush=True)
    genes = Q.all_set_genes(wfs)
    rows = []
    pso, w = H.gated_samples(out6, args.subsample, args.seed)
    rows += unit_rows("outpatient_skin", "skin_psoriasis", lambda k: "skin", pso, w, wfs)
    disc, _ = D.discovery_all(args.subsample, args.seed)
    ad = disc[disc.disease == "atopic_dermatitis"]
    s, w = Q.gated_cohorts(ad, pd.read_csv(out / "a53_skin_per_series.csv").query("disease == 'atopic_dermatitis'"))
    rows += unit_rows("outpatient_skin", "skin_atopic", lambda k: "skin", s, w, wfs)
    ibd, _ = D.ibd(args.subsample, args.seed)
    s, w = Q.gated_cohorts(ibd, pd.read_csv(out / "a55_ibd_per_series.csv"))
    rows += unit_rows("endoscopy", "ibd", lambda k: "colon", s, w, wfs)
    pc = cancer_table(run["quick"])
    if pc is not None:
        t, w = Q.tcga_table(pc, genes, args.subsample, args.seed)
        rows += unit_rows("surgery", "tcga", lambda k: k, t, w, wfs)
    u = pd.DataFrame(rows)
    ref = g[g.context == "gtex_fast"].set_index("unit").log_rho
    u["ref"] = u.tissue.map(lambda t: ref.reindex(MATCH[t]).mean() if t in MATCH else np.nan)
    u["delta"] = u.log_rho - u.ref
    u.to_csv(out / "a88_units.csv", index=False)
    x = u[(u.status == "uninvolved") & u.delta.notna()]
    sg = x[x.context.isin(["surgery", "endoscopy"])].delta.to_numpy()
    sk = x[x.context == "outpatient_skin"].delta.to_numpy()
    rng = np.random.default_rng(SEED)
    diff = sg.mean() - sk.mean()
    boot = np.array([rng.choice(sg, len(sg)).mean() - rng.choice(sk, len(sk)).mean() for _ in range(n_boot)])
    v = dict(surgical_endoscopy=sg.mean(), outpatient_skin=sk.mean(), diff=diff, lo=np.quantile(boot, 0.025),
             hi=np.quantile(boot, 0.975), p_lower=float((1 + np.sum(boot >= 0)) / (1 + n_boot)), n_surgical=len(sg),
             n_skin=len(sk))
    v["H88"] = "PASS" if (diff < 0 and v["p_lower"] < 0.05) else "FAIL"
    dig = x[x.tissue.isin(["COAD", "LIHC", "colon"]) & x.context.isin(["surgery", "endoscopy"])].delta
    oth = x[~x.tissue.isin(["COAD", "LIHC", "colon"]) & (x.context == "surgery")].delta
    v.update(digestive=dig.mean(), other_surgical=oth.mean())
    pd.Series(v).to_csv(out / "a88_verdict.csv", header=False)
    print(u.round(3).to_string(), "\n", pd.Series(v).to_string())


if __name__ == "__main__":
    main()
