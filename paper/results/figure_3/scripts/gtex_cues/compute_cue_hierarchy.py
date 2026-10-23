"""A89 and A90: from body-wide cues to organ-specific rhythms, and are p53 targets clock-locked or cue-like? (GTEx)

Pre-registered in docs/PLAN_internal_desync.md (archive) (A89, A90).
A89, descriptive plus one test (H89). Gene x tissue calls from the 26-tissue selection table (rhythm in all fast deaths,
lock from the ventilated selection half; clock-plane and D6 clock-output genes excluded): cue-like = p_fast < 0.01 and
lock < 0.5; clock-like = p_fast < 0.01 and lock >= 0.8. Levels per gene: body-wide (cue-like in >= 3 a priori tissue
groups), group-specific (>= 2 tissues of one group, none outside), organ-specific (one tissue, p_fast > 0.1 in >= 80 % of
the other tissues where expressed). Phase coherence R of acrophases per group (cue-like vs clock-like genes); Hallmark
enrichment (hypergeometric, BH within level) when a Hallmark gmt is present under $SCCC_DATA/reference; tissue clustering
(average linkage on the Jaccard overlap of cue-like sets, 6 clusters). H89: mean within-group R of cue-like genes minus
that of random tissue sets of the same sizes (1,000), one-sided.
A90: per p53 target and tissue: lock with 95 % intervals from 200 resamples of the ventilated selection donors (host fitted
once); locked if the interval lies above 0.6 in most rhythmic tissues (p_fast < 0.05), cue-like if below 0.6, else undetermined.
Hosts: `clocksync` `Dataset.fit(...).host`. Aggregates only.

    python compute_cue_hierarchy.py [--subsample N] [--jobs 2]
Outputs: a89_levels.csv, a89_tissues.csv, a89_group_coherence.csv, a89_enrichment.csv, a89_jaccard.csv, a89_verdict.csv,
a90_p53_lock.csv, a90_summary.csv (and the cached selection tables of compute_cue_reference.py).
Legacy: legacy/desync/validation/paired_damping/scripts/a89_cue_hierarchy.py, a90_p53_lock.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.cluster.hierarchy import fcluster, linkage  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _genes as G  # noqa: E402
import _shared as S  # noqa: E402
from _genes import FOLDER, W  # noqa: E402

GROUPS = {
    "digestive": ["Colon - Transverse", "Colon - Sigmoid", "Stomach", "Esophagus - Mucosa", "Esophagus - Muscularis",
                  "Esophagus - Gastroesophageal Junction", "Liver"],
    "barrier": ["Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)", "Esophagus - Mucosa", "Lung"],
    "endocrine": ["Thyroid", "Adrenal Gland"],
    "cardiovascular_muscle": ["Heart - Atrial Appendage", "Heart - Left Ventricle", "Artery - Aorta",
                              "Artery - Coronary", "Artery - Tibial", "Muscle - Skeletal"],
    "adipose_breast": ["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Breast - Mammary Tissue"],
    "other": ["Nerve - Tibial", "Prostate", "Testis", "Kidney - Cortex", "Uterus"],
}
P53 = ["TP53", "CDKN1A", "MDM2", "TP53I3", "TP53I11", "GDF15", "SESN1", "SESN2", "ZMAT3", "RRM2B", "BAX", "FAS",
       "TNFRSF10B", "PLK3", "BTG2", "DDB2", "XPC", "TP53INP1", "AEN", "PHLDA3"]
SEED_A89, SEED_A90, N_BOOT = 89, 90, 200


# ---- A89 ---------------------------------------------------------------------------------------------------------

def group_R(g: pd.DataFrame, tissues) -> float:
    """Mean over genes (cue-like in >= 2 of the tissues) of the phase coherence of their acrophases."""
    x = g[g.tissue.isin(tissues) & g.cue]
    r = x.groupby("gene").acro_h.agg(["size", G.circ_R])
    r = r[r["size"] >= 2]
    return float(r.circ_R.mean()) if len(r) else np.nan


def levels(g: pd.DataFrame) -> pd.DataFrame:
    grp_of = {}
    for k, ts in GROUPS.items():
        for t in ts:
            if k != "other":
                grp_of.setdefault(t, set()).add(k)
    rows = []
    for gene, x in g[g.cue].groupby("gene"):
        all_t = g[g.gene == gene]
        if len(set().union(*[grp_of.get(t, set()) for t in x.tissue])) >= 3:
            level = "body-wide"
        elif len(x) == 1:
            others = all_t[all_t.tissue != x.tissue.iloc[0]]
            level = "organ-specific" if len(others) and (others.p_fast > 0.1).mean() >= 0.8 else "unclassified"
        else:
            level = "unclassified"
            for k, ts in GROUPS.items():
                if k != "other" and x.tissue.isin(ts).sum() >= 2 and x.tissue.isin(ts).all():
                    level = f"group:{k}"
        rows.append(dict(gene=gene, level=level, tissues=";".join(sorted(x.tissue)), n=len(x),
                         R=G.circ_R(x.acro_h) if len(x) >= 2 else np.nan,
                         acro_h=float(np.angle(np.exp(1j * W * x.acro_h).mean()) / W % 24),
                         amp=x.amp_fast.median(), lock=x.lock.median()))
    return pd.DataFrame(rows)


def enrichment(lv: pd.DataFrame, universe: set, H: dict) -> pd.DataFrame:
    rows = []
    for level, x in lv.groupby("level"):
        s = set(x.gene) & universe
        for name, hs in H.items():
            hs = hs & universe
            k = len(s & hs)
            if k:
                rows.append(dict(level=level, set=name, k=k, n=len(s), K=len(hs),
                                 p=hypergeom.sf(k - 1, len(universe), len(hs), len(s))))
    e = pd.DataFrame(rows)
    if len(e):
        e["q"] = e.groupby("level").p.transform(lambda p: np.minimum(1, p * len(p) / p.rank()))
    return e.sort_values("p") if len(e) else e


def a89(g26: pd.DataFrame, out: Path, n_perm: int) -> None:
    g = g26[~g26.gene.isin(G.EXCLUDE)].copy()
    g["cue"] = (g.p_fast < 0.01) & (g.lock < 0.5)
    g["clock"] = (g.p_fast < 0.01) & (g.lock >= 0.8)
    tissues = sorted(g.tissue.unique())
    lv = levels(g)
    lv.to_csv(out / "a89_levels.csv", index=False)
    counts = lv.level.value_counts()
    per_tissue = g.groupby("tissue")[["cue", "clock"]].sum().join(g.groupby("tissue").size().rename("expressed"))
    coh = pd.DataFrame({k: dict(cue=group_R(g, ts), clock=group_R(g.assign(cue=g.clock), ts))
                        for k, ts in GROUPS.items()}).T
    sizes = [len(ts) for k, ts in GROUPS.items() if k != "other"]
    obs = np.nanmean([group_R(g, ts) for k, ts in GROUPS.items() if k != "other"])
    rng = np.random.default_rng(SEED_A89)
    null = np.array([np.nanmean([group_R(g, list(rng.choice(tissues, s, replace=False))) for s in sizes])
                     for _ in range(n_perm)])
    v = dict(within_group_R=obs, random_R=float(np.nanmean(null)), n_perm=n_perm,
             p_upper=cs.stats.perm_p(null, obs, "greater", tie=0.0))
    v["H89"] = "PASS" if (obs > np.nanmean(null) and v["p_upper"] < 0.05) else "FAIL"
    sets = {t: set(g[(g.tissue == t) & g.cue].gene) for t in tissues}
    J = np.array([[len(sets[a] & sets[b]) / max(1, len(sets[a] | sets[b])) for b in tissues] for a in tissues])
    clus = pd.Series(fcluster(linkage(squareform(1 - J, checks=False), "average"), 6, "maxclust"), index=tissues,
                     name="cluster")
    H = G.hallmarks()
    enr = enrichment(lv, set(g[g.p_fast < 0.01].gene), H) if H else pd.DataFrame()
    per_tissue.join(clus).to_csv(out / "a89_tissues.csv")
    coh.to_csv(out / "a89_group_coherence.csv")
    enr.to_csv(out / "a89_enrichment.csv", index=False)
    pd.DataFrame(J, index=tissues, columns=tissues).to_csv(out / "a89_jaccard.csv")
    pd.Series(v).to_csv(out / "a89_verdict.csv", header=False)
    print(counts.to_string(), "\n", coh.round(3).to_string(), "\n", pd.Series(v).to_string(), flush=True)


# ---- A90 ---------------------------------------------------------------------------------------------------------

def p53_one(tis, fv, fast_sub, vent_sub, seed) -> pd.DataFrame:
    genes = P53 + list(G.CORE_GENES)
    Xf, gf, obf = G.load_rows(tis, (1, 2), fast_sub, genes=genes)
    Xv, gv, obv = G.load_rows(tis, (0,), vent_sub, genes=genes)
    Yf, Yv = np.log(Xf + 0.5), np.log(Xv + 0.5)
    c_fast, p_fast = S.timed_cosinor(Yf, obf.tod.to_numpy(float) % 24, S.d2_covariates(obf))
    subj = obv.subject.to_numpy()
    m = np.isin(subj, fv.index)
    Yv, obv, f = Yv[m], obv[m], fv[["x", "y"]].loc[subj[m]].to_numpy()
    Zc = S.d2_covariates(obv)
    clock = np.isin(gf, list(G.CORE_GENES))

    def lock(idx):
        ff = f[idx] - f[idx].mean(0)
        Zi = Zc[idx][:, Zc[idx].std(0) > 0]
        B, _ = S.lstsq_fit(Yv[idx], ff, Zi)
        r = np.abs(B[-2] + 1j * B[-1]) / np.abs(c_fast)
        ref = clock & (p_fast < 0.01)
        return r / np.median(r[ref]) if ref.any() else np.full(len(r), np.nan)

    n = len(Yv)
    obs = lock(np.arange(n))
    rng = np.random.default_rng(seed)
    bs = np.array([lock(rng.integers(0, n, n)) for _ in range(N_BOOT)])
    return pd.DataFrame(dict(tissue=tis, gene=gf, p_fast=p_fast, amp_fast=np.abs(c_fast), acro_h=np.angle(c_fast) / W % 24,
                             lock=obs, lo=np.nanquantile(bs, 0.025, axis=0), hi=np.nanquantile(bs, 0.975, axis=0)))


def a90(out: Path, H, sub, jobs: int) -> None:
    fast = H["fast"]["sel"] | H["fast"]["eval"]
    ts = G.A81_TISSUES + G.A84_ADDED + G.A89_EXTRA  # as the legacy run: A81 tissues incl. whole blood (lock undefined)
    hosts = G.hosts_for(ts, (0,), donors=H["vent"]["sel"])
    res = pd.concat(Parallel(n_jobs=jobs)(delayed(p53_one)(t, hosts[t], fast, H["vent"]["sel"], SEED_A90 + i)
                                          for i, t in enumerate(ts)), ignore_index=True)
    res.to_csv(out / "a90_p53_lock.csv", index=False)
    r = res[res.gene.isin(P53) & (res.p_fast < 0.05)].copy()
    r["call"] = np.where(r.lo > 0.6, "locked", np.where(r.hi < 0.6, "cue-like", "undetermined"))
    summ = r.groupby("gene").agg(rhythmic_tissues=("tissue", "size"), median_lock=("lock", "median"),
                                 locked=("call", lambda x: (x == "locked").sum()),
                                 cue_like=("call", lambda x: (x == "cue-like").sum()),
                                 undetermined=("call", lambda x: (x == "undetermined").sum()))
    summ["class"] = np.where(summ.locked > summ.rhythmic_tissues / 2, "locked",
                             np.where(summ.cue_like > summ.rhythmic_tissues / 2, "cue-like", "undetermined"))
    summ.to_csv(out / "a90_summary.csv")
    print(summ.sort_values("median_lock").round(2).to_string(), flush=True)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    o = G.obs_table()
    sub = G.pick_donors(o, args.subsample, args.seed)
    H = G.halves(o, sub)
    g26 = G.selection_tables(out, args.subsample, args.seed, args.jobs)
    a89(g26, out, 50 if run["quick"] else 1000)
    a90(out, H, sub, args.jobs)


if __name__ == "__main__":
    main()
