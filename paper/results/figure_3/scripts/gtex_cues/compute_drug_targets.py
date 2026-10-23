"""A100: are rhythmic drug targets clock-locked or cue-driven, and how much of their wall-time rhythm survives in hospital?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A100). Targets: Human Protein Atlas class "FDA approved drug targets"
($SCCC_DATA/drug_targets/raw/hpa_fda_targets.tsv, sha256 in its PROVENANCE.md). Gene x tissue table: the 26-tissue selection
table (rhythm in fast deaths, ventilated lock) and the A91 classes; tissues that pass the A89 noise check (>= 50 rhythmic
calls above the 1 % chance count); rhythmic pair = p_fast < 0.01 (clock-plane and D6 clock-output genes excluded). New
quantity: wall-time retention = timed cosinor amplitude on time of death (D2 covariates) in ventilated donors (and Hardy 3 + 4)
over fast deaths. H100: share of cue-like pairs (lock < 0.5) among rhythmic pairs, targets minus other genes, pooled over
tissues; target labels permuted among genes (a gene keeps one label across tissues, 2,000); two-sided. Descriptive: shares
locked (lock >= 0.8) / cue-like / intermediate, median retention per class, A91 class shares.
Needs compute_internal_external.py (A91 classes). Aggregates only.

    python compute_drug_targets.py [--subsample N] [--jobs 2]
Outputs: a100_summary.csv, a100_targets.csv (rhythmic targets with tissues, lock, retention), a100_rhythmic_pairs.csv.gz
(gene x tissue aggregates, reused by compute_covariate_sensitivity.py).
Legacy: legacy/desync/validation/paired_damping/scripts/a100_drug_targets.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

import _genes as G  # noqa: E402
import _shared as S  # noqa: E402
from _genes import FOLDER, W  # noqa: E402

SEED = 100
TARGETS = DATA / "drug_targets" / "raw" / "hpa_fda_targets.tsv"
GROUPS = {"vent": (0,), "H34": (3, 4), "fast": (1, 2)}


def timed_amp(Y, ob):
    t = ob.tod.to_numpy(float) % 24
    Zc = S.d2_covariates(ob)
    B, _ = S.lstsq_fit(Y, np.c_[np.cos(W * t), np.sin(W * t)], Zc[:, Zc.std(0) > 0])
    return np.abs(B[-2] + 1j * B[-1])


def retention(tis, genes, sub):
    amp = {}
    for grp, hardy in GROUPS.items():
        X, g, ob = G.load_rows(tis, hardy, sub, genes=genes)
        amp[grp] = (timed_amp(np.log(X + 0.5), ob) if len(X) >= 40 else np.full(len(g), np.nan), g, len(X))
    gf = amp["fast"][1]
    return pd.DataFrame(dict(tissue=tis, gene=gf, n_vent=amp["vent"][2], n_h34=amp["H34"][2],
                             ret_vent=amp["vent"][0] / amp["fast"][0], ret_h34=amp["H34"][0] / amp["fast"][0]))


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = 200 if run["quick"] else 2000
    o = G.obs_table()
    sub = G.pick_donors(o, args.subsample, args.seed)
    targets = set(pd.read_csv(TARGETS, sep="\t").Gene.astype(str))
    g = G.selection_tables(out, args.subsample, args.seed, args.jobs)
    g = g[~g.gene.isin(G.EXCLUDE)]
    excess = g.groupby("tissue").apply(lambda d: (d.p_fast < 0.01).sum() - 0.01 * len(d), include_groups=False)
    keep_t = excess.index[excess >= 50]
    r = g[g.tissue.isin(keep_t) & (g.p_fast < 0.01)].copy()
    a91 = pd.read_csv(out / "a91_genes.csv.gz", usecols=["tissue", "gene", "cls"])
    r = r.merge(a91, on=["tissue", "gene"], how="left")
    r["target"] = r.gene.isin(targets)
    r["lock_class"] = np.select([r.lock >= 0.8, r.lock < 0.5], ["locked", "cue-like"], "intermediate")
    jobs = [(t, x.gene.unique(), sub) for t, x in r.groupby("tissue")]
    ret = pd.concat(Parallel(n_jobs=args.jobs)(delayed(retention)(*j) for j in jobs), ignore_index=True)
    r = r.merge(ret, on=["tissue", "gene"], how="left")
    cue = (r.lock_class == "cue-like").to_numpy()
    genes = r.gene.to_numpy()
    uniq = np.unique(genes)
    is_t = np.isin(uniq, list(targets))
    pos = np.searchsorted(uniq, genes)
    stat = lambda lab: cue[lab[pos]].mean() - cue[~lab[pos]].mean()
    obs = stat(is_t)
    rng = np.random.default_rng(SEED)
    null = np.array([stat(rng.permutation(is_t)) for _ in range(n_perm)])
    v = dict(tissues=len(keep_t), rhythmic_pairs=len(r), target_pairs=int(r.target.sum()),
             target_genes=int(r[r.target].gene.nunique()), cue_share_targets=float(cue[r.target].mean()),
             cue_share_others=float(cue[~r.target].mean()), H100_diff=float(obs), n_perm=n_perm,
             H100_p=cs.stats.perm_p(null, obs, "two-sided", tie=0.0))
    v["H100"] = "PASS" if v["H100_p"] < 0.05 else "FAIL"
    for who, x in (("targets", r[r.target]), ("others", r[~r.target])):
        for c in ("locked", "intermediate", "cue-like"):
            y = x[x.lock_class == c]
            v[f"{who}_{c}_share"] = len(y) / len(x)
            v[f"{who}_{c}_ret_vent_median"] = float(y.ret_vent.median())
            v[f"{who}_{c}_ret_h34_median"] = float(y.ret_h34.median())
        v[f"{who}_ret_vent_median"] = float(x.ret_vent.median())
        v[f"{who}_ret_h34_median"] = float(x.ret_h34.median())
        for c in ("clock", "time", "both"):
            v[f"{who}_a91_{c}_share"] = float((x.cls == c).mean())
    r.to_csv(out / "a100_rhythmic_pairs.csv.gz", index=False)
    tab = r[r.target].groupby("gene").agg(tissues=("tissue", lambda s: ";".join(sorted(s))), n=("tissue", "size"),
                                          median_lock=("lock", "median"), median_amp=("amp_fast", "median"),
                                          ret_vent=("ret_vent", "median"), ret_h34=("ret_h34", "median"))
    tab.sort_values("n", ascending=False).to_csv(out / "a100_targets.csv")
    pd.Series(v).to_csv(out / "a100_summary.csv", header=False)
    print(pd.Series(v).to_string(), flush=True)


if __name__ == "__main__":
    main()
