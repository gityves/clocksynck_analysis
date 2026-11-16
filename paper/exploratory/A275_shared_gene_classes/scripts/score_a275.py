"""A275: which genes are shared between the cell types of a donor? Per-gene sharing of clock genes, glucocorticoid
targets, immediate-early genes and expression-matched genes (the A272 atlases).

Registered in paper/docs/registered_experiments.md (A275; descriptive decomposition of the responsive companion of
A272). For a gene g and a pair of compartments (a, b) of an atlas: the Pearson correlation over donors of
log((count + 0.5) / depth), centred within stratum and pooled over strata (the A272 pairs). Sharing of a gene in a pair
class = median over the pairs of the class; of a gene class = median over its genes of the per-atlas value, then the
mean over atlases.

    python score_a275.py [--n-matched 200]
Outputs: paper/exploratory/A275_shared_gene_classes/output/{genes,classes}.csv
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync._core.controls import GC_TARGETS, IMMEDIATE_EARLY, pseudo_clock_genes  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

spec = importlib.util.spec_from_file_location("a275_a272", PAPER / "exploratory/A272_backbone_pairwise/scripts/score_a272.py")
A272 = importlib.util.module_from_spec(spec)
sys.modules["a275_a272"] = A272
spec.loader.exec_module(A272)
A150 = A272.A150
AMENDMENT, SEED = "A275_shared_gene_classes", 275


def gene_corr(Y, pr):
    """Pooled within-stratum correlation of every gene (columns of Y) between the two compartments of a pair."""
    s = pr["s"]
    ns = np.bincount(s)
    def centre(M):
        m = np.zeros((s.max() + 1, M.shape[1]))
        np.add.at(m, s, M)
        return M - (m / ns[:, None])[s]
    a, b = centre(Y[pr["ia"]]), centre(Y[pr["ib"]])
    den = np.sqrt((a ** 2).sum(0) * (b ** 2).sum(0))
    return np.divide((a * b).sum(0), den, out=np.full(Y.shape[1], np.nan), where=den > 0)


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--n-matched", type=int, default=200)
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(AMENDMENT)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for k, n in enumerate(A272.ATLASES):
        A = A150.load_atlas(n, None, 0)
        pb = A["pb"].reset_index(drop=True)
        pairs = A272.pair_index(pb)
        genes = A["genes"]
        X = pb[genes].to_numpy(float)
        frac = (X / pb.depth.to_numpy(float)[:, None]).mean(0)
        keep = (X.mean(0) >= 1) | np.isin(genes, A272.CLOCK)
        lm = pd.Series(np.log(frac + 1e-12), index=genes)[keep]
        sets = pseudo_clock_genes(lm, args.n_matched, np.random.default_rng([SEED, k]), genes=A272.CLOCK)
        cls = {g: "clock" for g in A272.CLOCK}
        cls.update({g: "gc_target" for g in GC_TARGETS if g in lm.index and g not in cls})
        cls.update({g: "immediate_early" for g in IMMEDIATE_EARLY if g in lm.index and g not in cls})
        for g in {g for s in sets for g in s}:
            cls.setdefault(g, "matched")
        use = [g for g in cls if g in pb.columns]
        Y = np.log((pb[use].to_numpy(float) + 0.5) / pb.depth.to_numpy(float)[:, None])
        kl = [A272.klass(A272.group(A150.FAM.family(p["a"])), A272.group(A150.FAM.family(p["b"]))) for p in pairs]
        C = np.array([gene_corr(Y, p) for p in pairs])  # pairs x genes
        for c in sorted(set(kl)):
            m = np.array(kl) == c
            med = np.nanmedian(C[m], 0)
            rows += [dict(atlas=n, klass=c, gene=g, gene_class=cls[g], pairs=int(m.sum()), sharing=float(v),
                          mean_log_fraction=float(lm.get(g, np.nan))) for g, v in zip(use, med)]
        A272.log(n, "genes", len(use), "pairs", len(pairs))
    G = pd.DataFrame(rows)
    G.to_csv(out / "genes.csv", index=False)
    per = G.groupby(["klass", "gene_class", "atlas"]).sharing.median().reset_index()
    S = per.groupby(["klass", "gene_class"]).agg(sharing=("sharing", "mean"), atlases=("atlas", "nunique")).reset_index()
    S["genes"] = [G[(G.klass == r.klass) & (G.gene_class == r.gene_class)].gene.nunique() for r in S.itertuples()]
    S.to_csv(out / "classes.csv", index=False)
    print(S.pivot(index="gene_class", columns="klass", values="sharing").round(3).to_string())
    top = G[G.klass == "SV_X"].groupby(["gene_class", "gene"]).sharing.mean().reset_index()
    print(top[top.gene_class != "matched"].sort_values("sharing", ascending=False).round(3).to_string())


if __name__ == "__main__":
    main()
