"""Post hoc descriptive check (after A342 / A343, user request 2026-11-22): why do subtraction and direct intron counts
give different pre-mRNA leads? Per clock gene, in timed skin (GSE326530) and in the dexamethasone set (GSE320593):
total subtraction and direct intron counts, the direct / subtraction ratio, the annotation's overlap flags and
overlapping genes, and intron-only length (stranded vs unstranded regions). In skin, per gene and site, the timed lag of
each intron definition over exon (first-harmonic acrophase of log CPM, volunteer means removed; positive = intron
ahead). No hypothesis, no verdict.

    python intron_definitions.py
Outputs (`../output/intron_definitions/`): genes.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import power_a342 as A  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

ROOT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/sccc_mapping")
OUT = HERE.parents[1] / "output" / "intron_definitions"
HOURS = 24 / (2 * np.pi)


def load(ds, name):
    x = pd.read_csv(ROOT / ds / f"counts_{name}.tsv.gz", sep="\t")
    x["symbol"] = x.symbol.astype(str).map(canonical)
    return x.drop_duplicates("symbol").set_index("symbol").drop(columns="gene_id")


def acro(y, t, vol):
    y = y - pd.Series(y).groupby(vol).transform("mean").to_numpy()
    b = np.mean(y * np.exp(-1j * t / HOURS))
    return np.angle(b) * HOURS, 2 * abs(b)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    d = A.design()
    for ds in ("gse326530_skin", "gse320593_liver_fat_dex"):
        ex, sub, dr = load(ds, "exon"), load(ds, "intron"), load(ds, "intron_direct")
        ann = pd.read_csv(ROOT / ds / "gene_annotation.tsv.gz", sep="\t")
        ann["symbol"] = ann.symbol.astype(str).map(canonical)
        ann = ann.drop_duplicates("symbol").set_index("symbol")
        for g in A.GENES:
            r = dict(dataset=ds, gene=g, exon=int(ex.loc[g].sum()), intron_sub=int(sub.loc[g].sum()),
                     intron_direct=int(dr.loc[g].sum()))
            r["direct_over_sub"] = r["intron_direct"] / max(r["intron_sub"], 1)
            a = ann.loc[g]
            r.update(genebody_kb=a.genebody_len_bp / 1e3, intron_kb_stranded=a.intron_only_len_bp_stranded / 1e3,
                     intron_kb_unstranded=a.intron_only_len_bp_unstranded / 1e3,
                     overlaps_same_strand=a.overlaps_same_strand, overlaps_any_strand=a.overlaps_any_strand,
                     overlapping_any=str(a.overlapping_genes_any_strand)[:120])
            if ds == "gse326530_skin":
                for site in ("P", "E"):
                    dd = d[d.site == site]
                    t, vol = dd.t.to_numpy(float), dd.vol.to_numpy()
                    dep = ex.loc[:, dd.run].sum().to_numpy(float)
                    ph = {}
                    for k, X in (("exon", ex), ("sub", sub), ("direct", dr)):
                        y = np.log((X.loc[g, dd.run].to_numpy(float) + 0.5) / dep * 1e6)
                        ph[k], r[f"{site}_amp_{k}"] = acro(y, t, vol)
                    r[f"{site}_lag_sub_h"] = A.P.score_a337.wrap(ph["sub"] - ph["exon"])
                    r[f"{site}_lag_direct_h"] = A.P.score_a337.wrap(ph["direct"] - ph["exon"])
            rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "genes.csv", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 60)
    print(t.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
