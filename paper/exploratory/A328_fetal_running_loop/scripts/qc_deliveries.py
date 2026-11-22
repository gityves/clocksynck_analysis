"""Blind-safe QC of the Jed mapping deliveries ($SCCC_DATA/sccc_mapping), before the A328 power stage.

Reads no clock-gene value: per dataset it reports the design resolved from the metadata, mapping and intron depth
(qc.tsv), whether the 14 clock genes and the 11 core genes have a row (presence only), and, where a processed matrix of
the same libraries is already on disk, the all-gene agreement of the delivered exon counts with it (Spearman of log
CPM over shared genes; each delivered library matched to its best processed column). del Valle has no processed
matrix on disk and stays blind (A328 part a).

`--cells` (A328 note 2, task 1): the 10x deliveries only. Spliced sums over the A244 cells per fetus against every
A242 whole-adrenal pseudobulk (Spearman over all genes except the clock and core genes; identity mismatch = best match
is another fetus; flag = own r < 0.9), and the genome-wide unspliced fraction per GSM and class (totals over all genes).

Output: paper/exploratory/A328_fetal_running_loop/output/qc/{summary.csv, design.csv, agreement.csv};
`--cells`: {agreement_10x_a242.csv, unspliced_fraction_10x.csv}.
Usage: python qc_deliveries.py [--cells]
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from clocksync.genes import CLOCK_GENES, canonical

DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
SRC = DATA / "sccc_mapping"
OUT = Path(__file__).resolve().parents[1] / "output" / "qc"
CORE = ("ARNTL", "CLOCK", "NPAS2", "PER3", "CRY1", "CRY2", "NR1D2", "DBP", "RORA", "TEF", "CIART")
MOUSE = {"ARNTL": "ARNTL", "BMAL1": "ARNTL", "CIART": "CIART", "GM129": "CIART"}


def read_counts(d: Path, kind: str) -> pd.DataFrame:
    t = pd.read_csv(d / f"counts_{kind}.tsv.gz", sep="\t", index_col=0)
    return t


def log_cpm(m: pd.DataFrame) -> pd.DataFrame:
    return np.log2(m / m.sum() * 1e6 + 1)


def agreement(new: pd.DataFrame, old: pd.DataFrame, expected: dict | None) -> pd.DataFrame:
    """Best processed column per delivered library, Spearman over shared genes expressed in both."""
    g = new.index.intersection(old.index)
    a, b = log_cpm(new.loc[g]), np.log2(old.loc[g].clip(lower=0) + 1)
    keep = (a.mean(1) > 1) & (b.mean(1) > 1)
    ra = np.apply_along_axis(rankdata, 0, a[keep].to_numpy())
    rb = np.apply_along_axis(rankdata, 0, b[keep].to_numpy())
    ra, rb = (ra - ra.mean(0)) / ra.std(0), (rb - rb.mean(0)) / rb.std(0)
    r = ra.T @ rb / len(ra)
    rows = []
    for i, run in enumerate(new.columns):
        o = np.argsort(r[i])[::-1]
        best = old.columns[o[0]]
        rows.append({"run": run, "best": best, "r_best": r[i, o[0]], "r_second": r[i, o[1]],
                     "expected": expected.get(run) if expected else None,
                     "match": (best == expected.get(run)) if expected else None, "genes": int(keep.sum())})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summary, designs, agree = [], [], []
    for d in sorted(p for p in SRC.iterdir() if p.is_dir() and not p.name.startswith(".")):
        s = pd.read_csv(d / "samples.tsv", sep="\t")
        q = pd.read_csv(d / "qc.tsv", sep="\t")
        exon = read_counts(d, "exon")
        sym = exon["symbol"].astype(str).map(canonical)
        present = set(sym)
        row = {"dataset": d.name, "runs": len(s), "genes": len(exon),
               "clock14_present": sum(g in present for g in CLOCK_GENES),
               "core11_present": sum(g in present for g in CORE),
               "clock_missing": ";".join(g for g in CLOCK_GENES if g not in present)}
        for c in ("input_reads", "uniquely_mapped_pct", "assigned_exon_pct", "intronic_fraction",
                  "intron_direct_fraction"):
            row[f"{c}_median"] = q[c].median()
            row[f"{c}_min"] = q[c].min()
            row[f"{c}_max"] = q[c].max()
        row["strand"] = ";".join(f"{k}:{v}" for k, v in q["strand_used"].value_counts().items())
        row["ambiguous_direct_pct_median"] = (100 * q["intron_direct_ambiguous_n"]
                                              / (q["intron_direct_assigned_n"] + q["intron_direct_ambiguous_n"])).median()
        summary.append(row)

        counts = exon.drop(columns="symbol")
        if d.name.startswith("delvalle"):
            g = s.groupby(["sdrf_Characteristics_organism_part", "sdrf_Characteristics_age",
                           "sdrf_Characteristics_sex"]).size()
            designs += [{"dataset": d.name, "group": " | ".join(map(str, k)), "runs": v} for k, v in g.items()]
            designs.append({"dataset": d.name, "group": "distinct sample names / runs",
                            "runs": f"{s['sample_title'].nunique()} / {len(s)}"})
        elif d.name.startswith("cardoso"):
            g = s.groupby(["sdrf_Characteristics_organism_part", "sdrf_Characteristics_age"]).size()
            designs += [{"dataset": d.name, "group": " | ".join(map(str, k)), "runs": v} for k, v in g.items()]
            designs.append({"dataset": d.name, "group": "individuals",
                            "runs": s["sdrf_Characteristics_individual"].nunique()})
            cm = DATA / "natexp/cardoso_moreira2019"
            old = pd.read_parquet(cm / "expr_cpm.parquet")
            new = counts.copy()
            new.index = new.index.str.split(".").str[0]
            exp = {r: f"{o.capitalize()}.{a}wpc" for r, o, a in zip(s["run_accession"],
                   s["sdrf_Characteristics_organism_part"], s["sdrf_Characteristics_age"])}
            t = agreement(new, old, None)
            t["expected"] = t["run"].map(exp)
            t["match"] = [b.rsplit(".", 1)[0] == e for b, e in zip(t["best"], t["expected"])]
            agree.append(t.assign(dataset=d.name))
        elif d.name.startswith("mure"):
            g = s.groupby("char_tissue").size()
            designs += [{"dataset": d.name, "group": k, "runs": v} for k, v in g.items()]
            old = pd.read_csv(DATA / "baboon/GSE98965/GSE98965_baboon_tissue_expression_FPKM.csv.gz", index_col=0)
            old = old.drop(columns="Symbol")
            new = counts.copy()
            new.index = new.index.str.split(".").str[0]
            exp = dict(zip(s["run_accession"], s["geo_title"].str.replace("_", ".")))
            agree.append(agreement(new, old, exp).assign(dataset=d.name))
        elif d.name.startswith("greenwell"):
            g = s.groupby("geo_source_name_ch1").size()
            designs += [{"dataset": d.name, "group": k, "runs": v} for k, v in g.items()]
            old = pd.read_csv(DATA / "greenwell2019/GSE118967_RAW_Table.txt.gz", sep="\t", index_col=0)
            old.index = old.index.str.split(".").str[0]
            new = counts.copy()
            new.index = new.index.str.split(".").str[0]
            # titles AR1ZT02 / DA1ZT02 / ... vs table columns A1ZT02 / D1ZT02 / ...: first letter + replicate + ZT
            exp = {r: t[0] + t[-5:] for r, t in zip(s["run_accession"], s["geo_title"])}
            agree.append(agreement(new, old, exp).assign(dataset=d.name))

    pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)
    pd.DataFrame(designs).to_csv(OUT / "design.csv", index=False)
    a = pd.concat(agree, ignore_index=True)
    a.to_csv(OUT / "agreement.csv", index=False)
    pd.set_option("display.width", 200)
    print(pd.DataFrame(summary).T.to_string())
    print(pd.DataFrame(designs).to_string(index=False))
    print(a.groupby("dataset").agg(r_best_median=("r_best", "median"), r_best_min=("r_best", "min"),
                                   gap_min=("r_second", lambda x: np.nan), matched=("match", "mean"),
                                   genes=("genes", "first")).to_string())
    for ds, t in a.groupby("dataset"):
        print(ds, "min r_best - r_second:", round((t["r_best"] - t["r_second"]).min(), 4),
              "| best column used twice:", int(t["best"].duplicated().sum()))
        if t["match"].notna().any():
            print(t.loc[t["match"] == False, ["run", "best", "expected", "r_best"]].head(10).to_string(index=False))


ADRENAL_10X = ("kameneva2021_fetal_adrenal", "wang2022_fetal_adrenal", "cheng2022_fetal_adrenal",
               "dong2020_fetal_adrenal")
ORGAN_10X = ("knightschrijver2023_fetal_heart", "hochane2019_fetal_kidney")


def spearman_matrix(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    """Spearman of log CPM between the rows of a and b (samples x genes) over genes expressed in both."""
    la, lb = (np.log2(m.div(m.sum(1), axis=0) * 1e6 + 1) for m in (a, b))
    keep = (la.mean() > 1) & (lb.mean() > 1)
    ra = np.apply_along_axis(rankdata, 1, la.loc[:, keep].to_numpy())
    rb = np.apply_along_axis(rankdata, 1, lb.loc[:, keep].to_numpy())
    ra = (ra - ra.mean(1, keepdims=True)) / ra.std(1, keepdims=True)
    rb = (rb - rb.mean(1, keepdims=True)) / rb.std(1, keepdims=True)
    return pd.DataFrame(ra @ rb.T / keep.sum(), index=a.index, columns=b.index).assign(genes=int(keep.sum()))


def main_cells():
    """A328 note 2 QC (task 1): identity agreement of the delivered spliced sums over the A244 cells with the A242
    whole-adrenal pseudobulks (clock genes and the core genes excluded from the genes compared), and the genome-wide
    unspliced fraction per GSM and class from the derived cells tables (totals over all genes; no gene-level value)."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from common.mapping import load

    OUT.mkdir(parents=True, exist_ok=True)
    drop = set(CLOCK_GENES) | set(CORE)
    old = pd.read_parquet(DATA / "natexp/fetal_adrenal_open/pseudobulk.parquet")
    rows, frac = [], []
    for ds in ADRENAL_10X:
        # 'a244_spliced': the note 2 rule as fixed (spliced layer, A244 cells); 'gene_all': like for like with A242
        # (Gene counts = X, every non-blood barcode >= 500 UMI), added after the first comparison did not discriminate
        for kind, (new, sheet) in (("a244_spliced", load.pseudobulk(ds, "spliced", min_cells=0)),
                                   ("gene_all", load.pseudobulk(ds, "X", group=None, min_cells=0))):
            new = new.groupby(sheet.donor.to_numpy()).sum()
            genes = new.columns.intersection(old.columns).difference(list(drop))
            ref = old.loc[:, genes].dropna(axis=1, how="any")
            r = spearman_matrix(new[ref.columns], ref)
            n_genes = int(r.pop("genes").iloc[0])
            for d in new.index:
                o = r.loc[d].sort_values(ascending=False)
                rows.append({"dataset": ds, "kind": kind, "donor": d, "in_a242": d in ref.index,
                             "r_own": r.loc[d, d] if d in ref.index else np.nan, "best": o.index[0],
                             "r_best": o.iloc[0], "second": o.index[1], "r_second": o.iloc[1], "genes": n_genes})
    for ds in ADRENAL_10X + ORGAN_10X:
        c = load.cells(ds)
        if ds in ADRENAL_10X:
            c = c[c.a244_class.astype(str).ne("") & (c.total_counts >= 500)]
            cls = c.a244_class.astype(str)
        else:
            c = c[~c.blood.astype(bool) & (c.total_counts >= 500)]
            cls = pd.Series("organ", index=c.index)
        g = c.assign(cls=cls).groupby(["gsm", "donor", "cls"])[["spliced", "unspliced", "ambiguous"]].sum()
        g["cells"] = c.assign(cls=cls).groupby(["gsm", "donor", "cls"]).size()
        g["unspliced_frac"] = g.unspliced / (g.spliced + g.unspliced)
        frac.append(g.reset_index().assign(dataset=ds))
    a = pd.DataFrame(rows)
    a["mismatch"] = a.in_a242 & (a.best != a.donor)
    a["flag_low"] = a.in_a242 & (a.r_own < 0.9)
    f = pd.concat(frac, ignore_index=True)
    a.to_csv(OUT / "agreement_10x_a242.csv", index=False)
    f.to_csv(OUT / "unspliced_fraction_10x.csv", index=False)
    pd.set_option("display.width", 200)
    print(a.round(3).to_string(index=False))
    print(f.groupby(["dataset", "cls"]).unspliced_frac.describe()[["count", "min", "50%", "max"]].round(3).to_string())


if __name__ == "__main__":
    import sys
    main_cells() if "--cells" in sys.argv else main()
