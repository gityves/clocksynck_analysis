"""A214 scout: inventory (counts and UMI only, no clock statistic) of the single-cell tumour cohorts on disk.

Question: per study, how many patients can compare malignant cells with normal epithelium of the same patient?
Eligible patient = >= 100 malignant cells (`tumour:malignant`, copy-number call) AND >= 100 normal epithelial cells, with
two definitions of normal epithelium: (a) `adjacent_normal:epithelial`, (b) `tumour:epithelial_nonmalignant` (non-malignant
epithelium inside the tumour sample, copy-number call "normal").

Sources (nothing rebuilt)
- 43 `cancer_sc` cohorts: per-cell tables `cells_cancer.load(study, "dispersion")` (legacy `outputs/celldisp` files in
  practice; 38 exist), copy-number calls `$SCCC_DATA/cnv/<study>.csv.gz` (29 cohorts), split by `cell_groups` of
  `figure_2/scripts/backbone/compute_cancer_split.py`. Without barcodes no split is possible: (b) and `tumour:malignant`
  are then NA; (a) is then NA as well because malignant cells cannot be told from normal ones.
- A149 new cohorts (17): only pseudobulks (`A149_paired_tumour_normal/output/pb/*.parquet`, donor x compartment, `n_cells`,
  14 clock counts, `depth` = summed genome UMI). Patient numbers from the pseudobulk table; clock UMI per cell is the group
  mean (no medians), and compartments of < 50 cells are absent by construction.
- Wu GSE176078 (author labels, `metadata.csv`): cell counts and `nCount_RNA` only (no clock counts read).
- CoPhaser fits: judgement of `C4_single_cell/fits/FITS.md` (hand-copied, `FIT_VALID`), cross-checked with `fits/output/qc.csv`.
Cycling annotation possible = valid fit or median genome UMI per malignant cell > 1,500.

    python inventory_disk.py
Outputs: ../inventory/inventory_patients.csv, inventory_studies.csv
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
REPO = HERE.parents[4]
PAPER = REPO / "paper"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(PAPER))
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import CLOCK_GENES  # noqa: E402
from common import data  # noqa: E402
from common.paths import ARCHIVE, DATA  # noqa: E402

limit_threads(4)
OUT = HERE.parents[1] / "inventory"
GENES = list(CLOCK_GENES)
MIN_CELLS = 100
MAL, NMAL, ADJ = "tumour:malignant", "tumour:epithelial_nonmalignant", "adjacent_normal:epithelial"
spec = importlib.util.spec_from_file_location("ccs", PAPER / "results/figure_2/scripts/backbone/compute_cancer_split.py")
ccs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ccs)  # cell_groups (read only)
cc, scb = data._builder("cells_cancer"), data._builder("cancer_sc")
LEGACY = ARCHIVE / "outputs" / "celldisp"
A149 = PAPER / "exploratory/A149_paired_tumour_normal"
units = pd.read_csv(A149 / "output/units.csv")
CANCER = dict(zip(units.study, units.cancer))
# FITS.md validity table + pilot paragraph; values: y, y_low (one seed / low confidence), suspect, weak, n
FIT_VALID = {"GSE161529": "y", "ERP149257": "y_low", "GSE166555": "suspect", "GSE181919": "y", "GSE200997": "weak",
             "GSE294300": "y", "GSE315534": "suspect", "GSE335070": "weak", "SRP416180": "y", "GSE159115": "n",
             "GSE212966": "n", "GSE270680": "n", "GSE178341": "y", "GSE132465": "y", "GSE183904": "y",
             "SRP327758": "y", "GSE176078": "y", "GSE131907": "y"}


def stats(g: pd.DataFrame) -> dict:
    u = g[GENES].to_numpy(float).sum(1)
    return dict(n_cells=len(g), clock_umi_median=float(np.median(u)), clock_umi_mean=float(u.mean()),
                depth_median=float(g.depth.median()))


def cohort(study: str):
    """Per-(patient, compartment) rows, pooled per-cell arrays of eligible patients, notes."""
    f = LEGACY / f"{study}_cells.csv.gz"
    cells = cc.load(study, "dispersion")
    has_bc = "barcode" in cells
    has_calls = (scb.S.CNV / f"{study}.csv.gz").exists()
    calls = None
    note = []
    if has_bc and has_calls:
        calls = scb.S.read_calls(study).rename("cnv_call").reset_index()
    elif not has_bc and has_calls:
        note.append("calls exist but legacy table has no barcodes: no malignant/normal split")
    elif not has_calls:
        note.append("no copy-number calls: tumour epithelium not split")
    rows, grp = [], {}
    for comp, r, g in ccs.cell_groups(cells, calls):
        rows.append(dict(study=study, patient=r, compartment=comp, **stats(g)))
        grp[(r, comp)] = g
    return pd.DataFrame(rows), grp, note, f.stat().st_size / 1e6, cells


def study_row(study, P, grp, note, mb, cells, can_split):
    d = P.pivot(index="patient", columns="compartment", values="n_cells") if len(P) else pd.DataFrame()
    n = lambda c: d[c].fillna(0) if c in d else pd.Series(0, index=d.index)
    out = dict(study=study, cancer=CANCER.get(study, ""), assay="droplet (UMI)", n_patients=cells.resection.nunique(),
               table_mb=round(mb, 2), source="per-cell table")
    chem = cells.stratum.replace("", np.nan).dropna().astype(str)
    out["chemistry"] = ";".join(sorted(chem.unique())[:4]) if len(chem) else ""
    for tag, norm in (("a", ADJ), ("b", NMAL)):
        if not can_split:
            out.update({f"n_eligible_{tag}": np.nan, f"mal_cells_{tag}": np.nan, f"norm_cells_{tag}": np.nan,
                        f"mal_clock_umi_median_{tag}": np.nan, f"norm_clock_umi_median_{tag}": np.nan})
            continue
        ok = list(d.index[(n(MAL) >= MIN_CELLS) & (n(norm) >= MIN_CELLS)])
        out[f"n_eligible_{tag}"] = len(ok)
        out[f"mal_cells_{tag}"] = int(sum(n(MAL)[ok])) if ok else 0
        out[f"norm_cells_{tag}"] = int(sum(n(norm)[ok])) if ok else 0
        for k, comp in (("mal", MAL), ("norm", norm)):
            u = np.concatenate([grp[(r, comp)][GENES].to_numpy(float).sum(1) for r in ok]) if ok else np.array([])
            out[f"{k}_clock_umi_median_{tag}"] = float(np.median(u)) if len(u) else np.nan
        out[f"_ok_{tag}"] = ok
    mal = [grp[k].depth for k in grp if k[1] == MAL]
    out["depth_median_malignant"] = float(np.median(np.concatenate(mal))) if mal else np.nan
    out["notes"] = "; ".join(note)
    return out


def a149_pb():
    """A149 cohorts without per-cell tables: patient rows from the pseudobulks."""
    rows, srows = [], []
    for f in sorted((A149 / "output/pb").glob("*.parquet")):
        s = f.stem
        pb = pd.read_parquet(f)
        u = pb[GENES].sum(axis=1) / pb.n_cells
        pb = pb.assign(cu=u, dm=pb.depth / pb.n_cells)
        for r in pb.itertuples():
            rows.append(dict(study=s, patient=r.donor, compartment=r.compartment, n_cells=int(r.n_cells), clock_umi_median=np.nan,
                             clock_umi_mean=float(r.cu), depth_median=np.nan, depth_mean=float(r.dm), source="pseudobulk"))
        d = pb.pivot_table(index="donor", columns="compartment", values="n_cells", aggfunc="sum")
        g = lambda c: d[c].fillna(0) if c in d else pd.Series(0, index=d.index)
        malc = MAL if MAL in d else ("tumour:epithelial" if "tumour:epithelial" in d else None)
        row = dict(study=s, cancer=CANCER.get(s, ""), assay="droplet/other (see A149 candidates.md)",
                   n_patients=pb.donor.nunique(), table_mb=np.nan, source="pseudobulk only (no per-cell table)", chemistry="")
        notes = [] if malc == MAL else ["no malignant call: 'tumour:epithelial' used as the tumour group (not copy-number split)"]
        for tag, norm in (("a", ADJ), ("b", NMAL)):
            if tag == "b":  # the pseudobulks carry no non-malignant tumour epithelium compartment
                row.update({f"n_eligible_{tag}": np.nan, f"mal_cells_{tag}": np.nan, f"norm_cells_{tag}": np.nan})
                continue
            if malc is None or norm not in d:
                row.update({f"n_eligible_{tag}": 0, f"mal_cells_{tag}": 0, f"norm_cells_{tag}": 0})
                continue
            ok = list(d.index[(g(malc) >= MIN_CELLS) & (g(norm) >= MIN_CELLS)])
            row.update({f"n_eligible_{tag}": len(ok), f"mal_cells_{tag}": int(g(malc)[ok].sum()), f"norm_cells_{tag}": int(g(norm)[ok].sum())})
            for k, c in (("mal", malc), ("norm", norm)):
                x = pb[(pb.compartment == c) & pb.donor.isin(ok)].cu
                row[f"{k}_clock_umi_median_{tag}"] = float(x.median()) if len(x) else np.nan  # median over patients of group means
        row["depth_median_malignant"] = float(pb[pb.compartment == malc].dm.median()) if malc else np.nan
        row["notes"] = "; ".join(notes + ["(b) not available in pseudobulks", "clock UMI = median over patients of group means"])
        srows.append(row)
    return pd.DataFrame(rows), srows


def wu_row():
    import glob
    b = glob.glob(str(DATA / "sc_breast_er/GSE176078/extract/*"))[0]
    m = pd.read_csv(b + "/metadata.csv", index_col=0)
    m = m[m.celltype_major.isin(["Cancer Epithelial", "Normal Epithelial"])]
    t = m.groupby(["orig.ident", "celltype_major"]).size().unstack(fill_value=0)
    ok = t.index[(t["Cancer Epithelial"] >= MIN_CELLS) & (t["Normal Epithelial"] >= MIN_CELLS)]
    return dict(study="GSE176078", cancer="breast", assay="droplet (UMI)", n_patients=int(pd.read_csv(b + "/metadata.csv", usecols=["orig.ident"])["orig.ident"].nunique()),
                source="author labels (metadata.csv)", n_eligible_a=0, n_eligible_b=len(ok), mal_cells_b=int(t.loc[ok, "Cancer Epithelial"].sum()),
                norm_cells_b=int(t.loc[ok, "Normal Epithelial"].sum()), depth_median_malignant=float(
                    pd.read_csv(b + "/metadata.csv").query("celltype_major == 'Cancer Epithelial'")["nCount_RNA"].median()),
                notes="Wu: tumour samples only (no adjacent normal); (b) = author 'Normal Epithelial' vs 'Cancer Epithelial' (not copy-number); clock UMI not read")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows, srows = [], []
    for study in list(units.study[~units.new]):
        if not (LEGACY / f"{study}_cells.csv.gz").exists():
            srows.append(dict(study=study, cancer=CANCER.get(study, ""), source="no per-cell table on disk", notes="in cancer_sc (pseudobulk) only"))
            continue
        P, grp, note, mb, cells = cohort(study)
        can = (scb.S.CNV / f"{study}.csv.gz").exists() and "barcode" in cells
        srows.append(study_row(study, P, grp, note, mb, cells, can))
        P["source"] = "per-cell table"
        all_rows.append(P)
        print(study, flush=True)
    pp = pd.concat(all_rows, ignore_index=True)
    pb_rows, a149 = a149_pb()
    srows += a149
    srows.append(wu_row())
    pp = pd.concat([pp, pb_rows], ignore_index=True)
    S = pd.DataFrame(srows)
    S["cophaser_fit_valid"] = S.study.map(FIT_VALID).fillna("not fitted")
    S["cycling_annotation_possible"] = (S.cophaser_fit_valid.isin(["y", "y_low"]) | (S.depth_median_malignant > 1500)).map({True: "y", False: "n"})
    S = S.drop(columns=[c for c in S if c.startswith("_ok")])
    S["n_eligible_either"] = np.nan
    # either definition: patients eligible under (a) or (b) (per-cell tables only; same patient counted once)
    for s in S.study:
        sub = pp[(pp.study == s) & (pp.source == "per-cell table")] if (S[S.study == s].source == "per-cell table").any() else pp[pp.study == s]
        d = sub.pivot_table(index="patient", columns="compartment", values="n_cells", aggfunc="sum")
        if MAL not in d:
            continue
        m = d[MAL].fillna(0) >= MIN_CELLS
        e = m & ((d[ADJ].fillna(0) >= MIN_CELLS if ADJ in d else False) | (d[NMAL].fillna(0) >= MIN_CELLS if NMAL in d else False))
        S.loc[S.study == s, "n_eligible_either"] = int(e.sum())
    pp.to_csv(OUT / "inventory_patients.csv", index=False)
    S.to_csv(OUT / "inventory_studies.csv", index=False)
    print(S.drop(columns="notes").round(1).to_string())


if __name__ == "__main__":
    main()
