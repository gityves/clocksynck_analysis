"""Prepare two public single-cell collections of the human fetal gonad: donor x cell-class pseudobulks and design tables.

Sources (CELLxGENE Discover h5ad, in `$SCCC_DATA/natexp/<folder>/raw/`; raw integer counts in `raw/X`, Ensembl ids in the
index of `raw/var`, symbols in `raw/var.feature_name`):
  garcia_alonso2022  collection 661a402a  "Human Male" 535e9336, "Human Female" 18e2a8c5 (Garcia-Alonso et al. 2022 Nature)
  gonad_tract2024    collection 222e800b  "male gonadal tissue" 43dc52cb, "female gonadal tissue" 7178ec0a (2024)

Design
- unit (donor) = `donor_id`; libraries of a donor (Garcia-Alonso `sample`, 2024 `batch`) are summed;
- cells flagged `is_doublet` (Garcia-Alonso, scrublet) are dropped; nothing else is filtered;
- class = `CLASS` rules on the CELLxGENE `cell_type` and the authors' finer label (written to `class_map_<key>.csv`);
- stage in weeks post fertilisation from `development_stage` (Carnegie stages converted with `CS_WPF`, approximate);
- 2024 collection: donors sample one tissue each (ovary / testis = gonad; mesonephros / epididymis dropped from the main
  tables, kept in `*_alltissue` class tables); Garcia-Alonso: all donors are gonad;
- outputs per dataset (`<key>` = garcia_female, garcia_male, gonad24_female, gonad24_male): `<key>_pseudobulk_class.h5ad`
  (donor x class, all genes, raw counts) and `<key>_pseudobulk_organ.h5ad` (donor, all gonad cells) in the data folder;
  design tables and clock-gene detection in `output/`.

Computed: sums of raw counts, cell counts, design summaries, median raw count of the clock genes per pseudobulk
(detection level), donor-identity checks between the collections (cell-level count tuples, whole-organ profile
similarity of same-sex donors, id overlap with Suo 2022).
NOT computed: any clock statistic (no gene-gene correlation, phase, coherence, gate, clocksync fit).

Usage: RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/scouts_2026-11-10/prep_fetal_gonad/scripts/prep_gonad.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

REPO = Path(__file__).resolve().parents[5]
sys.path[:0] = [str(REPO), str(REPO / "paper")]
import clocksync as cs  # noqa: E402
from clocksync.genes import CLOCK_GENES  # noqa: E402
from common.paths import DATA  # noqa: E402

NAT = DATA / "natexp"
OUT = Path(__file__).resolve().parents[1] / "output"
DATASETS = {  # key: (folder, h5ad, sex, authors' fine label, library column)
    "garcia_female": ("garcia_alonso2022", "002a804c", "female", "lineageSomatic", "sample"),
    "garcia_male": ("garcia_alonso2022", "62227d79", "male", "lineageSomatic", "sample"),
    "gonad24_female": ("gonad_tract2024", "485304df", "female", "cell_type_authors", "batch"),
    "gonad24_male": ("gonad_tract2024", "ac7034f8", "male", "cell_type_authors", "batch"),
}
GONAD = {"gonad", "ovary", "testis"}
CS_WPF = {13: 4.6, 14: 4.9, 15: 5.3, 16: 5.6, 17: 5.9, 18: 6.3, 19: 6.7, 20: 7.1, 21: 7.4, 22: 7.7, 23: 8.1}
CLASSES = ["endothelial", "mesenchymal_stromal", "mural", "supporting", "steroidogenic", "epithelial", "germ",
           "immune_blood", "neural", "other"]
CLASS = [  # (regex on "cell_type | author label", class); first match wins
    (r"^Leydig", "steroidogenic"), (r"^mesonephric podocyte", "other"), (r"Mixed population of germ", "other"),
    (r"^endothelial", "endothelial"), (r"^(mural|smooth muscle)", "mural"),
    (r"^(Sertoli|granulosa|pre-granulosa|supporting)", "supporting"),
    (r"^(mesenchymal|interstitial|stromal|cycling stromal)", "mesenchymal_stromal"),
    (r"^(epithelial|mesothelial|ovarian surface|secretory)", "epithelial"), (r"^germ", "germ"),
    (r"^(hematopoietic|erythrocyte)", "immune_blood"), (r"^neural", "neural"),
]
ALIAS = {"ARNTL": ["ARNTL", "BMAL1"], "CIART": ["CIART", "C1orf51"]}
AGE_BINS = [(0, 10, "<=9"), (10, 13, "10-12"), (13, 17, "13-16"), (17, 99, ">=17")]
MIN200, MIN50 = 200, 50
CHUNK = 20000


def to_class(label: str) -> str:
    return next((c for p, c in CLASS if re.search(p, label)), "other")


def stage_wpf(s: str) -> float:
    m = re.match(r"(\d+)(st|nd|rd|th) week", s)
    if m:
        return float(m.group(1))
    m = re.match(r"Carnegie stage (\d+)", s)
    return CS_WPF.get(int(m.group(1)), np.nan) if m else np.nan


def age_bin(w: float) -> str:
    return next((n for lo, hi, n in AGE_BINS if lo <= w < hi), "NA")


def cell_table(key: str) -> tuple[pd.DataFrame, Path]:
    folder, h5, sex, fine, lib = DATASETS[key]
    p = NAT / folder / "raw" / f"{h5}.h5ad"
    a = ad.read_h5ad(p, backed="r")
    o = a.obs
    cols = {"donor": "donor_id", "library": lib, "cell_type": "cell_type", "author": fine, "stage": "development_stage",
            "tissue": "tissue", "assay": "assay", "sex": "sex"}
    t = pd.DataFrame({k: o[v].astype(str).to_numpy() for k, v in cols.items()}, index=o.index.astype(str))
    t["doublet"] = o["is_doublet"].astype(str).eq("True").to_numpy() if "is_doublet" in o else False
    a.file.close()
    t["label"] = t.cell_type + " | " + t.author
    t["cls"] = [to_class(l) for l in t.label]
    t["wpf"] = t.stage.map(stage_wpf)
    t["gonad"] = t.tissue.isin(GONAD)
    return t, p


def read_raw(p: Path, t: pd.DataFrame, groupings: dict):
    """One pass over raw/X: group sums (dict name -> (labels, n_groups) -> dense float64) and per-cell (total, nnz)."""
    with h5py.File(p, "r") as f:
        var = read_elem(f["raw/var"])
        g = f["raw/X"]
        ip = g["indptr"][:]
        n, ng = len(ip) - 1, len(var)
        assert n == len(t), (n, len(t))
        samp = g["data"][:200000]
        integer = bool(np.all(samp == np.round(samp)))
        accs = {k: np.zeros((len(u), ng)) for k, (lab, u) in groupings.items()}
        tot, nnz = np.zeros(n), np.zeros(n, int)
        for s in range(0, n, CHUNK):
            e = min(s + CHUNK, n)
            a, b = ip[s], ip[e]
            X = sp.csr_matrix((g["data"][a:b], g["indices"][a:b], ip[s:e + 1] - a), shape=(e - s, ng))
            tot[s:e], nnz[s:e] = np.asarray(X.sum(1)).ravel(), np.diff(X.indptr)
            for k, (lab, u) in groupings.items():
                l = lab[s:e]
                ok = l >= 0
                G = sp.csr_matrix((np.ones(ok.sum()), (l[ok], np.where(ok)[0])), shape=(len(u), e - s))
                accs[k] += (G @ X).toarray()
    return var, accs, tot, nnz, integer


def codes(keys: pd.Series, mask: np.ndarray):
    """Integer group code per cell (-1 where not mask) and the sorted unique keys."""
    u = sorted(keys[mask].unique())
    m = {k: i for i, k in enumerate(u)}
    return np.where(mask, keys.map(m).fillna(-1).to_numpy(), -1).astype(int), u


def write_pb(path, X, rows: pd.DataFrame, var: pd.DataFrame):
    v = pd.DataFrame({"ensembl_id": var.index.astype(str), "symbol": var["feature_name"].astype(str).to_numpy()},
                     index=var.index.astype(str))
    a = ad.AnnData(X.astype(np.float32), obs=rows.set_index(rows.index.astype(str)), var=v)
    a.write_h5ad(path, compression="gzip")


def donor_info(t: pd.DataFrame) -> pd.DataFrame:
    g = t.groupby("donor")
    j = lambda c: g[c].agg(lambda x: ";".join(sorted(set(x))))  # noqa: E731
    return pd.DataFrame({"sex": j("sex"), "stage": j("stage"), "wpf": g.wpf.first(), "tissue": j("tissue"),
                         "assay": g.assay.agg(lambda x: x.value_counts().index[0]), "assays_all": j("assay"),
                         "n_libraries": g.library.nunique(), "libraries": j("library"), "n_cells": g.size()})


def build(key: str):
    t, p = cell_table(key)
    use = (~t.doublet).to_numpy()
    main = use & t.gonad.to_numpy()
    cl_all, ua = codes(t.donor + "||" + t.cls + "||" + t.tissue, use)
    cl_main = np.where(main, cl_all, -1)
    og, uo = codes(t.donor, main)
    var, accs, tot, nnz, integer = read_raw(p, t, {"cls": (cl_all, ua), "org": (og, uo)})
    t["total"], t["nnz"] = tot, nnz
    info = donor_info(t[use])
    rows = pd.DataFrame([k.split("||") for k in ua], columns=["donor", "compartment", "tissue"])
    cells = pd.Series(np.bincount(cl_all[cl_all >= 0], minlength=len(ua)))
    rows["n_cells"] = cells.to_numpy()
    rows["depth"] = accs["cls"].sum(1)
    for c in ("sex", "stage", "wpf", "assay"):
        rows[c] = rows.donor.map(info[c])
    rows["dataset"], rows["gonad"] = key, rows.tissue.isin(GONAD)
    rows.index = rows.donor + "|" + rows.compartment
    sel = rows.gonad.to_numpy()
    write_pb(NAT / DATASETS[key][0] / f"{key}_pseudobulk_class.h5ad", accs["cls"][sel], rows[sel], var)
    if not sel.all():
        write_pb(NAT / DATASETS[key][0] / f"{key}_pseudobulk_class_alltissue.h5ad", accs["cls"], rows, var)
    orows = pd.DataFrame({"donor": uo, "compartment": "organ"})
    orows["n_cells"] = np.bincount(og[og >= 0], minlength=len(uo))
    orows["depth"] = accs["org"].sum(1)
    for c in ("sex", "stage", "wpf", "assay", "tissue"):
        orows[c] = orows.donor.map(info[c])
    orows["dataset"] = key
    orows.index = orows.donor
    write_pb(NAT / DATASETS[key][0] / f"{key}_pseudobulk_organ.h5ad", accs["org"], orows, var)
    # class map for review
    m = t[use].groupby(["cell_type", "author", "cls"]).size().rename("n_cells").reset_index()
    m.columns = ["cell_type", "author_label", "class", "n_cells"]
    m.to_csv(OUT / f"class_map_{key}.csv", index=False)
    info.to_csv(OUT / f"donors_{key}.csv")
    # report
    print(f"[{key}] cells {len(t)}, doublets dropped {int(t.doublet.sum())}, non-gonad cells {int((~t.gonad & ~t.doublet).sum())}, "
          f"raw/X integer {integer}, genes {len(var)}, donors {t.donor.nunique()}, gonad donors {t[t.gonad].donor.nunique()}",
          flush=True)
    print(f"  libraries per donor (max) {info.n_libraries.max()}, donors with >1 library {int((info.n_libraries > 1).sum())}, "
          f"donors with >1 assay {int((info.assays_all.str.contains(';')).sum())}", flush=True)
    pbc = rows[sel].reset_index(drop=True)
    cellmat = (var, accs["cls"][sel], pbc, accs["org"], orows)
    return t, info, cellmat


def design(rows: pd.DataFrame, info: pd.DataFrame, name: str):
    """Design tables for one dataset or collection: `rows` = donor x class pseudobulk rows (gonad)."""
    r = rows.copy()
    r["bin"] = r.wpf.map(age_bin)
    r["core"] = ~r.compartment.isin(["immune_blood", "germ", "other"]) & (r.n_cells >= MIN200)
    strata = [("all", "all", r)] + [("sex", s, g) for s, g in r.groupby("sex")] + \
             [("age", b, g) for b, g in r.groupby("bin")] + \
             [("assay", a, g) for a, g in r.groupby("assay")]
    out, multi = [], []
    for typ, s, g in strata:
        for c in sorted(set(r.compartment)):
            x = g[g.compartment == c]
            out.append(dict(dataset=name, stratum_type=typ, stratum=s, cls=c, donors_ge200=int((x.n_cells >= MIN200).sum()),
                            donors_ge50=int((x.n_cells >= MIN50).sum()), donors_any=len(x),
                            median_n_cells=float(x.n_cells.median()) if len(x) else np.nan,
                            median_depth=float(x.depth.median()) if len(x) else np.nan))
        k = g[g.core].groupby("donor").size()
        multi.append(dict(dataset=name, stratum_type=typ, stratum=s, donors=g.donor.nunique(),
                          ge2_classes=int((k >= 2).sum()), ge3_classes=int((k >= 3).sum()), ge4_classes=int((k >= 4).sum())))
    return pd.DataFrame(out), pd.DataFrame(multi)


def clock_detection(key, var, X, rows):
    sym = var["feature_name"].astype(str).to_numpy()
    genes = list(cs.gtex().genes)
    assert tuple(genes) == tuple(CLOCK_GENES)
    out = []
    for g in genes:
        idx = [i for a in ALIAS.get(g, [g]) for i in np.where(sym == a)[0]]
        for c in sorted(set(rows.compartment)):
            m = (rows.compartment == c).to_numpy() & (rows.n_cells.to_numpy() >= MIN50)
            v = X[m][:, idx].sum(1) if idx and m.any() else np.array([])
            out.append(dict(dataset=key, cls=c, gene=g, symbol_found=";".join(sym[idx]) if idx else "ABSENT", n_pseudobulks=int(m.sum()),
                            median_count=float(np.median(v)) if len(v) else np.nan,
                            frac_ge1=float((v >= 1).mean()) if len(v) else np.nan))
    pd.DataFrame(out).to_csv(OUT / f"clock_detection_{key}.csv", index=False)


def overlap(T: dict, organ: dict):
    """Donor identity between the collections and with Suo 2022 (ids, library names, per-cell count tuples, profiles)."""
    res = []
    for sx in ("female", "male"):
        g, n = T[f"garcia_{sx}"], T[f"gonad24_{sx}"]
        gs = g[~g.doublet].assign(k=lambda d: list(zip(d.total.round().astype(int), d.nnz)))
        ns = n[~n.doublet].assign(k=lambda d: list(zip(d.total.round().astype(int), d.nnz)))
        gset = {}
        for lib, k in zip(gs.library, gs.k):
            gset.setdefault(k, set()).add(lib)
        hits = pd.Series([len(gset.get(k, ())) > 0 for k in ns.k], index=ns.index).groupby(ns.donor.to_numpy()).mean()
        hits.rename("frac_identical_tuple").to_csv(OUT / f"donor_overlap_cell_tuples_{sx}.csv")
        print(f"[overlap {sx}] fraction of 2024 cells with an identical (total counts, n genes) tuple in a Garcia cell, per donor:",
              hits.round(3).to_dict(), flush=True)
        # whole-organ log-CPM profiles on the shared Ensembl genes
        va, Xa, ra = organ[f"garcia_{sx}"]
        vb, Xb, rb = organ[f"gonad24_{sx}"]
        sh = va.index.intersection(vb.index)
        la = np.log1p(1e6 * Xa[:, va.index.get_indexer(sh)] / Xa.sum(1, keepdims=True))
        lb = np.log1p(1e6 * Xb[:, vb.index.get_indexer(sh)] / Xb.sum(1, keepdims=True))
        la, lb = la - la.mean(0), lb - lb.mean(0)
        C = np.corrcoef(np.vstack([lb, la]))[: len(lb), len(lb):]
        best = pd.DataFrame({"donor24": rb.donor.to_numpy(), "best_garcia": ra.donor.to_numpy()[C.argmax(1)],
                             "r_best": C.max(1), "r_median": np.median(C, 1),
                             "wpf24": rb.wpf.to_numpy(), "wpf_garcia": ra.wpf.to_numpy()[C.argmax(1)]})
        res.append(best.assign(sex=sx))
        print(best.round(3).to_string(), flush=True)
    pd.concat(res).to_csv(OUT / "donor_overlap_profile_similarity.csv", index=False)
    # id overlap with Suo 2022 and between collections
    gids = set(pd.concat([T["garcia_female"], T["garcia_male"]]).donor)
    nids = set(pd.concat([T["gonad24_female"], T["gonad24_male"]]).donor)
    suo = ad.read_h5ad(NAT / "suo2022" / "pseudobulk" / "pseudobulk_organ_all" / "pseudobulk.h5ad", backed="r").obs
    sids = set(suo["donor_id"].astype(str))
    print("Garcia donors", len(gids), "2024 donors", len(nids), "Suo donors", len(sids), flush=True)
    print("Garcia & 2024 ids:", sorted(gids & nids), "| Garcia & Suo ids:", sorted(gids & sids), flush=True)
    print("Suo ids:", sorted(sids), flush=True)
    glib = set(pd.concat([T["garcia_female"], T["garcia_male"]]).library)
    print("Garcia libraries (sample ids) shared with 2024 batch labels:", sorted(glib & set(pd.concat(
        [T["gonad24_female"], T["gonad24_male"]]).library)), flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    T, INFO, CM = {}, {}, {}
    for key in DATASETS:
        T[key], INFO[key], CM[key] = build(key)
        var, Xc, rows, Xo, orows = CM[key]
        clock_detection(key, var, np.vstack([Xc, Xo]), pd.concat([rows[["compartment", "n_cells"]],
                                                              orows[["compartment", "n_cells"]]], ignore_index=True))
    d, m = [], []
    for coll, keys in {"garcia_alonso2022": ["garcia_female", "garcia_male"], "gonad_tract2024": ["gonad24_female", "gonad24_male"]}.items():
        rows = pd.concat([CM[k][2] for k in keys], ignore_index=True)
        info = pd.concat([INFO[k] for k in keys])
        a, b = design(rows, info, coll)
        d.append(a), m.append(b)
        da = info[info.tissue.str.contains("|".join(GONAD))].groupby("assay").size().rename("donors").reset_index().assign(dataset=coll)
        da.to_csv(OUT / f"design_donors_per_assay_{coll}.csv", index=False)
        print(f"[design {coll}] donors per assay\n{da.to_string(index=False)}", flush=True)
    D, M = pd.concat(d), pd.concat(m)
    D.to_csv(OUT / "design_donors_per_class.csv", index=False)
    M.to_csv(OUT / "design_multiclass_donors.csv", index=False)
    pd.set_option("display.width", 250, "display.max_rows", 500)
    print(D[D.stratum_type.isin(["all", "sex", "age"])].to_string(index=False), flush=True)
    print(M.to_string(index=False), flush=True)
    overlap(T, {k: (CM[k][0], CM[k][3], CM[k][4]) for k in CM})
    print("done; no clock statistic computed", flush=True)


if __name__ == "__main__":
    main()
