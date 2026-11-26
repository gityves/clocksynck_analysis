"""A401: after build_census_tumour.py: dedup_candidates.csv, merged cell-type mapping, per-TCGA counts, PROVENANCE.md (sccc python).

Counts only (donors with tumour and >= 50 cells in k stromal compartments); no clock statistic.
"""
import hashlib
import json
import re
import time
import urllib.request
from pathlib import Path

import pandas as pd

OUT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/cellxgene_tumour_a401")
DATA = OUT.parent
HERE = Path(__file__).resolve().parent
FEAS = HERE.parent / "output" / "feasibility"
INV = HERE.parents[1] / "A388_malignant_tumour_host" / "output" / "feasibility" / "inventory.csv"
API = "https://api.cellxgene.cziscience.com/curation/v1/collections/"
ACC = re.compile(r"GSE\d{4,7}|SRP\d{5,7}|ERP\d{5,7}|PRJ[NE][AB]\d+|EGAS\d+|phs\d{6}|HRA\d+|syn\d{6,9}", re.I)
STROMA = ["fibroblast", "endothelial", "mural"]
NOTES = {  # known overlaps stated by the project (user brief 2026-11-26)
    "de5416ef": "Global Atlas: overlaps the other breast datasets of its collection; may include Wu 2021 / Pal 2021",
    "1e6a6ef9": "LuCA: may include on-disk lung cohorts",
    "5af90777": "already on disk as RCC_5af90777 (cellxgene_5af90777, Li 2022)",
    "b617ee1b": "pan-cancer atlas: partly on disk (A226 Qian 2020 donors, CRC and LUAD)",
}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def collection(cid):
    try:
        with urllib.request.urlopen(API + cid, timeout=60) as r:
            return json.load(r)
    except Exception as e:  # noqa: BLE001
        return {"error": repr(e)}


def main():
    t0 = time.time()
    stats = [json.loads(p.read_text()) for p in sorted(OUT.glob("*_stats.json"))]
    counts = pd.read_csv(FEAS / "census_counts.csv").drop_duplicates("dataset_id").set_index("dataset_id")
    inv = pd.read_csv(INV)
    prov = {}
    for p in list(DATA.glob("*/PROVENANCE.md")) + list(DATA.glob("*/*/PROVENANCE.md")) + list(DATA.glob("*/README*.md")):
        try:
            prov[str(p.relative_to(DATA))] = p.read_text(errors="ignore")
        except OSError:
            pass
    folders = [p.name for p in DATA.iterdir() if p.is_dir()]
    rows = []
    for s in stats:
        ds = s["dataset_id"]
        c = counts.loc[ds]
        col = collection(c.collection_id)
        pm = col.get("publisher_metadata") or {}
        authors = pm.get("authors") or []
        first = authors[0].get("family", "") if authors else ""
        year = pm.get("published_year", "")
        text = " ".join([col.get("description", "") or "", " ".join(l.get("link_url", "") for l in col.get("links", []))])
        accs = sorted(set(a.upper() for a in ACC.findall(text)))
        doi = col.get("doi", "") or ""
        hits = set()
        for f, t in prov.items():
            if ds in t or ds[:8] in t or c.collection_id in t or (doi and doi.lower() in t.lower()):
                hits.add(f.split("/PROVENANCE")[0])
            elif any(a in t.upper() for a in accs):
                hits.add(f.split("/PROVENANCE")[0] + " (accession)")
        hits |= {f for f in folders if ds[:8] in f or any(a.lower() in f.lower() for a in accs)}
        name_hits = {f for f in folders if first and len(first) > 3 and first.lower() in f.lower()}
        inv_hits = inv[inv.study.astype(str).str.upper().isin(accs)]
        rows.append(dict(
            dataset=s["short"], dataset_id=ds, dataset_title=c.dataset_title, collection_id=c.collection_id,
            collection_name=c.collection_name, doi=doi, journal=pm.get("journal", ""), year=year, first_author=first,
            is_preprint=pm.get("is_preprint", ""), links=";".join(f"{l.get('link_type')}:{l.get('link_url')}" for l in col.get("links", [])),
            geo_accessions_cited=";".join(accs), on_disk_provenance_hits=";".join(sorted(hits)),
            on_disk_author_name_hits=";".join(sorted(name_hits)), a388_inventory_hits=";".join(inv_hits.study.astype(str)),
            same_collection_donor_overlap=c.donor_overlap_datasets if isinstance(c.donor_overlap_datasets, str) else "",
            scout_on_disk=c.on_disk if isinstance(c.on_disk, str) else "", note=NOTES.get(ds[:8], ""),
            api_error=col.get("error", "")))
    pd.DataFrame(rows).to_csv(OUT / "dedup_candidates.csv", index=False)
    # merged mapping
    pd.concat([pd.read_csv(p) for p in sorted(OUT.glob("*_celltype_mapping.csv"))]).to_csv(OUT / "celltype_mapping.csv", index=False)
    # counts per TCGA code (tumour malignant/epithelial + k stromal at >= 50 cells), per donor x site
    import pyarrow.parquet as pq
    pb = pd.concat([pd.read_parquet(p, columns=["donor", "compartment", "n_cells", "tcga", "tissue_general", "disease"]
                                    + (["timepoint"] if "timepoint" in pq.read_schema(p).names else [])).assign(dataset=p.stem)
                    for p in sorted(OUT.glob("*.parquet")) if "_merged" not in p.stem])
    pb["timepoint"] = pb["timepoint"].fillna("") if "timepoint" in pb else ""
    t = pb[pb.compartment.str.startswith("tumour:") & (pb.n_cells >= 50)]
    t = t.assign(c=t.compartment.str.split(":").str[1])
    unit = ["dataset", "donor", "tcga", "tissue_general", "timepoint"]
    g = t.groupby(unit).apply(lambda x: pd.Series({
        "tumour": int(x.c.isin(["malignant", "epithelial"]).any()),
        "n_stroma": int(x.c.isin(STROMA).sum())}), include_groups=False).reset_index()
    g["q1"] = (g.tumour == 1) & (g.n_stroma >= 1)
    g["q2"] = (g.tumour == 1) & (g.n_stroma >= 2)
    tc = g.groupby("tcga").agg(donor_sites_q1=("q1", "sum"), donor_sites_q2=("q2", "sum"),
                               unique_donors_q1=("donor", lambda x: x[g.loc[x.index, "q1"]].nunique()),
                               unique_donors_q2=("donor", lambda x: x[g.loc[x.index, "q2"]].nunique()),
                               datasets=("dataset", "nunique")).reset_index()
    tc.to_csv(OUT / "counts_per_tcga.csv", index=False)
    g.groupby(["dataset", "tcga"]).agg(donor_sites_q1=("q1", "sum"), donor_sites_q2=("q2", "sum")).reset_index().to_csv(
        OUT / "counts_per_dataset_tcga.csv", index=False)
    # PROVENANCE
    files = sorted(p for p in OUT.rglob("*") if p.is_file() and p != OUT / "PROVENANCE.md")
    lines = ["# CELLxGENE Census tumour pseudobulks (A401)", "",
             f"- Built {time.strftime('%Y-%m-%d')} by `paper/exploratory/A401_cancer_timekeeping/scripts/build_census_tumour.py` "
             f"(+ `finalize_census_tumour.py`), Census 2025-11-08 (`cellxgene_census` 1.17.0, `census_env`), `is_primary_data == True`, "
             "`tissue_type == 'tissue'`. Raw counts summed to pseudobulks (donor x compartment x tissue site x assay x disease label) "
             "on the fly; no cell-level data kept. Pseudobulks with >= 20 cells kept. Columns: donor, compartment, stratum, n_cells, depth, "
             "median_umi (median obs raw_sum of the sample's cells: donor x state x tissue_general x assay), n_cells_sample, metadata, "
             "then gene symbols (duplicate symbols summed; all-zero genes dropped per dataset).",
             "- Selection: datasets/donors from `output/feasibility/census_counts.csv` (donor with >= 50 tumour cells and >= 50 cells in "
             "fibroblast / endothelial (non-lymphatic) / mural). For these donors all tumour-labelled and all disease 'normal' cells.",
             "- Compartment mapping from Census cell_type: `celltype_mapping.csv` (regexes in the build script). Tumour class: malignant "
             "labels when the dataset has any, else epithelial; adjacent_normal epithelial cells are `adjacent_normal:epithelial`.",
             "- Licence: CELLxGENE datasets are normally CC BY 4.0 (inferred; check the collection pages in `dedup_candidates.csv`).",
             "", "## Datasets", "", "| dataset | dataset id | collection | DOI | tumour label | donors | pseudobulks | genes |", "|---|---|---|---|---|---|---|---|"]
    dd = pd.DataFrame(rows).set_index("dataset")
    for s in stats:
        d = dd.loc[s["short"]]
        lines.append(f"| {s['short']} | {s['dataset_id']} | {d.collection_name} ({d.collection_id}) | {d.doi} | {s['tumour_label']} | "
                     f"{s['selected_donors']} | {s['pseudobulks']} | {s['genes']} |")
    lines += ["", "## Time-point split", "",
              "`integratedcancerce_7b20c613.parquet` is split by sample / time point (columns `timepoint` Pre / Post, `sample` = Sample|timepoint, outcome, treatment, "
              "study, ...); the earlier build keyed rows by donor only and mixed time points: kept as `integratedcancerce_7b20c613_merged.parquet` for reference. "
              "`integratedcancerce_7b20c613_obs_columns.csv` lists all source obs columns. BIOKEY labels: `biokey_labels/` (own PROVENANCE.md).",
              "", "## sha256", ""] + [f"- `{p.relative_to(OUT)}`  {sha256(p)}" for p in files]
    (OUT / "PROVENANCE.md").write_text("\n".join(lines) + "\n")
    print(tc.to_string())
    print("done", round(time.time() - t0), "s")


if __name__ == "__main__":
    main()
