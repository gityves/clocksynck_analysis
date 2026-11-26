"""A401: genome-wide donor x compartment pseudobulks of CELLxGENE Census tumour datasets (counts only, no statistic).

Run in `census_env` (cellxgene_census, tiledbsoma). Census 2025-11-08, primary data, tissue samples. Datasets and donors are
selected from the scout (`output/feasibility/census_counts.csv`): datasets with >= 1 donor with >= 50 tumour cells and >= 50
cells in >= 1 of fibroblast / endothelial / mural. For the selected donors every tumour-labelled cell and every
disease == 'normal' cell is streamed from the Census raw counts and summed into pseudobulks on the fly (no cell-level matrix
is kept). Output per dataset in $SCCC_DATA/cellxgene_tumour_a401/: <short>.parquet, <short>_donors.csv,
<short>_celltype_mapping.csv, <short>_stats.json; `finalize.py` writes PROVENANCE.md and dedup_candidates.csv.

    build_census_tumour.py [--datasets id8 id8 ...] [--workers 4] [--threads 2]
"""
import argparse
import json
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import census_scout as S  # noqa: E402  (regexes and tcga(); main() is guarded)

VERSION = S.VERSION
OUTDIR = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/cellxgene_tumour_a401")
COUNTS = HERE.parent / "output" / "feasibility" / "census_counts.csv"
MIN_CELLS = 20
MYELOID = re.compile(r"macrophage|monocyte|dendritic|myeloid|microglia|kupffer", re.I)
TNK = re.compile(r"\bT cell|T-regulatory|regulatory T|T helper|gamma-delta|T-cell|T lymphocyte|natural killer|\bNK\b|\bNKT\b|\bMAIT|mucosal invariant", re.I)
BPLASMA = re.compile(r"\bB cell|B-cell|plasma cell|plasmablast|plasmacyte", re.I)
OBS_COLS = ["soma_joinid", "donor_id", "cell_type", "assay", "tissue", "tissue_general", "disease", "sex",
            "development_stage", "suspension_type", "raw_sum", "tissue_type", "dataset_id", "observation_joinid"]
EXTRA_KEYS = re.compile(r"treat|therap|respon|chemo|drug|time|day|week|stage|grade|subtype|histolog|status|surviv|recur|"
                        r"prior|batch|sample|site|location|age|sex|gender|mutation|braak|pathol|biops|condition|diagnos|"
                        r"cohort|patient|tnm|molecular|hpv|msi|stage|risk", re.I)


def short_name(ds_id, title):
    slug = re.sub(r"[^a-z0-9]+", "", re.sub(r"\(.*?\)", "", title).lower())[:18]
    return f"{slug}_{ds_id[:8]}"


def compartment(ct, state, has_mal):
    """Compartment of a Census cell_type (None = not used). state: tumour | adjacent_normal."""
    if S.MURAL.search(ct):
        return "mural"
    if S.ENDO_LYMPH.search(ct):
        return None
    if S.ENDO.search(ct):
        return "endothelial"
    if S.FIB.search(ct) and not S.MALIG.search(ct):
        return "fibroblast"
    if MYELOID.search(ct):
        return "myeloid"
    if BPLASMA.search(ct):
        return "B_plasma"
    if TNK.search(ct):
        return "T_NK"
    mal = bool(S.MALIG.search(ct))
    if state == "tumour":
        if has_mal:
            return "malignant" if mal else None
        return "epithelial" if (S.EPI.search(ct) and not mal) else None
    return "epithelial" if (S.EPI.search(ct) and not mal) else None


# datasets whose pseudobulks are split by time point (source-obs columns merged on observation_joinid)
SPLIT_DATASETS = {"7b20c613-9add-43d1-87e9-defd3d9b9f8c": ["Sample", "pre_post", "outcome", "Treatment.or.Mode.of.Action", "Study_name",
                                                              "Combined_outcome", "Primary_or_met", "Cancer_type_update"]}


def read_source_extras(dataset_id, obs_joinids_needed, threads=2, only=None):
    """Donor-relevant columns of the source h5ad obs over anonymous S3 (treatment / response / time / stage ... fields)."""
    import cellxgene_census as cc
    import fsspec
    import h5py
    uri = cc.get_source_h5ad_uri(dataset_id, census_version=VERSION)["uri"]
    fs = fsspec.filesystem("s3", anon=True, default_block_size=2**22)
    out = {}
    with fs.open(uri, "rb") as f, h5py.File(f, "r") as h:
        obs = h["obs"]
        for c in obs.keys():
            if only is not None and c not in only and c != "observation_joinid":
                continue
            if only is None and (c.startswith("_") or re.search(r"ontology|dataset_id|is_primary|tissue_type", c) or not (c == "observation_joinid" or EXTRA_KEYS.search(c)) or c in ("donor_id", "sex", "tissue", "disease", "assay",
                                                                      "development_stage", "suspension_type")):
                continue
            try:
                g = obs[c]
                if isinstance(g, h5py.Group) and "categories" in g:
                    cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][()]], dtype=object)
                    codes = g["codes"][()]
                    out[c] = np.where(codes < 0, "NA", cats[np.clip(codes, 0, None)]).astype(str)
                elif isinstance(g, h5py.Dataset):
                    v = g[()]
                    out[c] = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v], dtype=object) \
                        if v.dtype.kind in "SO" else v.astype(str)
            except Exception as e:  # noqa: BLE001
                print("   skip source column", c, e, flush=True)
        oj = out.pop('observation_joinid', None)
    return pd.DataFrame(out), oj


def build_one(args):
    try:
        return _build_one(args)
    except Exception as e:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        return args[0][:8], "FAILED " + repr(e)[:200]


def _build_one(args):
    ds_id, threads = args
    import cellxgene_census as cc
    import tiledbsoma as soma
    t0 = time.time()
    meta = pd.read_csv(COUNTS)
    meta = meta[meta.dataset_id == ds_id].iloc[0]
    short = short_name(ds_id, meta.dataset_title)
    out = OUTDIR / f"{short}.parquet"
    if out.exists():
        return short, "exists"
    log = lambda *a: print(f"[{short} {time.time() - t0:5.0f}s]", *a, flush=True)  # noqa: E731
    ctx = soma.SOMATileDBContext(tiledb_config={"soma.init_buffer_bytes": str(128 * 2**20),
                                                "sm.compute_concurrency_level": str(threads),
                                                "sm.io_concurrency_level": str(threads)})
    census = cc.open_soma(census_version=VERSION, context=ctx)
    try:
        exp = census["census_data"]["homo_sapiens"]
        base_f = f"dataset_id == '{ds_id}' and is_primary_data == True and tissue_type == 'tissue'"
        obs = exp.obs.read(value_filter=base_f, column_names=OBS_COLS).concat().to_pandas()
        for c in obs.columns:
            if str(obs[c].dtype) == "category" or obs[c].dtype == object:
                obs[c] = obs[c].astype(str)
        xcols = []
        if ds_id in SPLIT_DATASETS:
            src, oj = read_source_extras(ds_id, None, only=SPLIT_DATASETS[ds_id])
            src["observation_joinid"] = oj
            obs = obs.merge(src, on="observation_joinid", how="left", validate="one_to_one")
            obs["timepoint"] = obs.pre_post.astype(str)
            obs["sample"] = obs.Sample.astype(str) + "|" + obs.timepoint
            xcols = ["timepoint", "sample", "outcome", "Treatment.or.Mode.of.Action", "Study_name", "Combined_outcome",
                     "Primary_or_met", "Cancer_type_update"]
            log(f"split by time point: {obs.timepoint.value_counts().to_dict()}")
        tp = (lambda x: "|" + x) if xcols else (lambda x: "")
        obs["is_tumour"] = obs.disease.map(lambda d: bool(S.TUMOUR.search(d)) and not S.EXCLUDE.search(d))
        obs["state"] = np.where(obs.is_tumour, "tumour", np.where(obs.disease == "normal", "adjacent_normal", "other"))
        has_mal = bool(obs[obs.is_tumour].cell_type.str.contains(S.MALIG).any())
        obs["comp"] = [compartment(c, s, has_mal) if s != "other" else None for c, s in zip(obs.cell_type, obs.state)]
        # qualifying donors (as the scout): per disease x tissue_general x donor
        t = obs[obs.is_tumour]
        rows = []
        for (dis, tg, don), g in t.groupby(["disease", "tissue_general", "donor_id"]):
            n = g.comp.value_counts()
            rows.append(dict(donor_id=don, disease=dis, tissue_general=tg, tum=n.get("malignant", 0) + n.get("epithelial", 0),
                             stroma=max(n.get("fibroblast", 0), n.get("endothelial", 0), n.get("mural", 0))))
        q = pd.DataFrame(rows)
        sel = sorted(q[(q.tum >= 50) & (q.stroma >= 50)].donor_id.unique())
        log(f"has_malignant_label={has_mal}; {obs.donor_id.nunique()} donors, {len(sel)} qualifying (q1); {len(obs)} cells")
        if not sel:
            return short, "no qualifying donors"
        o = obs[obs.donor_id.isin(sel) & (obs.state != "other")].copy()
        mapping = o.groupby(["state", "cell_type"]).size().rename("n_cells").reset_index()
        mapping["compartment"] = [compartment(c, s, has_mal) or "" for s, c in zip(mapping.state, mapping.cell_type)]
        mapping.insert(0, "dataset", short)
        mapping["tumour_label"] = ("malignant" if has_mal else "epithelial")
        o["sample_key"] = o.donor_id + "|" + o.state + "|" + o.tissue_general + "|" + o.assay + (tp(o.timepoint) if xcols else "")
        med = o.groupby("sample_key").raw_sum.median()
        n_sample_cells = o.groupby("sample_key").size()
        u = o[o.comp.notna()].copy()
        u["group"] = u.donor_id + "|" + u.state + "|" + u.comp + "|" + u.tissue_general + "|" + u.assay + "|" + u.disease + (tp(u.timepoint) if xcols else "")
        codes, uniq = pd.factorize(u.group)
        u["gidx"] = codes
        base = int(o.soma_joinid.min())
        lut = np.full(int(o.soma_joinid.max()) - base + 1, -1, dtype=np.int32)
        lut[u.soma_joinid.to_numpy() - base] = codes
        # genes
        var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_name"]).concat().to_pandas()
        vcode, vsym = pd.factorize(var.feature_name.astype(str))
        vbase = int(var.soma_joinid.min())
        vlut = np.full(int(var.soma_joinid.max()) - vbase + 1, -1, dtype=np.int32)
        vlut[var.soma_joinid.to_numpy() - vbase] = vcode
        ng, nsym = len(uniq), len(vsym)
        flt = (f"dataset_id == '{ds_id}' and is_primary_data == True and tissue_type == 'tissue' and donor_id in "
               f"[{', '.join(repr(d) for d in sel)}] and disease in "
               f"[{', '.join(repr(d) for d in sorted(set(o.disease)))}]")
        log(f"{ng} pseudobulks from {len(u)} cells ({len(o)} cells of selected donors); streaming X")
        acc, parts, nb, tl, nnz = None, [], 0, time.time(), 0
        qy = exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=flt))
        for tbl in qy.X("raw").tables():
            r = lut[tbl["soma_dim_0"].to_numpy() - base]
            ok = r >= 0
            if ok.any():
                c = vlut[tbl["soma_dim_1"].to_numpy()[ok] - vbase]
                d = np.rint(tbl["soma_data"].to_numpy()[ok])
                parts.append(sp.coo_matrix((d, (r[ok], c)), shape=(ng, nsym)).tocsr())
                nnz += int(ok.sum())
            nb += 1
            if len(parts) >= 40:
                s = sum(parts[1:], parts[0])
                acc = s if acc is None else acc + s
                parts = []
            if time.time() - tl > 300:
                tl = time.time()
                log(f"batch {nb}, {nnz / 1e6:.0f}M nnz used")
        qy.close()
        if parts:
            s = sum(parts[1:], parts[0])
            acc = s if acc is None else acc + s
    finally:
        census.close()
    dense = np.asarray(acc.todense())
    keep = dense.sum(0) > 0
    dense, genes = dense[:, keep], np.asarray(vsym)[keep]
    dense = dense.astype(np.int32 if dense.max() < 2**31 - 1 else np.int64)
    meta_rows = u.groupby("group", sort=False).agg(
        donor=("donor_id", "first"), state=("state", "first"), comp=("comp", "first"), tissue=("tissue", lambda x: ";".join(sorted(set(x)))),
        tissue_general=("tissue_general", "first"), assay=("assay", "first"), disease=("disease", "first"),
        suspension=("suspension_type", "first"), n_cells=("soma_joinid", "size"), obs_raw_sum=("raw_sum", "sum"),
        sample_key=("sample_key", "first"), **{c: (c, "first") for c in xcols}).reindex(uniq)
    df = pd.DataFrame({
        "donor": f"{short}|" + meta_rows.donor.to_numpy(),
        "compartment": np.where(meta_rows.state == "tumour", "tumour:", "adjacent_normal:") + meta_rows.comp.to_numpy(),
        "stratum": short + "|" + meta_rows.assay.to_numpy(),
        "n_cells": meta_rows.n_cells.to_numpy(), "depth": dense.sum(1).astype(np.int64),
        "median_umi": med.reindex(meta_rows.sample_key).to_numpy(), "n_cells_sample": n_sample_cells.reindex(meta_rows.sample_key).to_numpy(),
        "disease": meta_rows.disease.to_numpy(), "tissue": meta_rows.tissue.to_numpy(),
        "tissue_general": meta_rows.tissue_general.to_numpy(), "suspension": meta_rows.suspension.to_numpy(),
        "tcga": [S.tcga(d, tg) for d, tg in zip(meta_rows.disease, meta_rows.tissue_general)],
        "obs_raw_sum": meta_rows.obs_raw_sum.to_numpy().astype(np.int64),
        **{c.replace(".", "_"): meta_rows[c].to_numpy() for c in xcols}})
    df = df[meta_rows.n_cells.to_numpy() >= MIN_CELLS]
    dense = dense[(meta_rows.n_cells.to_numpy() >= MIN_CELLS)]
    cols = {c: pa.array(df[c].to_numpy()) for c in df.columns}
    cols.update({g: pa.array(dense[:, i]) for i, g in enumerate(genes)})
    # sum of duplicated symbols already done by the symbol factorisation
    tmp = out.with_suffix(".parquet.tmp")
    pq.write_table(pa.table(cols), tmp, compression="zstd")
    tmp.rename(out)
    # donors
    dn = o.groupby("donor_id").agg(
        disease=("disease", lambda x: ";".join(sorted(set(x) - {"normal"}))),
        tissue=("tissue", lambda x: ";".join(sorted(set(x)))), sex=("sex", lambda x: ";".join(sorted(set(x)))),
        development_stage=("development_stage", lambda x: ";".join(sorted(set(x)))),
        assay=("assay", lambda x: ";".join(sorted(set(x)))), n_cells_total=("soma_joinid", "size"),
        has_normal=("state", lambda x: bool((x == "adjacent_normal").any())),
        **{c: (c, lambda x: ";".join(sorted(set(x)))) for c in xcols if c != "sample"}).reset_index()
    dn["tcga"] = [S.tcga(d.split(";")[0], t.split(";")[0]) if d else "" for d, t in zip(dn.disease, dn.tissue)]
    pv = df.pivot_table(index=df.donor.str.split("|").str[-1], columns="compartment", values="n_cells", aggfunc="sum")
    dn = dn.merge(pv.add_prefix("n_").reset_index().rename(columns={"donor": "donor_id", "index": "donor_id"}), on="donor_id", how="left")
    dn.insert(0, "donor", f"{short}|" + dn.donor_id)
    try:  # source-obs fields (treatment, response, time, stage, ...), kept when constant within donor
        ex, oj = read_source_extras(ds_id, None)
        census2 = cc.open_soma(census_version=VERSION)
        try:
            oi = census2["census_data"]["homo_sapiens"].obs.read(
                value_filter=base_f, column_names=["soma_joinid", "observation_joinid", "donor_id"]).concat().to_pandas()
        finally:
            census2.close()
        if oj is not None and len(ex) == len(oj) and len(ex.columns):
            ex["observation_joinid"] = oj
            m = oi.merge(ex, on="observation_joinid")
            m["donor_id"] = m.donor_id.astype(str)
            m = m[m.donor_id.isin(sel)]
            agg = m.drop(columns=["soma_joinid", "observation_joinid"]).groupby("donor_id").agg(
                lambda x: ";".join(sorted(set(map(str, x))))[:200]).reset_index()
            dn = dn.merge(agg.add_prefix("src_").rename(columns={"src_donor_id": "donor_id"}), on="donor_id", how="left")
        else:
            log("source obs not aligned or no extra columns")
    except Exception as e:  # noqa: BLE001
        log("source-obs extras failed:", repr(e)[:200])
    dn.to_csv(OUTDIR / f"{short}_donors.csv", index=False)
    mapping.to_csv(OUTDIR / f"{short}_celltype_mapping.csv", index=False)
    comp = df.groupby("compartment").agg(pseudobulks=("n_cells", "size"), ge50=("n_cells", lambda x: int((x >= 50).sum())))
    stats = dict(dataset_id=ds_id, short=short, tumour_label="malignant" if has_mal else "epithelial", selected_donors=len(sel),
                 pseudobulks=len(df), genes=int(len(genes)), seconds=round(time.time() - t0),
                 compartments=comp.to_dict("index"), file_bytes=out.stat().st_size,
                 timepoints=(df.groupby("timepoint")["sample"].nunique().to_dict() if xcols else None))
    (OUTDIR / f"{short}_stats.json").write_text(json.dumps(stats, indent=1, default=int))
    log(f"done: {len(df)} pseudobulks, {len(genes)} genes, {out.stat().st_size / 1e6:.0f} MB")
    return short, "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*", help="dataset id prefixes (default: all with a qualifying donor)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    r = pd.read_csv(COUNTS)
    q = r[r.donors_qual1 > 0].drop_duplicates("dataset_id")
    nnz = pd.read_csv(COUNTS.parent / "census_nnz_by_dataset.csv").set_index("dataset_id").nnz
    ids = list(q.dataset_id)
    if a.datasets:
        ids = [i for i in ids if any(i.startswith(p) for p in a.datasets)]
    ids = sorted(ids, key=lambda i: -nnz.get(i, 0))  # biggest first
    print(len(ids), "datasets", flush=True)
    with ProcessPoolExecutor(a.workers) as ex:
        for res in ex.map(build_one, [(i, a.threads) for i in ids]):
            print("RESULT", res, flush=True)


if __name__ == "__main__":
    main()
