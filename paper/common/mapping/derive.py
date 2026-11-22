"""Derived tables of the Jed deliveries: depth, size factors, QC and harmonised design, beside the untouched counts.

The delivered files (`$SCCC_DATA/sccc_mapping/<dataset>/`, sha256-verified) are raw counts and stay as delivered;
clocksync reads counts with their depth, so nothing is normalised in place. Written to
`$SCCC_DATA/sccc_mapping_derived/<dataset>/`:
- bulk, `sheet.csv`, one row per run: library (the unit to sum runs over), donor, compartment, age / stage, sex,
  time_h (clock time or ZT, timed sets), group, used (authors' library list where one exists), depth_exon,
  depth_intron, depth_intron_direct, depth_genebody (genome-wide sums), size_factor (depth_exon / dataset median),
  uniquely mapped % and strand;
- single cell, `cells.parquet`, one row per barcode of `<dataset>_cells200` (>= 200 Gene UMIs; the STARsolo
  CellRanger2.2 call is `starsolo_filtered`), or of the `filtered_` files if not delivered: gsm, file, barcode, donor, age, sex, region, total_counts
  (Gene), n_genes, spliced / unspliced / ambiguous totals, unspliced_frac, pct_mito, size_factor (total_counts /
  dataset median), the A244 marker rule applied to the delivered counts (`class_rule`, '' = two classes, none or
  blood; `blood`), and for the adrenal cohorts the A244 cell table: `a244_cell` (barcode in it) and `a244_class` (the
  rule on the authors' counts of that cell, as A244 scored).

Blind-safe: only genome-wide sums, mito and marker genes are read; no clock-gene value is printed or written.

    python -I derive.py [dataset ...]        (default: every delivered dataset)
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[2] / "exploratory" / "A244_fetal_adrenal_classes" / "scripts"))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync.genes import canonical  # noqa: E402
from common.paths import DATA  # noqa: E402
from prep_a244 import CELLS, MARKERS, classify  # noqa: E402

SRC, DEST = DATA / "sccc_mapping", DATA / "sccc_mapping_derived"
CARDOSO = DATA / "natexp" / "cardoso_moreira2019" / "library_table.csv"
A242 = DATA / "natexp" / "fetal_adrenal_open" / "samples.csv"
A244_COHORT = {"kameneva2021_fetal_adrenal": "kameneva", "wang2022_fetal_adrenal": "wang",
               "cheng2022_fetal_adrenal": "sasaki", "dong2020_fetal_adrenal": "dong"}
REGIMEN = {"ad-libitum": "AL", "night-restricted": "NR", "did not adapt": "DR", "arrhythmic": "AR"}
# Delivered `filtered_` files whose GSM label is swapped (2026-11-22; the swap is in the GEO / SRA deposit, the
# `_cells200` files are already relabelled on Jed): the barcodes of filtered_GSM5114030 are those of the
# authors' GSM5114032 (Sample3, 7GW_female1) and vice versa (99.8-100 % shared), and the Y-gene counts follow the
# authors' files and the GEO sex (authors' GSM5114030: RPS4Y1/DDX3Y/KDM5D/UTY/EIF1AY 126,920 UMI; delivered: 4).
SWAP = {"wang2022_fetal_adrenal": {"GSM5114030": "GSM5114032", "GSM5114032": "GSM5114030"}}


def hours(s) -> pd.Series:
    """'12:00' / '7:00' / 'ZT06' / 'T18' -> hours."""
    s = pd.Series(s).astype(str)
    hm = s.str.extract(r"(\d+):(\d+)").astype(float)
    return (hm[0] + hm[1] / 60).fillna(s.str.extract(r"T(\d+)", expand=False).astype(float))


# ---------------------------------------------------------------- bulk design, one parser per dataset
def d_delvalle(s):
    return pd.DataFrame(dict(donor=s.sample_title.str.replace(" ", "_"), compartment=s.sdrf_Characteristics_organism_part,
                             age=s.sdrf_Characteristics_age, age_unit="weeks post conception",
                             sex=s.sdrf_Characteristics_sex))


def d_cardoso(s):
    a = s.sample_alias.str.split(":").str[-1].str.split(".", expand=True)  # <id>.Human.<organ>.<stage>.<sex>
    lib = pd.read_csv(CARDOSO).set_index("library")
    key = a[0].str.extract(r"^(\d+)", expand=False).astype(int)
    t = lib.reindex(key)
    # postnatal libraries have no individual in the SDRF nor in the authors' table: each library is its own donor
    donor = t.individual.astype("Int64").astype(str).where(t.individual.notna(), "lib" + key.astype(str).to_numpy())
    n = a[3].str.extract(r"^(\d+)(w|dpb|mpb|ypb)$")  # 10w (post conception); 226dpb / 6mpb / 13ypb (post birth)
    age = n[0].astype(float) / n[1].map({"w": 1, "dpb": 365.25, "mpb": 12, "ypb": 1})
    unit = n[1].map({"w": "weeks post conception", "dpb": "years post birth", "mpb": "years post birth",
                     "ypb": "years post birth"})
    return pd.DataFrame(dict(library=a[0].to_numpy(), donor=donor.to_numpy(), compartment=a[2].str.lower().to_numpy(),
                             stage=a[3].to_numpy(), age=age.to_numpy(), age_unit=unit.to_numpy(),
                             sex=a[4].str.lower().to_numpy(), used=t.used.eq("Yes").to_numpy()))


def d_greenwell(s):
    reg = s.geo_source_name_ch1.str.lower().map(lambda x: next(v for k, v in REGIMEN.items() if k in x))
    return pd.DataFrame(dict(donor=s.geo_title, compartment="liver", time_h=hours(s.geo_title.str.extract(
        r"(ZT\d+)", expand=False).str.replace("ZT", "T")).to_numpy(), group=reg, sex=s.char_Sex))


def d_mure(s):
    return pd.DataFrame(dict(donor=s.char_time, compartment=s.char_tissue,
                             time_h=hours(s.char_time.str.replace("ZT", "T")).to_numpy()))


def d_perrin(s):
    return pd.DataFrame(dict(library=s.sample_accession, donor=s.geo_title.str.split("_").str[0],
                             compartment="muscle", time_h=hours(s.char_time).to_numpy()))


def d_lundell(s):
    return pd.DataFrame(dict(library=s.geo_description.str.extract(r"Sample ID: (\S+)", expand=False),
                             donor=s.geo_title.str.extract(r"Participant: (\d+)", expand=False), compartment="muscle",
                             time_h=hours(s["char_time point"]).to_numpy(), group=s["char_sample type"]))


def d_dex(s):
    t = s.geo_title.str.extract(r"^[LA](\d+)_(.+?)_(surgical|non-surgical)_(NO_DEX|DEX)$")
    return pd.DataFrame(dict(donor=t[0], compartment=t[1], group=t[2] + "|" + t[3], batch=s.char_batch))


def d_skin(s):
    t = s.geo_title.str.extract(r"Volunteer (\d+), (\w+) skin, .*\((T\d+)\)")
    return pd.DataFrame(dict(donor=t[0], compartment=t[1], time_h=hours(t[2]).to_numpy()))


def d_oral(s):
    return pd.DataFrame(dict(donor=s.char_subject.astype(str), compartment="oral mucosa",
                             time_point=s["char_time point"], age=s.char_age, age_unit="years",
                             sex=s.char_gender.str.lower()))


BULK = {"delvalle2023_fetal_adrenal": d_delvalle, "cardosomoreira2019_prenatal": d_cardoso,
        "cardosomoreira2019_prenatal_v2": d_cardoso, "cardosomoreira2019_postnatal": d_cardoso,
        "greenwell2019_mouse_liver": d_greenwell, "mure2018_baboon": d_mure, "perrin2018_muscle": d_perrin,
        "perrin2018_muscle_trimmed": d_perrin,
        "lundell2020_muscle": d_lundell, "gse320593_liver_fat_dex": d_dex, "gse326530_skin": d_skin,
        "gse186805_oral_mucosa": d_oral}


def bulk(name: str) -> pd.DataFrame:
    d = SRC / name
    s = pd.read_csv(d / "samples.tsv", sep="\t")
    q = pd.read_csv(d / "qc.tsv", sep="\t").set_index("run")
    t = BULK[name](s).reset_index(drop=True)
    t.insert(0, "run", s.run_accession.to_numpy())
    if "library" not in t:
        t.insert(1, "library", t.run)
    for kind in ("exon", "intron", "intron_direct", "genebody"):
        c = pd.read_csv(d / f"counts_{kind}.tsv.gz", sep="\t", index_col=0)
        t[f"depth_{kind}"] = c[t.run].select_dtypes("number").sum().to_numpy()
    lib = t.groupby("library").depth_exon.transform("sum")
    t["size_factor"] = lib / lib.drop_duplicates().median() if t.library.nunique() < len(t) else \
        t.depth_exon / t.depth_exon.median()
    for col in ("uniquely_mapped_pct", "strand_used"):
        if col in q:
            t[col] = q[col].reindex(t.run).to_numpy()
    assert t.donor.notna().all() and t.compartment.notna().all(), f"{name}: unparsed design rows"
    return t


# ---------------------------------------------------------------- single cell
def sc_design(name: str) -> pd.DataFrame:
    s = pd.read_csv(SRC / name / "samples.tsv", sep="\t").set_index("gsm")
    out = pd.DataFrame(index=s.index)
    out["donor"] = name.split("_")[0] + ":" + s.donor_id.astype(str)
    sex = s.filter(regex="(?i)char_sex")
    out["sex"] = sex.iloc[:, 0].astype(str).str.lower().str[0].map({"m": "male", "f": "female"}) if sex.shape[1] \
        else np.nan
    if name in A244_COHORT:  # fetus and age exactly as A242 / A244 (age convention per cohort)
        a = pd.read_csv(A242)
        a = a[a.cohort == A244_COHORT[name]].assign(gsm=lambda x: x.samples.str.split("+")).explode("gsm")
        a["gsm"] = a.gsm.str.split("_").str[0]
        a = a.set_index("gsm").reindex(out.index)
        out["donor"] = a.id.fillna(out.donor)
        out["age"], out["age_unit"] = a.age, a.convention
    elif name.startswith("knightschrijver"):
        out["age"], out["age_unit"] = s["char_developmental stage_(days)"] / 7, "weeks (days / 7)"
        out["region"] = s.char_region
    elif name.startswith("hochane"):
        out["age"], out["age_unit"] = s.donor_id.astype(float), "weeks"
    return out


def sources(name: str) -> list[Path]:
    """The cell matrices of a 10x dataset: `<name>_cells200/cells200_<GSM>.h5ad` (every barcode with >= 200 Gene
    UMIs, labels per GEO with the Wang correction applied on Jed) when delivered, else `filtered_<GSM>.h5ad`."""
    big = sorted((SRC / f"{name}_cells200").glob("cells200_GSM*.h5ad"))
    return big or sorted((SRC / name).glob("filtered_GSM*.h5ad"))


def single_cell(name: str) -> pd.DataFrame:
    design, rows = sc_design(name), []
    a244 = None
    if name in A244_COHORT:
        c = pd.read_parquet(CELLS)
        c = c[c.cohort == A244_COHORT[name]].copy()
        c["a244_class"] = classify(c)[0].to_numpy()
        c["gsm"], c["barcode"] = c["sample"].str.split("_").str[0], c.barcode.str.replace(r"[-.]\d+$", "", regex=True)
        a244 = c.set_index(["gsm", "barcode"]).a244_class
    for f in sources(name):
        gsm = f.stem.split("_")[1]
        if f.stem.startswith("filtered_"):
            gsm = SWAP.get(name, {}).get(gsm, gsm)
        x = ad.read_h5ad(f)
        sym = pd.Index([canonical(g) for g in x.var.symbols.astype(str)])
        X = x.X.tocsr()
        t = pd.DataFrame(dict(gsm=gsm, file=f"{f.parent.name}/{f.name}", barcode=x.obs_names.to_numpy()))
        t["starsolo_filtered"] = x.obs["starsolo_filtered"].to_numpy(bool) if "starsolo_filtered" in x.obs else True
        t["total_counts"] = np.asarray(X.sum(1)).ravel()
        t["n_genes"] = np.diff(X.indptr)
        for layer in ("spliced", "unspliced", "ambiguous"):
            t[layer] = np.asarray(x.layers[layer].sum(1)).ravel()
        t["unspliced_frac"] = t.unspliced / (t.spliced + t.unspliced).replace(0, np.nan)
        t["pct_mito"] = 100 * np.asarray(X[:, sym.str.startswith("MT-")].sum(1)).ravel() / t.total_counts
        m = pd.DataFrame({g: np.asarray(X[:, sym == g].sum(1)).ravel() for g in sum(MARKERS.values(), [])})
        m["total"] = t.total_counts
        cls, met = classify(m)
        t["class_rule"], t["blood"] = cls.to_numpy(), met.blood.to_numpy()
        if a244 is not None:
            hit = a244.reindex(pd.MultiIndex.from_arrays([t.gsm, t.barcode]))
            t["a244_cell"], t["a244_class"] = hit.notna().to_numpy(), hit.fillna("").to_numpy()
        rows.append(t)
        print(f"  {gsm}: {len(t)} cells", flush=True)
    t = pd.concat(rows, ignore_index=True)
    t = t.join(design, on="gsm")
    t["size_factor"] = t.total_counts / t.total_counts.median()
    return t


def main():
    names = sys.argv[1:] or sorted(p.name for p in SRC.iterdir() if (p / "SHA256SUMS").exists())
    for name in names:
        out = DEST / name
        out.mkdir(parents=True, exist_ok=True)
        if name in BULK:
            t = bulk(name)
            t.to_csv(out / "sheet.csv", index=False)
            print(f"{name}: {len(t)} runs, {t.library.nunique()} libraries, {t.donor.nunique()} donors", flush=True)
        elif name.endswith("_cells200"):
            continue  # read through its base dataset
        elif sources(name):
            print(name, flush=True)
            t = single_cell(name)
            t.to_parquet(out / "cells.parquet", index=False)
            msg = f"{name}: {len(t)} cells, {t.gsm.nunique()} GSM, {t.donor.nunique()} donors"
            if "a244_cell" in t:
                n = pd.read_parquet(CELLS, columns=["cohort", "sample"])
                n = n[(n.cohort == A244_COHORT[name]) & n["sample"].str.split("_").str[0].isin(t.gsm)]
                msg += f", A244 cells found {t.a244_cell.sum()} of {len(n)}"
            print(msg, flush=True)
        else:
            print(f"{name}: no parser, skipped", flush=True)


if __name__ == "__main__":
    main()
