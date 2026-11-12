"""Tidy expression matrices (samples x HGNC symbols) and sample sheets for A254 (paired lesion / non-lesion biopsies)
and A255 (two tissues of living persons). Data wrangling only: no clock statistic, no group or tissue comparison.
The only expression-derived checks are shapes, medians of whole matrices, missing shares and a sex-gene concordance
of the person links (XIST vs Y-linked genes).

Conventions
  * rows = samples (index `sample`), columns = gene symbols; ARNTL (not BMAL1) and CIART (not CHRONO / C1orf51);
    first symbol of '///' lists; of several probes / ids of one symbol the one with the highest median is kept.
  * RNA-seq: raw counts as deposited (float32). Arrays: deposited normalised values, log2-transformed when the
    matrix median is > 30 (stated in `scale_note`).
Loaders (GEO series matrix, platform annotation) are reused from the scout in ../../scouts_2026-11-12/scripts.
Downloads are untrusted data: parsed only. Run: python -I prep_lesions_tissues.py [dataset ...]
"""
import sys, os, re, gzip, io, csv, zipfile, hashlib, urllib.request
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../scouts_2026-11-12/scripts"))
from scoutlib import read_sheet, read_matrix, read_annot, CLOCK, ALIAS, ROOT   # noqa: E402

OUT = os.path.normpath(os.path.join(HERE, "../output/prep"))
os.makedirs(OUT, exist_ok=True)
PL, LV = f"{ROOT}/paired_lesion", f"{ROOT}/living_multitissue"
SUMMARY, SEX = [], []


# ----------------------------------------------------------------------------- helpers
def collapse(mat, sym):
    """mat: ids x samples; sym: symbols aligned to the rows of mat (array) or a Series indexed by unique ids.
    Returns samples x symbols (best-median id per symbol)."""
    s = pd.Series(sym.reindex(mat.index).values if isinstance(sym, pd.Series) else np.asarray(sym), dtype=object)
    s = s.where(s.notna(), "").astype(str).str.strip().replace(ALIAS)
    ok = np.flatnonzero(~s.isin(["", "---", "nan", "None", "NA"]).values)
    s, mat = s.iloc[ok], mat.iloc[ok]
    med = np.nanmedian(mat.values, axis=1)
    order = np.argsort(-np.nan_to_num(med, nan=-np.inf), kind="stable")
    syms = s.values[order]
    _, first = np.unique(syms, return_index=True)
    keep = order[first]
    out = mat.iloc[keep].T
    out.columns = s.values[keep]
    out = out.sort_index(axis=1)
    out.index.name = "sample"
    return out


def first_symbols(annot):
    """read_annot output (probe -> list of symbols) -> probe -> first symbol."""
    return annot.apply(lambda L: L[0] if len(L) else np.nan)


def to_log2(expr):
    """Arrays: log2 if matrix median > 30 (non-positive values become missing). Returns (expr, note)."""
    med = float(np.nanmedian(expr.values))
    if med > 30:
        v = expr.values.astype("float64"); v[v <= 0] = np.nan
        return pd.DataFrame(np.log2(v), index=expr.index, columns=expr.columns), f"deposited values on linear scale (median {med:.0f}); log2-transformed"
    return expr, f"deposited normalised values already log scale (median {med:.2f}); unchanged"


def array_expr(g, topic, mat_file, annot_file, cols):
    """Series-matrix array -> (samples x symbols, scale note) restricted to sample ids `cols`."""
    m = read_matrix(f"{ROOT}/{topic}/{g}/{mat_file}")
    m.index = m.index.astype(str); m = m[cols].apply(pd.to_numeric, errors="coerce").astype("float32")
    sym = first_symbols(read_annot(f"{ROOT}/{topic}/{g}/{annot_file}")); sym.index = sym.index.astype(str)
    e = collapse(m, sym)
    return to_log2(e)


def sym_of_ensembl():
    t = pd.read_csv(f"{ROOT}/repeat_biopsy/GSE231509/GSE231509_Read_counts.txt.gz", sep="\t", usecols=["gene_id", "gene_name"])
    t["gene_id"] = t.gene_id.str.split(".").str[0]
    return t.drop_duplicates("gene_id").set_index("gene_id").gene_name


def counts_expr(cnt, sym):
    """cnt: ids x samples raw counts; sym: Series id -> symbol."""
    return collapse(cnt.apply(pd.to_numeric, errors="coerce").astype("float32"), sym)


def sex_check(name, expr, sheet, pcol, sexcol=None):
    """Concordance of a sex call (mean z of Y-linked genes minus z(XIST) > 0 = male; Y genes only when XIST is absent) among the samples of one person, and with the
    sheet's sex field when present. Only informative when both sexes are in the cohort."""
    ycols = [g for g in ("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY") if g in expr.columns]
    xcols = [g for g in ("XIST",) if g in expr.columns]
    if not ycols:
        SEX.append(dict(dataset=name, note="no Y-linked gene in matrix")); return
    x = expr[ycols + xcols].astype("float64")
    if expr.values.max() > 100:   # counts: log2 CPM-like on the retained genes
        x = np.log2(x / expr.sum(axis=1).values[:, None] * 1e6 + 1)
    z = (x - x.mean()) / x.std()
    score = z[ycols].mean(axis=1) - (z["XIST"] if xcols else 0)
    call = pd.Series(np.where(score > 0, "M", "F"), index=expr.index)
    s = sheet.set_index("sample")
    per = call.groupby(s.loc[call.index, pcol]).agg(lambda v: len(v) > 1 and v.nunique() == 1)
    n_multi = call.groupby(s.loc[call.index, pcol]).size().gt(1).sum()
    row = dict(dataset=name, persons_checked=int(n_multi), within_person_concordance=round(float(per.sum()) / max(n_multi, 1), 4),
               share_male_calls=round(float((call == "M").mean()), 3))
    if sexcol and sexcol in s and s[sexcol].notna().any():
        ref = s.loc[call.index, sexcol].astype(str).str[0].str.upper()
        ok = ref.isin(["M", "F"])
        row.update(sheet_sex_concordance=round(float((ref[ok] == call[ok]).mean()), 4), n_sex_in_sheet=int(ok.sum()),
                   cohort_both_sexes=bool(ref[ok].nunique() > 1))
    SEX.append(row)


def finish(name, expr, sheet, kind, scale_note, pair_col=None, pair_vals=None, pcol="person", sexcol=None, extra_note=""):
    """Write parquet + sheet and append the summary row. Pairs = persons with every value of `pair_vals`."""
    sheet = sheet.reset_index(drop=True)
    expr = expr.loc[sheet["sample"]]
    expr.astype("float32").to_parquet(f"{OUT}/{name}_expr.parquet")
    sheet.to_csv(f"{OUT}/{name}_sheet.csv", index=False)
    complete = None
    if pair_col:
        have = sheet.groupby(pcol)[pair_col].apply(set)
        complete = int(sum(set(pair_vals) <= v for v in have))
    hit = [c for c in CLOCK if c in expr.columns]
    SUMMARY.append(dict(dataset=name, samples=len(expr), persons=sheet[pcol].nunique(), genes=expr.shape[1], kind=kind,
                        matrix_median=round(float(np.nanmedian(expr.values)), 3), missing_share=round(float(expr.isna().values.mean()), 5),
                        clock_genes_present=len(hit), persons_complete_pairs=complete, scale_note=scale_note, note=extra_note))
    if sexcol is not None or name.startswith("GSE"):
        sex_check(name, expr, sheet, pcol, sexcol)
    print(SUMMARY[-1], flush=True)


def scout(topic, g):
    return pd.read_csv(os.path.join(HERE, f"../../scouts_2026-11-12/output/{topic}_{g}_samples.csv"), dtype=str)


def keep_pairs(sh, lesion_vals, non_vals):
    ok = sh.groupby("person").condition.apply(lambda c: bool(set(c) & set(lesion_vals)) and bool(set(c) & set(non_vals)))
    return sh[sh.person.map(ok)]


def lesion_sheet(g, cond, person, extra=("age", "sex", "group", "batch")):
    s = scout("paired_lesion", g)
    sh = pd.DataFrame(dict(sample=s.sample_id, matrix_column=s.matrix_column, person=person(s), condition=cond(s)))
    for c in extra:
        if c in s and s[c].notna().any(): sh[c] = s[c]
    return sh


# ----------------------------------------------------------------------------- part 1: paired lesions
def gingiva(g):
    sh = lesion_sheet(g, lambda s: np.where(s.status.str.startswith("lesion"), "lesion", "nonlesion"), lambda s: s.patient_id, ("group",))
    sh = keep_pairs(sh, ["lesion"], ["nonlesion"])
    e, note = array_expr(g, "paired_lesion", f"{g}_series_matrix.txt.gz", "GPL570.annot.gz", list(sh["sample"]))
    sh["kind"], sh["scale_note"] = "array", note
    finish(g, e, sh, "array", note, "condition", ["lesion", "nonlesion"], extra_note="Affymetrix HG-U133 Plus 2; several sites per person and condition")


def d_GSE43292():
    sh = lesion_sheet("GSE43292", lambda s: np.where(s.status.str.startswith("lesion"), "lesion", "nonlesion"), lambda s: s.patient_id, ("group",))
    e, note = array_expr("GSE43292", "paired_lesion", "GSE43292_series_matrix.txt.gz", "GPL6244.annot.gz", list(sh["sample"]))
    sh["kind"], sh["scale_note"] = "array", note
    finish("GSE43292", e, sh, "array", note, "condition", ["lesion", "nonlesion"], extra_note="Affymetrix HuGene 1.0 ST, carotid plaque vs intact")


def d_GSE57218():
    s = scout("paired_lesion", "GSE57218")
    sh = lesion_sheet("GSE57218", lambda s: s.status.map({"lesion (OA)": "lesion", "non-lesion (preserved, same joint)": "nonlesion"}),
                      lambda s: s.patient_id)
    sh = sh[sh.condition.notna()]
    sh = keep_pairs(sh, ["lesion"], ["nonlesion"])
    e, note = array_expr("GSE57218", "paired_lesion", "GSE57218_series_matrix.txt.gz", "GPL6947.annot.gz", list(sh["sample"]))
    sh["kind"], sh["scale_note"] = "array", note
    finish("GSE57218", e, sh, "array", note, "condition", ["lesion", "nonlesion"], sexcol="sex", extra_note="Illumina HT-12 v4; unpaired healthy donors dropped")


def d_GSE39491():
    sh = lesion_sheet("GSE39491", lambda s: s.status.map({"lesion (Barrett metaplasia, premalignant)": "lesion",
        "non-lesion (squamous oesophagus, same patient)": "nonlesion_squamous", "non-lesion (gastric cardia, same patient)": "nonlesion_cardia"}),
        lambda s: s.patient_id, ())
    sh = keep_pairs(sh, ["lesion"], ["nonlesion_squamous", "nonlesion_cardia"])
    e, note = array_expr("GSE39491", "paired_lesion", "GSE39491_series_matrix.txt.gz", "GPL571.annot.gz", list(sh["sample"]))
    sh["kind"], sh["scale_note"] = "array", note
    finish("GSE39491", e, sh, "array", note, "condition", ["lesion", "nonlesion_squamous"], extra_note="pairs counted as lesion + squamous; cardia as second control")


def d_GSE192352():
    f = f"{PL}/GSE192352/GSE192352_RAW_counts.csv.gz"
    hdr = pd.read_csv(f, sep=";", nrows=0).columns
    sc = [c for c in hdr if re.fullmatch(r"[CM]E\d+(_S\d+)?", c)]
    cnt = pd.read_csv(f, sep=";", usecols=["Geneid"] + sc).set_index("Geneid")
    sh = lesion_sheet("GSE192352", lambda s: np.where(s.status.str.startswith("lesion"), "lesion", "nonlesion"), lambda s: s.patient_id, ("age", "bmi", "group"))
    # sheet titles ME10_2_S17 / ME10_3_S22 are the count columns ME10_S17 / ME10_S22
    sh["matrix_column"] = sh.matrix_column.str.replace(r"_\d(?=_S\d+$)", "", regex=True)
    assert set(sh.matrix_column) <= set(cnt.columns), set(sh.matrix_column) - set(cnt.columns)
    sh = keep_pairs(sh, ["lesion"], ["nonlesion"])
    cnt = cnt[sh.matrix_column.values]; cnt.columns = sh["sample"].values
    e = counts_expr(cnt, cnt.index.values)
    sh["kind"], sh["scale_note"] = "counts", "raw gene counts as deposited (symbols)"
    finish("GSE192352", e, sh, "counts", sh.scale_note.iloc[0], "condition", ["lesion", "nonlesion"],
           extra_note="leiomyomas of one patient are separate lesion samples (patient 10: three)")


def d_GSE151243():
    f = f"{PL}/GSE151243/GSE151243_Raw_gene_counts_matrix.txt.gz"
    cnt = pd.read_csv(f, sep="\t", index_col=0)
    sh = lesion_sheet("GSE151243", lambda s: np.where(s.status.str.startswith("lesion"), "lesion", "nonlesion"), lambda s: s.patient_id, ("age", "sex", "group", "batch"))
    sh = keep_pairs(sh, ["lesion"], ["nonlesion"])
    cnt = cnt[sh.matrix_column.values]; cnt.columns = sh["sample"].values
    sy = sym_of_ensembl().reindex(cnt.index.astype(str).str.split(".").str[0]).values
    e = counts_expr(cnt, sy)
    sh["kind"], sh["scale_note"] = "counts", "raw gene counts as deposited; Ensembl ids mapped with the GSE231509 gene table"
    finish("GSE151243", e, sh, "counts", sh.scale_note.iloc[0], "condition", ["lesion", "nonlesion"], sexcol="sex",
           extra_note=f"Ensembl ids mapped: {int(pd.notna(sy).sum())} of {len(sy)}; perilesional skin")


# ----------------------------------------------------------------------------- part 2: living multi-tissue
def d_GSE73439():
    g = "GSE73439"; s = scout("living_multitissue", g)
    sh = pd.DataFrame(dict(sample=s.sample_id, person=s.patient_id, tissue=s.tissue.map({"subcutaneous adipose": "sc_fat", "visceral (omental) adipose": "visceral_fat"}), group=s.group, sex=s.sex))
    sh = sh[sh.groupby("person").tissue.transform("nunique") == 2]
    m = read_matrix(f"{LV}/{g}/GSE73439-GPL5175_series_matrix.txt.gz")
    m.index = m.index.astype(str); m = m[list(sh["sample"])].apply(pd.to_numeric, errors="coerce").astype("float32")
    # platform table of the 4.5 GB family file: stream to '!platform_table_end' only, keep ID and gene_assignment
    sym = {}
    with io.TextIOWrapper(gzip.open(f"{LV}/{g}/GPL5175_family.soft.gz", "rb"), errors="replace") as fh:
        for l in fh:
            if l.startswith("!platform_table_begin"): break
        cols = next(fh).rstrip("\n").split("\t"); ga = cols.index("gene_assignment")
        for l in fh:
            if l.startswith("!platform_table_end"): break
            p = l.rstrip("\n").split("\t")
            if len(p) > ga and p[0] in m.index:
                ent = [e.split("//") for e in p[ga].split("///") if "//" in e]
                ent = [x[1].strip() for x in ent if len(x) > 1 and x[1].strip() not in ("", "---")]
                if ent: sym[p[0]] = ent[0]
    e, note = to_log2(collapse(m, pd.Series(sym)))
    sh["kind"], sh["scale_note"] = "array", note
    finish(g, e, sh, "array", note, "tissue", ["sc_fat", "visceral_fat"], extra_note="Affymetrix HuEx 1.0 ST, transcript-cluster level; symbol = first of gene_assignment; pregnant and non-pregnant women")


def d_GSE24883():
    g = "GSE24883"; s = scout("living_multitissue", g)
    sh = pd.DataFrame(dict(sample=s.sample_id, person=s.patient_id, tissue=s.tissue.map({"subcutaneous adipose tissue": "sc_fat", "visceral adipose tissue": "visceral_fat"}), group=s.group, sex="F"))
    sh = sh[sh.groupby("person").tissue.transform("nunique") == 2]
    m = read_matrix(f"{LV}/{g}/GSE24883_series_matrix.txt.gz")
    m.index = m.index.astype(str); m = m[list(sh["sample"])].apply(pd.to_numeric, errors="coerce").astype("float32")
    sym = first_symbols(read_annot(f"{LV}/{g}/GPL4133.annot.gz")); sym.index = sym.index.astype(str)
    e, note = to_log2(collapse(m, sym))
    note = "two-colour log2 ratio vs reference pool; " + note
    sh["kind"], sh["scale_note"] = "array", note
    finish(g, e, sh, "array", note, "tissue", ["sc_fat", "visceral_fat"], extra_note="Agilent 4x44K; person = group + running index in the title (UNVERIFIED, all women so not checkable by sex genes)")


def d_GSE133979():
    g = "GSE133979"; s = scout("living_multitissue", g)
    sh = pd.DataFrame(dict(sample=s.sample_id, person=s.patient_id, tissue=s.tissue.map({"skeletal muscle (rectus abdominis)": "muscle", "subcutaneous adipose": "sc_fat"}),
                           group=s.group, sex=s.sex, col=s.matrix_column))
    files = {tis: pd.read_csv(f"{LV}/{g}/GSE133979_{fn}_variance_stabilized_counts.txt.gz", sep="\t", index_col=0)
             for fn, tis in [("Subcutaneous_adipose", "sc_fat"), ("Rectus_Abdominis_Muscle", "muscle")]}
    sh = sh[[c in files[t_].columns for c, t_ in zip(sh.col, sh.tissue)]]     # a few sheet samples have no matrix column
    sh = sh[sh.groupby("person").tissue.transform("nunique") == 2]
    parts = []
    for tis, t in files.items():
        ss = sh[sh.tissue == tis]; t = t[ss.col.values]; t.columns = ss["sample"].values
        parts.append(collapse(t.astype("float32"), t.index.astype(str).values))
    e = pd.concat(parts)
    note = "variance-stabilised counts (DESeq2-style vst, log2-like); not raw counts (only these deposited); one file per tissue, genes = union"
    sh = sh.drop(columns="col"); sh["kind"], sh["scale_note"] = "vst", note
    finish(g, e, sh, "vst", note, "tissue", ["sc_fat", "muscle"], sexcol="sex", extra_note="cancer-surgery patients; tissue-specific files with different gene sets")


def d_muther():
    """MuTHER (E-TABM-1140): per-person columns of the three processed zips; probes -> symbols with ADF A-MEXP-1171."""
    d = f"{LV}/E-TABM-1140"; adf = f"{d}/A-MEXP-1171.adf.txt"
    if not os.path.exists(adf):
        urllib.request.urlretrieve("https://ftp.ebi.ac.uk/pub/databases/microarray/data/array/MEXP/A-MEXP-1171/A-MEXP-1171.adf.txt", adf)
    sha = hashlib.sha256(open(adf, "rb").read()).hexdigest()
    pv = open(f"{d}/PROVENANCE.md").read()
    if "A-MEXP-1171" not in pv:
        with open(f"{d}/PROVENANCE.md", "a") as f:
            f.write(f"| A-MEXP-1171.adf.txt | https://ftp.ebi.ac.uk/pub/databases/microarray/data/array/MEXP/A-MEXP-1171/A-MEXP-1171.adf.txt (probe annotation, column `Reporter Database Entry[hugo]`; downloaded 2026-11-12) | {os.path.getsize(adf)} | {sha} |\n")
    with open(adf, errors="replace") as fh:
        for i, l in enumerate(fh):
            if l.startswith("Reporter Name"): skip = i; break
    ann = pd.read_csv(adf, sep="\t", skiprows=skip, usecols=["Reporter Name", "Reporter Database Entry[hugo]"], dtype=str, quoting=csv.QUOTE_NONE)
    sym = ann.set_index("Reporter Name")["Reporter Database Entry[hugo]"].str.split("///").str[0]
    persons = {}
    for n, tis in [(2, "fat"), (4, "skin"), (3, "lcl")]:
        with zipfile.ZipFile(f"{d}/E-TABM-1140.processed.{n}.zip") as z, z.open(z.infolist()[0]) as fh:
            m = pd.read_csv(fh, sep="\t", index_col=0, skiprows=[1], low_memory=False)
        pid = pd.Series([c.split(";")[0].split("_")[0] for c in m.columns])
        assert (pid.str.startswith("TWPID")).all()
        m.columns = pid.values; m = m.T.groupby(level=0).mean().T   # one column per person (merged replicates averaged)
        m = m.astype("float32")
        e = collapse(m, sym); e.index.name = "person"
        e.to_parquet(f"{OUT}/muther_{tis}_expr.parquet"); persons[tis] = e.index
        hit = [c for c in CLOCK if c in e.columns]
        SUMMARY.append(dict(dataset=f"muther_{tis}", samples=len(e), persons=len(e), genes=e.shape[1], kind="array",
                            matrix_median=round(float(np.nanmedian(e.values)), 3), missing_share=round(float(e.isna().values.mean()), 5),
                            clock_genes_present=len(hit), persons_complete_pairs=None,
                            scale_note="deposited log2 normalised intensity (as is; normalisation differs between the three files)", note="Illumina HT-12 v3; ADF A-MEXP-1171 hugo symbols"))
        print(SUMMARY[-1], flush=True)
        del m, e
    R = list(csv.reader(open(f"{d}/E-TABM-1140.sdrf.txt"), delimiter="\t")); h = R[0]
    rec = {}
    for r in R[1:]:
        rec.setdefault(r[h.index("Characteristics[Individual]")], (r[h.index("Characteristics[co-twin]")], r[h.index("Characteristics[twin zygosity]")]))
    zy = {"monozygotic": "MZ", "dizygotic": "DZ"}
    allp = sorted(set().union(*[set(v) for v in persons.values()]))
    pe = pd.DataFrame(dict(person=allp, cotwin=[rec.get(p, ("", ""))[0] or "" for p in allp],
                           zygosity=[zy.get(rec.get(p, ("", ""))[1].strip().lower(), "unknown") for p in allp]))
    for t in persons: pe[f"has_{t}"] = pe.person.isin(persons[t])
    pe.to_csv(f"{OUT}/muther_persons.csv", index=False)
    fs, fl, sl = (pe.has_fat & pe.has_skin).sum(), (pe.has_fat & pe.has_lcl).sum(), (pe.has_skin & pe.has_lcl).sum()
    for r in [r for r in SUMMARY if r["dataset"].startswith("muther")]:
        if r["dataset"] in ("muther_fat", "muther_skin"): r["persons_complete_pairs"] = int(fs)
        if r["dataset"] == "muther_lcl": r["persons_complete_pairs"] = int(fl); r["note"] += f"; pairs = lcl+fat; skin+lcl {int(sl)}"
        if r["dataset"] != "muther_lcl": r["note"] += f"; pairs = fat+skin; fat+lcl {int(fl)}"
    print(f"MuTHER persons {len(pe)}; zygosity {pe.zygosity.value_counts().to_dict()}; cotwin in cohort {int(pe.cotwin.isin(pe.person).sum())}", flush=True)


FUNCS = {"GSE16134": lambda: gingiva("GSE16134"), "GSE10334": lambda: gingiva("GSE10334"), "GSE57218": d_GSE57218, "GSE43292": d_GSE43292,
         "GSE192352": d_GSE192352, "GSE151243": d_GSE151243, "GSE39491": d_GSE39491,
         "muther": d_muther, "GSE73439": d_GSE73439, "GSE133979": d_GSE133979, "GSE24883": d_GSE24883}

if __name__ == "__main__":
    todo = sys.argv[1:] or list(FUNCS)
    for w in todo:
        FUNCS[w]()
    tag = "" if len(todo) == len(FUNCS) else "_" + "_".join(todo)
    pd.DataFrame(SUMMARY).to_csv(f"{OUT}/prep_summary{tag}.csv", index=False)
    pd.DataFrame(SEX).to_csv(f"{OUT}/prep_sexcheck{tag}.csv", index=False)
