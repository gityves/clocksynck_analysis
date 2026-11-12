"""Sample sheets + clock-gene detection summaries for the downloaded repeat-biopsy datasets.
Only sample metadata and the presence / median-expression percentile of the 14 clock genes are computed
(no rhythm, phase or before/after statistic). Run: python -I build_sheets.py"""
import glob, gzip, io, os, re, sys, tarfile, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from gsm_meta import read_meta

R = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/repeat_biopsy"
OUT = "/home/maxine/Documents/paychere/scCircadianOverview/paper/exploratory/scouts_2026-11-12/output"
SC = "/tmp/claude-1000/-home-maxine-Documents-paychere-scCircadianOverview/a7711f70-ea70-4615-a0ec-0b08050a8353/scratchpad"
CLOCK = "ARNTL CLOCK NPAS2 PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 DBP RORA TEF CIART".split()
ENS = {l.split("\t")[1].split(".")[0]: l.split("\t")[0] for l in open(f"{SC}/clock_ens.tsv").read().strip().split("\n")}
ENTREZ = dict(zip("406 9575 4862 5187 8864 8863 1407 1408 9572 7376 1628 6095 7008 400970".split(), CLOCK))
COLS = ["sample_id", "matrix_column", "subject_id", "subject_id_derivation", "tissue", "timepoint_label", "elapsed_h",
        "elapsed_basis", "intervention_between", "first_biopsy_clock_time", "group_arm", "age", "sex", "platform", "data_type"]

def read_annot(gpl, col="Gene symbol"):
    lines = gzip.open(f"{R}/platforms/{gpl}.annot.gz", "rt", errors="replace").read().split("\n")
    i = next(k for k, l in enumerate(lines) if l.startswith("ID\t"))
    a = pd.read_csv(io.StringIO("\n".join(lines[i:])), sep="\t", dtype=str, low_memory=False)
    return a.set_index("ID")[col]

def detect(expr, kind):
    """expr: genes x samples DataFrame indexed by gene symbol (or Ensembl / Entrez as per kind). Return clock-gene table."""
    if kind == "ens":
        expr = expr.copy(); expr.index = expr.index.astype(str).str.split(".").str[0]
        names = pd.Series(expr.index).map(ENS).values
    elif kind == "entrez":
        names = pd.Series(expr.index.astype(str).str.replace("_at", "")).map(ENTREZ).values
    else:
        names = pd.Series(expr.index.astype(str)).replace({"BMAL1": "ARNTL", "CHRONO": "CIART"}).values
    med = expr.median(axis=1).values
    s = pd.Series(med, index=pd.Index(names, dtype=object))
    allmed = pd.Series(med)
    rows = []
    for g in CLOCK:
        m = s[s.index == g]
        if len(m) == 0: rows.append((g, False, np.nan, np.nan)); continue
        v = m.max()
        rows.append((g, True, v, 100 * (allmed <= v).mean()))
    return pd.DataFrame(rows, columns=["gene", "present", "median_expr", "median_percentile_among_all_genes"]), len(med)

def collapse_probes(expr, gpl):
    sym = read_annot(gpl).reindex(expr.index)
    ex = expr.assign(_s=sym.values).dropna(subset=["_s"])
    ex["_s"] = ex["_s"].str.split("///")
    ex = ex.explode("_s"); ex = ex[ex["_s"].str.strip() != ""]
    ex["_m"] = ex.drop(columns="_s").median(axis=1)
    ex = ex.sort_values("_m", ascending=False).drop_duplicates("_s").drop(columns="_m")
    return ex.set_index("_s")

def sheet(df, **kw):
    out = pd.DataFrame({c: kw.get(c, np.nan) for c in COLS}); return out[COLS]

def read_table(path, **kw):
    return pd.read_csv(path, sep=kw.pop("sep", "\t"), index_col=0, **kw)

results = {}
def finish(acc, sh, expr, kind, sameday_max_h=48):
    sh.to_csv(f"{OUT}/repeat_{acc}_samples.csv", index=False)
    inmat = sh[sh.matrix_column.isin(expr.columns)] if expr is not None else sh
    sd = inmat[inmat.elapsed_h.fillna(0) <= sameday_max_h]
    # complete pairs: subject (within arm/visit unit) with >= 2 same-window samples
    multi = acc in ("GSE142401", "GSE118280", "GSE120862")   # several visits/arms/legs per subject
    unit = sd.subject_id.astype(str) + "|" + sd.group_arm.fillna("").astype(str) if multi else sd.subject_id.astype(str)
    vc = unit.value_counts(); ok = set(vc[vc >= 2].index)
    n_pairs = int(sd.loc[unit.isin(ok), "subject_id"].nunique()); n_units = len(ok)
    det, ngenes = detect(expr, kind)
    det.to_csv(f"{OUT}/repeat_{acc}_clockgenes.csv", index=False)
    results[acc] = dict(n_samples_sheet=len(sh), n_samples_in_matrix=len(inmat), n_subjects=int(sh.subject_id.nunique()),
                        n_subjects_with_ge2_samples_within_window=n_pairs, n_subject_arm_units_with_ge2=n_units, n_genes=ngenes,
                        clock_present=int(det.present.sum()), pct_range=[round(det.median_percentile_among_all_genes.min(), 1), round(det.median_percentile_among_all_genes.max(), 1)])
    print(acc, results[acc])

# ---- GSE252357
d, _ = read_meta(f"{R}/GSE252357/GSE252357_series_matrix.txt.gz")
ex = pd.read_csv(f"{R}/GSE252357/GSE252357_normalized_counts_protein_coding.csv.gz", index_col=0)
el = {"PRE": 0, "30m": 1.75, "3h": 4.25, "8h": 9.25, "24h": 25.25}
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d["c_subject id"] + "_" + d.c_timepoint, subject_id=d["c_subject id"],
           subject_id_derivation="characteristic 'subject id'", tissue="skeletal muscle (vastus lateralis)", timepoint_label=d.c_timepoint,
           elapsed_h=d.c_timepoint.map(el), elapsed_basis="protocol: PRE biopsy 1 h after breakfast; exercise/rest starts 30 min later, lasts 45 min; later biopsies timed from end of exercise (+30 min, 3 h, 8 h, 24 h)",
           intervention_between=np.where(d.c_group2 == "RE", "resistance exercise 45 min (leg press + extension)", "none (45 min seated rest; time-matched biopsies)"),
           first_biopsy_clock_time="not stated (1 h after standard breakfast, overnight fast)", group_arm=d.c_group2, age=d["c_age years"], sex=d.c_Sex,
           platform="GPL24676 NovaSeq 6000", data_type="DESeq-type normalised counts, protein-coding genes only")
finish("GSE252357", sh, ex, "sym", 9.5)

# ---- GSE118280
d, ex = read_meta(f"{R}/GSE118280/GSE118280_series_matrix.txt.gz")
el = {"8:40": 0.0, "12:20": 11 / 3, "19:00": 31 / 3}
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d.geo_accession, subject_id=d["c_patient no"], subject_id_derivation="characteristic 'patient no'",
           tissue="subcutaneous adipose", timepoint_label=d["c_time point"], elapsed_h=d["c_time point"].map(el),
           elapsed_basis="clock times of the three biopsies stated in GEO design (08:40, 12:20, 19:00)",
           intervention_between="isocaloric diet-day meals (HC breakfast/lunch + HF snack/dinner, or reverse); biopsies at end of 4-wk diet period; two diet visits per subject",
           first_biopsy_clock_time="08:40", group_arm=d.c_diet, age=np.nan, sex=d.c_gender, platform="GPL20103 Affymetrix HuGene 2.0 ST (Brainarray Entrez CDF)",
           data_type="RMA/quantile-normalised log2 array (series matrix)")
finish("GSE118280", sh, ex, "entrez")

# ---- GSE142401
d, ex = read_meta(f"{R}/GSE142401/GSE142401_series_matrix.txt.gz")
t = d.title
sub = t.str.extract(r"^(ctrl|MetS|mets|MetS)\s*participant\s*(\d+)", flags=re.I)
sid = np.where(t.str.lower().str.startswith("ctrl"), "ctrl", "MetS") + "_" + t.str.extract(r"participant\s*(\d+)")[0]
meal = t.str.extract(r",\s*(\w+) meal intervention")[0]
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d.geo_accession, subject_id=sid, subject_id_derivation="parsed from sample title 'ctrl/MetS participant N'",
           tissue="subcutaneous adipose (abdominal; per paper)", timepoint_label=np.where(t.str.contains("pre meal 0"), "0 h fasting", "4 h post meal"),
           elapsed_h=np.where(t.str.contains("pre meal 0"), 0.0, 4.0), elapsed_basis="sample title (0 h / 4 h post meal)",
           intervention_between="high-fat breakfast meal (" + meal + "-based), crossover; both meals in same subject on separate days", first_biopsy_clock_time="not stated (fasting, breakfast study)",
           group_arm=d.c_disease + " / " + meal, age=np.nan, sex=d.c_Sex, platform="GPL10558 Illumina HumanHT-12 v4", data_type="quantile-normalised log2 array (series matrix)")
finish("GSE142401", sh, collapse_probes(ex, "GPL10558"), "sym")

# ---- GSE4386 (characteristics are sparse; sex taken from non-empty values per sample)
d, ex = read_meta(f"{R}/GSE4386/GSE4386_series_matrix.txt.gz")
cc = [c for c in d.columns if c.startswith("c_")]
sex = d[cc].apply(lambda r: next((x for x in r if x in ("male", "female")), np.nan), axis=1)
p = d.title.str.extract(r"^(H\d+)_([SP])_B([12])")
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d.geo_accession, subject_id=p[0], subject_id_derivation="sample title H<NN>_<S|P>_B<1|2>",
           tissue="right atrial tissue (heart)", timepoint_label=np.where(p[2] == "1", "B1 baseline, beginning of off-pump CABG", "B2 end of CABG"),
           elapsed_h=np.nan, elapsed_basis="NOT stated in GEO (duration of off-pump CABG not given)", intervention_between="off-pump CABG under sevoflurane (S) or propofol (P) anaesthesia",
           first_biopsy_clock_time="not stated", group_arm=np.where(p[1] == "S", "sevoflurane", "propofol"), age=np.nan, sex=sex, platform="GPL570 Affymetrix HG-U133 Plus 2", data_type="RMA log2 array (series matrix)")
finish("GSE4386", sh, collapse_probes(ex, "GPL570"), "sym", 1e9)

# ---- GSE151648
d, _ = read_meta(f"{R}/GSE151648/GSE151648_series_matrix.txt.gz")
ex = read_table(f"{R}/GSE151648/GSE151648_liver-iri-counts.txt.gz")
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d.title, subject_id=d.title.str.extract(r"^(\d+)-Bx")[0], subject_id_derivation="sample title '<patient>-Bx1/Bx2'",
           tissue="liver (needle core; allograft)", timepoint_label=d.c_timepoint, elapsed_h=np.nan, elapsed_basis="NOT stated in GEO (pre-implantation vs post-reperfusion, same transplant operation)",
           intervention_between="orthotopic liver transplantation: cold ischaemia, implantation, reperfusion (TRANSPLANT; donor graft)", first_biopsy_clock_time="not stated",
           group_arm=d.c_iri, age=np.nan, sex=np.nan, platform="GPL21290 HiSeq 3000", data_type="raw counts (featureCounts)")
finish("GSE151648", sh, ex, "ens", 1e9)

# ---- GSE75432
d, _ = read_meta(f"{R}/GSE75432/GSE75432_series_matrix.txt.gz")
cols = {}
with tarfile.open(f"{R}/GSE75432/GSE75432_RAW.tar") as tf:
    for m in tf.getmembers():
        g = m.name.split("_")[0]
        cols[g] = pd.read_csv(io.BytesIO(gzip.decompress(tf.extractfile(m).read())), sep="\t", header=None, index_col=0)[1]
ex = pd.DataFrame(cols)
post = d.c_tissue.str.contains("2hr")
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d.geo_accession, subject_id=d.title.str.extract(r"^(\d+)-")[0], subject_id_derivation="sample title '<patient>-1 (pre) / -3 (post)'",
           tissue="skeletal muscle (vastus lateralis)", timepoint_label=np.where(post, "2 h after surgery with tourniquet", "immediately before surgery"),
           elapsed_h=np.where(post, 2.0, 0.0), elapsed_basis="label says 2 h after surgery+tourniquet; true interval = surgery duration + 2 h is not given, so 2.0 is a lower bound",
           intervention_between="total knee arthroplasty with tourniquet (ischaemia), general or spinal anaesthesia not given in GEO", first_biopsy_clock_time="not stated", group_arm="TKA tourniquet",
           age=np.nan, sex=np.nan, platform="GPL18573 NextSeq 500", data_type="raw counts (per-sample count files in RAW.tar)")
finish("GSE75432", sh, ex, "ens", 1e9)

# ---- GSE231509
d, _ = read_meta(f"{R}/GSE231509/GSE231509_series_matrix.txt.gz")
ex = pd.read_csv(f"{R}/GSE231509/GSE231509_Read_counts.txt.gz", sep="\t", index_col=1).drop(columns=["gene_id", "gene_biotype"])
ex = ex[~ex.index.duplicated()]
aft = d["c_time point"].str.contains("after")
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d.title, subject_id=d.c_diagnosis.str.replace(r"\W+", "", regex=True) + "_" + d.c_individual,
           subject_id_derivation="characteristic 'individual' within diagnosis group (individual numbers repeat across groups)", tissue="skeletal muscle (vastus lateralis)",
           timepoint_label=d["c_time point"], elapsed_h=np.where(aft, 1.0, 0.0), elapsed_basis="characteristic 'time point' (1 hour after MMTT)",
           intervention_between="mixed-meal tolerance test (Resource 2.0, 3 ml/kg)", first_biopsy_clock_time="arrival 09:00 after 12 h overnight fast; biopsy clock time not stated",
           group_arm=d.c_diagnosis, age=d.c_age, sex=d.c_Sex, platform="GPL21697 NextSeq 550", data_type="raw read counts")
finish("GSE231509", sh, ex, "sym", 1e9)

# ---- GSE273041
d, _ = read_meta(f"{R}/GSE273041/GSE273041_series_matrix.txt.gz")
ex = pd.read_csv(f"{R}/GSE273041/GSE273041_raw_counts.txt.gz", sep="\t", index_col=0)
mc = {c.split("_", 1)[0]: c for c in ex.columns}
tp = d.c_timepoint
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d.title.map(mc), subject_id=d["c_participant id"].str.split("_").str[0], subject_id_derivation="characteristic 'participant id' (<n>_<time>)",
           tissue="skeletal muscle (vastus lateralis; per paper)", timepoint_label=tp.map({"Base": "BL rest", "IP": "H0 immediately post 3 h run", "IP120": "H2 two hours into recovery"}),
           elapsed_h=tp.map({"Base": 0.0, "IP": 3.0, "IP120": 5.0}), elapsed_basis="protocol: 3 h treadmill run then +2 h; assumes BL biopsy immediately before the run (lower bound)",
           intervention_between="3 h submaximal treadmill run", first_biopsy_clock_time="not stated", group_arm=d.c_diet, age=np.nan, sex="male (all)", platform="GPL16791 HiSeq 2500",
           data_type="counts (non-integer; normalised, 'raw_counts' file)")
finish("GSE273041", sh, ex, "ens", 1e9)

# ---- GSE151066
d, _ = read_meta(f"{R}/GSE151066/GSE151066_series_matrix.txt.gz")
ex = pd.read_csv(f"{R}/GSE151066/GSE151066_rsem_genes_count.txt.gz", sep=r"\s+", quotechar='"')
tp = d["c_time.point"]
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d.title, subject_id=d["c_subject.id"], subject_id_derivation="characteristic 'subject.id'",
           tissue="skeletal muscle (vastus lateralis)", timepoint_label=tp.map({"Pre": "pre", "Post": "immediately post", "3hrPost": "3 h post"}),
           elapsed_h=tp.map({"Pre": 0.0, "Post": 40 / 60, "3hrPost": 40 / 60 + 3.0}), elapsed_basis="protocol: 40 min cycling at 60-70% HRR; 3 h after exercise end (biopsy delay ignored)",
           intervention_between="40 min endurance cycling", first_biopsy_clock_time="not stated", group_arm=d.c_cohort, age=np.nan, sex=d.c_gender, platform="GPL24676 NovaSeq 6000",
           data_type="RSEM expected counts")
finish("GSE151066", sh, ex, "ens", 1e9)

# ---- GSE303983
d, _ = read_meta(f"{R}/GSE303983/GSE303983_series_matrix.txt.gz")
ex = pd.read_csv(f"{R}/GSE303983/GSE303983_genes.counts.csv.gz", index_col=0)
ex.index = ex.index.str.split("|").str[-1]
ex = ex[~ex.index.duplicated()]
mc = {}
for c in ex.columns:
    m = re.match(r"(3X-\d+)_(PRE|POST)_", c); mc[(m.group(1), m.group(2))] = c
tt = d.title.str.extract(r"^(3X-\d+)_(PRE|POST)")
post = tt[1] == "POST"
sh = sheet(d, sample_id=d.geo_accession, matrix_column=[mc.get((a, b)) for a, b in zip(tt[0], tt[1])], subject_id=tt[0], subject_id_derivation="sample title '3X-<NN>_PRE/POST'",
           tissue="abdominal subcutaneous adipose", timepoint_label=d.c_time, elapsed_h=np.where(post, 1.5, 0.0),
           elapsed_basis="1.5 h after END of exercise (exercise 20-60 min by arm, so interval from first biopsy is longer); arm of each subject not in GEO metadata",
           intervention_between="one exercise session: LOW 60 min 30% VO2max, MOD 45 min 65%, or HIGH 10x1 min 90% (n=15 each)", first_biopsy_clock_time="not stated", group_arm="not in GEO metadata",
           age=np.nan, sex=np.nan, platform="GPL24676 NovaSeq 6000", data_type="raw gene counts (STAR/featureCounts-like)")
finish("GSE303983", sh, ex, "sym", 1e9)

# ---- GSE120862
ds = [read_meta(f"{R}/GSE120862/GSE120862-{g}_series_matrix.txt.gz")[0] for g in ("GPL16791", "GPL18573")]
d = pd.concat(ds, ignore_index=True)
ex = pd.read_csv(f"{R}/GSE120862/GSE120862_DESeq_normalized_counts.txt.gz", sep="\t", index_col=0)
ex = ex[~ex.index.duplicated()]
tp = d["c_time point"].str.extract(r"(baseline|\+1 h|\+4 h)")[0]
leg = d["c_exercised/non-exercised muscle"].str.contains("non").map({True: "NE (non-exercised leg, control)", False: "Ex (exercised leg)"})
sh = sheet(d, sample_id=d.geo_accession, matrix_column=d.title, subject_id=d.c_individual, subject_id_derivation="characteristic 'individual'",
           tissue="skeletal muscle (vastus lateralis)", timepoint_label=tp, elapsed_h=tp.map({"baseline": 0.0, "+1 h": 2.0, "+4 h": 5.0}),
           elapsed_basis="protocol: 60 min one-legged exercise, biopsies 1 h and 4 h after its end; assumes baseline biopsy just before exercise (otherwise 1 h / 4 h from exercise end)",
           intervention_between="60 min one-legged knee extension (other leg = within-person no-exercise control); repeated after 2 months training ('trained' rows, a different day)",
           first_biopsy_clock_time="not stated", group_arm=d["c_training status"] + " / " + leg, age=d.c_age, sex=d.c_Sex, platform="GPL16791 HiSeq 2500 + GPL18573 NextSeq 500",
           data_type="DESeq-normalised counts")
finish("GSE120862", sh, ex, "sym", 1e9)

json.dump(results, open(f"{OUT}/repeat_summary.json", "w"), indent=1)
