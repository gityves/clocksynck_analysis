"""Sample sheets, pair counts and clock-gene detection summaries for the downloaded lesion / multi-tissue datasets.
Only sample metadata and detection (present / median percentile among all genes) are computed; no clock statistic.
Run: python -I build_samples.py [GSE ...]"""
import sys, re, os
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scoutlib import *
OUT = "/home/maxine/Documents/paychere/scCircadianOverview/paper/exploratory/scouts_2026-11-12/output"
COLS = ["dataset","sample_id","title","patient_id","patient_id_source","tissue","status","group","age","sex","platform","data_type","batch","matrix_column"]
SUMMARY = []

def P(g, topic, f): return f"{ROOT}/{topic}/{g}/{f}"
def num(x): return pd.to_numeric(x, errors="coerce")

def finish(g, topic, d, pair_col, pair_vals, extra_note, clock=None, nprobes_genes=None, clock_note=""):
    d = d.copy()
    for c in COLS:
        if c not in d: d[c] = ""
    d["dataset"] = g
    d[COLS].to_csv(f"{OUT}/{topic}_{g}_samples.csv", index=False)
    pat = d.groupby("patient_id")[pair_col].apply(lambda s: set(s.dropna()))
    complete = int(sum(set(pair_vals) <= v for v in pat.values))
    row = dict(dataset=g, topic=topic, n_samples=len(d), n_patients=d.patient_id.nunique(), pair_def=" + ".join(pair_vals),
               n_patients_complete_pairs=complete, note=extra_note)
    if clock is not None:
        clock.to_csv(f"{OUT}/{topic}_{g}_clockgenes.csv", index=False)
        row.update(clock_genes_present=int(clock.present.sum()), median_pct_over_genes=float(clock.median_percentile.median()) if clock.median_percentile.notna().any() else np.nan,
                   n_genes_ranked=nprobes_genes, clock_note=clock_note)
    SUMMARY.append(row); print(g, row, flush=True)

def array_clock(g, topic, smat, annot, ids=None):
    m = read_matrix(P(g, topic, smat)); sym = read_annot(P(g, topic, annot))
    m.index = m.index.astype(str); sym.index = sym.index.astype(str)
    if ids is not None: m = m[[c for c in m.columns if c in set(ids)]]
    return clock_report(m, sym)

def counts_clock(df):
    """df genes x samples raw counts (symbols as index): log2 CPM, rank among all genes (zeros included)."""
    df = df.apply(pd.to_numeric, errors="coerce"); df = df.groupby(df.index.astype(str)).sum()
    lc = np.log2(df / df.sum(axis=0) * 1e6 + 1)
    return clock_report(lc, None)

# ---------------------------------------------------------------- paired_lesion
def d_GSE43292():
    t = "paired_lesion"; s = read_sheet(P("GSE43292", t, "GSE43292_series_matrix.txt.gz"))
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=s["patient"], patient_id_source="characteristic 'patient' (also in title)",
        tissue="carotid artery", status=s.tissue.map({"Atheroma plaque":"lesion (plaque)","Macroscopically intact tissue":"non-lesion (intact adjacent)"}),
        group=s["disease state"], platform=s.platform_id, data_type="normalised log2 array (series matrix)", matrix_column=s.geo))
    c, n = array_clock("GSE43292", t, "GSE43292_series_matrix.txt.gz", "GPL6244.annot.gz")
    finish("GSE43292", t, d, "status", ["lesion (plaque)","non-lesion (intact adjacent)"], "same endarterectomy; age and sex not recorded", c, n)

def d_GSE57218():
    t = "paired_lesion"; s = read_sheet(P("GSE57218", t, "GSE57218_series_matrix.txt.gz"))
    pid = s.title.str.extract(r"RAAK_(\d+)_")[0]
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=pid, patient_id_source="RAAK number parsed from title",
        tissue="articular cartilage (knee/hip)", status=s["disease state"].map({"OA":"lesion (OA)","Preserved":"non-lesion (preserved, same joint)","Healthy":"healthy donor (unpaired)"}),
        group="end-stage OA arthroplasty", age=s["age (yrs)"], sex=s["Sex"], platform=s.platform_id, data_type="normalised log2 array (series matrix)", matrix_column=s.geo))
    c, n = array_clock("GSE57218", t, "GSE57218_series_matrix.txt.gz", "GPL6947.annot.gz")
    finish("GSE57218", t, d, "status", ["lesion (OA)","non-lesion (preserved, same joint)"], "RAAK; 7 healthy donors unpaired; array batch not in sheet", c, n)

def d_gingiva(g):
    t = "paired_lesion"; s = read_sheet(P(g, t, f"{g}_series_matrix.txt.gz"))
    ph = s[[c for c in s.columns if "phenotype" in c.lower()][0]]
    pid = s.source_name_ch1.str.extract(r"Patient (\d+)")[0]
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=pid, patient_id_source="'Patient N' in source_name; numbering is per series (overlap between GSE10334 and GSE16134 not documented)",
        tissue="gingival papilla", status=ph.str.startswith("Diseased").map({True:"lesion (periodontitis site)",False:"non-lesion (healthy site, same mouth)"}),
        group="periodontitis", platform=s.platform_id, data_type="normalised array (series matrix)", matrix_column=s.geo))
    c, n = array_clock(g, t, f"{g}_series_matrix.txt.gz", "GPL570.annot.gz")
    nrep = d.groupby(["patient_id","status"]).size()
    finish(g, t, d, "status", ["lesion (periodontitis site)","non-lesion (healthy site, same mouth)"],
           f"several sites per patient and status (max {nrep.max()}); age, sex, smoking not in sheet; sites sampled at same visit", c, n)

def d_GSE26051():
    t = "paired_lesion"; s = read_sheet(P("GSE26051", t, "GSE26051_series_matrix.txt.gz"))
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=s.donor, patient_id_source="characteristic 'donor'",
        tissue=s.tissue, status=s["disease state"].map({"Lesional":"lesion (tendinopathy)","Non-lesional":"non-lesion"}), group="tendinopathy",
        age=s.age, sex=s.gender, platform=s.platform_id, data_type="normalised array (series matrix)", matrix_column=s.geo))
    c, n = array_clock("GSE26051", t, "GSE26051_series_matrix.txt.gz", "GPL570.annot.gz")
    loc = d.groupby("patient_id").tissue.nunique()
    finish("GSE26051", t, d, "status", ["lesion (tendinopathy)","non-lesion"], f"{int((loc>1).sum())} donors have different tendon labels in the two samples", c, n)

def d_GSE151243():
    t = "paired_lesion"; s = read_sheet(P("GSE151243", t, "GSE151243_series_matrix.txt.gz"))
    cnt = pd.read_csv(P("GSE151243", t, "GSE151243_Raw_gene_counts_matrix.txt.gz"), sep="\t", index_col=0)
    colmap = {c: re.match(r"HidradenitisSuppurativa_(\d+[LN])", c).group(1) for c in cnt.columns}
    s["matrix_column"] = [next(c for c, k in colmap.items() if k == re.search(r"_(\d+[LN])$", ti).group(1)) for ti in s.title]
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=s.patient_id, patient_id_source="characteristic 'patient_id'",
        tissue="skin", status=s.tissue_type.map({"Lesion":"lesion (HS)","Perilesion":"non-lesion (perilesional)"}), group="hidradenitis suppurativa",
        age=s.age, sex=s.gender, platform=s.platform_id, data_type="raw counts (Ensembl ids)", batch="library suffix S/L in column name", matrix_column=s.matrix_column))
    em = ensembl_map(); cnt.index = [em.get(i.split(".")[0], "NA") for i in cnt.index]; cnt = cnt[cnt.index != "NA"]
    c, n = counts_clock(cnt)
    finish("GSE151243", t, d, "status", ["lesion (HS)","non-lesion (perilesional)"], "perilesional (not distant) skin; skin, lower priority", c, n, "log2 CPM, Ensembl ids mapped with local GENCODE v29")

def d_GSE192352():
    t = "paired_lesion"; s = read_sheet(P("GSE192352", t, "GSE192352_series_matrix.txt.gz"))
    cnt = pd.read_csv(P("GSE192352", t, "GSE192352_RAW_counts.csv.gz"), sep=";", quotechar='"', low_memory=False)
    sc = [c for c in cnt.columns if re.fullmatch(r"[CM]E\d+(_.*)?", c)]
    # matrix columns: CE = leiomyoma?, ME = myometrium?; resolve through the sheet titles / descriptions
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=s.patient, patient_id_source="characteristic 'patient'",
        tissue="uterus", status=s.tissue.map({"Uterine leiomyoma":"lesion (leiomyoma, benign tumour)","Uterine myometrium":"non-lesion (myometrium)"}),
        group="Caucasian women (leiomyoma)", age=s.age, sex="F", platform=s.platform_id, data_type="raw counts (gene symbols)",
        batch="", matrix_column=s.title))
    d["bmi"] = s.get("bmi")
    cnt2 = cnt.set_index("Geneid")[sc]
    miss = [x for x in d.matrix_column if x not in cnt2.columns]
    c, n = counts_clock(cnt2)
    finish("GSE192352", t, d, "status", ["lesion (leiomyoma, benign tumour)","non-lesion (myometrium)"],
           f"benign tumour; {len(miss)} sheet titles not found as count columns; count columns {len(sc)}", c, n, "log2 CPM")

def d_GSE119717():
    t = "paired_lesion"; s = read_sheet(P("GSE119717", t, "GSE119717_series_matrix.txt.gz"))
    def dep(ti):
        if ti.startswith("dilated"): return "lesion (dilated AAA wall PVAT)"
        if ti.startswith("non-dilated"): return "non-lesion (non-dilated neck PVAT)"
        if ti.startswith("omental"): return "omental visceral fat"
        return "subcutaneous abdominal fat"
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=s["subject id"], patient_id_source="characteristic 'subject id'",
        tissue="perivascular / omental / subcutaneous adipose", status=s.title.map(dep), group="abdominal aortic aneurysm, all male, Italy",
        age=s["age (years)"], sex="M", platform=s.platform_id, data_type="normalised log2 array (series matrix)", matrix_column=s.geo))
    c, n = array_clock("GSE119717", t, "GSE119717_series_matrix.txt.gz", "GPL10558.annot.gz")
    finish("GSE119717", t, d, "status", ["lesion (dilated AAA wall PVAT)","non-lesion (non-dilated neck PVAT)"],
           "perivascular fat, not aortic wall; arrays duplicated for subjects 50 and 51; also 4-depot multi-tissue", c, n)

def d_GSE39491():
    t = "paired_lesion"; s = read_sheet(P("GSE39491", t, "GSE39491_series_matrix.txt.gz"))
    pid = s.title.str.extract(r"Case \w+ (\d+)")[0]; ct = s["cell type"]
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=pid, patient_id_source="case number parsed from title",
        tissue=ct.map({"BE":"Barrett oesophagus","NE":"squamous oesophagus","NC":"gastric cardia"}), status=ct.map({"BE":"lesion (Barrett metaplasia, premalignant)","NE":"non-lesion (squamous oesophagus, same patient)","NC":"non-lesion (gastric cardia, same patient)"}),
        group="Barrett patients", platform=s.platform_id, data_type="normalised array (series matrix)", matrix_column=s.geo))
    c, n = array_clock("GSE39491", t, "GSE39491_series_matrix.txt.gz", "GPL571.annot.gz")
    finish("GSE39491", t, d, "status", ["lesion (Barrett metaplasia, premalignant)","non-lesion (squamous oesophagus, same patient)"], "premalignant lesion; three tissues per patient; age and sex not in sheet", c, n)

# ---------------------------------------------------------------- living_multitissue
def d_GSE73439():
    t = "living_multitissue"; s = read_sheet(P("GSE73439", t, "GSE73439-GPL5175_series_matrix.txt.gz"))
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=s["subject id"], patient_id_source="characteristic 'subject id'",
        tissue=s.region.map({"Subcutaneous":"subcutaneous adipose","Visceral":"visceral (omental) adipose"}), status=s["subject status"], group=s["subject status"],
        sex="F", platform=s.platform_id, data_type="log2 array gene-level (HuEx, series matrix)", matrix_column=s.geo))
    m = read_matrix(P("GSE73439", t, "GSE73439-GPL5175_series_matrix.txt.gz")); m.index = m.index.astype(str)
    sym = read_annot(P("GSE73439", t, "GPL5175_family.soft.gz")); sym.index = sym.index.astype(str)
    c, n = clock_report(m, sym)
    finish("GSE73439", t, d, "tissue", ["subcutaneous adipose","visceral (omental) adipose"],
           "pregnant (elective c-section +/- labour) and non-pregnant women; not an obesity/T2D cohort; BMI, age not recorded", c, n)

def d_GSE24883():
    t = "living_multitissue"; s = read_sheet(P("GSE24883", t, "GSE24883_series_matrix.txt.gz"))
    idx = s.title.str.extract(r"_(\d+)$")[0]; grp = s["agent"].str.lower()
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=grp.str.replace(" ", "_") + "_" + idx,
        patient_id_source="UNVERIFIED: group + running index in title (GEO says paired SAT/VAT, no explicit id)",
        tissue=s.tissue, status=grp, group=grp, sex="F", platform=s.platform_id, data_type="two-colour log ratio vs reference pool (series matrix)", matrix_column=s.geo))
    sym = read_annot(P("GSE24883", t, "GPL4133.annot.gz")); c, n = presence_only(sym)
    finish("GSE24883", t, d, "tissue", ["subcutaneous adipose tissue","visceral adipose tissue"],
           "32 women, 4 groups x 8; patient link by title index only (cannot be checked by sex genes: all women)", c, n, "presence only; log ratio vs pool has no meaningful percentile")

def d_GSE320593():
    t = "living_multitissue"; s = read_sheet(P("GSE320593", t, "GSE320593_series_matrix.txt.gz"))
    pat = s.title.str.split("_").str[0]; surg = s.title.str.extract(r"(non-surgical|surgical)")[0]
    lib = s.description.str.replace("Library name: ", "")
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=pat, patient_id_source="title prefix (code shared by the two biopsies of one patient)",
        tissue=s.tissue, status=surg, group=s.treatment, platform=s.platform_id, data_type="raw counts (gene symbols)", batch=s.batch, matrix_column=lib))
    at = pd.read_csv(P("GSE320593", t, "GSE320593_AT.counts.csv.gz"), index_col=0); lv = pd.read_csv(P("GSE320593", t, "GSE320593_liver.counts.csv.gz"), index_col=0)
    c, n = counts_clock(at)
    finish("GSE320593", t, d, "status", ["surgical","non-surgical"], "NOT multi-tissue: no patient has both liver and SAT; pairs are surgical vs non-surgical biopsy of one tissue (DEX confound; batch A-D)", c, n, "log2 CPM, adipose file")
    c2, n2 = counts_clock(lv); c2.to_csv(f"{OUT}/{t}_GSE320593_clockgenes_liver.csv", index=False)

def d_GSE133979():
    t = "living_multitissue"; s = read_sheet(P("GSE133979", t, "GSE133979_series_matrix.txt.gz"))
    pat = s.title.str.extract(r"^([A-Z]\d+)")[0]; tis = s.tissue.map({"Rectus abdominis muscle":"skeletal muscle (rectus abdominis)","Subcutaneous adipose":"subcutaneous adipose"})
    col = pat + np.where(tis.str.startswith("skeletal"), "_Muscle", "_Adipose")
    d = pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id=pat, patient_id_source="title prefix (S/M/C + number)", tissue=tis, status=s["diagnosis"],
        group=s["diagnosis"], sex=s["Sex"], platform=s.platform_id, data_type="variance-stabilised counts (symbols)", matrix_column=col))
    ad = pd.read_csv(P("GSE133979", t, "GSE133979_Subcutaneous_adipose_variance_stabilized_counts.txt.gz"), sep="\t", index_col=0)
    c, n = clock_report(ad, None)
    finish("GSE133979", t, d, "tissue", ["skeletal muscle (rectus abdominis)","subcutaneous adipose"], "cancer-surgery patients (PDAC vs benign controls); cachexia study; VST not counts", c, n, "VST, adipose file")

def d_greenawalt():
    t = "living_multitissue"; parts = []
    sx = {}
    for g, tis in [("GSE24293","liver"),("GSE24294","omental adipose"),("GSE24297","subcutaneous adipose")]:
        s = read_sheet(P(g, t, f"{g}_series_matrix.txt.gz"))
        idc = [c for c in s.columns if c == "ID"]; idx = s["ID"].str.extract(r"_(\d+)$")[0]
        parts.append(pd.DataFrame(dict(sample_id=s.geo, title=s.title, patient_id="MGH_" + idx, patient_id_source="NOT A PERSON ID: running index per tissue series (MGH_<tissue>_<k>); sex-gene check across tissues fails (see identity_check csv); no linking file in GEO",
            tissue=tis, status="RYGB bariatric", group="morbidly obese, RYGB, MGH 2000-2007", platform=s.platform_id, data_type="two-colour log ratio vs tissue pool (series matrix)", matrix_column=s.geo, series=g)))
    d = pd.concat(parts, ignore_index=True)
    sym = read_annot(P("GSE24293", t, "GPL4372.annot.gz")); c, n = presence_only(sym)
    # identity check with sex genes (no clock genes)
    sexg = {"XIST":["XIST"], "Y":["RPS4Y1","KDM5D","DDX3Y","UTY","EIF1AY"]}
    probes = {}
    for k, gl in sexg.items():
        probes[k] = [p for p, L in sym.items() if any(x in gl for x in L)]
    keep = set(map(str, probes["XIST"] + probes["Y"]))
    scores = {}
    for g, tis in [("GSE24293","liver"),("GSE24294","omental adipose"),("GSE24297","subcutaneous adipose")]:
        m = read_matrix_rows(P(g, t, f"{g}_series_matrix.txt.gz"), keep)
        x = m.loc[[str(p) for p in probes["XIST"] if str(p) in m.index]].mean(); y = m.loc[[str(p) for p in probes["Y"] if str(p) in m.index]].mean()
        sc = (x - y); gs = d[d.series == g]; sc.index = gs.sample_id.values if len(sc) == len(gs) else sc.index
        scores[tis] = pd.Series(sc.values, index=gs.patient_id.values)
    rng = np.random.default_rng(0); chk = []
    for a, b in [("liver","omental adipose"),("liver","subcutaneous adipose"),("omental adipose","subcutaneous adipose")]:
        A, B = scores[a], scores[b]; com = A.index.intersection(B.index)
        r = np.corrcoef(A[com], B[com])[0, 1]
        rp = np.mean([np.corrcoef(A[com], rng.permutation(B[com].values))[0, 1] for _ in range(200)])
        chk.append(dict(pair=f"{a} vs {b}", n_common_index=len(com), r_same_index=r, r_permuted_index_mean=rp, n_probes_xist=len(probes["XIST"]), n_probes_y=len(probes["Y"])))
    pd.DataFrame(chk).to_csv(f"{OUT}/{t}_GSE24293-4-7_identity_check_sexgenes.csv", index=False); print(pd.DataFrame(chk).to_string(), flush=True)
    d = d.drop(columns="series")
    finish("GSE24293-4-7", t, d, "tissue", ["liver","omental adipose","subcutaneous adipose"], "pair count is index-based and NOT valid: same-index samples of different tissues do not share sex (sex-score r ~ 0.05); GSE24335 is the 2200-sample union of the three series (no extra fields)", c, n, "presence only; two-colour ratios")

FUNCS = {"GSE43292":d_GSE43292,"GSE57218":d_GSE57218,"GSE10334":lambda:d_gingiva("GSE10334"),"GSE16134":lambda:d_gingiva("GSE16134"),
         "GSE26051":d_GSE26051,"GSE151243":d_GSE151243,"GSE192352":d_GSE192352,"GSE119717":d_GSE119717,"GSE39491":d_GSE39491,
         "GSE73439":d_GSE73439,"GSE24883":d_GSE24883,"GSE320593":d_GSE320593,"GSE133979":d_GSE133979,"greenawalt":d_greenawalt}
if __name__ == "__main__":
    for w in (sys.argv[1:] or FUNCS):
        try: FUNCS[w]()
        except Exception as e:
            import traceback; traceback.print_exc(); print("FAILED", w, flush=True)
    if len(sys.argv) == 1: pd.DataFrame(SUMMARY).to_csv(f"{OUT}/lesion_multitissue_summary.csv", index=False)
    else: pd.DataFrame(SUMMARY).to_csv(f"{OUT}/lesion_multitissue_summary_partial_{'_'.join(sys.argv[1:])}.csv", index=False)
