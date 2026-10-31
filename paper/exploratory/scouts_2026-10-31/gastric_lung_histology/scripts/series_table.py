"""Per-series scout table: counts from the header-derived sample tables (metadata only), platform coverage from platform_coverage.csv, scale from scale_diagnostics.csv."""
from pathlib import Path
import sys
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent)); import parse_hdr as P
R = Path(__file__).resolve().parents[1]; D = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_gastric_lung"
cov = pd.read_csv(R / "platform_coverage.csv").set_index("platform")
sc = pd.read_csv(R / "scale_diagnostics.csv").set_index("series")
def T(name): return P.table(P.parse(R / "hdr" / f"{name}_series_matrix.hdr"))
def n(t, col, vals): return int(t[col].isin(vals).sum())
rows = []
def add(series, group, histology, nT, nN, plat, scale, ann, overlap, dl, cover=None):
    if plat in cov.index:
        c = cov.loc[plat]; cl = f"{c.n_clock}/{c.n_clock_total}"; miss = c.clock_missing if isinstance(c.clock_missing, str) else ""; g = f"{c.n_g2m}/{c.n_g2m_total}"; pu = f"{c.n_purity}/{c.n_purity_total}"
    else: cl = miss = g = pu = ""; miss = cover or "not assessed"
    if dl and series in sc.index and pd.notna(sc.loc[series, "global_max"]):
        r = sc.loc[series]; scale = f"{scale} [whole matrix min {r.global_min:.3g}, median {r.global_median:.3g}, max {r.global_max:.3g}, NaN {r.frac_nan:.2%}]"
    rows.append(dict(series=series, group=group, tumour_histology_counts=histology, n_tumours=nT, n_normal=nN, platform=plat, scale=scale, n_clock_genes=cl, clock_missing=miss,
                     n_g2m_genes=g, n_purity_genes=pu, annotations=ann, overlap_notes=overlap, downloaded="yes" if dl else "no", path=f"{D}/{series}/" if dl else ""))
# ---- gastric
add("GSE62254", "gastric", "gastric adenocarcinoma (ACRG), all primary", 300, 0, "GPL570", "NOT log2-like as deposited: range fits log10 of MAS5 (inference); to be established before use",
    "Only tissue and patient ID in GEO fields. ACRG molecular subtype (MSS/TP53+, MSS/TP53-, MSS/EMT, MSI), Lauren, stage, EBV, site are NOT in the series matrix (PMID 25894828 / 29725014 supplement or authors; unverified)",
    "ACRG (Samsung Medical Center cohort, Cristescu 2015). GSE66229 contains the same 300 tumours + 100 normals (identical patient IDs): use one, never both. Possible patient overlap with GSE26253 (Samsung) not verifiable", True)
add("GSE66229", "gastric", "same 300 ACRG tumours as GSE62254 + 100 normal gastric", 300, 100, "GPL570", "not downloaded (same arrays as GSE62254 expected)",
    "tissue (Gastric tumor / normal gastric tissue) and patient ID only; normals paired for a subset (patient IDs shared)", "OVERLAPS GSE62254 (300 tumours); only adds 100 normals", False)
add("GSE15459", "gastric", "primary gastric tumours, Singapore: 200 arrays, 192 gastric adenocarcinoma after the authors' exclusions (4 failed QC, 4 not adenocarcinoma, listed in series design)", 200, 0, "GPL570",
    "linear MAS5-like", "tissue only in GEO fields; Lauren/stage/survival in GSE15459_outcome.xls (supplementary, not downloaded, not checked); no normals",
    "National Cancer Centre Singapore (Ooi/Tan 2009; Ooi, Tan lab). GSE34942 is the Singapore batch B (56 other tumours, not overlapping by sample names GC-2000xxx vs GC-0xx but not provable)", True)
add("GSE84437", "gastric", "gastric cancer, Yonsei: 483 arrays; 433 with clinical fields, 50 without", 483, 0, "GPL6947", "linear (Illumina HT-12 v3, raw-like; 1.6% NaN)",
    "age, sex, pT, pN, death, OS for 433; no Lauren, subtype, MSI, EBV or site; SuperSeries of Yonsei sub-series; no normals", "Yonsei University Seoul (Yoon/Kim); independent of ACRG by institution, not verifiable", True)
add("GSE26253", "gastric", "FFPE gastric adenocarcinoma stage 1B-4 (stage 2 167, 3A 111, 1B 68, 4 67, 3B 19)", 432, 0, "GPL8432", "log2-like (DASL)",
    "pathological stage, recurrence status and time; no Lauren, subtype, MSI, EBV, site; FFPE degraded RNA; no normals",
    "Samsung (curative resection + adjuvant chemo(radio)therapy, Cho 2011); possible patient overlap with ACRG (Samsung Medical Center) not verifiable", True)
add("GSE26899", "gastric", "gastric adenocarcinoma 93 (+3 GIST), Korea University", 93, 12, "GPL6947", "not downloaded (Illumina HT-12 v3)",
    "Lauren (intestinal 59, diffuse 31, mixed 2), AJCC stage (1 11, 2 18, 3 27, 4 36), location (antrum 55, body 29, cardia 7, fundus 2), age, sex, adjuvant chemo; 12 surrounding normals",
    "Korea University Guro Hospital; <150 tumours", False)
add("GSE29272", "gastric", "tumours 134 (cardia 62, non-cardia 72)", 134, 134, "GPL96", "not downloaded",
    "cardia vs non-cardia site only; 134 paired adjacent normals", "NCI (Hu 2012); <150 tumours; 268 arrays", False)
add("GSE27342", "gastric", "80 tumours", 80, 80, "GPL5175", "not downloaded (Exon ST)", "stage, grade, age, sex; paired adjacent normal for all 80", "<150 tumours", False, cover="not assessed (Exon 1.0 ST)")
add("GSE57303", "gastric", "70 tumours", 70, 0, "GPL570", "not downloaded", "none (phenotype: gastric cancer)", "<150 tumours; Chinese cohort (AstraZeneca Shanghai)", False)
add("GSE34942", "gastric", "56 tumours (Singapore batch B)", 56, 0, "GPL570", "not downloaded", "none", "<150 tumours; Singapore, sibling of GSE15459", False)
# ---- lung
t = pd.read_csv(R / "meta/GSE30219_samples.csv")
add("GSE30219", "lung", f"ADC {n(t,'histology',['ADC'])}, squamous (SQC) {n(t,'histology',['SQC'])}, large-cell NE {n(t,'histology',['LCNE'])}, basaloid {n(t,'histology',['BAS'])}, carcinoid {n(t,'histology',['CARCI'])}, 'SCC' {n(t,'histology',['SCC'])} (very likely small-cell, not squamous: inferred from the CIT publication, not from GEO fields), LCC 3, other 4",
    293, 14, "GPL570", "log2-like (RMA)", "histology, pT/pN/pM, relapse, survival, age, sex; 14 non-tumoural lung; no smoking, no mutations",
    "CIT programme, Grenoble (Rousseaux 2013). Squamous n = 61 passes the 60 threshold barely; basaloid (39) is a separate class", True)
t = pd.read_csv(R / "meta/GSE37745_samples.csv")
add("GSE37745", "lung", f"adeno {n(t,'histology',['adeno'])}, squamous {n(t,'histology',['squamous'])}, large {n(t,'histology',['large'])}", 196, 0, "GPL570", "log2-like (RMA)",
    "histology, stage, WHO performance status, age, sex, adjuvant, recurrence, survival (days); no smoking, no normals", "Uppsala University (Botling 2013), operated 1995-2005; GSE81089 is the Uppsala RNA-seq cohort 2006-2010 (different patients by dates)", True)
t = pd.read_csv(R / "meta/GSE41271_samples.csv")
add("GSE41271", "lung", f"Adenocarcinoma {n(t,'histology',['Adenocarcinoma'])}, squamous {n(t,'histology',['Squamous'])}, other 12", 275, 0, "GPL6884", "log2-like (Illumina WG6 v3)",
    "histology, stage, tobacco history, sex, race, recurrence; no normals", "MD Anderson (Wistuba); copied from scout_luad_lihc (not re-downloaded)", True)
add("GSE81089", "lung", "squamous 67 (code 1), adenocarcinoma 108 (code 2), large cell/NOS 24 (code 3) + 19 non-malignant (coding from the series design text)", 199, 19, "GPL16791", "RNA-seq (Illumina HiSeq 2500); series matrix has NO expression table, FPKM is only in a supplementary file (GSE81089_FPKM_cufflinks.tsv.gz), not downloaded (rule: matrices and platform tables only)",
    "histology (coded), TNM stage (coded), age, sex, smoking (1 current 96, 2 ex 84, 3 never 19: coding as in design text), WHO PS, vital status", "Uppsala (Djureinovic 2016, 2006-2010); the only RNA-seq lung cohort with both histologies", True, cover="RNA-seq, whole-transcriptome (no probe table)")
add("GSE50081", "lung", "adenocarcinoma 127 (+1 favor ADC), squamous 42 (+1), large 7, adenosquamous 2, mixed 1", 181, 0, "GPL570", "not downloaded", "stage I-II, smoking, recurrence, survival", "UHN Toronto; squamous 42-43 < 60", False)
add("GSE19188", "lung", "ADC 45, SCC 27, LCC 19", 91, 65, "GPL570", "not downloaded", "histology, survival, sex; 65 adjacent healthy", "Erasmus MC; squamous 27 < 60", False)
add("GSE33532", "lung", "adeno 40, squamous 16, mixed 24", 80, 20, "GPL570", "not downloaded", "multi-region sampling (several arrays per tumour), stage", "squamous 16 < 60; replicate regions", False)
add("GSE3141", "lung", "adenocarcinoma 58, squamous 53", 111, 0, "GPL570", "not downloaded", "cell type, survival only", "Duke; squamous 53 < 60; also in many meta-analyses", False)
# ---- oesophagus
add("GSE13898", "oesophagus", "EAC 75, no ESCC", 75, 28, "GPL6102", "not downloaded (Illumina HumanWG-6 v2)", "pathology, Barrett/dysplasia", "no squamous carcinoma: fails the both-histologies rule", False, cover="not assessed")
add("GSE26886", "oesophagus", "EAC 21, ESCC 9, Barrett 20, normal squamous 19", 30, 19, "GPL570", "not downloaded", "tissue class only", "ESCC 9 < 50", False)
add("GSE53624", "oesophagus", "ESCC only (119 tumours + 119 paired normal)", 119, 119, "GPL18109", "not downloaded (Agilent lncRNA+mRNA)", "tobacco, alcohol, stage, grade, survival", "ESCC only", False, cover="not assessed")
add("GSE23400", "oesophagus", "ESCC only (53 + 51 tumours, two arrays GPL96 + GPL97)", 104, 104, "GPL96", "not downloaded", "none", "ESCC only; GPL97 is the second array of the U133 set", False)
add("GSE20347", "oesophagus", "ESCC only", 17, 17, "GPL571", "not downloaded", "none", "ESCC only, small", False, cover="not assessed")
df = pd.DataFrame(rows); df.to_csv(R / "series_table.csv", index=False)
pd.set_option("display.width", 250); print(df.drop(columns=["annotations", "overlap_notes", "path", "scale", "clock_missing", "tumour_histology_counts"]).to_string())
