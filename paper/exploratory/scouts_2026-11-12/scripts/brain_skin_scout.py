"""Scout: sample sheets, design tables and clock-gene probe detection for GSE45642, GSE71620, GSE139300/GSE139301.
Parses public series matrices only (no clock statistic). Run: python -I brain_skin_scout.py"""
import gzip, io, re, collections
import numpy as np, pandas as pd

R = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/"
OUT = "/home/maxine/Documents/paychere/scCircadianOverview/paper/exploratory/scouts_2026-11-12/output/"
GENES = "ARNTL CLOCK NPAS2 PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 DBP RORA TEF CIART".split()
ALIAS = {"ARNTL": {"ARNTL", "BMAL1"}, "CIART": {"CIART", "C1ORF51"}}
ENTREZ = dict(ARNTL=406, CLOCK=9575, NPAS2=4862, PER1=5187, PER2=8864, PER3=8863, CRY1=1407, CRY2=1408,
              NR1D1=9572, NR1D2=9975, DBP=1628, RORA=6095, TEF=7008, CIART=148523)


def read_matrix(path):
    meta = collections.defaultdict(list); lines = []; inmat = False
    with gzip.open(path, "rt") as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"): inmat = True; continue
            if l.startswith("!series_matrix_table_end"): break
            if inmat: lines.append(l)
            elif l.startswith("!Sample_"):
                p = l.rstrip("\n").split("\t"); meta[p[0]].append([x.strip('"') for x in p[1:]])
    X = pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0)
    X.index = X.index.astype(str).str.strip('"'); X.columns = [c.strip('"') for c in X.columns]
    return meta, X


def chars(meta):
    d = {}
    for row in meta["!Sample_characteristics_ch1"]:
        k = row[0].split(":")[0].strip().lower(); d[k] = [x.split(":", 1)[1].strip() for x in row]
    return d


def detect(X, comp, probe2gene, name):
    """median per probe per compartment; percentile of each clock probe among all probes' medians."""
    rows = []
    for c in sorted(comp.unique()):
        med = X.loc[:, comp.index[comp == c]].median(axis=1)
        for g in GENES:
            ps = [p for p, s in probe2gene.items() if g in s and p in med.index]
            if not ps: rows.append(dict(dataset=name, compartment=c, gene=g, probe="absent", median_log2=np.nan, percentile=np.nan)); continue
            for p in ps:
                rows.append(dict(dataset=name, compartment=c, gene=g, probe=p, median_log2=med[p],
                                 percentile=100 * (med < med[p]).mean(), n_probes=len(med)))
    return pd.DataFrame(rows)


def hist2h(t):
    return pd.cut(t, np.arange(0, 26, 2), right=False).value_counts().sort_index()


# ------------------------------------------------------------------ GSE45642
meta, X = read_matrix(R + "brain_tod/GSE45642/GSE45642_series_matrix.txt.gz")
ch = chars(meta)
s = pd.DataFrame(dict(sample_id=meta["!Sample_geo_accession"][0], title=meta["!Sample_title"][0],
                      subject_id=ch["subject id"], region=ch["brain region"], group=ch.get("disease status", ["control"] * X.shape[1]) if False else "control",
                      assay_site=ch["genechip assay site"], original_chip=ch["original genechip type"]))
s["subject_id"] = "S" + s["subject_id"].astype(str)
s["replicate_flag"] = s["title"].str.extract(r"_(\d)_[A-Z]$")[0]
for c in "age sex tod_zt tod_unit tod_type ph pmi rin".split(): s[c] = np.nan
s["platform"] = "GPL17027 (Affymetrix HG-U133A/Plus2 chips, custom CDF Entrez-gene probesets)"
s["subject_key_note"] = "GEO characteristic 'subject id' (integer, shared across regions); age/sex/TOD/pH/PMI/RIN not in GEO (PNAS SI only)"
s.to_csv(OUT + "GSE45642_samples.csv", index=False)
reg = s.region.values
u = s.drop_duplicates(["subject_id", "region"])
ct = pd.crosstab(u.subject_id, u.region); ct.to_csv(OUT + "GSE45642_subject_by_region.csv")
n_all = int((ct.sum(axis=1) == 6).sum())
d = [dict(dataset="GSE45642", item="subjects", value=s.subject_id.nunique()),
     dict(dataset="GSE45642", item="chips", value=len(s)),
     dict(dataset="GSE45642", item="subjects_with_all_6_regions", value=n_all)]
for k, v in ct.sum(axis=1).value_counts().sort_index().items(): d.append(dict(dataset="GSE45642", item=f"subjects_with_{k}_regions", value=v))
for r, v in ct.sum().items(): d.append(dict(dataset="GSE45642", item=f"subjects_{r}", value=v))
for r, v in s.region.value_counts().items(): d.append(dict(dataset="GSE45642", item=f"chips_{r}", value=v))
for k, v in s.assay_site.value_counts().items(): d.append(dict(dataset="GSE45642", item=f"chips_site_{k}", value=v))
pd.crosstab(s.region, s.assay_site).to_csv(OUT + "GSE45642_region_by_assay_site.csv")
# within-subject site mixing: subjects whose regions come from >1 assay site
sites = s.groupby(["subject_id", "region"]).assay_site.agg(lambda x: ",".join(sorted(set(x)))).unstack()
subj_nsites = s.groupby("subject_id").assay_site.nunique()
d.append(dict(dataset="GSE45642", item="subjects_with_chips_from_gt1_assay_site", value=int((subj_nsites > 1).sum())))
chip_mult = s.groupby(["subject_id", "region"]).size()
d.append(dict(dataset="GSE45642", item="subject-region cells with >1 chip", value=int((chip_mult > 1).sum())))
d.append(dict(dataset="GSE45642", item="subject-region cells total", value=len(chip_mult)))
pd.DataFrame(d).to_csv(OUT + "GSE45642_design.csv", index=False)
ann = pd.read_csv(R + "brain_tod/GSE45642/GPL17027_probe_annotation.tsv", sep="\t", dtype=str)
inv = {str(v): g for g, v in ENTREZ.items()}
p2g = {p: {inv[p.split("_")[0]]} for p in X.index if p.split("_")[0] in inv}
det45 = detect(X, pd.Series(reg, index=X.columns), p2g, "GSE45642")
rng = (float(X.values.min()), float(np.nanmedian(X.values)), float(X.values.max()))
info45 = dict(n_probes=X.shape[0], min=rng[0], median=rng[1], max=rng[2])

# ------------------------------------------------------------------ GSE71620
meta, X2 = read_matrix(R + "brain_tod/GSE71620/GSE71620_series_matrix.txt.gz")
ch = chars(meta)
t = pd.DataFrame(dict(sample_id=meta["!Sample_geo_accession"][0], title=meta["!Sample_title"][0],
                      age=pd.to_numeric(ch["age"], errors="coerce"), pmi=pd.to_numeric(ch["pmi"], errors="coerce"),
                      ph=pd.to_numeric(ch["ph"], errors="coerce"), rin=pd.to_numeric(ch["rin"], errors="coerce"),
                      sex=ch["sex"], race=ch["race"], tod_zt=pd.to_numeric(ch["tod"], errors="coerce")))
t["region"] = t.title.str.extract(r"^(BA\d+)")[0]
t["sample_number"] = t.title.str.extract(r"#(\d+)")[0].astype(int)
t["group"] = "control"
# subject key: same demographic/covariate tuple in BA11 and BA47 (as in the lab's earlier derivation)
t["key"] = t[["age", "sex", "race", "pmi", "ph", "tod_zt"]].round(4).apply(lambda r: "|".join(str(x) for x in r), axis=1)
kc = t.groupby("key").region.agg(lambda x: tuple(sorted(x)))
t["subject_id"] = t.key.map({k: f"P{i+1:03d}" for i, k in enumerate(sorted(kc.index))})
t["tod_zt_0_24"] = t.tod_zt % 24
t["tod_unit"] = "hours; ZT, ZT0 = local sunrise on day of death (per Chen 2016; not stated in GEO); negative = before ZT0"
t["platform"] = "GPL11532 (Affymetrix Human Gene 1.1 ST)"
t["subject_key_note"] = "derived: identical (age, sex, race, pmi, pH, TOD) tuple in BA11 and BA47; no donor id in GEO"
t.drop(columns="key").to_csv(OUT + "GSE71620_samples.csv", index=False)
npr = t.groupby("subject_id").region.nunique()
tt = t.drop_duplicates("subject_id")
d = [dict(dataset="GSE71620", item="chips", value=len(t)), dict(dataset="GSE71620", item="subjects_by_tuple", value=t.subject_id.nunique()),
     dict(dataset="GSE71620", item="subjects_with_both_BA11_BA47", value=int((npr == 2).sum())),
     dict(dataset="GSE71620", item="subjects_with_one_region_only", value=int((npr == 1).sum())),
     dict(dataset="GSE71620", item="subjects_with_TOD", value=int(tt.tod_zt.notna().sum())),
     dict(dataset="GSE71620", item="subjects_without_TOD", value=int(tt.tod_zt.isna().sum())),
     dict(dataset="GSE71620", item="sex_F_subjects", value=int((tt.sex == "F").sum())), dict(dataset="GSE71620", item="sex_M_subjects", value=int((tt.sex == "M").sum())),
     dict(dataset="GSE71620", item="age_min/median/max", value=f"{tt.age.min():.0f}/{tt.age.median():.0f}/{tt.age.max():.0f}")]
pd.DataFrame(d).to_csv(OUT + "GSE71620_design.csv", index=False)
bins = pd.cut(tt.age, [0, 29, 39, 49, 59, 69, 79, 100]).value_counts().sort_index()
bins.to_csv(OUT + "GSE71620_age_bins.csv", header=["subjects"])
h = hist2h(tt.tod_zt.dropna() % 24); h.to_csv(OUT + "GSE71620_tod_hist_2h.csv", header=["subjects"])
g = pd.read_csv(io.StringIO("".join(l for l in gzip.open(R + "brain_tod/GSE71620/GPL11532.annot.gz", "rt") if not l.startswith(("#", "^", "!")))), sep="\t", dtype=str)
p2g = {r.ID: {x.strip().upper() for x in str(r["Gene symbol"]).split("///")} for _, r in g.iterrows() if isinstance(r["Gene symbol"], str)}
p2g = {p: {gg for gg in GENES if (ALIAS.get(gg, {gg}) & s_)} for p, s_ in p2g.items()}
p2g = {p: s_ for p, s_ in p2g.items() if s_}
det71 = detect(X2, t.set_index("sample_id").region.reindex(X2.columns), p2g, "GSE71620")
info71 = dict(n_probes=X2.shape[0], min=float(X2.values.min()), median=float(np.nanmedian(X2.values)), max=float(X2.values.max()))

# ------------------------------------------------------------------ skin
rows = []; dets = []; info_sk = {}
plat = pd.read_csv(io.StringIO("".join(l for l in open(R + "skin_wu2020_gse139301/GPL13667_full.txt", errors="replace") if not l.startswith(("!", "^", "#")))), sep="\t", dtype=str, usecols=["ID", "Gene Symbol"])
p2gs = {}
for _, r in plat.iterrows():
    if isinstance(r["Gene Symbol"], str):
        ss = {x.strip().upper() for x in r["Gene Symbol"].split("///")}
        m = {gg for gg in GENES if ALIAS.get(gg, {gg}) & ss}
        if m: p2gs[r.ID] = m
layers = {"GSE139301": "epidermis", "GSE139300": "dermis"}
Xs = {}
for acc, layer in layers.items():
    meta, Xk = read_matrix(R + f"skin_wu2020_gse139301/{acc}_series_matrix.txt.gz"); Xs[acc] = Xk
    ch = chars(meta); title = meta["!Sample_title"][0]; src = meta["!Sample_source_name_ch1"][0]
    df = pd.DataFrame(dict(sample_id=meta["!Sample_geo_accession"][0], title=title, source=src, collection_time_hhmm=ch["collection time"], sex=ch["gender"]))
    df["layer"] = layer; df["series"] = acc
    df["site"] = df.source.str.extract(r"human (\w+)$")[0].replace({"arm": "forearm"})
    df["batch"] = df.title.str.extract(r"^(GSS\d+)_")[0]
    df["subject_num"] = df.title.str.extract(r"^GSS\d+_(\d+)")[0]
    df["arm"] = np.where(df.sex == "male", "timed_men", "untimed_women")
    df["subject_id"] = df.arm.map({"timed_men": "M", "untimed_women": "F"}) + df.subject_num
    df["tod_clock_hhmm"] = pd.to_numeric(df.collection_time_hhmm.replace("NA", np.nan), errors="coerce")
    df["tod_note"] = "GEO 'collection time' hhmm clock time for the timed men (0, 600, 1200, 1800); NA for women"
    df["age"] = np.nan; df["group"] = "healthy"; df["platform"] = "GPL13667 (Affymetrix HG-U219 GeneTitan)"
    rows.append(df)
    info_sk[acc] = dict(n_probes=Xk.shape[0], min=float(Xk.values.min()), median=float(np.nanmedian(Xk.values)), max=float(Xk.values.max()))
sk = pd.concat(rows, ignore_index=True)
sk["compartment"] = sk.layer + "_" + sk.site
sk["subject_key_note"] = "title field 2 (subject number) + sex arm; 'GSS' prefix = hybridisation batch; no age in GEO"
sk.to_csv(OUT + "wu2020_GSE139300_GSE139301_samples.csv", index=False)
for acc in layers:
    sub = sk[sk.series == acc]
    dets.append(detect(Xs[acc], sub.set_index("sample_id").compartment.reindex(Xs[acc].columns), p2gs, acc))
w = sk[sk.arm == "untimed_women"]
ct = pd.crosstab(w.subject_id, w.compartment); ct.to_csv(OUT + "wu2020_women_subject_by_compartment.csv")
cnt = (ct > 0).sum(axis=1)
d = [dict(item="women_subjects", value=w.subject_id.nunique()), dict(item="women_samples", value=len(w)),
     dict(item="women_with_all_6_compartments", value=int((cnt == 6).sum()))]
for k, v in cnt.value_counts().sort_index().items(): d.append(dict(item=f"women_with_{k}_compartments", value=v))
for c, v in (ct > 0).sum().items(): d.append(dict(item=f"women_{c}", value=v))
ep = set(w[w.layer == "epidermis"].subject_id); de = set(w[w.layer == "dermis"].subject_id)
d += [dict(item="women_in_epidermis_series", value=len(ep)), dict(item="women_in_dermis_series", value=len(de)), dict(item="women_in_both_series", value=len(ep & de))]
bw = w.groupby("subject_id").batch.nunique(); d.append(dict(item="women_with_samples_from_gt1_batch", value=int((bw > 1).sum())))
bl = w.groupby(["batch", "layer"]).size().unstack(fill_value=0); bl.to_csv(OUT + "wu2020_women_batch_by_layer.csv")
bsl = w.groupby(["batch", "site"]).size().unstack(fill_value=0); bsl.to_csv(OUT + "wu2020_women_batch_by_site.csv")
dup = w.groupby(["subject_id", "compartment"]).size(); d.append(dict(item="women_subject_compartment_cells_with_gt1_sample", value=int((dup > 1).sum())))
m = sk[sk.arm == "timed_men"]
d += [dict(item="men_subjects", value=m.subject_id.nunique()), dict(item="men_samples", value=len(m)),
      dict(item="men_samples_epidermis", value=int((m.layer == "epidermis").sum())), dict(item="men_samples_dermis", value=int((m.layer == "dermis").sum()))]
for tt_, v in m.collection_time_hhmm.value_counts().sort_index().items(): d.append(dict(item=f"men_collection_time_{tt_}", value=v))
# subject-number overlap between arms
d.append(dict(item="subject_numbers_shared_between_men_and_women", value=len(set(m.subject_num) & set(w.subject_num))))
pd.DataFrame(d).to_csv(OUT + "wu2020_design.csv", index=False)
pd.DataFrame(info_sk).T.to_csv(OUT + "wu2020_matrix_ranges.csv")
pd.DataFrame([info45, info71], index=["GSE45642", "GSE71620"]).to_csv(OUT + "brain_matrix_ranges.csv")
pd.concat([det45, det71] + dets).to_csv(OUT + "brain_skin_clock_probe_detection.csv", index=False)
print("done")
