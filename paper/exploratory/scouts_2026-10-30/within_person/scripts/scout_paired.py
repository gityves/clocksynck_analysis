"""Metadata-only scout for within-person disease designs (roadmap 6.6).
Reads sample metadata and cell counts only (no expression values). Run: python -I scripts/scout_paired.py
Writes CSVs to ../ (this folder) and prints summary tables."""
import os, sys, gzip, glob
import numpy as np, pandas as pd, anndata as ad

D = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian"
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MIN = 50
summ = []

def add(name, fam, n_pat, n_samples, note=""):
    summ.append(dict(candidate=name, family=fam, n_paired_patients=n_pat, n_samples=n_samples, note=note))

def pair_table(df, pat, cond, fams, a, b, name):
    """df: one row per (patient, cond, [site/visit]) unit with columns per family n_cells."""
    ok = {f: (df[f] >= MIN) for f in fams}
    df = df.assign(**{f"ok_{f}": ok[f] for f in fams})
    df["n_elig"] = sum(ok[f].astype(int) for f in fams)
    res = {}
    for f in fams:
        pa = set(df.loc[df[cond] == a, pat][ok[f][df[cond] == a]])
        pb = set(df.loc[df[cond] == b, pat][ok[f][df[cond] == b]])
        res[f] = len(pa & pb)
    both = set(df.loc[df[cond] == a, pat]) & set(df.loc[df[cond] == b, pat])
    return df, res, both

# ---------- 1. TAURUS (Thomas 2024) ----------
T = f"{D}/ibd/TAURUS/pseudobulk"
s = pd.read_csv(f"{T}/samples.csv"); cc = pd.read_csv(f"{T}/class_counts.csv").set_index("sample_id")
fibc = ["LP_fibroblast", "Epi_fibroblast", "C3pos_fibroblast", "Fibroblast", "THY1pos_FAPpos_PDPNpos_fibroblast"]
cc2 = pd.DataFrame({"fibroblast": cc[fibc].sum(1), "myofibroblast": cc["Myofibroblast"], "endothelial": cc["Endothelium"],
                    "mural": cc["Pericyte"], "epithelial": cc[["Non_ileal_epithelium", "Ileal_epithelium"]].sum(1),
                    "T": cc[["CD4_T", "CD8_T", "Unconventional_T"]].sum(1), "B_plasma": cc[["B", "Plasma"]].sum(1),
                    "myeloid": cc[["Mono_macro", "DC", "Cycling_MNP", "Mast"]].sum(1)})
t = s.merge(cc2, left_on="sample", right_index=True, how="left")
t.to_csv(f"{OUT}/taurus_sample_class_cells.csv", index=False)
print("TAURUS samples", len(t), "patients", t.patient.nunique()); print(t.disease.value_counts().to_dict(), t.inflammation.value_counts().to_dict())
print(t.groupby("patient").inflammation.nunique().value_counts().to_dict(), "patients by n inflammation levels")
print("timepoint", t.timepoint.value_counts().to_dict(), "Library", t.LibraryType.value_counts().to_dict(), "site", t.site.value_counts().to_dict())
fams = ["fibroblast", "endothelial"]
# unit 1: sample
d, res, both = pair_table(t, "patient", "inflammation", fams, "Inflamed", "Non_Inflamed", "taurus")
print("TAURUS patients with both:", len(both), "per-family (>=50 cells in both states, any sample):", res)
# unit 2: patient x condition pooled over sites
pc = t.groupby(["patient", "inflammation"])[["fibroblast", "myofibroblast", "endothelial", "mural", "epithelial", "T", "B_plasma", "myeloid"]].sum().reset_index()
d2, res2, both2 = pair_table(pc, "patient", "inflammation", fams, "Inflamed", "Non_Inflamed", "taurus_pooled")
both_both = len(set(d2[(d2.inflammation=="Inflamed")&d2.ok_fibroblast&d2.ok_endothelial].patient) & set(d2[(d2.inflammation=="Non_Inflamed")&d2.ok_fibroblast&d2.ok_endothelial].patient))
print("pooled over sites:", res2, "both families both states:", both_both)
pc.to_csv(f"{OUT}/taurus_patient_condition_cells.csv", index=False)
# same-site pairs (inflamed and non-inflamed at the same site, rare)
ss = t.groupby(["patient", "site"]).inflammation.nunique(); print("patient-sites with both states", int((ss == 2).sum()))
elig = t.assign(n_elig=(t.fibroblast >= MIN).astype(int) + (t.endothelial >= MIN).astype(int), n_elig3=(t.fibroblast >= MIN).astype(int) + (t.endothelial >= MIN) + (t.mural >= MIN))
print("samples with n_elig (fib,endo) =", elig.n_elig.value_counts().sort_index().to_dict(), "by inflammation:", elig.groupby("inflammation").n_elig.value_counts().unstack(fill_value=0).to_dict("index"))
# median cells
for f in ["fibroblast", "endothelial", "mural", "epithelial", "T", "myeloid"]:
    print(f, "median cells/sample by inflammation", t.groupby("inflammation")[f].median().to_dict())
print("pretreatment", t.pretreatment.value_counts().to_dict(), "Remission", t.Remission_status.value_counts().to_dict())
# paired with >= 2 eligible families at pooled level, by disease
pat_ok = set(d2[(d2.inflammation=="Inflamed")&d2.ok_fibroblast&d2.ok_endothelial].patient)&set(d2[(d2.inflammation=="Non_Inflamed")&d2.ok_fibroblast&d2.ok_endothelial].patient)
print("pooled both-families patients by disease", t[t.patient.isin(pat_ok)].drop_duplicates("patient").disease.value_counts().to_dict())
print("sample-level: patients with >=1 inflamed sample and >=1 non-inflamed sample where each has fib & endo >=50:",
      len(set(elig[(elig.inflammation=='Inflamed')&(elig.n_elig==2)].patient)&set(elig[(elig.inflammation=='Non_Inflamed')&(elig.n_elig==2)].patient)))
add("TAURUS (Thomas 2024)", "fib / endo (pooled sites)", f"{res2['fibroblast']} / {res2['endothelial']}; both {both_both}", len(t))

# ---------- 2. Psoriasis sc (GSE173706, GSE228421, GSE278330) ----------
p = pd.read_csv(f"{D}/psoriasis_sc/pseudobulk.csv", usecols=["study", "patient", "visit", "status", "sample", "compartment", "n_cells"])
w = p.pivot_table(index=["study", "patient", "visit", "status"], columns="compartment", values="n_cells", aggfunc="sum", fill_value=0).reset_index()
w.to_csv(f"{OUT}/psoriasis_sample_compartment_cells.csv", index=False)
w["pat"] = w.study + "_" + w.patient.astype(str)
print("\nPSORIASIS units", len(w), "per study patients", w.groupby("study").pat.nunique().to_dict())
for scope, sub in [("all visits pooled", w.groupby(["study", "pat", "status"])[[c for c in ["fibroblast","endothelial","mural","keratinocyte","T","myeloid","melanocyte"]]].sum().reset_index()),
                   ("same visit", w.assign(pat=w.pat + "_" + w.visit.astype(str)))]:
    sub = sub.copy()
    if "pat" not in sub: pass
    dd, r, both = pair_table(sub, "pat", "status", fams, "lesional", "nonlesional", "pso")
    bf = set(dd[(dd.status=="lesional")&dd.ok_fibroblast&dd.ok_endothelial].pat) & set(dd[(dd.status=="nonlesional")&dd.ok_fibroblast&dd.ok_endothelial].pat)
    st = dd.drop_duplicates("pat").set_index("pat").study if "study" in dd else None
    print(scope, "patients with both:", len(both), r, "both fams:", len(bf), "by study", pd.Series([x.split("_")[0] for x in bf]).value_counts().to_dict())
    if scope == "same visit": add("Psoriasis sc (3 GEO series), same visit", "fib / endo", f"{r['fibroblast']} / {r['endothelial']}; both {len(bf)}", len(dd))
    else: add("Psoriasis sc (3 GEO series), visits pooled", "fib / endo", f"{r['fibroblast']} / {r['endothelial']}; both {len(bf)}", len(dd))
print("median cells fib/endo/mural by status", w.groupby("status")[["fibroblast","endothelial","mural"]].median().to_dict())

# ---------- 3. Kong atlases (CELLxGENE Census pseudobulk, Crohn vs normal tissue, same donor) ----------
print()
kong = []
for k, label in [("kong_ti_str", "Kong TI stroma"), ("kong_co_str", "Kong colon stroma"), ("ileum_kong", "Kong ileum immune"), ("kong_ti_epi","Kong TI epithelium"), ("kong_co_epi","Kong colon epithelium"), ("kong_co_imm","Kong colon immune")]:
    f = f"{D}/atlas/{k}/pseudobulk.h5ad"
    if not os.path.exists(f): continue
    o = ad.read_h5ad(f, backed="r").obs.reset_index(drop=True)
    o["cls"] = o.cell_type.astype(str)
    o["dis"] = o.disease.astype(str).map(lambda x: "CD" if "Crohn" in x else "normal")
    g = o.groupby(["donor_id", "dis", "cls"]).n_cells.sum().unstack("cls", fill_value=0).reset_index()
    donors_both = set(g[g.dis == "CD"].donor_id) & set(g[g.dis == "normal"].donor_id)
    fcols = [c for c in ["fibroblast", "endothelial cell"] if c in g]
    row = dict(dataset=label, donors=g.donor_id.nunique(), donors_with_CD_and_normal=len(donors_both), tissues=",".join(sorted(o.tissue.astype(str).unique())))
    for c in g.columns[2:]:
        if c in ("fibroblast", "endothelial cell", "myofibroblast cell", "pericyte"):
            a = set(g[(g.dis=="CD")&(g[c]>=MIN)].donor_id); b = set(g[(g.dis=="normal")&(g[c]>=MIN)].donor_id)
            row[f"paired_{c}>=50"] = len(a & b & donors_both)
    if len(fcols) == 2:
        a = set(g[(g.dis=="CD")&(g[fcols[0]]>=MIN)&(g[fcols[1]]>=MIN)].donor_id); b = set(g[(g.dis=="normal")&(g[fcols[0]]>=MIN)&(g[fcols[1]]>=MIN)].donor_id)
        row["paired_both_families"] = len(a & b)
    kong.append(row); g.to_csv(f"{OUT}/{k}_donor_disease_class_cells.csv", index=False)
kg = pd.DataFrame(kong); kg.to_csv(f"{OUT}/kong_summary.csv", index=False); print(kg.T.to_string())
r = kg[kg.dataset == "Kong TI stroma"]
if len(r): add("Kong 2025 TI stroma (Census)", "fib / endo", f"{int(r['paired_fibroblast>=50'].iloc[0])} / {int(r['paired_endothelial cell>=50'].iloc[0])}; both {int(r.paired_both_families.iloc[0])}", "")
r = kg[kg.dataset == "Kong colon stroma"]
if len(r): add("Kong 2025 colon stroma (Census)", "fib / endo", f"{int(r['paired_fibroblast>=50'].iloc[0])} / {int(r['paired_endothelial cell>=50'].iloc[0])}; both {int(r.paired_both_families.iloc[0])}", "")

# ---------- 4. GSE266546 (annotation only on disk) ----------
print()
parts = []
for f in glob.glob(f"{D}/inflammation_scout/GSE266546/*.txt.gz"):
    x = pd.read_csv(f, sep="\t", usecols=["SampleID", "PatientID", "Tissue", "Histology", "Source", "Celltypes"]); x["file"] = os.path.basename(f); parts.append(x)
x = pd.concat(parts)
fib = x.Celltypes.str.contains("Fibro"); endo = x.Celltypes.str.startswith("Endo")
x["cls"] = np.where(fib, "fibroblast", np.where(endo, "endothelial", "other"))
g = x[x.cls != "other"].groupby(["PatientID", "Tissue", "Histology", "Source", "SampleID", "cls"]).size().unstack("cls", fill_value=0).reset_index()
g.to_csv(f"{OUT}/gse266546_sample_cells.csv", index=False)
print("GSE266546 cells", len(x), "samples", x.SampleID.nunique(), "patients", x.PatientID.nunique()); print(x.drop_duplicates("SampleID").groupby(["Tissue", "Source", "Histology"]).size())
pp = x.drop_duplicates("SampleID").groupby("PatientID").Histology.nunique(); print("patients with >1 histology category:", int((pp > 1).sum()))
pp2 = x.drop_duplicates("SampleID").groupby(["PatientID", "Tissue"]).Histology.nunique(); print("patient x tissue with >1 histology:", int((pp2 > 1).sum()))
print("fibroblast/endothelial cells per sample: median", g[["fibroblast", "endothelial"]].median().to_dict(), "samples with fib>=50", int((g.fibroblast >= MIN).sum()), "endo>=50", int((g.endothelial >= MIN).sum()), "of", x.SampleID.nunique())
mix = x.drop_duplicates("SampleID").assign(inflamed=lambda d: d.Histology.map(lambda h: "inflamed" if h in ("Severe","Moderate","Mild","Inflamed") else ("normal" if h in ("Control","Normal_CD","Normal") else "other")))
print(mix.Histology.value_counts().to_dict())

# ---------- 5. IBDverse (E-MTAB-16986/16999) ----------
print()
for e in ["E-MTAB-16986", "E-MTAB-16999"]:
    im = pd.read_csv(f"{D}/inflammation_scout/{e}/individual_metadata.tsv", sep="\t"); sm = pd.read_csv(f"{D}/inflammation_scout/{e}/sample_metadata.tsv", sep="\t")
    print(e, "individuals", len(im), im.disease_status.value_counts().to_dict(), im.ti_ses_inflamed.value_counts().to_dict(), "samples", len(sm), sm.biopsy_type.value_counts().to_dict())
    if os.path.exists(f"{D}/inflammation_scout/{e}/sample_individual_mapping.tsv"):
        mp = pd.read_csv(f"{D}/inflammation_scout/{e}/sample_individual_mapping.tsv", sep="\t"); print(" samples per individual:", mp.groupby("individual_id").size().value_counts().to_dict())
pd.DataFrame(summ).to_csv(f"{OUT}/summary_raw.csv", index=False)

# ---------- addenda: TAURUS design structure, psoriasis per study ----------
print("\nTAURUS timepoint x inflammation\n", pd.crosstab(t.timepoint, t.inflammation), "\nMatch\n", t.Match.value_counts().to_dict())
pp = t[t.inflammation.isin(["Inflamed","Non_Inflamed"])].groupby("patient").apply(lambda d: pd.Series(dict(n_pre=(d.timepoint=="Pre").sum(), n_post=(d.timepoint=="Post").sum(), sites=d.site.nunique())), include_groups=False)
print("TAURUS patient structure: patients with both Pre and Post", int(((pp.n_pre>0)&(pp.n_post>0)).sum()))
both_same_tp = [pt for pt, d in t[t.inflammation.isin(["Inflamed","Non_Inflamed"])].groupby("patient") if any(dd.inflammation.nunique()==2 for _, dd in d.groupby("timepoint"))]
print("patients with both states at the same timepoint:", len(both_same_tp))
print("TAURUS healthy donors:", t[t.disease=="Healthy"].patient.nunique())
print("\nPSORIASIS per study fibroblast >=50 samples by status\n", w.assign(f=w.fibroblast>=MIN, e=w.endothelial>=MIN).groupby(["study","status"])[["f","e"]].sum())
print(w.groupby(["study","status"]).size())
