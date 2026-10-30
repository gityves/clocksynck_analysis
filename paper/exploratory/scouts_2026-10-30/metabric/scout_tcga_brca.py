"""TCGA-BRCA feasibility scout: ER (IHC, recount3 XML metadata) and PAM50 (PanCanAtlas via cBioPortal). Metadata only. Run with python -I."""
import pandas as pd
R = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/tcga/recount3/BRCA/"
S = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/tcga/pancan_subtypes/brca_tcga_pan_can_atlas_2018/"
m = pd.read_csv(R+"tcga.tcga.BRCA.MD.gz", sep="\t", low_memory=False)
m["patient"] = m["gdc_cases.submitter_id"]
print("recount3 BRCA samples", len(m), m.cgc_sample_sample_type.value_counts().to_dict())
print("histology (xml_histological_type) among primaries:")
pr = m[m.cgc_sample_sample_type == "Primary Tumor"].copy()
print(pr["xml_histological_type"].value_counts(dropna=False).to_dict() if "xml_histological_type" in pr else "NA col")
print("primary rows", len(pr), "patients", pr.patient.nunique(), "patients with >1 primary aliquot", int(pr.patient.duplicated().sum()))
p1 = pr.sort_values("tcga_barcode").drop_duplicates("patient").copy()  # one primary per patient (first by barcode; scout only)
st = pd.read_csv(S+"data_clinical_patient.txt", sep="\t", comment="#")[["PATIENT_ID","SUBTYPE"]]
p1 = p1.merge(st, left_on="patient", right_on="PATIENT_ID", how="left")
p1["ER"] = p1["xml_breast_carcinoma_estrogen_receptor_status"].fillna("NA")
print("\nER (IHC) one primary per patient:", p1.ER.value_counts().to_dict())
print("PAM50 (PanCanAtlas SUBTYPE) coverage:", p1.SUBTYPE.fillna("NA").value_counts().to_dict())
print("\nPAM50 x ER IHC"); print(pd.crosstab(p1.SUBTYPE.fillna("NA"), p1.ER, margins=True))
print("\nER detection method:", p1["xml_er_detection_method_text"].value_counts(dropna=False).to_dict())
