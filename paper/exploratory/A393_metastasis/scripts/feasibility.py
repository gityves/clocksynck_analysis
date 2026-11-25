"""A393 feasibility (counts only; no clock statistic). Task 2: metastasis samples with single-cell stroma. Per cohort of the
existing builders (A149 / A226 tables; >= 50 cells per pseudobulk, sample median UMI >= 1,500, compartment kept if in >= 60 %
of patients or >= 6): patients with a metastasis malignant / epithelial pseudobulk and >= 2 metastasis stromal pseudobulks
(fibroblast, endothelial, mural), and those that also have a primary-tumour or adjacent-normal sample. Also reads the obs of the
Moorman 2025 CELLxGENE object (counts of cells by sample type and cell class; no expression read).
Writes output/feasibility/counts.csv and moorman_counts.csv."""
import sys, warnings
from pathlib import Path
import pandas as pd

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
for p in ("A149_paired_tumour_normal", "A226_stroma_pairing_confirm", "A227_tme_compartment_pairing"):
    sys.path.insert(0, str(PAPER / "exploratory" / p / "scripts"))
import a149_lib as L  # noqa: E402
import compute_a226 as K  # noqa: E402
import compute_a227 as M  # noqa: E402

OUT = HERE.parent / "output" / "feasibility"
ANCH = ("fibroblast", "endothelial", "mural")
MOOR = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/cellxgene_829a3cd1_moorman2025/829a3cd1.h5ad")
K.install()
A227 = set(M.cohort_list())


def one(study):
    pb = K.table(study) if study in K.COHORTS else L.table(study)[0]
    comp = pb.compartment
    mets = pb[comp.str.startswith("metastasis:")]
    if mets.empty:
        return None
    S = lambda c: set(pb[comp == c].donor)  # noqa: E731
    mal = S("metastasis:epithelial") | S("metastasis:malignant")
    anc = {k: S(f"metastasis:{k}") for k in ANCH}
    n_anc = {d: sum(d in anc[k] for k in ANCH) for d in set(mets.donor)}
    ok = {d for d in mal if n_anc.get(d, 0) >= 2}
    prim = set(pb[comp.str.startswith("tumour:")].donor)
    norm = set(pb[comp.str.startswith("adjacent_normal:")].donor)
    return dict(study=study, in_A227=study in A227, patients_with_metastasis=len(set(mets.donor)),
                met_malignant_or_epithelial=len(mal), met_fib=len(anc["fibroblast"]), met_endo=len(anc["endothelial"]),
                met_mural=len(anc["mural"]), met_mal_and_ge2_anchors=len(ok), of_which_with_primary_tumour=len(ok & prim),
                of_which_with_adjacent_normal=len(ok & norm), of_which_with_either=len(ok & (prim | norm)))


def moorman():
    import anndata as ad
    o = ad.read_h5ad(MOOR, backed="r").obs[["Patient", "Sample Type", "Site", "Treatment", "cell_type"]].astype(str).fillna("none")
    g = o.groupby(["Patient", "Sample Type", "Site", "Treatment"], dropna=False).size().rename("n_cells").reset_index()
    c = o.groupby(["Sample Type", "cell_type"]).size().rename("n_cells").reset_index()
    return g, c


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in list(dict.fromkeys(list(L.cohorts()) + list(K.COHORTS) + sorted(A227))):
        try:
            r = one(s)
        except Exception as e:  # noqa: BLE001
            r = dict(study=s, in_A227=s in A227, error=repr(e)[:100])
        if r:
            rows.append(r)
    pd.DataFrame(rows).to_csv(OUT / "counts.csv", index=False)
    g, c = moorman()
    g.to_csv(OUT / "moorman_patient_sample_counts.csv", index=False)
    c.to_csv(OUT / "moorman_celltype_counts.csv", index=False)


if __name__ == "__main__":
    main()
