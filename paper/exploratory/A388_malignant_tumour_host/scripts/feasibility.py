"""A388 feasibility inventory (counts only; no clock statistic is computed).

Per cohort of the A149 / A226 / A227 / H138.1 lists: patients with a tumour malignant (or tumour-epithelium) pseudobulk
and how many of them have >= 2 / >= 3 tumour stromal compartments among fibroblast, endothelial, mural (cell / UMI /
presence rules of the existing table builders: >= 50 cells per pseudobulk, sample median UMI >= 1,500, compartment kept
if in >= 60 % of patients or >= 6). Writes output/feasibility/cohorts.csv (inventory.csv is assembled by `inventory.py`
step in this file's `main`, with the breast atlases added by `breast_atlas.py`).
"""
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
STR = ("fibroblast", "endothelial", "mural")
H138 = set(pd.read_csv(PAPER / "exploratory/A138_cancer_estimability/output/score/per_unit.csv").query("hypothesis=='H138.1'").unit)
A149_ESTIMABLE = set(pd.read_csv(M.A149_SCORE).query("kind=='s1' and estimable").cohort)
K.install()


def a83(c):
    return "keeping" if c in L.KEEPING else "at null" if c in L.AT_NULL else "other"


def one(study):
    src = [n for n, s in (("A149", set(L.cohorts())), ("A226", set(K.COHORTS)),
                          ("A227", set(M.cohort_list()))) if study in s]
    if study in H138:
        src.append("H138.1")
    if study in K.COHORTS:
        pb, cancer = K.table(study), K.COHORTS[study]
        mal = L.MAL if L.MAL in set(pb.compartment) else L.EPI
    else:
        pb, mal = L.table(study)
        cancer = M.cancer(study)
    m = pb[pb.compartment == mal]
    pts = set(m.donor)
    n_anchor = {d: 0 for d in pts}
    per = {}
    for k in STR:
        dk = set(pb[pb.compartment == f"tumour:{k}"].donor) & pts
        per[k] = len(dk)
        for d in dk:
            n_anchor[d] += 1
    na = pd.Series(n_anchor)
    return dict(study=study, cancer=cancer, a83_class=a83(cancer), source="+".join(src), malignant_label=mal,
                chemistry_strata=";".join(f"{k}:{v}" for k, v in m.drop_duplicates("donor").stratum.astype(str).value_counts().items()),
                n_patients_total=pb.donor.nunique(), n_malignant=len(pts), n_fibroblast=per["fibroblast"],
                n_endothelial=per["endothelial"], n_mural=per["mural"], n_ge2_anchors=int((na >= 2).sum()),
                n_ge3_anchors=int((na >= 3).sum()), in_H138_1=study in H138, error="")


def main():
    studies = list(dict.fromkeys(L.cohorts() + list(K.COHORTS)))
    rows = []
    for s in studies:
        try:
            rows.append(one(s))
        except Exception as e:  # noqa: BLE001
            rows.append(dict(study=s, cancer=M.cancer(s), error=repr(e)))
            print("FAIL", s, repr(e), flush=True)
    pd.DataFrame(rows).to_csv(OUT / "cohorts.csv", index=False)
    print(len(rows), "cohorts")


if __name__ == "__main__" and "--assemble" not in sys.argv:
    main()


def assemble():
    """cohorts.csv + breast_patients_private_ids.csv (breast_atlas.py) -> inventory.csv, per_cancer.csv, breast_counts.csv."""
    c = pd.read_csv(OUT / "cohorts.csv")
    b = pd.read_csv(OUT / "breast_patients_private_ids.csv")
    rows = []
    for coh, study, lab in (("Wu", "GSE176078", "authors' copy-number-based Cancer Epithelial"),
                            ("Pal", "GSE161529", "authors' tumour clusters (inferCNV-annotated)")):
        d = b[b.cohort == coh]
        rows.append(dict(study=study, cancer="breast", a83_class="keeping", source="breast atlas (A170 / A187)",
                         malignant_label=lab, chemistry_strata="not read per sample (10x)",
                         n_patients_total=len(d), n_malignant=int(d.has_malignant.sum()),
                         n_fibroblast=int(d.has_fibroblast.sum()), n_endothelial=int(d.has_endothelial.sum()),
                         n_mural=int(d.has_mural.sum()), n_ge2_anchors=int(d.mal_ge2.sum()), n_ge3_anchors=int(d.mal_ge3.sum()),
                         in_H138_1=False, error="stromal classes: authors' CAF / Endothelial / PVL" if coh == "Wu"
                         else "stromal classes by cancer_common marker rule on non-tumour-cluster cells (authors give none)"))
    inv = pd.concat([c, pd.DataFrame(rows)], ignore_index=True)
    inv["new_vs_H138_1"] = ~inv.in_H138_1.fillna(False).astype(bool)
    inv.to_csv(OUT / "inventory.csv", index=False)
    ok = inv[inv.error.isna() | inv.error.astype(str).str.startswith("stromal")]
    g = ok.groupby(["cancer", "a83_class"]).agg(cohorts=("study", "size"), cohorts_with_ge2=("n_ge2_anchors", lambda s: int((s > 0).sum())),
                                                  patients_malignant=("n_malignant", "sum"), patients_ge2=("n_ge2_anchors", "sum"),
                                                  patients_ge3=("n_ge3_anchors", "sum")).reset_index().sort_values("patients_ge2", ascending=False)
    g.to_csv(OUT / "per_cancer.csv", index=False)
    rows = []
    for coh, d in b.groupby("cohort"):
        for var in ("er", "subtype"):
            for k, x in d.groupby(var):
                rows.append(dict(cohort=coh, variable=var, level=k, n_patients=len(x), n_malignant=int(x.has_malignant.sum()),
                                 n_mal_ge2=int(x.mal_ge2.sum()), n_mal_ge3=int(x.mal_ge3.sum())))
    pd.DataFrame(rows).to_csv(OUT / "breast_counts.csv", index=False)


if __name__ == "__main__" and "--assemble" in sys.argv:
    assemble()
