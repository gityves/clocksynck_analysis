"""A392 feasibility (counts only; no clock statistic). Task 1: tumour + patient-matched adjacent-normal single-cell cohorts
for a replication of A387. Per cohort of the existing builders (A149 / A226 tables, A149 rules: >= 50 cells per pseudobulk,
sample median UMI >= 1,500, compartment kept if in >= 60 % of patients or >= 6): patients with >= 2 normal anchors
(fibroblast / endothelial / mural) and a tumour fibroblast or endothelial pseudobulk, and whether the cohort is in the A227 list.
Writes output/feasibility/counts.csv (candidates.csv and report.md are written by hand from it and from the file inventory)."""
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
K.install()
A227 = set(M.cohort_list())


def one(study):
    if study in K.COHORTS:
        pb, cancer = K.table(study), K.COHORTS[study]
    else:
        pb, cancer = L.table(study)[0], M.cancer(study)
    sets = {(t, k): set(pb[pb.compartment == f"{t}:{k}"].donor) for t in ("tumour", "adjacent_normal") for k in ANCH}
    tum_any = set().union(*[sets["tumour", k] for k in ANCH]) | set(pb[pb.compartment.str.startswith("tumour:")].donor)
    nor_any = set(pb[pb.compartment.str.startswith("adjacent_normal:")].donor)
    both = tum_any & nor_any
    n_norm = {d: sum(d in sets["adjacent_normal", k] for k in ANCH) for d in both}
    tf = sets["tumour", "fibroblast"] | sets["tumour", "endothelial"]
    ok = [d for d in both if n_norm[d] >= 2 and d in tf]
    ok_fib = [d for d in both if n_norm[d] >= 2 and d in sets["tumour", "fibroblast"]]
    ok_end = [d for d in both if n_norm[d] >= 2 and d in sets["tumour", "endothelial"]]
    return dict(study=study, cancer=cancer, in_A227=study in A227, patients_total=pb.donor.nunique(),
                patients_tumour_any=len(tum_any), patients_normal_any=len(nor_any), patients_tumour_and_normal=len(both),
                both_with_ge2_normal_anchors=sum(v >= 2 for v in n_norm.values()),
                n_norm_fib=len(sets["adjacent_normal", "fibroblast"]), n_norm_endo=len(sets["adjacent_normal", "endothelial"]),
                n_norm_mural=len(sets["adjacent_normal", "mural"]),
                n_tum_fib=len(sets["tumour", "fibroblast"]), n_tum_endo=len(sets["tumour", "endothelial"]),
                usable_ge2norm_and_tumfib_or_endo=len(ok), usable_with_tum_fib=len(ok_fib), usable_with_tum_endo=len(ok_end),
                error="")


def main():
    studies = list(dict.fromkeys(list(L.cohorts()) + list(K.COHORTS) + sorted(A227)))
    rows = []
    for s in studies:
        try:
            rows.append(one(s))
        except Exception as e:  # noqa: BLE001
            rows.append(dict(study=s, in_A227=s in A227, error=repr(e)[:120]))
        print(s, flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / "counts.csv", index=False)


if __name__ == "__main__":
    main()
