"""A369 prep: STAGE (GSE40231) RMA values -> per-tissue matrices (persons x genes, log2) and a sample sheet.

Values: `rma_stage.R` on the CEL files (the deposited series-matrix values carry wrong probe identities: tissue markers
do not separate the tissues, checked 2026-11-24); the series matrix gives the sample sheet only. Parsing only
(downloads are untrusted): run with `python -I prep_stage.py`. Probe -> gene symbol of the GEO annotation (a clock gene
if any symbol of the probe is one, e.g. PER1 listed after MIR6883; otherwise the first), best-median probe per symbol
(A254 `collapse`).
Checks: tissue markers (ALB liver, MYH2 / ACTA1 muscle, ADIPOQ fat, MYH11 arteries: mean log2 per tissue); identity
by a sex call (Y-linked genes minus XIST, split at the largest gap within each tissue), concordant across a person's
tissues.
No clock statistic is computed here.
Outputs (`../output/prep/`): stage_expr.parquet, stage_sheet.csv, prep_summary.csv.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.normpath(os.path.join(HERE, "../.."))
sys.path.insert(0, os.path.join(EXP, "scouts_2026-11-12/scripts"))
sys.path.insert(0, os.path.join(EXP, "A254_paired_lesions/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scoutlib import ROOT, read_annot, read_sheet  # noqa: E402
import prep_lesions_tissues as A254  # noqa: E402

D = f"{ROOT}/living_multitissue/GSE40231"
OUT = os.path.normpath(os.path.join(HERE, "../output/prep"))
TISSUE = {"Visceral fat": "visceral_fat", "Liver": "liver", "Skeletal muscle": "muscle",
          "Atherosclerotic aortic wall": "aortic_wall", "Internal mammary artery": "mammary_artery"}


def main():
    os.makedirs(OUT, exist_ok=True)
    s = read_sheet(f"{D}/GSE40231_series_matrix.txt.gz")
    sh = pd.DataFrame(dict(sample=s.geo.astype(str), title=s.title.astype(str), tissue=s.tissue.map(TISSUE)))
    sh["person"] = sh.title.str.extract(r"Patient ID\s+(\d+)")[0]
    assert sh.person.notna().all() and sh.tissue.notna().all()
    m = pd.read_csv(f"{D}/GSE40231_rma_probesets.tsv.gz", sep="\t", index_col=0)
    m.index = m.index.astype(str)
    m = m[list(sh["sample"])].astype("float32")
    ann = read_annot(f"{D}/GPL570.annot.gz")
    clock = set(A254.CLOCK)
    sym = ann.apply(lambda L: next((g for g in L if g in clock), L[0]) if len(L) else np.nan)
    sym.index = sym.index.astype(str)
    e, note = A254.collapse(m, sym), "RMA log2 (rma_stage.R, affy 1.88)"
    e = e.loc[sh["sample"]]
    e.astype("float32").to_parquet(f"{OUT}/stage_expr.parquet")
    sh.to_csv(f"{OUT}/stage_sheet.csv", index=False)
    t = sh.set_index("sample").tissue
    markers = e[[g for g in ("ALB", "APOA1", "MYH2", "ACTA1", "ADIPOQ", "LEP", "MYH11", "CNN1") if g in e]].groupby(t.values).mean()
    markers.round(2).to_csv(f"{OUT}/markers.csv")
    print(markers.round(1).to_string())
    y = [g for g in ("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY") if g in e.columns]
    calls = []
    for tis, i in t.groupby(t.values).groups.items():
        sc = (e.loc[i, y].mean(axis=1) - e.loc[i, "XIST"]).sort_values()
        k = int(np.argmax(np.diff(sc.to_numpy())))
        calls.append(pd.Series(np.where(np.arange(len(sc)) > k, "M", "F"), index=sc.index))
    call = pd.concat(calls)
    per = call.groupby(t.index.map(sh.set_index("sample").person).to_series(index=t.index).loc[call.index]).nunique()
    sex = dict(persons_sex_concordant=int((per == 1).sum()), persons_checked=len(per), share_male_calls=float((call == "M").mean()))
    clock = [g for g in A254.CLOCK if g in e.columns]
    summ = dict(samples=len(e), persons=sh.person.nunique(), genes=e.shape[1], scale_note=note,
                matrix_median=float(np.nanmedian(e.values)), missing_share=float(e.isna().values.mean()),
                clock_genes_present=len(clock), clock_missing=",".join(sorted(set(A254.CLOCK) - set(clock))),
                **{f"n_{t}": int((sh.tissue == t).sum()) for t in TISSUE.values()},
                persons_3plus=int((sh.groupby("person").tissue.nunique() >= 3).sum()),
                persons_5=int((sh.groupby("person").tissue.nunique() == 5).sum()), **sex)
    pd.Series(summ).to_csv(f"{OUT}/prep_summary.csv", header=["value"])
    print(pd.Series(summ).to_string())


if __name__ == "__main__":
    main()
