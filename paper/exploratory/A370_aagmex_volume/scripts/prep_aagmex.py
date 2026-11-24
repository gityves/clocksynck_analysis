"""A370 prep: AAGMEx (GSE95674 adipose, GSE95675 muscle) -> per-tissue matrices (persons x genes, log2) and a sheet.

Parsing only (downloads are untrusted): `python -I prep_aagmex.py`. Rows are gene symbols (aliases mapped as A254);
deposited values are quantile-normalised on the linear scale, log2-transformed (A254 `to_log2`). Person = adipose title;
muscle title minus "-M" / "-RM" (a "-RM" title replaces a plain "-M" one of the same person if both exist).
Checks: tissue markers; a sex call (Y-linked genes minus XIST, split at the largest gap within each tissue) concordant
across a person's two tissues. No clock statistic is computed here.
Outputs (`../output/prep/`): aagmex_expr.parquet, aagmex_sheet.csv, markers.csv, prep_summary.csv.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.normpath(os.path.join(HERE, "../.."))
sys.path.insert(0, os.path.join(EXP, "scouts_2026-11-12/scripts"))
sys.path.insert(0, os.path.join(EXP, "A254_paired_lesions/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scoutlib import ALIAS, CLOCK, ROOT, read_matrix, read_sheet  # noqa: E402
import prep_lesions_tissues as A254  # noqa: E402

D = f"{ROOT}/living_multitissue/AAGMEx_GSE95676"
OUT = os.path.normpath(os.path.join(HERE, "../output/prep"))
SERIES = {"adipose": "GSE95674", "muscle": "GSE95675"}


def main():
    os.makedirs(OUT, exist_ok=True)
    sheets, mats = [], []
    for tis, g in SERIES.items():
        f = f"{D}/{g}_series_matrix.txt.gz"
        s = read_sheet(f)
        sh = pd.DataFrame(dict(sample=s.geo.astype(str), title=s.title.astype(str), tissue=tis))
        sh["person"] = sh.title.map(lambda t: re.sub(r"-R?M$", "", t))
        sh["rerun"] = sh.title.str.endswith("-RM")
        sh = sh.sort_values("rerun").drop_duplicates("person", keep="last")
        m = read_matrix(f)
        m.index = m.index.astype(str)
        m = m[list(sh["sample"])].apply(pd.to_numeric, errors="coerce").astype("float32")
        sym = pd.Series(m.index.map(lambda x: ALIAS.get(x, x)), index=m.index)
        mats.append(A254.collapse(m, sym))
        sheets.append(sh)
    sh = pd.concat(sheets, ignore_index=True)
    common = sorted(set(mats[0].columns) & set(mats[1].columns))
    e, note = A254.to_log2(pd.concat([x[common] for x in mats]))
    e = e.loc[sh["sample"]]
    both = set(sh[sh.tissue == "adipose"].person) & set(sh[sh.tissue == "muscle"].person)
    e.astype("float32").to_parquet(f"{OUT}/aagmex_expr.parquet")
    sh.to_csv(f"{OUT}/aagmex_sheet.csv", index=False)
    t = sh.set_index("sample").tissue
    markers = e[[g for g in ("ADIPOQ", "LEP", "PLIN1", "MYH2", "ACTA1", "CKM", "MYH7") if g in e]].groupby(t.values).mean()
    markers.round(2).to_csv(f"{OUT}/markers.csv")
    print(markers.round(1).to_string())
    y = [g for g in ("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY") if g in e.columns]
    calls = []
    for tis, i in t.groupby(t.values).groups.items():
        sc = (e.loc[i, y].mean(axis=1) - e.loc[i, "XIST"]).sort_values()
        k = int(np.argmax(np.diff(sc.to_numpy())))
        calls.append(pd.Series(np.where(np.arange(len(sc)) > k, "M", "F"), index=sc.index))
    call = pd.concat(calls)
    pers = sh.set_index("sample").person
    pb = call.groupby(pers.loc[call.index].to_numpy()).agg(list)
    pb = pb[pb.map(len) == 2]
    clock = [g for g in CLOCK if g in e.columns]
    summ = dict(samples=len(e), persons_adipose=int((sh.tissue == "adipose").sum()), persons_muscle=int((sh.tissue == "muscle").sum()),
                persons_both=len(both), reruns_used=int(sh.rerun.sum()), genes=e.shape[1], scale_note=note,
                matrix_median=float(np.nanmedian(e.values)), missing_share=float(e.isna().values.mean()),
                clock_genes_present=len(clock), clock_missing=",".join(sorted(set(CLOCK) - set(clock))),
                sex_y_genes=",".join(y), persons_sex_concordant=int(pb.map(lambda v: v[0] == v[1]).sum()),
                persons_sex_checked=len(pb), share_male_calls=float((call == "M").mean()))
    pd.Series(summ).to_csv(f"{OUT}/prep_summary.csv", header=["value"])
    print(pd.Series(summ).to_string())


if __name__ == "__main__":
    main()
