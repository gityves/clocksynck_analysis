"""NHS breast tumour / adjacent-normal pairs (GSE93601, GPL22920): data preparation for the A175 H175.1 replication, blind to every
pairing statistic. Only expression, the A97 G2M score of the paired tumours, its tertiles (rank within cohort) and ER (IHC) are built.

Pairs: patients (`patient id`) with exactly one tumour and one tumour-adjacent normal sample (445). Values: deposited RMA log2 (APT),
probes averaged per symbol (A171 `expression_file`; platform table `Symbol` column, '---' dropped). Reuses A177 `finish` (G2M, tertiles,
parquet) and the A171 parser by import.
The 467 HTA 2.0 (GPL17585) samples of GSE115577 have no expression in the downloaded files (the 1,110-sample csv are re-analyses of the
GSE93601 CELs; the GPL17585 series matrix has an empty table; values are only in RAW.tar): they are out (`hta_check` in prep_log.csv).
Run: python -I prep_nhs.py   (parses downloaded text only; nothing downloaded is executed)"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "A177_c1_replication/scripts"))
import prep_a177 as B  # noqa: E402  (also puts A168 / A171 scripts on sys.path)

Q = B.Q
OUT = HERE.parent / "output"
SCCC = B.SCCC
B.OUT = OUT  # `finish` writes here, not into A177/output
KEY = "nhs_gpl22920"


def probe_symbols(gpl):
    """Platform table of the full GEO SOFT (`nhs_platforms/GPL22920_full.txt`): column `Symbol`, '---' = no gene."""
    with open(SCCC / "nhs_platforms" / f"{gpl}_full.txt", errors="replace") as fh:
        lines = [l for l in fh if not l.startswith(("^", "!", "#"))]
    a = pd.read_csv(pd.io.common.StringIO("".join(lines)), sep="\t", dtype=str, usecols=["ID", "Symbol"]).dropna()
    a = a[a.Symbol.str.strip() != "---"]
    a["sym"] = a.Symbol.str.split("///")
    a = a.explode("sym")
    a["sym"] = a.sym.str.strip().map(Q.P.canon)
    a = a[a.sym != ""]
    Q.P.LOG += [(gpl, "probes_with_symbol", a.ID.nunique()), (gpl, "probes_with_several_symbols", int((a.groupby("ID").size() > 1).sum()))]
    return a.set_index("ID")["sym"]


Q.probe_symbols = probe_symbols
Q.P.D = SCCC


def sample_meta(path):
    """GSM -> tissue, patient id, ER status (GSE93601 `estrogen receptor status`)."""
    rows, cur = [], None
    for l in open(path, encoding="utf-8", errors="replace"):
        l = l.rstrip("\n")
        if l.startswith("^SAMPLE"):
            cur = {"gsm": l.split("=")[1].strip()}
            rows.append(cur)
        elif cur is not None and l.startswith("!Sample_characteristics_ch1"):
            k, _, v = l.split("=", 1)[1].partition(":")
            cur[k.strip()] = v.strip()
    return pd.DataFrame(rows)


def er_ihc_from_115577():
    """ER (IHC) of GSE115577 (same patients, `sample id` = `patient id`): per patient, 1 positive / 0 negative / NaN."""
    rows, cur = [], None
    for l in open(SCCC / "nhs_gse115577/GSE115577_gsm_brief.txt", encoding="utf-8", errors="replace"):
        l = l.rstrip("\n")
        if l.startswith("^SAMPLE"):
            cur = {}
            rows.append(cur)
        elif cur is not None and l.startswith("!Sample_characteristics_ch1"):
            k, _, v = l.split("=", 1)[1].partition(":")
            cur[k.strip()] = v.strip()
    d = pd.DataFrame(rows)
    d["er"] = pd.to_numeric(d["er ihc"], errors="coerce")
    return d.groupby("sample id").er.agg(lambda x: x.dropna().iloc[0] if x.notna().any() else np.nan)


def hta_check():
    """Does any downloaded file hold expression for the GPL17585 samples? Series-matrix table row count and csv column overlap."""
    gpl = [x for x in open(SCCC / "nhs_gse115577/GSE115577_gsm_brief.txt", errors="replace") if x.startswith("!Sample_platform_id")]
    n_hta = sum("GPL17585" in x for x in gpl)
    with gzip.open(SCCC / "nhs_gse115577/GSE115577-GPL17585_series_matrix.txt.gz", "rt") as fh:
        t = [l for l in fh if not l.startswith("!")]
    rows_in_table = len(t) - 1  # minus the ID_REF header
    with gzip.open(SCCC / "nhs_gse115577/GSE115577_data_1110_25334probes_GEO.csv.gz", "rt") as fh:
        ncol = len(fh.readline().split(",")) - 1
    return [("nhs", "hta_check_GPL17585_samples", n_hta), ("nhs", "hta_check_GPL17585_series_matrix_probe_rows", rows_in_table),
            ("nhs", "hta_check_1110_csv_columns", ncol), ("nhs", "hta_check_csv_vs_GPL22920_samples", "csv columns = GPL22920 re-analyses only")]


def main():
    m = sample_meta(SCCC / "nhs_gse93601/GSE93601_gsm_brief.txt")
    m["tis"] = m.tissue.map({"breast tumor": "tumour", "breast tumor-adjacent normal": "normal"})
    assert m.tis.notna().all()
    cnt = m.groupby(["patient id", "tis"]).size().unstack(fill_value=0)
    both = sorted(cnt[(cnt.tumour == 1) & (cnt.normal == 1)].index)
    t = m[m.tis == "tumour"].set_index("patient id").loc[both]
    n = m[m.tis == "normal"].set_index("patient id").loc[both]
    ex = Q.expression_file("nhs_gse93601", "GSE93601_series_matrix.txt.gz", "GPL22920")
    mx = [r for r in Q.P.LOG if r[0] == "nhs_gse93601" and r[1] == "deposited_max"][0][2]
    assert mx < 30, mx  # RMA log2 (a linear scale would exceed 100 and be converted by expression_file)
    ex = ex[ex.isna().mean(axis=1) <= 0.2]
    X = ex.T.astype(np.float32)
    X = X.fillna(X.mean())
    B.finish(KEY, both, X.loc[t.gsm.to_numpy()].reset_index(drop=True), X.loc[n.gsm.to_numpy()].reset_index(drop=True))
    # ER: GSE93601 `estrogen receptor status` (primary) and GSE115577 `er ihc` (cross-check), on the paired tumours
    er_a = t["estrogen receptor status"].map({"Positive": 1, "Negative": 0}).to_numpy(float)
    er_b = pd.Series(both).map(er_ihc_from_115577()).to_numpy(float)
    both_ok = np.isfinite(er_a) & np.isfinite(er_b)
    p = pd.read_csv(OUT / f"{KEY}_pairs.csv")
    p["er"] = np.where(np.isfinite(er_a), er_a, -1).astype(int)  # -1 = missing (excluded from the ER contrast)
    p["er_ihc_115577"] = np.where(np.isfinite(er_b), er_b, -1).astype(int)
    p.to_csv(OUT / f"{KEY}_pairs.csv", index=False)
    log = [(KEY, "er_positive", int((p.er == 1).sum())), (KEY, "er_negative", int((p.er == 0).sum())), (KEY, "er_missing", int((p.er < 0).sum())),
           (KEY, "er_source_agreement_n", int(both_ok.sum())), (KEY, "er_source_disagreements", int((er_a[both_ok] != er_b[both_ok]).sum())),
           (KEY, "tumour_tertile_sizes", ";".join(str(int((p.tertile == k).sum())) for k in range(3))),
           (KEY, "n_patients_gse93601", int(m["patient id"].nunique())), (KEY, "n_patients_exactly_1_plus_1", len(both))]
    pd.DataFrame(B.LOG + Q.P.LOG + log + hta_check(), columns=["series", "item", "value"]).to_csv(OUT / "prep_log.csv", index=False)


if __name__ == "__main__":
    main()
