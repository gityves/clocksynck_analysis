"""NHS breast (GSE93601, GSE115577): sample metadata, deposit overlap, unique tumour + adjacent patients, ER, and clock/G2M gene
coverage on the platform annotations. Metadata and annotation only; no expression value is read or used.
Run: python -I meta_nhs.py   (parses downloaded text only)"""
import re
from pathlib import Path
import pandas as pd

SCCC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
OUT = Path(__file__).resolve().parents[1] / "output"
CLOCK = ["ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "TEF", "HLF", "CIART"]
import sys
sys.path.insert(0, "/home/maxine/Documents/paychere/scCircadianOverview/paper/exploratory/A168_endocrine/scripts")
ALIAS = {"ARNTL": ["BMAL1"], "CIART": ["CHRONO", "C1orf51"], "NR1D1": ["REV-ERBA"], "NR1D2": ["REV-ERBB"]}


def _g2m():
    for line in open(SCCC / "reference/h.all.v2024.1.Hs.symbols.gmt").read().splitlines():
        f = line.split("\t")
        if f[0] == "HALLMARK_G2M_CHECKPOINT":
            return {g for g in f[2:] if g not in CLOCK}


HALLMARK_G2M = _g2m()


def parse(path):
    rows, cur = [], None
    for l in open(path, encoding="utf-8", errors="replace"):
        l = l.rstrip("\n")
        if l.startswith("^SAMPLE"):
            cur = {"gsm": l.split("=")[1].strip()}
            rows.append(cur)
        elif cur is not None and l.startswith("!Sample_characteristics_ch1"):
            k, _, v = l.split("=", 1)[1].partition(":")
            cur[k.strip()] = v.strip()
        elif cur is not None and l.startswith("!Sample_title"):
            cur["title"] = l.split("=", 1)[1].strip()
        elif cur is not None and l.startswith("!Sample_platform_id"):
            cur["platform"] = l.split("=")[1].strip()
        elif cur is not None and l.startswith("!Sample_relation"):
            cur["relation"] = l.split("=", 1)[1].strip()
    return pd.DataFrame(rows)


def symbols(path, col_candidates):
    """Gene symbols of the platform table."""
    syms, hdr, on = set(), None, False
    for l in open(path, encoding="utf-8", errors="replace"):
        if l.startswith("!platform_table_begin"):
            on = True
            continue
        if l.startswith("!platform_table_end"):
            break
        if on:
            f = l.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                idx = [hdr.index(c) for c in col_candidates if c in hdr][0]
                continue
            if len(f) > idx:
                syms.update(re.split(r"\s*///\s*|\s*;\s*|\s*,\s*", f[idx].strip()))
    syms.discard("---")
    return syms


def main():
    a = parse(SCCC / "nhs_gse93601/GSE93601_gsm_brief.txt")
    b = parse(SCCC / "nhs_gse115577/GSE115577_gsm_brief.txt")
    b["tissue2"] = b["tumor tumoradj"].map({"tumor": "tumour", "tumor-adjacent normal": "normal"})
    a["tissue2"] = a["tissue"].map({"breast tumor": "tumour", "breast tumor-adjacent normal": "normal"})
    # deposit overlap: GSE115577 'Reanalysis of: GSM...' points at the GSE93601 sample (same CEL file)
    r = b.relation.dropna().str.extract(r"(GSM\d+)")[0].dropna()
    rows = [dict(item="GSE93601 samples", value=len(a)), dict(item="GSE115577 samples", value=len(b)),
            dict(item="GSE115577 samples on GPL22920 (hGlue 3.0, = GSE93601 CELs)", value=int((b.platform == "GPL22920").sum())),
            dict(item="GSE115577 samples on GPL17585 (HTA 2.0, new)", value=int((b.platform == "GPL17585").sum())),
            dict(item="GSE115577 'Reanalysis of' a GSE93601 GSM (in GSE93601)", value=int(r.isin(set(a.gsm)).sum())),
            dict(item="GSE93601 patients", value=a["patient id"].nunique()),
            dict(item="GSE115577 patients (sample id)", value=b["sample id"].nunique()),
            dict(item="patients in both deposits", value=len(set(a["patient id"]) & set(b["sample id"])))]
    # pairs per platform: patients (sample id) with exactly one tumour and one adjacent on the same platform
    pairs = []
    for plat, d in b.groupby("platform"):
        cnt = d.groupby(["sample id", "tissue2"]).size().unstack(fill_value=0)
        for k in ("tumour", "normal"):
            if k not in cnt:
                cnt[k] = 0
        one = cnt[(cnt.tumour == 1) & (cnt.normal == 1)].index
        anyp = cnt[(cnt.tumour > 0) & (cnt.normal > 0)].index
        rows += [dict(item=f"{plat}: patients with >= 1 tumour and >= 1 adjacent", value=len(anyp)),
                 dict(item=f"{plat}: patients with exactly 1 + 1", value=len(one))]
        pairs.append(d[d["sample id"].isin(one)].assign(pair_platform=plat))
    pr = pd.concat(pairs)
    allp = b.groupby("sample id").tissue2.agg(lambda x: set(x))
    rows.append(dict(item="unique patients with tumour + adjacent in the union (any platform)", value=int(sum(s == {"tumour", "normal"} for s in allp))))
    both_plat = b.groupby("sample id").platform.nunique()
    rows.append(dict(item="patients with samples on both platforms", value=int((both_plat > 1).sum())))
    t = pr[pr.tissue2 == "tumour"].copy()
    for col, lab in (("er ihc", "ER (IHC; 1 = positive)"), ("molecularsubtype pam50 modifiedmedianmethod", "PAM50 (modified median)"), ("ki67 ihc", "Ki67 IHC")):
        for plat, d in t.groupby("pair_platform"):
            for k, v in d[col].value_counts(dropna=False).items():
                rows.append(dict(item=f"{plat} paired tumours, {lab}: {k}", value=int(v)))
    rows.append(dict(item="proliferation score available (paired tumours, non-missing)",
                     value=int(pd.to_numeric(t["proliferationscore modifiedmedianmethod"], errors="coerce").notna().sum())))
    rows.append(dict(item="metadata fields (GSE115577)", value=";".join(b.columns)))
    cov = []
    for gpl, cols in (("GPL22920", ["Symbol"]), ("GPL17585", ["gene_symbols"])):
        s = symbols(SCCC / f"nhs_platforms/{gpl}_full.txt", cols)
        pres = [g for g in CLOCK if g in s or any(x in s for x in ALIAS.get(g, []))]
        g2m = sorted(HALLMARK_G2M & s)
        cov.append(dict(platform=gpl, n_symbols=len(s), n_clock=len(pres), present=";".join(pres), missing=";".join(g for g in CLOCK if g not in pres),
                        g2m_covered=len(g2m), g2m_total=len(HALLMARK_G2M)))
    pd.DataFrame(rows).to_csv(OUT / "nhs_counts.csv", index=False)
    pd.DataFrame(cov).to_csv(OUT / "nhs_clock_coverage.csv", index=False)
    pr.drop(columns=["title", "relation"]).to_csv(OUT / "nhs_pairs_samples.csv", index=False)  # sample-level metadata only
    print(pd.DataFrame(rows).to_string()); print(pd.DataFrame(cov).to_string())


main()
