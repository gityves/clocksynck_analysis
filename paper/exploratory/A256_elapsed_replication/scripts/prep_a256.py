"""A256 data preparation: per-sample matrices (HGNC symbols) and sample sheets of the fresh repeat-biopsy studies.

Wrangling only. No clock statistic, no before / after comparison and no expression value of a clock gene is computed or
printed: the script reports shapes, missing shares, whole-matrix medians and counts of persons, units and time points.
The design (arms, class, elapsed hours, clock time of the first biopsy) is fixed in `DESIGN` below from the GEO pages
and the papers' methods (quotes in `design_addendum_draft.md`).

    python -I prep_a256.py            # -> ../output/prep/{<GSE>_expr.parquet, <GSE>_sheet.csv, design.csv, prep_summary.csv}

Conventions
- Unit = person x arm (x visit / leg / bout / session): the first biopsy of the unit has elapsed_h = 0, later ones their
  hours from it. Samples that are not on the same day as the unit's first biopsy are dropped (GSE59088 'Pre',
  GSE202295 'basal').
- Symbols: ids mapped with the platform table, first symbol of `A /// B`, `canonical` aliases (BMAL1 -> ARNTL,
  CHRONO -> CIART); several ids per gene: the one with the highest median; empty symbols dropped.
- Arrays as deposited (log2 if the matrix median > 30); counts as deposited (GSE163356 deposits non-integer values).
- Person-link check (sex genes XIST / RPS4Y1, mixed-sex or unknown-sex cohorts only): per sample the sign of
  RPS4Y1 - XIST against the cohort's midpoint of two means; number of persons whose samples disagree is reported.
"""

from __future__ import annotations

import gzip
import io
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
OUT = HERE.parent / "output" / "prep"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/scouts_2026-11-12/scripts"))
from clocksync.genes import canonical  # noqa: E402
from gsm_meta import read_meta  # noqa: E402

R = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/repeat_biopsy")
CLOCK = "ARNTL CLOCK NPAS2 PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 DBP RORA TEF CIART".split()

# ---- fixed design -------------------------------------------------------------------------------------------------
NOTHING, INTERV = "nothing_done", "intervention"
MORNING = 'assumed (protocol: "{}")'
DESIGN = {
    "GSE40589": dict(tissue="subcutaneous adipose", kind="array", t0=9.0,
                     t0_source=MORNING.format("All studies were done in the morning in the overnight fasted state (>10h since last meal)")),
    "GSE9105": dict(tissue="skeletal muscle", kind="array", t0=9.0,
                    t0_source=MORNING.format("After a 10-h overnight fast, all subjects reported to the General Clinical Research Center (GCRC) at 0800") + "; first biopsy 60 min before insulin start"),
    "GSE59088": dict(tissue="skeletal muscle", kind="array", t0=13.0,
                     t0_source='stated: "After exercise, all subjects rested under fasting conditions ... until 3.30 pm (corresponding to 5 h post exercise)"; first biopsy (2.5 h) = 13:00'),
    "GSE220899": dict(tissue="skeletal muscle", kind="array", t0=9.0,
                      t0_source=MORNING.format("participants completed a resistance exercise bout in the morning hours (between 7:00 AM-10:00 AM) following an overnight fast")),
    "GSE202295": dict(tissue="skeletal muscle", kind="counts", t0=9.0,
                      t0_source=MORNING.format("All participants reported to the laboratory after an overnight fast") + "; no clock time given for the exercise day"),
    "GSE163356": dict(tissue="skeletal muscle", kind="counts", t0=9.0, t0_source="assumed (no time of day given)"),
    "GSE126296": dict(tissue="skeletal muscle", kind="array", t0=9.0,
                      t0_source=MORNING.format("The subjects reported to the laboratory in the morning after an overnight fast")),
}
REJECTED = {  # dataset: (samples, persons, reason); counts from the GEO sample sheets (no expression read)
    "GSE7146": (24, 12, "two separate protocols of 6 persons each (2 h on GPL96, 3 h on GPL80), < 8 persons per protocol and platform; insulin clamp only (no nothing-done arm); person overlap between protocols unknown"),
    "GSE164081": (50, 10, "no per-sample matrix (deposited: TPM per time point, pooled over persons; 3 series on 3 platforms)"),
    "GSE43219": (28, 14, "rest biopsy taken on a separate day (>= 1 week before the exercise day, before an insulin clamp): not a same-day pair"),
}


def say(*a):
    print(*a, flush=True)


# ---- platform tables ----------------------------------------------------------------------------------------------
def plat_table(name: str, annot: bool = False) -> pd.DataFrame:
    """Platform table of a file in repeat_biopsy/platforms (SOFT read from !platform_table_begin to the end of the table)."""
    p = R / "platforms" / name
    lines = []
    with gzip.open(p, "rt", errors="replace") as fh:
        started = not annot and False
        for l in fh:
            if l.startswith("!platform_table_begin"):
                started = True
                continue
            if l.startswith("!platform_table_end"):
                break
            if annot and not started and l.startswith("ID\t"):
                started = True
            if started:
                lines.append(l)
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str, low_memory=False, quoting=3)


def probe_symbols(gpl: str) -> pd.Series:
    if gpl == "GPL96":
        t = plat_table("GPL96.table.txt.gz")
        return t.set_index("ID")["Gene Symbol"]
    if gpl == "GPL6848":
        t = plat_table("GPL6848.table.txt.gz")
        return t[t.CONTROL_TYPE == "FALSE"].set_index("ID")["GENE_SYMBOL"]
    if gpl == "GPL6244":
        return plat_table("GPL6244.annot.gz", annot=True).set_index("ID")["Gene symbol"]
    if gpl == "GPL24324":    # Clariom S HT: the symbol is only in the RefSeq description "... (SYMBOL), transcript variant 1, mRNA."
        t = plat_table("GPL24324.table.txt.gz")
        txt = t.iloc[:, -1].fillna("")
        sym = txt.str.extract(r"\(([A-Za-z0-9._@\-]+)\)(?:, transcript variant \w+)?, (?:mRNA|transcript|non-coding|long|partial|microRNA|antisense|small)")[0]
        return pd.Series(sym.to_numpy(), index=t["ID"])
    raise KeyError(gpl)


def collapse(expr: pd.DataFrame, sym: pd.Series) -> pd.DataFrame:
    """genes x samples with HGNC symbols as index; per symbol the id with the highest median; empty symbols dropped."""
    s = sym.reindex(expr.index)
    ok = s.notna().to_numpy() & (s.astype(str).str.strip() != "").to_numpy()
    e = expr[ok].copy()
    e.index = [canonical(str(x).split("///")[0].strip()) for x in s[ok]]
    e = e.iloc[np.argsort(-e.median(axis=1).to_numpy(), kind="stable")]
    return e[~e.index.duplicated()]


def ens_map() -> pd.Series:
    t = pd.read_csv(R / "GSE231509/GSE231509_Read_counts.txt.gz", sep="\t", usecols=["gene_id", "gene_name"])
    return pd.Series(t.gene_name.to_numpy(), index=t.gene_id.str.split(".").str[0])


def series(acc: str):
    meta, e = read_meta(str(R / acc / f"{acc}_series_matrix.txt.gz"))
    meta = meta.set_index("geo_accession")
    if e is not None:
        e.index = e.index.astype(str).str.strip('"')
        e = e.apply(pd.to_numeric, errors="coerce")
    return meta, e


# ---- sheets -------------------------------------------------------------------------------------------------------
def rec(sample, person, unit, arm, cls, interv, elapsed):
    return dict(sample=sample, person=person, unit=unit, arm=arm, arm_class=cls, intervention=interv, elapsed_h=elapsed)


def sheet_40589(meta):
    drug = {"GIP": "GIP infusion 2 pmol/kg/min for 240 min", "NACL": "saline"}
    clamp = {"INF": "", "EU": "euglycemic hyperinsulinemic clamp (80 mg/dl)", "HC": "hyperglycemic hyperinsulinemic clamp (140 mg/dl) with diazoxide"}
    rows = []
    for g, r in meta.iterrows():
        m = re.match(r"SCAT_(\d+)_(\w+)_(\w+)_(\d+)", r["title"])
        t, d, c = int(m[1]), m[2], m[3]
        person = f"P{int(re.search(r'replicate(\d+)', r['description'])[1])}"
        parts = [x for x in (clamp[c], drug[d] if d == "GIP" else ("saline infusion" if c == "INF" else "saline")) if x]
        cls = NOTHING if (d == "NACL" and c == "INF") else INTERV
        rows.append(rec(g, person, f"{person}|{c}_{d}", f"{c}_{d}", cls, "" if cls == NOTHING else " + ".join(parts), t / 60))
    return rows


def sheet_9105(meta):
    el = {"0min": 0.0, "30min": 1.5, "240min": 5.0}     # biopsies at -60, +30, +240 min of the insulin infusion (paper)
    rows = []
    for g, r in meta.iterrows():
        m = re.match(r"(S\d+)_Muscle_\w+_(\d+min)", r["title"])
        rows.append(rec(g, m[1], m[1], "insulin_clamp", INTERV, "euglycemic hyperinsulinemic clamp 80 mU/m2/min (4 h)", el[m[2]]))
    return rows


def sheet_59088(meta):
    arms = {"C": (NOTHING, ""), "E": (INTERV, "2 h cycling at 60% VO2peak (08:30-10:30) before the first biopsy; fasting rest after"),
            "S": (INTERV, "30 min resistance exercise (ends 10:30) before the first biopsy; fasting rest after")}
    rows = []
    for g, r in meta.iterrows():
        m = re.match(r"(\d+)([ABC])_([ESC])_", r["title"])
        if m[2] == "A":        # 'Pre': biopsy before the training phase (>= 10 weeks earlier), not the same day
            continue
        cls, iv = arms[m[3]]
        rows.append(rec(g, f"P{m[1]}", f"P{m[1]}", {"C": "control", "E": "endurance", "S": "resistance"}[m[3]], cls, iv, 0.0 if m[2] == "B" else 2.5))
    return rows


def sheet_220899(meta):
    rows = []
    for g, r in meta.iterrows():
        m = re.search(r"_(X\d+)_(30|80)_(PRE|3H|6H)", r["title"])
        iv = f"resistance exercise to failure, 4 sets squat + 4 sets leg extension at {m[2]}% 1RM"
        rows.append(rec(g, m[1], f"{m[1]}|{m[2]}", f"{m[2]}%1RM", INTERV, iv, {"PRE": 0.0, "3H": 4.0, "6H": 7.0}[m[3]]))   # 1 h bout assumed
    return rows


def sheet_202295(meta):
    rows = []
    for g, r in meta.iterrows():
        m = re.match(r"(NGT|T2D)_(\d+)_(basal|post|rec)", r["title"])
        if m[3] == "basal":    # rest biopsy of a baseline visit 5-10 days before the exercise day
            continue
        rows.append(rec(g, f"{m[1]}_{m[2]}", f"{m[1]}_{m[2]}", m[1], INTERV,
                        "30 min cycling at 85% max heart rate before the first (post) biopsy; rest after", 0.0 if m[3] == "post" else 3.0))
    return rows


def sheet_163356(meta):
    rows = []
    for g, r in meta.iterrows():
        s, tcol = r["c_subject"], re.search(r"(T[19])_([12])$", r["title"])
        rows.append(rec(g, s, f"{s}|{tcol[1]}", f"session{tcol[1][1:]}", INTERV, "HIIT bout 10 x 4 min at 91% HRmax",
                        0.0 if tcol[2] == "1" else 4.0))      # 3 h after the bout; bout of about 1 h assumed
    return rows


def sheet_126296(meta):
    rows = []
    for g, r in meta.iterrows():
        p = r["c_individual"]
        rows.append(rec(g, p, p, "sprint", INTERV, "3 x 30 s all-out cycle sprints, 20 min rest between", 0.0 if "pre" in r["title"] else 220 / 60))
    return rows


SHEETS = {"GSE40589": sheet_40589, "GSE9105": sheet_9105, "GSE59088": sheet_59088, "GSE220899": sheet_220899,
          "GSE202295": sheet_202295, "GSE163356": sheet_163356, "GSE126296": sheet_126296}
PLATFORM = {"GSE40589": "GPL6848", "GSE9105": "GPL96", "GSE59088": "GPL6244", "GSE220899": "GPL24324", "GSE126296": "GPL6244"}


def counts_matrix(acc, meta):
    """genes x GSM (symbols) for the two count tables."""
    if acc == "GSE202295":
        e = pd.read_csv(R / acc / "GSE202295_gene_counts.txt.gz", sep=" ", quotechar='"')
        sym = pd.Series(e["gene_name"].to_numpy(), index=e.index)
        e = e.drop(columns=["Geneid", "gene_name"]) if "Geneid" in e.columns else e.drop(columns=["gene_name"])
        e.columns = [re.match(r"(P7401_\d+)", c)[1] for c in e.columns]
        by = {re.search(r"\[(P7401_\d+)\]", t)[1]: g for g, t in meta["title"].items()}
        e = e.rename(columns=by)
        return collapse(e, sym)
    e = pd.read_csv(R / acc / "GSE163356_counts.tsv.gz", sep="\t", index_col=0)
    by = {re.search(r"(\w_T[19]_[12])$", t)[1]: g for g, t in meta["title"].items()}
    e = e.rename(columns=by)
    e.index = e.index.str.split(".").str[0]
    return collapse(e[~e.index.duplicated()], ens_map())


def sex_check(e: pd.DataFrame, sheet: pd.DataFrame) -> str:
    """Persons whose samples disagree on sex. Score = mean of the Y genes present (RPS4Y1 KDM5D DDX3Y EIF1AY UTY) minus
    XIST when present; two-means split; 'uninformative' when the two means are < 4 within-person SDs apart."""
    y = [g for g in ("RPS4Y1", "KDM5D", "DDX3Y", "EIF1AY", "UTY") if g in e.index]
    if not y:
        return "sex genes absent"
    ids = sheet["sample"]
    d = e.loc[y, ids].astype(float).mean(axis=0) - (e.loc["XIST", ids].astype(float) if "XIST" in e.index else 0.0)
    lo, hi = d.min(), d.max()
    for _ in range(50):
        thr = (lo + hi) / 2
        lo, hi = d[d <= thr].mean(), d[d > thr].mean()
    within = float(np.sqrt((d.groupby(sheet["person"].to_numpy()).transform(lambda v: v - v.mean()) ** 2).mean()))
    if hi - lo < 4 * within:
        return "uninformative (no sex separation: single-sex cohort or ratio-scale data)"
    call = (d > (lo + hi) / 2).to_numpy()
    s = sheet.assign(male=call).groupby("person").male.agg(["min", "max"])
    return f"{int((s['min'] != s['max']).sum())}/{len(s)} persons discordant (called male {int(call.sum())}/{len(call)} samples)"


# ---- main ---------------------------------------------------------------------------------------------------------
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summ, design = [], []
    for acc, D in DESIGN.items():
        say(acc)
        meta, e = series(acc)
        sh = pd.DataFrame(SHEETS[acc](meta))
        if D["kind"] == "counts":
            e = counts_matrix(acc, meta)
            med, miss = float(np.nanmedian(e.to_numpy())), float(np.isnan(e.to_numpy()).mean())
            deposited_median = med
            note = "counts as deposited" + ("; non-integer values (normalised / expected counts)" if not np.allclose(e.iloc[:200].to_numpy(), np.round(e.iloc[:200].to_numpy())) else "")
            e = e[[c for c in sh["sample"] if c in e.columns]]
        else:
            e = e[[c for c in sh["sample"] if c in e.columns]]
            deposited_median = float(np.nanmedian(e.to_numpy()))
            miss = float(np.isnan(e.to_numpy()).mean())
            e = collapse(e, probe_symbols(PLATFORM[acc]))
            note = f"array as deposited ({PLATFORM[acc]}); median {deposited_median:.2f}"
            if deposited_median > 30:
                e = np.log2(e.clip(lower=1.0))
                note += "; linear scale, log2-transformed"
        sh = sh[sh["sample"].isin(e.columns)]
        # units: one first biopsy and >= 1 later, same day
        sh = sh[sh.elapsed_h <= 12]
        ok = sh.groupby("unit").elapsed_h.transform(lambda v: (v == 0).sum() == 1 and (v > 0).sum() >= 1 and v.nunique() == len(v))
        sh = sh[ok].copy()
        sh["tissue"], sh["kind"], sh["scale_note"] = D["tissue"], D["kind"], note
        sh["t0_clock_h"], sh["t0_source"] = D["t0"], D["t0_source"]
        sh = sh[["sample", "person", "unit", "arm", "arm_class", "intervention", "elapsed_h", "t0_clock_h", "t0_source", "tissue", "kind", "scale_note"]]
        sh.to_csv(OUT / f"{acc}_sheet.csv", index=False)
        ex = e[sh["sample"]].T.astype(np.float32)
        ex.index.name = "sample"
        ex.to_parquet(OUT / f"{acc}_expr.parquet")
        chk = "male-only (metadata)" if acc in ("GSE59088", "GSE220899", "GSE202295", "GSE163356") else sex_check(e, sh)
        later = sh[sh.elapsed_h > 0]
        for (arm, el), g in later.groupby(["arm", "elapsed_h"]):
            a = sh[sh.arm == arm].iloc[0]
            design.append(dict(dataset=acc, tissue=D["tissue"], arm=arm, arm_class=a.arm_class, intervention=a.intervention, elapsed_h=el,
                               t0_clock_h=D["t0"], t0_source=D["t0_source"], units=g.unit.nunique(), persons=g.person.nunique()))
        npres = sum(g in ex.columns for g in CLOCK)
        nunit, npers = sh.unit.nunique(), sh.person.nunique()
        keep = npers >= 8
        summ.append(dict(dataset=acc, samples=len(sh), persons=npers, units=nunit, genes=ex.shape[1], kind=D["kind"],
                         matrix_median=round(deposited_median, 3), missing_share=round(miss, 6), clock_genes_present=npres,
                         person_link_check=chk, retained="yes" if keep and npres >= 12 else "no",
                         why_not="" if keep and npres >= 12 else "< 8 persons or clock genes missing"))
        say(f"  samples {len(sh)} persons {npers} units {nunit} genes {ex.shape[1]} clock {npres}/14  {chk}")
    for acc, (n, p, why) in REJECTED.items():
        summ.append(dict(dataset=acc, samples=n, persons=p, retained="no", why_not=why))
    pd.DataFrame(summ).to_csv(OUT / "prep_summary.csv", index=False)
    pd.DataFrame(design).to_csv(OUT / "design.csv", index=False)
    say(pd.DataFrame(design).drop(columns=["t0_source", "intervention"]).to_string())


if __name__ == "__main__":
    main()
