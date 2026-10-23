"""Single-cell cancer cohorts: one pseudobulk table per cohort (patient x sample type x compartment), and the pooled
table of the cohorts with copy-number calls (malignant vs normal epithelium).

Source:   scBaseCount h5ad per sample + curated inventory, and GEO cohort h5ads (see `cancer_common.studies`), copy-number
          calls `$SCCC_DATA/cnv/<study>.csv.gz` (recomputable with `cancer_common.studies.calls`).
Replaces: legacy/scripts/run_study.py (`load`, `_counts`), legacy/scripts/cnv_calls.py, and the `load` functions of
          legacy/desync/validation/paired_damping/scripts/a86_malignant_gate.py and a87_decomposition.py (pooled table,
          malignant / normal epithelium and stromal host compartments).
Used by:  the cancer pipeline (M4 studies), A83-A87 (tumour gate, malignant gate, decomposition), A97 (delay).
Studies:  43 cohorts (`STUDIES`; 29 with copy-number calls, `CNV_STUDIES`; 24 GEO + 19 scBaseCount).

Choices (as run_study.py; docs/cancer/ANALYSIS_PLAN.md amendments of 2026-10-03 and 2026-10-06)
- donor = patient (all eligible samples of the patient, majority 10x chemistry only; sorted fractions pooled);
  samples: biopsy or resection, not treated, patient known, sample type in tumour / adjacent_normal / healthy / metastasis,
  median UMI >= 1,500 (GEO cohorts also drop CD45+-sorted samples);
- compartment = `<sample_type>:<class>` with classes from marker genes (fibroblast, endothelial, myeloid, T_NK, mural,
  B_plasma, epithelial, mast ...; cells that clearly match none are `unassigned` and dropped). With `cnv=True` the epithelial
  cells of tumour samples become `tumour:malignant` (call malignant) and `tumour:epithelial_nonmalignant` (call normal);
  ambiguous cells and tumour epithelium of cohorts without calls are dropped;
- pseudobulks of fewer than 50 cells are dropped; a compartment is kept if it is present in >= 60 % of the patients or
  in >= 6 patients; stratum = 10x chemistry (empty where the cohort does not record it);
- the eligibility of a cohort for an analysis (>= 6 patients, >= 3 compartments, gate) belongs to the analysis.

Donor table (patient): study, patient, cancer type, sample types present, number of compartments, `has_malignant`.
The pooled table (`pooled`) has donor = `study|patient`, stratum = `study|stratum`, and the columns `study`; its donor
table adds `cancer`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(Path(__file__).resolve().parents[3])]
import clocksync as cs  # noqa: E402
from common.paths import DATA  # noqa: E402
from natexp_common import GENES, clock_adata, compare, standard, sum_tables  # noqa: E402

from cancer_common import studies as S  # noqa: E402

MIN_CELLS, MIN_PRESENT, MIN_PATIENTS = 50, 0.6, 6
MAL, NORM = "tumour:malignant", "adjacent_normal:epithelial"
STROMA = ("fibroblast", "endothelial", "mural")
STUDIES = """ERP130387 ERP149257 ERP153272 ERP156718 GSE131907 GSE132465 GSE144735 GSE149614 GSE159115 GSE166555 GSE178341
GSE181294 GSE181919 GSE183904 GSE184362 GSE200997 GSE203067 GSE212966 GSE222315 GSE236696 GSE270680 GSE271690 GSE281736
GSE294300 GSE315534 GSE326201 GSE335070 SRP273482 SRP301444 SRP309698 SRP327758 SRP329496 SRP335254 SRP353200 SRP407590
SRP414286 SRP416180 SRP429995 SRP467578 SRP514981 SRP539998 SRP567887 SRP572743""".split()
CNV_STUDIES = sorted(p.name.removesuffix(".csv.gz") for p in S.CNV.glob("*.csv.gz"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # paper/
from common.paths import ARCHIVE  # noqa: E402
LEGACY = ARCHIVE / "outputs" / "studies"


def _labels(comp, kind, sample, barcodes, calls):
    lab = pd.Series([f"{kind}:{c}" for c in comp], dtype=object)
    if calls is not None and kind == "tumour":
        c = calls.reindex(pd.MultiIndex.from_arrays([[sample] * len(barcodes), list(barcodes)])).to_numpy()
        name = {"malignant": "tumour:malignant", "normal": "tumour:epithelial_nonmalignant",
                "ambiguous": "tumour:unassigned"}
        e = (comp == "epithelial")
        lab[e] = [name.get(x, "tumour:unassigned") for x in c[e]]
    return lab.to_numpy()


def _filter(pb: pd.DataFrame, compartments: bool = True) -> pd.DataFrame:
    pb = pb[(pb.n_cells >= MIN_CELLS) & ~pb.compartment.str.endswith(":unassigned")]
    if not compartments:
        return pb.reset_index(drop=True)
    present = pb.groupby("compartment").donor.nunique()
    keep = present[(present >= MIN_PRESENT * pb.donor.nunique()) | (present >= MIN_PATIENTS)].index
    return pb[pb.compartment.isin(keep)].reset_index(drop=True)


def donor_table(pb: pd.DataFrame, study: str, cancer: str) -> pd.DataFrame:
    g = pb.groupby("donor")
    don = pd.DataFrame({"study": study, "patient": g.size().index, "cancer": cancer,
                        "sample_types": g.compartment.agg(lambda c: ";".join(sorted({x.split(":")[0] for x in c}))),
                        "n_compartments": g.compartment.nunique(),
                        "has_malignant": g.compartment.agg(lambda c: bool((c == MAL).any()))}, index=g.size().index)
    return don.rename_axis("donor")


def build(study: str, cnv: bool = True, donors: list | None = None, where: str | None = None,
          calls: pd.Series | None = None, compartment_rule: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pseudobulk, donors) of one cohort from the counts. `cnv`: split the tumour epithelium by the stored copy-number
    calls (`calls=` to pass others, e.g. `studies.calls(study)` indexed by (sample, barcode)); `donors`: restrict to
    these patients (checks); `where`: GEO only, 'obs_column=value' subset (sensitivity analyses);
    `compartment_rule=False` skips the cohort-level rule on compartments (meaningless on a slice of patients)."""
    us, cancer = S.units(study, where)
    if donors is not None:
        us = [u for u in us if u[1] in set(donors)]
    us = S.majority_chemistry(us)
    if cnv and calls is None:
        calls = S.read_calls(study)
    parts = []
    for sample, patient, kind, load, chem, _ in us:
        X, _genes, umi, comp, barcodes, sym = load()
        lab = _labels(comp, kind, sample, barcodes, calls if cnv else None)
        obs = pd.DataFrame({"donor": patient, "compartment": lab}, index=pd.Index(barcodes).astype(str))
        t = cs.pseudobulk(clock_adata(X, sym, umi, obs), donor="donor", compartment="compartment", depth="depth",
                          min_cells=1)
        parts.append(t.assign(stratum=chem))
    pb = _filter(sum_tables(parts, keys=("donor", "compartment", "stratum")), compartment_rule)
    return standard(pb), donor_table(pb, study, cancer)


def _read_legacy(study: str, cnv: bool) -> pd.DataFrame:
    f = LEGACY / f"{study}{'_cnv' if cnv else ''}_pseudobulk.csv"
    t = pd.read_csv(f, dtype={"resection": str}, keep_default_na=False)
    t = t.rename(columns={"resection": "donor", "BMAL1": "ARNTL"})
    t["stratum"] = t["stratum"].astype(str) if "stratum" in t else ""
    return t[t.depth > 0]


def from_legacy(study: str | None = None, cnv: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Standard table from the legacy `outputs/studies/<study>[_cnv]_pseudobulk.csv`. `study=None`: the pooled table of
    all cohorts with copy-number calls (`pooled`), as loaded by A86 / A87."""
    if study is None:
        return pooled(cnv=True)
    pb = _read_legacy(study, cnv)
    cancer = S.GEO_CANCER.get(study, "") if study.startswith("GSE") else _inventory_cancer(study)
    return standard(pb), donor_table(pb, study, cancer)


def _inventory_cancer(study: str) -> str:
    inv = pd.read_csv(S.INVENTORY, usecols=["study", "cancer_type"], low_memory=False)
    ct = inv[inv.study == study].cancer_type.dropna()
    return ct.mode().iloc[0] if len(ct) else ""


def pooled(studies=None, cnv: bool = True, source=from_legacy) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pooled table of cohorts (default: all with copy-number calls): donor = `study|patient`, stratum = `study|stratum`,
    column `study`; the input of A86 (gate of malignant vs normal epithelium centred within study x stratum) and A87.
    `source`: `from_legacy` (tables of the legacy run) or `build` (from the counts)."""
    studies = studies or (CNV_STUDIES if cnv else STUDIES)
    pbs, dons = [], []
    for s in studies:
        pb, don = source(s, cnv) if source is from_legacy else source(s, cnv=cnv)
        pb = pb.assign(study=s, donor=s + "|" + pb.donor, stratum=s + "|" + pb.stratum)
        pbs.append(pb)
        dons.append(don.assign(patient=don.index).rename(index=lambda d, s=s: f"{s}|{d}"))
    pb = pd.concat(pbs, ignore_index=True)
    return standard(pb, extra=("study",)), pd.concat(dons).rename_axis("donor")


def paired_epithelium(pb: pd.DataFrame, host=STROMA) -> list:
    """Donors (of a pooled or single table) with malignant epithelium, adjacent-normal epithelium and a stromal host
    compartment: the patients of A87 (alignment against the stromal host) and A97."""
    comp = pb.compartment
    cls = comp.str.split(":").str[1]
    sets = [set(pb.donor[comp == MAL]), set(pb.donor[comp == NORM]), set(pb.donor[cls.isin(host)])]
    return sorted(set.intersection(*sets))


if __name__ == "__main__":
    # checks (see paper/common/datasets/INDEX_natexp_cancer.md): a slice of patients per cohort against the legacy tables
    pooled_pb, pooled_don = from_legacy()
    print("pooled legacy table:", pooled_pb.donor.nunique(), "patients,", pooled_pb.study.nunique(), "studies,",
          pooled_pb.compartment.nunique(), "compartments;", len(paired_epithelium(pooled_pb)), "paired patients", flush=True)
    for study, cnv, n in (("GSE144735", True, 3), ("SRP273482", True, 3), ("GSE203067", True, 2), ("SRP309698", False, 3),
                          ("SRP327758", True, 3), ("GSE270680", True, 3), ("GSE178341", True, 3), ("GSE131907", False, 3)):
        leg, _ = from_legacy(study, cnv)
        sl = sorted(leg.donor.unique())[:n]
        new, don = build(study, cnv=cnv, donors=sl, compartment_rule=False)
        new = new[new.compartment.isin(leg.compartment)]
        print(study, "cnv" if cnv else "no cnv", sl, compare(new.drop(columns="stratum"), leg.drop(columns="stratum"), sl),
              "strata built:", sorted(new.stratum.unique()), flush=True)
    # copy-number calls recomputed for one small cohort against the stored calls
    from cancer_common import studies as S
    new = S.calls("GSE144735")
    old = pd.read_csv(S.CNV / "GSE144735.csv.gz")
    m = new.merge(old, on=["sample", "barcode"], suffixes=("", "_old"))
    print("GSE144735 calls recomputed:", len(new), "cells vs stored", len(old), "| same call:", float((m.cnv_call == m.cnv_call_old).mean()),
          "| max |clonality diff|:", float((m.clonality - m.clonality_old).abs().max()))
