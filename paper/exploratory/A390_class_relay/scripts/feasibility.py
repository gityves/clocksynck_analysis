"""A390 feasibility (counts only; no clock statistic): healthy atlases with genome-wide builders not read by A17 / A389,
donors with endothelial, fibroblast-stromal and mural pseudobulks (classes by `clocksync.anchors.family`), per atlas and
site (compartment prefix), at >= 20 and >= 50 cells per pseudobulk. Disease atlases: control donors only when a
condition column exists. Output: output/feasibility/{compartments,counts}.csv.

    python feasibility.py [--jobs 2]
"""
import argparse
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import importlib.util  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync.anchors import family  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output" / "feasibility"
ATLASES = ["breast_hbca_genome", "skin_ganier_genome", "hlca_fibrotic_genome", "lung_natri_genome", "colon_chrono_genome"]
USED = {"breast_hbca_genome": "organ_breast_genome", "hlca_fibrotic_genome": "lung_hlca_genome",
        "lung_natri_genome": "lung_natri_control_genome"}  # donor overlap with atlases A389 used
CLASSES = {"endothelial": "endothelial", "fibroblast_stromal": "fibroblast", "mural": "mural"}
CONTROL = {"normal", "control", "healthy"}


def run(name):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    spec = importlib.util.spec_from_file_location("grs", PAPER / "results/figure_2/scripts/gene_rhythms/_shared.py")
    grs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(grs)
    try:
        pb, donors = grs.genome_table(name)
    except Exception as e:  # noqa: BLE001
        return dict(atlas=name, error=f"{type(e).__name__}: {e}"[:200]), None
    note = ""
    cond = [c for c in donors.columns if any(k in c.lower() for k in ("disease", "condition", "status", "diagnosis"))]
    if cond:
        c = cond[0]
        ok = donors[c].astype(str).str.lower().str.strip().isin(CONTROL)
        note = f"controls by {c}: {int(ok.sum())} of {len(donors)} ({sorted(donors[c].astype(str).unique())[:6]})"
        pb = pb[pb.donor.isin(donors.index[ok])]
    comp = pb[["compartment"]].drop_duplicates()
    comp["site"] = comp.compartment.str.split(":").str[0]
    comp["cls"] = [CLASSES.get(str(family(c.split(":")[-1]))) for c in comp.compartment]
    rows = []
    for site, g in comp.groupby("site"):
        m = dict(zip(g.compartment, g.cls))
        for cut in (20, 50):
            d = pb[pb.compartment.isin([k for k, v in m.items() if v]) & (pb.n_cells >= cut)].copy()
            d["cls"] = d.compartment.map(m)
            has = d.groupby("donor").cls.agg(set)
            rows.append(dict(atlas=name, site=site, min_cells=cut, donors=int(pb.donor.nunique()),
                             endo=int(sum("endothelial" in s for s in has)), fib=int(sum("fibroblast" in s for s in has)),
                             mural=int(sum("mural" in s for s in has)),
                             endo_fib=int(sum({"endothelial", "fibroblast"} <= s for s in has)),
                             all_three=int(sum({"endothelial", "fibroblast", "mural"} <= s for s in has)), note=note))
    comp["atlas"] = name
    if name in USED:
        try:
            other, _ = grs.genome_table(USED[name] if USED[name] != "lung_natri_control_genome" else "lung_hlca_genome")
            ov = set(pb.donor.astype(str)) & set(other.donor.astype(str))
            for r in rows:
                r["overlap_with"] = USED[name]
                r["overlap_donors"] = len(ov)
        except Exception as e:  # noqa: BLE001
            for r in rows:
                r["overlap_with"] = f"{USED[name]}: {type(e).__name__}"
    return rows, comp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=2)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    res = Parallel(n_jobs=args.jobs)(delayed(run)(a) for a in ATLASES)
    rows, comps = [], []
    for r, c in res:
        rows += r if isinstance(r, list) else [r]
        if c is not None:
            comps.append(c)
    pd.DataFrame(rows).to_csv(OUT / "counts.csv", index=False)
    pd.concat(comps).to_csv(OUT / "compartments.csv", index=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
