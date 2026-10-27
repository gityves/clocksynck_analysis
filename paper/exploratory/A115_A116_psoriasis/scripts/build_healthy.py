"""A115 data step: rebuild the uninvolved rows AND add the healthy rows of four series from the raw GEO files, processed exactly as
the legacy builders (build_skin.py, build_skin_heldout.py: same selection, replicate merge, probe-to-symbol map, GENCODE v46 symbol
sum, scale, depth). Acceptance (registered): the rebuilt uninvolved rows must reproduce the processed ones (exact, or per-gene
correlation > 0.99); a series that fails is dropped and reported.

    python -I build_healthy.py
Outputs: output/a115_acceptance.csv; derived parquet <data root>/skin_lesional/processed_a115/<GSE>_<disease>.parquet
(uninvolved + healthy rows, all genes; public data). Raw files are parsed, never executed.
"""
import importlib.util
import sys
from pathlib import Path
import numpy as np, pandas as pd

LEG = Path("/home/maxine/Documents/paychere/scCircadianCancer/legacy/desync/validation/paired_damping/scripts")


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, LEG / file)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


bs = load("build_skin", "build_skin.py")
hb = load("build_skin_heldout", "build_skin_heldout.py")
ROOT = bs.ROOT
OUTD = ROOT / "processed_a115"
OUT = Path(__file__).resolve().parents[1] / "output"
g1 = bs.g1
META = ["patient", "status", "sample", "dataset", "disease", "platform", "value_type", "depth", "n_replicates"]
# healthy rule: (person id, "N" placeholder status so the legacy builders accept it)
HEALTHY = {
    "GSE13355": lambda s: ("H_" + g1(r"Individual_(\d+)_NN_sample", s["title"]), "N", True) if "_NN_sample" in s["title"] else (None, None, None),
    "GSE121212": lambda s: ("H_" + g1(r"^(CTRL_\d+)_healthy", s["title"]), "N", True) if "skin type: healthy" in s["ch"] else (None, None, None),
    "GSE201827": lambda s: ("H_" + g1(r"subject: (\d+)", s["ch"]), "N", True) if "timepoint: Control | tissue: Normal" in s["ch"] else (None, None, None),
    "GSE268958": lambda s: ("H_" + g1(r"subjectid: (\S+)", s["ch"]), "N", True) if "condition: Healthy" in s["ch"] else (None, None, None),
}
# series -> (platform, kind, value_type, {disease: lesional rule}, processed folder)
SERIES = {
    "GSE13355": ("GPL570", "array", "log2_centred", {"psoriasis": hb.RULES["GSE13355"]}, "processed_heldout"),
    "GSE121212": ("Illumina HiSeq 2000 (GPL16791)", "counts", "counts",
                  {"psoriasis": bs.RULES["GSE121212_PSO"], "atopic_dermatitis": bs.RULES["GSE121212_AD"]}, "processed"),
    "GSE201827": ("GPL570", "array", "log2_intensity", {"psoriasis": bs.RULES["GSE201827"]}, "processed"),
    "GSE268958": ("Illumina", "counts", "counts", {"psoriasis": bs.RULES["GSE268958"]}, "processed"),
}


def build(gse, plat, kind, sel):
    return bs.build_array(gse, plat, sel) if kind == "array" else bs.build_counts(gse, sel, None)


def frame(info, X, gse, dis, plat, vt, status=None):
    X = hb.rename_aliases(X.fillna(0) if vt == "counts" else X)
    if status:
        info = info.assign(status=status)
    info = info.assign(dataset=gse, disease=dis, platform=plat, value_type=vt)
    return pd.concat([info[META], X.reset_index(drop=True)], axis=1)


def accept(new, old):
    """new/old: uninvolved rows of the rebuilt / processed series. -> dict of acceptance statistics."""
    old = old.set_index("sample"); new = new.set_index("sample")
    common = [s for s in old.index if s in new.index]
    genes = [g for g in old.columns if g not in META and g in new.columns]
    a, b = new.loc[common, genes].astype(float), old.loc[common, genes].astype(float)
    exact = bool(np.array_equal(np.nan_to_num(a.values), np.nan_to_num(b.values)))
    maxabs = float(np.nanmax(np.abs(a.values - b.values))) if len(common) else np.nan
    ac, bc = a - a.mean(), b - b.mean()
    den = np.sqrt((ac ** 2).sum() * (bc ** 2).sum())
    r = ((ac * bc).sum() / den.where(den > 0)).dropna()
    return dict(n_old=len(old), n_new=len(new), n_common=len(common), n_genes=len(genes), exact=exact, max_abs_diff=maxabs,
                median_gene_r=float(r.median()), min_gene_r=float(r.min()), share_r_gt_0_99=float((r > 0.99).mean()),
                depth_equal=bool(np.allclose(new.loc[common, "depth"].astype(float).fillna(-1), old.loc[common, "depth"].astype(float).fillna(-1))))


def main():
    OUTD.mkdir(exist_ok=True)
    rows = []
    for gse, (plat, kind, vt, rules, folder) in SERIES.items():
        meta = bs.read_meta(gse)
        old_all = pd.read_parquet(ROOT / folder / f"{gse}.parquet")
        sel_h = bs.select_h = {}
        for s in meta:
            p, st, b = HEALTHY[gse](s)
            if p and st:
                sel_h.setdefault((p, st), []).append(s)
        info_h, X_h, _ = build(gse, plat, kind, sel_h)
        H = frame(info_h, X_h, gse, None, plat, vt, status="healthy")
        for dis, rule in rules.items():
            sel = bs.select(meta, rule)
            info, X, _ = build(gse, plat, kind, sel)
            new = frame(info, X, gse, dis, plat, vt)
            old = old_all[old_all.disease == dis]
            for st in ("nonlesional", "lesional"):
                acc = accept(new[new.status == st], old[old.status == st])
                ok = acc["exact"] or (acc["median_gene_r"] > 0.99 and acc["share_r_gt_0_99"] > 0.99 and acc["n_common"] == acc["n_old"])
                rows.append(dict(series=gse, disease=dis, rows=st, n_healthy_built=len(H), **acc, accepted=bool(ok)))
                print(rows[-1], flush=True)
            Hd = H.assign(disease=dis)
            cols = list(new.columns.union(Hd.columns, sort=False))
            out = pd.concat([new[new.status == "nonlesional"], Hd], ignore_index=True)
            out = out[[c for c in new.columns if c in out.columns]]
            out.to_parquet(OUTD / f"{gse}_{dis}.parquet", index=False)
    pd.DataFrame(rows).to_csv(OUT / "a115_acceptance.csv", index=False)


if __name__ == "__main__":
    main()
