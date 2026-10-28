"""A152 post hoc diagnostic (descriptive): the A17 lymphocyte classes without the mast cells their regexes capture.

heart: `heart:lymphocyte` rebuilt from lymphocyte + natural killer cells, mast cells as `heart:mast`;
breast: `breast:T_other` rebuilt without mast cells, `breast:mast` separate. Same A17 setup (donors, strata, learned
anchors, majority batch, left ventricle for heart, >= 20 cells), point estimates (n_boot = 0). Per-cell parquets of
extract_cells.py. Output output/mast_diagnostic.csv; mast share = mast cells / cells of the original row (median over donors).

    python mast_diagnostic.py
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a152_lib as L  # noqa: E402
from _shared import Target  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

TARGET = {"heart": ("heart", "heart:lymphocyte"), "breast": ("organ_breast", "breast:T_other")}


def main():
    warnings.filterwarnings("ignore")
    rows = []
    for tissue, (tgt, orig) in TARGET.items():
        cfg = L.CFG[tissue]
        cells = pd.read_parquet(DATA / "cells" / cfg["cells"] / "cells.parquet")
        cells["donor"] = cfg["atlas"] + ":" + cells.donor_id.astype(str)
        batch = L.donor_batch_cached(tissue)
        cb = cells.assay.astype(str) + (("|" + cells.suspension_type.astype(str)) if tissue == "heart" else "")
        keep = (cb == cells.donor.map(batch)).to_numpy().copy()
        if tissue == "heart":
            keep = keep & cells.tissue.eq("heart left ventricle").to_numpy()
        cells = cells[keep]
        cls = L.classes_of(tissue, cells.cell_type)
        mast = cells.cell_type.str.contains(L.NOT_LYMPH).to_numpy()
        lab_new = cls.where(~mast, f"{tissue}:mast")
        T = Target(tgt)
        pb0 = T.pb
        base = pb0[~pb0.compartment.eq(orig)]
        new = L.pseudobulk(cells[cls.eq(orig).to_numpy()], lab_new[cls.eq(orig).to_numpy()], tissue, 20)
        pb = pd.concat([base, new], ignore_index=True)
        anc = L.anchors(tgt)
        t0 = L.fit_table(pb0, T.donors, anc, tgt)
        t1 = L.fit_table(pb, T.donors, anc, tgt)
        # mast share of the cells of the original row, per donor with an original row (>= 20 cells)
        o = cells[cls.eq(orig).to_numpy()]
        sh = o.assign(m=mast[cls.eq(orig).to_numpy()]).groupby("donor").m.agg(["mean", "size"])
        sh = sh[sh["size"] >= 20]["mean"]
        for name, tab, p in ((f"{orig} (A17 original)", t0, pb0[pb0.compartment == orig]),
                             (f"{orig} (rebuilt without mast)", t1, pb[pb.compartment == orig]),
                             (f"{tissue}:mast", t1, pb[pb.compartment == f"{tissue}:mast"])):
            key = orig if name.startswith(orig) else name
            n_don = int(p.donor.nunique())
            rows.append(dict(tissue=tissue, row=name, n_donors=n_don, enough_donors=n_don >= 15,
                             coupling=tab.loc[key, "coupling"] if key in tab.index and n_don >= 15 else np.nan,
                             coupling_if_fitted=tab.loc[key, "coupling"] if key in tab.index else np.nan,
                             boundary=tab.loc[key, "boundary"] if key in tab.index else np.nan,
                             median_cells=p.n_cells.median(), median_depth=p.depth.median(),
                             mast_share_median=float(sh.median()), mast_share_q25=float(sh.quantile(.25)),
                             mast_share_q75=float(sh.quantile(.75)), donors_in_share=len(sh)))
    out = pd.DataFrame(rows)
    out.to_csv(L.OUT / "mast_diagnostic.csv", index=False)
    print(out.round(4).to_string())


if __name__ == "__main__":
    main()
