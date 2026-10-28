"""A152 part 2 (heart, breast): reproduction check of the per-cell clock counts against the atlas pseudobulk, A6 residency
state per cell, state pseudobulks on the A17 grouping.

  1. check: per-cell sums over every A17 lymphocyte-class (donor, class) vs the stored A17 table (counts, n_cells, depth);
  2. score: score = mean x(CD69, ITGAE, CXCR6) - mean x(SELL, CCR7, KLF2, S1PR1); tertiles within block (assay x
     suspension type; one study per atlas), ties at random (seed stream 152), on all lymphocyte cells pulled (mast cells,
     which the A17 regexes put in the class, are not scored and not in the state pseudobulks);
  3. state pseudobulks (>= 20 cells) of the cells of the donor's majority assay (and nuclei, heart; left ventricle, heart).
Outputs (output/): part2_<tissue>_pb.parquet (A17 table with the lymphocyte classes replaced by the two states; also the
`_cd8` and `_blockLV` variants), part2_repro_check.csv, part2_cells_summary.csv.

    python prepare_part2.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a152_lib as L  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402


def main():
    chk, summ = [], []
    for ti, tissue in enumerate(("heart", "breast")):
        cfg = L.CFG[tissue]
        cells = pd.read_parquet(DATA / "cells" / cfg["cells"] / "cells.parquet")
        cells["donor"] = cfg["atlas"] + ":" + cells.donor_id.astype(str)
        cells["cls"] = L.classes_of(tissue, cells.cell_type).to_numpy()
        batch = L.donor_batch_cached(tissue)
        cells["cbatch"] = cells.assay.astype(str) + (("|" + cells.suspension_type.astype(str)) if tissue == "heart" else "")
        keep = (cells.cbatch == cells.donor.map(batch)).to_numpy().copy()  # donor's majority assay (heart: nuclei only)
        if tissue == "heart":
            keep = keep & cells.tissue.eq("heart left ventricle").to_numpy()
        cells["keep"] = keep
        # 1. reproduction on all cells of the replaced classes (mast cells included)
        pb0, _ = data.load(cfg["target"] if tissue == "breast" else "heart_global")
        a = L.pseudobulk(cells[cells.keep], cells[cells.keep].cls, tissue, min_cells=20)
        ref = pb0[pb0.compartment.isin(cfg["replaced"])].set_index(["donor", "compartment"]).sort_index()
        mine = a.set_index(["donor", "compartment"]).sort_index()
        com = ref.index.intersection(mine.index)
        num = ["n_cells", "depth"] + L.GENES
        d = np.abs(ref.loc[com, num].to_numpy(float) - mine.loc[com, num].to_numpy(float))
        for cl in cfg["replaced"]:
            m = com.get_level_values(1) == cl
            dd = d[m]
            chk.append(dict(tissue=tissue, cls=cl, rows_atlas=int((ref.index.get_level_values(1) == cl).sum()),
                            rows_mine=int((mine.index.get_level_values(1) == cl).sum()), rows_common=int(m.sum()),
                            max_diff_counts=float(dd[:, 2:].max()) if m.any() else np.nan,
                            max_diff_n_cells=float(dd[:, 0].max()) if m.any() else np.nan,
                            max_diff_depth=float(dd[:, 1].max()) if m.any() else np.nan))
        # 2. scores on lymphocyte cells (all pulled, not only the analysed ones; A6 scored every T cell of the dataset)
        ly = cells[~cells.cell_type.str.contains(L.NOT_LYMPH)].copy()
        ly["score"] = L.score_cells(ly)
        for variant in ("all", "analysed"):
            sub = ly if variant == "all" else ly[ly.keep]
            blk = (sub.assay.astype(str) + "|" + sub.suspension_type.astype(str)).to_numpy()
            st = L.tertile(sub.score.to_numpy(), blk, np.random.default_rng([L.SEED, ti, 0]))
            ly.loc[sub.index, f"state_{variant}"] = st
        cells.loc[ly.index, ["score", "state_all", "state_analysed"]] = ly[["score", "state_all", "state_analysed"]]
        summ.append(ly.groupby(["assay", "suspension_type", "state_all"]).size().rename("cells").reset_index().assign(tissue=tissue))
        # 3. state pseudobulks replace the lymphocyte classes
        rest = pb0[~pb0.compartment.isin(cfg["replaced"])]
        for variant, col, sel in (("", "state_all", None), ("_blockLV", "state_analysed", None), ("_cd8", "state_all", cfg["cd8"])):
            if variant == "_cd8" and sel is None:
                continue
            c = cells[cells.keep & cells.score.notna()]
            if sel:
                c = c[c.cls == sel]
            lab = c[col].map({"resident_like": f"{tissue}:lymph_resident_like", "circulating_like": f"{tissue}:lymph_circulating_like"})
            if variant == "_cd8":
                lab = lab.str.replace("lymph_", "CD8_")
            st = L.pseudobulk(c, lab, tissue, min_cells=20)
            out = pd.concat([rest, st], ignore_index=True)
            out.to_parquet(L.OUT / f"part2_{tissue}{variant}_pb.parquet", index=False)
    pd.DataFrame(chk).to_csv(L.OUT / "part2_repro_check.csv", index=False)
    pd.concat(summ).to_csv(L.OUT / "part2_cells_summary.csv", index=False)
    print(pd.DataFrame(chk).to_string())


if __name__ == "__main__":
    main()
