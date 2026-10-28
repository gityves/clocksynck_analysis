"""A152 part 1 (lung, HLCA healthy) and D152.0: lung T under the A17 setup.

The A17 table of `lung_hlca` (healthy donors, the frozen compartments) with the T class replaced by the T rows of the
`hlca_tstates` builder (modes "T", "state", "cd8_state"; raw input $SCCC_DATA/atlas/hlca_tcell, A6 states folded in).
Every other row is the A17 row; rows of the builder for the other classes are compared with it (check).
Outputs (output/): lung_<variant>_pb.parquet for variants T50 (mode T, >= 50 cells: the A17 rule), T20 (>= 20 cells),
state, cd8_state; lung_repro_check.csv.

    python prepare_lung.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a152_lib as L  # noqa: E402
from _shared import Target  # noqa: E402

sys.path.insert(0, str(L.PAPER / "common" / "datasets" / "hlca_tstates"))


def main():
    import tstates as M
    T = Target("lung_hlca")
    donors = T.donors
    A = T.pb
    ids = [d.split(":", 1)[1] for d in donors.index]
    chk = []
    for variant, mode, tmin in (("T50", "T", 50), ("T20", "T", 20), ("state", "state", 20), ("cd8_state", "cd8_state", 20)):
        pb, _ = M.build("hlca_tstates", donors=ids, mode=mode, t_min=tmin)
        pb = pb[pb.donor.isin(donors.index)]
        if variant == "T50":
            pb = pb[~(pb.compartment.eq("lung:T") & (pb.n_cells < 50))]
        t = pb[pb.compartment.str.startswith("lung:T")]
        t = t.assign(stratum=t.donor.map(A.drop_duplicates("donor").set_index("donor").stratum))
        rest_b = pb[~pb.compartment.str.startswith("lung:T")].set_index(["donor", "compartment"]).sort_index()
        rest_a = A[~A.compartment.eq("lung:T")].set_index(["donor", "compartment"]).sort_index()
        com = rest_a.index.intersection(rest_b.index)
        num = ["n_cells", "depth"] + L.GENES
        d = np.abs(rest_a.loc[com, num].to_numpy(float) - rest_b.loc[com, num].to_numpy(float)).max()
        chk.append(dict(variant=variant, rows_a17_other=len(rest_a), rows_builder_other=len(rest_b), rows_common=len(com),
                        max_abs_diff_other=float(d)))  # d: max over rows and columns
        keepc = set(L.compartments_of("lung_hlca")) - {"lung:T"}
        # "b": every row from the builder (hlca_tcell: all groups, n_cells >= 1 before the donor sum); "m": A17 rows for the
        # other classes (hlca_full: groups of < 10 cells dropped), builder rows for T only
        b = pb[pb.compartment.isin(keepc) | pb.compartment.str.startswith("lung:T")]
        b = b.assign(stratum=b.donor.map(A.drop_duplicates("donor").set_index("donor").stratum))
        b[A.columns].reset_index(drop=True).to_parquet(L.OUT / f"lung_b_{variant}_pb.parquet", index=False)
        m = pd.concat([A[~A.compartment.eq("lung:T")], t[A.columns]], ignore_index=True)
        m.to_parquet(L.OUT / f"lung_m_{variant}_pb.parquet", index=False)
        chk[-1].update(rows_builder_in_a17_set=int(b.compartment.isin(keepc).sum()),
                       n_cells_diff_rows=int((rest_a.loc[com, "n_cells"] != rest_b.loc[com, "n_cells"]).sum()))
    pd.DataFrame(chk).to_csv(L.OUT / "lung_repro_check.csv", index=False)
    print(pd.DataFrame(chk).to_string())


if __name__ == "__main__":
    main()
