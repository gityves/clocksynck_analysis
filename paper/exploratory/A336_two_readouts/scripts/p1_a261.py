"""A336 part 1, A261 H261.3 (night-shift programs, Figure 1d): the program statistic is waveform-free; the tissue-clock
hypothesis takes the core-loop constants of A257 (s and R from `part1/a257.csv` core values) in place of the 14-gene
ones (2.68 h, 0.233). The all-rhythmic program's class is recomputed under both. No drop or trio readout (registered).

    python p1_a261.py
Output: ../output/part1/a261.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(PAPER / "exploratory/A261_nightshift_programs/scripts"))
import a336_common as C  # noqa: E402
import score_a261 as S  # noqa: E402
from compute_a257 import load  # noqa: E402


def main():
    ref = S.reference().set_index("gene")
    d = S.arrays()
    al, mis = d["ctl"], d["mis"]
    genes = [g for g in ref.index if g in al.columns and g not in S.CLOCK]
    ref = ref.loc[genes]
    rhy = list(ref.index[ref.p < 0.01])
    Cc = 2 * ref.amp * np.cos(S.W * (19 - ref.acro_h))
    Ss = 2 * ref.amp * np.sin(S.W * (19 - ref.acro_h))
    zA = S.per_subject_z(al[rhy].to_numpy(float), Cc[rhy].to_numpy(), Ss[rhy].to_numpy())
    zM = S.per_subject_z(mis[rhy].to_numpy(float), Cc[rhy].to_numpy(), Ss[rhy].to_numpy())
    stored = pd.read_csv(PAPER / "exploratory/A261_nightshift_programs/output/tests.csv").set_index("test").loc["H261.3"]
    r14 = S.read(zA, zM)
    ok = r14["cls"] == stored.value and np.isclose(r14["shift_h"], float(stored.shift_h), atol=1e-6)
    if not ok:
        raise SystemExit("A261 H261.3 not reproduced")
    a257 = pd.read_csv(HERE.parents[1] / "output/part1/a257.csv").set_index("id")
    s_core, R_core = float(a257.loc["A257_H257.2", "value_core"]), 1 - float(a257.loc["A257_H257.3", "value_core"])
    S.HYP["tissue_clock"] = (s_core, R_core)
    rc = S.read(zA, zM)
    status = "holds" if rc["cls"] == r14["cls"] else "falls"
    row = dict(id="A261_H261.3", figure_panel="Fig 1d", statistic="all-rhythmic program class (behaviour / tissue clock / central)",
               value_14_stored=stored.value, value_14=r14["cls"], reproduced=True, value_core=rc["cls"],
               core_criterion="class of the all-rhythmic program with the tissue-clock hypothesis at the core constants",
               core_status=status, reading=C.reading(status, None) if status == "holds" else "falls: class changes",
               blind="core constants not blind (A335 1a)", gene_drop="none",
               notes=(f"program shift {rc['shift_h']:.2f} h, ratio {rc['ratio']:.3f} (waveform-free); tissue-clock hypothesis "
                      f"14 genes (2.68 h, 0.233): p {r14['p_tissue_clock']:.2e}; core ({s_core:.2f} h, {R_core:.3f}): p {rc['p_tissue_clock']:.2e}; "
                      f"behaviour p {rc['p_behaviour']:.3f}; central p {rc['p_central']:.2e}"),
               script="exploratory/A336_two_readouts/scripts/p1_a261.py")
    C.write_rows("a261", [row])
    print(row)


if __name__ == "__main__":
    main()
