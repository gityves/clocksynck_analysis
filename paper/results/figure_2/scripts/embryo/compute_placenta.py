"""A106 (primary): does the first-trimester placenta share the mother's clock state? (Arutyunyan 2023)

Pre-registered in docs/PLAN_internal_desync.md (archive) (A106: H106.1, H106.2 and its ambient check).
Units: pregnancies (>= 20 cells per pseudobulk, compartments with >= 8 pregnancies; stratum = assay).
  H106.1  prerequisite: the gate over the maternal compartments (decidua and placenta libraries summed) passes.
  H106.2  cross gate (`Dataset.gate(targets=, hosts=)`, 2000 permutations): fetal VCT and SCT x the maternal
          compartments; the pairing null permutes the donors of the fetal compartments only. Shared if p_pairing and
          p_acrophase < 0.05 and the ambient check does not explain it.
  Ambient check (pre-specified): synthetic fetal pseudobulks without a clock of their own (mean fetal profile mixed
          with the pregnancy's maternal soup, a = 10 x the compartment's XIST soup share), same cross gate; the pass
          stands only if the synthetic data do not pass (p_pairing >= 0.05) and the observed clock_cov exceeds the
          synthetic by > 5x (column `cov_ratio`).

    python compute_placenta.py [--subsample N] [--seed S]
Outputs: paper/results/figure_2/output/embryo/a106_gates.csv (H106.1, H106.2, ambient check), a106_soup.csv
Legacy:  legacy/desync/validation/natural_experiments/scripts/a106_placenta.py and natexp.py (cross_gate).
         local: soup share and the synthetic mixture (_shared.soup_share, _shared.synthetic_ambient).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402

from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

FOLDER = "embryo"
N_CROSS = 2000


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_cross = 100 if run["quick"] else N_CROSS
    pb, donors = data.load("placenta_arutyunyan", subsample=args.subsample, seed=args.seed)
    min_don = 5 if run["quick"] else 8  # a subsample has fewer pregnancies
    t = S.placenta_table(pb, min_donors=min_don)
    comps = sorted(t.compartment.unique())
    mat = [c for c in comps if c in S.MAT]
    prim = [c for c in S.PRIMARY_FET if c in comps]
    share = S.soup_share(pb, donors)
    share.rename("xist_soup_share").to_csv(out / "a106_soup.csv")
    print(t.groupby("compartment").donor.nunique().to_dict(), flush=True)
    rows = [S.run_gate(t[t.compartment.isin(mat)].reset_index(drop=True), "H106.1 maternal gate", 1, run["n_perm"])]
    ds = S.dataset(t, donors)
    g = ds.gate(n_perm=n_cross, seed=2, targets=prim, hosts=mat)
    rows.append(S.gate_row(g, test="H106.2 cross gate VCT+SCT x maternal", donors=g.n_donors))
    syn = S.dataset(S.synthetic_ambient(t, prim, share), donors).gate(n_perm=n_cross, seed=2, targets=prim, hosts=mat)
    rows.append(S.gate_row(syn, test="H106.2 ambient check (synthetic)", donors=syn.n_donors))
    res = pd.DataFrame(rows)
    res["cov_ratio"] = float("nan")
    res.loc[res.test.str.startswith("H106.2 cross"), "cov_ratio"] = res.clock_cov.iloc[1] / res.clock_cov.iloc[2]
    res.to_csv(out / "a106_gates.csv", index=False)
    print(res.drop(columns=["compartments"]).round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
