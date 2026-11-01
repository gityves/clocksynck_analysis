"""A181 step 0: reproduce the A107b gates of H107b.1 and H107b.2 (Cao 2020) before the responsive-gene gate (H181.2)."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import pandas as pd  # noqa: E402

from common import cli, data  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402
from results.figure_2.scripts.embryo.compute_fetal_atlases import CAO, MIN_DON  # noqa: E402

OUT = HERE.parents[1] / "output"


def main():
    cli.setup(cli.parser(__doc__).parse_args([]))
    ref = pd.read_csv(HERE.parents[3] / "results/figure_2/output/embryo/a107b_gates.csv").set_index("test")
    rows = []
    for j, (test, (lev, comps, mc)) in enumerate(list(CAO.items())[:2]):
        pb, _ = data.load("fetal_cao", level=lev)
        r = S.run_gate(S.eligible(pb, mc, MIN_DON, compartments=comps), test, 40 + j, 1000)
        for k in ("p_pairing", "p_acrophase", "clock_share", "clock_cov", "donors", "n_pairs"):
            rows.append(dict(test=test, quantity=k, new=r[k], a107b=ref.loc[test, k],
                             match=bool(abs(r[k] - ref.loc[test, k]) < 1e-9)))
    t = pd.DataFrame(rows)
    OUT.mkdir(exist_ok=True)
    t.to_csv(OUT / "repro_a107b.csv", index=False)
    print(t.to_string(index=False))


if __name__ == "__main__":
    main()
