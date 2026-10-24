"""WP3 scout (A111, A112, A129): design table of the candidate datasets. No hypothesis statistic is computed.

One row per dataset x arm: mice, mice per time point (min-max over the sampled times), tissues, mice that carry
every tissue / at least two tissues (hosts need >= 2), count type, clock genes present (of 14) and with median
count >= 5 per tissue (`tissue_qc` rule). Reads the standard tables built by `common/datasets/<name>/build.py`.

    paper/common/run.sh paper/exploratory/WP3_scout/scout_wp3.py      -> output/scout_wp3.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import data  # noqa: E402

OUT = Path(__file__).resolve().parent / "output"
GENES = list(cs.gtex().genes)
# name -> (accession, count type, arm column, time column, role)
SETS = {
    "manella2021": ("GSE159135", "raw 3' MARS-seq reads (genes pre-filtered, ~11.6k)", "arm", "zt", "A111 primary, H111.3, A129"),
    "xin2021": ("GSE150380/150381/151221/151228", "TPM only (pseudo-counts)", "arm", "zt", "A111 replication"),
    "greenwell2019": ("GSE118967", "raw 3' counts + TPM", "arm", "zt", "A112"),
    "jouffe2022_bmal1ko": ("GSE150102", "raw counts", "stratum", "zt", "A129 candidate (whole-body KO, one ZT)"),
    "vipko_gse296764": ("GSE296764", "raw counts", "stratum", "ct", "A129 candidate (central VIP-Cre KO, CT0/CT12)"),
}


def main() -> None:
    rows = []
    for name, (acc, ctype, arm_col, t_col, role) in SETS.items():
        pb, donors = data.load(name)
        pb = pb.merge(donors[[arm_col, t_col]], left_on="donor", right_index=True) if arm_col not in pb else pb
        pb["_t"] = pb.donor.map(donors[t_col])
        pb["_arm"] = pb.donor.map(donors[arm_col]) if arm_col != "stratum" else pb["stratum"]
        for arm, g in pb.groupby("_arm"):
            nt = g.groupby("donor").compartment.nunique()
            ntis = g.compartment.nunique()
            per_t = g.drop_duplicates("donor").groupby("_t").size()
            med = g.groupby("compartment")[GENES].median()
            rows.append(dict(
                dataset=name, accession=acc, arm=arm, role=role, count_type=ctype,
                n_mice=g.donor.nunique(), time_col=t_col, n_times=per_t.size,
                mice_per_time_min=int(per_t.min()), mice_per_time_max=int(per_t.max()),
                mice_per_time="; ".join(f"{int(k)}:{int(v)}" for k, v in per_t.items()),
                tissues=",".join(sorted(g.compartment.unique())), n_tissues=ntis,
                mice_all_tissues=int((nt == ntis).sum()), mice_ge2_tissues=int((nt >= 2).sum()),
                clock_genes_present=int(sum(c in pb.columns for c in GENES)),
                clock_genes_median_ge5_all_tissues=int((med >= 5).all(axis=0).sum()),
                weakest_tissue_gene_median="; ".join(f"{t}:{med.loc[t].min():.0f}" for t in med.index)))
    OUT.mkdir(exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "scout_wp3.csv", index=False)
    print(df.drop(columns=["accession", "role", "count_type", "mice_per_time", "time_col"]).to_string())


if __name__ == "__main__":
    main()
