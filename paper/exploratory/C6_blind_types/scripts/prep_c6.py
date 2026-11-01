"""A182 data preparation (blind to the clock statistic): the C6 GEO series matrices (colon GSE39582 / GSE17536, prostate
GSE116918, lung squamous GSE157009 / GSE157010, head and neck GSE65858 descriptive) parsed with the A171 parser into
../output/{key}_meta.parquet and {key}_expr.parquet (tumours only per the registered filters of header_summary.csv; log2 if max > 100;
probes averaged per symbol; G2M / purity-proxy scores = mean gene z over the kept tumours). Writes ../output/prep_log.csv
(scale decision, probe counts, n per cohort). No coherence here. GSE53625 (no symbol annotation) is dropped (registered).
Run: python -I prep_c6.py   (reads downloaded text files only; nothing downloaded is executed)"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A171_luad_lihc" / "scripts"))
sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts"))
sys.path.insert(0, str(HERE))
import a168_lib as L  # noqa: E402
import header_c6 as H  # noqa: E402
import prep_a171 as Q  # noqa: E402

OUT = HERE.parent / "output"
L.OUT = OUT
Q.P.D = H.D


def probe_symbols(gpl):
    a = H.platform_symbols(gpl)
    Q.P.LOG.append((gpl, "probes_with_symbol", a.ID.nunique()))
    Q.P.LOG.append((gpl, "probes_with_several_symbols", int((a.groupby("ID").size() > 1).sum())))
    return a.set_index("ID")["sym"]


Q.probe_symbols = probe_symbols  # expression_file resolves probe_symbols from prep_a171's namespace

if __name__ == "__main__":
    for typ, key, gse, folder in H.COHORTS:
        if key == "gse53625":
            continue
        f = next((H.D / folder).glob("*_series_matrix.txt.gz"))
        head, _ = H.read_header(f)
        s, gpl, _, _ = H.samples(head)
        keep, _ = H.tumour_rule(key)[1:]
        m = pd.DataFrame([d for d in s if keep(d)])
        extra = {"dataset": "dataset"} if key == "gse39582" else {}
        m = m.set_index("gsm")[["title", *extra]]
        # as prep_a171.expression_file (linear -> log2), reading the cohort file directly
        ex = pd.read_csv(f, sep="\t", comment="!", index_col=0, quotechar='"', na_values=["null", "NA", ""]).astype(np.float32)
        ex.index = ex.index.astype(str)
        mx, mn = float(np.nanmax(ex.values)), float(np.nanmin(ex.values))
        if mx > 100:
            ex = np.log2(ex.clip(lower=float(ex.values[ex.values > 0].min())))
        Q.P.LOG += [(gse, "deposited_max", mx), (gse, "deposited_min", mn),
                    (gse, "deposited_scale", "linear -> log2" if mx > 100 else "log2 (kept)"), (gse, "n_samples_matrix", ex.shape[1])]
        sym = probe_symbols(gpl[0]); sym = sym[sym.index.isin(ex.index)]
        ex = ex.loc[sym.index].groupby(sym.to_numpy()).mean()
        Q.P.LOG += [(gse, "probe_symbol_rows", len(sym)), (gse, "genes", len(ex)), (gse, "n_tumours_kept", len(m))]
        Q.P.finish(gse, key, m, ex)
    pd.DataFrame(Q.P.LOG, columns=["series", "item", "value"]).to_csv(OUT / "prep_log.csv", index=False)
    for k in ("gse39582", "gse17536", "gse116918", "gse157009", "gse157010", "gse65858"):
        m = pd.read_parquet(OUT / f"{k}_meta.parquet")
        print(k, len(m), m.dataset.value_counts().to_dict() if "dataset" in m else "", flush=True)
