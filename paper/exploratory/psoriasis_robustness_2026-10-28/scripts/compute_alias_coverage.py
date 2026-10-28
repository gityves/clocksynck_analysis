"""Post hoc, descriptive: does mapping platform names with clocksync.genes.canonical change the core / cue coverage and the
headline offsets? Columns whose canonical symbol is a core or cue gene absent under its own name are renamed to it (first
column wins). Same estimator, series and weights as compute_robustness.py."""
from __future__ import annotations
import os
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import compute_robustness as C

H, Q, PD = C.H, C.Q, C.PD


def main():
    wfs = Q.cue_waveforms(H.reference_sets()["systemic"])
    cue = wfs["cue"]
    targets = set(H.GENES) | set(cue.genes)
    wanted = set(targets)
    orig = H.T.read_all_genes

    def slim(path, meta):
        names = pq.ParquetFile(path).schema.names
        al = {c for c in names if C.canonical(c) in targets}
        d = pd.read_parquet(path, columns=[c for c in names if c in set(meta) | wanted | al])
        return d.rename(columns={"BMAL1": "ARNTL"}) if "ARNTL" not in names else d
    H.T.read_all_genes = slim
    samples, w = H.gated_samples(C.OLD)
    S = C.Series(samples, w)
    rows, ren = [], []

    def run(S):
        full = lambda k, c: H.WF.subset([g for g in H.GENES if g in c.columns])
        fc = lambda k, c: H.restrict(cue, [g for g in cue.genes if g in c.columns])
        co, cu = S.offset(full), S.offset(fc)
        return co, cu, C.wrap(co - cu)
    base = run(S)
    for k in S.names:
        for i in (0, 1):
            x = S.pairs[k][i]
            m = {}
            for c in x.columns:
                t = C.canonical(c)
                if t in targets and t not in x.columns and t not in m.values() and c != t:
                    m[c] = t
            if i == 0:
                ren += [dict(series=k, platform_name=a, mapped_to=b) for a, b in m.items()]
            S.pairs[k] = tuple(x.rename(columns=m) if j == i else S.pairs[k][j] for j in (0, 1))
    new = run(S)
    cov = [dict(series=k, core=sum(g in S.pairs[k][0].columns for g in H.GENES), cue=sum(g in S.pairs[k][0].columns for g in cue.genes)) for k in S.names]
    out = pd.DataFrame(dict(estimate=["core", "cue", "core_minus_cue"], raw_names_h=base, canonical_names_h=new))
    out["shift_h"] = out.canonical_names_h - out.raw_names_h
    out.to_csv(C.OUT / "alias_mapping_effect.csv", index=False)
    pd.DataFrame(ren).to_csv(C.OUT / "alias_mapping_renamed_columns.csv", index=False)
    pd.DataFrame(cov).to_csv(C.OUT / "coverage_after_canonical_mapping.csv", index=False)
    print(out.round(4).to_string(index=False), "\n", pd.DataFrame(ren).to_string(index=False))


if __name__ == "__main__":
    main()
