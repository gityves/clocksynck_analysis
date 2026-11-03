"""Post hoc, descriptive (2026-11-03, user question after A201 was scored): is there a pre-existing SMC - fibroblast clock phase offset
under aligned conditions, and does the misaligned group keep it? Offset = arg sum_g z_g^SMC conj(z_g^fib) (h), same cosinor as A201;
also per gene. Aligned controls split into the 4 Tempo libraries (99990-99993, male) and the other 10 as a consistency check.
SD by library jackknife within each group. Output: output/posthoc_baseline_offset.csv, output/posthoc_gene_offsets.csv.
"""
import numpy as np
import pandas as pd

import score_a201 as S

TEMPO = {"99990", "99991", "99992", "99993"}


def offset(lib, C, cond, keep=None):
    zs, zf = S.fit(lib, C, cond, "SMC", keep), S.fit(lib, C, cond, "Fibroblast", keep)
    return float(np.angle(np.sum(zs * np.conj(zf))) / S.W), zs, zf


def jack(lib, C, cond, keep):
    libs = lib.loc[(lib.cond == cond).to_numpy() & keep, "sample_id"].astype(str).unique()
    v = np.array([offset(lib, C, cond, keep & (lib.sample_id.astype(str) != s).to_numpy())[0] for s in libs])
    n = len(v)
    return float(np.sqrt((n - 1) / n * np.sum(S.wrap(v - np.angle(np.exp(1j * S.W * v).mean()) / S.W) ** 2)))


def main():
    o, X, genes = S.load()
    keep = o.cell_type.isin(S.MAIN).to_numpy()
    lib, C = S.pseudobulk(o[keep].reset_index(drop=True), X[keep])
    allk = np.ones(len(lib), bool)
    tempo = lib.sample_id.astype(str).isin(TEMPO).to_numpy()
    rows, genes_rows = [], []
    for lab, cond, k in (("aligned, all 14", "ali", allk), ("aligned, Tempo 4 (male)", "ali", tempo),
                         ("aligned, other 10", "ali", ~tempo), ("misaligned, all 8", "mis", allk),
                         ("aligned, male", "ali", (lib.sex == "male").to_numpy()), ("misaligned, male", "mis", (lib.sex == "male").to_numpy()),
                         ("aligned, female", "ali", (lib.sex == "female").to_numpy()), ("misaligned, female", "mis", (lib.sex == "female").to_numpy())):
        try:
            off, zs, zf = offset(lib, C, cond, k)
            sd = jack(lib, C, cond, k)
        except np.linalg.LinAlgError:
            continue
        rows.append(dict(group=lab, smc_minus_fib_h=off, sd_jack_h=sd))
        if lab in ("aligned, all 14", "misaligned, all 8"):
            for g, a, b in zip(genes, zs, zf):
                genes_rows.append(dict(group=lab, gene=g, amp_smc=abs(a), amp_fib=abs(b),
                                       smc_minus_fib_h=float(np.angle(a * np.conj(b)) / S.W)))
    r = pd.DataFrame(rows)
    r.to_csv(S.OUT / "posthoc_baseline_offset.csv", index=False)
    gtab = pd.DataFrame(genes_rows).pivot(index="gene", columns="group", values="smc_minus_fib_h")
    gtab.to_csv(S.OUT / "posthoc_gene_offsets.csv")
    print(r.round(3).to_string(index=False)); print(gtab.round(2).to_string())


if __name__ == "__main__":
    main()
