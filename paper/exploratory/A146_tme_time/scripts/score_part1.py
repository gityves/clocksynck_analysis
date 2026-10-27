"""A146 part 1: time-of-day programs of the mouse TNBC atlas in GSE233242 (bulk TPM), internal vs recorded time.

Registered in paper/docs/registered_experiments.md (A146). Outputs: paper/exploratory/A146_tme_time/output/part1_*.csv.
    python -I score_part1.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "4"
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from common.paths import DATA, DATASETS  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "output"
TPM = DATA / "validation" / "GSE233242" / "GSE233242_Processed_file_tpm.csv.gz"
TRUTH = DATASETS / "breast_timed_gse233242" / "resection_times.csv"
V25 = ROOT / "results" / "figure_1" / "output" / "windows" / "v25_5_patients.csv"
SEED, NPERM, NBOOT = 146, 10_000, 2_000

# name: (genes, mouse CT, residual-on)
PROG = {
    "P1_proliferation": ("MKI67 CCNB1 CENPA TOP2A AURKB PCNA MCM2 BUB1B BRCA2 MRE11", 1.5, None),
    "P2_MHC_I": ("B2M PSMB8 TAPBP TAP1 HLA-A HLA-B HLA-C", 22, None),
    "P3_interferon": ("STAT1 IRF1 IFNGR1 IFNGR2 IFNAR2 STAT3", 22, None),
    "P4_Tcell_activation": ("TCF7 LEF1 RUNX3 BCL6", 18, "P7_Tcell_abundance"),
    "P5_myeloid_AP": ("CIITA HLA-DRA HLA-DRB1 CD74", 18, "P8_myeloid_abundance"),
    "P6_phagocytosis": ("DNM2 CDC42 SYK FCGR1A FCGR2A", 18, "P8_myeloid_abundance"),
    "P7_Tcell_abundance": ("CD3D CD3E CD2 TRAC", 16, None),
    "P8_myeloid_abundance": ("CD68 CD14 CSF1R LYZ", 4, None),
    "P9_PDL1": ("CD274", 22, None),
    "P10_PD1": ("PDCD1", 4, "P7_Tcell_abundance"),
}
H = {k: (v[1] - 5) % 24 for k, v in PROG.items()}


def z(a):
    return (a - a.mean(0)) / a.std(0, ddof=1)


def scores(x, cols, patients, log):
    """Program scores (n patients x K) for one tissue; log collects found/dropped genes."""
    S = {}
    for k, (g, _, res) in PROG.items():
        gl = g.split()
        found = [q for q in gl if q in x.index]
        log.append(dict(program=k, n_genes=len(gl), n_found=len(found), found=" ".join(found),
                        dropped=" ".join(q for q in gl if q not in found),
                        scored=len(found) >= len(gl) / 2))
        if len(found) < len(gl) / 2:
            continue
        L = np.log2(x.loc[found, cols].T.to_numpy() + 1)
        S[k] = z(L).mean(1)
    for k, (_, _, res) in PROG.items():  # residual after raw scores are all available
        if res and k in S and res in S:
            a = np.c_[np.ones(len(patients)), S[res]]
            S[k] = S[k] - a @ np.linalg.lstsq(a, S[k], rcond=None)[0]
    return pd.DataFrame({k: z(S[k][:, None])[:, 0] for k in PROG if k in S}, index=patients)


def tstat(Y, C):
    """t of slope of each column of Y (n,K) on C (P,n,K) or (n,K): returns (P,K) or (K,)."""
    n = Y.shape[0]
    Yc = Y - Y.mean(0)
    Cc = C - C.mean(-2, keepdims=True)
    r = (Cc * Yc).sum(-2) / np.sqrt((Cc ** 2).sum(-2) * (Yc ** 2).sum(0))
    return r * np.sqrt((n - 2) / (1 - r ** 2))


def design(t, hs):
    return np.cos(2 * np.pi * (np.asarray(t)[..., None] - np.asarray(hs)) / 24)


def omnibus(Y, t, hs):
    return tstat(Y, design(t, hs)).mean(-1)


def holm(p):
    o = np.argsort(p)
    adj = np.empty(len(p))
    run = 0
    for i, j in enumerate(o):
        run = max(run, (len(p) - i) * p[j])
        adj[j] = min(1, run)
    return adj


def main():
    rng = np.random.default_rng(SEED)
    tr = pd.read_csv(TRUTH)
    tr = tr[tr.note.fillna("") == ""].set_index("patient_id")
    v = pd.read_csv(V25).set_index("patient")
    d = tr.join(v[["est_day"]], how="inner")
    assert len(d) == 42, len(d)
    assert np.allclose(d.resection_h, v.loc[d.index, "true_h"])
    pts = d.index.to_numpy()
    times = {"internal": d.est_day.to_numpy(), "recorded": d.resection_h.to_numpy()}
    x = pd.read_csv(TPM).drop_duplicates("gene_name").set_index("gene_name")

    log, Y = [], {}
    for tis in ("tumour", "normal"):
        lg = []
        Y[tis] = scores(x, d[f"tpm_column_{tis}"].to_numpy(), pts, lg)
        log += [dict(tissue=tis, **r) for r in lg]
    pd.DataFrame(log).to_csv(OUT / "part1_programs.csv", index=False)

    perms = np.array([rng.permutation(len(pts)) for _ in range(NPERM)])  # same for every program and time
    prog = {t: list(Y[t].columns) for t in Y}
    hk = {t: np.array([H[k] for k in prog[t]]) for t in Y}
    prog_n = [k for k in prog["normal"] if k != "P1_proliferation"]
    hn = np.array([H[k] for k in prog_n])

    tests, per, info = [], [], {}
    for tname, t in times.items():
        for tis in ("tumour", "normal"):
            cols = prog[tis] if tis == "tumour" else prog_n
            hs = hk[tis] if tis == "tumour" else hn
            Yk = Y[tis][cols].to_numpy()
            tk = tstat(Yk, design(t, hs))
            tn = tstat(Yk, design(t[perms], hs))  # (P,n,K) via fancy index -> (P,K)
            Sobs, Sn = tk.mean(), tn.mean(1)
            p_s = (1 + (Sn >= Sobs).sum()) / (NPERM + 1)
            pk = (1 + (tn >= tk).sum(0)) / (NPERM + 1)
            for k, a, b in zip(cols, tk, pk):
                per.append(dict(tissue=tis, time=tname, program=k, human_h=H[k], t_k=a, p_perm=b))
            info[(tname, tis)] = (Sobs, p_s, Sn, len(cols))
            if tis == "tumour":
                tests.append(dict(hypothesis="H146.1" if tname == "internal" else "H146.2", time=tname,
                                  n_patients=len(pts), n_programs=len(cols), S=Sobs, p=p_s, null_mean=Sn.mean(),
                                  null_sd=Sn.std(ddof=1), null_q95=np.quantile(Sn, .95),
                                  MDE=np.quantile(Sn, .95) + 0.84 * Sn.std(ddof=1)))
            else:
                per.append(dict(tissue="normal", time=tname, program="OMNIBUS_P2-P10", human_h=np.nan,
                                t_k=Sobs, p_perm=p_s))
    T = pd.DataFrame(tests)
    T["p_holm"] = holm(T.p.to_numpy())
    T["verdict"] = np.where((T.p_holm < 0.05) & (T.S > 0), "pass", "fail (not detected at this design)")
    T.to_csv(OUT / "part1_tests.csv", index=False)
    pd.DataFrame(per).to_csv(OUT / "part1_per_program.csv", index=False)

    # offset sweep (shown, not tested)
    cols = prog["tumour"]
    Yk = Y["tumour"][cols].to_numpy()
    base = np.array([(PROG[k][1] - 5) for k in cols], float)
    sw = []
    for off in (-6, -4, -2, 0, 2, 4, 6):
        for tname, t in times.items():
            sw.append(dict(offset_h=off, time=tname, S=omnibus(Yk, t, (base + off) % 24)))
    pd.DataFrame(sw).to_csv(OUT / "part1_offset_sweep.csv", index=False)

    # internal vs recorded, patient bootstrap
    hs = hk["tumour"]
    diff = np.empty(NBOOT)
    for b in range(NBOOT):
        i = rng.integers(0, len(pts), len(pts))
        diff[b] = omnibus(Yk[i], times["internal"][i], hs) - omnibus(Yk[i], times["recorded"][i], hs)
    Si, Sr = info[("internal", "tumour")][0], info[("recorded", "tumour")][0]
    pd.DataFrame([dict(S_internal=Si, S_recorded=Sr, diff=Si - Sr, boot_mean=diff.mean(), boot_sd=diff.std(ddof=1),
                       ci_lo=np.quantile(diff, .025), ci_hi=np.quantile(diff, .975),
                       p_boot_diff_le0=(1 + (diff <= 0).sum()) / (NBOOT + 1), n_boot=NBOOT)]
                 ).to_csv(OUT / "part1_internal_vs_recorded.csv", index=False)
    pd.set_option("display.width", 250)
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
