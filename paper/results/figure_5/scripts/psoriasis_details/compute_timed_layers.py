"""A67 and A76: the healthy epidermis-dermis clock shift in timed skin (GSE205155, 11 subjects x 7 times, both layers).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A67 H67.1, A76 H76.1). Known sampling times, no host.
A67: per gene and layer a cosinor on clock time with subject intercepts (log2 intensity; F-test p, amplitude, acrophase).
Shift = wrap(acrophase dermis - acrophase epidermis) (> 0: the epidermis peaks earlier). H67.1: among the 14 clock genes
rhythmic in both layers (p < 0.01 each) the circular mean shift > 0 with the subject-bootstrap 95 % interval > 0 and the
circular SD of the per-gene shifts < 2 h. Descriptive: all genes rhythmic in both layers.
A76: per subject, the clock-plane state (GTEx waveform, one projector for both layers), complex regression on clock time
per layer, offset = arg(a_epidermis / a_dermis) (h, > 0: epidermis ahead). Epidermal indices (subject means, log2):
differentiation (differentiated minus basal set), proliferation. H76.1: Spearman(differentiation, offset) > 0, one-sided
permutation p < 0.05 (10,000).

    python compute_timed_layers.py [--subsample N]      (--subsample: fewer bootstrap draws and permutations)
Outputs: paper/results/figure_5/output/psoriasis_details/a67_clock_genes.csv, a67_genomewide_rhythmic.csv, a67_verdict.csv,
  a76_verdict.csv, a76_summary.csv (public GEO data; per-subject offsets are not written)
Legacy: legacy/desync/validation/paired_damping/scripts/a67_timed_layers.py, a76_subject_offsets.py
Uses a local per-gene cosinor with an F test (local: not in clocksync; `clocksync.timed.cosinor` has no test).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402

W, wrap = H.W, H.S.wrap


def circ(h):
    z = np.exp(1j * W * np.asarray(h))
    R = abs(z.mean())
    return float(np.angle(z.mean()) / W), float(np.sqrt(-2 * np.log(max(R, 1e-12))) / W)


def a67(s, m, out, n_boot, seed=67):
    clock = [g for g in H.GENES if g in m.index]
    r = H.per_layer(s, m, clock)
    d = pd.DataFrame({"acro_epi": r["epidermis"][0], "amp_epi": r["epidermis"][1], "p_epi": r["epidermis"][2],
                      "acro_derm": r["dermis"][0], "amp_derm": r["dermis"][1], "p_derm": r["dermis"][2]}, index=clock)
    d["shift_h"] = wrap(d.acro_derm - d.acro_epi)
    d["log2_amp_ratio"] = np.log2(d.amp_epi / d.amp_derm)
    rhy = (d.p_epi < 0.01) & (d.p_derm < 0.01)
    mean, sd = circ(d.shift_h[rhy])
    rng = np.random.default_rng(seed)
    subj = s.subject.unique()
    boot = []
    for _ in range(n_boot):
        pick = rng.choice(subj, len(subj))
        sb = pd.concat([s[s.subject == p].assign(subject=f"{p}#{j}") for j, p in enumerate(pick)], ignore_index=True)
        rb = H.per_layer(sb, m, list(d.index[rhy]))
        boot.append(circ(wrap(rb["dermis"][0] - rb["epidermis"][0]))[0])
    dd = wrap(np.asarray(boot) - mean)
    lo, hi = mean + np.quantile(dd, 0.025), mean + np.quantile(dd, 0.975)
    d["rhythmic_both"] = rhy
    d.rename_axis("gene").to_csv(out / "a67_clock_genes.csv")
    genes = list(m.index)
    g = H.per_layer(s, m, genes)
    gw = pd.DataFrame({"acro_epi": g["epidermis"][0], "p_epi": g["epidermis"][2], "amp_epi": g["epidermis"][1],
                       "acro_derm": g["dermis"][0], "p_derm": g["dermis"][2], "amp_derm": g["dermis"][1]}, index=genes)
    gr = gw[(gw.p_epi < 0.01) & (gw.p_derm < 0.01)].copy()
    gr["shift_h"] = wrap(gr.acro_derm - gr.acro_epi)
    gr.rename_axis("gene").to_csv(out / "a67_genomewide_rhythmic.csv")
    gm, gsd = circ(gr.shift_h)
    v = dict(n_clock_rhythmic_both=int(rhy.sum()), clock_mean_shift_h=mean, clock_shift_lo=lo, clock_shift_hi=hi,
             clock_circ_sd_h=sd, n_genome_rhythmic_both=len(gr), n_rhythmic_epi=int((gw.p_epi < 0.01).sum()),
             n_rhythmic_derm=int((gw.p_derm < 0.01).sum()), genome_mean_shift_h=gm, genome_circ_sd_h=gsd,
             genome_median_shift_h=float(np.median(gr.shift_h)),
             genome_frac_epidermis_ahead=float(np.mean(gr.shift_h > 0)),
             H67_1="PASS" if (mean > 0 and lo > 0 and sd < 2) else "FAIL")
    pd.Series(v).to_csv(out / "a67_verdict.csv", header=False)
    print(d.round(2).to_string(), "\n", pd.Series(v).to_string(), flush=True)


def a76(s, m, out, n_perm, seed=76):
    y = m.loc[H.GENES, s.gsm].T * np.log(2.0)
    y.index = s.index
    L = H.PD.projector(y, H.WF, counts=False)
    beta = (y - y.mean(1).to_numpy()[:, None]).to_numpy() @ L.T
    s = s.assign(c=beta[:, 0] + 1j * beta[:, 1])
    idx = lambda genes: m.loc[[g for g in genes if g in m.index], s.gsm].mean(0).to_numpy()
    s["diff"], s["prolif"] = idx(H.DIFF) - idx(H.BASAL), idx(H.PROLIF)
    rows = []
    for subj, d in s.groupby("subject"):
        a = {}
        for layer, x in d.groupby("layer"):
            e = np.exp(1j * W * x.clock_h.to_numpy(float))
            X = np.column_stack([np.ones_like(e), e, np.conj(e)])
            a[layer] = np.linalg.lstsq(X, x.c.to_numpy(), rcond=None)[0][1]
        ep = d[d.layer == "epidermis"]
        rows.append(dict(subject=subj, offset_h=float(np.angle(a["epidermis"] / a["dermis"]) / W),
                         amp_ratio=float(abs(a["epidermis"] / a["dermis"])), diff=ep["diff"].mean(),
                         prolif=ep.prolif.mean()))
    d = pd.DataFrame(rows)
    rng = np.random.default_rng(seed)
    o = {}
    for k in ("diff", "prolif"):
        rho = spearmanr(d[k], d.offset_h).statistic
        null = np.array([spearmanr(rng.permutation(d[k].to_numpy()), d.offset_h).statistic for _ in range(n_perm)])
        o[f"rho_{k}"] = rho
        o[f"p_upper_{k}"] = float((1 + np.sum(null >= rho)) / (1 + n_perm))
        o[f"p_lower_{k}"] = float((1 + np.sum(null <= rho)) / (1 + n_perm))
    o["mean_offset_h"] = float(np.angle(np.mean(np.exp(1j * W * d.offset_h))) / W)
    o["H76_1"] = "PASS" if (o["rho_diff"] > 0 and o["p_upper_diff"] < 0.05) else "FAIL"
    pd.Series(o).to_csv(out / "a76_verdict.csv", header=False)
    d[["offset_h", "amp_ratio", "diff", "prolif"]].describe().to_csv(out / "a76_summary.csv")
    print(pd.Series(o).to_string(), flush=True)


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    s, m = H.timed_skin()
    sz = H.sizes(args, run, 10000, 1000)
    a67(s, m, out, sz["n_boot"])
    a76(s, m, out, sz["n_perm"])


if __name__ == "__main__":
    main()
