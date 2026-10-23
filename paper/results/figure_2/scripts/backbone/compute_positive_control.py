"""A18 part 2 positive control (V10): the per-cell split on mouse dermis with known ZT.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A18, part 2, "Positive control"): no public dataset has the same cells deep-
sequenced both synchronised and unsynchronised; the control is the mouse dermis ZT time course (GSE223109, 18 samples, ZT 2-22
every 4 h) and mixtures of its ZT groups. Expected readings: fibroblasts above immune cells in R; R falling with mixing.
Estimand per cell type: D (within-sample dispersion of each (sample, cell type) group of >= 100 cells, `cells.vonmises.group_stats`,
minus the mean of 10 acrophase-permuted planes, cell-weighted), rho^2 = s / (1 - R_W^2) with s the between-sample variance of the
cell type's pseudobulk clock coefficients minus the mean measurement variance and R_W = |mean exp(i omega ZT)|; lambda, R from
lambda R = rho and the von Mises model. Mixtures: two samples Delta h apart (equal cell numbers), D of the mixed group with
rho |cos(omega Delta / 2)|, against the predicted R |cos(omega Delta / 2)|. Every row also after binomial thinning of the clock counts
to 2.4 clock UMI per cell (human T/NK depth). Null: acrophases permuted across genes. No bootstrap.

    python compute_positive_control.py [--subsample N]       (N: cells per group, quick check)
Outputs: paper/results/figure_2/output/backbone/a18_positive_control.csv, a18_positive_control_groups.csv
Legacy: legacy/desync/common/positive_control.py (+ `load` of legacy/scripts/v10_cell_dispersion.py, builder data/cells_mouse_skin).
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.cells import vonmises as VM  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

from compute_cancer_split import split_lr  # noqa: E402  (local: not in clocksync)

FOLDER = "backbone"
OMEGA = 2 * np.pi / 24
MIN_GROUP, TARGET, N_PERM = 100, 2.4, 10


def coherent_variance(by, counts, G, ct, thin, genes) -> float:
    """s = mean |beta - mean beta|^2 - mean tr V over the sample pseudobulks of one cell type (clock-plane coefficients
    of `Dataset.study`; local: reads two attributes of the study object)."""
    pb = []
    for s, idx in by.items():
        Cc, N = counts(idx, thin)
        pb.append(dict(donor=s, compartment=ct, stratum="", depth=N.sum(), n_cells=len(idx), **dict(zip(genes, Cc.sum(0)))))
    pb = pd.DataFrame(pb)
    # the Dataset needs two compartments; the second is a copy (projection is per compartment, so ct is unaffected)
    ds = cs.Dataset(pd.concat([pb, pb.assign(compartment=ct + "_copy")], ignore_index=True), anchors=[ct, ct + "_copy"],
                    min_donors=3)
    st = ds.study
    i = st.compartment == ct
    b = st.beta[i] - st.beta[i].mean(0)
    return float(np.mean(np.sum(b ** 2, 1)) - np.mean(np.trace(st.V[i], axis1=1, axis2=2)))


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    rng = np.random.default_rng(0)
    cells = data._builder("cells_mouse_skin").load()
    wf = cs.gtex()
    genes = list(wf.genes)
    C0, N0 = cells[genes].to_numpy(float), cells.umi.to_numpy(float)
    obs = cells[["sample", "zt", "ct"]]
    c = wf.amp * np.exp(-1j * wf.psi)
    LA = VM._split_projectors(c, 30, rng)  # local: private constructor of the split projectors
    LAp = [VM._split_projectors(p.amp * np.exp(-1j * p.psi), 30, rng) for p in (wf.permuted(rng) for _ in range(N_PERM))]

    def counts(idx, thin):
        Cc, N = C0[idx], N0[idx]
        if thin:
            f = min(1.0, TARGET / (Cc.sum() / len(idx)))
            Cc, N = rng.binomial(Cc.astype(np.int64), f).astype(float), N * f
        return Cc, N

    def D_of(idx, thin):
        Cc, N = counts(idx, thin)
        st = VM.group_stats(Cc, N, LA, rng=rng)
        null = np.mean([VM.group_stats(Cc, N, L, rng=rng)["D"] for L in LAp])
        return st["D"] - null, Cc.sum() / len(idx)

    rows, summ = [], []
    for thin in (False, True):
        for ct, sub in obs.groupby("ct"):
            by = {s: g.index.to_numpy() for s, g in sub.groupby("sample") if len(g) >= MIN_GROUP}
            if len(by) < 6:
                continue
            if run["quick"]:
                by = {s: v[:args.subsample or 150] if len(v) > (args.subsample or 150) else v for s, v in by.items()}
            zt = sub.groupby("sample").zt.first()
            s2 = coherent_variance(by, counts, None, ct, thin, genes)
            rw = float(abs(np.mean(np.exp(1j * OMEGA * zt[list(by)].to_numpy()))))
            rho = np.sqrt(s2 / (1 - rw ** 2)) if s2 > 0 else np.nan
            Ds = {s: D_of(idx, thin) for s, idx in by.items()}
            w = np.array([len(by[s]) for s in by], float)
            D = float(np.sum(w * np.array([Ds[s][0] for s in by])) / w.sum())
            lam, R = split_lr(D, rho, c, LA)
            umi = float(np.median([Ds[s][1] for s in by]))
            summ.append(dict(thin=thin, ct=ct, kind="real", delta_h=0.0, groups=len(by), clock_umi=umi, D=D, rho=rho, rw=rw,
                             lam=lam, R=R, R_pred=np.nan))
            rows += [dict(thin=thin, ct=ct, kind="real", sample=s, delta_h=0.0, D=Ds[s][0], cells=len(by[s])) for s in by]
            mix = {}
            for s1, s2_ in itertools.combinations(sorted(by), 2):
                dh = (zt[s2_] - zt[s1]) % 24
                dh = min(dh, 24 - dh)
                if dh == 0:
                    continue
                n = min(len(by[s1]), len(by[s2_]))
                idx = np.concatenate([rng.choice(by[x], n, replace=False) for x in (s1, s2_)])
                Dm = D_of(idx, thin)[0]
                mix.setdefault(dh, []).append((Dm, 2 * n))
                rows.append(dict(thin=thin, ct=ct, kind="mixed", sample=f"{s1}+{s2_}", delta_h=dh, D=Dm, cells=2 * n))
            for dh, v in sorted(mix.items()):
                Dm = float(np.average([x[0] for x in v], weights=[x[1] for x in v]))
                lm, Rm = split_lr(Dm, rho * abs(np.cos(OMEGA * dh / 2)), c, LA)
                summ.append(dict(thin=thin, ct=ct, kind="mixed", delta_h=dh, groups=len(v), clock_umi=umi, D=Dm,
                                 rho=rho * abs(np.cos(OMEGA * dh / 2)), rw=rw, lam=lm, R=Rm,
                                 R_pred=R * abs(np.cos(OMEGA * dh / 2)) if np.isfinite(R) else np.nan))
    pd.DataFrame(rows).to_csv(out / "a18_positive_control_groups.csv", index=False)
    s = pd.DataFrame(summ)
    s.to_csv(out / "a18_positive_control.csv", index=False)
    print(s.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
