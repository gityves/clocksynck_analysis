"""A68: are differentiated-keratinocyte-rich skin samples phase-advanced? GTEx skin, fast deaths, time of death.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A68; scores from A58 / A61). GTEx v10 skin, sun exposed (lower leg) and not sun
exposed (suprapubic), fast deaths (Hardy 1-2) with a time of death t (controlled access: aggregates only). Per sample: clock-plane
state beta (14-gene GTEx projector per tissue, Poisson-corrected, `clocksync.paired`), z = (beta_1 + i beta_2) exp(-i w t).
Scores from log(1 + CPM) (genes z-scored within tissue, set means, score z-scored within tissue): x_diff = differentiated minus
basal keratinocyte set (A61); x_epi = keratinocyte (A58) minus fibroblast + endothelial (A58) sets. Within each tissue u = z /
mean(z); u = 1 + c x by least squares (x centred), pooled over tissues with one c; k = Im(c) / w (h per SD; > 0: advanced at a
higher score), Re(c) = relative amplitude slope. H68.1 (x_diff) and H68.2 (x_epi): k > 0, one-sided p < 0.05 against x permuted
within tissue (n_perm 2,000). Secondary: per tissue; scores residualised on sex and age bracket (`_adj`).

    python compute_gtex_skin.py [--subsample N] [--n-perm 2000]
Outputs (paper/results/figure_5/output/psoriasis_details/): a68_gtex_skin.csv, a68_verdict.csv (aggregates only)
Legacy: legacy/desync/validation/paired_damping/scripts/a68_gtex_skin.py
"""

from __future__ import annotations

import sys
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("psoriasis")))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _cues as Q  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402

TISSUES = ["Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"]
W, SEED = H.W, 68
BASAL = "KRT5 KRT14 KRT15 COL17A1 TP63 ITGA6 ITGB4 CXCL14".split()
DIFF = "KRT1 KRT10 KRT2 DSG1 DSC1 SBSN KRTDAP LORICRIN FLG IVL CNFN LCE1B LCE3D SPRR1A SPRR1B SPRR2A".split()
KER = ("KRT5 KRT14 KRT1 KRT10 KRT15 KRT17 KRT6A KRT16 KRT2 SBSN LORICRIN SFN PERP DSP JUP KRT6B DMKN LCE3D LCE1B CNFN "
       "CRCT1 SERPINB3 SPRR2A S100A2 KLK7").split()
DERM = "DCN LUM COL1A1 COL1A2 COL3A1 PDGFRA CCL19 APOE PECAM1 VWF CLDN5 PLVAP ACKR1 SELE LYVE1 PROX1 CCL21".split()
ALIAS = {"LORICRIN": "LOR"}


def score(lc: pd.DataFrame, pos, neg) -> np.ndarray:
    f = lambda gs: [ALIAS.get(g, g) if g not in lc else g for g in gs]
    z = (lc - lc.mean()) / lc.std().replace(0, 1)
    s = z[[g for g in f(pos) if g in z]].mean(1) - z[[g for g in f(neg) if g in z]].mean(1)
    return ((s - s.mean()) / s.std()).to_numpy()


def residual(x, ob) -> np.ndarray:
    Dm = pd.get_dummies(ob[["sex", "age_bracket"]].astype(str), drop_first=True).to_numpy(float)
    Dm = np.column_stack([np.ones(len(x)), Dm])
    r = x - Dm @ np.linalg.lstsq(Dm, x, rcond=None)[0]
    return r / r.std()


def slope(us, xs) -> complex:
    return sum(np.sum((x - x.mean()) * (u - 1)) for u, x in zip(us, xs)) / sum(np.sum((x - x.mean()) ** 2) for x in xs)


def test(us, xs, rng, label, n_perm):
    c = slope(us, xs)
    k = float(np.imag(c) / W)
    null = np.array([np.imag(slope(us, [rng.permutation(x) for x in xs])) / W for _ in range(n_perm)])
    return dict(set=label, n=int(sum(len(x) for x in xs)), k_h_per_sd=k, amp_slope=float(np.real(c)),
                p_upper=float((1 + np.sum(null >= k)) / (1 + n_perm)))


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    n_perm = H.sizes(args, run, 2000)["n_perm"]
    rng = np.random.default_rng(SEED)
    score_genes = sorted(set(DIFF + BASAL + KER + DERM))
    data = {}
    for t in TISSUES:
        d = Q.gtex_rows(t, (1, 2), None, H.GENES + score_genes, depth="sum")
        if args.subsample:
            d = d.sample(min(args.subsample, len(d)), random_state=args.seed).reset_index(drop=True)
        cl = d[H.GENES]
        L = H.PD.projector(cl, H.WF, counts=True)
        beta, _ = H.PD.states(cl, L, H.WF, depth=d.depth.to_numpy(float))
        z = (beta[:, 0] + 1j * beta[:, 1]) * np.exp(-1j * W * d.tod.to_numpy(float))
        u = z / z.mean()
        lc = np.log1p(d[[g for g in score_genes if g in d]].div(d.depth, axis=0) * 1e6)
        sc = {"diff": score(lc, DIFF, BASAL), "epi": score(lc, KER, DERM)}
        sc.update({f"{k}_adj": residual(v, d) for k, v in list(sc.items())})
        data[t] = (u, sc, float(np.angle(z.mean()) / W) % 24, float(abs(z.mean()) / np.mean(abs(z))))
        print(t, len(u), "corr(diff, epi) %.2f" % np.corrcoef(sc["diff"], sc["epi"])[0, 1], flush=True)
    rows = []
    for key in ("diff", "epi", "diff_adj", "epi_adj"):
        rows.append(dict(score=key, **test([d[0] for d in data.values()], [d[1][key] for d in data.values()], rng,
                                           "pooled", n_perm)))
        for t, d in data.items():
            rows.append(dict(score=key, **test([d[0]], [d[1][key]], rng, t, n_perm)))
    res = pd.DataFrame(rows)
    pr = res[res.set == "pooled"].set_index("score")
    v = dict(H68_1="PASS" if (pr.k_h_per_sd["diff"] > 0 and pr.p_upper["diff"] < 0.05) else "FAIL",
             H68_2="PASS" if (pr.k_h_per_sd["epi"] > 0 and pr.p_upper["epi"] < 0.05) else "FAIL",
             **{f"timing_consistency_{t}": d[3] for t, d in data.items()})
    res.to_csv(out / "a68_gtex_skin.csv", index=False)
    pd.Series(v).to_csv(out / "a68_verdict.csv", header=False)
    print(res.round(3).to_string(), "\n", pd.Series(v).to_string())


if __name__ == "__main__":
    main()
