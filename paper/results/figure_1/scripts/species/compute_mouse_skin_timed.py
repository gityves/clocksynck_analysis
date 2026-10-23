"""A72: epidermal stem cells vs whole epidermis across timed mouse datasets (ZT known, LD 12:12).

Pre-registered as A72 (2026-10-18; docs/PLAN_internal_desync.md (archive); data rule and implementation note before scoring; run note: H72.1 passes). Design
(weaker than A71, stated there): datasets are compared, not populations within animals; confounded by lab, age, platform; supporting evidence only. Data
(`mouse_skin_timed` builder): stem cells GSE84510 / GSE84511 / GSE84521 (Affymetrix arrays), whole epidermis GSE190032 / GSE294336 / GSE294339; strict
sensitivity = young adults (GSE84510, GSE84511, GSE190032, GSE294336). Per dataset: clock-plane state beta of the 14 clock genes (GTEx waveform, mouse
orthologues; counts -> log fractions, arrays -> log2 x ln 2; one projector per dataset), complex state c = beta_1 + i beta_2; regression
c = m + a e^{iwt} + b e^{-iwt} on ZT; clock offset = arg(a) / w (h, GTEx frame minus ZT; larger = ahead); bootstrap of samples within ZT (1,000, seed 72).
H72.1: the circular mean offset of whole-epidermis datasets minus that of stem-cell datasets > 0 with a bootstrap 95 % interval > 0 (infeasible without one
population). Descriptive: per-gene acrophases (cosinor on ZT) per dataset. Post hoc (run note): dataset-level label permutation over the 6 datasets (20 labellings, one-sided; mean-difference and rank-sum statistics).

    python compute_mouse_skin_timed.py [--subsample N] [--n-boot B]      (--subsample: 100 bootstraps)
Outputs: paper/results/figure_1/output/species/a72_{primary,strict}_{datasets,verdict,gene_acrophases}.csv, a72_dataset_permutation.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a72_mouse_timed.py -> output/a72_*
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import log_expr, projector, states  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "species"
W = 2 * np.pi / 24
WF = cs.gtex()
GENES = list(WF.genes)
PRIMARY = ["GSE84510", "GSE84511", "GSE84521", "GSE190032", "GSE294336", "GSE294339"]
STRICT = ["GSE84510", "GSE84511", "GSE190032", "GSE294336"]


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def cmean(h):
    return float(np.angle(np.mean(np.exp(1j * W * np.asarray(h)))) / W)


def dataset_state(d: pd.DataFrame) -> np.ndarray:
    counts = d.value_type.iloc[0] == "counts"
    x = d[GENES].astype(float)
    L = projector(x, WF, counts=counts)
    beta, _ = states(x, L, WF, depth=d.depth.to_numpy(float) if counts else None)
    return beta[:, 0] + 1j * beta[:, 1]


def offset(c, zt) -> float:
    e = np.exp(1j * W * zt)
    coef = np.linalg.lstsq(np.column_stack([np.ones_like(e), e, np.conj(e)]), c, rcond=None)[0]
    return float(np.angle(coef[1]) / W)


def gene_acro(d) -> pd.Series:
    zt = d.zt.to_numpy(float)
    X = np.column_stack([np.ones_like(zt), np.cos(W * zt), np.sin(W * zt)])
    counts = d.value_type.iloc[0] == "counts"
    y = log_expr(d[GENES].astype(float), d.depth.to_numpy(float) if counts else None).to_numpy()
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    return pd.Series(np.mod(np.arctan2(b[2], b[1]) / W, 24), index=GENES)


def run(T, names, tag, out, n_boot):
    rng = np.random.default_rng(72)
    rows, acro, boots = [], {}, {}
    for name in names:
        d = T[T.cohort == name].reset_index(drop=True)
        c, zt = dataset_state(d), d.zt.to_numpy(float)
        o = offset(c, zt)
        groups = [np.flatnonzero(zt == z) for z in np.unique(zt)]
        bs = np.array([offset(c[idx], zt[idx]) for idx in (np.concatenate([rng.choice(g, len(g)) for g in groups]) for _ in range(n_boot))])
        dd = wrap(bs - o)
        pop = d.condition.iloc[0]
        rows.append(dict(dataset=name, population=pop, n=len(d), n_zt=len(groups), offset_h=o, lo=o + np.quantile(dd, 0.025), hi=o + np.quantile(dd, 0.975)))
        acro[name], boots[name] = gene_acro(d), (pop, bs)
    r = pd.DataFrame(rows)
    pops = {p: [k for k, (q, _) in boots.items() if q == p] for p in ("stem", "whole_epidermis")}
    v = dict(n_stem=len(pops["stem"]), n_whole=len(pops["whole_epidermis"]))
    off = r.set_index("dataset").offset_h
    if all(pops.values()):
        diff = float(wrap(cmean(off[pops["whole_epidermis"]]) - cmean(off[pops["stem"]])))
        bd = np.array([wrap(cmean([boots[k][1][i] for k in pops["whole_epidermis"]]) - cmean([boots[k][1][i] for k in pops["stem"]]))
                       for i in range(n_boot)])
        dd = wrap(bd - diff)
        v.update(whole_minus_stem_h=diff, lo=diff + np.quantile(dd, 0.025), hi=diff + np.quantile(dd, 0.975))
        v["H72_1"] = "PASS" if (diff > 0 and v["lo"] > 0) else "FAIL"
        if tag == "primary":  # post hoc dataset-level permutation of the population labels (one-sided)
            labs = np.array([boots[k][0] for k in names])
            ids = range(len(names))
            nw = int((labs == "whole_epidermis").sum())
            sets = [(sel, [i for i in ids if i not in sel]) for sel in combinations(ids, nw)]
            names_ = lambda ix: [names[i] for i in ix]
            wsel = tuple(i for i in ids if labs[i] == "whole_epidermis")
            md = lambda w, st: float(wrap(cmean(off[names_(w)]) - cmean(off[names_(st)])))
            rk = off[names].rank().to_numpy()
            v["perm_p_mean_diff"] = float(np.mean([md(w, st) >= md(wsel, [i for i in ids if i not in wsel]) - 1e-12 for w, st in sets]))
            v["perm_p_rank_sum"] = float(np.mean([rk[list(w)].sum() >= rk[list(wsel)].sum() - 1e-12 for w, st in sets]))
    else:
        v["H72_1"] = "INFEASIBLE"
    r.to_csv(out / f"a72_{tag}_datasets.csv", index=False)
    pd.DataFrame(acro).to_csv(out / f"a72_{tag}_gene_acrophases.csv")
    pd.Series(v).to_csv(out / f"a72_{tag}_verdict.csv", header=False)
    print(tag, "\n", r.round(2).to_string(index=False), "\n", pd.Series(v).to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run_ = cli.setup(args)
    out = out_dir(FOLDER, run_["quick"])
    T, _ = data._builder("mouse_skin_timed").build()
    n_boot = args.n_boot or (100 if run_["quick"] else 1000)
    run(T, PRIMARY, "primary", out, n_boot)
    run(T, STRICT, "strict", out, n_boot)


if __name__ == "__main__":
    main()
