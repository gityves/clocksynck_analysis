"""A276, registered diagnostics after V276 failed: why is the clock's null |r| above its companions'?

Per atlas, compartment and gene set (clock, 200 matched, 200 responsive; states as `calibrate_a276.load`), on the states
centred within stratum: anisotropy = share of the plane variance on its first axis (0.5 isotropic, 1 one-dimensional),
kurtosis of |z|^2 / mean |z|^2, and the effective number of genes carrying the projection (1 / sum of squared normalised
loadings). Per pair (`output/pair_null.csv`): the mean null |r| over the 500 draws against the pair's anisotropy, for the
clock and the companion means. Expectation for Gaussian states: E|r|^2 = 1 / df whatever the anisotropy, E|r| larger for
isotropic states (0.886 / sqrt(df) against 0.798 / sqrt(df) for one-dimensional ones).

    python diagnose_a276.py
Outputs (../output/): diag_compartments.csv, diag_pairs.csv, diag_summary.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import calibrate_a276 as C  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402


def shape(z):
    """(anisotropy, kurtosis ratio) of centred complex states."""
    X = np.stack([z.real, z.imag])
    ev = np.linalg.eigvalsh(np.cov(X))
    m2 = np.abs(z) ** 2
    return float(ev[-1] / ev.sum()), float(np.mean(m2 ** 2) / np.mean(m2) ** 2)


def n_eff(pb, w):
    L = P.projector(pb[list(w.genes)], w, counts=True)
    q = (L ** 2).sum(0)
    q = q / q.sum()
    return float(1 / np.sum(q ** 2))


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--reuse", action="store_true", help="summarise an existing diag_compartments.csv")
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(C.AMENDMENT)
    C.A150.ATLASES.update(C.A274.ALL)
    rows = []
    for study, names in () if args.reuse else (("A272", C.A272_ATLASES), ("A274", C.A274_ATLASES)):
        for n in names:
            A = C.A150.load_atlas(n, None, 0)
            pb = A["pb"].reset_index(drop=True)
            wfs = [("clock", cs.gtex())]
            for kind in C.KINDS:
                ws = (C.A272.pseudo_sets(A, C.A272_ATLASES.index(n), 200, kind) if study == "A272"
                      else C.A274.pseudo_sets(A, list(C.A274.ALL).index(n), 200, kind))
                wfs += [(kind, w) for w in ws]
            for j, (kind, w) in enumerate(wfs):
                z = C.A272.atlas_states(pb[C.A272.META + list(w.genes)], w)[0]
                ne = n_eff(pb, w)
                for comp, g in pb.groupby("compartment"):
                    idx = g.index.to_numpy()
                    s = pd.factorize(g.stratum)[0]
                    zc = C.A272._centre(z[idx], s)
                    a, k = shape(zc)
                    rows.append(dict(study=study, atlas=n, compartment=comp, set=j, kind=kind, n_eff=ne, anisotropy=a,
                                     kurtosis=k, rows=len(idx)))
            C.log(study, n)
    if args.reuse:
        D = pd.read_csv(out / "diag_compartments.csv")
    else:
        D = pd.DataFrame(rows)
        D.to_csv(out / "diag_compartments.csv", index=False)

    comp = D.groupby(["study", "atlas", "compartment", "kind"]).agg(anisotropy=("anisotropy", "mean"),
                                                                    kurtosis=("kurtosis", "mean"),
                                                                    n_eff=("n_eff", "mean")).reset_index()
    wide = comp.pivot_table(index=["study", "atlas", "compartment"], columns="kind",
                            values=["anisotropy", "kurtosis", "n_eff"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    pn = pd.read_csv(out / "pair_null.csv")
    for side in ("a", "b"):
        pn = pn.merge(wide.rename(columns={c: f"{c}_{side}" for c in wide.columns if c not in ("study", "atlas")}
                                  | {"compartment": side}), on=["study", "atlas", side], how="left")
    for kind in ("clock", "matched", "responsive"):
        pn[f"iso_{kind}"] = 1 - (pn[f"anisotropy_{kind}_a"] + pn[f"anisotropy_{kind}_b"]) / 2
    pn.to_csv(out / "diag_pairs.csv", index=False)
    summ = []
    for kind in ("matched", "responsive"):
        d_r = pn.clock_null - pn[f"{kind}_null"]
        d_iso = pn.iso_clock - pn[f"iso_{kind}"]
        summ.append(dict(family=kind, pairs=len(pn), null_r_clock_median=float(pn.clock_null.median()),
                         null_r_family_median=float(pn[f"{kind}_null"].median()),
                         share_pairs_clock_above=float((d_r > 0).mean()),
                         anisotropy_clock_median=float(comp[comp.kind == "clock"].anisotropy.median()),
                         anisotropy_family_median=float(comp[comp.kind == kind].anisotropy.median()),
                         kurtosis_clock_median=float(comp[comp.kind == "clock"]["kurtosis"].median()),
                         kurtosis_family_median=float(comp[comp.kind == kind]["kurtosis"].median()),
                         n_eff_clock_median=float(comp[comp.kind == "clock"].n_eff.median()),
                         n_eff_family_median=float(comp[comp.kind == kind].n_eff.median()),
                         rho_dnull_diso=float(spearmanr(d_r, d_iso).statistic),
                         rho_dnull_dkurt=float(spearmanr(d_r, (pn.kurtosis_clock_a + pn.kurtosis_clock_b)
                                                         - (pn[f"kurtosis_{kind}_a"] + pn[f"kurtosis_{kind}_b"])).statistic)))
    S = pd.DataFrame(summ)
    S.to_csv(out / "diag_summary.csv", index=False)
    print(S.T.to_string())


if __name__ == "__main__":
    main()
