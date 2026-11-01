"""A181 H181.2: responsive-gene shared state in the A107b whole-organ (H107b.1) and stromal (H107b.2) designs (Cao 2020).

Genes: FOS FOSB JUN JUNB EGR1 EGR3 IER2 ATF3 NR4A1 DUSP1 FKBP5 TSC22D3 ZBTB16 DDIT4 KLF9 (from the 68-gene panel h5ad).
`Dataset.gate` needs the 14 clock genes and a clock-plane waveform, so the pairing null is computed here with clocksync's
internals: statistic = sum over compartment pairs of ||C_kl||^2 (cross-covariance of the centred log fractions across
fetuses, identity projector: a flat waveform has no plane), null = fetus pairing shuffled within each organ
(clocksync shuffle_pairing), 1000 permutations, one-sided add-one p. Holm over the two designs. Descriptive companion:
the same statistic on the 14 clock genes (identity projector); the registered companion is the A107b clock-plane p_pairing.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core.controls import shuffle_pairing  # noqa: E402
from clocksync._core.study import centred_log_fraction, cross_covariances  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import data  # noqa: E402
from results.figure_2.scripts.embryo.compute_fetal_atlases import CAO, MIN_DON  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

OUT = HERE.parents[1] / "output"
CAO_DIR = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/natexp/cao2020")
RESP = "FOS FOSB JUN JUNB EGR1 EGR3 IER2 ATF3 NR4A1 DUSP1 FKBP5 TSC22D3 ZBTB16 DDIT4 KLF9".split()
N_PERM, SEED = 1000, 1812


def pairing(counts, depth, res, comp, strat, seed):
    y = centred_log_fraction(counts, depth)
    stat = lambda r: sum(np.sum(C ** 2) for *_, C in cross_covariances(y, r, comp, strat, 5))  # noqa: E731
    obs = stat(res)
    rng = np.random.default_rng(seed)
    grp = np.char.add(np.char.add(comp.astype(str), "|"), strat.astype(str))
    null = np.array([stat(shuffle_pairing(res, grp, rng)) for _ in range(N_PERM)])
    return obs, (1 + np.sum(null >= obs)) / (1 + N_PERM), float(np.median(null))


def main():
    limit_threads(4)
    clock = list(cs.gtex().genes)
    rows = []
    for j, (test, (lev, comps, mc)) in enumerate(list(CAO.items())[:2]):
        pb, _ = data.load("fetal_cao", level=lev)
        el = S.eligible(pb, mc, MIN_DON, compartments=comps)
        a = ad.read_h5ad(CAO_DIR / f"pseudobulk_{lev}" / "pseudobulk.h5ad")
        key = a.obs.donor_id.astype(str) + "|" + (a.obs.organ.astype(str) if lev == "organ" else a.obs.cell_type.astype(str))
        names = [canonical(g) for g in a.var_names]
        X = pd.DataFrame(np.asarray(a.X, float), columns=names, index=key.to_numpy())
        k = (el.donor + "|" + el.compartment).to_numpy()
        assert X.index.is_unique and set(k) <= set(X.index)
        X = X.loc[k]
        depth = el.depth.to_numpy(float)
        assert np.allclose(depth, a.obs.depth.to_numpy()[pd.Index(key).get_indexer(k)])
        res, comp, strat = el.donor.to_numpy(str), el.compartment.to_numpy(str), el.stratum.to_numpy(str)
        for label, genes in (("responsive", RESP), ("clock_identity", clock)):
            cols = [canonical(g) for g in genes]
            obs, p, nm = pairing(X[cols].to_numpy(), depth, res, comp, strat, SEED + j)
            rows.append(dict(test=test, gene_set=label, n_genes=len(cols), donors=el.donor.nunique(),
                             compartments=";".join(sorted(set(comp))), stat=obs, null_median=nm, p_pairing=p))
    t = pd.DataFrame(rows)
    r = t.gene_set == "responsive"
    t["p_holm"] = np.nan
    t.loc[r, "p_holm"] = cs.stats.holm(t.loc[r, "p_pairing"].to_numpy())
    t.loc[r, "passed_holm"] = t.loc[r, "p_holm"] < 0.05
    OUT.mkdir(exist_ok=True)
    t.to_csv(OUT / "h181_2_responsive.csv", index=False)
    print(t.drop(columns="compartments").round(4).to_string(index=False))


if __name__ == "__main__":
    main()
