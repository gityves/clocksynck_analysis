"""A17 step 7: cross-organ confirmation of the anchor rule (descriptive).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17 step 7 "cross-organ confirmation (descriptive)"; run note 2026-10-12). No hypothesis id.
  tabula  Tabula Sapiens: pair excess coherence e_kl of compartments from different organs of the same donors, anchor x anchor pairs (learned
          and fixed rule, frozen table) against the others; interval from compartment bootstraps (a pair is kept when both compartments are
          drawn; the registered donor bootstrap is not possible because the pair table holds pairs, not donors).
  gtex    GTEx bulk, fast deaths: Spearman correlation over the 8 tissues of the Eraslan snRNA atlas between the tissue's stromal-vascular
          nucleus fraction (fibroblast, myofibroblast, endothelial, pericyte; also with smooth-muscle myocytes) and (a) its agreement with
          recorded time (`phase_dep_sd`, lower = closer) and (b) its untimed coupling to the other tissues, also partial given the tissue's
          clock amplitude. Aggregates over tissues only; composition from other donors (Eraslan, ventilator cases) than the GTEx bulk.
Eraslan per-donor phases against time of death are not run (12 timed donors, all Hardy 0: a ventilator check only; as legacy).
Reads the pair table of `host_identity` (`a17_coherence_pairs.csv`) and `organs/tissues.csv` (the GTEx tissue table), both
documented outputs; run those first.

    python compute_step7_confirmations.py [--subsample N] [--n-boot B]
Outputs: paper/results/figure_2/output/atlas_fits/a17_step7_tabula.csv, a17_step7_gtex_tissues.csv, a17_step7_gtex.csv
Legacy: legacy/desync/common/step7_confirm.py -> common_output/a17_step7_{tabula,gtex}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

from _shared import FOLDER_4, FOLDER_5, FOLDER_S04, frozen  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

ERASLAN = [DATA.parent / "CoPhaser_data" / "data" / "gtex" / "snrna" / "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad",
           DATA.parent / "data" / "gtex" / "snrna" / "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"]


sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_out  # noqa: E402
def tabula(quick: bool, n_boot: int) -> pd.DataFrame:
    a = frozen()
    a = a[a.target == "tabula_sapiens"].set_index("compartment")
    sub = "subsample" if quick else ""
    f = group_out(FOLDER_S04) / sub / "a17_coherence_pairs.csv"
    p = pd.read_csv(f)
    p = p[p.name == "tabula_sapiens"]
    p = p[p.k.str.split(":").str[0] != p.l.str.split(":").str[0]].copy()
    comps = np.unique(np.r_[p.k, p.l])
    rows = []
    for rule in ("anchor", "anchor_fixed"):
        A = a[rule].astype(bool)
        p["aa"] = [bool(A.get(k, False) and A.get(l, False)) for k, l in zip(p.k, p.l)]
        stat = lambda q: q[q.aa].excess.mean() - q[~q.aa].excess.mean()
        rng = np.random.default_rng(0)
        reps = []
        for _ in range(n_boot):
            s = rng.choice(comps, len(comps))
            q = p[p.k.isin(s) & p.l.isin(s)]
            if q.aa.any() and (~q.aa).any():
                reps.append(stat(q))
        rows.append(dict(rule=rule, pairs=len(p), anchor_pairs=int(p.aa.sum()), excess_anchor=p[p.aa].excess.mean(),
                         excess_other=p[~p.aa].excess.mean(), diff=stat(p), lo=np.quantile(reps, 0.025) if reps else np.nan,
                         hi=np.quantile(reps, 0.975) if reps else np.nan, n_boot=len(reps)))
    return pd.DataFrame(rows)


def partial(x, y, z) -> float:
    """Partial Spearman correlation of x and y given z (the tissue's clock amplitude |a_k|)."""
    ok = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    Z = np.c_[np.ones(ok.sum()), rankdata(z[ok])]
    res = lambda v: rankdata(v[ok]) - Z @ np.linalg.lstsq(Z, rankdata(v[ok]), rcond=None)[0]
    return float(np.corrcoef(res(x), res(y))[0, 1])


def gtex(quick: bool):
    path = next(p for p in ERASLAN if p.exists())
    o = ad.read_h5ad(path, backed="r").obs[["Tissue Site Detail", "Broad cell type"]].astype(str)
    sv = o["Broad cell type"].str.match(r"Fibroblast|Myofibroblast|Endothelial|Pericyte")
    sm = o["Broad cell type"].str.contains("smooth muscle")
    tis = o["Tissue Site Detail"]
    f = pd.DataFrame(dict(sv=sv.groupby(tis).mean(), sv_smooth=(sv | sm).groupby(tis).mean(), nuclei=tis.value_counts()))
    t = pd.read_csv(group_out("organs") / ("subsample" if quick else "") / "tissues.csv")
    t = t[t.cohort == "fast"].set_index("tissue")  # the organs table holds both cohorts; step 7 reads fast deaths
    f = f.join(t[["donors", "phase_dep_sd", "coupling", "abs_amp"]], how="left")
    rows = []
    for c in ("sv", "sv_smooth"):
        for y in ("phase_dep_sd", "coupling"):
            d = f[[c, y]].dropna()
            r = spearmanr(d[c], d[y])
            rows.append(dict(fraction=c, measure=y, tissues=len(d), spearman=r.statistic, p=r.pvalue,
                             partial_given_amp=partial(f[c].to_numpy(float), f[y].to_numpy(float), f.abs_amp.to_numpy(float))))
    return f, pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER_5, run["quick"])
    pd.set_option("display.width", 250)
    t = tabula(run["quick"], args.n_boot or (200 if run["quick"] else 2000))
    t.to_csv(out / "a17_step7_tabula.csv", index=False)
    print(t.round(3).to_string(index=False))
    f, r = gtex(run["quick"])
    f.to_csv(out / "a17_step7_gtex_tissues.csv")
    r.to_csv(out / "a17_step7_gtex.csv", index=False)
    print(r.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
