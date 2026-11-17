"""A293 Part 3 (registered in `paper/docs/registered_experiments.md`): does each cell class carry the same gain relative to
endothelium in every organ?

A17 healthy atlas fits (rule "learned", status ok), body atlases. Per atlas and family: mean log_amp_rel of its
compartments minus that of the atlas's endothelial compartments; SE from the 95 % bootstrap intervals, combined in
quadrature. Per family over atlases: random-effects pooled value, tau, Cochran's Q (`clocksync.stats.meta`).
H293.4: tau <= 0.35 in every family with >= 4 atlases.

    python part3_relative_gain.py
Outputs (`../output/part3/`): atlases.csv, families.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2  # noqa: E402

from clocksync import stats as S  # noqa: E402
from clocksync.anchors import family  # noqa: E402

FITS = HERE.parents[3] / "results" / "figure_2" / "output" / "atlas_fits" / "a17_healthy_fits.csv"
OUT = HERE.parents[1] / "output" / "part3"
EXCLUDE = {"blood_aida", "blood_onek1k", "organ_bone_marrow", "ia_mesenteric_ln", "ia_thoracic_ln", "ia_spleen",
           "psychad_hbcc", "psychad_mssm", "seaad_dlpfc", "seaad_mtg", "pd_control", "retina", "tabula_sapiens"}
FAMILIES = ("fibroblast_stromal", "mural", "epithelial", "parenchyma", "myeloid")
TAU_MAX, MIN_ATLASES = 0.35, 4


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f = pd.read_csv(FITS)
    f = f[(f.rule == "learned") & (f.status == "ok") & ~f.target.isin(EXCLUDE)].copy()
    f["family"] = f.compartment.map(family)
    f["se"] = (f.log_amp_hi - f.log_amp_lo) / 3.92
    rows = []
    for atlas, g in f.groupby("target"):
        ec = g[g.family == "endothelial"]
        if not len(ec) or ec.log_amp_rel.isna().all():
            continue
        e_mu, e_se = ec.log_amp_rel.mean(), np.sqrt(np.nansum(ec.se ** 2)) / len(ec)
        for fam in FAMILIES:
            x = g[g.family == fam]
            if not len(x) or x.log_amp_rel.isna().all():
                continue
            se = np.sqrt(np.nansum(x.se ** 2)) / len(x)
            rows.append(dict(atlas=atlas, family=fam, compartments=len(x), endothelial=len(ec),
                             g=float(x.log_amp_rel.mean() - e_mu), se=float(np.sqrt(se ** 2 + e_se ** 2))))
    at = pd.DataFrame(rows)
    at.to_csv(OUT / "atlases.csv", index=False)
    fam_rows = []
    for fam, g in at.groupby("family"):
        g = g[np.isfinite(g.se) & (g.se > 0)]
        if len(g) < 2:
            continue
        m = S.meta(g.g.to_numpy(), g.se.to_numpy())
        q = float(np.sum(((g.g - m["fixed"]) / g.se) ** 2)) if "fixed" in m else np.nan
        fam_rows.append(dict(family=fam, atlases=len(g), pooled=m["random"], se_pooled=m["se_random"],
                             tau=float(np.sqrt(max(m["tau2"], 0))), Q=q, Q_p=float(chi2.sf(q, len(g) - 1)),
                             min_g=float(g.g.min()), max_g=float(g.g.max()), tested=len(g) >= MIN_ATLASES))
    fm = pd.DataFrame(fam_rows)
    fm["within_bound"] = fm.tau <= TAU_MAX
    fm.to_csv(OUT / "families.csv", index=False)
    print(at.round(3).to_string(index=False))
    print(fm.round(3).to_string(index=False))
    t = fm[fm.tested]
    print("H293.4", "PASS" if len(t) and t.within_bound.all() else "FAIL", f"({len(t)} families tested)")


if __name__ == "__main__":
    main()
