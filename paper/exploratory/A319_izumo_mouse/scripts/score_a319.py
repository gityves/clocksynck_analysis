"""A319 (registered): forebrain Bmal1 KO (Izumo et al. 2014 eLife): rearranged or dissolved organ clocks?

From the published source data (organ mean phase and circular variance across mice; per-mouse circular variance across
the six organs): order variance V_order, order share V_order / mean V_obs, an independence null (von Mises per organ),
the excess within-mouse coherence E and Rayleigh tests per organ.

    python score_a319.py
Outputs (`../output/`): conditions.csv (one row per condition x genotype), organs.csv, tests.csv.
"""
import os
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
OUT = HERE.parents[1] / "output"
SRC = Path(os.environ.get("SCCC_DATA", "/mnt/max8tb/data/CoPhaser_data/cancer_circadian")) / "izumo2014_elife"
ORGANS = ["Pituitary", "Liver", "Kidney", "Heart", "Lung", "Spleen"]
N_SIM, SEED = 20000, 319


def rows(fig):
    x = zipfile.ZipFile(SRC / f"elife-04617-fig{fig}-data1-v2.docx").read("word/document.xml").decode()
    cell = lambda tc: "".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", tc)).strip()
    return [[cell(tc) for tc in re.findall(r"<w:tc>.*?</w:tc>", tr, re.S)] for tr in re.findall(r"<w:tr[ >].*?</w:tr>", x, re.S)]


def num(s):
    m = re.match(r"\s*([-\d.]+)", s)
    return float(m.group(1)) if m else np.nan


def organ_tables():
    """{(condition, genotype): DataFrame(organ, mu_deg, V, n)} from tables (a) of figures 5 and 8."""
    out = {}
    r5, r8 = rows(5), rows(8)
    for i, r in enumerate(r5):
        if r[0] in ("LD", "DD") and r[1] == "Pituitary" and len(out) < 4:
            for j, g in ((i + 1, "WT"), (i + 4, "BKO")):
                out[(r[0], g)] = (r5[j][2:8], r5[j + 2][2:8])
    lab = {"Fx/Fxin ZT": ("FR_ZT", "WT"), "Fx/Fxin CT": ("FR", "WT"), "BKO": ("FR", "BKO")}
    for i, r in enumerate(r8):
        if r[0] in lab and len(r) > 2 and "Mean Phase" in r[1] and lab[r[0]] not in out:
            out[lab[r[0]]] = (r[2:8], r8[i + 2][2:8])
    tabs = {}
    for k, (mu, v) in out.items():
        tabs[k] = pd.DataFrame(dict(organ=ORGANS, mu_deg=[num(s) for s in mu], V=[num(s) for s in v],
                                    n=[int(re.search(r"n\s*=\s*(\d+)", s).group(1)) for s in mu]))
    return tabs


def mouse_tables():
    """{(condition, genotype): per-mouse circular variances across organs} from tables (b), (c)."""
    out = {}
    for fig, conds in ((5, ("LD", "DD")), (8, ("FR",))):
        rr = rows(fig)
        for i, r in enumerate(rr):
            if r and r[0] in conds and len(r) == 8 and r[1].startswith("Mean Vector"):
                wt, ko = [], []
                for q in rr[i + 1:]:
                    if not q or not re.match(r"WT\d+", q[0]):
                        break
                    wt.append(num(q[3]))
                    ko.append(num(q[7]))
                out[(r[0], "WT")] = np.array([v for v in wt if np.isfinite(v)])
                out[(r[0], "BKO")] = np.array([v for v in ko if np.isfinite(v)])
    return out


def a1inv(R):
    """Inverse of A1 (Best & Fisher 1981): von Mises concentration from mean resultant length."""
    R = np.asarray(R, float)
    return np.where(R < 0.53, 2 * R + R ** 3 + 5 * R ** 5 / 6,
                    np.where(R < 0.85, -0.4 + 1.39 * R + 0.43 / (1 - R), 1 / np.maximum(R ** 3 - 4 * R ** 2 + 3 * R, 1e-9)))


def rayleigh_p(R, n):
    Rn = R * n
    return float(np.exp(np.sqrt(1 + 4 * n + 4 * (n ** 2 - Rn ** 2)) - (1 + 2 * n)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tabs, mice = organ_tables(), mouse_tables()
    rng = np.random.default_rng(SEED)
    cond_rows, organ_rows = [], []
    for (cond, g), t in tabs.items():
        mu = np.deg2rad(t.mu_deg.to_numpy())
        R = 1 - t.V.to_numpy()
        for o, r_, n_, m_ in zip(t.organ, R, t.n, t.mu_deg):
            organ_rows.append(dict(condition=cond, genotype=g, organ=o, mu_deg=m_, R=r_, n=n_, rayleigh_p=rayleigh_p(r_, n_)))
        row = dict(condition=cond, genotype=g, V_order=float(1 - abs(np.mean(np.exp(1j * mu)))),
                   organs_rayleigh_05=int(sum(rayleigh_p(r_, n_) < 0.05 for r_, n_ in zip(R, t.n))))
        key = ("FR", g) if cond == "FR_ZT" else (cond, g)
        if key in mice and cond != "FR_ZT":
            v = mice[key]
            th = rng.vonmises(mu[None, None, :], a1inv(R)[None, None, :], size=(N_SIM, len(v), len(mu)))
            vsim = (1 - np.abs(np.exp(1j * th).mean(2))).mean(1)
            row.update(n_mice=len(v), V_obs_mean=float(v.mean()), V_indep_mean=float(vsim.mean()),
                       E=float(vsim.mean() - v.mean()), p_excess=float((1 + np.sum(vsim <= v.mean())) / (1 + N_SIM)),
                       order_share=row["V_order"] / float(v.mean()))
        cond_rows.append(row)
    C, O = pd.DataFrame(cond_rows), pd.DataFrame(organ_rows)
    # descriptive: BKO vs WT organ means per condition
    for cond in ("LD", "DD", "FR"):
        a, b = tabs[(cond, "BKO")], tabs[(cond, "WT")]
        d = np.rad2deg(np.angle(np.exp(1j * np.deg2rad(a.mu_deg.to_numpy() - b.mu_deg.to_numpy()))))
        for o, x in zip(ORGANS, d):
            O.loc[(O.condition == cond) & (O.genotype == "BKO") & (O.organ == o), "shift_vs_WT_deg"] = x
        ca, cb = np.deg2rad(a.mu_deg.to_numpy()), np.deg2rad(b.mu_deg.to_numpy())
        sa, sb = np.sin(ca - np.angle(np.mean(np.exp(1j * ca)))), np.sin(cb - np.angle(np.mean(np.exp(1j * cb))))
        C.loc[(C.condition == cond) & (C.genotype == "BKO"), "circ_corr_vs_WT"] = float(np.sum(sa * sb) / np.sqrt(np.sum(sa ** 2) * np.sum(sb ** 2)))
    dd, fr = C.set_index(["condition", "genotype"]).loc[("DD", "BKO")], C.set_index(["condition", "genotype"]).loc[("FR", "BKO")]
    tests = pd.DataFrame([dict(test="H319.1 BKO DD no excess", E=dd.E, p=dd.p_excess, passes=bool(dd.p_excess >= 0.05)),
                          dict(test="H319.2 BKO FR rearranged", organs_rayleigh_05=fr.organs_rayleigh_05, order_share=fr.order_share,
                               passes=bool(fr.organs_rayleigh_05 >= 5 and fr.order_share >= 0.5))])
    C.to_csv(OUT / "conditions.csv", index=False)
    O.to_csv(OUT / "organs.csv", index=False)
    tests.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 200)
    print(C.round(3).to_string(index=False))
    print(O.round(3).to_string(index=False))
    print(tests.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
