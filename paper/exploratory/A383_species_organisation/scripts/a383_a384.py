"""A383 (amplitude-coupling rule and own-term scaling, three species) and A384 (organ gain ranks across species),
registered 2026-11-25.

One estimator in all three species (the A263 R2 setup): `Dataset.fit(gate=False)`, peripheral tissues as anchors read
leave-one-out. Human: the reference map (`presentation/output/organs/couplings.csv`, GTEx v10 fast deaths, 200 donor
bootstraps); sheep: the A378 fit (same setup, 200 animal bootstraps, `A378_core_shell_gain/output/sheep_fit_unread.csv`);
baboon: fitted here (200 animal bootstraps). Coupling c = S / (S + 2 psi), S = z2 V.

Stages:
  power  A383: peripheral tissues per species, analytic power of part 1, tissue-bootstrap SE of the slope b of log 2 psi
         on log S (b not printed), one companion fit timed per species. A384: organ map, bootstrap SE of each log gain,
         simulated power of H384.1 and H384.2 at a latent correlation of 0.6. No cross-species statistic, no companion
         statistic.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A383_species_organisation/scripts/a383_a384.py --stage power
Outputs (`../output/`): baboon_fit.csv, power_a383.csv, power_a384.csv, organs_a384.csv. Aggregates only.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER / "exploratory" / "A263_companions" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm, spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
import score_a263 as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

OUT = HERE.parents[1] / "output"
HUMAN = PAPER / "presentation" / "output" / "organs" / "couplings.csv"
SHEEP = PAPER / "exploratory" / "A378_core_shell_gain" / "output" / "sheep_fit_unread.csv"
SEED = 383
SPECIES = ("human", "baboon", "sheep")
GTEX_COMPANION_MEDIAN = {"matched": -0.160, "responsive": 0.109}  # A271 H271.2 (periphery), planning values
# A384 organ map: (GTEx, baboon, sheep); None = absent in that species
ORGANS = [
    ("Stomach", "STF", "Abomasum"), ("Adipose - Visceral (Omentum)", "OMF", "Caul adipose"),
    ("Adipose - Subcutaneous", "WAS", "Subcutaneous adipose"), ("Adrenal Gland", "ADC", "Adrenal"),
    ("Artery - Aorta", "AOR", "Artery"), ("Colon - Transverse", "ASC", "Colon"),
    ("Heart - Left Ventricle", "HEA", "Heart"), ("Kidney - Cortex", "KIC", "Kidney"), ("Liver", "LIV", "Liver"),
    ("Lung", "LUN", "Lung"), ("Muscle - Skeletal", "MUG", "Muscle"),
    ("Skin - Not Sun Exposed (Suprapubic)", "SKI", "Skin"), ("Thyroid", "THR", "Thyroid"),
    ("Esophagus - Muscularis", "OES", None), ("Colon - Sigmoid", "DEC", None), ("Prostate", "PRO", None),
    ("Testis", "TES", None),
    (None, "WAR", "Perirenal adipose"), (None, "DUO", "Duodenum"), (None, "ILE", "Ileum"), (None, "CEC", "Cecum"),
    (None, "SPL", "Spleen"), (None, "BOM", "Bone marrow"), (None, "MEL", "Lymph node"),
]
PAIRS = (("human", "baboon"), ("human", "sheep"), ("baboon", "sheep"))


def setup(species, pb):
    brain, per, excl = A.groups(species, sorted(pb.compartment.unique()))
    cols = [c for c in ("donor", "compartment", "depth", "stratum") if c in pb]
    return brain, per, excl, cols


def fit(species, pb, wf, n_boot, seed=SEED):
    brain, per, excl, cols = setup(species, pb)
    ds = cs.Dataset(pb[cols + list(wf.genes)], anchors=per, waveform=wf, mixed_strata=species == "sheep",
                    name=f"{species}_a383", min_donors=5)
    return ds.fit(n_boot=n_boot, gate=False, targets=brain + per + excl, seed=seed).table, per


def fits():
    """Clock fits of the three species (table, peripheral tissues); the baboon fit is computed and written once."""
    h = pd.read_csv(HUMAN, index_col="tissue")
    out = {"human": (h, list(h.index[h.group == "periphery"]))}
    pb = A.table("baboon", A.WF.genes)
    f = OUT / "baboon_fit.csv"
    if not f.exists():
        tab, _ = fit("baboon", pb, A.WF, 200)
        tab.rename_axis("tissue").to_csv(f)
    out["baboon"] = (pd.read_csv(f, index_col="tissue"), setup("baboon", pb)[1])
    s = pd.read_csv(SHEEP, index_col="tissue")
    out["sheep"] = (s, setup("sheep", A.table("sheep", A.WF.genes))[1])
    return out


def part2_rows(tab, per):
    t = tab.reindex(per)
    S = t.z2 * t.V
    keep = (S > 0) & (t.psi > 0) & ~t.boundary.astype(bool)
    return np.log(S[keep].to_numpy(float)), np.log(2 * t.psi[keep].to_numpy(float))


def slope(x, y):
    return float(np.polyfit(x, y, 1)[0])


def boot_slopes(x, y, rng, n=2000):
    out = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        if np.ptp(x[i]) > 0:
            out.append(slope(x[i], y[i]))
    return np.array(out)


def spearman_power(n, rho_true, threshold):
    se = np.sqrt(1.06 / (n - 3))
    return float(1 - norm.cdf((np.arctanh(threshold) - np.arctanh(rho_true)) / se))


def time_companion(species, rng):
    pb = A.table(species)
    genes = [c for c in pb.columns if c not in ("donor", "compartment", "depth", "stratum")]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    w = pseudo_waveforms(lm, 1, rng, A.WF, "matched")[0]
    t0 = time.time()
    fit(species, pb[[c for c in pb.columns if c in ("donor", "compartment", "depth", "stratum")] + list(w.genes)], w, 0)
    return time.time() - t0


def power_a383(F):
    rows = []
    for sp, (tab, per) in F.items():
        t = tab.reindex(per)
        n1 = int(t[["own_amplitude", "coupling"]].notna().all(axis=1).sum())
        x, y = part2_rows(tab, per)
        se_b = float(np.std(boot_slopes(x, y, np.random.default_rng([SEED, SPECIES.index(sp)])), ddof=1))
        r = dict(species=sp, tissues_part1=n1, tissues_part2=len(x), se_b=se_b, mde_1_minus_b=2.49 * se_b,
                 part2_decidable=bool(2.49 * se_b <= 1))
        for k, m in GTEX_COMPANION_MEDIAN.items():
            thr = m + 1.645 / np.sqrt(n1 - 1)
            r[f"threshold_{k}"] = thr
            r[f"power_{k}_rho06"] = spearman_power(n1, 0.6, thr)
        r["part1_decidable"] = bool(r["power_responsive_rho06"] >= 0.8)
        if sp != "human":
            r["companion_fit_s"] = time_companion(sp, np.random.default_rng([SEED, 9]))
        rows.append(r)
        print(sp, "done", flush=True)
    return pd.DataFrame(rows)


def organ_table(F):
    rows = []
    for o in ORGANS:
        r = {}
        for sp, k in zip(SPECIES, o):
            tab = F[sp][0]
            r[sp] = k
            r[f"se_{sp}"] = (float((tab.log_amp_hi[k] - tab.log_amp_lo[k]) / 3.92)
                             if k is not None and k in tab.index else np.nan)
            if k is not None and k not in tab.index:
                r[sp] = None
        rows.append(r)
    return pd.DataFrame(rows)


def pair_rhos(vals):
    """vals: species -> array over the organ universe (nan if absent); Spearman per pair over shared organs."""
    out = {}
    for a, b in PAIRS:
        m = ~np.isnan(vals[a]) & ~np.isnan(vals[b])
        out[(a, b)] = spearmanr(vals[a][m], vals[b][m]).statistic if m.sum() >= 4 else np.nan
    return out


def perm_null(vals, rng, n):
    T, Tbs = [], []
    for _ in range(n):
        v = {}
        for sp, x in vals.items():
            y = x.copy()
            m = ~np.isnan(y)
            y[m] = rng.permutation(y[m])
            v[sp] = y
        r = pair_rhos(v)
        T.append(np.mean(list(r.values())))
        Tbs.append(r[("baboon", "sheep")])
    return np.array(T), np.array(Tbs)


def power_a384(o, sd_true, rng, n_sim=2000, n_perm=500, rho=0.6):
    present = {sp: o[sp].notna().to_numpy() for sp in SPECIES}
    se = {sp: o[f"se_{sp}"].to_numpy(float) for sp in SPECIES}
    hit_T = hit_bs = 0
    for _ in range(n_sim):
        f = rng.normal(0, 1, len(o))
        vals = {}
        for sp in SPECIES:
            v = sd_true * (np.sqrt(rho) * f + np.sqrt(1 - rho) * rng.normal(0, 1, len(o)))
            v = v + rng.normal(0, np.nan_to_num(se[sp]))
            vals[sp] = np.where(present[sp], v, np.nan)
        r = pair_rhos(vals)
        T, Tbs = np.mean(list(r.values())), r[("baboon", "sheep")]
        nT, nbs = perm_null(vals, rng, n_perm)
        hit_T += (1 + np.sum(nT >= T)) / (n_perm + 1) < 0.05
        hit_bs += (1 + np.sum(nbs >= Tbs)) / (n_perm + 1) < 0.05
    return hit_T / n_sim, hit_bs / n_sim


def power():
    limit_threads(4)
    t0 = time.time()
    F = fits()
    print(f"fits ({time.time() - t0:.0f}s)", flush=True)
    p383 = power_a383(F)
    p383.to_csv(OUT / "power_a383.csv", index=False)
    o = organ_table(F)
    o.drop(columns=[]).to_csv(OUT / "organs_a384.csv", index=False)
    h = F["human"][0]
    hk = o.human.dropna()
    sd_true = float(np.std(h.log_amp_rel.reindex(hk).to_numpy(float), ddof=1))  # GTEx: not blind
    pT, pbs = power_a384(o, sd_true, np.random.default_rng([SEED, 384]))
    n = {f"organs_{a}_{b}": int((o[a].notna() & o[b].notna()).sum()) for a, b in PAIRS}
    p384 = dict(**n, sd_true_gtex=sd_true, median_se_human=float(o.se_human.median()),
                median_se_baboon=float(o.se_baboon.median()), median_se_sheep=float(o.se_sheep.median()),
                power_H384_1_rho06=pT, power_H384_2_rho06=pbs, gate_H384_1=bool(pT >= 0.8), gate_H384_2=bool(pbs >= 0.8))
    pd.Series(p384).to_csv(OUT / "power_a384.csv", header=["value"])
    pd.set_option("display.width", 220)
    print(p383.round(3).to_string(index=False))
    print(o.drop(columns=["human", "baboon", "sheep"]).round(3).assign(organ=o.human.fillna(o.sheep)).to_string(index=False))
    print(pd.Series(p384).to_string(), f"\n({time.time() - t0:.0f}s)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power"], required=True)
    ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    power()


if __name__ == "__main__":
    main()
