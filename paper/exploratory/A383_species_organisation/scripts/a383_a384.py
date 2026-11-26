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
  score  (written after the power note) A383: the A263 companion sets of each species (same genes, from
         `A263_companions/output/sets_<species>.csv`) and the core loop fitted as the clock (n_boot=0); H383.1 (rho,
         10,000 tissue permutations, companion q95) and H383.2 (b, 2000 tissue bootstraps, same seeds as the power
         stage); descriptives. A384: 200 expression-matched sets on the genes present in all three tables, matched on
         GTEx fast-death expression (seed 384), fitted in all three species; H384.1 (T, 10,000 permutations, matched
         q95); H384.2 reported only (power 0.75 < 0.8, power note).

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A383_species_organisation/scripts/a383_a384.py --stage power
    RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A383_species_organisation/scripts/a383_a384.py --stage score [--jobs 2]
Outputs (`../output/`): power: baboon_fit.csv, power_a383.csv, power_a384.csv, organs_a384.csv; score:
companion_fits.csv (per set and tissue), sets_a384.csv, tests_a383.csv, companions_a383.csv, own_term_a383.csv,
tests_a384.csv, gains_a384.csv. Tissue- and organ-level aggregates only.
"""
import argparse
import sys
import time
from multiprocessing import Pool
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
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
HUMAN = PAPER / "presentation" / "output" / "organs" / "couplings.csv"
SHEEP = PAPER / "exploratory" / "A378_core_shell_gain" / "output" / "sheep_fit_unread.csv"
SEED = 383
A263_SETS = PAPER / "exploratory" / "A263_companions" / "output"
TRIO = ("PER1", "PER2", "NR1D1")
KEYS = ("donor", "compartment", "depth", "stratum")
FIT_COLS = ["coupling", "own_amplitude", "z2", "psi", "V", "boundary", "log_amp_rel"]
N_PERM, N_BOOT_B = 10000, 2000
W = {}
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


# ---- score stage ---------------------------------------------------------------------------------------------

def waveform(genes, name):
    """A companion waveform rebuilt from its gene list (A.WF amplitudes and phases, as `pseudo_waveforms`)."""
    return Waveform(tuple(genes), A.WF.amp, A.WF.psi, None, None, name)


def a263_sets(species):
    t = pd.read_csv(A263_SETS / f"sets_{species}.csv")
    t = t[t.kind != "clock"]
    return [(k, int(i), g.split(";")) for k, i, g in zip(t.kind, t.set, t.genes)]


def a384_sets(rng, n):
    """Matched sets on the genes present in all three tables, matched on GTEx fast-death mean log fraction."""
    import pyarrow.parquet as pq
    from common.paths import DATA
    other = [set(pq.read_schema(DATA / "derived" / "A263" / f"{sp}_genomewide.parquet").names) - set(KEYS)
             for sp in ("baboon", "sheep")]
    from common import gtex
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = sorted(set(a.var_names) & other[0] & other[1])
    a.file.close()
    pb = A.table("human", genes)
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    assert all(g in lm.index for g in A.WF.genes), "clock genes missing from the common table"
    return [w.genes for w in pseudo_waveforms(lm, n, rng, A.WF, "matched")], len(genes), len(lm)


def task(job):
    species, family, i, genes = job
    pb = W[species]
    wf = A.WF.drop([g for g in TRIO if g in A.WF.genes]) if family == "core" else waveform(genes, f"{family}{i}")
    try:
        tab, _ = fit(species, pb[[c for c in KEYS if c in pb] + list(wf.genes)], wf, 0)
        out = tab[FIT_COLS].rename_axis("tissue").reset_index()
    except Exception as e:  # a set whose fit fails is counted, not dropped silently
        out = pd.DataFrame([dict(tissue=None, error=f"{type(e).__name__}: {e}"[:200])])
    return out.assign(species=species, family=family, set=i)


def rho(tab, per):
    t = tab.reindex(per)[["own_amplitude", "coupling"]].dropna()
    return float(spearmanr(t.own_amplitude, t.coupling).statistic) if len(t) >= 4 else np.nan, t


def b_of(tab, per):
    x, y = part2_rows(tab, per)
    return (slope(x, y) if len(x) >= 4 else np.nan), x, y


def gains(F, which=None):
    """species -> log gain over the A384 organ universe (nan where absent or not fitted)."""
    out = {}
    for j, sp in enumerate(SPECIES):
        tab = F[sp] if which is None else which[sp]
        out[sp] = np.array([float(tab.log_amp_rel.get(o[j], np.nan)) if o[j] is not None else np.nan for o in ORGANS])
    return out


def score(n_sets, jobs, check):
    limit_threads(1)
    t0 = time.time()
    sfx = "_check" if check else ""
    F = fits()
    a384, n_common, n_matchable = a384_sets(np.random.default_rng([SEED, 384]), n_sets)
    pd.DataFrame(dict(set=range(len(a384)), genes=[";".join(g) for g in a384])).to_csv(OUT / f"sets_a384{sfx}.csv",
                                                                                        index=False)
    print(f"fits and A384 sets ({n_common} common genes, {n_matchable} matchable; {time.time() - t0:.0f}s)", flush=True)
    jobs_ = []
    for sp in SPECIES:
        s = a263_sets(sp)
        if check:
            s = [x for x in s if x[0] == "matched"][:n_sets] + [x for x in s if x[0] == "responsive"][:n_sets]
        jobs_ += [(sp, k, i, g) for k, i, g in s] + [(sp, "core", 0, None)]
        jobs_ += [(sp, "a384", i, list(g)) for i, g in enumerate(a384)]
    for sp in SPECIES:
        genes = sorted({g for j in jobs_ if j[0] == sp and j[3] is not None for g in j[3]} | set(A.WF.genes))
        W[sp] = A.table(sp, genes)
    print(f"{len(jobs_)} fits; tables loaded ({time.time() - t0:.0f}s)", flush=True)
    parts = []
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, jobs_, chunksize=1)):
            parts.append(r)
            if check or (k + 1) % 50 == 0 or k + 1 == len(jobs_):
                print(f"{k + 1}/{len(jobs_)} {r.species[0]} {r.family[0]} ({time.time() - t0:.0f}s)", flush=True)
    C = pd.concat(parts, ignore_index=True)
    C.to_csv(OUT / f"companion_fits{sfx}.csv", index=False)
    fitted = {key: g.dropna(subset=["tissue"]).set_index("tissue") for key, g in C.groupby(["species", "family", "set"])}
    n_err = int(C.tissue.isna().sum())
    if check:  # fits only: no statistic before the scoring run
        print(C.groupby(["species", "family"]).agg(fits=("set", "nunique"), tissues=("tissue", "count")).to_string(),
              f"failed fits: {n_err}", sep="\n")
        return

    # A383
    rows, comp, own = [], [], []
    for j, sp in enumerate(SPECIES):
        tab, per = F[sp]
        r, t = rho(tab, per)
        rng = np.random.default_rng([SEED, 1, j])
        null = np.array([spearmanr(t.own_amplitude, rng.permutation(t.coupling.to_numpy())).statistic
                         for _ in range(N_PERM)])
        b, x, y = b_of(tab, per)
        bs = boot_slopes(x, y, np.random.default_rng([SEED, j]), N_BOOT_B)  # the power stage's seed
        lo, hi = np.quantile(bs, [0.05, 0.95])
        core = fitted.get((sp, "core", 0))
        bc, xc, yc = b_of(core, per) if core is not None else (np.nan, [], [])
        bcs = boot_slopes(np.asarray(xc), np.asarray(yc), np.random.default_rng([SEED, 2, j]), N_BOOT_B) \
            if len(xc) >= 4 else np.array([np.nan])
        row = dict(species=sp, tissues_part1=len(t), rho=r, p_perm=float((1 + np.sum(null >= r)) / (1 + N_PERM)),
                   tissues_part2=len(x), b=b, b_lo90=lo, b_hi90=hi, b_excludes_1=bool(hi < 1 or lo > 1),
                   b_below_1=bool(hi < 1), b_excludes_0=bool(lo > 0 or hi < 0),
                   b_core=bc, b_core_lo90=float(np.nanquantile(bcs, 0.05)), b_core_hi90=float(np.nanquantile(bcs, 0.95)),
                   tissues_core=len(xc))
        h1 = bool(r > 0 and row["p_perm"] < 0.05)
        m2psi = []
        for fam in ("matched", "responsive"):
            rs, bs_ = [], []
            for (s_, f_, i_), tb in fitted.items():
                if s_ == sp and f_ == fam:
                    rs.append(rho(tb, per)[0])
                    bs_.append(b_of(tb, per)[0])
                    if fam == "matched":
                        m2psi.append(2 * tb.psi.reindex(per))
            rs, bs_ = np.array(rs, float), np.array(bs_, float)
            row[f"{fam}_valid"] = int(np.isfinite(rs).sum())
            row[f"{fam}_rho_q95"] = float(np.nanquantile(rs, 0.95))
            row[f"{fam}_pct_below_clock"] = float(100 * np.mean(rs[np.isfinite(rs)] < r))
            h1 &= bool(r > row[f"{fam}_rho_q95"])
            comp.append(dict(species=sp, family=fam, rho_median=float(np.nanmedian(rs)),
                             rho_q05=float(np.nanquantile(rs, 0.05)), rho_q95=float(np.nanquantile(rs, 0.95)),
                             b_valid=int(np.isfinite(bs_).sum()), b_median=float(np.nanmedian(bs_)),
                             b_q05=float(np.nanquantile(bs_, 0.05)), b_q95=float(np.nanquantile(bs_, 0.95))))
        row["H383_1"], row["H383_2"] = h1, row["b_excludes_1"]
        rows.append(row)
        med = pd.concat(m2psi, axis=1).median(axis=1)
        own.append(pd.DataFrame(dict(species=sp, tissue=per, own_clock=(2 * tab.psi.reindex(per)).to_numpy(float),
                                     own_matched_median=med.to_numpy(float),
                                     log_ratio=np.log((2 * tab.psi.reindex(per)) / med).to_numpy(float))))
    T383, comp, own = pd.DataFrame(rows), pd.DataFrame(comp), pd.concat(own, ignore_index=True)

    # A384
    G = gains({sp: F[sp][0] for sp in SPECIES})
    r_obs = pair_rhos(G)
    T = float(np.mean(list(r_obs.values())))
    nT, nbs = perm_null(G, np.random.default_rng([SEED, 384, 1]), N_PERM)
    Tm, bsm = [], []
    for i in range(len(a384)):
        tabs = {sp: fitted.get((sp, "a384", i)) for sp in SPECIES}
        if any(v is None for v in tabs.values()):
            Tm.append(np.nan), bsm.append(np.nan)
            continue
        rr = pair_rhos(gains(None, tabs))
        Tm.append(np.mean(list(rr.values()))), bsm.append(rr[("baboon", "sheep")])
    Tm, bsm = np.array(Tm, float), np.array(bsm, float)
    t384 = dict(T=T, p_perm=float((1 + np.sum(nT >= T)) / (1 + N_PERM)), matched_valid=int(np.isfinite(Tm).sum()),
                matched_T_median=float(np.nanmedian(Tm)), matched_T_q95=float(np.nanquantile(Tm, 0.95)),
                matched_pct_below=float(100 * np.mean(Tm[np.isfinite(Tm)] < T)))
    t384["H384_1"] = bool(T > 0 and t384["p_perm"] < 0.05 and T > t384["matched_T_q95"])
    for a, b in PAIRS:
        m = ~np.isnan(G[a]) & ~np.isnan(G[b])
        t384[f"rho_{a}_{b}"], t384[f"organs_{a}_{b}"] = float(r_obs[(a, b)]), int(m.sum())
    t384.update(H384_2_reported_only=True, p_perm_baboon_sheep=float((1 + np.sum(nbs >= r_obs[("baboon", "sheep")]))
                                                                     / (1 + N_PERM)),
                matched_baboon_sheep_q95=float(np.nanquantile(bsm, 0.95)), a384_common_genes=n_common,
                a384_matchable_genes=n_matchable, failed_fits=n_err)
    for (a, b) in PAIRS[:2]:  # per-pair permutation p from independent permutations (reported)
        rng = np.random.default_rng([SEED, 384, 2, SPECIES.index(b)])
        m = ~np.isnan(G[a]) & ~np.isnan(G[b])
        nl = np.array([spearmanr(G[a][m], rng.permutation(G[b][m])).statistic for _ in range(N_PERM)])
        t384[f"p_perm_{a}_{b}"] = float((1 + np.sum(nl >= r_obs[(a, b)])) / (1 + N_PERM))
    gt = pd.DataFrame({**{sp: [o[j] for o in ORGANS] for j, sp in enumerate(SPECIES)},
                       **{f"log_gain_{sp}": G[sp] for sp in SPECIES}})
    for sp in SPECIES:
        gt[f"rank_{sp}"] = gt[f"log_gain_{sp}"].rank(ascending=False)

    T383.to_csv(OUT / f"tests_a383{sfx}.csv", index=False)
    comp.to_csv(OUT / f"companions_a383{sfx}.csv", index=False)
    own.to_csv(OUT / f"own_term_a383{sfx}.csv", index=False)
    pd.Series(t384).to_csv(OUT / f"tests_a384{sfx}.csv", header=["value"])
    gt.to_csv(OUT / f"gains_a384{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(T383.round(3).T.to_string(), comp.round(3).to_string(index=False), sep="\n")
    print(own.groupby("species").log_ratio.describe().round(3).to_string())
    print(pd.Series(t384).to_string(), gt.round(2).to_string(index=False), f"({time.time() - t0:.0f}s)", sep="\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--check", action="store_true", help="3 sets per family, fits only (*_check outputs, no statistic)")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    power() if a.stage == "power" else score(3 if a.check else 200, a.jobs, a.check)


if __name__ == "__main__":
    main()
