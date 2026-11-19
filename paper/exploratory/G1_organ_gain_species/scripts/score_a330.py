"""A330 (registered, check G1): do GTEx organ gains (sqrt S) rank like the timed clock amplitudes of the same organs in
the baboon atlas (Mure 2018)?

    python score_a330.py --stage power   # precondition (|g| vs ZT permutations) and power; no cross-species statistic
    python score_a330.py --stage score   # rho, permutation p, matched companion
Outputs (`../output/`): organs.csv, power.csv | tests.csv, sets.csv. Aggregates only.
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_1" / "scripts" / "species"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import _shared as SH  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import data  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/derived/A263/baboon_genomewide.parquet")
GAINS = PAPER / "exploratory" / "A294_floor" / "output" / "tissues.csv"
MAP = {"ADC": "Adrenal Gland", "AOR": "Artery - Aorta", "ASC": "Colon - Transverse", "DEC": "Colon - Sigmoid",
       "HEA": "Heart - Left Ventricle", "KIC": "Kidney - Cortex", "LIV": "Liver", "LUN": "Lung",
       "MUG": "Muscle - Skeletal", "OES": "Esophagus - Muscularis", "OMF": "Adipose - Visceral (Omentum)",
       "WAS": "Adipose - Subcutaneous", "PIT": "Pituitary", "PRO": "Prostate",
       "SKI": "Skin - Not Sun Exposed (Suprapubic)", "STF": "Stomach", "TES": "Testis", "THR": "Thyroid"}
W = 2 * np.pi / 24


def load():
    pb = pd.read_parquet(DATA)
    _, donors = data.load("baboon")
    zt = donors.zt.astype(float)
    return pb[pb.compartment.isin(MAP)].reset_index(drop=True), zt


def gain(sub, zt, wf):
    return SH.timed_gain(sub, zt.reindex(sub.donor).to_numpy(float), wf)[0]


def perm_q95(sub, zt, wf, rng, n=2000):
    z = zt.reindex(sub.donor).to_numpy(float)
    g = list(wf.genes)
    cpm = sub[g].to_numpy(float) / sub.depth.to_numpy(float)[:, None] * 1e6
    ex = np.median(cpm, 0) >= 5
    y = np.log(cpm[:, ex] + 0.5)
    x = (wf.amp * np.exp(1j * wf.psi))[ex]
    X = np.column_stack([np.ones_like(z), np.cos(W * z), np.sin(W * z)])

    def g_of(t):
        Xt = np.column_stack([np.ones_like(t), np.cos(W * t), np.sin(W * t)])
        b = np.linalg.lstsq(Xt, y, rcond=None)[0]
        c = b[1] - 1j * b[2]  # amplitude * exp(i W acrophase): a cos + b sin = A cos(W t - W phi)
        c = np.abs(c) * np.exp(1j * np.arctan2(b[2], b[1]))
        return complex(np.conj(x) @ c / (np.conj(x) @ x))

    obs = abs(g_of(z))
    null = np.array([abs(g_of(rng.permutation(z))) for _ in range(n)])
    # residual bootstrap SE of log |g|
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    fit, res = X @ b, y - X @ b
    boot = []
    for _ in range(500):
        yb = fit + res[rng.integers(0, len(z), len(z))]
        bb = np.linalg.lstsq(X, yb, rcond=None)[0]
        c = np.abs(bb[1] - 1j * bb[2]) * np.exp(1j * np.arctan2(bb[2], bb[1]))
        boot.append(abs(complex(np.conj(x) @ c / (np.conj(x) @ x))))
    return obs, float(np.quantile(null, 0.95)), float(np.std(np.log(boot), ddof=1))


def power(logS, se, rho_true, rng, n_sim=2000, n_perm=500):
    n = len(logS)
    hits = 0
    for _ in range(n_sim):
        if rho_true >= 1:
            s = 0.0
        else:  # species deviation with SD s so that corr(logS, logS + dev) ~ rho_true (Pearson on the latent scale)
            s = np.std(logS) * np.sqrt(1 / rho_true ** 2 - 1)
        y = logS + rng.normal(0, s, n) + rng.normal(0, se)
        r0 = spearmanr(y, logS).statistic
        null = np.array([spearmanr(rng.permutation(y), logS).statistic for _ in range(n_perm)])
        hits += (1 + np.sum(null >= r0)) / (n_perm + 1) < 0.05
    return hits / n_sim


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    pb, zt = load()
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    wf = cs.gtex()
    S = pd.read_csv(GAINS, index_col=0).S
    if a.stage == "power":
        rng = np.random.default_rng(330)
        rows = []
        for b, gt in MAP.items():
            sub = pb[pb.compartment == b]
            obs, q95, se = perm_q95(sub, zt, wf, rng)
            rows.append(dict(baboon=b, gtex=gt, n_animals=len(sub), g_abs=obs, perm_q95=q95, passes=obs > q95,
                             se_log_g=se, sqrt_S=float(np.sqrt(S[gt]))))
        o = pd.DataFrame(rows)
        # cross-check with A36's g_abs (same estimator)
        a36 = pd.read_csv(PAPER / "results/figure_1/output/species/a36_tissues.csv", index_col=0).g_abs
        o["a36_g_abs"] = a36.reindex(o.baboon).to_numpy()
        o.to_csv(OUT / "organs.csv", index=False)
        p = o[o.passes]
        logS = np.log(np.sqrt(p.sqrt_S.to_numpy() ** 2))
        pw = dict(organs_matched=len(o), organs_passing=len(p),
                  power_rho1=power(logS, p.se_log_g.to_numpy(), 1.0, rng),
                  power_rho06=power(logS, p.se_log_g.to_numpy(), 0.6, rng))
        k = rng.choice(len(p), 8, replace=False) if len(p) >= 8 else np.arange(len(p))
        pw["power_rho06_n8"] = power(logS[k], p.se_log_g.to_numpy()[k], 0.6, rng)
        pw["precondition"] = bool(len(p) >= 8)
        pw["gate"] = bool(pw["precondition"] and pw["power_rho06"] >= 0.8)
        pd.Series(pw).to_csv(OUT / "power.csv", header=["value"])
        print(o.drop(columns=["sqrt_S"]).round(3).to_string())
        print(pd.Series(pw).to_string())
        return
    o = pd.read_csv(OUT / "organs.csv")
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("power gate not passed: not decidable, not scored")
    p = o[o.passes].reset_index(drop=True)
    rng = np.random.default_rng(330)
    rho = spearmanr(p.g_abs, p.sqrt_S).statistic
    null = np.array([spearmanr(rng.permutation(p.g_abs.to_numpy()), p.sqrt_S).statistic for _ in range(10_000)])
    pval = (1 + np.sum(null >= rho)) / 10_001
    genes = [c for c in pb.columns if c not in ("donor", "compartment", "depth")]
    x = pb[genes].to_numpy(float) / pb.depth.to_numpy(float)[:, None] * 1e6
    lm = pd.Series(np.log(x.mean(0) + 1e-12), index=genes)[np.median(x, 0) >= 1]
    sets = pseudo_waveforms(lm, 200, np.random.default_rng([330, 1]), wf, kind="matched")
    rows = []
    for j, w in enumerate(sets):
        g = {b: abs(gain(pb[pb.compartment == b], zt, w)) for b in p.baboon}
        gs = np.array([g[b] for b in p.baboon])
        rows.append(dict(set=j, rho=spearmanr(gs, p.sqrt_S).statistic, **{f"g_{b}": g[b] for b in p.baboon}))
    st = pd.DataFrame(rows)
    st.to_csv(OUT / "sets.csv", index=False)
    med = st[[f"g_{b}" for b in p.baboon]].median().to_numpy()
    rel = p.g_abs.to_numpy() / med
    q95 = float(st.rho.quantile(0.95))
    t = dict(test="H330.1", organs=len(p), rho=float(rho), p_one_sided=float(pval), matched_median=float(st.rho.median()),
             matched_q95=q95, pct_matched=float((st.rho < rho).mean() * 100),
             rho_relative_to_matched=float(spearmanr(rel, p.sqrt_S).statistic))
    t["verdict"] = ("PASS" if rho > 0 and pval < 0.05 and rho > q95 else
                    "not clock-specific" if rho > 0 and pval < 0.05 else "FAIL")
    pd.DataFrame([t]).to_csv(OUT / "tests.csv", index=False)
    print(pd.Series(t).to_string())


if __name__ == "__main__":
    main()
