"""A378 (registered 2026-11-25): core against shell organ gain (H378.1) and conservation of the GTEx same-class gain
ratios (H378.2) in GTEx, baboon and sheep.

Stages:
  power  standard error of every contrast and its MDE (2.49 SE); no gain ratio is computed or printed.
         GTEx: per-tissue SE of the log relative amplitude from the reference-map bootstrap interval (95 %, 200 draws,
         `presentation/output/organs/couplings.csv`), SE = width / 3.92. Baboon: residual-bootstrap SD of log |g| per
         organ (500 draws, the A330 estimator). Sheep: one clocksync fit of the 45 compartments (A263 / A289 setup:
         peripheral anchors read leave-one-out, mixed strata, 200 animal bootstraps), SE from its 95 % intervals; the
         fit table is written to `sheep_fit_unread.csv` and not read at this stage.
         Choices fixed here, before scoring: (i) clocksync exposes bootstrap replicates of couplings only, so contrast
         SEs combine per-tissue SEs as independent (paired resampling as registered is not available without refitting;
         positively correlated tissue errors make this conservative); the scoring interval uses the same rule;
         (ii) a side with several tissues is the mean of their log gains; (iii) a species' part-1 value is the mean of
         its pair contrasts (GTEx: arteries and skin; baboon: muscle and fat; sheep: fat), and its MDE is computed on it.
  score  (written after the power note)

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A378_core_shell_gain/scripts/a378.py --stage power
Outputs (`../output/`): power_tissues.csv, power_contrasts.csv, sheep_fit_unread.csv. Aggregates only.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_2" / "scripts" / "organs"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_brain_periphery_species as SP  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import data  # noqa: E402

OUT = HERE.parents[1] / "output"
GTEX = PAPER / "presentation" / "output" / "organs" / "couplings.csv"
W = 2 * np.pi / 24
SEED = 378
# contrasts: (part, species, name, core tissues, shell tissues)
CONTRASTS = [
    ("1", "human", "arteries", ["Artery - Aorta", "Artery - Coronary"], ["Artery - Tibial"]),
    ("1", "human", "skin", ["Skin - Not Sun Exposed (Suprapubic)"], ["Skin - Sun Exposed (Lower leg)"]),
    ("1", "baboon", "muscle", ["MUA"], ["MUG"]),
    ("1", "baboon", "fat", ["OMF", "WAM", "WAP", "WAR"], ["WAS"]),
    ("1", "sheep", "fat", ["Perirenal adipose", "Caul adipose"],
     ["Subcutaneous adipose", "Tail adipose", "Intermuscular adipose"]),
    ("1", "sheep", "fat, tail apart", ["Perirenal adipose", "Caul adipose"], ["Tail adipose"]),
    ("2", "human", "squamous", ["Esophagus - Mucosa"], ["Skin - Sun Exposed (Lower leg)"]),
    ("2", "baboon", "squamous", ["OES"], ["SKI"]),
    ("2", "sheep", "squamous", ["Rumen", "Reticulum", "Omasum"], ["Skin"]),
    ("2", "human", "fat", ["Adipose - Visceral (Omentum)"], ["Adipose - Subcutaneous"]),
    ("2", "baboon", "fat", ["OMF"], ["WAS"]),
    ("2", "sheep", "fat", ["Caul adipose"], ["Subcutaneous adipose"]),
]
SPECIES_PART1 = {"human": ["arteries", "skin"], "baboon": ["muscle", "fat"], "sheep": ["fat"]}


def baboon_se(codes, wf, rng, n=500):
    """Residual-bootstrap SD of log |g| per organ (A330 estimator); the point gain is not returned."""
    pb, donors = data.load("baboon")
    zt = donors.zt.astype(float)
    out = {}
    for k in codes:
        sub = pb[pb.compartment == k]
        z = zt.reindex(sub.donor).to_numpy(float)
        g = list(wf.genes)
        cpm = sub[g].to_numpy(float) / sub.depth.to_numpy(float)[:, None] * 1e6
        ex = np.median(cpm, 0) >= 5
        y = np.log(cpm[:, ex] + 0.5)
        x = (wf.amp * np.exp(1j * wf.psi))[ex]
        X = np.column_stack([np.ones_like(z), np.cos(W * z), np.sin(W * z)])
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        fit, res = X @ b, y - X @ b
        boot = []
        for _ in range(n):
            bb = np.linalg.lstsq(X, fit + res[rng.integers(0, len(z), len(z))], rcond=None)[0]
            c = np.abs(bb[1] - 1j * bb[2]) * np.exp(1j * np.arctan2(bb[2], bb[1]))
            boot.append(abs(complex(np.conj(x) @ c / (np.conj(x) @ x))))
        out[k] = dict(se=float(np.std(np.log(boot), ddof=1)), n=len(z))
    return out


def sheep_se(wf):
    pb, _ = data.load("sheep")
    brain, per, excl = SP.groups("sheep", sorted(pb.compartment.unique()))
    cols = ["donor", "compartment", "depth", "stratum"] + list(wf.genes)
    ds = cs.Dataset(pb[cols], anchors=per, waveform=wf, mixed_strata=True, name="sheep_a378", min_donors=5)
    tab = ds.fit(n_boot=200, gate=False, targets=brain + per + excl, seed=SEED).table
    tab.rename_axis("tissue").to_csv(OUT / "sheep_fit_unread.csv")  # point estimates: not read before scoring
    se = (tab.log_amp_hi - tab.log_amp_lo) / 3.92
    return {t: dict(se=float(se[t]), n=int(tab.n[t])) for t in tab.index}


def power():
    limit_threads(4)
    t0 = time.time()
    wf = cs.gtex()
    g = pd.read_csv(GTEX, index_col="tissue")
    tissues = {"human": {t: dict(se=float((g.log_amp_hi[t] - g.log_amp_lo[t]) / 3.92), n=int(g.n[t])) for t in g.index}}
    codes = sorted({t for _, sp, _, a, b in CONTRASTS if sp == "baboon" for t in a + b})
    tissues["baboon"] = baboon_se(codes, wf, np.random.default_rng(SEED))
    print(f"baboon SEs ({time.time() - t0:.0f}s)", flush=True)
    tissues["sheep"] = sheep_se(wf)
    print(f"sheep fit ({time.time() - t0:.0f}s)", flush=True)
    rows = []
    for part, sp, name, a, b in CONTRASTS:
        se = np.sqrt(sum(tissues[sp][t]["se"] ** 2 for t in a) / len(a) ** 2
                     + sum(tissues[sp][t]["se"] ** 2 for t in b) / len(b) ** 2)
        rows.append(dict(part=part, species=sp, contrast=name, core=";".join(a), shell=";".join(b), se=se,
                         mde=2.49 * se, n_min=min(tissues[sp][t]["n"] for t in a + b)))
    con = pd.DataFrame(rows)
    for sp, names in SPECIES_PART1.items():
        s = con[(con.part == "1") & (con.species == sp) & con.contrast.isin(names)].se
        se = float(np.sqrt((s ** 2).sum()) / len(s))
        con.loc[len(con)] = dict(part="1", species=sp, contrast="species mean", core="", shell="", se=se,
                                 mde=2.49 * se, n_min=np.nan)
    con["decidable_part1"] = np.where(con.contrast == "species mean", con.mde <= np.log(1.5), np.nan)
    tis = pd.DataFrame([dict(species=sp, tissue=t, **v) for sp, d in tissues.items() for t, v in d.items()
                        if any(t in a + b for _, s2, _, a, b in CONTRASTS if s2 == sp)])
    tis.to_csv(OUT / "power_tissues.csv", index=False)
    con.to_csv(OUT / "power_contrasts.csv", index=False)
    pd.set_option("display.width", 200)
    print(tis.round(3).to_string(index=False))
    print(con.round(3).to_string(index=False), f"\n({time.time() - t0:.0f}s)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power"], required=True)
    ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    power()


if __name__ == "__main__":
    main()
