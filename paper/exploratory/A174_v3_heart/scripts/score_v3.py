"""A174 (V3 port) scoring: gate, PM - AM host displacement and its angle on the current estimator (registered 2026-11-01).

Primary: the 16 patients, compartments cardiomyocyte, fibroblast, endothelial, neuron (floor 20 cells), learned A17 heart
anchors present here (fibroblast = reference, cardiomyocyte, neuron), default waveform. Descriptive sensitivities: each
surgical arm, the archive floor of 50 cells, the stromal-vascular anchors (fibroblast, endothelial), myeloid added, the
archive left-ventricle waveform.
    python score_v3.py [--n-gate 2000] [--n-perm 5000]
Output: output/v3_scores.csv (one row per analysis), output/v3_host.csv (primary host state per patient).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync.config import limit_threads

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from common.paths import exploratory_dir  # noqa: E402

HOURS = 24 / (2 * np.pi)
FORWARD_H = 18.0  # mean window midpoint (8 h, 16 h) + 6 h: where forward time moves the host
PRIMARY = ["heart:cardiomyocyte", "heart:fibroblast", "heart:endothelial", "heart:neuron"]
LEARNED = ["heart:fibroblast", "heart:cardiomyocyte", "heart:neuron"]  # first = reference (ties go to the first)
FIXED = ["heart:fibroblast", "heart:endothelial"]


def displacement(host: pd.DataFrame, window: pd.Series, n_perm: int, seed: int = 0) -> dict:
    """Mean PM minus mean AM host state: length, window-label permutation p, angle (h) and error to forward time."""
    f, w = host[["x", "y"]].to_numpy(), window.reindex(host.index).to_numpy()
    d = lambda lab: f[lab == "PM"].mean(0) - f[lab == "AM"].mean(0)
    obs = d(w)
    rng = np.random.default_rng(seed)
    null = np.array([np.linalg.norm(d(rng.permutation(w))) for _ in range(n_perm)])
    ang = np.mod(np.arctan2(obs[1], obs[0]) * HOURS, 24)
    return {"length": np.linalg.norm(obs), "p_separation": (1 + np.sum(null >= np.linalg.norm(obs))) / (1 + n_perm),
            "angle_h": ang, "error_h": (ang - FORWARD_H + 12) % 24 - 12, "n_am": int(np.sum(w == "AM")),
            "n_pm": int(np.sum(w == "PM"))}


def score(name, pb, donors, comps, anchors, args, waveform="human"):
    ds = cs.Dataset(pb[pb.compartment.isin(comps)], anchors=anchors, donors=donors, waveform=waveform, name=name)
    g = ds.gate(n_perm=args.n_gate, seed=0)
    host = ds.host()
    row = {"analysis": name, "n_patients": host.shape[0], "compartments": ";".join(comps), "anchors": ";".join(anchors),
           "gate_passes": g.passes, "p_pairing": g.p_pairing, "p_acrophase": g.p_acrophase,
           **displacement(host, donors.window, args.n_perm)}
    print(pd.Series(row).to_string(), flush=True)
    return row, host


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-gate", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=5000)
    args = ap.parse_args()
    limit_threads(4)
    out = exploratory_dir("A174_v3_heart")
    pb, donors = pd.read_parquet(out / "v3_heart.parquet"), pd.read_parquet(out / "v3_heart_donors.parquet")
    pb20, pb50 = pb[pb.n_cells >= 20], pb[pb.n_cells >= 50]
    rows = []
    r, host = score("primary", pb20, donors, PRIMARY, LEARNED, args)
    rows.append(r)
    host.join(donors[["group", "window"]]).to_csv(out / "v3_host.csv")
    for arm in ("NCR", "ACR"):
        keep = donors.index[donors.group == arm]
        rows.append(score(f"arm_{arm}", pb20[pb20.donor.isin(keep)], donors.loc[keep], PRIMARY, LEARNED, args)[0])
    rows.append(score("floor_50", pb50, donors, PRIMARY, LEARNED, args)[0])
    rows.append(score("anchors_fixed", pb20, donors, PRIMARY, FIXED, args)[0])
    rows.append(score("with_myeloid", pb20, donors, PRIMARY + ["heart:myeloid"], LEARNED, args)[0])
    rows.append(score("waveform_lv", pb20, donors, PRIMARY, LEARNED, args, waveform="Heart - Left Ventricle")[0])
    pd.DataFrame(rows).to_csv(out / "v3_scores.csv", index=False)


if __name__ == "__main__":
    main()
