"""A332 D2 (descriptive, registered): brain death per organ, the arrows of Figure 5. No test.

Per A317 body tissue: phase change = circular mean over the tissue's pairs (oriented tissue minus partner) of the pair
phase in organ donors minus that in ventilated post-mortem donors (A317 states; organ donors at the post-mortem n, 20
draws, circular mean over draws), with a 1,000-draw donor bootstrap (both groups resampled; one organ-donor draw each);
loudness = A315 `tissues.csv` (clock log variance ratio L minus its matched median, amplitude exp(diff / 2)).
Aggregates only (GTEx controlled access).

    python d2_per_organ.py
Output: `../output/d2_per_organ.csv`.
"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(EXP / "A109_coupling_nulls" / "scripts"))
sys.path.insert(0, str(EXP / "A132_brain_death" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
import compute_brain_death as B  # noqa: E402
import compute_coupling_nulls as N  # noqa: E402
import score_a317 as S317  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
HOURS = 24 / (2 * np.pi)


def per_tissue(Z, M, od, pm, mask):
    """Per tissue: circular mean over its masked pairs of (phase_od - phase_pm), oriented tissue minus partner."""
    d = np.angle(np.exp(1j * (H.pair_phase(N._H_of(Z, M, od)) - H.pair_phase(N._H_of(Z, M, pm)))))
    out = np.full(len(H.BODY), np.nan + 0j)
    for j in range(len(H.BODY)):
        v = [(d[k] if a == j else -d[k]) for k, (a, b) in enumerate(H.PAIRS) if mask[k] and j in (a, b)]
        if v:
            out[j] = np.mean(np.exp(1j * np.array(v)))
    return out


def main():
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    rng = np.random.default_rng(S317.SEED)  # the A317 draws
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(S317.N_DRAWS)]
    Z, M = S317.raw(T, "clock", 0, T.R)
    mask = B.pair_mask(M, fast, pm, *od_draws[:1]) & B.pair_mask(M, od)
    mask = mask.astype(bool)
    print(f"states ({time.time() - t0:.0f}s): pairs {int(mask.sum())}", flush=True)
    obs = np.array([per_tissue(Z, M, d, pm, mask) for d in od_draws])
    est = np.angle(np.nanmean(obs / np.abs(obs), 0)) * HOURS
    rb = np.random.default_rng([332, 2])
    boot = []
    for i in range(1000):
        p_ = rb.choice(pm, len(pm))
        o_ = rb.choice(od, len(pm))
        boot.append(np.angle(per_tissue(Z, M, o_, p_, mask)) * HOURS)
        if i % 200 == 0:
            print(f"boot {i} ({time.time() - t0:.0f}s)", flush=True)
    boot = np.array(boot)
    dev = np.angle(np.exp(1j * (boot - est) / HOURS)) * HOURS  # deviations from the estimate, wrapped
    lo, hi = est + np.nanquantile(dev, 0.05, 0), est + np.nanquantile(dev, 0.95, 0)
    t315 = pd.read_csv(EXP / "A315_braindeath_gain/output/tissues.csv", index_col=0)
    rows = []
    for j, t in enumerate(H.BODY):
        npair = int(sum(mask[k] for k, (a, b) in enumerate(H.PAIRS) if j in (a, b)))
        L = t315.loc[t] if t in t315.index else None
        rows.append(dict(tissue=t, pairs=npair, phase_change_h=est[j], phase_lo90_h=lo[j], phase_hi90_h=hi[j],
                         L=None if L is None else L.L, L_matched_median=None if L is None else L.matched_median,
                         amplitude_ratio_vs_matched=None if L is None else float(np.exp((L.L - L.matched_median) / 2)),
                         amplitude_ratio_raw=None if L is None else float(np.exp(L.L / 2)),
                         L_lo90=None if L is None else L.lo90, L_hi90=None if L is None else L.hi90))
    o = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    o.to_csv(OUT / "d2_per_organ.csv", index=False)
    print(o.round(3).to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
