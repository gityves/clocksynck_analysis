"""A349 D349.4 (registered, descriptive): the A132b compression ratio after removing the A109 non-temporal axis.

1. Reproduce A132b (`compute_posthoc_shrinkage`, 14 genes, fast vs ventilated post-mortem): median ratio of pole deviations
   (stored 0.585); stop if it does not reproduce.
2. Non-temporal axis: A109 `nontemporal` (leading axis of the time-of-death-residual cross-covariance of fast deaths).
3. Projection (the registration's "host states projected off the leading displacement axis"): per group (fast and
   post-mortem) the host state w_d of each donor (A109 `host_frame` / `host_states` on the group's Hermitian covariance)
   is split along the axis, s_d = Re(w_d e^{-i theta}), and the shared displacement s_d e^{i theta} is subtracted from
   every tissue state of that donor. Pair phases and the same readouts are recomputed on the projected states, the fast
   reference too. Aggregates only.

    python a349_d4.py
Output (`../output/`): d4_compression.csv
"""
import sys

import numpy as np
import pandas as pd

import a349_common as C
from a349_common import B, H
import compute_coupling_nulls as N
import compute_posthoc_shrinkage as PS


def projected_H(Zraw, M, V, gd, theta):
    Zc, m = N._centre(Zraw, M, gd)
    Hm = H.hermitian(Zc, m, V[gd])[0]
    u, _ = H.host_frame(Hm)
    w = H.host_states(Zc, m, u)
    s = np.real(w * np.exp(-1j * theta))
    Zp = (Zc - (s * np.exp(1j * theta))[:, None]) * m
    return H.cross_moments(Zp, m)[0], Zp, m, Zc


def main():
    T, fast, pm = C.tissues()
    Zraw, M, V = N._raw_matrix(T)
    mask = B.pair_mask(M, fast, pm)
    th_f = B.phase(Zraw, M, fast, mask)
    th_p = H.pair_phase(N._H_of(Zraw, M, pm))[mask]
    r0 = PS.readouts(th_f, th_p)
    print("reproduction", r0, flush=True)
    if abs(r0["ratio_median"] - 0.5848837067825516) > 1e-6:
        sys.exit("A132b ratio not reproduced; stopping before the projection")
    nt = N.nontemporal(T)
    nt_h = float(nt[nt.group == "fast"].nt_axis_h.iloc[0])
    th = nt_h / H.HOURS
    Hf, Zf, mf, Zcf = projected_H(Zraw, M, V, fast, th)
    Hp, Zp, mp, Zcp = projected_H(Zraw, M, V, pm, th)
    pf, pp = H.pair_phase(Hf)[mask], H.pair_phase(Hp)[mask]
    r1 = PS.readouts(pf, pp)
    sh = lambda Z, m: H.share_along(H.cross_cov_matrix(Z, m), nt_h)  # noqa: E731
    rows = [dict(stage="A132b reproduction (no projection)", n_pairs=int(mask.sum()), **r0, stored_ratio=0.5848837067825516,
                 share_along_nt_fast=sh(Zcf, mf), share_along_nt_pm=sh(Zcp, mp)),
            dict(stage="axis removed (host component along A109 axis)", n_pairs=int(mask.sum()), **r1,
                 share_along_nt_fast=sh(Zf, mf), share_along_nt_pm=sh(Zp, mp), nt_axis_h=nt_h,
                 nt_eig_share=float(nt[nt.group == "fast"].nt_eig_share.iloc[0]))]
    out = pd.DataFrame(rows)
    C.OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(C.OUT / "d4_compression.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.round(4).T.to_string())


if __name__ == "__main__":
    main()
