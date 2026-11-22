"""A344 scoring (registered in `paper/docs/registered_experiments.md`, A344): untimed readouts of timed skin
(GSE326530) on time subsets against the full-cycle timed truth of A342. One run.

Subsets: D3 = {6, 12, 18} (primary), D2 = {6, 12}, N2 = {0, 6}, full. Spliced `counts_exon`, unspliced
`counts_intron_direct` (A342 erratum). One projector over both sites per subset; the GTEx not-sun-exposed skin
waveform; 14 genes and core loop. Bootstraps (2,000) resample the 20 volunteers jointly for the subset statistic and the
full-cycle truth (paired intervals).
- Offset E-P: arg sum z_E conj(z_P), states centred within site over the subset. Truth: full-cycle timed offset.
- Lead per site: matched estimator (20 replicates), volunteer covariate within the subset. Truth: full-cycle timed
  lead (plain states, A342).
- Log relative amplitude E/P: half `log_var_ratio`. Truth: full-cycle timed log ratio.
- Per-visit phase: angle of (z_P + z_E) / 2 in hours minus clock time, one common rotation removed (circular mean);
  concentration R = |mean exp(i (angle - t))|, median absolute error (h), p of R against times permuted within
  volunteer (2,000).
H344.1: D3 offset difference interval within +-1 h. H344.2: D3 lead difference intervals within +-1.5 h, per site
(both interval ends inside the bound). H344.3: the same on D2. Undecidable when the interval is wider
than the bound (half-width > bound): reported as such, not as a fail.

    python score_a344.py
Outputs (`../output/score_a344/`): subsets.csv, verdicts.csv.
"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import power_a342 as A  # noqa: E402
import score_a342 as S  # noqa: E402
from power_a342 import P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import log_var_ratio, projector  # noqa: E402
from common import development as dev  # noqa: E402

OUT = HERE.parents[1] / "output" / "score_a344"
SEED, N_BOOT, N_PERM, REPS = 3441, 2000, 2000, 20
HOURS = A.HOURS
wrap = P.score_a337.wrap
SUBSETS = {"full": (0, 6, 12, 18), "D3": (6, 12, 18), "D2": (6, 12), "N2": (0, 6)}


def rows_for(dd, vols):
    """Row indices of dd for a resampled list of volunteers (each volunteer's rows in time order)."""
    pos = {v: np.flatnonzero(dd.vol.to_numpy() == v) for v in np.unique(dd.vol)}
    return np.concatenate([pos[v] for v in vols])


def site_lead(s, u, dd, wf, rng, VB):
    """Matched untimed lead (point, bootstrap draws) and plain timed lead on dd."""
    vol, t = dd.vol.to_numpy(), dd.t.to_numpy(float)
    ds, du = dd.depth_s.to_numpy(float)[:, None], dd.depth_u.to_numpy(float)[:, None]
    Z = dev.design(dd.assign(v=dd.vol.astype(str)), age=None, factors=("v",))
    idx = [A.GENES.index(g) for g in wf.genes]
    L = projector(pd.DataFrame(s[:, idx], columns=list(wf.genes)), wf)
    zs, zu = S.zfun(s, ds, Z, wf, L), S.zfun(u, du, Z, wf, L)
    zm = np.zeros(len(s), complex)
    for _ in range(REPS):
        sm, dsm = P.depth_match(s, u, ds[:, 0], du[:, 0], np.zeros(len(s), int), rng)
        zm += S.zfun(sm, dsm, Z, wf, L)
    zm /= REPS
    dU = float(np.angle(np.sum(zu * np.conj(zm))) * HOURS)
    dT = float(np.angle(A.timed(zu, t, vol) * np.conj(A.timed(zs, t, vol))) * HOURS)
    dU_b, dT_b = [], []
    for vb in VB:
        b = rows_for(dd, vb)
        dU_b.append(np.angle(np.sum(zu[b] * np.conj(zm[b]))) * HOURS)
        dT_b.append(np.angle(A.timed(zu[b], t[b], vol[b]) * np.conj(A.timed(zs[b], t[b], vol[b]))) * HOURS)
    return dU, np.array(dU_b), dT, np.array(dT_b)


def pair(sP, sE, dP, dE, wf, VB, timed: bool):
    zP, zE, st = S.pair_states(sP, sE, dP, dE, wf)
    vol, t = dP.vol.to_numpy(), dP.t.to_numpy(float)

    def stats(b):
        off_u = np.angle(np.sum(zE[b] * np.conj(zP[b]))) * HOURS
        amp_u = 0.5 * log_var_ratio(st["E"][0][b], st["E"][1][b], st["P"][0][b], st["P"][1][b])
        if not timed:
            return off_u, amp_u, np.nan, np.nan
        tP, tE = A.timed(zP[b], t[b], vol[b]), A.timed(zE[b], t[b], vol[b])
        return off_u, amp_u, np.angle(tE * np.conj(tP)) * HOURS, np.log(abs(tE) / abs(tP))
    point = stats(np.arange(len(zP)))
    draws = np.array([stats(rows_for(dP, vb)) for vb in VB])
    return point, draws, zP, zE


def visit_phase(zP, zE, dP, rng):
    z = (zP + zE) / 2
    t = dP.t.to_numpy(float)
    vol = dP.vol.to_numpy()
    d = np.angle(z) - t / HOURS
    R = abs(np.mean(np.exp(1j * d)))
    rot = np.angle(np.mean(np.exp(1j * d)))
    err = np.abs(wrap((d - rot) * HOURS))
    perm = P.score_a337.within(vol, rng, N_PERM)
    Rn = np.array([abs(np.mean(np.exp(1j * (np.angle(z) - t[p] / HOURS)))) for p in perm])
    return dict(phase_R=R, phase_median_err_h=float(np.median(err)), phase_p=float((1 + np.sum(Rn >= R)) / (1 + N_PERM)))


def interval(x):
    return float(np.quantile(x, .05)), float(np.quantile(x, .95))


def main():
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    d = A.design()
    ex, itd = S.load("exon"), S.load("intron_direct")
    s_all, u_all = ex[d.run].T, itd[d.run].T
    S.GENOME = s_all
    vols = np.unique(d.vol)
    VB = [rng.choice(vols, len(vols)) for _ in range(N_BOOT)]
    rows = []
    for w, wf in A.WAVEFORMS.items():
        truth = {}
        for name, times in SUBSETS.items():
            dd = d[d.t.isin(times)].reset_index(drop=True)
            data = {}
            for site in ("P", "E"):
                x = dd[dd.site == site].reset_index(drop=True)
                data[site] = (x, s_all.loc[x.run, A.GENES].to_numpy(np.int64), u_all.loc[x.run, A.GENES].to_numpy(np.int64))
            dP, sP, _ = data["P"]
            dE, sE, _ = data["E"]
            assert (dP.vol.to_numpy() == dE.vol.to_numpy()).all() and (dP.t.to_numpy() == dE.t.to_numpy()).all()
            (off_u, amp_u, off_t, amp_t), dr, zP, zE = pair(sP, sE, dP, dE, wf, VB, timed=name == "full")
            leads = {site: site_lead(data[site][1], data[site][2], data[site][0], wf, rng, VB) for site in ("P", "E")}
            if name == "full":
                truth = dict(off=(off_t, dr[:, 2]), amp=(amp_t, dr[:, 3]),
                             **{f"lead_{site}": (leads[site][2], leads[site][3]) for site in ("P", "E")})
            r = dict(waveform=w, subset=name, visits=len(dP), off_u=off_u, off_truth=truth["off"][0],
                     amp_u=amp_u, amp_truth=truth["amp"][0], **visit_phase(zP, zE, dP, rng))
            r["off_diff"] = float(wrap(off_u - truth["off"][0]))
            r["off_diff_lo"], r["off_diff_hi"] = interval(wrap(dr[:, 0] - truth["off"][1]))
            r["amp_diff"] = amp_u - truth["amp"][0]
            r["amp_diff_lo"], r["amp_diff_hi"] = interval(dr[:, 1] - truth["amp"][1])
            for site in ("P", "E"):
                dU, dU_b, _, _ = leads[site]
                tv, tb = truth[f"lead_{site}"]
                r[f"lead_{site}"], r[f"lead_{site}_truth"] = dU, tv
                r[f"lead_{site}_diff"] = float(wrap(dU - tv))
                r[f"lead_{site}_diff_lo"], r[f"lead_{site}_diff_hi"] = interval(wrap(dU_b - tb))
            rows.append(r)
            print(w, name, f"{time.time() - t0:.0f}s", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "subsets.csv", index=False)
    v = []
    for w in A.WAVEFORMS:
        for name, test in (("D3", "H344"), ("D2", "H344.3")):
            r = t[(t.waveform == w) & (t.subset == name)].iloc[0]
            hw = (r.off_diff_hi - r.off_diff_lo) / 2
            v.append(dict(test=f"{test}.1 offset" if name == "D3" else f"{test} offset", waveform=w, subset=name,
                          diff=r.off_diff, lo=r.off_diff_lo, hi=r.off_diff_hi,
                          verdict="undecidable" if hw > 1 else "pass" if (r.off_diff_lo > -1 and r.off_diff_hi < 1)
                          else "fail"))
            for site in ("P", "E"):
                lo, hi = r[f"lead_{site}_diff_lo"], r[f"lead_{site}_diff_hi"]
                v.append(dict(test=f"{test}.2 lead {site}" if name == "D3" else f"{test} lead {site}", waveform=w,
                              subset=name, diff=r[f"lead_{site}_diff"], lo=lo, hi=hi,
                              verdict="undecidable" if (hi - lo) / 2 > 1.5 else "pass" if (lo > -1.5 and hi < 1.5)
                              else "fail"))
    v = pd.DataFrame(v)
    v.to_csv(OUT / "verdicts.csv", index=False)
    pd.set_option("display.width", 250)
    print(v.round(3).to_string(index=False), t.round(3).T.to_string(), f"runtime {time.time() - t0:.0f}s", sep="\n")


if __name__ == "__main__":
    main()
