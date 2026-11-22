"""A342 scoring (registered in `paper/docs/registered_experiments.md`, A342 and its power note): timed living human skin
(GSE326530), untimed readouts against timed truth, and the sun-exposure contrast. One run.

Units: 20 volunteers x 2 sites (P photoprotected, E photoexposed) x 4 clock times. Spliced = `counts_exon`, unspliced
= `counts_intron` (`counts_intron_direct` as sensitivity); symbols through `clocksync.genes.canonical`. One waveform
(GTEx not-sun-exposed skin) for every statistic; site-specific waveforms as sensitivity (power note). Gene sets: 14
genes ("all") and the core loop. Bootstraps (2,000) resample volunteers with all their samples; permutations 10,000.
- P342.0: per site, |beta_s| (beta = sum z exp(-i 2 pi t / 24), states residualised on volunteer) against times
  permuted within volunteer.
- H342.1: per site, untimed lead (A328 part b, matched estimator: spliced thinned gene by gene to the unspliced level,
  20 replicates averaged; plain and per-gene alongside), volunteer as covariate; gate = |c| / norms against unspliced
  rows permuted within volunteer; pass if gate p < 0.05 (Holm over the two sites) and the 90 % interval in (0, 6 h).
- H342.2: delta_U - delta_T (timed lead on the plain states) within +-1.5 h (90 % paired interval), where H342.1
  passes.
- H342.3: E-P pairing on one projector, states centred within site (no time, no volunteer covariate): gate against P
  rows permuted (10,000), clock-specific above the q95 of 200 matched and 200 responsive companions; offset arg c
  against the timed offset arg(beta_E conj beta_P): pass if the untimed offset is inside the timed 90 % interval or the
  two are < 1 h apart.
- H342.4 (descriptive): untimed log relative amplitude (half `log_var_ratio`) vs timed log(|beta_E| / |beta_P|).
- H342.5 (descriptive): timed log amplitude ratio, its interval and the 5-95 % range of the same ratio on the 200
  matched companions; the timed offset against +-1 h.
Descriptive: unspliced timed amplitude ratio; per-gene timed acrophases; QC sensitivity (runs < 30 % assigned to exons
removed, with their pair partner for the pairing); direct introns (lead); site waveforms; volunteers 1-11 vs 12-20.

    python score_a342.py [--check]
`--check`: unspliced := spliced thinned to 10 %, times permuted within volunteer, P rows permuted across volunteers
(true lead 0, no rhythm, no pairing). Outputs (`../output/score[_check]/`): verdicts.csv, leads.csv, pairing.csv,
genes.csv, sensitivity.csv.
"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import power_a342 as A  # noqa: E402  (sets the import paths)
from power_a342 import P  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import log_var_ratio, projector, pseudo_waveforms, states  # noqa: E402
from common import development as dev  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, N_PERM, N_BOOT, N_COMP, REPS = 3421, 10000, 2000, 200, 20
HOURS = A.HOURS
wrap = P.score_a337.wrap
CHECK = False


def load(name: str) -> pd.DataFrame:
    x = pd.read_csv(A.DATA / f"counts_{name}.tsv.gz", sep="\t")
    x["symbol"] = x.symbol.astype(str).map(canonical)
    return x.drop_duplicates("symbol").set_index("symbol").drop(columns="gene_id")


def cluster_boot(vol, rng, n):
    """n x N bootstrap indices resampling volunteers (requires equal samples per volunteer)."""
    ids = np.unique(vol)
    rows = np.array([np.flatnonzero(vol == v) for v in ids])
    return rows[rng.integers(0, len(ids), (n, len(ids)))].reshape(n, -1)


def zfun(x, depth, Z, wf, L):
    idx = [A.GENES.index(g) for g in wf.genes]
    y = dev.residual(np.log((x[:, idx] + 0.5) / np.broadcast_to(depth, x.shape)[:, idx]), Z)
    b = (y - y.mean(0)) @ L.T
    return b[:, 0] + 1j * b[:, 1]


def timed_b(z, t, vol):
    return A.timed(z, t, vol)


def site_stats(s, u, dd, wf, rng, B) -> dict:
    """Leads (matched, plain, per-gene), timed lead and rhythm of one site."""
    vol, t = dd.vol.to_numpy(), dd.t.to_numpy(float)
    ds, du = dd.depth_s.to_numpy(float)[:, None], dd.depth_u.to_numpy(float)[:, None]
    Z = dev.design(dd.assign(v=dd.vol.astype(str)), age=None, factors=("v",))
    idx = [A.GENES.index(g) for g in wf.genes]
    L = projector(pd.DataFrame(s[:, idx], columns=list(wf.genes)), wf)
    zs, zu = zfun(s, ds, Z, wf, L), zfun(u, du, Z, wf, L)
    zm = np.zeros(len(s), complex)
    zero = np.zeros(len(s), int)
    for _ in range(REPS):
        sm, dsm = P.depth_match(s, u, ds[:, 0], du[:, 0], zero, rng)
        zm += zfun(sm, dsm, Z, wf, L)
    zm /= REPS
    perm = P.score_a337.within(vol, rng, N_PERM)
    bs_t = np.array([timed_b(zs[b], t[b], vol[b]) for b in B])
    bu_t = np.array([timed_b(zu[b], t[b], vol[b]) for b in B])
    beta_s, beta_u = timed_b(zs, t, vol), timed_b(zu, t, vol)
    dT = float(np.angle(beta_u * np.conj(beta_s)) * HOURS)
    dT_draws = np.angle(bu_t * np.conj(bs_t)) * HOURS
    null = np.array([abs(timed_b(zs, t[p], vol)) for p in perm[:N_PERM]])
    out = dict(rhythm_p=float((1 + np.sum(null >= abs(beta_s))) / (1 + N_PERM)), amp_s=abs(beta_s),
               amp_u=abs(beta_u), dT=dT, dT_lo=P.score_a337.interval(dT, dT_draws)[0],
               dT_hi=P.score_a337.interval(dT, dT_draws)[1])
    for est, z_s in (("matched", zm), ("plain", zs)):
        c = np.sum(zu * np.conj(z_s))
        k = abs(c) / np.sqrt(np.sum(abs(zu) ** 2) * np.sum(abs(z_s) ** 2))
        kp = np.abs((zu[perm] * np.conj(z_s)[None]).sum(1)) / np.sqrt(np.sum(abs(zu) ** 2) * np.sum(abs(z_s) ** 2))
        d = float(np.angle(c) * HOURS)
        dr = np.angle((zu[B] * np.conj(z_s[B])).sum(1)) * HOURS
        lo, hi = P.score_a337.interval(d, dr)
        dlo, dhi = P.score_a337.interval(wrap(d - dT), wrap(dr - dT_draws))
        out.update({f"{est}_kappa": k, f"{est}_gate_p": float((1 + np.sum(kp >= k)) / (1 + N_PERM)),
                    f"{est}_delta": d, f"{est}_lo": lo, f"{est}_hi": hi, f"{est}_diff_lo": dlo, f"{est}_diff_hi": dhi})
    # per-gene estimator (A328 part b companion), volunteer-cluster bootstrap
    P.boot = lambda groups, r, n: cluster_boot(groups, r, n)
    P.N_PERM, P.N_BOOT = 2000, N_BOOT

    def Y(x, depth):
        y = dev.residual(np.log((x[:, idx] + 0.5) / np.broadcast_to(depth, x.shape)[:, idx]), Z)
        return y - y.mean(0)
    d, kap, gp, draws = P.pergene(Y(s, ds), Y(u, du), s, idx, wf, vol, rng)
    lo, hi = P.score_a337.interval(d, draws)
    out.update(pergene_kappa=kap, pergene_gate_p=gp, pergene_delta=d, pergene_lo=lo, pergene_hi=hi)
    return out


def pair_states(sP, sE, dP, dE, wf):
    idx = [A.GENES.index(g) for g in wf.genes]
    g = list(wf.genes)
    L = projector(pd.DataFrame(np.vstack([sP, sE])[:, idx], columns=g), wf)
    st = {}
    for site, x, dd in (("P", sP, dP), ("E", sE, dE)):
        b, v = states(pd.DataFrame(x[:, idx], columns=g), L, wf, dd.depth_s.to_numpy(float))
        st[site] = (b - b.mean(0), v)
    zP = st["P"][0][:, 0] + 1j * st["P"][0][:, 1]
    zE = st["E"][0][:, 0] + 1j * st["E"][0][:, 1]
    return zP, zE, st


def coherence(zE, zP):
    return abs(np.sum(zE * np.conj(zP))) / np.sqrt(np.sum(abs(zE) ** 2) * np.sum(abs(zP) ** 2))


def pairing(sP, sE, dP, dE, wf, rng, B, comps) -> dict:
    vol, t = dP.vol.to_numpy(), dP.t.to_numpy(float)
    zP, zE, st = pair_states(sP, sE, dP, dE, wf)
    k = coherence(zE, zP)
    kp = np.array([coherence(zE, zP[rng.permutation(len(zP))]) for _ in range(N_PERM)])
    c = np.sum(zE * np.conj(zP))
    off_u = float(np.angle(c) * HOURS)
    amp_u = 0.5 * log_var_ratio(st["E"][0], st["E"][1], st["P"][0], st["P"][1])
    bP, bE = timed_b(zP, t, vol), timed_b(zE, t, vol)
    off_t, amp_t = float(np.angle(bE * np.conj(bP)) * HOURS), float(np.log(abs(bE) / abs(bP)))
    bo = []
    for b in B:
        tP, tE = timed_b(zP[b], t[b], vol[b]), timed_b(zE[b], t[b], vol[b])
        bo.append((np.angle(np.sum(zE[b] * np.conj(zP[b]))) * HOURS,
                   0.5 * log_var_ratio(st["E"][0][b], st["E"][1][b], st["P"][0][b], st["P"][1][b]),
                   np.angle(tE * np.conj(tP)) * HOURS, np.log(abs(tE) / abs(tP))))
    bo = np.array(bo)
    ci = P.score_a337.interval
    out = dict(coherence=k, gate_p=float((1 + np.sum(kp >= k)) / (1 + N_PERM)), off_u=off_u, off_t=off_t,
               amp_u=amp_u, amp_t=amp_t)
    out["off_u_lo"], out["off_u_hi"] = ci(off_u, bo[:, 0])
    out["off_t_lo"], out["off_t_hi"] = ci(off_t, bo[:, 2])
    out["amp_u_lo"], out["amp_u_hi"] = np.quantile(bo[:, 1], [.05, .95])
    out["amp_t_lo"], out["amp_t_hi"] = np.quantile(bo[:, 3], [.05, .95])
    out["amp_diff_lo"], out["amp_diff_hi"] = np.quantile(bo[:, 1] - bo[:, 3], [.05, .95])
    for kind, ws in comps.items():
        kc, rc = [], []
        for w in ws:
            gi = list(w.genes)
            Lc = projector(pd.DataFrame(np.vstack([CP[gi].to_numpy(), CE[gi].to_numpy()]), columns=gi), w)
            bPc, _ = states(CP[gi], Lc, w, dP.depth_s.to_numpy(float))
            bEc, _ = states(CE[gi], Lc, w, dE.depth_s.to_numpy(float))
            zPc, zEc = (bPc - bPc.mean(0)) @ [1, 1j], (bEc - bEc.mean(0)) @ [1, 1j]
            kc.append(coherence(zEc, zPc))
            rc.append(np.log(abs(timed_b(zEc, t, vol)) / abs(timed_b(zPc, t, vol))))
        out[f"{kind}_q95"] = float(np.quantile(kc, .95))
        if kind == "matched":
            out["amp_t_comp_q05"], out["amp_t_comp_q95"] = np.quantile(rc, [.05, .95])
    out["clock_specific"] = bool(k > out["matched_q95"] and k > out["responsive_q95"])
    return out


CP = CE = GENOME = LOGMEAN = None


def run(variant, s_all, u_all, d, rng, waveform_for_site=None, pair=True) -> tuple[list, list]:
    """All statistics of one variant. s_all, u_all: run x gene (14 genes, A.GENES order) DataFrames."""
    global CP, CE
    leads, pairs = [], []
    data = {}
    for site in ("P", "E"):
        dd = d[d.site == site].reset_index(drop=True)
        data[site] = (dd, s_all.loc[dd.run, A.GENES].to_numpy(np.int64), u_all.loc[dd.run, A.GENES].to_numpy(np.int64))
    vol = data["P"][0].vol.to_numpy()
    equal = len(set(np.bincount(vol)[np.unique(vol)])) == 1
    B = cluster_boot(vol, rng, N_BOOT) if equal else P.score_a337.boot(vol, rng, N_BOOT)
    for w, wf0 in A.WAVEFORMS.items():
        for site in ("P", "E"):
            wf = wf0 if waveform_for_site is None else waveform_for_site[site][w]
            dd, s, u = data[site]
            leads.append(dict(variant=variant, site=site, waveform=w, n=len(dd), **site_stats(s, u, dd, wf, rng, B)))
        if pair:
            dP, sP, _ = data["P"]
            dE, sE, _ = data["E"]
            assert (dP.vol.to_numpy() == dE.vol.to_numpy()).all() and (dP.t.to_numpy() == dE.t.to_numpy()).all()
            CP, CE = GENOME.loc[dP.run].reset_index(drop=True), GENOME.loc[dE.run].reset_index(drop=True)
            comps = {kind: pseudo_waveforms(LOGMEAN, N_COMP, np.random.default_rng([SEED, 7, len(w)]), wf0, kind)
                     for kind in ("matched", "responsive")}
            pairs.append(dict(variant=variant, waveform=w, n=len(dP), **pairing(sP, sE, dP, dE, wf0, rng, B, comps)))
    return leads, pairs


def main():
    global CHECK, GENOME, LOGMEAN
    CHECK = "--check" in sys.argv
    out = OUT / ("score_check" if CHECK else "score")
    out.mkdir(parents=True, exist_ok=True)
    limit_threads(4)
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    d = A.design()
    ex, it, itd = load("exon"), load("intron"), load("intron_direct")
    q = pd.read_csv(A.DATA / "qc.tsv", sep="\t").set_index("run")
    s_all, u_all, ud_all = ex[d.run].T, it[d.run].T, itd[d.run].T
    if CHECK:
        u_all = pd.DataFrame(rng.binomial(s_all.to_numpy(np.int64), 0.1), index=s_all.index, columns=s_all.columns)
        ud_all = u_all
        for v in d.vol.unique():  # one permutation of times per volunteer, applied to both sites
            perm = dict(zip(np.array(A.TIMES), rng.permutation(np.array(A.TIMES))))
            i = d.index[d.vol == v]
            d.loc[i, "t"] = d.loc[i, "t"].map(perm)
        ip = d.index[d.site == "P"]
        d.loc[ip, "run"] = rng.permutation(d.loc[ip, "run"].to_numpy())
        d = d.sort_values(["site", "vol", "t"]).reset_index(drop=True)
    keep = s_all.sum(0) > 0
    GENOME = s_all.loc[:, keep]
    LOGMEAN = np.log(GENOME.div(GENOME.sum(1), axis=0).mean(0) + 1e-12)
    LOGMEAN = LOGMEAN[GENOME.mean(0) >= 5]
    print(f"load {time.time() - t0:.0f}s", flush=True)
    leads, pairs = run("main", s_all, u_all, d, rng)
    print(f"main {time.time() - t0:.0f}s", flush=True)
    lt, pt = pd.DataFrame(leads), pd.DataFrame(pairs)
    lt.to_csv(out / "leads.csv", index=False)
    pt.to_csv(out / "pairing.csv", index=False)
    # verdicts
    v = []
    m = lt[lt.waveform == "all"].set_index("site")
    holm = dict(zip(m.index, dev.holm(m.matched_gate_p.to_numpy())))
    for w in ("all", "core"):
        for site in ("P", "E"):
            r = lt[(lt.waveform == w) & (lt.site == site)].iloc[0]
            h1 = (holm[site] if w == "all" else r.matched_gate_p) < .05 and r.matched_lo > 0 and r.matched_hi < 6
            v.append(dict(test="P342.0", waveform=w, site=site, value=r.rhythm_p, passed=r.rhythm_p < .05))
            v.append(dict(test="H342.1", waveform=w, site=site, value=r.matched_delta, passed=bool(h1)))
            v.append(dict(test="H342.2", waveform=w, site=site, value=wrap(r.matched_delta - r.dT),
                          passed=bool(h1 and r.matched_diff_lo > -1.5 and r.matched_diff_hi < 1.5)))
        r = pt[pt.waveform == w].iloc[0]
        inside = r.off_t_lo <= r.off_u <= r.off_t_hi if r.off_t_lo <= r.off_t_hi else True
        h3 = r.gate_p < .05 and r.clock_specific and (inside or abs(wrap(r.off_u - r.off_t)) < 1)
        v.append(dict(test="H342.3", waveform=w, site="E-P", value=wrap(r.off_u - r.off_t), passed=bool(h3)))
        v.append(dict(test="H342.4 (descr.)", waveform=w, site="E-P", value=r.amp_u - r.amp_t,
                      passed=None))
        v.append(dict(test="H342.5 (descr.)", waveform=w, site="E-P", value=r.amp_t, passed=None))
    pd.DataFrame(v).to_csv(out / "verdicts.csv", index=False)
    # per-gene timed acrophases (descriptive)
    g = []
    for site in ("P", "E"):
        dd = d[d.site == site].reset_index(drop=True)
        x = s_all.loc[dd.run, A.GENES].to_numpy(float)
        y = np.log((x + 0.5) / dd.depth_s.to_numpy(float)[:, None])
        y = y - pd.DataFrame(y).groupby(dd.vol.to_numpy()).transform("mean").to_numpy()
        e = np.exp(-1j * dd.t.to_numpy(float) / HOURS)
        for j, gene in enumerate(A.GENES):
            b = 2 * np.mean(y[:, j] * e)
            g.append(dict(site=site, gene=gene, amplitude=abs(b), acrophase_h=(-np.angle(b) * HOURS) % 24))
    pd.DataFrame(g).to_csv(out / "genes.csv", index=False)
    # sensitivities
    sens = []
    if not CHECK:
        low = set(q.index[q.assigned_exon_pct < 30])
        bad = set(d.loc[d.run.isin(low), ["vol", "t"]].itertuples(index=False, name=None))
        dq = d[[(a, b) not in bad for a, b in zip(d.vol, d.t)]].reset_index(drop=True)
        for name, args in (("qc_exon30", (s_all, u_all, dq)), ("direct_intron", (s_all, ud_all, d))):
            l2, p2 = run(name, *args, rng, pair=name == "qc_exon30")
            sens += l2 + p2
        site_wf = {site: {"all": cs.gtex(A.GTEX_SITE[site]),
                          "core": cs.gtex(A.GTEX_SITE[site]).drop(["PER1", "PER2", "NR1D1"])} for site in ("P", "E")}
        l2, _ = run("site_waveforms", s_all, u_all, d, rng, site_wf, pair=False)
        sens += l2
        for name, vols in (("vol_1_11", range(1, 12)), ("vol_12_20", range(12, 21))):
            _, p2 = run(name, s_all, u_all, d[d.vol.isin(list(vols))].reset_index(drop=True), rng)
            sens += p2
    pd.DataFrame(sens).to_csv(out / "sensitivity.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(v).round(3).to_string(index=False))
    print(lt.round(3).T.to_string(), pt.round(3).T.to_string(), f"runtime {time.time() - t0:.0f}s", sep="\n",
          flush=True)


if __name__ == "__main__":
    main()
