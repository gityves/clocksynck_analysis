"""A342 power stage (registered in `paper/docs/registered_experiments.md`, A342), before any clock-gene count of
GSE326530 is read.

Read from the delivery: per run the exon and intron depths (`qc.tsv`: assigned_exon_n; intronic_fraction x
assigned_genebody_n, the genome-wide subtraction-intron total) and the design (volunteer, site, time; `samples.tsv`). No gene count of the set is read.

Stand-ins: GTEx donors with both skin tissues and a time of death within +-2 h (fast deaths, Hardy 1-2) of 0, 6, 12 or 18 h
(time of death in memory only). One draw = 20 donors per time bin, grouped into 20 "volunteers" of four (one per bin);
not-sun-exposed skin stands in for the photoprotected site (P), sun-exposed for the photoexposed (E). Clock times are
the bin centres. Spliced counts = GTEx clock counts (clock at relative amplitude a, `simulate.scale_clock`); E may be
turned and scaled (`simulate.clock_change`); unspliced = spliced with the state turned by the lead
(`simulate.thin_clock`). Both thinned to the delivery's per-run targets: GTEx median skin fraction of each gene x the
run's exon depth (spliced) or intron depth (unspliced; same per-gene fractions, i.e. the clock genes are assumed to carry
the genome-wide intron share), matched to runs by volunteer x site x time.
Settings: main (a = 1, lead 1.27 h, E unchanged), a0.75, lead0 (no turn), shifted (lead 1.27 h, E turned +1 h and
scaled 0.8). One waveform and one projector for both sites (not-sun-exposed skin; see the power note).

Statistics per draw (14 genes and core loop): per site the A328 part b lead (matched estimator for all draws, plain and
per-gene on the first 200), the timed lead and the rhythm p (500 permutations of time within volunteer); the untimed
E-P offset, gate and log relative amplitude (`clocksync.paired.log_var_ratio`, halved) and their timed counterparts,
with 90 % intervals from 500 volunteer bootstraps (paired for the differences).

    python power_a342.py [--n-sub 1000] [--jobs 4] [--check]
Outputs (`../output/power/`): draws.csv, rates.csv (`_check` with --check). Controlled access: aggregates only.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "A328_fetal_running_loop" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import power_a328b as P  # noqa: E402  (sets the import paths of the analysis repo)
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import log_var_ratio, projector, states  # noqa: E402
from clocksync.simulate import clock_change, scale_clock, thin_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402

DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/sccc_mapping/gse326530_skin")
OUT = HERE.parents[1] / "output" / "power"
SEED, N_PERM, N_BOOT = 3420, 500, 500
HOURS = 24 / (2 * np.pi)
GTEX_SITE = {"P": "Skin - Not Sun Exposed (Suprapubic)", "E": "Skin - Sun Exposed (Lower leg)"}
WF = cs.gtex(GTEX_SITE["P"])
GENES = list(WF.genes)
assert GENES == P.GENES
WAVEFORMS = {"all": WF, "core": WF.drop(["PER1", "PER2", "NR1D1"])}
TIMES = (0, 6, 12, 18)
SETTINGS = {"main": (1.0, 1.27, 0.0, 1.0), "a0.75": (0.75, 1.27, 0.0, 1.0), "lead0": (1.0, 0.0, 0.0, 1.0),
            "shifted": (1.0, 1.27, 1.0, 0.8)}  # (clock scale, lead h, E turn h, E scale)
S = {}


def design() -> pd.DataFrame:
    s = pd.read_csv(DATA / "samples.tsv", sep="\t")
    q = pd.read_csv(DATA / "qc.tsv", sep="\t").set_index("run")
    t = s.geo_title.str.extract(r"Volunteer (\d+), (\w+) skin, .*\(T(\d\d)\)")
    d = pd.DataFrame(dict(run=s.run_accession, vol=t[0].astype(int), site=t[1].map({"photoprotected": "P",
                     "photoexposed": "E"}), t=t[2].astype(int)))
    d["depth_s"] = q.loc[d.run, "assigned_exon_n"].to_numpy(float)
    d["depth_u"] = (q.loc[d.run, "intronic_fraction"] * q.loc[d.run, "assigned_genebody_n"]).to_numpy(float)
    assert d.notna().all().all() and len(d) == 160
    return d.sort_values(["site", "vol", "t"]).reset_index(drop=True)


def pools(rng):
    pb, don = gtex.load(list(GTEX_SITE.values()), hardy=gtex.FAST, genes="clock", depth="library_size")
    tod = don.tod.to_numpy(float) % 24
    dist = np.abs(((tod[:, None] - np.array(TIMES)[None]) + 12) % 24 - 12)
    bin_ = np.where(dist.min(1) <= 2.0, np.array(TIMES)[dist.argmin(1)], -1)
    tbin = pd.Series(bin_, index=don.index)
    both = pb.groupby("donor").compartment.nunique()
    keep = [d for d in both.index[both == 2] if tbin.get(d, -1) >= 0]
    out = {"tbin": tbin.loc[keep].to_numpy(), "donors": np.array(keep)}
    for site, tissue in GTEX_SITE.items():
        c = pb[pb.compartment == tissue].set_index("donor").loc[keep].reset_index()
        out[site] = dict(depth=c.depth.to_numpy(float), xs={}, xu={})
        for name, (a, lead, turn, sc) in SETTINGS.items():
            s = scale_clock(c[["depth", *GENES]], WF, a, depth=c.depth, rng=rng)
            beta, _ = states(s, projector(s, WF), WF, c.depth.to_numpy(float))
            x = np.round(s[GENES].to_numpy(float))
            if site == "E" and (turn or sc != 1):
                x = thin_clock(x, WF.X, clock_change(beta - beta.mean(0), sc, turn), rng)
                xd = pd.DataFrame(x, columns=GENES)
                beta, _ = states(xd, projector(xd, WF), WF, c.depth.to_numpy(float))
            out[site]["xs"][name] = x.astype(np.int64)
            out[site]["xu"][name] = thin_clock(x, WF.X, clock_change(beta - beta.mean(0), 1.0, lead), rng).astype(np.int64)
    counts = pd.Series(out["tbin"]).value_counts().sort_index()
    return out, pd.DataFrame(dict(time=counts.index, donors=counts.to_numpy()))


def zstates(x, depth, Z, wf):
    idx = [GENES.index(g) for g in wf.genes]
    L = projector(pd.DataFrame(x[:, idx], columns=list(wf.genes)), wf)
    y = dev.residual(np.log((x[:, idx] + 0.5) / np.broadcast_to(depth, x.shape)[:, idx]), Z)
    b = (y - y.mean(0)) @ L.T
    return b[:, 0] + 1j * b[:, 1], L, idx


def vboot(vol, rng, n):
    """n volunteer-cluster bootstrap index arrays (rows of resampled volunteers, all their samples)."""
    ids = np.unique(vol)
    rows = {v: np.flatnonzero(vol == v) for v in ids}
    return [np.concatenate([rows[v] for v in rng.choice(ids, len(ids))]) for _ in range(n)]


def timed(z, t, vol):
    """beta = sum z exp(-i 2 pi t / 24), z centred within volunteer."""
    zc = z - pd.Series(z).groupby(vol).transform("mean").to_numpy()
    return np.sum(zc * np.exp(-1j * t / HOURS))


def interval(est, draws, circ=True):
    d = (P.score_a337.wrap(np.asarray(draws) - est) + est) if circ else np.asarray(draws)
    return float(np.quantile(d, .05)), float(np.quantile(d, .95))


def one(task) -> list:
    i, = task
    rng = np.random.default_rng([SEED, i])
    d, pool = S["design"], S["pools"]
    pick = {T: rng.choice(np.flatnonzero(pool["tbin"] == T), 20, replace=False) for T in TIMES}
    out = []
    for name in SETTINGS:
        data = {}
        for site in ("P", "E"):
            dd = d[d.site == site].reset_index(drop=True)
            rows = np.array([pick[T][v - 1] for v, T in zip(dd.vol, dd.t)])
            depth = pool[site]["depth"][rows]
            ts = np.outer(dd.depth_s.to_numpy(), S["frac"])
            tu = np.outer(dd.depth_u.to_numpy(), S["frac"])
            s, _ = P.thin(pool[site]["xs"][name][rows], ts, depth, rng)
            u, _ = P.thin(pool[site]["xu"][name][rows], tu, depth, rng)
            data[site] = (dd, s, u)
        vol, t = data["P"][0].vol.to_numpy(), data["P"][0].t.to_numpy(float)
        B = vboot(vol, rng, N_BOOT)
        for w, wf in WAVEFORMS.items():
            row = dict(sub=i, setting=name, waveform=w)
            zs = {}
            for site in ("P", "E"):
                dd, s, u = data[site]
                obs = dd.assign(cohort="all")
                Z = dev.design(obs.assign(v=obs.vol.astype(str)), age=None, factors=("v",))
                ests = ("matched", "plain", "pergene") if i < 200 else ("matched",)
                for est in ests:
                    r = P.stats(s, u, obs, Z, np.zeros(len(obs), int), wf, rng, est)
                    row.update({f"{site}_{est}_{k}": r[k] for k in ("gate_p", "delta", "lo", "hi")})
                z_s, _, _ = zstates(s, dd.depth_s.to_numpy(float)[:, None], Z, wf)
                z_u, _, _ = zstates(u, dd.depth_u.to_numpy(float)[:, None], Z, wf)
                bs, bu = timed(z_s, t, vol), timed(z_u, t, vol)
                row[f"{site}_dT"] = float(np.angle(bu * np.conj(bs)) * HOURS)
                perm = P.score_a337.within(vol, rng, N_PERM)
                null = np.array([abs(timed(z_s, t[p], vol)) for p in perm])
                row[f"{site}_rhythm_p"] = float((1 + np.sum(null >= abs(bs))) / (1 + N_PERM))
                zs[site] = z_s
            # pairing E-P on one projector, states centred within site (no volunteer covariate)
            sP, sE = data["P"][1], data["E"][1]
            idx = [GENES.index(g) for g in wf.genes]
            L = projector(pd.DataFrame(np.vstack([sP, sE])[:, idx], columns=list(wf.genes)), wf)
            st = {}
            for site, x in (("P", sP), ("E", sE)):
                dep = data[site][0].depth_s.to_numpy(float)
                b, v = states(pd.DataFrame(x[:, idx], columns=list(wf.genes)), L, wf, dep)
                st[site] = (b - b.mean(0), v)
            zP = st["P"][0][:, 0] + 1j * st["P"][0][:, 1]
            zE = st["E"][0][:, 0] + 1j * st["E"][0][:, 1]
            c = np.sum(zE * np.conj(zP))
            k = abs(c) / np.sqrt(np.sum(abs(zE) ** 2) * np.sum(abs(zP) ** 2))
            kp = [abs(np.sum(zE * np.conj(zP[rng.permutation(len(zP))]))) / np.sqrt(np.sum(abs(zE) ** 2)
                  * np.sum(abs(zP) ** 2)) for _ in range(N_PERM)]
            off_u = float(np.angle(c) * HOURS)
            amp_u = 0.5 * log_var_ratio(st["E"][0], st["E"][1], st["P"][0], st["P"][1])
            bP, bE = timed(zP, t, vol), timed(zE, t, vol)
            off_t = float(np.angle(bE * np.conj(bP)) * HOURS)
            amp_t = float(np.log(abs(bE) / abs(bP)))
            bo = []
            for b in B:
                cb = np.sum(zE[b] * np.conj(zP[b]))
                au = 0.5 * log_var_ratio(st["E"][0][b], st["E"][1][b], st["P"][0][b], st["P"][1][b])
                tP, tE = timed(zP[b], t[b], vol[b]), timed(zE[b], t[b], vol[b])
                bo.append((np.angle(cb) * HOURS, au, np.angle(tE * np.conj(tP)) * HOURS, np.log(abs(tE) / abs(tP))))
            bo = np.array(bo)
            row.update(pair_gate_p=float((1 + np.sum(np.array(kp) >= k)) / (1 + N_PERM)), off_u=off_u, off_t=off_t,
                       amp_u=amp_u, amp_t=amp_t)
            row["off_u_lo"], row["off_u_hi"] = interval(off_u, bo[:, 0])
            row["off_t_lo"], row["off_t_hi"] = interval(off_t, bo[:, 2])
            row["off_diff_lo"], row["off_diff_hi"] = interval(P.score_a337.wrap(off_u - off_t),
                                                              P.score_a337.wrap(bo[:, 0] - bo[:, 2]))
            row["amp_u_lo"], row["amp_u_hi"] = interval(amp_u, bo[:, 1], False)
            row["amp_diff_lo"], row["amp_diff_hi"] = interval(amp_u - amp_t, bo[:, 1] - bo[:, 3], False)
            row["amp_t_se"] = float(np.std(bo[:, 3]))
            out.append(row)
    return out


def rates(res) -> pd.DataFrame:
    rows = []
    for (st, w), r in res.groupby(["setting", "waveform"], sort=False):
        row = dict(setting=st, waveform=w, draws=len(r))
        for site in ("P", "E"):
            row[f"{site}_rhythm"] = (r[f"{site}_rhythm_p"] < .05).mean()
            for est in ("matched", "plain", "pergene"):
                g = r.dropna(subset=[f"{site}_{est}_delta"]) if f"{site}_{est}_delta" in r else r.iloc[:0]
                if len(g):
                    ok = (g[f"{site}_{est}_gate_p"] < .05) & (g[f"{site}_{est}_lo"] > 0) & (g[f"{site}_{est}_hi"] < 6)
                    row[f"{site}_{est}_H1"] = ok.mean()
                    row[f"{site}_{est}_width"] = (g[f"{site}_{est}_hi"] - g[f"{site}_{est}_lo"]).median()
                    row[f"{site}_{est}_agree"] = (abs(P.score_a337.wrap(g[f"{site}_{est}_delta"] - g[f"{site}_dT"]))
                                                  < 1.5).mean()
        row.update(pair_gate=(r.pair_gate_p < .05).mean(),
                   off_halfwidth=((r.off_u_hi - r.off_u_lo) / 2).median(),
                   off_agree=((r.off_diff_lo > -1) & (r.off_diff_hi < 1)).mean(),
                   amp_halfwidth=((r.amp_u_hi - r.amp_u_lo) / 2).median(),
                   amp_agree=((r.amp_diff_lo > -.15) & (r.amp_diff_hi < .15)).mean(),
                   amp_t_se=r.amp_t_se.median(), off_t_median=r.off_t.median(), amp_t_median=r.amp_t.median())
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0, sfx = time.time(), "_check" if a.check else ""
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng([SEED, 0])
    pool, counts = pools(rng)
    print("GTEx donors with both skins per time bin:", counts.to_dict("list"), flush=True)
    frac = np.median(pool["P"]["xs"]["main"] / pool["P"]["depth"][:, None], axis=0)
    S.update(design=design(), pools=pool, frac=frac)
    n = 3 if a.check else a.n_sub
    res = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(one, [(i,) for i in range(n)])):
            res.extend(r)
            if (k + 1) % 25 == 0 or k + 1 == n:
                print(f"{k + 1}/{n} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res)
    res.to_csv(OUT / f"draws{sfx}.csv", index=False)
    rt = rates(res)
    rt.to_csv(OUT / f"rates{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(rt.round(3).T.to_string(), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
