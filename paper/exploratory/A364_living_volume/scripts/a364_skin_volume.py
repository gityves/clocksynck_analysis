"""A364 (registered): do living persons differ in a clock amplitude shared by two skin sites, measured from a timed day?

GSE326530 (A342 delivery): 20 volunteers x photoprotected (P) / photoexposed (E) skin x 4 clock times (0, 6, 12, 18 h).
One projector on all 160 samples (GTEx not-sun-exposed skin waveform, as A342); states centred within volunteer x site
(removes each person's offsets); per volunteer and site, the first-harmonic amplitude a = |sum_t z_t exp(-i 2 pi t / 24)|
(independent of the person's phase and of the biopsy time). Statistic r = Pearson correlation over volunteers of
log a_P and log a_E (core loop primary, 14 genes reported). Null: volunteers permuted (10,000). Companions: 200
expression-matched and 200 responsive pseudo-clock waveforms (seeds [364, 1] / [364, 2]). Core readout: the 364
three-gene drops of the 14 genes (CLAUDE.md gene-drop rule).

    python a364_skin_volume.py --stage power   # residual noise per site, simulation, companions; no clock r
    python a364_skin_volume.py --stage score   # registered test (after the power note)
Outputs (`../output/`): power.csv, companions.csv, sim.csv | tests.csv, drops.csv.
"""
import argparse
import itertools
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "A342_skin_timed" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import power_a342 as A  # noqa: E402  (sets the import paths)
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import projector, pseudo_waveforms, states  # noqa: E402

OUT = HERE.parents[1] / "output"
WF = A.WF
TRIO = ["PER1", "PER2", "NR1D1"]
CORE = WF.drop(TRIO)
CORE_IDX = [k for k, g in enumerate(WF.genes) if g in CORE.genes]
N_SETS, N_PERM, N_SIM = 200, 10_000, 2000
SIGMA_A = 0.352   # GTEx person spread of log amplitude (A331 power stage)
T = np.array(A.TIMES, float)


def load():
    d = A.design().sort_values(["vol", "site", "t"]).reset_index(drop=True)
    x = pd.read_csv(A.DATA / "counts_exon.tsv.gz", sep="\t")
    x["symbol"] = x.symbol.astype(str).map(canonical)
    x = x.drop_duplicates("symbol").set_index("symbol").drop(columns="gene_id")[d.run].T
    x.index = range(len(d))
    return d, x


def site_states(x, d, w):
    """Per sample complex state on one projector over all samples, centred within volunteer x site."""
    g = list(w.genes)
    X = x[g]
    L = projector(X, w)
    b, _ = states(X, L, w, d.depth_s.to_numpy(float))
    z = b[:, 0] + 1j * b[:, 1]
    key = d.vol.astype(str) + d.site
    return z - pd.Series(z).groupby(key.to_numpy()).transform("mean").to_numpy()


def amplitudes(z, d):
    """a[v, s] = |sum_t z exp(-i 2 pi t / 24)| (4 times; 4 m for a clock radius m); per-sample noise variance per site
    from the residual after the first harmonic."""
    out, res = {}, {"P": [], "E": []}
    for (v, s), i in d.groupby(["vol", "site"]).groups.items():
        i = np.asarray(i)
        e = np.exp(-1j * 2 * np.pi * d.t.to_numpy(float)[i] / 24)
        beta = np.sum(z[i] * e)
        out[(v, s)] = abs(beta)
        res[s].append(np.abs(z[i] - beta / len(i) * np.conj(e)) ** 2)
    a = pd.Series(out).unstack()
    # centred 4-point series: the residual after the first harmonic keeps 2 of 4 complex degrees of freedom
    return a[["P", "E"]], {s: 2 * float(np.mean(np.concatenate(r))) for s, r in res.items()}


def r_of(a):
    la = np.log(a.to_numpy(float))
    return float(np.corrcoef(la[:, 0], la[:, 1])[0, 1])


def core_of(w, j):
    from clocksync.waveform import Waveform
    return Waveform(tuple(w.genes[k] for k in CORE_IDX), CORE.amp, CORE.psi, None, None, f"core{j}")


def simulate(m, s2, n, sigma, rng):
    """r under a person amplitude of log SD sigma shared by the sites; site mean amplitude m[s] (per-sample clock-plane
    radius) and residual noise s2[s] (per-sample complex variance), 4 times, uniform phases."""
    A_ = np.exp(rng.normal(0, sigma, n))
    ph = rng.uniform(0, 2 * np.pi, n)
    la = []
    e = np.exp(-1j * 2 * np.pi * T / 24)
    for s in ("P", "E"):
        z = m[s] * A_[:, None] * np.exp(1j * (2 * np.pi * T[None, :] / 24 + ph[:, None]))
        z = z + np.sqrt(s2[s] / 2) * (rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape))
        z = z - z.mean(1, keepdims=True)
        la.append(np.log(np.abs((z * e[None, :]).sum(1))))
    return float(np.corrcoef(la[0], la[1])[0, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    d, x = load()
    n = d.vol.nunique()
    genome = x.loc[:, x.sum(0) > 0]
    lm = np.log(genome.div(genome.sum(1), axis=0).mean(0) + 1e-12)
    lm = lm[genome.mean(0) >= 5]
    lm = lm[~lm.index.duplicated()]
    print(f"{len(d)} samples, {n} volunteers ({time.time() - t0:.0f}s)", flush=True)

    if a.stage == "power":
        rows = []
        for kind, tag in (("matched", 1), ("responsive", 2)):
            for j, w in enumerate(pseudo_waveforms(lm, N_SETS, np.random.default_rng([364, tag]), WF, kind=kind)):
                for gs, ww in (("14 genes", w), ("core", core_of(w, j))):
                    am, _ = amplitudes(site_states(x, d, ww), d)
                    rows.append(dict(kind=kind, set=j, gene_set=gs, r=r_of(am)))
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "companions.csv", index=False)
        rng = np.random.default_rng([364, 3])
        res, sims = {}, {}
        for gs, w in (("core", CORE), ("14 genes", WF)):
            am, s2 = amplitudes(site_states(x, d, w), d)
            # per-sample clock radius m: E|beta|^2 = 16 m^2 + 4 noise
            m = {s: float(np.sqrt(max(np.mean(am[s].to_numpy() ** 2) - 4 * s2[s], 1e-12)) / 4) for s in ("P", "E")}
            null = np.array([simulate(m, s2, n, 0.0, rng) for _ in range(N_SIM)])
            for lab, sg in (("gtex", SIGMA_A), ("half", SIGMA_A / 2)):
                sims[f"{gs}|{lab}"] = np.array([simulate(m, s2, n, sg, rng) for _ in range(N_SIM)])
            c = comp[comp.gene_set == gs]
            mq95 = float(c[c.kind == "matched"].r.quantile(0.95))
            thr = max(mq95, float(np.quantile(null, 0.95)))
            res.update({f"{gs}|noise_P": s2["P"], f"{gs}|noise_E": s2["E"], f"{gs}|m_P": m["P"], f"{gs}|m_E": m["E"],
                        f"{gs}|null_sim_q95": float(np.quantile(null, 0.95)), f"{gs}|matched_median": float(c[c.kind == "matched"].r.median()),
                        f"{gs}|matched_q95": mq95, f"{gs}|responsive_median": float(c[c.kind == "responsive"].r.median()),
                        f"{gs}|responsive_q95": float(c[c.kind == "responsive"].r.quantile(0.95)),
                        f"{gs}|expected_r_gtex": float(np.median(sims[f'{gs}|gtex'])),
                        f"{gs}|power_gtex": float(np.mean(sims[f"{gs}|gtex"] > thr)),
                        f"{gs}|power_half": float(np.mean(sims[f"{gs}|half"] > thr))})
        res["gate"] = bool(res["core|power_gtex"] >= 0.8)
        pd.DataFrame(sims).to_csv(OUT / "sim.csv", index=False)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("A364 not powered: not decidable at this design; not scored")
    comp = pd.read_csv(OUT / "companions.csv")
    rng = np.random.default_rng([364, 4])
    rows = []
    for gs, w in (("core", CORE), ("14 genes", WF)):
        am, _ = amplitudes(site_states(x, d, w), d)
        r = r_of(am)
        la = np.log(am.to_numpy(float))
        null = np.array([np.corrcoef(la[:, 0], la[rng.permutation(n), 1])[0, 1] for _ in range(N_PERM)])
        c = comp[comp.gene_set == gs]
        mt, rs = c[c.kind == "matched"].r, c[c.kind == "responsive"].r
        row = dict(test="H364.1" if gs == "core" else "H364.1 (14 genes, reported)", volunteers=n, r=r,
                   p_perm=float((1 + np.sum(null >= r)) / (1 + N_PERM)), matched_median=float(mt.median()),
                   matched_q95=float(mt.quantile(.95)), pct_matched=float((mt < r).mean() * 100),
                   responsive_median=float(rs.median()), responsive_q95=float(rs.quantile(.95)),
                   pct_responsive=float((rs < r).mean() * 100), expected_r_gtex=float(pw[f"{gs}|expected_r_gtex"]))
        row["passes"] = bool(r > row["matched_q95"] and row["p_perm"] < 0.05)
        rows.append(row)
    dr = []
    for c3 in itertools.combinations(WF.genes, 3):
        am, _ = amplitudes(site_states(x, d, WF.drop(list(c3))), d)
        dr.append(dict(drop="+".join(c3), trio=list(c3) == TRIO, r=r_of(am)))
    dr = pd.DataFrame(dr)
    dr.to_csv(OUT / "drops.csv", index=False)
    tr = float(dr[dr.trio].r.iloc[0])
    rows[0].update(drops_median=float(dr.r.median()), drops_q10=float(dr.r.quantile(.10)),
                   trio_drop_pct=float((dr.r < tr).mean() * 100))
    rows[0]["reading"] = ("living persons share a clock amplitude across two skin sites, beyond matched genes"
                          if rows[0]["passes"] else "not shown in living timed skin at this power")
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "tests.csv", index=False)
    print(t.round(4).T.to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
