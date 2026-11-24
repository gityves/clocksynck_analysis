"""Post hoc (2026-11-24; descriptive, underpowered; plan written in registered_experiments.md before computing): is the
person-level rhythm amplitude shared by two skin samples carried also by random genes?

Part 1: our clock statistic (A364 two sites, A365 two layers: r of log first-harmonic amplitude, core loop and 14
genes) against the companions of their power stages.
Part 2: Saint-Antoine & Anafi's global amplitude (re-implemented from their `cosinor_pipeline.R` / `fig3.R`): per layer,
collective cosinor with subject intercepts (BH q <= 0.05), per-subject amplitude / mesor on the linear scale,
percentile across subjects, mean over the layer's rhythmic genes; r between layers. Control: 200 random non-rhythmic
gene sets per layer matched to the rhythmic set's mean-expression deciles (seed [372, 1]).

    python posthoc_skin_shared.py
Outputs (`../output/`): posthoc_clock.csv, posthoc_global.csv, posthoc_global_random.csv.
"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import a364_skin_volume as S4  # noqa: E402
import a365_delolmo_volume as S5  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
N_RAND = 200


def part1():
    rows = []
    d, x = S5.load()
    comp = pd.read_csv(OUT / "a365_companions.csv")
    for gs, w in (("core", S5.CORE), ("14 genes", S5.WF)):
        rows.append(("del Olmo layers (A365)", gs, S5.r_of(S5.amplitudes(S5.site_states(x, d, w), d)[0]), comp[comp.gene_set == gs]))
    d4, x4 = S4.load()
    comp4 = pd.read_csv(OUT / "companions.csv")
    for gs, w in (("core", S4.CORE), ("14 genes", S4.WF)):
        rows.append(("GSE326530 sites (A364)", gs, S4.r_of(S4.amplitudes(S4.site_states(x4, d4, w), d4)[0]), comp4[comp4.gene_set == gs]))
    out = []
    for ds, gs, r, c in rows:
        mt, rs = c[c.kind == "matched"].r, c[c.kind == "responsive"].r
        out.append(dict(dataset=ds, gene_set=gs, r=r, matched_median=mt.median(), matched_q95=mt.quantile(.95),
                        pct_matched=(mt < r).mean() * 100, responsive_median=rs.median(),
                        responsive_q95=rs.quantile(.95), pct_responsive=(rs < r).mean() * 100))
    return pd.DataFrame(out)


class Global:
    """Saint-Antoine & Anafi global amplitude for one layer: Y samples x genes (linear), subjects, times (h)."""

    def __init__(self, Y, subj, t):
        keep = Y.mean(0) >= Y.mean(0).median()                      # expressed: top half by mean (our choice)
        self.Y, self.genes = Y.loc[:, keep].to_numpy(float), Y.columns[keep]
        self.subj, self.t = np.asarray(subj), np.asarray(t, float)
        self.subjects = sorted(set(self.subj))
        D = pd.get_dummies(pd.Series(self.subj)).to_numpy(float)
        w = 2 * np.pi * self.t / 24
        Xf, Xr = np.column_stack([D, np.cos(w), np.sin(w)]), D
        rss = lambda X: ((self.Y - X @ np.linalg.lstsq(X, self.Y, rcond=None)[0]) ** 2).sum(0)  # noqa: E731
        rf, rr = rss(Xf), rss(Xr)
        dfe = len(self.Y) - Xf.shape[1]
        F = ((rr - rf) / 2) / (rf / dfe)
        p = stats.f.sf(F, 2, dfe)
        self.q = pd.Series(stats.false_discovery_control(p, method="bh"), index=self.genes)
        rel = np.empty((len(self.subjects), self.Y.shape[1]))
        for k, s in enumerate(self.subjects):
            i = self.subj == s
            X = np.column_stack([np.ones(i.sum()), np.cos(w[i]), np.sin(w[i])])
            b = np.linalg.lstsq(X, self.Y[i], rcond=None)[0]
            rel[k] = np.hypot(b[1], b[2]) / b[0]
        n = len(self.subjects)
        self.pct = pd.DataFrame((stats.rankdata(rel, axis=0) - 0.5) / n, index=self.subjects, columns=self.genes)
        self.logmean = pd.Series(np.log(self.Y.mean(0)), index=self.genes)

    def rhythmic(self):
        return list(self.q.index[self.q <= 0.05])

    def global_amp(self, genes):
        return self.pct[genes].mean(1)

    def matched_random(self, genes, rng):
        pool = self.q.index[self.q > 0.5]
        bins = pd.qcut(self.logmean, 10, labels=False)
        out = []
        for b, k in bins[genes].value_counts().items():
            cand = pool[bins[pool] == b]
            out += list(rng.choice(cand, k, replace=len(cand) < k))
        return out


def part2():
    rng = np.random.default_rng([372, 1])
    rows, rrows = [], []
    d, x = S5.load()
    d4, x4 = S4.load()
    cpm = x4.div(x4.sum(1), axis=0) * 1e6
    sets = {"del Olmo layers": (np.exp(x), d.subject.to_numpy(), d.t.to_numpy(), d.layer.to_numpy(), S5.LAYERS),
            "GSE326530 sites": (cpm, d4.vol.astype(str).to_numpy(), d4.t.to_numpy(), d4.site.to_numpy(), ("P", "E"))}
    for ds, (Y, subj, t, lay, labs) in sets.items():
        G = {L: Global(Y[lay == L].reset_index(drop=True), subj[lay == L], t[lay == L]) for L in labs}
        R = {L: G[L].rhythmic() for L in labs}
        ga = [G[L].global_amp(R[L]) for L in labs]
        common = ga[0].index.intersection(ga[1].index)
        r = float(np.corrcoef(ga[0][common], ga[1][common])[0, 1])
        rand = []
        for j in range(N_RAND):
            gr = [G[L].global_amp(G[L].matched_random(R[L], rng)) for L in labs]
            rand.append(float(np.corrcoef(gr[0][common], gr[1][common])[0, 1]))
            rrows.append(dict(dataset=ds, set=j, r=rand[-1]))
        rand = np.array(rand)
        rows.append(dict(dataset=ds, subjects=len(common), **{f"rhythmic_{L}": len(R[L]) for L in labs},
                         **{f"expressed_{L}": len(G[L].genes) for L in labs}, r=r, R2=r ** 2,
                         random_median=float(np.median(rand)), random_q95=float(np.quantile(rand, .95)),
                         pct_random=float((rand < r).mean() * 100)))
        print(rows[-1], flush=True)
    return pd.DataFrame(rows), pd.DataFrame(rrows)


def main():
    limit_threads(4)
    t0 = time.time()
    p1 = part1()
    p1.to_csv(OUT / "posthoc_clock.csv", index=False)
    print(p1.round(3).to_string(index=False), flush=True)
    p2, rr = part2()
    p2.to_csv(OUT / "posthoc_global.csv", index=False)
    rr.to_csv(OUT / "posthoc_global_random.csv", index=False)
    print(p2.round(3).T.to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
