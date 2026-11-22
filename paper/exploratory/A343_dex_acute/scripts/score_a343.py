"""A343 scoring (registered in `paper/docs/registered_experiments.md`, A343 and its identity / power note): does an acute
dexamethasone input raise pre-mRNA before mature mRNA in living human fat and liver? One run.

Units: per tissue, patients with a non-surgical (N) and a surgical (S) biopsy; primary = N and S share strandedness.
Per run and gene, r = log((u + 0.5) / (s + 0.5)) (u = `counts_intron`, subtraction; `counts_intron_direct` as the
sensitivity), s = `counts_exon`. GR set: FKBP5, TSC22D3, ZBTB16, KLF9. Controls: 50 genes per set gene drawn from its
bin (deciles of N-run mean exon CPM x deciles of intron-only length, over genes outside the clock panel, `CLOCK_LIKE`,
the responsive families and the GR set, CPM >= 10, intron-only length >= 1 kb; a set gene outside the pool's range takes
the edge bin). R = mean r(set) - mean r(controls) per run, D = R(S) - R(N) per patient. E = the same with log CPM in
place of r. Companions: 200 four-gene sets, each gene drawn from the bin of the GR gene it replaces, with their own
controls. Sign-flip tests over patients (10,000).
- P343.0: DEX patients, mean E > 0 (one-sided) and above the companions' q95.
- H343.1: DEX patients, mean D > 0 (one-sided) and above the companions' q95. Scored in fat (decidable); liver
  descriptive (power note).
- H343.2 (descriptive by the power note): DEX minus NO_DEX mean D, permutation of arm labels (10,000), primary units;
  all units with the library-type pair as covariate.
- H343.3 (descriptive by the power note): DEX minus NO_DEX in the PER1 S - N log CPM change and in the S - N
  displacement of the core-loop state (GTEx tissue waveform, one projector per tissue), 90 % bootstrap intervals.
Descriptive: per-gene D and E (GR set, PER1); the NO_DEX S - N change; two-tissue patients (code 53 excluded by the
identity check): liver-fat complex correlation and offset of the 14-gene state at N and at S.

    python score_a343.py [--check]
`--check`: the GR set replaced by four random control-pool genes; PER1, the core loop and the two-tissue clock pairing
are skipped (no GR or clock count read). Outputs (`../output/score[_check]/`): verdicts.csv, genes.csv, h343_2.csv, h343_3.csv, two_tissue.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import power_a343 as PW  # noqa: E402
from clocksync import waveform as W  # noqa: E402
from clocksync.paired import ccorr, projector, states  # noqa: E402

OUT = HERE.parents[1] / "output" / "score"
SEED, N_COMP, N_CTRL, N_FLIP, N_PERM, N_BOOT = 3431, 200, 50, 10000, 10000, 2000
HOURS = 24 / (2 * np.pi)
EXCLUDE_TWO_TISSUE = {"53"}
SCORED = {"A": "H343.1 scored", "L": "descriptive (power note)"}
CHECK = False


def edges(v, k=10):
    return np.quantile(v, np.linspace(0, 1, k + 1)[1:-1])


class Matcher:
    """Expression x intron-length bins over the control pool of one tissue."""

    def __init__(self, ex, ann, runs_n):
        cpm = ex[runs_n].div(ex[runs_n].sum(), axis=1) * 1e6
        self.lc = np.log(cpm.mean(axis=1) + 1e-3)
        self.ll = np.log(ann.intron_only_len_bp_unstranded.reindex(ex.index).astype(float).clip(lower=1))
        pool = [g for g in ex.index if g not in PW.TARGETS and np.exp(self.lc[g]) >= 10 and np.exp(self.ll.get(g, 0)) >= 1000]
        self.ec, self.el = edges(self.lc[pool]), edges(self.ll[pool])
        key = pd.Series([self.key(g) for g in pool], index=pool)
        self.by = {k: v.index.to_numpy() for k, v in key.groupby(key)}

    def key(self, g):
        return (int(np.searchsorted(self.ec, self.lc[g])), int(np.searchsorted(self.el, self.ll[g])))

    def draw(self, g, n, rng):
        return rng.choice(self.by[self.key(g)], n)


def patient_stat(X, genes, ctrl, pts):
    """Per patient S - N of (mean X over genes - mean X over controls); X genes x runs."""
    def R(runs):
        return X.loc[genes, runs].mean().to_numpy() - X.loc[ctrl, runs].mean().to_numpy()
    return R(pts.run_S.tolist()) - R(pts.run_N.tolist())


def signflip(d, rng, n=N_FLIP):
    null = rng.choice([-1, 1], (n, len(d))) @ d / len(d)
    return float((1 + np.sum(null >= d.mean())) / (1 + n))


def tissue_scores(tis, ex, it, itd, ann, t, rng) -> tuple[list, list]:
    pts_all = PW.primary(t, tis)
    m = Matcher(ex, ann, pts_all.run_N.tolist())
    runs = pts_all.run_N.tolist() + pts_all.run_S.tolist()
    cpm = np.log((ex[runs] + 0.5) / ex[runs].sum() * 1e6)
    Rr = {"intron": np.log((it[runs] + 0.5) / (ex[runs] + 0.5)), "direct": np.log((itd[runs] + 0.5) / (ex[runs] + 0.5))}
    ctrl = np.concatenate([m.draw(g, N_CTRL, rng) for g in PW.GR])
    comps = []
    for _ in range(N_COMP):
        genes = [str(m.draw(g, 1, rng)[0]) for g in PW.GR]
        comps.append((genes, np.concatenate([m.draw(g, N_CTRL, rng) for g in genes])))
    rows, genes_rows = [], []
    for units in ("primary", "all"):
        for arm in ("DEX", "NO_DEX"):
            p = pts_all[(pts_all.dex == arm) & (pts_all.same | (units == "all"))]
            for readout, X in (("E", cpm), ("D", Rr["intron"]), ("D_direct", Rr["direct"])):
                d = patient_stat(X, list(PW.GR), ctrl, p)
                cv = np.array([patient_stat(X, g, c, p).mean() for g, c in comps])
                rows.append(dict(tissue=tis, units=units, arm=arm, readout=readout, patients=len(p), mean=d.mean(),
                                 p_signflip=signflip(d, rng), comp_q95=float(np.quantile(cv, .95)),
                                 comp_median=float(np.median(cv)), above_comp=bool(d.mean() > np.quantile(cv, .95))))
            if units == "primary":
                for g in [*PW.GR, *([] if CHECK else ["PER1"])]:
                    c = m.draw(g, N_CTRL, rng)
                    genes_rows.append(dict(tissue=tis, arm=arm, gene=g, patients=len(p),
                                           E=patient_stat(cpm, [g], c, p).mean(),
                                           D=patient_stat(Rr["intron"], [g], c, p).mean(),
                                           D_direct=patient_stat(Rr["direct"], [g], c, p).mean()))
    return rows, genes_rows, (pts_all, cpm, Rr["intron"], ctrl)


def h2(tis, pts, X, ctrl, rng) -> list:
    out = []
    for units in ("primary", "all"):
        p = pts[pts.same | (units == "all")].reset_index(drop=True)
        d = patient_stat(X, list(PW.GR), ctrl, p)
        arm = (p.dex == "DEX").to_numpy()
        lib = (p.strand_N + "|" + p.strand_S).to_numpy()
        Z = pd.get_dummies(lib, drop_first=True).to_numpy(float) if units == "all" else np.zeros((len(p), 0))

        def diff(a):
            M = np.column_stack([np.ones(len(p)), a.astype(float), Z])
            return np.linalg.lstsq(M, d, rcond=None)[0][1]
        obs = diff(arm)
        null = np.array([diff(rng.permutation(arm)) for _ in range(N_PERM)])
        out.append(dict(tissue=tis, units=units, dex=int(arm.sum()), no_dex=int((~arm).sum()), diff=obs,
                        p_perm=float((1 + np.sum(null >= obs)) / (1 + N_PERM))))
    return out


def h3(tis, ex, pts, cpm, rng) -> list:
    wf = W.gtex(PW.TISSUE[tis]).drop(PW.CORE_DROP)
    g = list(wf.genes)
    runs = pts.run_N.tolist() + pts.run_S.tolist()
    x = ex.loc[g, runs].T
    b, _ = states(x, projector(x, wf), wf, depth=ex[runs].sum().to_numpy(float))
    n = len(pts)
    z = (b[n:, 0] - b[:n, 0]) + 1j * (b[n:, 1] - b[:n, 1])
    per1 = (cpm.loc["PER1", pts.run_S.tolist()].to_numpy() - cpm.loc["PER1", pts.run_N.tolist()].to_numpy())
    out = []
    for units in ("primary", "all"):
        k = (pts.same | (units == "all")).to_numpy()
        arm = (pts.dex == "DEX").to_numpy()
        dx, nx = k & arm, k & ~arm
        if dx.sum() < 3 or nx.sum() < 3:
            continue
        iz, jz = np.flatnonzero(dx), np.flatnonzero(nx)
        bs = np.array([(per1[rng.choice(iz, len(iz))].mean() - per1[rng.choice(jz, len(jz))].mean(),
                        abs(z[rng.choice(iz, len(iz))].mean() - z[rng.choice(jz, len(jz))].mean()))
                       for _ in range(N_BOOT)])
        out.append(dict(tissue=tis, units=units, dex=int(dx.sum()), no_dex=int(nx.sum()),
                        per1_diff=per1[dx].mean() - per1[nx].mean(), per1_lo=np.quantile(bs[:, 0], .05),
                        per1_hi=np.quantile(bs[:, 0], .95), core_disp_diff=abs(z[dx].mean() - z[nx].mean()),
                        core_ub90=np.quantile(bs[:, 1], .95), core_disp_dex=abs(z[dx].mean()),
                        core_disp_nodex=abs(z[nx].mean())))
    return out


def two_tissue(ex, t) -> list:
    out = []
    pl, pa = PW.primary(t, "L"), PW.primary(t, "A")
    ids = sorted(set(pl.pid) & set(pa.pid) - EXCLUDE_TWO_TISSUE)
    pl, pa = pl.set_index("pid").loc[ids], pa.set_index("pid").loc[ids]
    for bio in ("N", "S"):
        st = {}
        for tis, p in (("L", pl), ("A", pa)):
            wf = W.gtex(PW.TISSUE[tis])
            runs = p[f"run_{bio}"].tolist()
            x = ex.loc[list(wf.genes), runs].T
            b, v = states(x, projector(x, wf), wf, depth=ex[runs].sum().to_numpy(float))
            st[tis] = (b - b.mean(0), v)
        c = ccorr(st["L"][0], st["L"][1], st["A"][0], st["A"][1])
        out.append(dict(biopsy=bio, patients=len(ids), coherence=abs(c), offset_liver_minus_fat_h=np.angle(c) * HOURS))
    return out


def main():
    global CHECK, OUT
    CHECK = "--check" in sys.argv
    if CHECK:
        OUT = OUT.with_name("score_check")
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t = PW.samples()
    ex, it, itd = PW.load("exon"), PW.load("intron"), PW.load("intron_direct")
    ann = pd.read_csv(PW.DATA / "gene_annotation.tsv.gz", sep="\t").drop_duplicates("symbol").set_index("symbol")
    if CHECK:
        pts = PW.primary(t, "A")
        m = Matcher(ex, ann, pts.run_N.tolist())
        PW.GR = tuple(rng.choice(np.concatenate(list(m.by.values())), 4, replace=False))
        print("check: GR set replaced by", PW.GR, flush=True)
    allrows, generows, h2rows, h3rows = [], [], [], []
    for tis in ("A", "L"):
        rows, gr, (pts, cpm, X, ctrl) = tissue_scores(tis, ex, it, itd, ann, t, rng)
        allrows += rows
        generows += gr
        h2rows += h2(tis, pts, X, ctrl, rng)
        if not CHECK:
            h3rows += h3(tis, ex, pts, cpm, rng)
    v = pd.DataFrame(allrows)
    v["status"] = v.tissue.map(SCORED)
    prim = (v.units == "primary") & (v.arm == "DEX")
    v["P343_0"] = np.where(prim & (v.readout == "E"), (v.p_signflip < .05) & v.above_comp, None)
    v["H343_1"] = np.where(prim & (v.readout == "D"), (v.p_signflip < .05) & v.above_comp, None)
    v.to_csv(OUT / "verdicts.csv", index=False)
    pd.DataFrame(generows).to_csv(OUT / "genes.csv", index=False)
    pd.DataFrame(h2rows).to_csv(OUT / "h343_2.csv", index=False)
    pd.DataFrame(h3rows).to_csv(OUT / "h343_3.csv", index=False)
    tt = pd.DataFrame([] if CHECK else two_tissue(ex, t))
    tt.to_csv(OUT / "two_tissue.csv", index=False)
    pd.set_option("display.width", 220)
    for f in ("verdicts", "genes", "h343_2") + (() if CHECK else ("h343_3", "two_tissue")):
        print(f"== {f}", pd.read_csv(OUT / f"{f}.csv").round(3).to_string(index=False), sep="\n")


if __name__ == "__main__":
    main()
