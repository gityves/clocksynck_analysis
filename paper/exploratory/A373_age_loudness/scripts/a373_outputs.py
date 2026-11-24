"""A373 part 3 (registered): are the rhythmic output genes of old donors quieter once death type is balanced?

Output programme per tissue (blind to age): all donors of both cohorts with a CHIRAL donor internal phase (DIP), all ages,
genes with median CPM >= 5 outside the clock-like and responsive lists, cosinor of log(CPM + 0.5) on the DIP, BH q <= 0.05,
top 100 by F (time axis changed from the time of death to the DIP before any value, user correction);
waveform = their cosinor amplitudes and acrophases. H373.3: A373's loudness statistic (debiased plane variance W, old /
young on the A231 matched pairs, median over tissues of both cohorts) on the programme. D373.3t (both cohorts, on the
DIP, the CHIRAL analogue): per-gene cosinor amplitude ratio old / young and the rhythmic share per arm.

    python a373_outputs.py --stage power   # programmes, label null of Q_out, companions; no observed age contrast
    python a373_outputs.py --stage score
Outputs (`../output/`): out_programmes.csv, out_power.csv, out_null.csv, out_companions.csv | out_tests.csv,
out_tissues.csv. Controlled access: ages and times of death stay in memory; only aggregates are written.
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import a373_age_loudness as P  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE, RESPONSIVE_GENES  # noqa: E402
from clocksync._core.responsive_broad import RESPONSIVE_BROAD  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, pseudo_waveforms, states  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
from common import gtex  # noqa: E402

OUT = P.OUT
EXCLUDE = set(CLOCK_LIKE) | set(RESPONSIVE_GENES) | set(RESPONSIVE_BROAD)
N_PROG, MIN_PROG, N_SETS, N_PERM, N_BOOT = 100, 20, 200, 1000, 2000
W24 = 2 * np.pi / 24
META = ["donor", "compartment", "depth"]


def cosinor(logx, h):
    """F, p, amplitude, acrophase per column of logx (n x G) on [1, cos, sin] of hours h."""
    t = np.asarray(h, float) * W24
    X = np.c_[np.ones(len(t)), np.cos(t), np.sin(t)]
    B = np.linalg.lstsq(X, logx, rcond=None)[0]
    rss = ((logx - X @ B) ** 2).sum(0)
    tss = ((logx - logx.mean(0)) ** 2).sum(0)
    n = len(t)
    F = ((tss - rss) / 2) / (rss / (n - 3))
    return F, stats.f.sf(F, 2, n - 3), np.hypot(B[1], B[2]), np.arctan2(B[2], B[1])


class Tissue:
    def __init__(self, x, genes, donors, depth):
        self.x, self.genes, self.donors, self.depth = x, genes, donors, depth
        self.col = {g: k for k, g in enumerate(genes)}
        self.pos = {d: k for k, d in enumerate(donors)}

    def frame(self, gl):
        k = [self.col[g] for g in gl]
        return pd.DataFrame(self.x[:, k], columns=list(gl)).assign(depth=self.depth)

    def zv(self, w):
        seen, idx = set(), []
        for i, x in enumerate(w.genes):     # a gene drawn twice in a companion set is used once
            if x in self.col and x not in seen:
                seen.add(x)
                idx.append(i)
        g = [w.genes[i] for i in idx]
        ww = Waveform(tuple(g), np.asarray(w.amp)[idx], np.asarray(w.psi)[idx], None, None, w.name)
        f = self.frame(g)
        L = projector(f[g], ww)
        b, v = states(f, L, ww, self.depth)
        return b[:, 0] + 1j * b[:, 1], v


def load(c, tissues):
    h = P.ML.COH[c]
    pb, don = gtex.load(tissues=tissues, hardy=h, genes="all", depth="library_size", require_tod=(c == "fast"))
    cols = pd.Index([g for g in pb.columns if g not in META])
    genes = list(cols[~cols.duplicated(keep=False)])     # symbols that occur once only (no ambiguous duplicates)
    out = {}
    for t in tissues:
        r = pb[pb.compartment == t]
        if len(r) == 0:
            continue
        out[t] = Tissue(r.loc[:, genes].to_numpy(np.float32), genes, list(r.donor.astype(str)), r.depth.to_numpy(float))
    cpm_mean = pd.Series(np.concatenate([v.x for v in out.values()]).mean(0), index=genes)
    del pb
    return out, don, cpm_mean


def programme(Ts, dip):
    """Ts: the tissue's data in each cohort (same gene order); pooled over cohorts, ages and donors with a DIP."""
    genes = Ts[0].genes
    cpm = np.concatenate([T.x / T.depth[:, None] * 1e6 for T in Ts])
    h = np.concatenate([dip.reindex(T.donors).to_numpy(float) for T in Ts])
    ok = np.median(cpm, 0) >= 5
    ok &= np.array([g not in EXCLUDE for g in genes])
    m = np.isfinite(h)
    F, p, amp, psi = cosinor(np.log(cpm[m][:, ok] + 0.5), h[m])
    q = stats.false_discovery_control(p, method="bh")
    g = np.array(genes)[ok]
    sel = np.flatnonzero(q <= 0.05)
    if len(sel) < MIN_PROG:
        return None, len(sel)
    top = sel[np.argsort(-F[sel])[:N_PROG]]
    return Waveform(tuple(g[top]), amp[top], psi[top], None, None, "outputs"), len(sel)


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    cohorts = P.ML.load_cohorts()
    plans, _, _, _, M, Bal = P.ML.build(P.C.DESIGN, cohorts)
    sys.path.insert(0, str(P.EXP.parent / "results" / "figure_3" / "scripts" / "gtex_cues"))
    import compute_dip_arm as DA
    dip = DA.read_dip()
    data = {}
    for c in ("fast", "vent"):
        data[c], _, cm = load(c, list(plans[c]))
        if c == "fast":
            lm = np.log(cm[cm >= 1] / 1.0)
            lm = lm[~lm.index.duplicated()]
        P.say(c, "loaded", f"({time.time() - t0:.0f}s)")
    progs, prow = {}, []
    for t, T in data["fast"].items():
        Ts = [T] + ([data["vent"][t]] if t in data["vent"] and data["vent"][t].genes == T.genes else [])
        w, nq = programme(Ts, dip)
        prow.append(dict(tissue=t, rhythmic_q05=nq, programme=0 if w is None else len(w.genes)))
        if w is not None:
            progs[t] = w
    pd.DataFrame(prow).to_csv(OUT / "out_programmes.csv", index=False)
    P.say("programmes:", len(progs), "tissues")
    keys = [(c, t) for c in plans for t in plans[c] if t in progs and t in data[c]]
    pairs, pooled = {}, {}
    for c, t in keys:
        T, q = data[c][t], cohorts[c][1]
        pr = [(T.pos[str(s[0])], T.pos[str(s[1])]) for s in plans[c][t] if str(s[0]) in T.pos and str(s[1]) in T.pos]
        pairs[(c, t)] = np.array(pr)
        op, yp = P.C.pools(q)
        pooled[(c, t)] = (np.array([T.pos[d] for d in T.donors if d in q.index and op.get(d, False)]),
                          np.array([T.pos[d] for d in T.donors if d in q.index and yp.get(d, False)]))
    arms = {k: (p[:, 0], p[:, 1]) for k, p in pairs.items()}
    sets = {t: {kind: pseudo_waveforms(lm, N_SETS, np.random.default_rng([373, tag, i]), progs[t], kind=kind)
                for kind, tag in (("matched", 5), ("responsive", 6))} for i, t in enumerate(sorted(progs))}
    clock = {k: data[k[0]][k[1]].zv(progs[k[1]]) for k in keys}
    P.say("states of the programmes", f"({time.time() - t0:.0f}s)")

    if a.stage == "power":
        rng = np.random.default_rng([373, 7])
        null = []
        for _ in range(N_PERM):
            sw = {}
            for k, p in pairs.items():
                f = rng.random(len(p)) < 0.5
                sw[k] = (np.where(f, p[:, 1], p[:, 0]), np.where(f, p[:, 0], p[:, 1]))
            null.append(P.Q_of(clock, sw)[0])
        null = np.array(null)
        pd.DataFrame(dict(Q=null)).to_csv(OUT / "out_null.csv", index=False)
        rows = []
        for kind in ("matched", "responsive"):
            for j in range(N_SETS):
                st = {k: data[k[0]][k[1]].zv(sets[k[1]][kind][j]) for k in keys}
                rows.append(dict(kind=kind, set=j, Q=P.Q_of(st, arms)[0], Q_pooled=P.Q_of(st, pooled)[0]))
            P.say(kind, "companions done", f"({time.time() - t0:.0f}s)")
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "out_companions.csv", index=False)
        sd = float(np.std(null))
        res = dict(tissues=len(keys), tissues_with_programme=len(progs), label_null_sd=sd,
                   label_null_q05=float(np.quantile(null, .05)), expected_half_width=1.645 * sd, margin=P.MARGIN,
                   matched_median=float(comp[comp.kind == "matched"].Q.median()),
                   matched_q05=float(comp[comp.kind == "matched"].Q.quantile(.05)),
                   responsive_median=float(comp[comp.kind == "responsive"].Q.median()),
                   responsive_q05=float(comp[comp.kind == "responsive"].Q.quantile(.05)))
        res["gate_H373.3_equivalence"] = bool(res["expected_half_width"] < P.MARGIN)
        pd.Series(res).to_csv(OUT / "out_power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "out_power.csv", index_col=0)["value"]
    comp = pd.read_csv(OUT / "out_companions.csv")
    null = pd.read_csv(OUT / "out_null.csv").Q.to_numpy()
    q, L = P.Q_of(clock, arms)
    qp, Lp = P.Q_of(clock, pooled)
    rng = np.random.default_rng([373, 8])
    bt = []
    for _ in range(N_BOOT):
        ab = {}
        for k, p in pairs.items():
            i = rng.integers(0, len(p), len(p))
            ab[k] = (p[i, 0], p[i, 1])
        bt.append(P.Q_of(clock, ab)[0])
    lo, hi = np.quantile(bt, [.05, .95])
    mt, rs = comp[comp.kind == "matched"].Q, comp[comp.kind == "responsive"].Q
    r = dict(test="H373.3", Q_out=q, amp_ratio=float(np.exp(q / 2)), lo90=float(lo), hi90=float(hi),
             amp_ratio_lo90=float(np.exp(lo / 2)), amp_ratio_hi90=float(np.exp(hi / 2)), Q_star=q - float(mt.median()),
             p_label=float((1 + np.sum(null <= q)) / (1 + len(null))), pct_matched=float((mt < q).mean() * 100),
             pct_responsive=float((rs < q).mean() * 100), Q_pooled=qp, amp_ratio_pooled=float(np.exp(qp / 2)),
             Q_fast=float(np.nanmedian([v for k, v in L.items() if k[0] == "fast"])),
             Q_vent=float(np.nanmedian([v for k, v in L.items() if k[0] == "vent"])))
    r["effect"] = bool(q < float(mt.quantile(.05)) and q < float(np.quantile(null, .05)) and hi < 0)
    r["equivalence_decidable"] = str(pw["gate_H373.3_equivalence"]) == "True"
    r["equivalence"] = bool(lo > -P.MARGIN and hi < P.MARGIN) if r["equivalence_decidable"] else None
    # D373.3t: timed, fast deaths, matched pairs and pooled arms
    trow = []
    for (c, t), p in pairs.items():
        T, w = data[c][t], progs[t]
        cpm = T.frame(list(w.genes))[list(w.genes)].to_numpy(float) / T.depth[:, None] * 1e6
        h = dip.reindex(T.donors).to_numpy(float)
        res_t = {}
        for lab, (o, y) in (("matched", (p[:, 0], p[:, 1])), ("pooled", pooled[(c, t)])):
            o, y = o[np.isfinite(h[o])], y[np.isfinite(h[y])]
            if len(o) < 10 or len(y) < 10:
                continue
            _, po, ao, _ = cosinor(np.log(cpm[o] + 0.5), h[o])
            _, py, ay, _ = cosinor(np.log(cpm[y] + 0.5), h[y])
            res_t[f"log_amp_ratio_{lab}"] = float(np.median(np.log(ao / ay)))
            res_t[f"rhythmic_share_diff_{lab}"] = float(np.mean(stats.false_discovery_control(po, method="bh") <= 0.05)
                                                         - np.mean(stats.false_discovery_control(py, method="bh") <= 0.05))
        trow.append(dict(cohort=c, tissue=t, pairs=len(p), L_matched=L[(c, t)], **res_t))
    tt = pd.DataFrame(trow)
    for col in [x for x in tt.columns if x.startswith(("log_amp", "rhythmic"))]:
        r[f"D373.3t_{col}_median"] = float(tt[col].median())
    r["D373.3t_amp_ratio_matched"] = float(np.exp(r["D373.3t_log_amp_ratio_matched_median"]))
    r["D373.3t_amp_ratio_pooled"] = float(np.exp(r["D373.3t_log_amp_ratio_pooled_median"]))
    pd.DataFrame([r]).to_csv(OUT / "out_tests.csv", index=False)
    tt.to_csv(OUT / "out_tissues.csv", index=False)
    print(pd.Series(r).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
