"""A373 (registered): how large an age fall of absolute clock loudness is excluded once death type is balanced?

GTEx, A231 design (old 60-70 vs young 21-49 y; within cohort, 1:1 pairs exact on sex and Hardy class, nearest
ischaemic time). Per cohort and tissue: clock-plane states on counts (projector per tissue), per arm the debiased
clock-plane variance W (A315); L_t = log(W_old / W_young); Q = median over tissues of both cohorts. Companions:
expression-matched and responsive pseudo-clock sets. Matched-set bootstrap; label swaps within pairs.

    python a373_age_loudness.py --stage power   # matching, balance, label null of Q, companions; no clock Q
    python a373_age_loudness.py --stage score   # H373.1, H373.2, descriptive
Outputs (`../output/`): power.csv, companions.csv, null.csv, balance.csv | tests.csv, tissues.csv.
Controlled access: ages, Hardy classes and times of death stay in memory; only aggregates are written.
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A231_age_coupling" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import compute_a231 as C  # noqa: E402
import matched_label as ML  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, pseudo_waveforms, states  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
WF = ML.cs.gtex()
TRIO = ["PER1", "PER2", "NR1D1"]
CORE = WF.drop(TRIO)
N_SETS, N_BOOT, N_PERM, MIN_N = 200, 2000, 1000, 10
MARGIN = 0.21


def say(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


class Cohort:
    """Per tissue: donor -> row, and for each waveform the states z and Poisson variance trace v."""

    def __init__(self, pb):
        self.pb = pb
        self.tissues = sorted(pb.compartment.unique())
        self.rows = {t: pb.index[pb.compartment == t] for t in self.tissues}
        self.pos = {t: {d: k for k, d in enumerate(pb.donor[self.rows[t]].astype(str))} for t in self.tissues}

    def zv(self, w, t):
        x = self.pb.loc[self.rows[t]]
        k = [i for i, g in enumerate(w.genes) if g in x.columns]
        ww = w if len(k) == len(w.genes) else Waveform(tuple(w.genes[i] for i in k), np.asarray(w.amp)[k],
                                                       np.asarray(w.psi)[k], None, None, w.name)
        L = projector(x[list(ww.genes)], ww)
        b, v = states(x, L, ww, x.depth.to_numpy(float))
        return b[:, 0] + 1j * b[:, 1], v


def W(z, v, idx):
    if len(idx) < MIN_N:
        return np.nan
    zz = z[idx]
    return float(np.mean(np.abs(zz - zz.mean()) ** 2) - np.mean(v[idx]))


def L_t(z, v, old, young):
    wo, wy = W(z, v, old), W(z, v, young)
    return np.log(wo / wy) if (wo > 0 and wy > 0) else np.nan


def Q_of(states_, arms):
    """states_: {(cohort, tissue): (z, v)}; arms: {(cohort, tissue): (old idx, young idx)} -> (Q, per-tissue L)."""
    L = {k: L_t(*states_[k], *arms[k]) for k in arms}
    return float(np.nanmedian(list(L.values()))), L


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    cohorts = ML.load_cohorts()
    plans, _, _, _, M, Bal = ML.build(C.DESIGN, cohorts)
    say("matched plans:", {c: len(p) for c, p in plans.items()})
    # genome per cohort (A231 filters), companion draws per cohort, then only the needed genes kept
    co, sets = {}, {}
    for c, h in ML.COH.items():
        if c not in plans:
            continue
        pb, _ = gtex.load(tissues=list(plans[c]), hardy=h, genes="all", depth="library_size", require_tod=(c == "fast"))
        genes = [g for g in pb.columns if g not in ("donor", "compartment", "depth")]
        x = pb[genes].to_numpy(np.float32)
        lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
        lm = lm[~lm.index.duplicated()]
        del x
        sets[c] = {kind: pseudo_waveforms(lm, N_SETS, np.random.default_rng([373, tag]), WF, kind=kind)
                   for kind, tag in (("matched", 1), ("responsive", 2))}
        keep = sorted({g for ws in sets[c].values() for w in ws for g in w.genes} | set(WF.genes))
        co[c] = Cohort(pb[["donor", "compartment", "depth"] + keep].reset_index(drop=True))
        del pb
        say(c, "genome loaded; kept", len(keep), "genes", f"({time.time() - t0:.0f}s)")
    # matched arms (row indices per tissue) and pooled arms
    arms, pooled, pairs = {}, {}, {}
    for c, plan in plans.items():
        q = cohorts[c][1]
        old_pool, young_pool = C.pools(q)
        for t, ss in plan.items():
            if t not in co[c].pos:
                continue
            pos = co[c].pos[t]
            pr = [(pos[str(s[0])], pos[str(s[1])]) for s in ss if str(s[0]) in pos and str(s[1]) in pos]
            pairs[(c, t)] = np.array(pr)
            arms[(c, t)] = (pairs[(c, t)][:, 0], pairs[(c, t)][:, 1])
            ids = np.array(list(pos))
            o = np.array([pos[d] for d in ids if d in q.index and old_pool.get(d, False)])
            y = np.array([pos[d] for d in ids if d in q.index and young_pool.get(d, False)])
            pooled[(c, t)] = (o, y)
    bal = Bal.assign() if isinstance(Bal, pd.DataFrame) else pd.DataFrame()
    say("tissues:", len(arms))

    def states_for(w_by_c):
        return {k: co[k[0]].zv(w_by_c[k[0]], k[1]) for k in arms}

    if a.stage == "power":
        ML.balance_summary(Bal, M).to_csv(OUT / "balance.csv", index=False)
        clock = states_for({c: WF for c in co})
        rng = np.random.default_rng([373, 4])
        null = []
        for _ in range(N_PERM):
            sw = {}
            for k, p in pairs.items():
                f = rng.random(len(p)) < 0.5
                o = np.where(f, p[:, 1], p[:, 0])
                y = np.where(f, p[:, 0], p[:, 1])
                sw[k] = (o, y)
            null.append(Q_of(clock, sw)[0])
        null = np.array(null)
        pd.DataFrame(dict(Q=null)).to_csv(OUT / "null.csv", index=False)
        rows = []
        for kind in ("matched", "responsive"):
            for j in range(N_SETS):
                st = states_for({c: sets[c][kind][j] for c in co})
                rows.append(dict(kind=kind, set=j, Q=Q_of(st, arms)[0], Q_pooled=Q_of(st, pooled)[0]))
            say(kind, "companions done", f"({time.time() - t0:.0f}s)")
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "companions.csv", index=False)
        sd = float(np.std(null))
        res = dict(tissues=len(arms), pairs_median=float(np.median([len(p) for p in pairs.values()])),
                   label_null_sd=sd, label_null_q05=float(np.quantile(null, 0.05)), expected_half_width=1.645 * sd,
                   margin=MARGIN, matched_median=float(comp[comp.kind == "matched"].Q.median()),
                   matched_q05=float(comp[comp.kind == "matched"].Q.quantile(0.05)),
                   responsive_median=float(comp[comp.kind == "responsive"].Q.median()),
                   responsive_q05=float(comp[comp.kind == "responsive"].Q.quantile(0.05)),
                   matched_median_pooled=float(comp[comp.kind == "matched"].Q_pooled.median()),
                   responsive_median_pooled=float(comp[comp.kind == "responsive"].Q_pooled.median()))
        res["gate_H373.2"] = bool(res["expected_half_width"] < MARGIN)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    comp = pd.read_csv(OUT / "companions.csv")
    null = pd.read_csv(OUT / "null.csv").Q.to_numpy()
    out, trows = [], []
    for gs, w in (("14 genes", WF), ("core", CORE)):
        st = states_for({c: w for c in co})
        q, L = Q_of(st, arms)
        qp, Lp = Q_of(st, pooled)
        rng = np.random.default_rng([373, 3])
        bt = []
        for _ in range(N_BOOT):
            arms_b = {}
            for k, p in pairs.items():
                i = rng.integers(0, len(p), len(p))
                arms_b[k] = (p[i, 0], p[i, 1])
            bt.append(Q_of(st, arms_b)[0])
        lo, hi = np.quantile(bt, [0.05, 0.95])
        mt, rs = comp[comp.kind == "matched"].Q, comp[comp.kind == "responsive"].Q
        r = dict(gene_set=gs, Q=q, amp_ratio=float(np.exp(q / 2)), lo90=float(lo), hi90=float(hi),
                 amp_ratio_lo90=float(np.exp(lo / 2)), amp_ratio_hi90=float(np.exp(hi / 2)),
                 Q_star=q - float(mt.median()), p_label=float((1 + np.sum(null <= q)) / (1 + len(null))),
                 pct_matched=float((mt < q).mean() * 100), pct_responsive=float((rs < q).mean() * 100),
                 Q_pooled=qp, amp_ratio_pooled=float(np.exp(qp / 2)),
                 Q_pooled_star=qp - float(comp[comp.kind == "matched"].Q_pooled.median()))
        for c in co:
            r[f"Q_{c}"] = float(np.nanmedian([v for k, v in L.items() if k[0] == c]))
            r[f"Q_pooled_{c}"] = float(np.nanmedian([v for k, v in Lp.items() if k[0] == c]))
        if gs == "14 genes":
            r["H373.1"] = bool(q < float(mt.quantile(.05)) and q < float(np.quantile(null, .05)) and hi < 0)
            r["H373.2_decidable"] = str(pw["gate_H373.2"]) == "True"
            r["H373.2"] = bool(lo > -MARGIN and hi < MARGIN) if r["H373.2_decidable"] else None
        out.append(r)
        trows += [dict(gene_set=gs, cohort=k[0], tissue=k[1], pairs=len(pairs[k]), L_matched=L[k], L_pooled=Lp.get(k))
                  for k in arms]
    pd.DataFrame(out).to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(trows).to_csv(OUT / "tissues.csv", index=False)
    print(pd.DataFrame(out).round(4).T.to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
