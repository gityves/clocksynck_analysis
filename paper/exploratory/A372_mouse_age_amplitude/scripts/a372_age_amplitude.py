"""A372 (registered): does the clock amplitude of all peripheral organs fall together with age beyond random genes?

Wolff 2023 (GSE201207; 216 libraries = pools of 2 mice; 6 tissues x 3 ages x 12 CT, constant darkness). Raw counts with
mouse symbols upper-cased (orthology by symbol, as the clock-gene mapping of A48). Timed gain per tissue and age: A48
`timed_gain` generalised to any waveform. D = mean over the five peripheral tissues of log(|gain old| / |gain young|).
Label null: young / old swapped within CT (all 4,096 assignments); companions: expression-matched and responsive
pseudo-clock sets.

    python a372_age_amplitude.py --stage power   # label null of the clock D (clock D not computed), companions, MDE
    python a372_age_amplitude.py --stage score   # registered test
Outputs (`../output/`): power.csv, companions.csv, null.csv | tests.csv, tissues.csv, drops.csv.
"""
import argparse
import itertools
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
SRC = DATA / "mouse_ageing" / "GSE201207" / "GSE201207_count_216.csv.gz"
PERIPH = ["adrenal", "heart", "kidney", "lung", "muscle"]
W = 2 * np.pi / 24
WF = cs.gtex()
TRIO = ["PER1", "PER2", "NR1D1"]
N_SETS = 200
ALT = np.log(0.8)


def load():
    c = pd.read_csv(SRC, index_col=0)
    c.index = c.index.astype(str).str.upper()
    c = c.groupby(level=0).sum()
    meta = []
    for s in c.columns:
        ct, age, tis = s.split("-")
        meta.append(dict(lib=s, ct=int(ct[2:]) % 24, ct_raw=int(ct[2:]), age=age, tissue=tis))
    meta = pd.DataFrame(meta).set_index("lib")
    cpm = c / c.sum(0) * 1e6
    return cpm.T, meta          # libraries x genes, meta


def design(h):
    t = np.asarray(h, float) * W
    return np.c_[np.ones(len(t)), np.cos(t), np.sin(t)]


def gain(cpm, h, w):
    """A48 timed_gain for waveform w: cosinor of log(CPM + 0.5) of w's genes with median CPM >= 5."""
    g = [x for x in w.genes if x in cpm.columns]
    k = [i for i, x in enumerate(w.genes) if x in cpm.columns]
    y = cpm[g].to_numpy(float)
    ex = np.median(y, 0) >= 5
    if ex.sum() < 3:
        return np.nan
    B = np.linalg.lstsq(design(h), np.log(y[:, ex] + 0.5), rcond=None)[0]
    x = (np.asarray(w.amp)[k] * np.exp(1j * np.asarray(w.psi)[k]))[ex]
    return float(abs((np.conj(x) @ (B[1] + 1j * B[2])) / (np.conj(x) @ x)))


def D_of(cpm, meta, w, swap=None, tissues=PERIPH, ages=("young", "old")):
    """mean over tissues of log(|gain ages[1]| / |gain ages[0]|); swap: CT -> bool (exchange the two ages' pools)."""
    out = []
    for t in tissues:
        gs = []
        for a in ages:
            m = meta[(meta.tissue == t) & meta.age.isin(ages)]
            lab = m.age.copy()
            if swap is not None:
                flip = m.ct_raw.map(swap).to_numpy(bool)
                lab[flip] = np.where(lab[flip] == ages[0], ages[1], ages[0])
            sel = lab.index[lab == a]
            gs.append(gain(cpm.loc[sel], meta.loc[sel, "ct"], w))
        out.append(np.log(gs[1] / gs[0]))
    return float(np.mean(out)), out


def all_swaps(meta):
    cts = sorted(meta.ct_raw.unique())   # the 12 sampling times (registered: 2^12)
    for bits in itertools.product([False, True], repeat=len(cts)):
        yield dict(zip(cts, bits))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    cpm, meta = load()
    missing = [g for g in WF.genes if g not in cpm.columns]
    print(f"{len(meta)} libraries, {cpm.shape[1]} genes; clock genes missing: {missing} ({time.time() - t0:.0f}s)", flush=True)
    lm = np.log(cpm.mean(0) + 1e-3)
    lm = lm[cpm.mean(0) >= 1]
    swaps = list(all_swaps(meta))

    if a.stage == "power":
        # label null of the clock D over the 4,096 assignments minus the identity and the full swap (both give the
        # observed D up to sign); the observed D itself is not summarised
        null = np.array([D_of(cpm, meta, WF, s)[0] for s in swaps if any(s.values()) and not all(s.values())])
        pd.DataFrame(dict(D=null)).to_csv(OUT / "null.csv", index=False)
        rows = []
        for kind, tag in (("matched", 1), ("responsive", 2)):
            for j, w in enumerate(pseudo_waveforms(lm, N_SETS, np.random.default_rng([372, 2, tag]), WF, kind=kind)):
                rows.append(dict(kind=kind, set=j, D=D_of(cpm, meta, w)[0]))
            print(f"{kind} done ({time.time() - t0:.0f}s)", flush=True)
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "companions.csv", index=False)
        sd = float(np.std(null))
        crit = float(np.quantile(null, 0.05))
        mq05 = float(comp[comp.kind == "matched"].D.quantile(0.05))
        thr = min(mq05, crit)
        from scipy.stats import norm
        res = dict(libraries=len(meta), label_assignments=len(null), label_null_sd=sd, label_null_q05=crit,
                   mde_label=2.49 * sd, matched_median=float(comp[comp.kind == "matched"].D.median()), matched_q05=mq05,
                   matched_sd=float(comp[comp.kind == "matched"].D.std()),
                   responsive_median=float(comp[comp.kind == "responsive"].D.median()),
                   responsive_q05=float(comp[comp.kind == "responsive"].D.quantile(0.05)), threshold=thr,
                   power_20pct=float(norm.cdf((thr - ALT) / sd)))
        res["gate"] = bool(res["power_20pct"] >= 0.8)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("A372 not powered: not decidable; not scored")
    comp = pd.read_csv(OUT / "companions.csv")
    null = pd.read_csv(OUT / "null.csv").D.to_numpy()
    d, per = D_of(cpm, meta, WF)
    d_aged, per_aged = D_of(cpm, meta, WF, ages=("young", "aged"))
    mt, rs = comp[comp.kind == "matched"].D, comp[comp.kind == "responsive"].D
    row = dict(test="H372.1", D=d, fold=float(np.exp(d)), p_label=float((1 + np.sum(null <= d)) / (1 + len(null))),
               threshold=float(pw["threshold"]), matched_q05=float(mt.quantile(.05)), pct_matched=float((mt < d).mean() * 100),
               responsive_q05=float(rs.quantile(.05)), pct_responsive=float((rs < d).mean() * 100), D_aged_vs_young=d_aged)
    row["passes"] = bool(d < row["matched_q05"] and row["p_label"] < 0.05)
    core = WF.drop(TRIO)
    row["D_core"] = D_of(cpm, meta, core)[0]
    dr = pd.DataFrame([dict(drop="+".join(c), trio=list(c) == TRIO, D=D_of(cpm, meta, WF.drop(list(c)))[0])
                       for c in itertools.combinations(WF.genes, 3)])
    dr.to_csv(OUT / "drops.csv", index=False)
    row.update(drops_median=float(dr.D.median()), trio_drop_pct=float((dr.D < row["D_core"]).mean() * 100))
    pd.DataFrame([row]).to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(dict(tissue=PERIPH, log_old_over_young=per, log_aged_over_young=per_aged)).to_csv(OUT / "tissues.csv", index=False)
    print(pd.Series(row).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
