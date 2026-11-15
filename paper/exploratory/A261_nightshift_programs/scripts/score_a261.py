"""A261: which daily programs of human muscle follow central time, the muscle clock or behaviour after a 12 h shift of
behaviour? (GSE106800, reference GSE108539)

Registered in paper/docs/registered_experiments.md (A261; not blind for the clock genes). Reference and array loading
are the scout's (`scouts_2026-11-15/T1_nightshift/scout_t1.py`: `reference`, `arrays`). Per subject and condition, z =
regression of the program's 19:00 minus 07:00 differences on 2 A (cos, sin)(w (19 - phi)); hypotheses central (0 h, 1),
tissue clock (2.68 h, 0.233), behaviour (12 h, 1) excluded at Hotelling p < 0.01; rhythm gate against 200 matched
non-rhythmic sets.

    python score_a261.py
Outputs (output/): programs.csv (one row per program), matched.csv (companion sets), tests.csv.

Choices not fixed by the registration
- Matched sets: for a program of G genes, G genes drawn without replacement from reference non-rhythmic genes on the
  array (p > 0.5), one per program gene from the same decile of mean array intensity (log2, all 44 arrays); seed 261.
- A decile with no unused non-rhythmic gene left (programs of hundreds of genes) draws from all unused non-rhythmic genes, and
  when none is left (a program larger than the non-rhythmic pool), from its decile with replacement.
- Class calls of the matched sets use the same tests; their distribution is reported per program size.
- Jackknife-t intervals: log ratio and angle of mean z_mis / mean z_al, leave-one-subject-out, t(n - 1).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/scouts_2026-11-15/T1_nightshift"))
import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from scout_t1 import HYP, arrays, reference  # noqa: E402

OUT = HERE.parent / "output"
W = 2 * np.pi / 24
SEED, N_SETS, ALPHA, MIN_G = 261, 200, 0.01, 20
GMT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/reference/h.all.v2024.1.Hs.symbols.gmt")
CLOCK = {canonical(g) for g in cs.gtex().genes}


def per_subject_z(D: np.ndarray, C: np.ndarray, S: np.ndarray) -> np.ndarray:
    X = np.c_[C, S]
    B, *_ = np.linalg.lstsq(X, D.T, rcond=None)
    return B[0] + 1j * B[1]


def read(zA: np.ndarray, zM: np.ndarray) -> dict:
    n = len(zA)
    q0 = np.log(zM.mean() / zA.mean())
    jk = np.array([np.log(np.delete(zM, i).mean() / np.delete(zA, i).mean()) for i in range(n)])
    jk = jk.real + 1j * (q0.imag + np.angle(np.exp(1j * (jk.imag - q0.imag))))
    se = lambda v: np.sqrt((n - 1) / n * ((v - v.mean()) ** 2).sum())
    t = stats.t.ppf(0.95, n - 1)
    out = dict(z_al_abs=float(abs(zA.mean())), z_al_arg_h=float(np.angle(zA.mean()) / W), shift_h=q0.imag / W,
               shift_lo=(q0.imag - t * se(jk.imag)) / W, shift_hi=(q0.imag + t * se(jk.imag)) / W,
               ratio=float(np.exp(q0.real)), ratio_lo=float(np.exp(q0.real - t * se(jk.real))),
               ratio_hi=float(np.exp(q0.real + t * se(jk.real))))
    for h, (s, r) in HYP.items():
        dd = zM - r * np.exp(1j * W * s) * zA
        X = np.c_[dd.real, dd.imag]
        m, Sg = X.mean(0), np.cov(X.T)
        T2 = n * m @ np.linalg.solve(Sg, m)
        out[f"p_{h}"] = float(stats.f.sf((n - 2) / (2 * (n - 1)) * T2, 2, n - 2))
    kept = [h for h in HYP if out[f"p_{h}"] >= ALPHA]
    out["cls"] = kept[0] if len(kept) == 1 else ("mixed" if not kept else "unresolved")
    out["not_excluded"] = "+".join(kept)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    ref = reference().set_index("gene")
    d = arrays()
    al, mis = d["ctl"], d["mis"]
    genes = [g for g in ref.index if g in al.columns and g not in CLOCK]
    ref = ref.loc[genes]
    rhy = ref.index[ref.p < 0.01]
    non = ref.index[ref.p > 0.5]
    from compute_a257 import load
    _, y = load()
    mean_int = (y / np.log(2.0)).mean(0)
    mean_int = mean_int[~mean_int.index.duplicated()].reindex(genes)
    dec = pd.qcut(mean_int, 10, labels=False)
    C = 2 * ref.amp * np.cos(W * (19 - ref.acro_h))
    S = 2 * ref.amp * np.sin(W * (19 - ref.acro_h))

    def program(gs):
        gs = list(gs)
        zA = per_subject_z(al[gs].to_numpy(float), C[gs].to_numpy(), S[gs].to_numpy())
        zM = per_subject_z(mis[gs].to_numpy(float), C[gs].to_numpy(), S[gs].to_numpy())
        return read(zA, zM)

    progs = {}
    for line in open(GMT):
        f = line.rstrip("\n").split("\t")
        gs = [g for g in dict.fromkeys(canonical(x) for x in f[2:]) if g in set(rhy)]
        if len(gs) >= MIN_G:
            progs[f[0]] = ("hallmark", gs)
    progs["ALL_RHYTHMIC"] = ("all", list(rhy))
    for b in range(6):
        gs = list(rhy[((ref.loc[rhy].acro_h // 4) == b).to_numpy()])
        progs[f"PHASE_{4 * b:02d}_{4 * b + 4:02d}h"] = ("phase (descriptive)", gs)

    by_dec = {k: list(non[(dec.reindex(non) == k).to_numpy()]) for k in range(10)}
    rows, mrows = [], []
    for name, (kind, gs) in progs.items():
        obs = program(gs)
        null = []
        for i in range(N_SETS):
            pick, used = [], set()
            for g in gs:
                pool = [x for x in by_dec[int(dec[g])] if x not in used]
                if not pool:  # decile exhausted (large programs): any unused non-rhythmic gene
                    pool = [x for x in non if x not in used] or list(by_dec[int(dec[g])])  # then with replacement
                x = pool[rng.integers(len(pool))]
                used.add(x)
                pick.append(x)
            r = program(pick)
            null.append(r)
            mrows.append(dict(program=name, set=i, G=len(gs), z_al_abs=r["z_al_abs"], cls=r["cls"]))
        q95 = float(np.quantile([r["z_al_abs"] for r in null], 0.95))
        rows.append(dict(program=name, kind=kind, G=len(gs), **obs, matched_q95=q95, gated=bool(obs["z_al_abs"] > q95)))
        print(name, len(gs), obs["cls"], round(obs["shift_h"], 2), round(obs["ratio"], 3), flush=True)
    P = pd.DataFrame(rows)
    P.to_csv(OUT / "programs.csv", index=False)
    pd.DataFrame(mrows).to_csv(OUT / "matched.csv", index=False)
    allr = P.set_index("program").loc["ALL_RHYTHMIC"]
    hm = P[(P.kind == "hallmark") & P.gated]
    tests = [dict(test="H261.1", readout="all rhythmic: gated and aligned argument within 3 h of 0",
                  value=float(allr.z_al_arg_h), gated=bool(allr.gated),
                  passed=bool(allr.gated and abs(allr.z_al_arg_h) <= 3))]
    for c in ("central", "tissue_clock", "behaviour", "mixed", "unresolved"):
        tests.append(dict(test="H261.2", readout=f"gated Hallmark programs: {c}", value=int((hm.cls == c).sum()),
                          share=float((hm.cls == c).mean()) if len(hm) else np.nan, n_gated=len(hm),
                          n_hallmark=int((P.kind == "hallmark").sum())))
    tests.append(dict(test="H261.3", readout="all rhythmic: class", value=allr.cls, shift_h=allr.shift_h,
                      ratio=allr.ratio))
    T = pd.DataFrame(tests)
    T.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(T.to_string(index=False))


if __name__ == "__main__":
    main()
