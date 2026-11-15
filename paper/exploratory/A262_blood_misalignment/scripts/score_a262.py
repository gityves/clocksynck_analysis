"""A262: in human blood, does the genome-wide rhythmic program follow central time or the displaced sleep schedule?
(GSE48113 forced desynchrony; GSE122725 dim-light arm = GSE107537 simulated night shift; reference GSE39445)

Registered in paper/docs/registered_experiments.md (A262). Loaders and the reference cosinor are the scout's
(`scouts_2026-11-15/T5_blood_misalignment/scout_blood.py`). Per subject and condition, z = regression of the program's
gene-centred values on A_g (cos, sin)(w (t - phi_g)); hypotheses excluded at Hotelling p < 0.01; rhythm gate against 200
intensity-matched non-rhythmic sets.

    python -I score_a262.py
Outputs (output/): programs.csv (one row per dataset x program x central hypothesis), matched.csv, tests.csv.

Choices not fixed by the registration
- Clock-time mapping of the reference: circular mean of (clock time - circadian phase) over the GSE39445 sleep-extension
  samples used in the cosinor.
- Intensity deciles: mean log2 over all arrays of the dataset, over the genes shared by reference and dataset.
- Matched draw per decile: a random permutation of the decile's non-rhythmic genes, completed with replacement when the
  program needs more genes than the decile holds.
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
sys.path.insert(0, str(PAPER / "exploratory/scouts_2026-11-15/T5_blood_misalignment"))
import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from scout_blood import DATA, gpl15331, load_archer, load_kervezee, loading, reference, series  # noqa: E402

OUT = HERE.parent / "output"
W = 2 * np.pi / 24
SEED, N_SETS, ALPHA = 262, 200, 0.01
CLOCK = {canonical(g) for g in cs.gtex().genes}
DESIGNS = {
    "D1 GSE48113": (load_archer, {"central": (0.0, 1.0), "behaviour": (12.0, 1.0), "lost": (0.0, 0.0)}, [0.0]),
    "D2 GSE122725 dim light": (load_kervezee, {"central": (1.5, 1.0), "behaviour": (10.0, 1.0), "lost": (0.0, 0.0)},
                               [1.5, 0.0, 2.0]),
}


def clock_offset() -> float:
    smp, _ = series(DATA / "blood_timed/GSE39445/GSE39445_series_matrix.txt.gz")
    smp = smp[smp.sleepprotocol == "Sleep Extension"].dropna(subset=["circadianphase"])
    t = np.array([int(h) + int(m) / 60 for h, m in (v.split(":") for v in smp.timesampletaken)])
    d = W * (t - smp.circadianphase.astype(float).to_numpy())
    return float(np.mod(np.angle(np.exp(1j * d).mean()), 2 * np.pi) / W)


def read(zA: np.ndarray, zM: np.ndarray, hyp: dict) -> dict:
    n = len(zA)
    q0 = np.log(zM.mean() / zA.mean())
    jk = np.array([np.log(np.delete(zM, i).mean() / np.delete(zA, i).mean()) for i in range(n)])
    jk = jk.real + 1j * (q0.imag + np.angle(np.exp(1j * (jk.imag - q0.imag))))
    se = lambda v: np.sqrt((n - 1) / n * ((v - v.mean()) ** 2).sum())
    t = stats.t.ppf(0.95, n - 1)
    out = dict(n=n, z_al_abs=float(abs(zA.mean())), z_al_arg_h=float(np.angle(zA.mean()) / W),
               z_mis_abs=float(abs(zM.mean())), shift_h=q0.imag / W, shift_lo=(q0.imag - t * se(jk.imag)) / W,
               shift_hi=(q0.imag + t * se(jk.imag)) / W, ratio=float(np.exp(q0.real)),
               ratio_lo=float(np.exp(q0.real - t * se(jk.real))), ratio_hi=float(np.exp(q0.real + t * se(jk.real))))
    for h, (s, r) in hyp.items():
        dd = zM - r * np.exp(1j * W * s) * zA
        X = np.c_[dd.real, dd.imag]
        m, Sg = X.mean(0), np.cov(X.T)
        T2 = n * m @ np.linalg.solve(Sg, m)
        out[f"p_{h}"] = float(stats.f.sf((n - 2) / (2 * (n - 1)) * T2, 2, n - 2))
    kept = [h for h in hyp if out[f"p_{h}"] >= ALPHA]
    out["cls"] = kept[0] if len(kept) == 1 else ("mixed" if not kept else "unresolved")
    out["not_excluded"] = "+".join(kept)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    ref = reference(gpl15331())
    off = clock_offset()
    ref = ref.assign(phi=np.mod(ref.acro_h + off, 24))
    ref = ref[~ref.gene.map(canonical).isin(CLOCK)].set_index("gene")
    print(f"reference: {len(ref)} genes, {(ref.p < 0.01).sum()} rhythmic; clock - phase offset {off:.2f} h", flush=True)
    prows, mrows, tests = [], [], []
    for name, (loader, hyp, centrals) in DESIGNS.items():
        smp, y = loader()
        smp = smp[smp.gsm.isin(y.columns)]
        units = {}
        for (s, c), g in smp.groupby(["subject", "cond"]):
            units.setdefault(s, {})[c] = g
        units = {s: u for s, u in units.items() if len(u) == 2 and min(len(v) for v in u.values()) >= 3}
        genes = [g for g in ref.index if g in y.index]
        r = ref.loc[genes]
        Y = {s: {c: (lambda v: v - v.mean(1).to_numpy()[:, None])(y.loc[genes, g.gsm]) for c, g in u.items()}
             for s, u in units.items()}
        T = {s: {c: g.t.to_numpy(float) for c, g in u.items()} for s, u in units.items()}
        dec = pd.qcut(y.loc[genes, smp.gsm].mean(1), 10, labels=False)
        rhy = r.index[r.p < 0.01]
        non = r.index[r.p > 0.5]
        by_dec = {k: np.array(non[(dec.reindex(non) == k).to_numpy()]) for k in range(10)}
        print(name, "subjects", len(units), "genes", len(genes), "rhythmic", len(rhy), "non-rhythmic", len(non),
              flush=True)

        def program(gs):
            gs = list(gs)
            amp, phi = r.loc[gs, "amp"].to_numpy(), r.loc[gs, "phi"].to_numpy()
            zA = np.array([loading(Y[s]["al"].loc[gs].to_numpy().T, T[s]["al"], amp, phi) for s in units])
            zM = np.array([loading(Y[s]["mis"].loc[gs].to_numpy().T, T[s]["mis"], amp, phi) for s in units])
            return zA, zM

        def matched(gs):
            need = dec.loc[gs].value_counts()
            pick = []
            for k, n_k in need.items():
                pool = by_dec[int(k)]
                p = rng.permutation(pool)[:n_k]
                pick += list(p) + list(rng.choice(pool, n_k - len(p)))
            return pick

        progs = {"ALL_RHYTHMIC": list(rhy)}
        for b in range(6):
            progs[f"PHASE_{4 * b:02d}_{4 * b + 4:02d}h"] = list(rhy[((r.loc[rhy].phi // 4) == b).to_numpy()])
        for pname, gs in progs.items():
            zA, zM = program(gs)
            null = []
            for i in range(N_SETS):
                nA, nM = program(matched(gs))
                rr = read(nA, nM, hyp)
                null.append(rr["z_al_abs"])
                mrows.append(dict(design=name, program=pname, set=i, G=len(gs), z_al_abs=rr["z_al_abs"], cls=rr["cls"]))
            q95 = float(np.quantile(null, 0.95))
            for sc in centrals:
                h = hyp | {"central": (sc, 1.0)}
                o = read(zA, zM, h)
                prows.append(dict(design=name, program=pname, central_h=sc, G=len(gs), **o, matched_q95=q95,
                                  gated=bool(o["z_al_abs"] > q95)))
                print(name, pname, sc, len(gs), o["cls"], round(o["shift_h"], 2), round(o["ratio"], 3),
                      "gated" if o["z_al_abs"] > q95 else "not gated", flush=True)
        P = pd.DataFrame(prows)
        a = P[(P.design == name) & (P.program == "ALL_RHYTHMIC") & (P.central_h == centrals[0])].iloc[0]
        tests.append(dict(test="H262.1", design=name, readout="all rhythmic gated, aligned argument within 3 h",
                          value=float(a.z_al_arg_h), gated=bool(a.gated),
                          passed=bool(a.gated and abs(a.z_al_arg_h) <= 3)))
        tests.append(dict(test="H262.2" if name.startswith("D1") else "H262.3", design=name,
                          readout="all rhythmic: class", value=a.cls, shift_h=a.shift_h, ratio=a.ratio,
                          not_excluded=a.not_excluded))
    pd.DataFrame(prows).to_csv(OUT / "programs.csv", index=False)
    pd.DataFrame(mrows).to_csv(OUT / "matched.csv", index=False)
    Tt = pd.DataFrame(tests)
    Tt.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(Tt.to_string(index=False))


if __name__ == "__main__":
    main()
