"""A334b (registered): after a 6 h advance, does the mouse aorta lose coherent clock loudness (rho^2), and does the
within-library dispersion D rise (cells out of phase) or fall (quieter cells)? Model-free sign test.

    python score_a334b.py --stage power   # aligned + Bmal1 KO: clock check, simulated D changes, placebo MDEs, gate
    python score_a334b.py --stage score   # misaligned vs aligned: H334b.1-2 (only if the gate passed)
Outputs (`../output/`): b_power.csv, b_sim.csv, b_placebo.csv | b_tests.csv, b_boot.csv. Public data, parsed only.
Per-group statistics (D, split-half rho^2) and the loaders are those of `score_a334.py`.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import optimize, special

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_a334 as S  # noqa: E402
from clocksync.cells import vonmises as VM  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = S.OUT
CLASSES = ("SMC", "Fibroblast")


def arm(recs, M, mes=None):
    """Cell-weighted D and rho^2 of a set of groups (mesor ZT-balanced within the set unless given)."""
    m1, m2 = S.mesor(recs) if mes is None else mes
    n = np.array([r["n"] for r in recs], float)
    D = float(np.sum(n * [r["D"] for r in recs]) / n.sum())
    r2 = float(np.sum(n * [(M.L @ (r["y1"] - m1)) @ (M.L @ (r["y2"] - m2)) for r in recs]) / n.sum())
    return D, r2


def null_D(raw, wf, seed):
    """Per group D on the 20 permuted-acrophase planes (same planes for every group): groups x 20."""
    rng = np.random.default_rng(seed)
    planes = []
    for _ in range(S.N_NULL):
        w = wf.permuted(rng)
        planes.append(VM._split_projectors(w.amp * np.exp(-1j * w.psi), S.N_SPLITS, rng))
    return np.array([[VM.group_stats(C, N, LA, rng=rng)["D"] for LA in planes] for C, N in raw])


def power(jobs):
    wf, c = S.waveform()
    M = S.Model(c, [334, 60])
    recs, raw = S.load_aorta(list(wf.genes), ["aligned", "ko"], M, [334, 61])
    rows, sims, plc = [], [], []
    for cls in CLASSES:
        al = [i for i, r in enumerate(recs) if r["cls"] == cls and r["arm"] == "aligned"]
        ko = [i for i, r in enumerate(recs) if r["cls"] == cls and r["arm"] == "ko"]
        ra = [recs[i] for i in al]
        D_al, r2_al = arm(ra, M)
        nul = S.d_null([raw[i] for i in al], wf, [334, 62])
        nko = np.array([recs[i]["n"] for i in ko], float)
        D_ko = float(np.sum(nko * [recs[i]["D"] for i in ko]) / nko.sum())
        clock_ok = bool(D_al > nul.max() and D_ko <= 0.5 * D_al)
        lam, kap, R = M.solve(D_al, r2_al)
        # simulated D changes at a rho ratio of 0.8 (pure dephasing: R x 0.8; pure damping: lambda x 0.8)
        C_all = np.vstack([raw[i][0] for i in al])
        N_all = np.concatenate([raw[i][1] for i in al])
        p = C_all.sum(0) / N_all.sum()
        ncell = int(np.median([len(raw[i][1]) for i in al]))
        zts = [recs[i]["zt"] for i in al]
        kap_d = optimize.brentq(lambda k: special.ive(1, k) / special.ive(0, k) - 0.8 * R, 1e-6, 500)

        def sim(rep):
            r = np.random.default_rng([3342, 1, rep])
            mus = np.array(zts) * 2 * np.pi / 24 + r.normal(0, 0.1, len(zts))
            out = {}
            for k, (l, kk) in dict(ref=(lam, kap), deph=(lam, kap_d), damp=(0.8 * lam, kap)).items():
                g = S.simulate(len(zts), ncell, c, p, l, kk, mus, N_all, r)
                out[k] = arm([S.record(C, N, M, r, zt=z) for (C, N), z in zip(g, zts)], M)
            return dict(cls=cls, rep=rep, **{f"dlogD_{k}": np.log(out[k][0] / out["ref"][0]) for k in ("deph", "damp")},
                        **{f"dlogrho2_{k}": np.log(max(out[k][1], 1e-12) / out["ref"][1]) for k in ("deph", "damp")})
        sm = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(sim)(k) for k in range(100)))
        sims.append(sm)
        d_deph, d_damp = float(sm.dlogD_deph.median()), float(sm.dlogD_damp.median())
        mes = S.mesor(ra)
        by = {}
        for i in al:
            by.setdefault(recs[i]["zt"], []).append(i)

        def split(s):
            r = np.random.default_rng([3342, 2, s])
            a8 = [int(i) for z, v in sorted(by.items()) for i in r.choice(v, min(2, len(v)), replace=False)]
            b6 = [i for i in al if i not in a8]
            (D8, q8), (D6, q6) = arm([recs[i] for i in a8], M, mes), arm([recs[i] for i in b6], M, mes)
            return dict(cls=cls, split=s, dlogD=np.log(D8 / D6), dlogrho2=np.log(q8 / q6))
        pl = pd.DataFrame([split(s) for s in range(200)])
        plc.append(pl)
        mde_D, mde_r2 = float(2.49 * pl.dlogD.std()), float(2.49 * pl.dlogrho2.std())
        powered = bool(mde_r2 <= np.log(1 / 0.64) and mde_D <= min(abs(d_deph), abs(d_damp)))
        rows.append(dict(cls=cls, libraries_aligned=len(al), libraries_ko=len(ko), cells_aligned=int(sum(recs[i]["n"] for i in al)),
                         clock_umi_median=float(np.median([recs[i]["clock_umi"] for i in al])), D_aligned=D_al,
                         D_null_max=float(nul.max()), D_ko=D_ko, ko_over_aligned=D_ko / D_al, rho2_aligned=r2_al,
                         sim_lam=lam, sim_kappa=kap, sim_R=R, d_deph=d_deph, d_damp=d_damp,
                         sign_right_deph=float((sm.dlogD_deph > 0).mean()), sign_right_damp=float((sm.dlogD_damp < 0).mean()),
                         MDE_dlogD=mde_D, MDE_dlogrho2=mde_r2, clock_check=clock_ok, powered=powered,
                         gate=bool(clock_ok and powered)))
        S.log(str(rows[-1]))
    pd.DataFrame(rows).to_csv(OUT / "b_power.csv", index=False)
    pd.concat(sims).to_csv(OUT / "b_sim.csv", index=False)
    pd.concat(plc).to_csv(OUT / "b_placebo.csv", index=False)


def score(jobs):
    pw = pd.read_csv(OUT / "b_power.csv").set_index("cls")
    if str(pw.loc["SMC", "gate"]) != "True":
        raise SystemExit("smooth-muscle gate not passed (clock check or power): not decidable at this design; not scored")
    wf, c = S.waveform()
    M = S.Model(c, [334, 60])
    recs, raw = S.load_aorta(list(wf.genes), ["aligned", "misaligned"], M, [334, 61])
    rows, boots = [], []
    for cls in CLASSES:
        if str(pw.loc[cls, "clock_check"]) != "True":
            rows.append(dict(cls=cls, verdict="D not clock-specific"))
            continue
        idx = [i for i, r in enumerate(recs) if r["cls"] == cls]
        al = [i for i in idx if recs[i]["arm"] == "aligned"]
        mi = [i for i in idx if recs[i]["arm"] == "misaligned"]
        zt = {i: recs[i]["zt"] for i in idx}
        zts = sorted(set(zt.values()))
        ND = null_D([raw[i] for i in idx], wf, [3342, 3])
        pos = {i: k for k, i in enumerate(idx)}
        n = {i: recs[i]["n"] for i in idx}

        def contrast(a, b):
            (Da, qa), (Db, qb) = arm([recs[i] for i in a], M), arm([recs[i] for i in b], M)
            wa, wb = np.array([n[i] for i in a], float), np.array([n[i] for i in b], float)
            na = (wa[:, None] * ND[[pos[i] for i in a]]).sum(0) / wa.sum()
            nb = (wb[:, None] * ND[[pos[i] for i in b]]).sum(0) / wb.sum()
            return dict(dlogD=np.log(Db / Da), dlogrho2=np.log(max(qb, 1e-12) / qa), D_a=Da, D_b=Db, rho2_a=qa, rho2_b=qb,
                        null_dlogD=np.log(nb / na))
        obs = contrast(al, mi)

        def perm(s):
            r = np.random.default_rng([3342, 4, s])
            a, b = [], []
            for z in zts:
                v = [i for i in idx if zt[i] == z]
                pick = set(int(x) for x in r.choice(v, sum(zt[i] == z for i in mi), replace=False))
                b += [i for i in v if i in pick]
                a += [i for i in v if i not in pick]
            d = contrast(a, b)
            return d["dlogD"], d["dlogrho2"]

        def boot(s):
            r = np.random.default_rng([3342, 5, s])
            a = [int(i) for z in zts for i in r.choice([j for j in al if zt[j] == z], sum(zt[j] == z for j in al))]
            b = [int(i) for z in zts for i in r.choice([j for j in mi if zt[j] == z], sum(zt[j] == z for j in mi))]
            d = contrast(a, b)
            return dict(cls=cls, boot=s, dlogD=d["dlogD"], dlogrho2=d["dlogrho2"])
        pm = np.array(Parallel(n_jobs=jobs)(delayed(perm)(s) for s in range(2000)))
        bt = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(boot)(s) for s in range(2000)))
        boots.append(bt)
        p_r2 = float((1 + np.sum(pm[:, 1] <= obs["dlogrho2"])) / 2001)
        p_up = float((1 + np.sum(pm[:, 0] >= obs["dlogD"])) / 2001)
        p_down = float((1 + np.sum(pm[:, 0] <= obs["dlogD"])) / 2001)
        q05, q95 = np.quantile(obs["null_dlogD"], [0.05, 0.95])
        t = dict(cls=cls, primary=cls == "SMC", powered=str(pw.loc[cls, "powered"]) == "True",
                 libraries_aligned=len(al), libraries_misaligned=len(mi), D_aligned=obs["D_a"], D_misaligned=obs["D_b"],
                 rho2_aligned=obs["rho2_a"], rho2_misaligned=obs["rho2_b"], dlogrho2=obs["dlogrho2"], p_rho2_falls=p_r2,
                 dlogD=obs["dlogD"], p_D_rises=p_up, p_D_falls=p_down, null_planes_dlogD_q05=float(q05),
                 null_planes_dlogD_q95=float(q95),
                 **{f"{k}_{q}": float(bt[k].quantile(v)) for k in ("dlogD", "dlogrho2") for q, v in (("lo90", 0.05), ("hi90", 0.95))})
        h1 = obs["dlogrho2"] < 0 and p_r2 < 0.05
        t["H334b.1"] = "PASS" if h1 else "FAIL"
        if not h1:
            t["H334b.2"] = "not read (H334b.1 failed)"
        elif obs["dlogD"] > 0 and p_up < 0.05 and obs["dlogD"] > q95:
            t["H334b.2"] = "dephasing"
        elif obs["dlogD"] < 0 and p_down < 0.05 and obs["dlogD"] < q05:
            t["H334b.2"] = "damping"
        else:
            t["H334b.2"] = "not separated"
        if not t["powered"]:
            t["note"] = "gate not passed: descriptive"
        rows.append(t)
        S.log(str(t))
    pd.DataFrame(rows).to_csv(OUT / "b_tests.csv", index=False)
    if boots:
        pd.concat(boots).to_csv(OUT / "b_boot.csv", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    {"power": power, "score": score}[a.stage](a.jobs)
    S.log("done")


if __name__ == "__main__":
    main()
