"""A371 (registered): do living organs keep the post-mortem organ arrangement?

Pair phase phi_ab = arg sum_p z_a(p) conj z_b(p) (clock state shared by two organs of the same persons). Living: STAGE
(liver, muscle, visceral fat, aortic wall; A369 prep) and AAGMEx (subcutaneous fat, muscle; A370 prep), A255 `Cohort`
states. Post-mortem: GTEx fast deaths (A336 rows), `clocksync.paired` states per tissue. Agreement A = mean over the six
STAGE pairs of cos(phi_STAGE - phi_GTEx); label null over the 24 assignments; matched and responsive companions.

    python a371_living_order.py --stage power [--jobs 4]   # GTEx phases, simulations, companions; no living phase
    python a371_living_order.py --stage score              # registered tests
Outputs (`../output/`): gtex_pairs.csv, power.csv, sim.csv, companions.csv | tests.csv, pairs.csv.
"""
import argparse
import itertools
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
for p in ("A366_muther_timing", "A369_stage_volume", "A370_aagmex_volume", "A336_two_readouts"):
    sys.path.insert(0, str(EXP / p / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a336_person_gtex as G  # noqa: E402
import a366_muther_timing as M  # noqa: E402
import a369_stage as S  # noqa: E402
import a370_aagmex as AG  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
STAGE4 = ["liver", "muscle", "visceral_fat", "aortic_wall"]
GTEX = {"liver": "Liver", "muscle": "Muscle - Skeletal", "visceral_fat": "Adipose - Visceral (Omentum)",
        "aortic_wall": "Artery - Aorta", "adipose": "Adipose - Subcutaneous", "mammary_coronary": "Artery - Coronary",
        "mammary_tibial": "Artery - Tibial"}
PAIRS = list(itertools.combinations(STAGE4, 2))
PERMS = list(itertools.permutations(range(4)))
N_SIM, N_SETS, N_BOOT, N_BOOT_SIM = 500, 200, 2000, 500
H = 24 / (2 * np.pi)


def pair_phase(za, zb):
    ok = np.isfinite(za) & np.isfinite(zb)
    return float(np.angle(np.sum(za[ok] * np.conj(zb[ok]))))


def cdist(a, b):
    return float(np.angle(np.exp(1j * (a - b))))


def agreement(ph_live, ph_ref, perm=range(4)):
    """mean cos over the six pairs; perm relabels the living tissues (live tissue perm[i] read as STAGE4[i])."""
    out = []
    for i, j in itertools.combinations(range(4), 2):
        pl = ph_live[(STAGE4[perm[i]], STAGE4[perm[j]])] if perm[i] < perm[j] else -ph_live[(STAGE4[perm[j]], STAGE4[perm[i]])]
        out.append(np.cos(pl - ph_ref[(STAGE4[i], STAGE4[j])]))
    return float(np.mean(out))


def label_rank(ph_live, ph_ref):
    vals = np.array([agreement(ph_live, ph_ref, p) for p in PERMS])
    return vals[0], int(np.sum(vals >= vals[0]))   # PERMS[0] is the identity; rank 1 = first, ties against


class GTExStates:
    def __init__(self):
        pb, _ = G.load()
        self.pb = pb
        miss = [v for v in GTEX.values() if v not in set(pb.compartment)]
        assert not miss, miss

    def Z(self, w, tissues):
        out = {}
        for t in tissues:
            x = self.pb[self.pb.compartment == GTEX[t]]
            k = [i for i, c in enumerate(w.genes) if c in x.columns]
            ww = w if len(k) == len(w.genes) else Waveform(tuple(w.genes[i] for i in k), np.asarray(w.amp)[k],
                                                           np.asarray(w.psi)[k], None, None, w.name)
            L = projector(x[list(ww.genes)], ww)
            b, _ = states(x, L, ww, x.depth.to_numpy(float))
            z = pd.Series(b[:, 0] + 1j * b[:, 1], index=x.donor.astype(str).to_numpy())
            out[t] = z - z.mean()
        return out


def gtex_phase(Zg, a, b):
    j = Zg[a].index.intersection(Zg[b].index)
    return pair_phase(Zg[a][j].to_numpy(), Zg[b][j].to_numpy()), len(j)


def sim_live(n, sd, f, phi, rng):
    """persons on a ring (phase N(0, sd)), tissue t = T exp(i phi_t) + noise at clock fraction f_t."""
    T = np.exp(1j * rng.normal(0, sd, n))
    T = T - T.mean()
    v = np.mean(np.abs(T) ** 2)
    return {t: T * np.exp(1j * phi[t]) + np.sqrt(v * (1 - f[t]) / f[t] / 2) * (rng.normal(size=n) + 1j * rng.normal(size=n))
            for t in f}


def live_phases(Z, idx=None, tissues=STAGE4):
    """pair phases over the persons with both tissues of each pair (idx unused; kept for call symmetry)."""
    return {(a, b): pair_phase(Z[a], Z[b]) for a, b in itertools.combinations(tissues, 2)}


def boot_phase(za, zb, rng, nb):
    n = len(za)
    i = rng.integers(0, n, (nb, n))
    return np.angle(np.sum(za[i] * np.conj(zb[i]), axis=1))


def interval_excludes_opposite(bt, ref):
    d = np.angle(np.exp(1j * (bt - ref)))
    lo, hi = np.quantile(d, [0.05, 0.95])
    return bool(lo > -np.pi + 1e-9 and hi < np.pi - 1e-9 and (hi - lo) < 2 * np.pi - 1e-9 and not (lo <= np.pi <= hi)), lo, hi


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    _, cs_ = S.load()
    _, ca = AG.load()
    Zs, Za = cs_.Z(), ca.Z()
    idx_s = np.all([np.isfinite(Zs[t]) for t in STAGE4], axis=0)
    idx_a = np.isfinite(Za["adipose"]) & np.isfinite(Za["muscle"])
    gt = GTExStates()
    print(f"STAGE {int(idx_s.sum())} with four tissues, AAGMEx {int(idx_a.sum())}; GTEx loaded ({time.time() - t0:.0f}s)", flush=True)

    # GTEx pair phases (clock; not blind)
    Zg = gt.Z(cs_.wf, STAGE4 + ["mammary_coronary", "mammary_tibial"])
    Zga = gt.Z(ca.wf, ["adipose", "muscle"])
    ref, rows = {}, []
    for x, y in PAIRS:
        ref[(x, y)], nd = gtex_phase(Zg, x, y)
        rows.append(dict(pair=f"{x}-{y}", gtex_phase_h=ref[(x, y)] * H, donors=nd))
    ref_a, nd = gtex_phase(Zga, "adipose", "muscle")
    rows.append(dict(pair="adipose-muscle (AAGMEx genes)", gtex_phase_h=ref_a * H, donors=nd))
    gp = pd.DataFrame(rows)
    gp.to_csv(OUT / "gtex_pairs.csv", index=False)

    if a.stage == "power":
        print(gp.round(3).to_string(index=False), flush=True)
        p369 = pd.read_csv(EXP / "A369_stage_volume/output/power.csv", index_col=0)["value"]
        p370 = pd.read_csv(EXP / "A370_aagmex_volume/output/power.csv", index_col=0)["value"]
        A255 = M.a255()
        f = {t: float(p369[f"fraction_{t}"]) for t in ("liver", "muscle", "visceral_fat")}
        ra = [A255.r_pair(Zs["aortic_wall"][idx_s], Zs[t][idx_s])[0] ** 2 / f[t] for t in f]   # moduli only
        f["aortic_wall"] = float(min(np.mean(ra), 0.99))
        phi = {"liver": 0.0, **{t: -ref[("liver", t)] for t in ("muscle", "visceral_fat", "aortic_wall")}}
        fa = {"adipose": float(p370["fraction_adipose"]), "muscle": float(p370["fraction_muscle"])}
        rng = np.random.default_rng([371, 1])
        ns, na = int(idx_s.sum()), int(idx_a.sum())
        tri = np.all([np.isfinite(Zs[t]) for t in ("liver", "muscle", "visceral_fat")], axis=0)
        n66, has_aorta = int(tri.sum()), np.isfinite(Zs["aortic_wall"])[tri]
        # companions (blind: no clock)
        crow = []
        for kind, tag in (("matched", 5), ("responsive", 6)):
            for j, w in enumerate(cs_.pseudo(kind, N_SETS, tag)):
                Zp, Zgp = cs_.Z(w), gt.Z(w, STAGE4)
                pl = live_phases(Zp, idx_s)
                pr = {(x, y): gtex_phase(Zgp, x, y)[0] for x, y in PAIRS}
                crow.append(dict(kind=kind, set=j, cohort="STAGE", A=agreement(pl, pr)))
            for j, w in enumerate(ca.pseudo(kind, N_SETS, tag)):
                Zp, Zgp = ca.Z(w), gt.Z(w, ["adipose", "muscle"])
                crow.append(dict(kind=kind, set=j, cohort="AAGMEx",
                                 A=float(np.cos(pair_phase(Zp["adipose"][idx_a], Zp["muscle"][idx_a]) - gtex_phase(Zgp, "adipose", "muscle")[0]))))
            print(f"{kind} companions done ({time.time() - t0:.0f}s)", flush=True)
        comp = pd.DataFrame(crow)
        comp.to_csv(OUT / "companions.csv", index=False)
        mq = comp[comp.kind == "matched"].groupby("cohort").A.quantile(0.95)
        sims = []
        for _ in range(N_SIM):
            Zsim = sim_live(n66, float(p369["phase_sd_rad"]), f, phi, rng)
            Zsim["aortic_wall"][~has_aorta] = np.nan
            A_, rank = label_rank(live_phases(Zsim, np.ones(ns, bool)), ref)
            Zsa = sim_live(na, float(p370["phase_sd_rad"]), fa, {"adipose": 0.0, "muscle": -ref_a}, rng)
            ph = pair_phase(Zsa["adipose"], Zsa["muscle"])
            ok_int, _, _ = interval_excludes_opposite(boot_phase(Zsa["adipose"], Zsa["muscle"], rng, N_BOOT_SIM), ref_a)
            sims.append(dict(A=A_, rank=rank, aag_dev_h=abs(cdist(ph, ref_a)) * H, aag_cos=float(np.cos(ph - ref_a)), aag_interval_ok=ok_int))
        sm = pd.DataFrame(sims)
        sm.to_csv(OUT / "sim.csv", index=False)
        phs = np.array([phi[t] for t in STAGE4])
        res = dict(stage_persons=ns, aagmex_persons=na, fraction_aortic_wall=f["aortic_wall"],
                   gtex_phase_spread_h=float(np.sqrt(-2 * np.log(abs(np.mean(np.exp(1j * phs)))))) * H,
                   matched_q95_stage=float(mq["STAGE"]), matched_q95_aagmex=float(mq["AAGMEx"]),
                   responsive_q95_stage=float(comp[(comp.kind == "responsive") & (comp.cohort == "STAGE")].A.quantile(.95)),
                   responsive_q95_aagmex=float(comp[(comp.kind == "responsive") & (comp.cohort == "AAGMEx")].A.quantile(.95)),
                   sim_A_median=float(sm.A.median()), sim_rank1_share=float((sm["rank"] == 1).mean()),
                   power_H371_1=float(((sm.A > mq["STAGE"]) & (sm["rank"] == 1)).mean()),
                   sim_aag_dev_h_median=float(sm.aag_dev_h.median()),
                   power_H371_2=float(((sm.aag_dev_h <= 3) & sm.aag_interval_ok & (sm.aag_cos > mq["AAGMEx"])).mean()))
        res["gate_H371.1"] = bool(res["power_H371_1"] >= 0.8)
        res["gate_H371.2"] = bool(res["power_H371_2"] >= 0.8)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    comp = pd.read_csv(OUT / "companions.csv")
    pl = live_phases(Zs, idx_s)
    A_, rank = label_rank(pl, ref)
    rng = np.random.default_rng([371, 3])
    prow = []
    for x, y in PAIRS:
        m2 = np.isfinite(Zs[x]) & np.isfinite(Zs[y])
        bt = boot_phase(Zs[x][m2], Zs[y][m2], rng, N_BOOT)
        d = np.angle(np.exp(1j * (bt - pl[(x, y)])))
        prow.append(dict(cohort="STAGE", pair=f"{x}-{y}", living_h=pl[(x, y)] * H, gtex_h=ref[(x, y)] * H,
                         diff_h=cdist(pl[(x, y)], ref[(x, y)]) * H, lo90_h=(pl[(x, y)] + np.quantile(d, .05)) * H,
                         hi90_h=(pl[(x, y)] + np.quantile(d, .95)) * H))
    for t in ("liver", "muscle", "visceral_fat", "aortic_wall"):    # mammary artery, descriptive
        m = np.isfinite(Zs["mammary_artery"]) & np.isfinite(Zs[t])
        ph = pair_phase(Zs["mammary_artery"][m], Zs[t][m])
        for gname in ("mammary_coronary", "mammary_tibial"):
            j = Zg[gname].index.intersection(Zg[t].index)
            prow.append(dict(cohort="STAGE", pair=f"mammary_artery-{t} vs {GTEX[gname]}", living_h=ph * H,
                             gtex_h=pair_phase(Zg[gname][j].to_numpy(), Zg[t][j].to_numpy()) * H))
    pa = pair_phase(Za["adipose"][idx_a], Za["muscle"][idx_a])
    bta = boot_phase(Za["adipose"][idx_a], Za["muscle"][idx_a], rng, N_BOOT)
    ok_int, lo, hi = interval_excludes_opposite(bta, ref_a)
    prow.append(dict(cohort="AAGMEx", pair="adipose-muscle", living_h=pa * H, gtex_h=ref_a * H, diff_h=cdist(pa, ref_a) * H,
                     lo90_h=(ref_a + lo) * H, hi90_h=(ref_a + hi) * H))
    pd.DataFrame(prow).to_csv(OUT / "pairs.csv", index=False)
    cm = comp[comp.cohort == "STAGE"]
    cmA = comp[comp.cohort == "AAGMEx"]
    t1 = dict(test="H371.1", A=A_, label_rank=rank, label_p=rank / 24, matched_q95=float(cm[cm.kind == "matched"].A.quantile(.95)),
              pct_matched=float((cm[cm.kind == "matched"].A < A_).mean() * 100),
              pct_responsive=float((cm[cm.kind == "responsive"].A < A_).mean() * 100),
              powered=str(pw["gate_H371.1"]) == "True")
    t1["passes"] = bool(A_ > t1["matched_q95"] and rank == 1) if t1["powered"] else None
    cos_a = float(np.cos(pa - ref_a))
    t2 = dict(test="H371.2", A=cos_a, dev_h=abs(cdist(pa, ref_a)) * H, interval_excludes_opposite=ok_int,
              matched_q95=float(cmA[cmA.kind == "matched"].A.quantile(.95)),
              pct_matched=float((cmA[cmA.kind == "matched"].A < cos_a).mean() * 100),
              pct_responsive=float((cmA[cmA.kind == "responsive"].A < cos_a).mean() * 100),
              powered=str(pw["gate_H371.2"]) == "True")
    t2["passes"] = bool(t2["dev_h"] <= 3 and ok_int and cos_a > t2["matched_q95"]) if t2["powered"] else None
    t = pd.DataFrame([t1, t2])
    t.to_csv(OUT / "tests.csv", index=False)
    print(t.round(4).T.to_string(), "\n", pd.DataFrame(prow).round(2).to_string(index=False), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
