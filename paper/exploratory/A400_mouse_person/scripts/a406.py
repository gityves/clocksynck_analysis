"""A406 (registered 2026-11-27): does a within-time person factor of mouse clock loudness survive a mouse waveform
learned in other mice, and can the A400c cross-fit produce one from noise?

External waveforms: E (Wolff 2023, GSE201207: adrenal, heart, kidney, lung, muscle; group tissue x age; hours = CT) and
E2 (Manella 2021, Alb-Cre controls, liver / WAT / lung; group tissue x arm; hours = ZT), `fit_waveform` on log
((count + 0.5) / depth) of the 14 genes. Statistic, companions, cross-fit and permutations as A400c (`a400c.py`).

Stages:
  sim    Part 1 (Deota design simulated with real depths, A394 tissue gains and offsets, A400c Mf gene waveform as truth,
         recorded ZT + N(0, 1 h) per mouse, person factor exp(N(0, 0.25)) or none; read with Mx (cross-fit refitted on
         the simulated mice), E and E2), the real companions read with E and E2 (no clock value), the gates G406.1 and
         G406.2 and the Mx-null row. No real clock gene is read with E or E2.
  score  H406.1 (real Deota clock genes read with E; E2 beside it), D406 (per feeding arm). BLIND until the user
         decides: not run by the sim stage.
  check  1 replicate per cell, 4 companions, 20 permutations / bootstraps; nothing written.

    RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A400_mouse_person/scripts/a406.py --stage sim --jobs 4
Outputs (`../output/`): d406_waveforms.csv, d406_sim.csv, d406_companions.csv, d406_gates.csv; score: d406_scores.csv,
d406_arms.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
for p in (PAPER, HERE.parent, EXP / "A385_sheep_person" / "scripts", EXP / "A394_mouse_organisation" / "scripts"):
    sys.path.insert(0, str(p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a385b as R  # noqa: E402
import a394 as M  # noqa: E402
import a400c as C  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.waveform import fit_waveform  # noqa: E402
from common.datasets.manella2021.build import build as build_manella  # noqa: E402
from common.datasets.mouse_ageing_wolff.build import PERIPHERAL, build as build_wolff  # noqa: E402

OUT = HERE.parents[1] / "output"
A394_OUT = EXP / "A394_mouse_organisation" / "output"
VARIANTS, V2 = C.VARIANTS, C.VARIANTS[2]
PERSON_SD = (0.0, 0.25)
READS = ("Mx", "E", "E2")
N_REP, N_SETS, N_PERM, N_BOOT = 20, 100, 1000, 1000
S = C.S  # shared with a400c (its states / cross_states / leverage read this dict)


# ---- external waveforms ----------------------------------------------------------------------------------------------
def log_frac(d):
    return np.log((d[list(R.WF.genes)].to_numpy(float) + 0.5) / d.depth.to_numpy(float)[:, None])


def external_waveforms():
    """E (Wolff, peripheral tissues, group tissue x age, hours CT) and E2 (Manella Alb-Cre, group tissue x arm, ZT)."""
    pb, donors = build_wolff()
    pb = pb[pb.compartment.isin(PERIPHERAL)]
    ct = donors.ct.reindex(pb.donor).to_numpy(float)
    age = donors.age.reindex(pb.donor).to_numpy()
    e = fit_waveform(log_frac(pb), ct, R.WF.genes, group=pb.compartment.astype(str).to_numpy() + "|" + age,
                     name="E")
    pb, donors = build_manella()
    pb = pb[pb.stratum.str.startswith("AlbCre_")]
    zt = donors.zt.reindex(pb.donor).to_numpy(float)
    e2 = fit_waveform(log_frac(pb), zt, R.WF.genes, group=pb.compartment.astype(str) + "|" + pb.stratum.astype(str),
                      name="E2")
    return {"E": e, "E2": e2}


def waveform_table():
    """Gene-level E, E2 and Mf waveforms (amplitudes, phases) and their differences from Mf (descriptive)."""
    def centred(w, ref):
        ra = np.log(np.asarray(w.amp) / np.mean(w.amp)) - np.log(np.asarray(ref.amp) / np.mean(ref.amp))
        dp = np.angle(np.exp(1j * (np.asarray(w.psi) - ref.psi)))
        dp = np.angle(np.exp(1j * (dp - np.angle(np.mean(np.exp(1j * dp))))))
        return ra, dp * 24 / (2 * np.pi)
    out = pd.DataFrame(dict(gene=R.WF.genes, core=[g not in C.TRIO for g in R.WF.genes]))
    for k, w in (("Mf", S["wf"]["Mf"]), ("E", S["wf"]["E"]), ("E2", S["wf"]["E2"])):
        out[f"{k}_amp"], out[f"{k}_psi"] = w.amp, w.psi
    for k in ("E", "E2"):
        out[f"{k}_vs_Mf_log_rel_amp"], out[f"{k}_vs_Mf_phase_h"] = centred(S["wf"][k], S["wf"]["Mf"])
    return out


# ---- reading ---------------------------------------------------------------------------------------------------------
def read_states(read, rd, name, pb=None):
    """States of the mice read with Mx (cross-fit on `pb`, as A400c), E or E2 (fixed external waveform)."""
    if pb is not None:
        S["mouse"] = pb
    if read == "Mx":
        S["wf_fold"] = [C.fit_on(S["folds"][f]) for f in (0, 1)]
        return C.cross_states(R.WF.genes, rd, name)
    st = C.states(C.as_set(S["wf"][read], R.WF.genes, rd, name), S["mice"])
    return [(j, a, b, np.column_stack([cov, np.zeros(len(j))])) for j, a, b, cov in st]


# ---- simulation (Part 1) ---------------------------------------------------------------------------------------------
def sim_setup(pb):
    """Real depths, mean log fraction and SD per tissue x arm, A394 gain and offset per tissue (gain: sample-weighted
    mean 1, as `fit_waveform`)."""
    d = pb[M.KEYS + list(R.WF.genes)].reset_index(drop=True)
    lf = pd.DataFrame(log_frac(d), columns=list(R.WF.genes))
    key = (d.compartment.astype(str) + "|" + d.stratum.astype(str)).to_numpy()
    fit = pd.read_csv(A394_OUT / "mouse_fit.csv").set_index("tissue")
    g = np.exp(fit.log_amp_rel.reindex(d.compartment).to_numpy(float))
    assert not np.isnan(g).any(), "tissue missing from the A394 fit"
    g = g / g.mean()
    phi = 2 * np.pi * fit.offset_h.reindex(d.compartment).to_numpy(float) / 24
    return dict(d=d, depth=d.depth.to_numpy(float), g=g, phi=phi,
                m=lf.groupby(key).transform("mean").to_numpy(),
                s=lf.groupby(key).transform("std").fillna(lf.groupby(d.compartment.to_numpy()).transform("std")).to_numpy())


def simulate(sd, cell, r):
    """One simulated replicate: Poisson counts, mean depth x exp(m + lam_s g_k A_j cos(t_s - psi_j - phi_k) + e)."""
    X = S["sim"]
    rng = np.random.default_rng([406, cell, r])
    wf = S["wf"]["Mf"]
    mice = np.asarray(S["mice"])
    dep = pd.Series(rng.normal(0, 1, len(mice)), index=mice)  # individual departure from the recorded ZT (h)
    lam = pd.Series(np.exp(rng.normal(0, sd, len(mice))) if sd > 0 else 1.0, index=mice)
    don = X["d"].donor.astype(str)
    t = 2 * np.pi * (S["zt"].reindex(don).to_numpy(float) + dep.reindex(don).to_numpy()) / 24
    amp = X["g"][:, None] * wf.amp[None, :]
    osc = lam.reindex(don).to_numpy()[:, None] * amp * np.cos(t[:, None] - wf.psi[None, :] - X["phi"][:, None])
    sde = np.sqrt(np.maximum(X["s"] ** 2 - amp ** 2 / 2, X["s"] ** 2 / 4))
    mu = X["m"] + osc + rng.normal(size=osc.shape) * sde
    out = X["d"][M.KEYS].copy()
    out[list(R.WF.genes)] = rng.poisson(X["depth"][:, None] * np.exp(mu))
    return out


def sim_task(job):
    warnings.filterwarnings("ignore")
    sd, cell, r = job
    rows = []
    try:
        pb = simulate(sd, cell, r)
        for read in READS:
            for rd in ("core", "14 genes"):
                st = read_states(read, rd, f"sim{cell}_{r}", pb)
                for v in VARIANTS:
                    o = dict(person_sd=sd, rep=r, read=read, readout=rd, variant=v, leverage=C.leverage(st, v))
                    if rd == "core" and v == V2:
                        rng = np.random.default_rng([406, 9, cell, r, READS.index(read)])
                        null = np.array([C.leverage(st, v, rng) for _ in range(S["n_perm"])])
                        o.update(p_perm=float((1 + np.sum(null >= o["leverage"])) / (1 + len(null))),
                                 null_q95=float(np.quantile(null, 0.95)))
                    rows.append(o)
    except Exception as e:  # counted, not dropped silently
        rows.append(dict(person_sd=sd, rep=r, error=f"{type(e).__name__}: {e}"[:200]))
    return rows


# ---- companions (real genes, no clock value) -------------------------------------------------------------------------
def companion_task(job):
    warnings.filterwarnings("ignore")
    kind, i, genes = job
    rows = []
    try:
        for read in ("E", "E2"):
            for rd in ("core", "14 genes"):
                st = C.states(C.as_set(S["wf"][read], genes, rd, f"{kind}{i}"), S["mice"])
                st = [(j, a, b, np.column_stack([cov, np.zeros(len(j))])) for j, a, b, cov in st]
                rows += [dict(kind=kind, set=i, read=read, readout=rd, variant=v, leverage=C.leverage(st, v))
                         for v in VARIANTS]
    except Exception as e:
        rows.append(dict(kind=kind, set=i, error=f"{type(e).__name__}: {e}"[:200]))
    return rows


def run_pool(task, jobs, nproc, t0, every):
    rows = []
    with Pool(nproc, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, jobs, chunksize=1)):
            rows += r
            if (k + 1) % every == 0 or k + 1 == len(jobs):
                print(f"{task.__name__} {k + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    return pd.DataFrame(rows)


# ---- setup -----------------------------------------------------------------------------------------------------------
def setup(n_perm, n_boot):
    pb, per, _ = M.mouse_table()
    pb = pb[pb.compartment.isin(per)].reset_index(drop=True)
    pb["donor"] = pb.donor.astype(str)
    lm = M.log_mean(pb)
    sets = {k: pseudo_waveforms(lm, N_SETS, np.random.default_rng([394, t]), R.WF, kind=k)
            for k, t in (("matched", 2), ("responsive", 3))}
    clock = set(R.WF.genes)
    assert not any(clock & set(w.genes) for ws in sets.values() for w in ws), "a companion set holds a clock gene"
    genes = sorted({g for ws in sets.values() for w in ws for g in w.genes} | clock)
    arm = pb.groupby("donor").stratum.first()
    zt = pd.Series(arm.index.str.extract(r"_ZT(\d+)_", expand=False).astype(float), index=arm.index)
    cell = pd.Series(arm.to_numpy() + "_" + zt.astype(int).astype(str).to_numpy(), index=arm.index)
    fold = pd.Series(0, index=arm.index)
    for ms in cell.groupby(cell).groups.values():
        fold[sorted(ms)[1::2]] = 1
    S.update(mouse=pb[M.KEYS + genes].copy(), per=per, arm=arm, zt=zt, cell=cell, mice=np.asarray(arm.index),
             folds=[np.asarray(fold.index[fold == f]) for f in (0, 1)], halvings=R.halves(per, [400, 1]),
             n_perm=n_perm, n_boot=n_boot)
    S["mouse_full"] = S["mouse"]
    S["wf"] = {"Mf": C.fit_on(S["mice"]), **external_waveforms()}
    S["sim"] = sim_setup(pb)
    return sets


def check_old_companions(sets):
    """The first 50 draws of the 100 are the A400c draws (same seeds)?"""
    old = {k: pseudo_waveforms(M.log_mean(S["mouse_full"]), 50, np.random.default_rng([394, t]), R.WF, kind=k)
           for k, t in (("matched", 2), ("responsive", 3))}
    return all(old[k][i].genes == sets[k][i].genes for k in old for i in range(50))


# ---- gates -----------------------------------------------------------------------------------------------------------
def companion_q95(comp, read, rd="core", v=V2):
    c = comp[(comp.read == read) & (comp.readout == rd) & (comp.variant == v)]
    return {k: float(np.nanquantile(c[c.kind == k].leverage, 0.95)) for k in ("matched", "responsive")}


def mx_companions():
    c = pd.read_csv(OUT / "c_leverages.csv")
    c = c[(c.kind != "clock") & (c.waveform == "Mx") & (c.readout == "core") & (c.variant == V2)]
    return {k: float(np.nanquantile(c[c.kind == k].leverage, 0.95)) for k in ("matched", "responsive")}


def gates(sim, comp):
    q = {"Mx": mx_companions(), "E": companion_q95(comp, "E"), "E2": companion_q95(comp, "E2")}
    rows = []
    for read in READS:
        s = sim[(sim.read == read) & (sim.readout == "core") & (sim.variant == V2)]
        qm, qr = q[read]["matched"], q[read]["responsive"]
        null, pw = s[s.person_sd == 0], s[s.person_sd > 0]
        med = float(null.leverage.median())
        above = (pw.leverage > max(qm, qr)) & (pw.leverage > pw.null_q95)
        rows.append(dict(read=read, matched_q95=qm, responsive_q95=qr, n_null=len(null), median_null_V2_core=med,
                         max_null=float(null.leverage.max()),
                         frac_null_above_both_q95=float((null.leverage > max(qm, qr)).mean()),
                         median_null_perm_q95=float(null.null_q95.median()),
                         G406_1_pass=bool(med < 0.2 and med < min(qm, qr)),
                         null_median_at_or_above_both_q95=bool(med >= max(qm, qr)),
                         n_person=len(pw), median_person_V2_core=float(pw.leverage.median()),
                         frac_person_above_all=float(above.mean()), G406_2_pass=bool(above.mean() >= 0.8)))
    return pd.DataFrame(rows)


# ---- score (blind until the user decides) ----------------------------------------------------------------------------
def lev_rows(st, variant, keep=None, rng=None, boot=None):
    """A400 leverage on a subset: keep = arms kept; rng = permute B within arm; boot = rng resampling mice within cell."""
    L = []
    for j, a, b, cov in st:
        arm, cell = S["arm"].reindex(j).to_numpy(), S["cell"].reindex(j).to_numpy()
        rows = np.flatnonzero(np.isin(arm, keep)) if keep is not None else np.arange(len(j))
        if boot is not None:
            pos = pd.Series(rows).groupby(cell[rows]).apply(lambda x: boot.choice(x.to_numpy(), len(x))).to_numpy()
            rows = np.concatenate(pos)
        j, a, b, cov, arm, cell = j[rows], a[rows], b[rows], cov[rows], arm[rows], cell[rows]
        if rng is not None:
            b = b.copy()
            for p in np.unique(arm):
                m = np.flatnonzero(arm == p)
                b[m] = b[rng.permutation(m)]
        X = np.column_stack([cov, C.harm(j, range(1, 5))]) if variant == VARIANTS[1] else cov
        L.append(R.pooled_leverage(a, b, X, cell if variant == V2 else arm)[0])
    return float(np.nanmedian(L))


def score_task(job):
    warnings.filterwarnings("ignore")
    kind, i, read, rd, genes, arms = job
    try:
        st = [(j, a, b, np.column_stack([cov, np.zeros(len(j))])) for j, a, b, cov in
              C.states(C.as_set(S["wf"][read], genes, rd, f"{kind}{i}"), S["mice"])]
        rows = []
        for arm in arms:
            keep = None if arm == "all" else [arm]
            for v in VARIANTS:
                o = dict(kind=kind, set=i, read=read, readout=rd, arm=arm, variant=v, leverage=lev_rows(st, v, keep))
                if kind == "clock":
                    rng = np.random.default_rng([406, 5])
                    null = np.array([lev_rows(st, v, keep, rng) for _ in range(S["n_perm"])])
                    o.update(p_perm=float((1 + np.sum(null >= o["leverage"])) / (1 + len(null))),
                             null_q95=float(np.quantile(null, 0.95)))
                rows.append(o)
        if kind == "clock" and len(arms) > 1:  # D406: ad libitum minus restricted, mouse bootstrap within cell
            a1, a2 = [x for x in arms if x != "all"][:2]
            rb = np.random.default_rng([406, 6])
            d = [lev_rows(st, V2, [a1], boot=rb) - lev_rows(st, V2, [a2], boot=rb) for _ in range(S["n_boot"])]
            rows.append(dict(kind=kind, set=i, read=read, readout=rd, arm=f"{a1} minus {a2}", variant=V2,
                             leverage=lev_rows(st, V2, [a1]) - lev_rows(st, V2, [a2]),
                             lo90=float(np.quantile(d, 0.05)), hi90=float(np.quantile(d, 0.95))))
        return rows
    except Exception as e:
        return [dict(kind=kind, set=i, read=read, readout=rd, error=f"{type(e).__name__}: {e}"[:200])]


def score(sets, nproc, t0):
    """H406.1 and D406. Reads the real clock genes with E and E2. Run only after the user releases it."""
    g = pd.read_csv(OUT / "d406_gates.csv").set_index("read")
    passed = bool(g.loc["E", "G406_1_pass"] and g.loc["E", "G406_2_pass"])
    comp = pd.read_csv(OUT / "d406_companions.csv")
    arms = ["all"] + sorted(S["arm"].unique())
    jobs = [("clock", 0, rd_read[0], rd_read[1], R.WF.genes, arms) for rd_read in
            [(w, rd) for w in ("E", "E2") for rd in ("core", "14 genes")]]
    # per-arm companions (no clock value) under E, core
    jobs += [(k, i, "E", "core", w.genes, arms[1:]) for k, ws in sets.items() for i, w in enumerate(ws)]
    res = run_pool(score_task, jobs, nproc, t0, 10)
    clk = res[res.kind == "clock"]
    rows = []
    for _, o in clk[clk.arm == "all"].iterrows():
        q = companion_q95(comp, o.read, o.readout, o.variant)
        rows.append(dict(read=o.read, readout=o.readout, variant=o.variant, leverage=o.leverage, p_perm=o.p_perm,
                         null_q95=o.null_q95, matched_q95=q["matched"], responsive_q95=q["responsive"],
                         above_both=bool(o.leverage > max(q.values())),
                         H406_1=bool(o.read == "E" and o.readout == "core" and o.variant == V2 and o.p_perm < 0.05
                                     and o.leverage > max(q.values())),
                         gates_passed=passed))
    arms_out = pd.concat([clk[clk.arm != "all"], res[res.kind != "clock"]])
    return pd.DataFrame(rows), arms_out


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("sim", "score"), default="sim")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    sets = setup(20 if a.check else N_PERM, 20 if a.check else N_BOOT)
    if a.check:
        sets = {k: ws[:2] for k, ws in sets.items()}
    print(f"{len(S['mice'])} mice, {len(S['per'])} tissues, {S['mouse_full'].shape[1] - 4} genes; "
          f"simulation table {S['sim']['d'].shape[0]} samples ({time.time() - t0:.0f}s)", flush=True)
    if a.stage == "score":
        if a.check:
            sys.exit("--check applies to the sim stage")
        sc, arms = score(sets, a.jobs, t0)
        OUT.mkdir(parents=True, exist_ok=True)
        sc.to_csv(OUT / "d406_scores.csv", index=False)
        arms.to_csv(OUT / "d406_arms.csv", index=False)
        print(sc.round(4).to_string(index=False), flush=True)
        return
    same = check_old_companions(sets) if not a.check else None
    print(f"first 50 companion draws equal the A400c draws: {same}", flush=True)
    reps = 1 if a.check else N_REP
    jobs = [(sd, c, r) for c, sd in enumerate(PERSON_SD) for r in range(reps)]
    sim = run_pool(sim_task, jobs, a.jobs, t0, 4)
    OUT.mkdir(parents=True, exist_ok=True)
    if not a.check:
        sim.to_csv(OUT / "d406_sim.csv", index=False)  # kept if the companion stage fails
    S["mouse"] = S["mouse_full"]
    cjobs = [(k, i, list(w.genes)) for k, ws in sets.items() for i, w in enumerate(ws)]
    comp = run_pool(companion_task, cjobs, a.jobs, t0, 8)
    failed = int(sim.error.notna().sum()) if "error" in sim else 0
    failed += int(comp.error.notna().sum()) if "error" in comp else 0
    if a.check:
        print(f"check ok: sim {sim.shape}, companions {comp.shape}, failed {failed}, "
              f"score-path test on a simulated table follows", flush=True)
        S["mouse"] = simulate(0.25, 1, 0)
        st = [(j, a_, b, np.column_stack([cov, np.zeros(len(j))])) for j, a_, b, cov in
              C.states(C.as_set(S["wf"]["E"], R.WF.genes, "core", "chk"), S["mice"])]
        arms = sorted(S["arm"].unique())
        print(f"arms {arms}: lev_rows all/arm/boot ok", round(lev_rows(st, V2), 3),
              round(lev_rows(st, V2, [arms[0]], np.random.default_rng(1)), 3),
              round(lev_rows(st, V2, [arms[0]], boot=np.random.default_rng(2)), 3),
              f"({time.time() - t0:.0f}s)", flush=True)
        return
    G = gates(sim.dropna(subset=["variant"]), comp)
    OUT.mkdir(parents=True, exist_ok=True)
    waveform_table().to_csv(OUT / "d406_waveforms.csv", index=False)
    sim.to_csv(OUT / "d406_sim.csv", index=False)
    comp.to_csv(OUT / "d406_companions.csv", index=False)
    G.to_csv(OUT / "d406_gates.csv", index=False)
    pd.set_option("display.width", 250)
    print(G.round(4).to_string(index=False), f"\nfailed: {failed} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
