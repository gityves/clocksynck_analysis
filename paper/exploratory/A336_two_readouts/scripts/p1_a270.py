"""A336 part 1, A270 (Figure 1c): H270.1 (baboon), H270.2 (GTEx fast deaths) and the clock R against recorded time.

Original: `A270_responsive_time/scripts/score_a270.py` (imported, not edited: `load`, `log_mean`, SEED 270). Same data, seeds
([270, 1] responsive, [270, 2] matched, 200 sets each, `Dataset.timed(hours, n_perm=1000, seed=270)`); only the waveform
changes. Companions: the original 14-gene draws restricted to the core positions. GTEx: aggregates only (hours stay in
memory, nothing per donor is printed or written).

Rows: A270_H270.1 and A270_H270.2 (statistic = share of responsive sets with p_R < 0.05; criterion share >= 0.5, the
registered verdict; the human 14-gene verdict is FAIL (0.17), so for H270.2 "holds" = the core share stays below 0.5 with
the 90 % upper of its set bootstrap below 0.5, crossing = "falls" (not decidable unless the trio drop is below q10);
for H270.1 "holds" = core share >= 0.5); A270_clockR_baboon / _human (descriptive: clock R above the q95 of the matched and
of the responsive companions; both: holds; one: at the line).
Trio contribution: H270.x share: 200 bootstraps of the pseudo-clock SETS (the unit of a companion-share statistic; a donor
bootstrap would need 400 refits per draw), same resamples for both waveforms; clock R: 200 donor (animal) bootstraps
(copies renamed, host refitted, `timed(n_perm=1)`), shared. Gene-drop null: the clock R (no companions) on every
three-gene drop (364); for H270.x the clock statistic is the clock R (the share is a companion statistic).
Trio-alone gate: `Dataset.gate(n_perm=1000, seed=336)`, p_pairing and p_acrophase < 0.05 (all tissues as anchors).

    python p1_a270.py [--check] [--jobs 4] [--species baboon human]
Output: ../output/part1/a270.csv
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A270_responsive_time/scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import score_a270 as A  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

ENTRY = "a270"
BASE = C.WF
META = A.META
SEED = A.SEED
OUTA = PAPER / "exploratory/A270_responsive_time/output"
G = {}
SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a270.py"
PANEL = "Fig 1c"


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def dataset(pb, w, name="x"):
    return cs.Dataset(pb[META + list(w.genes)], anchors=sorted(pb.compartment.unique()), waveform=w, name=name)


def task(i):
    """(R, p_R) of job i = (waveform, n_perm); NaN on failure."""
    warnings.filterwarnings("ignore")
    w, n_perm = G["jobs"][i]
    try:
        r = dataset(G["pb"], w).timed(G["hours"], n_perm=n_perm, seed=SEED)
        return i, (float(r.R), float(r.p_R))
    except Exception:  # noqa: BLE001
        return i, (np.nan, np.nan)


def run(jobs_list, n_jobs):
    G["jobs"] = jobs_list
    with Pool(n_jobs, initializer=limit_threads, initargs=(1,)) as p:
        out = dict(p.imap_unordered(task, range(len(jobs_list))))
    return [out[i] for i in range(len(jobs_list))]


def boot_task(args):
    """Donor bootstrap draw: (R of the 14-gene clock, R of the core clock) on the same resampled donors."""
    warnings.filterwarnings("ignore")
    b, draw = args
    pb, hours = G["pb"], G["hours"]
    by = G["by_donor"]
    parts, hh = [], {}
    for j, d in enumerate(draw):
        x = pb.iloc[by[d]].copy()
        x["donor"] = f"{d}@{j}"
        parts.append(x)
        hh[f"{d}@{j}"] = hours.loc[d]
    rep = pd.concat(parts, ignore_index=True)
    hs = pd.Series(hh)
    out = []
    for w in (BASE, G["core"]):
        try:
            out.append(float(dataset(rep, w).timed(hs, n_perm=1, seed=SEED).R))
        except Exception:  # noqa: BLE001
            out.append(np.nan)
    return b, out


def species_rows(species, jobs, chk):
    n = 3 if chk else 200
    t0 = time.time()
    pb, hours = A.load(species)
    limit_threads(1)
    lm = A.log_mean(pb)
    core, trio = C.core(BASE), C.trio(BASE)
    sets = {k: pseudo_waveforms(lm, n, np.random.default_rng([SEED, j]), BASE, k) for j, k in ((1, "responsive"), (2, "matched"))}
    used = sorted({g for ws in sets.values() for w in ws for g in w.genes} | set(BASE.genes))
    pbu = pb[META + used].copy()
    del pb
    G.update(pb=pbu, hours=hours, core=core)
    stored = pd.read_csv(OUTA / f"tests_{species}.csv").iloc[0]
    # ---- 14 genes: clock and companions, reproduce -----------------------------------------------------------------
    r = run([(BASE, 1000)] + [(w, 1000) for k in ("responsive", "matched") for w in sets[k]], jobs)
    R14, comp14 = r[0][0], {"responsive": r[1:1 + n], "matched": r[1 + n:1 + 2 * n]}
    share14 = float(np.mean([p < 0.05 for R, p in comp14["responsive"] if np.isfinite(R)]))
    q14 = {k: float(np.nanquantile([x[0] for x in v], 0.95)) for k, v in comp14.items()}
    rep = {"R": np.isclose(R14, stored.R_clock, atol=1e-6), "share": np.isclose(share14, stored.responsive_share_p05, atol=1e-6),
           "resp_q95": np.isclose(q14["responsive"], stored.responsive_R_q95, atol=1e-6),
           "matched_q95": np.isclose(q14["matched"], stored.matched_R_q95, atol=1e-6)}
    if not chk:
        log(species, "reproduction", {k: bool(v) for k, v in rep.items()})
        if not all(rep.values()):
            return [dict(id=i, figure_panel=PANEL, statistic=s, direction=1, value_14_stored=v, value_14=v14, reproduced=False,
                         notes=f"{species}: not reproduced {rep}", script=SCRIPT)
                    for i, s, v, v14 in ((f"A270_H270.{1 if species == 'baboon' else 2}", "share of responsive sets with p_R < 0.05",
                                          stored.responsive_share_p05, share14),
                                         (f"A270_clockR_{species}", "clock R", stored.R_clock, R14))]
    log(species, "14-gene stage done", f"{(time.time() - t0) / 60:.1f} min")
    # ---- core ------------------------------------------------------------------------------------------------------
    dl, tl = C.drops(BASE)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:2]
    cw = {k: [C.restrict(w, BASE, core) for w in sets[k]] for k in sets}
    r = run([(core, 1000)] + [(w, 1000) for k in ("responsive", "matched") for w in cw[k]] + [(w, 1) for _, w in dl], jobs)
    Rc = r[0][0]
    compc = {"responsive": r[1:1 + n], "matched": r[1 + n:1 + 2 * n]}
    dv = {lab: x[0] for (lab, _), x in zip(dl, r[1 + 2 * n:])}
    assert np.isclose(dv[tl], Rc, equal_nan=True)
    dwf = dict(dl)
    okr = np.array([np.isfinite(R) for R, _ in compc["responsive"]])
    pr_c = np.array([p for _, p in compc["responsive"]])
    pr_14 = np.array([p for _, p in comp14["responsive"]])
    ok14 = np.array([np.isfinite(R) for R, _ in comp14["responsive"]])
    sharec = float(np.mean(pr_c[okr] < 0.05))
    qc = {k: float(np.nanquantile([x[0] for x in v], 0.95)) for k, v in compc.items()}
    log(species, "core and drops done", f"{(time.time() - t0) / 60:.1f} min")
    # ---- bootstraps ------------------------------------------------------------------------------------------------
    rb = np.random.default_rng([336, 1])
    nb = 5 if chk else 200
    sidx = [rb.integers(0, n, n) for _ in range(nb)]
    ds14 = np.array([np.mean((pr_14[i][ok14[i]] < 0.05)) for i in sidx])
    dsc = np.array([np.mean((pr_c[i][okr[i]] < 0.05)) for i in sidx])
    donors = sorted(hours.index.intersection(G["pb"].donor.unique()))
    by = G["pb"].groupby("donor").indices
    G["by_donor"] = {d: by[d] for d in donors}
    draws = [(b, [donors[k] for k in rb.integers(0, len(donors), len(donors))]) for b in range(nb)]
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        bres = dict(p.imap_unordered(boot_task, draws))
    bR = np.array([bres[b] for b in range(nb)])
    # ---- trio alone ------------------------------------------------------------------------------------------------
    pbt = G["pb"]
    g = dataset(pbt, trio, "trio").gate(n_perm=50 if chk else 1000, seed=336)
    gate = bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05)
    Rt = float(dataset(pbt, trio).timed(hours, n_perm=1, seed=SEED).R) if gate else np.nan
    if chk:
        return []
    return dict(species=species, R14=R14, Rc=Rc, share14=share14, sharec=sharec, q14=q14, qc=qc, stored=stored, dv=dv, tl=tl,
                dwf=dwf, ds=(ds14 - dsc), dsc=dsc, bR=bR, gate=gate, g=g, Rt=Rt, comp14=comp14, compc=compc, n=n,
                med=dict(resp=float(np.nanmedian([x[0] for x in compc["responsive"]])), mat=float(np.nanmedian([x[0] for x in compc["matched"]]))))


def rows_of(x):
    sp, s = x["species"], x["stored"]
    out = []
    # share row
    hid = "A270_H270.1" if sp == "baboon" else "A270_H270.2"
    lo_s, hi_s = np.quantile(x["dsc"], [0.05, 0.95])
    pr = C.drop_percentiles(x["dv"], x["tl"], x["dwf"], direction=1)
    if sp == "baboon":
        status = "holds" if x["sharec"] >= 0.5 else "falls"
        crit = "share of core responsive sets with p_R < 0.05 >= 0.5 (responsive family carries time)"
    else:
        status = "holds" if (x["sharec"] < 0.5 and hi_s < 0.5) else "falls"
        crit = "registered FAIL (share < 0.5): core share < 0.5 with 90 % set-bootstrap upper < 0.5"
    out.append(dict(
        id=hid, figure_panel=PANEL, statistic=f"share of 200 responsive pseudo-clock sets with p_R < 0.05 against time ({sp})", direction=1,
        value_14_stored=s.responsive_share_p05, value_14=x["share14"], reproduced=True, value_core=x["sharec"], core_criterion=crit,
        core_status=status, diff_14_minus_core=x["share14"] - x["sharec"], diff_lo90=float(np.quantile(x["ds"], .05)),
        diff_hi90=float(np.quantile(x["ds"], .95)), diff_draws="200 bootstraps of the 200 pseudo-clock sets (rng [336, 1]), shared",
        **pr, trio_gate_passes=x["gate"], trio_alone_value=x["Rt"] if x["gate"] else np.nan,
        reading=C.reading(status, pr["trio_drop_below_q10"]), blind="core blind",
        notes=(f"drop null and trio-alone value are of the clock R (the share is a companion statistic; trio-alone value = trio clock R); "
               f"core share 90 % set-bootstrap {lo_s:.3f} to {hi_s:.3f}; core responsive R median {x['med']['resp']:.3f} q95 {x['qc']['responsive']:.3f}; "
               f"trio gate p_pairing {x['g'].p_pairing:.4f}, p_acrophase {x['g'].p_acrophase:.4f}"), script=SCRIPT))
    # clock R row
    d = x["bR"][:, 0] - x["bR"][:, 1]
    conds = [bool(x["Rc"] > x["qc"]["matched"]), bool(x["Rc"] > x["qc"]["responsive"])]
    st = "holds" if all(conds) else "at the line" if any(conds) else "falls"
    out.append(dict(
        id=f"A270_clockR_{sp}", figure_panel=PANEL, statistic=f"clock R with recorded time ({sp})", direction=1,
        value_14_stored=s.R_clock, value_14=x["R14"], reproduced=True, value_core=x["Rc"],
        core_criterion="descriptive: core R above q95 of matched and of responsive companions (restricted to core)",
        core_status=st, diff_14_minus_core=x["R14"] - x["Rc"], diff_lo90=float(np.nanquantile(d, .05)), diff_hi90=float(np.nanquantile(d, .95)),
        diff_draws="200 donor bootstraps (rng [336, 1]), shared", **pr, trio_gate_passes=x["gate"],
        trio_alone_value=x["Rt"] if x["gate"] else np.nan, reading=C.reading(st, pr["trio_drop_below_q10"]), blind="core blind",
        notes=(f"core companion q95: matched {x['qc']['matched']:.3f}, responsive {x['qc']['responsive']:.3f} (14 genes {x['q14']['matched']:.3f} / {x['q14']['responsive']:.3f}); "
               f"core R 90 % donor bootstrap {np.nanquantile(x['bR'][:, 1], .05):.3f} to {np.nanquantile(x['bR'][:, 1], .95):.3f}; "
               f"trio gate p_pairing {x['g'].p_pairing:.4f}, p_acrophase {x['g'].p_acrophase:.4f}"), script=SCRIPT))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--species", nargs="*", default=["baboon", "human"])
    a = ap.parse_args()
    rows = []
    for sp in a.species:
        x = species_rows(sp, a.jobs, a.check)
        if a.check:
            continue
        rows += x if isinstance(x, list) else rows_of(x)
        C.write_rows(ENTRY, rows)
    if a.check:
        print("check ok")
        return
    df = C.write_rows(ENTRY, rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws"]).round(4).T.to_string())
    for r in rows:
        print(r["id"], r.get("notes"))


if __name__ == "__main__":
    main()
