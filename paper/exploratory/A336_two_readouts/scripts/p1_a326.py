"""A336 part 1 wrapper: A326 H326.1 (design power of A281 / A287 / A321 to read a lymphocyte clock of mouse size; smallest r with
power >= 0.8) on the 14 genes, the core loop, the trio contribution and the gene-drop null. One row per design column:
A326_H326.1:A281, :A287, :A321.

Original: L1_lymphocyte_power/scripts/score_a326.py (imported, not edited; `build_inputs`, `scaled`, `RS`, seeds, `poibin_gt`).
The inserted clock is the 14-gene waveform scaled by r in every readout (the same `scale_clock` draws, same generator
[3260, design, atlas, r, k]); only the READOUT waveform changes, on identical data and identical permutations (the
generator state before the first permutation is copied for every waveform). A281 / A321 readout: signed coherence of the
waveform's genes (`A281.perm_test(wf=...)`); A287 readout: U, pi and p_perm on the waveform's genes with the matched
sets restricted to the positions of those genes (A287 `U`, `prep`, `permuted`, copied with a gene-set argument).
14 genes: replicates.csv (p, passed, pi of every A281 / A287 replicate) and design_power.csv are reproduced to 1e-9 / 1e-6
before the core is read. The A321 column is core in the original: its stored replicates are reproduced by the core readout
and the 14-gene A321 readout is computed here (not stored).

Statistic and criterion (operationalised here, the registered criterion is "power >= 0.8 at r = 0.25"): value = design
power at a reference r chosen from the STORED 14-gene table (A281 0.55, A287 0.70, A321 0.55: where the stored power lies
in (0.1, 0.9)); also reported: power at r = 0.25 and the bound r* (smallest grid r with power >= 0.8). core_status:
holds if core power at r = 0.25 stays < 0.8 (same registered outcome as the 14 genes) and r*_core <= r*_14 + 0.15 (one
grid step), at the line if one of the two holds, falls otherwise. Powered atlases (A281 / A287): the 14 genes use the
stored registered thinned rate >= 0.8 (as the original); core and drops use their own r = 1 pass rate >= 0.8 (20
replicates); A321: the stored informative set for every waveform.
Trio contribution: 200 atlas (A321: compartment) bootstraps with replacement (rng [336, 1]), power at the reference r with the own-rate
powered rule for both waveforms. Gene drops: 100 random + trio, power at the reference r (drops evaluated at r = 1 and the
reference r only; the trio drop equals the core). No companions, no gate (simulation design: trio alone "not
computable").

    python p1_a326.py [--check] [--jobs 2]        (long: paper/common/run.sh)
Output: ../output/part1/a326.csv, a326_replicates.csv (rows: design, key, r, rep, wf, p, passed, pi; core and 14 genes)
"""
import argparse
import copy
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE.parent))
import importlib.util  # noqa: E402

import cloudpickle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a326.py"
spec = importlib.util.spec_from_file_location("p1_a326_orig", PAPER / "exploratory/L1_lymphocyte_power/scripts/score_a326.py")
S = importlib.util.module_from_spec(spec)
sys.modules["p1_a326_orig"] = S
spec.loader.exec_module(S)
cloudpickle.register_pickle_by_value(S)
A281, A287, A321 = S.A281, S.A287, S.A321
BASE = C.WF
CORE, TRIO = C.core(BASE), C.trio(BASE)
RS, SEED, N_REP, N_PERM_THIN = S.RS, S.SEED, S.N_REP, S.N_PERM_THIN
R_REF = {"A281": 0.55, "A287": 0.70, "A321": 0.55}
STORED = PAPER / "exploratory/L1_lymphocyte_power/output"
CLOCK = list(BASE.genes)


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def restrict_sets(matched, w):
    keep = [k for k, g in enumerate(BASE.genes) if g in w.genes]
    return [[s[k] for k in keep] for s in matched]


def readout287(t, sets, n_perm, rng, clock):
    """A287 `readout` on the genes `clock` (U, pi, p_perm, structure)."""
    t = t.reset_index(drop=True)
    y = A287.prep(t, clock)
    u = A287.U(y)
    comp = np.array([A287.U(A287.prep(t, s)) for s in sets])
    groups = A281.groups_of(t)
    null = np.array([A287.U(A287.permuted(y, groups, rng)) for _ in range(n_perm)])
    pi = float(np.mean(comp[np.isfinite(comp)] < u))
    p = float((1 + np.sum(null >= u)) / (1 + n_perm))
    return dict(pi=pi, p_perm=p, structure=bool(pi >= 0.95 and p < 0.05))


def rep_task(design, key, r, reps, di, ai, ri, ts, target, extra, wfs):
    """Replicates `reps` of (design, key, r): rows for every waveform in `wfs` [(label, wf)] on the same data and permutations."""
    warnings.filterwarnings("ignore")
    out = []
    for k in reps:
        rng = np.random.default_rng([SEED, di, ai, ri, k])
        tt = S.scaled(ts, r, rng)
        if design == "A287":
            genes = extra["genes"]
            s = tt[CLOCK].sum(1).to_numpy(float)
            pr = np.minimum(1.0, target / np.maximum(s, 1))
            tt[genes] = rng.binomial(tt[genes].to_numpy(np.int64), pr[:, None])
            tt["depth"] = tt.depth.to_numpy(float) * pr
        else:
            tt = A281.thin(tt, target, rng)
        state = copy.deepcopy(rng)
        for lab, w in wfs:
            g = copy.deepcopy(state)
            if design == "A287":
                if lab == "14":
                    res, _ = A287.readout(tt, extra["matched"], N_PERM_THIN, g)
                else:
                    res = readout287(tt, restrict_sets(extra["matched"], w), N_PERM_THIN, g, list(w.genes))
                out.append(dict(design=design, key=key, r=r, rep=k, wf=lab, p=res["p_perm"], passed=bool(res["structure"]), pi=res["pi"]))
            else:
                p = A281.perm_test(A281.y_of(tt, list(w.genes)), A281.groups_of(tt), N_PERM_THIN, g, wf=w)[1]
                out.append(dict(design=design, key=key, r=r, rep=k, wf=lab, p=p, passed=p < 0.05, pi=np.nan))
    return out


def poibin_power(rates):
    return S.poibin_gt(np.asarray(rates, float), 2)


def power_of(df, design, keys, r, powered, wf):
    """Design power of waveform `wf` at r over `keys` (list, repeats allowed) with the `powered` keys of A281 / A287."""
    d = df[(df.design == design) & (df.wf == wf) & (df.r == r)]
    if design == "A321":
        P = {k: g.sort_values("rep").p.to_numpy() for k, g in d.groupby("key")}
        z = sum(norm.isf(P[k]) for k in keys) / np.sqrt(len(keys))
        return float(np.mean(norm.sf(z) < 0.05))
    rate = d.groupby("key").passed.mean()
    pk = [k for k in keys if k in powered]
    if design == "A281":
        return poibin_power([rate[k] for k in pk]) if pk else np.nan
    return float(1 - np.prod([1 - rate[k] for k in pk])) if pk else np.nan


def own_powered(df, design, wf):
    d = df[(df.design == design) & (df.wf == wf) & (df.r == 1.0)]
    rate = d.groupby("key").passed.mean()
    return set(rate.index[rate >= 0.8])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    nrep = 4 if chk else N_REP
    specs = S.build_inputs()
    if chk:
        keep = {}
        for s in specs:
            keep.setdefault(s["design"], []).append(s)
        specs = [x for v in keep.values() for x in v[:1]]
    dl, tl = C.drops(BASE, n=100)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:1]
    main_wfs = [("14", BASE), ("core", CORE)]
    drop_wfs = [(f"drop:{lab}", w) for lab, w in dl]
    rs_all = [1.0, 0.55, 0.25] if chk else RS
    tasks = []
    for s in specs:
        di = S.DESIGNS.index(s["design"])
        for r in rs_all:
            ri = RS.index(r)
            wfs = main_wfs + (drop_wfs if r in (1.0, R_REF[s["design"]]) else [])
            for c in range(0, nrep, S.CHUNK):
                tasks.append(delayed(rep_task)(s["design"], s["key"], r, range(c, min(c + S.CHUNK, nrep)), di, s["ai"], ri, s["ts"],
                                               s["target"], s["extra"], wfs))
    log(f"{len(tasks)} tasks, {a.jobs} workers")
    res, t0 = [], time.time()
    for i, x in enumerate(Parallel(n_jobs=a.jobs, return_as="generator_unordered")(tasks), 1):
        res += x
        if i % 10 == 0 or i == len(tasks):
            log(f"{i}/{len(tasks)} tasks, {(time.time() - t0) / 60:.1f} min")
    df = pd.DataFrame(res).sort_values(["design", "key", "r", "rep", "wf"]).reset_index(drop=True)
    reg = {(s["design"], s["key"]): s["registered"] for s in specs}
    keys = {d: [s["key"] for s in specs if s["design"] == d] for d in S.DESIGNS}
    info321 = {s["key"] for s in specs if s["design"] == "A321" and s["extra"]["informative"]} or (set(keys["A321"]) if chk else set())
    # ---- 14 genes (and the A321 core column) reproduced first --------------------------------------------------------
    stored_rep = pd.read_csv(STORED / "replicates.csv")
    stored_dp = pd.read_csv(STORED / "design_power.csv")
    bad = []
    if not chk:
        for design, wf in (("A281", "14"), ("A287", "14"), ("A321", "core")):
            m = df[(df.design == design) & (df.wf == wf)].merge(stored_rep, on=["design", "key", "r", "rep"], suffixes=("", "_s"))
            n_exp = int(((stored_rep.design == design)).sum())
            if len(m) != n_exp or not np.allclose(m.p, m.p_s, atol=1e-9, rtol=0) or (m.passed != m.passed_s).any():
                bad.append((design, "replicates", len(m), n_exp))
        pw14 = {}
        for r in RS:
            powered = {"A281": {k for k in keys["A281"] if reg[("A281", k)] >= 0.8}, "A287": {k for k in keys["A287"] if reg[("A287", k)] >= 0.8}}
            for design, wf in (("A281", "14"), ("A287", "14"), ("A321", "core")):
                v = power_of(df, design, sorted(info321) if design == "A321" else keys[design], r, powered.get(design, set()), wf)
                want = stored_dp[(stored_dp.design == design) & np.isclose(stored_dp.r, r)].design_power.iloc[0]
                if not np.isclose(v, want, atol=1e-6, rtol=0):
                    bad.append((design, "design_power", r, v, want))
    log("14-gene reproduction:", "OK" if not bad else bad)
    if bad:
        C.write_rows("a326", [dict(id=f"A326_H326.1:{d}", figure_panel="Fig 2b", reproduced=False, notes=f"not reproduced: {bad}", script=SCRIPT)
                              for d in S.DESIGNS])
        raise SystemExit("A326 not reproduced")
    # ---- power curves -----------------------------------------------------------------------------------------------
    def curve(design, wf, powered):
        kk = sorted(info321) if design == "A321" else keys[design]
        return {r: power_of(df, design, kk, r, powered, wf) for r in rs_all}

    def rstar(c):
        ok = [r for r in sorted(c) if np.isfinite(c[r]) and c[r] >= 0.8]
        return min(ok) if ok else np.inf

    reg_pow = {d: {k for k in keys[d] if reg[(d, k)] >= 0.8} for d in ("A281", "A287")}
    cur = {}
    for design in S.DESIGNS:
        if design == "A321":
            cur[design] = {"14": curve(design, "14", set()), "core": curve(design, "core", set())}
        else:
            cur[design] = {"14": curve(design, "14", reg_pow[design]), "core": curve(design, "core", own_powered(df, design, "core"))}
    # ---- trio contribution: atlas / compartment bootstrap ---------------------------------------------------------------
    rng = np.random.default_rng([336, 1])
    nb = 5 if chk else 200
    boot = {}
    for design in S.DESIGNS:
        kk = keys[design]
        rr = R_REF[design] if R_REF[design] in rs_all else rs_all[1]
        p14 = {"14": own_powered(df, design, "14"), "core": own_powered(df, design, "core")}
        d_ = []
        for _ in range(nb):
            b = [kk[i] for i in rng.integers(0, len(kk), len(kk))]
            if design == "A321":
                b = [k for k in b if k in info321] or sorted(info321)
            d_.append(power_of(df, design, b, rr, p14["14"], "14") - power_of(df, design, b, rr, p14["core"], "core"))
        boot[design] = np.array(d_)
    # ---- gene drops --------------------------------------------------------------------------------------------------
    dvals = {d: {} for d in S.DESIGNS}
    for design in S.DESIGNS:
        rr = R_REF[design] if R_REF[design] in rs_all else rs_all[1]
        kk = sorted(info321) if design == "A321" else keys[design]
        for lab, _ in dl:
            wf = f"drop:{lab}"
            pw = own_powered(df, design, wf) if design != "A321" else set()
            dvals[design][lab] = power_of(df, design, kk, rr, pw, wf)
        assert np.isclose(dvals[design][tl], cur[design]["core"][rr], equal_nan=True), (design, dvals[design][tl], cur[design]["core"][rr])
    if chk:
        print("check ok")
        return
    out = HERE.parents[1] / "output" / "part1"
    out.mkdir(parents=True, exist_ok=True)
    df[~df.wf.str.startswith("drop:")].to_csv(out / "a326_replicates.csv", index=False)
    dwf = dict(dl)
    rows = []
    for design in S.DESIGNS:
        rr = R_REF[design]
        c14, cc = cur[design]["14"], cur[design]["core"]
        stat_stored = float(stored_dp[(stored_dp.design == design) & np.isclose(stored_dp.r, rr)].design_power.iloc[0])
        pr = C.drop_percentiles(dvals[design], tl, dwf, direction=1)
        conds = [bool(cc[0.25] < 0.8), bool(rstar(cc) <= rstar(c14) + 0.15 + 1e-9)]
        status = "holds" if all(conds) else "at the line" if any(conds) else "falls"
        d_ = boot[design]
        tab = "; ".join(f"r={r:g}: 14 {c14[r]:.3f} / core {cc[r]:.3f}" for r in RS)
        rows.append(dict(
            id=f"A326_H326.1:{design}", figure_panel="Fig 2b", statistic=f"design power of {design} to read a lymphocyte clock at r = {rr:g} of fibroblast amplitude",
            direction=1, value_14_stored=(stat_stored if design != "A321" else np.nan), value_14=c14[rr], reproduced=True, value_core=cc[rr],
            core_criterion="power at r = 0.25 stays < 0.8 (same outcome as the 14 genes) and bound r* <= 14-gene r* + 0.15 (operationalised; registered: power >= 0.8 at r = 0.25)",
            core_status=status, diff_14_minus_core=c14[rr] - cc[rr], diff_lo90=float(np.nanquantile(d_, .05)), diff_hi90=float(np.nanquantile(d_, .95)),
            diff_draws="200 atlas / compartment bootstraps (rng [336, 1]), own-rate powered rule for both waveforms", **pr,
            trio_gate_passes="not computable", reading=C.reading(status, pr["trio_drop_below_q10"]),
            blind="core not blind (A329 / A335)" if design == "A321" else "core blind", script=SCRIPT,
            notes=(f"power curve ({tab}); r*: 14 genes {rstar(c14):g}, core {rstar(cc):g}; power at 0.25: 14 {c14[0.25]:.3f}, core {cc[0.25]:.3f}; "
                   f"powered atlases core (own r=1 rate >= 0.8): {len(own_powered(df, design, 'core')) if design != 'A321' else len(info321)} "
                   f"vs 14 genes {len(reg_pow[design]) if design != 'A321' else len(info321)}"
                   + ("; A321 column: the core readout reproduces the stored replicates, the 14-gene value is computed here (not stored)" if design == "A321" else "")
                   + "; drops evaluated at r = 1 and the reference r only")))
    d = C.write_rows("a326", rows)
    pd.set_option("display.width", 250)
    print(d.drop(columns=["script", "notes", "core_criterion", "diff_draws", "figure_panel", "statistic"]).round(4).T.to_string())
    for r in rows:
        print(r["notes"])


if __name__ == "__main__":
    main()
