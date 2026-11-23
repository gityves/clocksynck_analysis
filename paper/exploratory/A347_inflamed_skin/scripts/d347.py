"""A347 descriptive rows (D347), 14 genes and core loop. Needs the user's go-ahead on real labels (`--go`).

  a  GSE154200 atopic eczema arm: lesional minus non-lesional pooled offset, A183 form (batches with >= 10 atopic pairs and
     a non-lesional reference gate p < 0.05 in `clocksync.paired`, projector per batch built from the atopic persons only,
     weights 1 / se^2 of the batch's 14-gene variance ratio kept for both readouts as A308 / A336), one-sided swap p
     (10,000), MDE = 2.49 SD of the swap null, 90 % patient bootstrap (2,000), matched / narrow responsive / broad
     responsive companions (200 sets each per batch, drawn against the 14-gene waveform, restricted to the core positions
     for the core readout) with the A308 calibrated q (2,000 swap draws).
  b  hidradenitis GSE151243 (A254 `Pair`: one dataset, genes present): same statistics, weight 1; the A254 gate is
     re-run on the waveform read and reported (a failed gate is uninformative).
  c  A345 descriptive: the A342 lead statistic (matched estimator, 20 depth-matched replicates, projector on all rows of the
     disease) within lesional and within non-lesional Tsoi 2019 samples, psoriasis and AD separately; spliced
     `counts_exon`, unspliced `counts_intron_direct`. Covariate: patient indicators fitted on both samples of the
     patients, then the lead is read over the rows of each condition ("patient"; the registration does not say how a
     patient covariate acts inside one condition with one sample per patient); "plain" (no covariate) alongside.
     Gate: kappa against the unspliced rows permuted within condition (10,000); 90 % patient bootstrap (2,000).

    python d347.py --part a|b|c|all [--smoke] [--go]
`--smoke`: labels swapped at random within patient (a, b) or shuffled across samples (c), reduced sizes, outputs to
`../output/d347_smoke/`. Real-label runs write to `../output/d347/` and need `--go`.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A183_skin_gse154200/scripts"))
sys.path.insert(0, str(PAPER / "exploratory/A254_paired_lesions/scripts"))
sys.path.insert(0, str(PAPER / "exploratory/A345_tsoi_premrna/scripts"))
sys.path.insert(0, str(PAPER / "exploratory/A342_skin_timed/scripts"))
import a347 as M  # noqa: E402  (offsets, complex_pairs, boot_offset, wsd, sd0, swap helpers, a336_common as C)
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

C = M.C
FAMILIES = M.FAMILIES
HOURS = P.HOURS
OUT = HERE.parents[1] / "output"
SIZE = dict(n_swap=10000, n_cal=2000, n_sets=200, n_boot=2000, n_perm=10000, reps=20)
SMOKE = dict(n_swap=500, n_cal=200, n_sets=20, n_boot=100, n_perm=500, reps=3)


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def interval(est, draws):
    d = M.wrap(np.asarray(draws) - est) + est
    return float(np.quantile(d, .05)), float(np.quantile(d, .95))


# ---- parts a, b: one paired arm, both readouts ----------------------------------------------------------------------

def run_arm(name, arm, w, tag, sz, out):
    """arm: key -> dict(case, ref, counts, wf, lm); w: key -> weight. One row per readout and family."""
    keys = list(arm)
    sizes = {k: len(arm[k]["case"]) for k in keys}
    F = M.swap_matrix(tag, sizes, sz["n_swap"])
    Fid = {k: np.zeros((1, sizes[k]), bool) for k in keys}
    Fcal = {k: np.vstack([Fid[k], F[k][:sz["n_cal"]]]) for k in keys}
    rows = []
    for ro in ("14", "core"):
        wfk = {k: arm[k]["wf"] if ro == "14" else C.core(arm[k]["wf"]) for k in keys}
        pr = {k: M.complex_pairs(arm[k]["case"], arm[k]["ref"], arm[k]["counts"], wfk[k]) for k in keys}
        o = M.offsets(pr, w, Fid)[0]
        null = M.offsets(pr, w, F)
        sd = M.sd0(null)
        lo, hi = interval(o, M.boot_offset(pr, w, np.random.default_rng([347, tag, 5, int(ro == "core")]), sz["n_boot"]))
        clock = M.offsets(pr, w, Fcal)
        sc = M.wsd(clock[1:])
        base = dict(arm=name, readout=ro, n_series=len(keys), n_pairs=sum(sizes.values()), offset_h=o, lo90=lo, hi90=hi,
                    p_delay=S.perm_p(null, o, "less", tie=0.0), p_two_sided=S.perm_p(null, o, tie=0.0),
                    null_sd_h=sd, mde_h=2.49 * sd, n_swaps=sz["n_swap"], clock_placebo_sd_h=sc)
        for f, fam in enumerate(FAMILIES):
            crng = np.random.default_rng([347, tag, f])
            per = {k: P.pseudo_waveforms(arm[k]["lm"], sz["n_sets"], crng, arm[k]["wf"], fam) for k in keys}
            Mx = np.empty((sz["n_sets"], sz["n_cal"] + 1))
            for j in range(sz["n_sets"]):
                wj = {k: C.restrict(per[k][j], base=arm[k]["wf"]) if ro == "core" else per[k][j] for k in keys}
                Mx[j] = M.offsets({k: M.complex_pairs(arm[k]["case"], arm[k]["ref"], arm[k]["counts"], wj[k]) for k in keys}, w, Fcal)
            sj = np.array([M.wsd(m[1:]) for m in Mx])
            m_, tau2 = Mx.mean(0), np.maximum(0.0, Mx.var(0, ddof=1) - np.mean(sj ** 2))
            qd = norm.cdf((clock - m_) / np.sqrt(tau2 + sc ** 2))
            rows.append(dict(base, family=fam, sets_q05_h=float(np.quantile(Mx[:, 0], .05)), sets_median_h=float(np.median(Mx[:, 0])),
                             beyond_q05=bool(o < np.quantile(Mx[:, 0], .05)), tau2_obs=float(tau2[0]), q=float(qd[0]),
                             v_pass_rate=float(np.mean(qd[1:] < 0.05))))
        say(name, ro, "done")
    T = pd.DataFrame(rows)
    T.to_csv(out / f"{name}.csv", index=False)
    say("\n" + T.round(3).to_string(index=False))


def swap_within_patient(s, seed):
    """Smoke: lesional / non-lesional labels exchanged for a random half of the patients."""
    rng = np.random.default_rng(seed)
    u = s.unit.unique()
    pick = set(rng.choice(u, len(u) // 2, replace=False))
    m = s.unit.isin(pick)
    s = s.copy()
    s.loc[m, "condition"] = s.loc[m, "condition"].map({"lesional": "nonlesional", "nonlesional": "lesional"})
    return s


def part_a(sz, smoke, out):
    import score_a183 as A
    s = A.load()
    if smoke:
        s = swap_within_patient(s, 0)
    means = A.gene_means(s)
    ae = s[s.diagnosis == "Atopic eczema"]
    res = P.paired(ae, unit="unit", condition="condition", levels=("nonlesional", "lesional"), cohort="cohort", counts=True,
                   n_perm=10, n_boot=1000, seed=A.SEED, min_pairs=10)
    used = [k for k in res.cohorts[res.cohorts.gated].cohort]
    say("atopic eczema batches (n pairs, gated):", res.cohorts[["cohort", "n_pairs", "gated"]].to_dict("records"))
    if not used:
        say("no gated batch: part a not scorable")
        return
    se = res.cohorts.set_index("cohort").se
    w = {k: 1 / se[k] ** 2 for k in used}
    arm = {}
    for k in used:
        x = ae[ae.cohort == k]
        case = x[x.condition == "lesional"].set_index("unit")
        ref = x[x.condition == "nonlesional"].set_index("unit")
        u = sorted(set(case.index) & set(ref.index))
        arm[k] = dict(case=case.loc[u], ref=ref.loc[u], counts=True, wf=cs.gtex(), lm=means[k])
    run_arm("gse154200_atopic_eczema", arm, w, 154200, sz, out)


def part_b(sz, smoke, out):
    import compute_a254 as A
    d = A.Pair("GSE151243", placebo=smoke)
    arm = {"GSE151243": dict(case=d.case, ref=d.ref, counts=d.counts, wf=d.wf, lm=d.lm)}
    run_arm("hidradenitis_GSE151243", arm, {"GSE151243": 1.0}, 151243, sz, out)
    rows = []
    for ro, wf in (("14", d.wf), ("core", C.core(d.wf))):
        yr = P.log_expr(d.ref[list(wf.genes)], d.ref.depth.to_numpy(float) if d.counts else None).to_numpy()
        g = P.reference_gate(yr, wf, 200, np.random.default_rng([A.SEED, 1, sum(map(ord, d.key))]))
        rows.append(dict(readout=ro, pairs=d.n, clock_genes=len(wf.genes), gate_p=g["p"], gate_ok=bool(g["p"] < .05)))
    pd.DataFrame(rows).to_csv(out / "hidradenitis_gate.csv", index=False)
    say("\n" + pd.DataFrame(rows).to_string(index=False))


# ---- part c: A342 lead within condition, Tsoi 2019 -----------------------------------------------------------------

def lead_within(s, u, ds, du, patient, cond, wf, covariate, rng, sz):
    """Matched lead per condition: arg sum zu conj(zm) over the rows of the condition (point, kappa, gate p, patient
    bootstrap draws). s, u: n x 14 counts in `a342.GENES` order; ds, du: genome-wide depths."""
    import power_a342 as A
    import score_a342 as S342
    from common import development as dev
    from clocksync.paired import projector
    import power_a328b as PA
    idx = [A.GENES.index(g) for g in wf.genes]
    Z = dev.design(pd.DataFrame(dict(v=patient)), age=None, factors=("v",)) if covariate else np.empty((len(s), 0))
    L = projector(pd.DataFrame(s[:, idx], columns=list(wf.genes)), wf)
    dsc, duc = ds[:, None], du[:, None]
    zu = S342.zfun(u, duc, Z, wf, L)
    zm = np.zeros(len(s), complex)
    for _ in range(sz["reps"]):
        sm, dsm = PA.depth_match(s, u, ds, du, np.zeros(len(s), int), rng)
        zm += S342.zfun(sm, dsm, Z, wf, L)
    zm /= sz["reps"]
    pats = {p: np.flatnonzero(patient == p) for p in np.unique(patient)}
    VB = [rng.choice(list(pats), len(pats)) for _ in range(sz["n_boot"])]
    out = {}
    for c in ("lesional", "non-lesional"):
        i = np.flatnonzero(cond == c)
        a, b = zu[i], zm[i]
        norm_ = np.sqrt(np.sum(abs(a) ** 2) * np.sum(abs(b) ** 2))
        z = np.sum(a * np.conj(b))
        kap = abs(z) / norm_
        perm = rng.permuted(np.tile(a, (sz["n_perm"], 1)), axis=1)
        kp = np.abs((perm * np.conj(b)[None]).sum(1)) / norm_
        draws = []
        for vb in VB:
            r = np.concatenate([pats[p] for p in vb])
            r = r[cond[r] == c]
            draws.append(np.angle(np.sum(zu[r] * np.conj(zm[r]))) * HOURS)
        out[c] = (float(np.angle(z) * HOURS), float(kap), float((1 + np.sum(kp >= kap)) / (1 + len(kp))), np.array(draws))
    return out


def part_c(sz, smoke, out):
    import a345 as A5
    import power_a342 as A
    d = A5.design("tsoi2019_skin_trimmed")
    ex, itd = A5.load("tsoi2019_skin_trimmed", "exon"), A5.load("tsoi2019_skin_trimmed", "intron_direct")
    rng = np.random.default_rng([347, 345])
    rows = []
    for dis in A5.DISEASES:
        dd = d[d.disease == dis].reset_index(drop=True)
        if smoke:
            dd["condition"] = rng.permutation(dd.condition.to_numpy())
        ds, du = ex[dd.run].sum().to_numpy(float), itd[dd.run].sum().to_numpy(float)
        s = ex.loc[A.GENES, dd.run].T.to_numpy(np.int64)
        u = itd.loc[A.GENES, dd.run].T.to_numpy(np.int64)
        pat, cond = dd.patient.to_numpy(), dd.condition.to_numpy()
        for ro, wf in A5.WAVEFORMS.items():
            res = {cov: lead_within(s, u, ds, du, pat, cond, wf, cov == "patient", rng, sz) for cov in ("patient", "plain")}
            for cov, r in res.items():
                for c, (pt, kap, gp, dr) in r.items():
                    lo, hi = interval(pt, dr)
                    rows.append(dict(disease=dis, readout="14" if ro == "all" else ro, covariate=cov, condition=c,
                                     n=int((cond == c).sum()), lead_h=pt, lo90=lo, hi90=hi, kappa=kap, gate_p=gp))
                a, b = r["lesional"], r["non-lesional"]
                dd_ = M.wrap(a[0] - b[0])
                lo, hi = interval(dd_, M.wrap(a[3] - b[3]))
                rows.append(dict(disease=dis, readout="14" if ro == "all" else ro, covariate=cov, condition="lesional_minus_non-lesional",
                                 n=len(pat), lead_h=float(dd_), lo90=lo, hi90=hi))
            say(dis, ro, "done")
    T = pd.DataFrame(rows)
    T.to_csv(out / "tsoi_lead_within_condition.csv", index=False)
    say("\n" + T.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--part", choices=["a", "b", "c", "all"], required=True)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--go", action="store_true", help="real-label run: only with the user's go-ahead")
    a = ap.parse_args()
    if not (a.smoke or a.go):
        sys.exit("real labels need --go (user's go-ahead); use --smoke to test the code path")
    limit_threads(4)
    out = OUT / ("d347_smoke" if a.smoke else "d347")
    out.mkdir(parents=True, exist_ok=True)
    sz = SMOKE if a.smoke else SIZE
    t0 = time.time()
    for p, f in (("a", part_a), ("b", part_b), ("c", part_c)):
        if a.part in (p, "all"):
            f(sz, a.smoke, out)
    say(f"runtime {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
