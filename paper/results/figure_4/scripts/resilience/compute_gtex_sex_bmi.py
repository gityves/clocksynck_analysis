"""A99 (sex), A96 part 1 (BMI, T2D, liver disease) and the A94 cue-versus-core secondary: GTEx individual amplitude.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A94-A99 common definitions; A94 secondary; A96 part 1; A99). Same per-sample
clock states, radius, wall-time alignment and inter-organ coherence as compute_gtex_age.py (imported: one definition).
  H99 (sex): radius (phase-free individual amplitude), female minus male at equal n per tissue, median over the 14 body tissues
    (none is sex-specific), fast and ventilated cohorts each then their mean; sex permuted within Hardy class x age tertile
    (n_perm); two-sided. Descriptive: R, C, brain-body C by sex; timed log amplitude and phase by sex (fast deaths).
  H96.1 (liver, BMI): timed log amplitude per SD of BMI < 0 (fast deaths; complex regression with a BMI x e^{iwt} term;
    covariates sex, Hardy, ischaemic time, age); BMI permuted within Hardy class (2,000 in the full run); one-sided.
    Secondary: the same in the other metabolic tissues and two companions (tibial artery, nerve); the radius slope on BMI and the
    radius ratio of T2D and recorded liver disease (>= 10 cases), fast and ventilated.
  A94 secondary, cue vs core (fast deaths): timed amplitude of the A81b cue state relative to the core state, oldest minus
    youngest age tertile, and the cue-minus-core phase offset (h); median over body tissues. Cue waveform: the A81b cue
    reference of figure_S07 (`a81b_cue_reference.csv`, run compute_cue_reference.py first).
Null: label permutations as named; estimand: differences of the readouts above. Controlled access: exact age, BMI, T2D, liver
disease, Hardy class and time of death are held in memory, only aggregates are written.

    python compute_gtex_sex_bmi.py [--subsample N] [--seed S]      (--subsample N: N random donors; 50 permutations)
Outputs (paper/results/figure_4/output/resilience/): a99_a96_summary.csv, a99_groups.csv, a99_timed_sex.csv,
a96_gtex_metabolic.csv, a94_cue_core.csv
Legacy: legacy/desync/validation/gtex_covariates/scripts/a94_age_sex_bmi.py (A99, A96 part 1, cue-core), gtex_states.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import projector as _proj, states as _states  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

import compute_gtex_age as A  # noqa: E402

FOLDER = "resilience"
W, BODY, WF = A.W, A.BODY, A.WF
METAB = ["Liver", "Pancreas", "Adipose - Visceral (Omentum)", "Adipose - Subcutaneous", "Muscle - Skeletal"]
COMPANION = ["Artery - Tibial", "Nerve - Tibial"]
CUE_REF = "a81b_cue_reference.csv"


def cue_waveform():
    """The A81b cue set as a clocksync waveform (local: the file of figure_S07 -> Waveform, see its README.txt)."""
    ref = pd.read_csv(out_dir("gtex_cues", False) / CUE_REF, index_col=0)
    x = ref[ref.set == "cue"]
    return cs.waveform.Waveform(tuple(x.index), x.amp.to_numpy(), x.acro_h.to_numpy() * W, None, None, "cue")


def sample_ischaemic_h(x: pd.DataFrame) -> np.ndarray:
    """Ischaemic time (h) of the sample behind every (donor, tissue) row (local: `gtex.load` keeps one value per donor; the
    legacy covariate was the tissue's own sample)."""
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs[["subject", "tissue", "ischemic_time"]].copy()
    finally:
        a.file.close()
    o = o.assign(donor=o.subject.astype(str), compartment=o.tissue.astype(str)).drop_duplicates(["donor", "compartment"])
    o = o.set_index(["donor", "compartment"]).ischemic_time.astype(float) / 60
    return o.reindex(pd.MultiIndex.from_frame(x[["donor", "compartment"]])).to_numpy()


def load_states(tissues, subsample, seed, cue) -> pd.DataFrame:
    """compute_gtex_age.load_states plus the protected BMI, T2D and liver disease (0/1), and the cue state (bc, vc)."""
    genes = list(WF.genes) + (list(cue.genes) if cue is not None else [])
    pb, don = gtex.load(tissues=tissues, hardy=(0, 1, 2, 3, 4), genes=genes, depth="library_size")
    ph = gtex.phenotypes(["AGE", "BMI", "MHT2D", "MHLVRDIS"])
    don = don.join(ph.AGE.astype(float).rename("age_exact")).join(pd.to_numeric(ph.BMI, errors="coerce").rename("bmi"))
    for c, n in (("MHT2D", "t2d"), ("MHLVRDIS", "liver_dis")):
        don = don.join(pd.to_numeric(ph[c], errors="coerce").where(lambda v: v.isin([0, 1])).rename(n))
    if subsample:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb = pb[pb.donor.isin(keep)]
    rows = []
    for t, x in pb.groupby("compartment"):
        if len(x) < 20:
            continue
        x = x.reset_index(drop=True)
        depth = x.depth.to_numpy(float)
        b, v = _states(x, _proj(x, WF, counts=True), WF, depth=depth)
        r = dict(tissue=t, subject=x.donor.to_numpy(), b=b[:, 0] + 1j * b[:, 1], v=v)
        if cue is not None:
            bc, vc = _states(x, _proj(x, cue, counts=True), cue, depth=depth)
            r.update(bc=bc[:, 0] + 1j * bc[:, 1], vc=vc)
        d = don.reindex(x.donor)
        r.update(tod=d.tod.to_numpy(float), hardy=d.hardy.to_numpy(int), sex=1 - d.male.to_numpy(float),
                 isch=sample_ischaemic_h(x), age=d.age_exact.to_numpy(float), bmi=d.bmi.to_numpy(float),
                 t2d=d.t2d.to_numpy(float), liver_dis=d.liver_dis.to_numpy(float))
        rows.append(pd.DataFrame(r))
    s = pd.concat(rows, ignore_index=True)
    s["isch"] = s.isch.fillna(s.isch.median())
    return s


def timed_slope(x, col, extra_age=False):
    """Log-amplitude and phase (h) change per unit of `col` with SEs (as compute_gtex_age.timed_slope; optional age)."""
    Z = A.covariates(x)
    if extra_age:
        Z = np.column_stack([Z, (x.age - x.age.mean()) / x.age.std()])
    f = A.timed_fit(x.b.to_numpy(), x.tod.to_numpy(), Z, x[col].to_numpy(float))
    q = f["a1"] / f["a"]
    se = f["se_a1"] / abs(f["a"])
    return dict(log_amp=q.real, phase_h=q.imag / W, se_log_amp=se, se_phase_h=se / W, n=len(x))


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm, n_bmi = (50, 100) if run["quick"] else (1000, 2000)
    rng = np.random.default_rng(94 + args.seed)
    seed = 94
    brain = A._brain_tissues()
    s = load_states(BODY + ["Pancreas"] + brain, args.subsample, args.seed, cue_waveform()).dropna(subset=["tod"])
    fast_s, vent_s = s[s.hardy.isin([1, 2])], s[s.hardy == 0]
    C = {"fast": A.Cohort(fast_s, True), "vent": A.Cohort(vent_s, False)}
    brain = sorted(t for t in s.tissue.unique() if t.startswith("Brain"))
    res, rows, labs = {}, [], {}
    for k, c in C.items():
        labs[k] = A.tertiles(c.donors.dropna(subset=["age"]).age, rng)

    # ---- A94 secondary: cue vs core by age tertile (fast deaths)
    cv, lab = [], labs["fast"]
    for t, x in fast_s.groupby("tissue"):
        if t not in BODY:
            continue
        r = {}
        for g in (0, 2):
            y = x[x.subject.isin(lab.index[lab == g])]
            if len(y) < A.MIN_N:
                break
            Z = A.covariates(y)
            fc = A.timed_fit(y.b.to_numpy(), y.tod.to_numpy(), Z)["a"]
            fu = A.timed_fit(y.bc.to_numpy(), y.tod.to_numpy(), Z)["a"]
            r[g] = (np.log(abs(fu) / abs(fc)), np.angle(fu / fc) / W)
        if len(r) == 2:
            cv.append(dict(tissue=t, d_log_cue_over_core=r[2][0] - r[0][0],
                           d_cue_minus_core_h=(r[2][1] - r[0][1] + 12) % 24 - 12))
    cv = pd.DataFrame(cv)
    cv.to_csv(out / "a94_cue_core.csv", index=False)
    res["cue_core_d_log_ratio_median"] = float(cv.d_log_cue_over_core.median())
    res["cue_core_d_offset_h_median"] = float(cv.d_cue_minus_core_h.median())

    # ---- A99 sex
    tis99, nulls = BODY, {}
    for k, cc in C.items():
        dn = cc.donors.dropna(subset=["age"])
        sexlab = dn.sex.astype(int)
        strata = dn.hardy.astype(str) + "|" + labs[k].reindex(dn.index).astype(str)
        stat = lambda lab, cc=cc: A.contrast(cc, lab, 1, 0, A.radius, tis99, seed)
        o = stat(sexlab)
        nulls[k] = A.perm_null(stat, sexlab, strata, rng, n_perm)
        res[f"H99_{k}"], res[f"H99_{k}_p"] = o, A.p_two(o, nulls[k])
        for g in (0, 1):
            d = sexlab.index[sexlab == g].to_numpy()
            row = dict(test="A99", cohort=k, group="female" if g else "male", donors=len(d),
                       log_radius=np.nanmedian([A.radius(cc, t, d) for t in tis99 if t in cc.t]),
                       C_body=A.coherence(cc, d, tis99))
            if k == "fast":
                row.update(R_wall=np.nanmedian([A.walltime_R(cc, t, d) for t in tis99 if t in cc.t]),
                           C_brain_body=A.coherence(cc, d, brain, tis99))
            rows.append(row)
    m = min(len(nulls["fast"]), len(nulls["vent"]))
    o = (res["H99_fast"] + res["H99_vent"]) / 2
    res["H99"], res["H99_p"] = o, A.p_two(o, (nulls["fast"][:m] + nulls["vent"][:m]) / 2)
    sx = pd.DataFrame([dict(tissue=t, **timed_slope(x.assign(fem=x.sex - x.sex.mean()), "fem"))
                       for t, x in fast_s.groupby("tissue") if t in tis99 and len(x) >= (40 if not run["quick"] else 20)])
    sx.to_csv(out / "a99_timed_sex.csv", index=False)
    for col in ("log_amp", "phase_h"):
        r = cs.stats.meta(sx[col], sx[f"se_{col}"])
        res[f"timed_{col}_female_minus_male"], res[f"timed_{col}_sex_se"] = r["fixed"], r["se_fixed"]

    # ---- A96 part 1 (BMI, T2D, liver disease)
    mt = []
    for t in METAB + COMPANION:
        x = fast_s[fast_s.tissue == t].dropna(subset=["bmi", "age"])
        if len(x) < A.MIN_N:
            continue
        x = x.assign(bmi_z=(x.bmi - x.bmi.mean()) / x.bmi.std())
        f = timed_slope(x, "bmi_z", extra_age=True)
        row = dict(tissue=t, **f)
        if t == "Liver":
            hard, null = x.hardy.to_numpy(), []
            for _ in range(n_bmi):
                xb, z = x.copy(), x.bmi_z.to_numpy().copy()
                for h in np.unique(hard):
                    i = np.flatnonzero(hard == h)
                    z[i] = rng.permutation(z[i])
                xb["bmi_z"] = z
                null.append(timed_slope(xb, "bmi_z", extra_age=True)["log_amp"])
            res["H96_1_liver_log_amp_per_sd_bmi"], res["H96_1_n"] = f["log_amp"], len(x)
            null = np.asarray(null)
            res["H96_1_p_lower"] = float((1 + np.sum(null <= f["log_amp"])) / (1 + len(null)))
        for k, cc in C.items():
            if t not in cc.t:
                continue
            d = cc.t[t]
            info = cc.donors.reindex(d["subj"])
            ok = info.bmi.notna().to_numpy()
            if ok.sum() > 10:
                z = (info.bmi[ok] - info.bmi[ok].mean()) / info.bmi[ok].std()
                row[f"radius_slope_per_sd_bmi_{k}"] = float(np.polyfit(z, d["r2"][ok], 1)[0] / d["r2"][ok].mean())
            for g in ("t2d", "liver_dis"):
                y = info[g].to_numpy()
                a_, b_ = d["r2"][y == 1], d["r2"][y == 0]
                row[f"n_{g}_{k}"] = int(len(a_))
                row[f"log_radius_{g}_{k}"] = (float(np.log(a_.mean()) - np.log(b_.mean()))
                                              if len(a_) >= 10 and a_.mean() > 0 and b_.mean() > 0 else np.nan)
        mt.append(row)
    pd.DataFrame(mt).to_csv(out / "a96_gtex_metabolic.csv", index=False)
    pd.DataFrame(rows).to_csv(out / "a99_groups.csv", index=False)
    v = pd.Series(res)
    v["n_perm"] = n_perm
    v["H99_verdict"] = "PASS" if v.H99_p < 0.05 else "FAIL"
    if "H96_1_p_lower" in v:
        v["H96_1_verdict"] = "PASS" if (v.H96_1_liver_log_amp_per_sd_bmi < 0 and v.H96_1_p_lower < 0.05) else "FAIL"
    v.to_csv(out / "a99_a96_summary.csv", header=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).round(3).to_string(index=False), "\n", v.to_string(), "\n",
          pd.DataFrame(mt).round(3).to_string(index=False), "\n", cv.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
