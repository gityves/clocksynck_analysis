"""A146 part 2 (pan-cancer) blind preparation, step 3: inventory, power and omnibus null.

BLIND: no program score is related to a time or phase anywhere here. Observed times enter only as the per-cohort set of
values that is PERMUTED; the unpermuted pairing of score and time is never evaluated. Inputs: output/part2pc_cache/
(host_cohorts.csv, host_times.csv from part2pc_host.py; <study>.csv.gz from part2pc_extract.py).

Planned statistic (A146 omnibus, applied within cohort, pooled by stacking): per program k, scores z-scored within
cohort over the patients that have the program compartment and an internal time; stacked OLS of the score on
cos(2 pi (t - h_k) / 24), cosine centred within cohort; t_k = slope t statistic (df = n - cohorts - 1); S = mean t_k.
Null: the internal times of each cohort permuted among its timed patients (one permutation per cohort and replicate, shared
by all programs), 2,000 replicates, seed 146.

Programs (compartment): P1, P2, P3, P9 malignant (CNV call where >= 5 patients have it, else tumour epithelium);
P4, P10 T_NK, residual on the T_NK cell fraction; P5, P6 myeloid, residual on the myeloid fraction; P7 T_NK fraction; P8
myeloid fraction (fractions of ALL cells of the patient's tumour samples, incl. unassigned). The omnibus (10 programs) has
P9 on malignant cells; the sensitivity `omnibus11` adds P9 in myeloid cells.
CHOICES (not in the registration, listed in the report): gene score = log2(count / depth * 1e6 + 1), mean of within-cohort gene z,
re-z-scored; fraction score = log2((n + .5) / (N + 1)); a patient needs >= 100 tumour cells (fractions) and >= 50 cells in the
compartment (pseudobulk table rule); a cohort enters the pool if its gate passes and it has >= 6 timed patients, a
cohort x program unit if it has >= 6 such patients; a program is scored in a cohort if >= half of its genes are present.

    python part2pc_prep.py [--n-perm 2000]
Outputs: output/part2pc_prep_{cohorts,compartments,power,keeping_split,omnibus_summary,omnibus_null}.csv
"""
from __future__ import annotations
import sys, argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm

REPO = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(REPO / "paper/exploratory/A138_cancer_estimability/scripts"), str(Path(__file__).resolve().parent)]
import q3_lib as q  # noqa: E402
from part2pc_extract import PROGRAMS  # noqa: E402

BASE = Path(__file__).resolve().parents[1] / "output"
CACHE = BASE / "part2pc_cache"
SEED, MIN_UNIT, MIN_CELLS_ALL = 146, 6, 100
HUMAN_H = dict(P1=20.5, P2=17, P3=17, P4=13, P5=13, P6=13, P7=11, P8=23, P9=17, P10=23)
# (program, compartment kind, abundance for residual)
UNITS = [("P1", "mal", None), ("P2", "mal", None), ("P3", "mal", None), ("P9", "mal", None), ("P4", "T_NK", "T_NK"),
         ("P5", "myeloid", "myeloid"), ("P6", "myeloid", "myeloid"), ("P7", "T_NK", "frac"), ("P8", "myeloid", "frac"),
         ("P10", "T_NK", "T_NK")]
EXTRA = [("P9m", "myeloid", None)]  # sensitivity unit (P9 in myeloid cells, predicted 17 h)
HIGH = {"BLCA", "LIHC", "LUAD", "PRAD", "KIRC"}; LOW = {"UCEC", "COAD"}
# cohort cancer label -> TCGA project of the H138.4 relation (conf: direct / histology not recorded)
MAP = {"colorectal": ("COAD", "direct"), "bladder": ("BLCA", "direct"), "prostate": ("PRAD", "direct"),
       "liver": ("LIHC", "histology_unrecorded"), "liver/biliary": ("LIHC", "histology_unrecorded"),
       "lung": ("LUAD", "histology_unrecorded"), "kidney": ("KIRC", "histology_unrecorded"),
       "breast": ("BRCA", "direct"), "thyroid": ("THCA", "direct")}
MDE_Z = norm.ppf(.95) + norm.ppf(.8)


def mde(n):
    return float(np.tanh(MDE_Z / np.sqrt(n - 3))) if n > 3 else np.nan


def keeping(cancer: str):
    k, conf = MAP.get(cancer.replace(" cancer", "").strip(), (None, ""))
    cls = "keeping" if k in HIGH else "not_keeping" if k in LOW else "unmapped"
    return k, conf, cls


def cohort_inputs(study: str, host: pd.Series, mal: str):
    """Per-patient inputs of one cohort (tumour samples): cell counts, program gene log-CPM per compartment, fractions."""
    c = pd.read_csv(CACHE / f"{study}.csv.gz")
    c["patient"] = c.patient.astype(str)
    tot = c.groupby("patient").n_cells.sum()
    epi = c[c.compartment == "epithelial"]
    m = epi[epi.call == "malignant"] if mal == "tumour:malignant" else epi
    kinds = {"mal": m, "T_NK": c[c.compartment == "T_NK"], "myeloid": c[c.compartment == "myeloid"]}
    return c, tot, {k: v.groupby("patient").sum(numeric_only=True) for k, v in kinds.items()}


def gscore(df: pd.DataFrame, genes) -> pd.Series | None:
    g = [x for x in genes if x in df and df[x].notna().any()]
    if len(g) < len(genes) / 2 or len(df) < 2:
        return None
    lc = np.log2(df[g].div(df.depth, axis=0) * 1e6 + 1)
    z = (lc - lc.mean()) / lc.std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1)


def zs(x: pd.Series) -> pd.Series:
    return (x - x.mean()) / x.std(ddof=0)


def units_of(study, c, tot, comp, pbc, timed):
    """{unit: z-scored score Series over eligible patients} of a cohort. Eligible: internal time, compartment in the
    pseudobulk table (>= 50 cells, builder rule), tumour cells >= MIN_CELLS_ALL."""
    out = {}
    ok_all = set(tot[tot >= MIN_CELLS_ALL].index) & set(timed)
    frac = {k: np.log2((comp[k].n_cells.reindex(tot.index).fillna(0) + .5) / (tot + 1)) for k in ("T_NK", "myeloid")}
    for prog, kind, ab in UNITS + EXTRA:
        pts = sorted(ok_all & set(comp[kind].index) & pbc[kind])
        if len(pts) < MIN_UNIT:
            continue
        if prog in ("P7", "P8"):
            s = frac[kind].loc[pts]
        else:
            s = gscore(comp[kind].loc[pts], PROGRAMS[prog.rstrip("m")].split())
            if s is None:
                continue
            s = zs(s)
            if ab:
                a = zs(frac[ab].loc[pts]).to_numpy()
                X = np.column_stack([np.ones(len(a)), a])
                s = pd.Series(s.to_numpy() - X @ np.linalg.lstsq(X, s.to_numpy(), rcond=None)[0], index=pts)
        s = zs(s)
        if s.notna().all():
            out[prog] = s
    return out


def omnibus_null(cohorts: dict, progs: list, n_perm: int, seed: int) -> np.ndarray:
    """cohorts: {study: (times array over timed patients, {prog: (idx into times, z score array)})}. Returns S per perm."""
    rng = np.random.default_rng(seed)
    W = 2 * np.pi / 24
    S = np.zeros(n_perm)
    studies = list(cohorts)
    perms = {s: np.argsort(rng.random((n_perm, len(cohorts[s][0]))), axis=1) for s in studies}
    for k, pr in enumerate(progs):
        sxy = np.zeros(n_perm); sxx = np.zeros(n_perm); syy = 0.0; n = 0; nc = 0
        for s in studies:
            t, u = cohorts[s]
            if pr not in u:
                continue
            idx, y = u[pr]
            cos = np.cos(W * (t[perms[s]][:, idx] - HUMAN_H[pr.rstrip("m") if pr != "P9m" else "P9"]))
            cos = cos - cos.mean(1, keepdims=True)
            sxy += cos @ y; sxx += (cos ** 2).sum(1); syy += (y ** 2).sum(); n += len(y); nc += 1
        if n - nc - 1 < 3:
            continue
        b = sxy / sxx
        s2 = (syy - b * sxy) / (n - nc - 1)
        S += b / np.sqrt(s2 / sxx)
    return S / len(progs)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--n-perm", type=int, default=2000); a = ap.parse_args()
    H = pd.read_csv(CACHE / "host_cohorts.csv")
    T = pd.read_csv(CACHE / "host_times.csv", dtype={"patient": str})
    rows, comp_rows, pools = [], [], {}
    for r in H.itertuples():
        study = r.study
        pb, mal = q.study_table(study)
        tum = pb[pb.compartment.str.startswith("tumour:")]
        c, tot, comp = cohort_inputs(study, None, mal)
        cancer = c.cancer.iloc[0]
        tcga, conf, cls = keeping(cancer)
        pbc = {"mal": set(tum[tum.compartment == mal].donor.astype(str)), "T_NK": set(tum[tum.compartment == "tumour:T_NK"].donor.astype(str)),
               "myeloid": set(tum[tum.compartment == "tumour:myeloid"].donor.astype(str))}
        timed = T[T.study == study]
        tmap = dict(zip(timed.patient, timed.internal_h))
        npat = {k: len(v) for k, v in pbc.items()}
        npat_t = {k: len(v & set(tmap)) for k, v in pbc.items()}
        frac_ok = tot[tot >= MIN_CELLS_ALL]
        tf = comp["T_NK"].n_cells.reindex(tot.index).fillna(0) / tot
        mf = comp["myeloid"].n_cells.reindex(tot.index).fillna(0) / tot
        pbn = tum[tum.compartment == mal].set_index(tum[tum.compartment == mal].donor.astype(str)).n_cells
        agree = float((comp["mal"].n_cells.reindex(pbn.index) == pbn).mean()) if len(pbn) else np.nan
        rows.append(dict(
            study=study, cancer_type=cancer, breast_cohort=("breast" in cancer) or (study == "SRP329496"),
            breast_note=("10 breast patients inside an ovarian-labelled cohort" if study == "SRP329496" else ""),
            n_patients_tumour=r.n_patients_tumour, n_patients_stromal_host=r.n_patients_stromal_host,
            n_anchors=r.n_anchors, anchors=r.anchors, anchor_source=r.anchor_source, gate_clock_share=r.gate_clock_share,
            gate_share_null=r.gate_share_null, gate_p_clock_share=r.gate_p_clock_share, gate_p_pairing=r.gate_p_pairing,
            gate_passes=r.gate_passes, gate_error=getattr(r, "error", ""), n_patients_internal_time=r.n_patients_internal_time,
            time_circ_sd_h=r.time_circ_sd_h, time_window_w_h=r.window_w_h,
            malignant_compartment=mal, n_patients_malignant=npat["mal"], n_patients_TNK=npat["T_NK"], n_patients_myeloid=npat["myeloid"],
            n_timed_malignant=npat_t["mal"], n_timed_TNK=npat_t["T_NK"], n_timed_myeloid=npat_t["myeloid"],
            n_patients_frac_denominator=len(frac_ok), n_timed_frac_denominator=len(set(frac_ok.index) & set(tmap)),
            median_cells_per_tumour_patient=float(tot.median()), median_T_NK_fraction=float(tf[frac_ok.index].median()),
            median_myeloid_fraction=float(mf[frac_ok.index].median()), cache_vs_pb_ncells_agree=agree,
            tcga_project=tcga, map_confidence=conf, clock_class=cls))
        if r.gate_passes is True or r.gate_passes == True:  # noqa: E712
            if len(tmap) >= MIN_UNIT:
                pools[study] = (np.array([tmap[p] for p in sorted(tmap)]), sorted(tmap),
                                units_of(study, c, tot, comp, pbc, sorted(tmap)))
    D = pd.DataFrame(rows)
    D.to_csv(BASE / "part2pc_prep_cohorts.csv", index=False)
    D["eligible_pool"] = D.study.isin(pools)

    # compartments (counts) summary
    cols = ["n_patients_tumour", "n_patients_stromal_host", "n_patients_internal_time", "n_patients_malignant", "n_patients_TNK",
            "n_patients_myeloid", "n_timed_malignant", "n_timed_TNK", "n_timed_myeloid", "n_patients_frac_denominator",
            "n_timed_frac_denominator"]
    comp_rows = []
    for name, sub in {"all_cohorts": D, "gate_passing": D[D.gate_passes == True], "pool": D[D.eligible_pool],  # noqa: E712
                      "pool_keeping": D[D.eligible_pool & (D.clock_class == "keeping")],
                      "pool_not_keeping": D[D.eligible_pool & (D.clock_class == "not_keeping")],
                      "pool_unmapped": D[D.eligible_pool & (D.clock_class == "unmapped")]}.items():
        comp_rows.append(dict(set=name, n_cohorts=len(sub), **{c: int(sub[c].fillna(0).sum()) for c in cols}))
    pd.DataFrame(comp_rows).to_csv(BASE / "part2pc_prep_compartments.csv", index=False)

    # keeping split
    ks = []
    for name, sub in {"all_cohorts": D, "gate_passing": D[D.gate_passes == True], "pool": D[D.eligible_pool]}.items():  # noqa: E712
        for cl, g in sub.groupby("clock_class"):
            ks.append(dict(set=name, clock_class=cl, n_cohorts=len(g), n_patients_tumour=int(g.n_patients_tumour.sum()),
                           n_patients_internal_time=int(g.n_patients_internal_time.fillna(0).sum()),
                           cohorts=";".join(g.study + "(" + g.tcga_project.fillna("-") + ")")))
    pd.DataFrame(ks).to_csv(BASE / "part2pc_prep_keeping_split.csv", index=False)

    # power per program (pooled over the pool; n = patients with the program compartment, time and denominator)
    prow = []
    for sets, keepf in {"pool": lambda s: True,
                        "pool_keeping": lambda s: D.set_index("study").clock_class[s] == "keeping"}.items():
        for prog, kind, ab in UNITS + EXTRA:
            ns = {s: len(v[2][prog]) for s, v in pools.items() if prog in v[2] and keepf(s)}
            n = sum(ns.values())
            prow.append(dict(set=sets, program=prog, compartment={"mal": "malignant/tumour epithelium"}.get(kind, kind),
                             n_cohorts=len(ns), n_patients=n, mde_abs_r=mde(n), cohorts=";".join(f"{s}:{k}" for s, k in ns.items())))
    pd.DataFrame(prow).to_csv(BASE / "part2pc_prep_power.csv", index=False)

    # omnibus null (permuted times only)
    summ, nulls = [], {}
    for name, keepf, progs in [("omnibus10", lambda s: True, [u[0] for u in UNITS]),
                               ("omnibus11", lambda s: True, [u[0] for u in UNITS] + ["P9m"]),
                               ("omnibus10_keeping", lambda s: D.set_index("study").clock_class[s] == "keeping", [u[0] for u in UNITS])]:
        co = {}
        for s, (t, ids, un) in pools.items():
            if not keepf(s):
                continue
            ix = {p: i for i, p in enumerate(ids)}
            co[s] = (t, {pr: (np.array([ix[p] for p in sc.index]), sc.to_numpy()) for pr, sc in un.items() if pr in progs})
        used = [p for p in progs if any(p in v[1] for v in co.values())]
        S = omnibus_null(co, used, a.n_perm, SEED)
        nulls[name] = S
        summ.append(dict(statistic=name, n_cohorts=len(co), n_programs=len(used), programs=";".join(used), n_perm=a.n_perm,
                         null_mean=float(S.mean()), null_sd=float(S.std(ddof=1)), null_q95=float(np.quantile(S, .95)),
                         mde=float(np.quantile(S, .95) + 0.84 * S.std(ddof=1)), seed=SEED))
    pd.DataFrame(summ).to_csv(BASE / "part2pc_prep_omnibus_summary.csv", index=False)
    pd.DataFrame(nulls).rename_axis("perm").to_csv(BASE / "part2pc_prep_omnibus_null.csv")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(pd.DataFrame(comp_rows).to_string()); print(pd.DataFrame(prow).drop(columns="cohorts").to_string()); print(pd.DataFrame(summ).drop(columns="programs").to_string())
