"""A305 power stage: malignant vs own-tumour stromal-vascular pairwise clock sharing, PLACEBO DRAWS ONLY.

Registered in paper/docs/registered_experiments.md (A305). T_MS = median over cohorts of the cohort-median clock |r|
(uncorrected, `clocksync.pairwise`) of the (malignant, tumour SV) pairs; the clock's percentile pi among 200 matched and
200 responsive companion sets (A283 form). Cohorts: the A227 list (27 A149 cohorts with estimable S1 + 7 A226 cohorts),
tumour-sample compartments only. Strata = cohort x chemistry (table `stratum`).

NO observed statistic: the identity pairing (true donor labels) is never evaluated for |r|, T_MS or any companion value.
Only `pairwise.placebo` draws (donor labels permuted within compartment x stratum, seeds [305, d]) enter `abs_r`.
Design counts (patients, pairs, df) come from `pair_index` on the table keys only.

    python power_a305.py [--draws 2000] [--jobs 4] [--smoke]
Outputs (../output/power/): design.csv, placebo_T.npz (draws x cohorts x sets), placebo_draws.csv, power.csv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
for p in ("A149_paired_tumour_normal", "A151_fibroblast_pairing", "A226_stroma_pairing_confirm",
          "A227_tme_compartment_pairing"):
    sys.path.insert(0, str(PAPER / "exploratory" / p / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync import pairwise as PW  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.waveform import gtex  # noqa: E402
from common import cli  # noqa: E402

import a149_lib as L  # noqa: E402
import compute_a226 as K  # noqa: E402
import compute_a227 as A227  # noqa: E402

OUT = HERE.parent / "output" / "power"
SEED, N_SET = 305, 200
SV = ("fibroblast", "endothelial", "mural")
FAM = {"matched": slice(1, 1 + N_SET), "responsive": slice(1 + N_SET, 1 + 2 * N_SET)}
MIN_P, MAX_MDE = 0.05, 0.10
KEY = ["donor", "compartment", "stratum"]


def log(*a):
    print(pd.Timestamp.now().strftime("%H:%M:%S"), *a, flush=True)


def mal_of(study):
    return L.table(study)[1] if study not in K.COHORTS else L.EPI


def genome(study, comps, keys):
    """Genome-wide pseudobulks (rows of `keys` = donor x compartment kept in the cohort table), labels as `cancer_sc`."""
    K.install()
    from cancer_sc import build as B
    us, _ = L.ST.units(study)
    us = [u for u in L.ST.majority_chemistry(us) if u[1] in set(keys.donor)]
    calls = None
    if L.MAL in comps:
        calls = L.new_calls(study) if study in L.NEW else L.ST.read_calls(study)
    sums, depth, genes = {}, {}, None
    for sample, patient, kind, load, chem, _ in us:
        X, _g, umi, comp, barcodes, sym = load()
        lab = B._labels(comp, kind, sample, barcodes, calls)
        genes = pd.Index(sym) if genes is None else genes
        for c in comps:
            i = np.flatnonzero(lab == c)
            if len(i):
                k = (patient, c)
                sums[k] = sums.get(k, 0) + np.asarray(X[i].sum(0)).ravel()
                depth[k] = depth.get(k, 0) + float(umi[i].sum())
    ks = sorted(k for k in sums if (k[0], k[1]) in set(zip(keys.donor, keys.compartment)))
    G = pd.DataFrame(np.vstack([sums[k] for k in ks]), columns=[canonical(str(g)) for g in genes])
    G = G.T.groupby(level=0).sum().T
    st = keys.set_index(["donor", "compartment"]).stratum
    meta = pd.DataFrame({"donor": [k[0] for k in ks], "compartment": [k[1] for k in ks],
                         "stratum": [st.loc[k] for k in ks], "depth": [depth[k] for k in ks]})
    return pd.concat([meta, G.reset_index(drop=True)], axis=1)


def cohort(study, k):
    """Design rows and the state matrix Z (sets x rows) of one cohort; None if no (malignant, SV) pair passes."""
    from clocksync.config import limit_threads
    limit_threads(1)
    warnings.filterwarnings("ignore")
    pb = A227.table(study)
    mal = mal_of(study)
    comps = set(pb.compartment)
    sv = sorted(c for c in comps if c.startswith("tumour:") and c.split(":")[1] in SV)
    tnk = "tumour:T_NK" if "tumour:T_NK" in comps else None
    if mal not in comps or not sv:
        return study, None, dict(study=study, malignant=mal, note="no malignant or no SV compartment")
    keep = [mal] + sv + ([tnk] if tnk else [])
    keys = pb[pb.compartment.isin(keep)][KEY].reset_index(drop=True)
    if "stratum" not in pb:
        raise ValueError(study)
    pi_all = PW.pair_index(keys, min_stratum=3, min_df=15)
    rows = []
    for p in pi_all:
        cls = ("MS" if mal in (p["a"], p["b"]) and (p["a"] in sv or p["b"] in sv) else
               "MT" if mal in (p["a"], p["b"]) and tnk in (p["a"], p["b"]) else
               "SS" if p["a"] in sv and p["b"] in sv else None)
        if cls:
            rows.append(dict(study=study, cancer=A227.cancer(study), klass=cls, a=p["a"], b=p["b"], donors=p["n"], df=p["df"],
                             strata=int(p["s"].max()) + 1))
    d = pd.DataFrame(rows)
    pat = keys[keys.compartment == mal].donor
    sv_pat = keys[keys.compartment.isin(sv)].donor
    info = dict(study=study, malignant=mal, patients_total=pb.donor.nunique(), patients_malignant=pat.nunique(),
                patients_mal_and_sv=len(set(pat) & set(sv_pat)), sv_compartments=";".join(sv),
                strata=keys.stratum.nunique())
    if d.empty or not (d.klass == "MS").any():
        return study, None, dict(info, design=d, note="no (malignant, SV) pair passes min_stratum 3, min_df 15")
    # state matrix on the genome-wide table restricted to malignant + SV rows
    ms_keys = keys[keys.compartment.isin([mal] + sv)].reset_index(drop=True)
    G = genome(study, [mal] + sv, ms_keys)
    wf = gtex()
    wfs = [wf]
    for j, kind in enumerate(PW.FAMILIES):
        wfs += PW.companions(G, N_SET, kind, np.random.default_rng([SEED, 7, k, j]), wf)
    Z = np.array([PW.states(G, w) for w in wfs])
    key = G[KEY].reset_index(drop=True)
    return study, (key, Z.astype(np.complex64), mal, sv), dict(info, design=d, note="")


def draw(seed, data):
    """T per cohort x set for one placebo draw (donor labels permuted within compartment x stratum; never the true pairing)."""
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(seed)
    out = []
    for key, Z, mal, sv in data:
        pk = PW.placebo(key, rng)
        pairs = [p for p in PW.pair_index(pk, min_stratum=3, min_df=15) if (p["a"] == mal) != (p["b"] == mal)
                 and (p["a"] in sv or p["b"] in sv)]
        out.append(np.nanmedian(PW.abs_r(Z.astype(np.complex128), pairs), 1))
    return np.array(out, np.float32)


def pi_of(T):
    """T: draws x sets -> pi per family."""
    return {f: PW.percentile(T, sl) for f, sl in FAM.items()}


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    args.jobs = min(args.jobs, 4)
    cli.setup(args)
    OUT.mkdir(parents=True, exist_ok=True)
    cohorts = A227.cohort_list(args.smoke)
    res = Parallel(n_jobs=args.jobs)(delayed(cohort)(s, i) for i, s in enumerate(cohorts))
    drows, data, names = [], [], []
    for s, d, info in res:
        dsg = info.pop("design", None)
        if dsg is not None and len(dsg):
            drows.append(dsg)
        info["pairs_MS"] = int((dsg.klass == "MS").sum()) if dsg is not None and len(dsg) else 0
        info["pairs_SS"] = int((dsg.klass == "SS").sum()) if dsg is not None and len(dsg) else 0
        info["pairs_MT"] = int((dsg.klass == "MT").sum()) if dsg is not None and len(dsg) else 0
        info["df_MS_total"] = int(dsg[dsg.klass == "MS"].df.sum()) if info["pairs_MS"] else 0
        info["in_T_MS"] = d is not None
        drows_info = info
        names.append(drows_info)
        if d is not None:
            data.append(d)
    coh = pd.DataFrame(names)
    coh.to_csv(OUT / "design.csv", index=False)
    if drows:
        pd.concat(drows).to_csv(OUT / "design_pairs.csv", index=False)
    log("design: cohorts", len(coh), "in T_MS", int(coh.in_T_MS.sum()), "pairs MS", int(coh.pairs_MS.sum()))
    in_names = coh[coh.in_T_MS].study.tolist()
    if not data:
        log("no cohort passes the pair rule; stop")
        return

    T = np.empty((args.draws, len(data), 1 + 2 * N_SET), np.float32)
    chunk = 100
    for a in range(0, args.draws, chunk):
        ds = range(a, min(a + chunk, args.draws))
        r = Parallel(n_jobs=args.jobs)(delayed(draw)([SEED, d], data) for d in ds)
        T[a:a + len(r)] = np.array(r)
        log("draws", a + len(r), "of", args.draws)
    np.savez_compressed(OUT / "placebo_T.npz", T=T, cohorts=np.array(in_names))
    TMS = np.nanmedian(T, 1)                       # draws x sets: median over cohorts of the cohort median |r|
    P = pi_of(TMS)
    D = TMS[:, 0] - np.nanmedian(TMS[:, FAM["matched"]], 1)
    pd.DataFrame(dict(draw=np.arange(args.draws), T_clock=TMS[:, 0], T_matched_median=np.nanmedian(TMS[:, FAM["matched"]], 1),
                      D=D, pi_matched=P["matched"], pi_responsive=P["responsive"])).to_csv(OUT / "placebo_draws.csv", index=False)
    mu, sd, q95 = float(D.mean()), float(D.std(ddof=1)), float(np.quantile(D, 0.95))
    mde = q95 + 0.84 * sd - mu
    row = dict(draws=args.draws, cohorts=len(in_names), pairs_MS=int(coh.pairs_MS.sum()),
               patients_mal_and_sv=int(coh[coh.in_T_MS].patients_mal_and_sv.sum()),
               **{f"placebo_rate_pi1_{f}": float(np.mean(P[f] == 1.0)) for f in FAM},
               **{f"smallest_p_cal_{f}": float((1 + np.sum(P[f] == 1.0)) / (1 + args.draws)) for f in FAM},
               D_mean=mu, D_sd=sd, D_q95=q95, MDE=mde)
    row["flag_p_matched_lt_0.05"] = row["smallest_p_cal_matched"] < MIN_P
    row["flag_MDE_le_0.10"] = mde <= MAX_MDE
    row["scored"] = row["flag_p_matched_lt_0.05"] and row["flag_MDE_le_0.10"]
    pd.DataFrame([row]).to_csv(OUT / "power.csv", index=False)
    print(pd.DataFrame([row]).T.to_string())
    log("done")


if __name__ == "__main__":
    main()
