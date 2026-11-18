"""A317 (registered): does brain death remove the clock relation between organs, clock-specifically and beyond procurement
confounds? GTEx Hardy 0 organ donors vs ventilated post-mortem donors; fast deaths as the arrangement reference.

A109 `Tissues(200, seed 0)` with responsive sets; A132 groups; organ donors drawn to the post-mortem n (20 draws).
H317.1: D_rho = rho(post-mortem, fast) - rho(organ donors, fast) (Fisher-Lee of pair phases, A132 / A293 estimator), with
a donor-label permutation null and the same difference for matched and responsive sets. H317.2: D_C = excess C
(post-mortem) - excess C (organ donors), excess = C_clock - median matched C (A109 H109.1 estimator). Sensitivities:
steroid users excluded, common ischaemic-time range, RIN >= 6. Aggregates only (GTEx controlled access).

    python score_a317.py --stage power | score
Outputs (`../output/`): power.csv | tests.csv, sensitivity.csv, sets.csv.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A109_coupling_nulls" / "scripts"))
sys.path.insert(0, str(EXP / "A132_brain_death" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
import compute_brain_death as B  # noqa: E402
import compute_coupling_nulls as N  # noqa: E402
import compute_h109_responsive as RS  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
N_DRAWS, N_PERM, SEED = 20, 2000, 317


def raw(T, kind, k, R):
    Z = np.zeros((len(T.donors), len(H.BODY)), complex)
    M = np.zeros((len(T.donors), len(H.BODY)))
    for j, t in enumerate(H.BODY):
        r = R[t]
        m = r >= 0
        Z[m, j] = T.Z[kind][t][k][r[m]]
        M[m, j] = 1
    return Z, M


def rho(Z, M, ref, g, mask):
    return float(ST.fisher_lee(B.phase(Z, M, ref, mask), B.phase(Z, M, g, mask)))


def d_rho(Z, M, fast, pm, od_draws, mask):
    return rho(Z, M, fast, pm, mask) - float(np.median([rho(Z, M, fast, d, mask) for d in od_draws]))


def excess_C(T, pm, od_draws, R):
    """Excess C (clock minus median matched) per group, and the responsive excess."""
    out = {}
    for g, draws in (("pm", [pm]), ("od", od_draws)):
        c = H.coherence(T, "clock", draws, R)[0]
        cm = H.coherence(T, "pseudo", draws, R)
        cr = H.coherence(T, "resp", draws, R)
        out[g] = dict(C=float(c), excess=float(c - np.nanmedian(cm)), pct_matched=float(np.nanmean(cm < c)),
                      resp_excess=float(c - np.nanmedian(cr)), pct_resp=float(np.nanmean(cr < c)))
    return out


def restrict(T, keep_rows):
    """R with rows not in keep_rows[t] (boolean per row) set to absent."""
    R = {}
    for t in H.BODY:
        r = T.R[t].copy()
        bad = (r >= 0) & ~keep_rows[t][np.maximum(r, 0)]
        r[bad] = -1
        R[t] = r
    return R


def sample_attrs(T):
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    o = a.obs[["subject", "tissue", "ischemic_time", "rin"]].copy()
    a.file.close()
    o["subject"], o["tissue"] = o.subject.astype(str), o.tissue.astype(str)
    o = o.drop_duplicates(["subject", "tissue"])  # first sample of a donor, as A109 Tissues
    return {t: o[o.tissue == t].set_index("subject").reindex(T.row_donor[t]) for t in H.BODY}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    RS.add_sets(T, 200, 0)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    rng = np.random.default_rng(SEED)
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(N_DRAWS)]
    Zc, Mc = raw(T, "clock", 0, T.R)
    mask = B.pair_mask(Mc, fast, pm, *od_draws[:1]) & B.pair_mask(Mc, od)
    print(f"states built ({time.time() - t0:.0f}s): od {len(od)}, pm {len(pm)}, pairs {int(mask.sum())}", flush=True)

    def perm_null(n):
        rho_n, c_n = [], []
        for _ in range(n):
            p = rng.permutation(h0)
            gp, go = p[:len(pm)], p[len(pm):2 * len(pm)]
            rho_n.append(rho(Zc, Mc, fast, gp, mask) - rho(Zc, Mc, fast, go, mask))
            c_n.append(H.coherence(T, "clock", [gp])[0] - H.coherence(T, "clock", [go])[0])
        return np.array(rho_n), np.array(c_n)

    if a.stage == "power":
        rn, cn = perm_null(N_PERM)
        res = dict(pairs=int(mask.sum()), sd_D_rho=float(np.nanstd(rn, ddof=1)), mde_D_rho=float(2.49 * np.nanstd(rn, ddof=1)),
                   a132b_D_rho=0.517, sd_D_C=float(np.nanstd(cn, ddof=1)), mde_D_C=float(2.49 * np.nanstd(cn, ddof=1)))
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).to_string(), flush=True)
        return

    rn, cn = perm_null(N_PERM)
    D = d_rho(Zc, Mc, fast, pm, od_draws, mask)
    sets = []
    for kind, lab in (("pseudo", "matched"), ("resp", "responsive")):
        for k in range(200):
            Z, M = raw(T, kind, k, T.R)
            sets.append(dict(family=lab, set=k, D_rho=d_rho(Z, M, fast, pm, od_draws, mask)))
    sets = pd.DataFrame(sets)
    sets.to_csv(OUT / "sets.csv", index=False)
    q = sets.groupby("family").D_rho.quantile(0.95)
    ex = excess_C(T, pm, od_draws, T.R)
    D_C = ex["pm"]["excess"] - ex["od"]["excess"]
    p_rho = float((1 + np.sum(rn >= D)) / (1 + len(rn)))
    p_C = float((1 + np.sum(cn >= D_C)) / (1 + len(cn)))
    tests = [dict(test="H317.1", D_rho=D, p_labels=p_rho, matched_q95=float(q["matched"]), responsive_q95=float(q["responsive"]),
                  pct_matched=float(np.mean(sets[sets.family == "matched"].D_rho < D)),
                  pct_responsive=float(np.mean(sets[sets.family == "responsive"].D_rho < D)),
                  passes=bool(D > 0 and p_rho < 0.05 and D > q["matched"] and D > q["responsive"])),
             dict(test="H317.2", D_C=D_C, p_labels=p_C, passes=bool(D_C > 0 and p_C < 0.05),
                  **{f"{g}_{k}": v for g in ("pm", "od") for k, v in ex[g].items()})]
    pd.DataFrame(tests).to_csv(OUT / "tests.csv", index=False)
    # sensitivities
    attrs = sample_attrs(T)
    steroid = pd.to_numeric(ph.MHSTRDLT, errors="coerce").to_numpy() == 1
    isch_pm = np.concatenate([attrs[t].ischemic_time.to_numpy(float)[T.R[t][pm][T.R[t][pm] >= 0]] for t in H.BODY])
    lo, hi = np.nanquantile(isch_pm, [0.05, 0.95])
    sens = []
    for name in ("steroids_excluded", "ischaemic_common_range", "rin_ge_6"):
        if name == "steroids_excluded":
            keep = ~steroid
            pm_s, od_s, R = pm[keep[pm]], od[keep[od]], T.R
        else:
            col = "ischemic_time" if name.startswith("isch") else "rin"
            rows = {t: (attrs[t][col].to_numpy(float) >= lo) & (attrs[t][col].to_numpy(float) <= hi) if col == "ischemic_time"
                    else attrs[t][col].to_numpy(float) >= 6 for t in H.BODY}
            R, pm_s, od_s = restrict(T, rows), pm, od
        Z, M = raw(T, "clock", 0, R)
        dr = [rng.choice(od_s, min(len(pm_s), len(od_s)), replace=False) for _ in range(N_DRAWS)]
        m = mask & B.pair_mask(M, fast, pm_s, dr[0])
        row = dict(sensitivity=name, n_pm=len(pm_s), n_od=len(od_s), pairs=int(m.sum()), D_rho=d_rho(Z, M, fast, pm_s, dr, m))
        comp = []
        for kind in ("pseudo", "resp"):
            vals = [d_rho(*raw(T, kind, k, R), fast, pm_s, dr[:5], m) for k in range(0, 200, 10)]
            row[f"{kind}_median_D_rho"] = float(np.nanmedian(vals))
        e = excess_C(T, pm_s, dr, R)
        row.update(D_C=e["pm"]["excess"] - e["od"]["excess"], pm_excess=e["pm"]["excess"], od_excess=e["od"]["excess"])
        sens.append(row)
    pd.DataFrame(sens).to_csv(OUT / "sensitivity.csv", index=False)
    print(pd.DataFrame(tests).round(4).T.to_string())
    print(pd.DataFrame(sens).round(3).to_string(index=False))
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
