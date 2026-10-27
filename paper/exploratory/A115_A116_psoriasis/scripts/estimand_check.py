"""A115 estimand check (allowed: a known published result). Lesional minus uninvolved core-cue offset, per gated A54
psoriasis series and pooled with the A54 weights; within-group estimator (A115) and paired estimator (A78 / A81b).
Cue sets: A81b GTEx cue (primary), A78 skin systemic set (sensitivity). Healthy samples are never read.

    python -I estimand_check.py
Outputs: output/estimand_check_series.csv, output/estimand_check_pooled.csv, output/gene_presence.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import _lib as L

OUT = Path(__file__).resolve().parents[1] / "output"
sets = L.cue_sets()
genes = L.all_set_genes(sets)
w54 = L.weights54()
FILES = {s: L.ROOT / "processed" / f"{s}.parquet" for s in ["GSE30999", "GSE117468", "GSE201827", "GSE136757", "GSE85034", "GSE252029", "GSE197056", "GSE268958", "GSE121212"]}
FILES |= {s: L.ROOT / "processed_heldout" / f"{s}.parquet" for s in ["GSE13355", "GSE137218", "GSE153007", "GSE193128", "GSE283313", "GSE66511", "GSE69967"]}

rows, comp, pres = [], {}, []
for s, f in FILES.items():
    if s not in w54:
        continue
    d = L.read_series(f, genes)
    d = d[d.disease == "psoriasis"]
    counts = bool((d.value_type == "counts").all())
    pat = sorted(set(d[d.status == "lesional"].patient) & set(d[d.status == "nonlesional"].patient))
    d = d[d.patient.isin(pat)]
    les = d[d.status == "lesional"].drop_duplicates("patient").set_index("patient").loc[pat]
    non = d[d.status == "nonlesional"].drop_duplicates("patient").set_index("patient").loc[pat]
    x = pd.concat([les, non])
    iL = np.r_[np.ones(len(les), bool), np.zeros(len(non), bool)]
    st, n = {}, {}
    for k, wf in {"core": L.WF, "cue": sets["cue"], "systemic": sets["systemic"]}.items():
        st[k], n[k] = L.states_for(x, wf, counts)
    rec = dict(series=s, n_pairs=len(pat), weight=w54[s], core_genes=n["core"], cue_genes=n["cue"], systemic_genes=n["systemic"])
    for k in ("cue", "systemic"):
        c = {g: L.ccorr_sub(*st["core"], *st[k], m) for g, m in (("L", iL), ("U", ~iL))}
        rec[f"within_{k}_L_h"] = np.angle(c["L"]) * L.HOURS
        rec[f"within_{k}_U_h"] = np.angle(c["U"]) * L.HOURS
        rec[f"within_{k}_L_coh"], rec[f"within_{k}_U_coh"] = abs(c["L"]), abs(c["U"])
        rec[f"within_{k}_diff_h"] = L.wrap(rec[f"within_{k}_L_h"] - rec[f"within_{k}_U_h"])
        comp.setdefault(("within", k), {})[s] = c
        pa = {}
        for name in ("core", k):
            b, v = st[name]
            pa[name] = S_c = L.PD.ccorr(b[iL], v[iL], b[~iL], v[~iL])
        rec[f"paired_{k}_core_h"] = np.angle(pa["core"]) * L.HOURS
        rec[f"paired_{k}_set_h"] = np.angle(pa[k]) * L.HOURS
        rec[f"paired_{k}_diff_h"] = L.wrap(rec[f"paired_{k}_core_h"] - rec[f"paired_{k}_set_h"])
        comp.setdefault(("paired", k), {})[s] = pa
    rows.append(rec)
    print(s, len(pat), flush=True)
ser = pd.DataFrame(rows)
ser.to_csv(OUT / "estimand_check_series.csv", index=False)
w = ser.weight.to_numpy(); S_ = ser.series.tolist()
out = []
for k in ("cue", "systemic"):
    cL = L.pooled([comp[("within", k)][s]["L"] for s in S_], w); cU = L.pooled([comp[("within", k)][s]["U"] for s in S_], w)
    out.append(dict(estimator="within-group", set=k, lesional_h=cL, uninvolved_h=cU, diff_h=L.wrap(cL - cU),
                    series=len(S_), pairs=int(ser.n_pairs.sum()), n_series_negative=int((ser[f"within_{k}_diff_h"] < 0).sum())))
    pc = L.pooled([comp[("paired", k)][s]["core"] for s in S_], w); ps = L.pooled([comp[("paired", k)][s][k] for s in S_], w)
    out.append(dict(estimator="paired (A78/A81b)", set=k, lesional_h=pc, uninvolved_h=ps, diff_h=L.wrap(pc - ps),
                    series=len(S_), pairs=int(ser.n_pairs.sum()), n_series_negative=int((ser[f"paired_{k}_diff_h"] < 0).sum())))
for r in out:
    if r["estimator"] == "within-group":
        k = r["set"]
        r["coh_L"] = float(abs(np.sum(w * np.array([comp[("within", k)][s]["L"] for s in S_]))) / w.sum())
        r["coh_U"] = float(abs(np.sum(w * np.array([comp[("within", k)][s]["U"] for s in S_]))) / w.sum())
out = pd.DataFrame(out).rename(columns={"lesional_h": "col1_h", "uninvolved_h": "col2_h"})
out["note"] = np.where(out.estimator == "within-group", "col1 = lesional core-cue, col2 = uninvolved core-cue", "col1 = core L-vs-U offset, col2 = set L-vs-U offset")
out.to_csv(OUT / "estimand_check_pooled.csv", index=False)
pd.set_option("display.width", 250)
print(ser.round(2).to_string(), "\n", out.round(3).to_string())
