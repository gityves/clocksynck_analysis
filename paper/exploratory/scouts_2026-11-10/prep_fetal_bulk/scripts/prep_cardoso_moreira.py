"""Cardoso-Moreira et al. 2019 (Nature 571:505), human bulk RNA-seq time series of seven organs (ArrayExpress E-MTAB-6814): prep ONLY.
NO clock statistic (no correlation, phase, coherence, gate or fit): design tables and detection-level summaries.

Source (untrusted, parse only): sdrf (313 samples, 93 prenatal individuals 4-20 wpc, 120 postnatal rows with individual "not available"),
`Human.RPKM.txt` and `Human.CPM.txt` (Ensembl gene x sample, linear scale, 51,299 genes, 297 samples).

Sample mapping. Matrix columns are named `<Organ>.<stage>.<k>` (Brain, Cerebellum, Heart, Kidney, Liver, Ovary, Testis); they carry no sdrf name.
Brain = forebrain, Cerebellum = hindbrain (sdrf organ counts 55 and 59 equal the maximum k of each organ). k is the rank of the sample within its organ
in the full 313-row table, sorted by age; 16 ranks are absent from the matrix. What is therefore EXACT per column: organ, stage label (matrix), age_wpc (matrix label, integer),
sex (from Y-gene detection, RPS4Y1/DDX3Y/UTY/KDM5D mean RPKM > 1; agrees with the assigned sdrf sex by construction). What is NOT recoverable exactly: which sdrf row
(hence individual) a column is. Each column is given its candidate sdrf rows (same organ, sex, stage within 1 week or same postnatal class); a Hungarian assignment
minimising |k - rank(age, sex, source number)| gives `individual_guess`; `individual` is filled only when all candidates belong to ONE individual
(`individual_certain`), else NA. `n_candidates` and `n_candidate_individuals` are reported.

Library size: CPM alone is not counts; `lib_size_est` = 1e6 / smallest non-zero CPM of the sample (exact if some gene has a single read and CPM is plain, else an estimate).

    python -I prep_cardoso_moreira.py
"""
from __future__ import annotations

import gzip
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from clocksync.genes import CLOCK_GENES, canonical

ROOT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
DAT = ROOT / "natexp/cardoso_moreira2019"
RAW = DAT / "raw"
OUT = Path(__file__).resolve().parents[1] / "output"
ORGAN = {"Brain": "forebrain", "Cerebellum": "hindbrain", "Heart": "heart", "Kidney": "kidney", "Liver": "liver", "Ovary": "ovary", "Testis": "testis"}
POST_CLASS = {"neonate": "newborn", "infant": "infant", "toddler": "toddler", "school age child": "school", "adolescent": "teen",
              "young adult": "youngAdult", "middle adult": "mid", "elderly": "senior"}
MATRIX_CLASS = {"newborn": "newborn", "infant": "infant", "toddler": "toddler", "school": "school", "youngTeenager": "teen", "oldTeenager": "teen",
                "teenager": "teen", "youngAdult": "youngAdult", "youngMidAge": "youngAdult", "olderMidAge": "mid", "senior": "senior", "Senior": "senior"}
Y_GENES = ["ENSG00000129824", "ENSG00000067048", "ENSG00000183878", "ENSG00000012817"]  # RPS4Y1, DDX3Y, UTY, KDM5D
AGE_BINS = [(4, 8), (9, 12), (13, 20)]


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def read_matrix(path):
    return pd.read_csv(path, sep=" ", index_col=0, quotechar='"').astype("float32")


def sdrf_table():
    s = pd.read_csv(RAW / "E-MTAB-6814.sdrf.txt", sep="\t", dtype=str)
    pre = s["Characteristics[developmental stage]"].str.contains("post conception")
    a = s["Characteristics[age]"].astype(float)
    days = a * s["Unit [time unit]"].map({"week": 7, "day": 1, "month": 30.4, "year": 365.25})
    t = pd.DataFrame({
        "source_name": s["Source Name"], "ena_run": s["Comment[ENA_RUN]"], "ena_sample": s["Comment[ENA_SAMPLE]"],
        "organ_matrix": s["Characteristics[organism part]"].map({v: k for k, v in ORGAN.items()}),
        "organ": s["Characteristics[organism part]"], "stage_sdrf": s["Characteristics[developmental stage]"],
        "individual": s["Characteristics[individual]"].where(s["Characteristics[individual]"].ne("not available")),
        "prenatal": pre, "age_wpc": a.where(pre), "age_days_since_conception": days + np.where(pre, 0, 280),
        "sex": s["Characteristics[sex]"], "post_class": s["Characteristics[developmental stage]"].map(POST_CLASS),
    })
    t["srcnum"] = t.source_name.str.extract(r"^(\d+)")[0].astype(int)
    return t


def gene_symbols(ids):
    """Ensembl id -> HGNC symbol: GENCODE v46 table on disk first, GENCODE v29 gene lines for the rest."""
    g = pd.read_csv(ROOT / "reference/gene_positions_gencode46.tsv", sep="\t", header=None, usecols=[0, 1], names=["id", "symbol"]).drop_duplicates("id").set_index("id").symbol
    sym = pd.Series(ids, index=ids).map(g)
    miss = sym.isna()
    if miss.any():
        old = {}
        with gzip.open(ROOT / "gencode_v29/gencode.v29.basic.annotation.gtf.gz", "rt") as f:
            for l in f:
                if l.startswith("#"):
                    continue
                p = l.split("\t", 8)
                if p[2] == "gene":
                    old[p[8].split('gene_id "')[1].split(".")[0]] = p[8].split('gene_name "')[1].split('"')[0]
        sym[miss] = pd.Series(ids, index=ids)[miss].map(old)
    return sym


def assign(m, s):
    """Candidate sdrf rows and Hungarian assignment per organ (see module docstring)."""
    out = []
    for org, g in s.groupby("organ_matrix"):
        g = g.sort_values(["age_days_since_conception", "sex", "srcnum"]).reset_index()
        g["rank"] = np.arange(1, len(g) + 1)
        mm = m[m.organ_matrix == org].reset_index()
        C = np.full((len(mm), len(g)), 1e6)
        for i, r in mm.iterrows():
            ok = g.sex.eq(r.sex)
            ok &= (g.age_wpc - r.age_wpc).abs().le(1) if r.prenatal else g.post_class.eq(r.post_class) | (r.stage == "youngMidAge") & g.post_class.eq("mid")
            C[i, ok.to_numpy()] = np.abs(g["rank"][ok.to_numpy()] - r.k)
        ri, ci = linear_sum_assignment(C)
        for i, j in zip(ri, ci):
            cand = g[C[i] < 1e6]
            inds = cand.individual.fillna("NA_" + cand.source_name)  # unavailable individuals are all distinct
            out.append({"sample": mm["sample"][i], "n_candidates": len(cand), "n_candidate_individuals": inds.nunique(),
                        "individual_certain": bool(inds.nunique() == 1 and cand.individual.notna().all()),
                        "individual_guess": g.individual[j], "sdrf_source_guess": g.source_name[j], "ena_run_guess": g.ena_run[j],
                        "assignment_feasible": bool(C[i, j] < 1e6)})
    return pd.DataFrame(out)


def prep():
    s = sdrf_table()
    rpkm, cpm = read_matrix(RAW / "Human.RPKM.txt"), read_matrix(RAW / "Human.CPM.txt")
    assert list(rpkm.columns) == list(cpm.columns) and (rpkm.index == cpm.index).all()
    say("matrix", rpkm.shape, "sdrf rows", len(s), "distinct runs", s.ena_run.nunique())
    m = pd.DataFrame({"sample": rpkm.columns})
    m[["organ_matrix", "stage", "k"]] = m["sample"].str.split(".", expand=True)
    m["k"] = m.k.astype(int)
    m["organ"] = m.organ_matrix.map(ORGAN)
    m["age_wpc"] = pd.to_numeric(m.stage.str.extract(r"^(\d+)wpc$")[0])
    m["prenatal"] = m.age_wpc.notna()
    m["post_class"] = m.stage.map(MATRIX_CLASS)
    ymean = rpkm.loc[Y_GENES].mean()
    m["sex"] = np.where(ymean.to_numpy() > 1, "male", "female")
    m = m.merge(assign(m, s), on="sample")
    m["individual"] = m.individual_guess.where(m.individual_certain)
    m["age_bin"] = [next((f"{a}-{b}" for a, b in AGE_BINS if a <= w <= b), np.nan) if p else "postnatal" for w, p in zip(m.age_wpc, m.prenatal)]
    nz = cpm.where(cpm > 0).min()
    m["cpm_sum"] = cpm.sum().to_numpy()
    m["lib_size_est"] = (1e6 / nz).to_numpy()
    m["n_detected_rpkm_gt_0"] = (rpkm > 0).sum().to_numpy()
    cols = ["sample", "individual", "individual_certain", "individual_guess", "n_candidates", "n_candidate_individuals", "organ", "organ_matrix", "stage",
            "age_bin", "prenatal", "age_wpc", "post_class", "sex", "sdrf_source_guess", "ena_run_guess", "assignment_feasible", "cpm_sum", "lib_size_est", "n_detected_rpkm_gt_0"]
    m[cols].to_csv(DAT / "samples.csv", index=False)
    m[cols].to_csv(OUT / "cardoso_samples.csv", index=False)
    s.drop(columns="srcnum").to_csv(DAT / "sdrf_rows.csv", index=False)
    say("mapping: columns", len(m), "infeasible assignments", int((~m.assignment_feasible).sum()), "certain individual", int(m.individual_certain.sum()),
        "candidates median", float(m.n_candidates.median()))
    say("sdrf rows without a matrix column (rank gaps):", len(s) - len(m), "; cpm_sum range", float(m.cpm_sum.min()), float(m.cpm_sum.max()),
        "; lib_size_est quantiles", np.quantile(m.lib_size_est, [0, .25, .5, .75, 1]).round(0).tolist())
    # expression + genes
    sym = gene_symbols(list(rpkm.index))
    say("genes", len(sym), "without symbol", int(sym.isna().sum()), "duplicated symbols", int(sym.dropna().duplicated().sum()))
    rpkm.to_parquet(DAT / "expr_rpkm.parquet")
    cpm.to_parquet(DAT / "expr_cpm.parquet")
    pd.DataFrame({"gene_id": sym.index, "symbol": sym.to_numpy()}).to_csv(DAT / "genes.csv", index=False)
    say("value scales: RPKM min/median/max", float(rpkm.values.min()), float(np.median(rpkm.values)), float(rpkm.values.max()),
        "; CPM", float(cpm.values.min()), float(np.median(cpm.values)), float(cpm.values.max()))
    design(m, s)
    detection(m, rpkm, sym)


def design(m, s):
    ct = lambda df, r, c: pd.crosstab(df[r], df[c])
    t = ct(m, "organ", "age_bin")[["4-8", "9-12", "13-20", "postnatal"]]
    t.to_csv(OUT / "design_cardoso_organ_by_agebin.csv")
    print("\nmatrix samples per organ x age bin\n", t.to_string())
    pc = ct(m[~m.prenatal], "organ", "post_class"); pc.to_csv(OUT / "design_cardoso_postnatal_class.csv")
    print("\npostnatal class per organ (matrix labels)\n", pc.to_string())
    print("\nsex x prenatal (expression-inferred)\n", ct(m, "sex", "prenatal").to_string())
    ct(m, "sex", "prenatal").to_csv(OUT / "design_cardoso_sex.csv")
    # individual-level design from the sdrf (all 313 rows; authoritative)
    p = s[s.prenatal & s.individual.notna()]
    n_org = p.groupby("individual").organ.nunique()
    per = n_org.value_counts().sort_index().rename_axis("n_organs").rename("individuals")
    per.to_csv(OUT / "design_cardoso_individuals_per_n_organs_sdrf.csv")
    print(f"\nsdrf prenatal individuals {p.individual.nunique()}, samples {len(p)}, with >=2 organs {int((n_org >= 2).sum())}; per number of organs\n", per.to_string())
    rows = []
    for lab, (a, b) in [("all", (4, 20))] + [(f"{a}-{b}", (a, b)) for a, b in AGE_BINS]:
        q = p[p.age_wpc.between(a, b)]
        pairs = {}
        for _, g in q.groupby("individual"):
            for x, y in __import__("itertools").combinations(sorted(set(g.organ)), 2):
                pairs[(x, y)] = pairs.get((x, y), 0) + 1
        rows += [{"age_bin": lab, "organ_a": x, "organ_b": y, "n_individuals": n} for (x, y), n in pairs.items()]
        rows.append({"age_bin": lab, "organ_a": "ANY_PAIR", "organ_b": "", "n_individuals": int((q.groupby("individual").organ.nunique() >= 2).sum())})
    pr = pd.DataFrame(rows); pr.to_csv(OUT / "design_cardoso_organ_pairs_sdrf.csv", index=False)
    print("\nsdrf same-individual organ pairs (individuals with that pair)\n", pr.pivot_table(index=["organ_a", "organ_b"], columns="age_bin", values="n_individuals", fill_value=0).to_string())
    # the same on the matrix with the (uncertain) individual guess, certain-only count for transparency
    mp = m[m.prenatal]
    print("\nmatrix prenatal columns", len(mp), "with certain individual", int(mp.individual_certain.sum()), "; individuals in guess", mp.individual_guess.nunique())
    g2 = mp.groupby("individual_guess").organ.nunique().value_counts().sort_index()
    g2.to_csv(OUT / "design_cardoso_individuals_per_n_organs_matrix_guess.csv")
    print("matrix guess: individuals per number of organs\n", g2.to_string())


def detection(m, rpkm, sym):
    from clocksync import gtex
    ref = list(gtex().genes)
    print("\nclocksync CLOCK_GENES == gtex().genes:", list(CLOCK_GENES) == ref)
    sy = sym.dropna()
    sy = sy.map(canonical)
    rows = []
    for gname in CLOCK_GENES:
        ids = sy.index[sy == canonical(gname)]
        r = {"gene": gname, "present": len(ids) > 0, "ensembl_ids": ";".join(ids)}
        if len(ids):
            x = rpkm.loc[ids[0]]
            for stage, sel in (("prenatal", m.prenatal), ("postnatal", ~m.prenatal)):
                for org in ORGAN.values():
                    ss = m["sample"][sel & (m.organ == org)]
                    r[f"median_rpkm_{stage}_{org}"] = float(x[ss].median()) if len(ss) else np.nan
        rows.append(r)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "clock_detection_cardoso.csv", index=False)
    d.to_csv(DAT / "clock_detection.csv", index=False)
    print("clock genes present", int(d.present.sum()), "/", len(d), "; missing", d.gene[~d.present].tolist())
    print(d.drop(columns="ensembl_ids").set_index("gene").round(2).to_string())
    q = rpkm.median(axis=1)
    print("\nbackground: all-gene median RPKM per gene, percentiles 10/25/50/75/90:", np.quantile(q, [.1, .25, .5, .75, .9]).round(3).tolist())


if __name__ == "__main__":
    prep()
