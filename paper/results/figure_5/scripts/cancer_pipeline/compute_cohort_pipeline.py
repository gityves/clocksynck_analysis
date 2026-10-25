"""The single-cell cancer cohort pipeline (M4): gate, E1-E3, contrasts H1-H3, own rhythm (E5), gene loadings (E4) and the
microenvironment contrasts, per cohort and pooled over the cohorts that pass the gate.

Pre-registered in docs/cancer/ANALYSIS_PLAN.md (frozen 2026-10-03, amendments to 2026-10-08). Per cohort (patients =
resections, compartments = sample type x class; stratum = 10x chemistry; one global GTEx waveform):
- eligibility: >= 6 patients and >= 3 compartments; reference compartment = the first of tumour fibroblast, endothelial,
  myeloid, T_NK, mural, B_plasma, else the most frequent compartment;
- gate: model-free (`Dataset.gate`, 1,000 permutations); a cohort enters the pooled tests only if both p < 0.05;
- E1-E3: joint factor model with the reference compartment, resection bootstrap (500);
- H1-H3 (`_contrasts`): malignant against every other compartment (log amplitude ratio, offset difference);
  `main` malignant = tumour epithelium, `cnv` = copy-number-called malignant cells (`tumour:malignant`);
- microenvironment (`_micro`): every tumour compartment against the same compartment of adjacent normal tissue, and
  non-malignant tumour epithelium against adjacent-normal epithelium (BH over comparisons);
- E5 own rhythm of every compartment (clock-plane variance of the residual against permuted-acrophase planes, 500);
- E4 gene loadings of the malignant compartment (descriptive: per-gene timing deviation from the conserved waveform).
Pooling (m4_summary): per comparison, DerSimonian-Laird over gated cohorts with SE = bootstrap interval width / 3.92,
pan-cancer and per cancer type; BH-FDR over comparisons within a quantity and family.

    python compute_cohort_pipeline.py [--subsample N] [--variants main cnv] [--studies S ...] [--boot B] [--jobs J]
`--jobs J`: cohorts fitted in parallel (each cohort has its own seed, so the outputs do not depend on J).
`--subsample N`: max(N, 12) random patients per cohort, three small cohorts, fewer resamples.
Outputs (paper/results/figure_5/output/cancer_pipeline/studies/): <study>[_cnv]_{gate,estimates,e5,e4,contrasts,micro}.csv and
summary_{main,cnv}_{gate,e5,e4_genes}.csv, summary_{family}_{contrasts,meta,meta_by_cancer}.csv
Legacy: legacy/scripts/run_study.py, legacy/scripts/m4_summary.py (E1-E5 of the joint model: local, see _shared.py).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer_pipeline"
HOST_PREF = ("fibroblast", "endothelial", "myeloid", "T_NK", "mural", "B_plasma")
MIN_PATIENTS, MIN_COMPARTMENTS = 6, 3
QUICK = ("SRP273482", "GSE132465", "GSE181919")


def micro_pairs(comps) -> list[tuple[str, str]]:
    pairs = [(c, "adjacent_normal:" + c.split(":", 1)[1]) for c in comps if c.startswith("tumour:")]
    pairs += [("tumour:epithelial_nonmalignant", "adjacent_normal:epithelial")]
    return [(k, l) for k, l in pairs if k in comps and l in comps]


def safe_run(study: str, cnv: bool, out: Path, a) -> None:
    try:
        run_study(study, cnv, out, a)
    except Exception as e:  # a failing cohort is reported, the run goes on
        print(study, "cnv" if cnv else "main", "FAILED:", type(e).__name__, e, flush=True)


def run_study(study: str, cnv: bool, out: Path, a) -> None:
    tag = study + ("_cnv" if cnv else "")
    pb, don = data.load("cancer_sc", study=study, cnv=cnv)
    if a.subsample is not None:
        keep = np.random.default_rng(a.seed).permutation(sorted(pb.donor.unique()))[:max(a.subsample, 12)]
        pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
    cancer = str(don.cancer.iloc[0]) if don is not None and len(don) else ""
    comps = pb.groupby("compartment").donor.nunique()
    if pb.donor.nunique() < MIN_PATIENTS or len(comps) < MIN_COMPARTMENTS:
        pd.Series({"status": "ineligible", "n_resections": pb.donor.nunique(), "n_compartments": len(comps)}).to_csv(
            out / f"{tag}_gate.csv")
        print(tag, "ineligible", flush=True)
        return
    ref = next((f"tumour:{h}" for h in HOST_PREF if f"tumour:{h}" in comps.index), comps.idxmax())
    ds = sh.joint_dataset(pb, ref, tag)
    g = ds.gate(n_perm=a.n_gate, seed=a.seed)
    gate = pd.Series({**g.to_series().to_dict(), "p_pairing": g.p_pairing, "p_acrophase": g.p_acrophase,
                      "ref": ref, "cancer": cancer, "waveform": "global"})
    gate.to_csv(out / f"{tag}_gate.csv")
    print(tag, g.summary(), flush=True)
    sh.joint_estimates(ds, a.boot, a.seed).to_csv(out / f"{tag}_estimates.csv")
    k = list(ds.compartments)
    mal = sh.MAL if cnv else "tumour:epithelial"
    if len(k) >= 3:
        pd.DataFrame([sh.own_rhythm(ds, c, a.n_e5, a.seed) for c in k]).to_csv(out / f"{tag}_e5.csv", index=False)
    if mal in k:
        try:
            genes, summ = sh.gene_loadings(ds, mal, a.n_e5, a.seed)
            genes.assign(**summ.to_dict()).to_csv(out / f"{tag}_e4.csv")
        except ValueError as e:  # fewer than two host compartments besides the malignant one
            print(tag, "E4 skipped:", e)
    pairs = [(mal, c) for c in k if c != mal and mal in k]
    if pairs:
        sh.joint_contrasts(ds, pairs, a.boot, a.seed).to_csv(out / f"{tag}_contrasts.csv", index=False)
    mp = micro_pairs(k)
    if mp:
        sh.joint_contrasts(ds, mp, a.boot, a.seed).to_csv(out / f"{tag}_micro.csv", index=False)


# ---- pooling (m4_summary) ----------------------------------------------------------------------------------------

def dl(y, se) -> dict:
    m = cs.stats.meta(y, se)
    z = m["random"] / m["se_random"] if m["k"] else np.nan
    return dict(k=m["k"], estimate=m["random"], lo=m["lo"], hi=m["hi"], p=float(2 * norm.sf(abs(z))), tau2=m["tau2"])


def pool(con: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    rows = []
    for key, g in con.groupby(by):
        key = dict(zip(by, key if isinstance(key, tuple) else (key,)))
        rows.append(dict(**key, quantity="log_amp_ratio", **dl(g.log_amp_ratio, g.se_amp)))
        rows.append(dict(**key, quantity="offset_diff_h", **dl(g.offset_diff_h, g.se_off)))
    out = pd.DataFrame(rows)
    out["q_bh"] = np.nan
    for _, g in out.groupby("quantity"):
        out.loc[g.index, "q_bh"] = cs.stats.bh(g.p.to_numpy())
    return out


def summarise(out: Path, variant: str) -> None:
    name = "cnv" if variant else "main"
    rows = []
    for f in sorted(out.glob(f"*{variant}_gate.csv")):
        tag = f.name[:-len("_gate.csv")]
        if tag.endswith("_cnv") != bool(variant) or tag.startswith("summary"):
            continue
        g = pd.read_csv(f, index_col=0).iloc[:, 0]
        g["study"] = tag.replace("_cnv", "")
        if g.get("status") != "ineligible":
            rows.append(g)
    if not rows:
        return
    gate = pd.DataFrame(rows).set_index("study")
    for c in ("p_pairing", "p_acrophase", "clock_share", "share_null_mean", "n_donors", "n_compartments"):
        gate[c] = pd.to_numeric(gate[c])
    gate["pass"] = (gate.p_pairing < 0.05) & (gate.p_acrophase < 0.05)
    gate.to_csv(out / f"summary_{name}_gate.csv")
    print(name, "gate:", int(gate["pass"].sum()), "of", len(gate), "pass", flush=True)
    passed = list(gate.index[gate["pass"]])
    if not passed:
        return
    e5 = pd.concat([pd.read_csv(out / f"{s}{variant}_e5.csv").assign(study=s) for s in passed
                    if (out / f"{s}{variant}_e5.csv").exists()], ignore_index=True)
    e5[["tissue", "type"]] = e5.compartment.str.split(":", n=1, expand=True)
    e5.groupby(["tissue", "type"]).agg(studies=("p", "size"), own_rhythm=("p", lambda p: int((p < 0.05).sum())),
                                      coupling_clock=("coupling_clock", "median")).reset_index().to_csv(
        out / f"summary_{name}_e5.csv", index=False)
    e4 = [pd.read_csv(out / f"{s}{variant}_e4.csv").assign(study=s) for s in passed if (out / f"{s}{variant}_e4.csv").exists()]
    if e4:
        e4 = pd.concat(e4)
        pd.DataFrame([dict(gene=g, **dl(x.dev_h, x.dev_se_h)) for g, x in e4.groupby("gene", sort=False)]).to_csv(
            out / f"summary_{name}_e4_genes.csv", index=False)
    for suffix, fam in (("contrasts", name), ("micro", f"{name}_micro")):
        con = []
        for s in passed:
            f = out / f"{s}{variant}_{suffix}.csv"
            if f.exists():
                con.append(pd.read_csv(f).assign(study=s, cancer=gate.loc[s].get("cancer", "")))
        if not con:
            continue
        con = pd.concat(con, ignore_index=True).drop_duplicates(["study", "k", "l"])
        con["comparison"] = con.k + " vs " + con.l
        con["se_amp"] = (con.log_amp_ratio_hi - con.log_amp_ratio_lo) / 3.92
        con["se_off"] = (con.offset_diff_hi - con.offset_diff_lo) / 3.92
        con.to_csv(out / f"summary_{fam}_contrasts.csv", index=False)
        meta = pool(con, ["comparison"])
        meta.to_csv(out / f"summary_{fam}_meta.csv", index=False)
        pool(con, ["cancer", "comparison"]).to_csv(out / f"summary_{fam}_meta_by_cancer.csv", index=False)
        print(fam, "\n", meta[["comparison", "quantity", "k", "estimate", "lo", "hi", "p", "q_bh"]].round(3).to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--variants", nargs="+", default=["main", "cnv"], choices=["main", "cnv"])
    ap.add_argument("--studies", nargs="*", help="cohort accessions (default: all; with --subsample: three small ones)")
    ap.add_argument("--boot", type=int, default=None, help="bootstrap replicates (default 500; 20 with --subsample)")
    a = ap.parse_args(argv)
    run = cli.setup(a)
    out = out_dir(FOLDER, run["quick"]) / "studies"
    out.mkdir(exist_ok=True)
    a.boot = a.boot or (20 if run["quick"] else 500)
    a.n_gate, a.n_e5 = (50, 50) if run["quick"] else (1000, 500)
    b = sh.sc_builder()
    for variant in a.variants:
        cnv = variant == "cnv"
        studies = a.studies or (list(QUICK) if run["quick"] else (b.CNV_STUDIES if cnv else b.STUDIES))
        Parallel(n_jobs=a.jobs)(delayed(safe_run)(s, cnv, out, a) for s in studies if not cnv or s in b.CNV_STUDIES)
        summarise(out, "_cnv" if cnv else "")


if __name__ == "__main__":
    main()
