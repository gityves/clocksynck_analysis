"""A17 (cancer targets) and A18 part 3: is the coupling of tumour compartments lower than that of adjacent-normal ones?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17 cancer targets; A18 copy-number follow-up). Cohorts: the 21
single-cell cohorts whose gate passes (both p < 0.05, 1,000 permutations) on the tumour / adjacent-normal pseudobulk
table. Host: the frozen A17 anchors of the cohort (`common.anchors`, rule "learned" and "fixed"; at least 3), each anchor
read on the host of the others; coupling c with its donor-bootstrap interval (host refitted, 500).
- A17: tumour minus adjacent-normal coupling per cell class, from the same bootstrap replicates (patients paired, SE =
  std of the replicate difference, percentile interval), pooled over cohorts per class (DerSimonian-Laird).
- A18.3: copy-number-called malignant cells (`tumour:malignant`) minus adjacent-normal epithelium (and minus non-malignant
  tumour epithelium) with the learned anchors, per cohort and pooled (cohorts with calls, no gate refilter).
Coupling scales with host spread: these are within-patient differences on one host fit, read as in the pre-registration.

    python compute_anchored_cancer.py [--subsample N] [--jobs J] [--seed S]
`--subsample N`: three cohorts with max(N, 12) patients, 10 bootstrap replicates.
Outputs (paper/results/figure_5/output/cancer/): a17_cancer.csv, a17_cancer_tumour_minus_normal.csv, a17_cancer_pooled.csv,
a18_cnv.csv, a18_cnv_pooled.csv
Legacy: legacy/desync/common/run_anchored.py (`cancer`), legacy/desync/common/followups.py (`cnv_one`).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.anchors import TABLE, anchors  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer"
QUICK = ("SRP273482", "GSE132465", "GSE181919")
CNV_PAIRS = ((sh.MAL, sh.NORM), ("tumour:epithelial", sh.NORM), (sh.MAL, "tumour:epithelial_nonmalignant"))


def table(study: str, cnv: bool, subsample, seed: int):
    pb, _ = data.load("cancer_sc", study=study, cnv=cnv)
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:max(subsample, 12)]
        pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
    return pb


def anchor_set(study: str, rule: str, present) -> list | None:
    try:
        a = [c for c in anchors(f"cancer_{study}", rule) if c in set(present)]
    except KeyError:
        return None
    return a if len(a) >= 3 else None


def differences(rep: pd.DataFrame, tab: pd.DataFrame, pairs) -> list[dict]:
    """Coupling difference (a minus b) of compartment pairs from the same bootstrap replicates."""
    out = []
    for a, b in pairs:
        if a in rep and b in rep:
            x = (rep[a] - rep[b]).to_numpy()
            x = x[np.isfinite(x)]
            out.append(dict(a=a, b=b, c_a=tab.coupling[a], c_b=tab.coupling[b], diff=tab.coupling[a] - tab.coupling[b],
                            se=float(np.std(x, ddof=1)) if len(x) > 2 else np.nan,
                            lo=np.quantile(x, 0.025) if len(x) else np.nan, hi=np.quantile(x, 0.975) if len(x) else np.nan,
                            n_boot=len(x)))
    return out


def cancer_one(study: str, rule: str, a) -> dict:
    """A17: one cohort, one anchor rule."""
    limit_threads(1)
    target = f"cancer_{study}"
    pb = table(study, False, a.subsample, a.seed)
    anc = anchor_set(study, rule, pb.compartment.unique())
    ds = cs.Dataset(pb, anchors=anc or "auto", name=study)
    g = ds.gate(n_perm=a.n_gate, seed=a.seed)
    head = dict(target=target, rule=rule, gate_p_pairing=g.p_pairing, gate_p_acrophase=g.p_acrophase, gate_passed=g.passes)
    if not g.passes or anc is None:
        return dict(rows=pd.DataFrame([head | dict(status="gate failed" if not g.passes else "no anchor")]))
    f = ds.fit(n_boot=a.n_boot, seed=a.seed, gate=g)
    tab = f.table.reset_index().assign(**head, status="ok")
    tab["site"] = tab.compartment.str.split(":").str[0]
    tab["cls"] = tab.compartment.str.split(":", n=1).str[1]
    diffs = []
    for c in sorted(set(tab.cls)):
        d = differences(f.replicates, f.table, [(f"tumour:{c}", f"adjacent_normal:{c}")])
        diffs += [dict(target=target, rule=rule, cls=c, **{k: v for k, v in r.items() if k not in ("a", "b")}) for r in d]
    return dict(rows=tab, diffs=pd.DataFrame(diffs))


def cnv_one(study: str, a) -> list[dict]:
    """A18.3: one cohort, copy-number-called malignant cells, learned anchors."""
    limit_threads(1)
    pb = table(study, True, a.subsample, a.seed)
    anc = anchor_set(study, "learned", pb.compartment.unique())
    if anc is None:
        return [dict(study=study, status="no anchor")]
    ds = cs.Dataset(pb, anchors=anc, name=study)
    f = ds.fit(n_boot=a.n_boot, seed=a.seed, gate=False)
    return [dict(study=study, status="ok", **r) for r in differences(f.replicates, f.table, CNV_PAIRS)]


def pooled(d: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    rows = []
    for key, g in d.groupby(by):
        m = cs.stats.meta(g["diff"], g.se)
        rows.append(dict(zip(by, key if isinstance(key, tuple) else (key,))) | dict(
            est=m["random"], lo=m["lo"], hi=m["hi"], tau2=m["tau2"], studies=m["k"], negative=int((g["diff"] < 0).sum())))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    a = ap.parse_args(argv)
    run = cli.setup(a)
    out = out_dir(FOLDER, run["quick"])
    a.n_boot, a.n_gate = (10, 50) if run["quick"] else (500, 1000)
    reg = pd.read_csv(TABLE)
    targets = sorted(t.removeprefix("cancer_") for t in reg.target.unique() if t.startswith("cancer_"))
    studies = [s for s in targets if s in QUICK] if run["quick"] else targets
    jobs = [(s, r) for s in studies for r in ("learned", "fixed")]
    res = Parallel(n_jobs=a.jobs)(delayed(cancer_one)(s, r, a) for s, r in jobs)
    rows = pd.concat([r["rows"] for r in res], ignore_index=True)
    rows.to_csv(out / "a17_cancer.csv", index=False)
    diffs = pd.concat([r["diffs"] for r in res if "diffs" in r and len(r["diffs"])], ignore_index=True)
    diffs.to_csv(out / "a17_cancer_tumour_minus_normal.csv", index=False)
    p = pooled(diffs, ["rule", "cls"])
    p.to_csv(out / "a17_cancer_pooled.csv", index=False)
    print(p.round(3).to_string(), flush=True)

    passed = sorted(set(rows[(rows.rule == "learned") & rows.gate_passed.fillna(False).astype(bool)].target.str.removeprefix("cancer_")))
    cnv_studies = [s for s in passed if s in sh.sc_builder().CNV_STUDIES]
    res = Parallel(n_jobs=a.jobs)(delayed(cnv_one)(s, a) for s in cnv_studies)
    c = pd.DataFrame([r for rr in res for r in rr])
    c.to_csv(out / "a18_cnv.csv", index=False)
    ok = c[c.status == "ok"] if len(c) else c
    if len(ok):
        pc = pooled(ok.assign(pair=ok.a + " - " + ok.b), ["pair"])
        pc.to_csv(out / "a18_cnv_pooled.csv", index=False)
        print(pc.round(3).to_string())


if __name__ == "__main__":
    main()
