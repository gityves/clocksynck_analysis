"""A17 step 8 (contrasts): every disease contrast of A1-A16 on one anchored host, at its original rules and labels,
with the contrast at a common host spread, a specificity companion and the minimum detectable effect.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17 steps 3, 5, 8; A18 addendum 3 for d*; reading rules of A17). Per
contrast (36 targets) and anchor rule (learned, fixed): model-free gate of the pooled table (1,000 permutations; a failed
gate withholds the contrast); one anchored host fitted on the pooled anchors (an anchor on its leave-one-out host);
per compartment c_control, c_case, d = c_case - c_control (donor bootstrap, host fixed), label-permutation p (within the
contrast's blocks, one permutation shared by all compartments), d* and p* at the common host spread, log ratios of V, z2
and psi, phase offset. Holm over the contrast's pre-specified primary compartments (column `primary` of the registry).
Specificity companion (CLAUDE.md rule; same permutations): the target's d minus the mean d of the fixed-rule
stromal-vascular anchors (neurons in brain) of its table. With `--mde`: minimum detectable coupling drop at 80 % power
per estimable compartment of the learned rule (A17 step 5; 100 data sets x 7 damping levels x 200 permutations).

    python compute_contrasts_a17.py [--subsample N] [--seed S] [--jobs 4] [--mde] [--only NAME ...]
Outputs (paper/results/figure_4/output/resilience/): a17_gates.csv, a17_contrast.csv, a17_contrast_companion.csv,
a17_mde.csv (with --mde)
Legacy: legacy/desync/common/run_anchored.py (steps contrast, mde), unified.py, common_v.py.
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
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "resilience"


def _dataset(row, rule, pb, donors, lab, name):
    anc = sh.anchor_set(row.name, rule)
    if anc is None:
        return None, None
    anc = [a for a in anc if a in set(pb.compartment)]
    if len(anc) < 3:
        return None, None
    return cs.Dataset(pb, anchors=anc, donors=donors.assign(**{"lab": lab}), name=name), anc


def one_target(row, run, subsample, seed):
    """Gate, and the contrast under both anchor rules, of one target."""
    limit_threads(1)
    pb, donors, lab, blocks = sh.load_contrast(row, subsample, seed)
    gate_row = dict(target=row.name, donors=lab.size, donors_control=int((lab == "control").sum()),
                    donors_case=int((lab == "case").sum()))
    tabs, comps = [], []
    fixed = sh.anchor_set(row.name, "fixed")
    companion = [a for a in (fixed or []) if a in set(pb.compartment)]
    prim = [p.removesuffix("_coupling") for p in str(row.primary).split(";") if p and p != "nan"]
    for rule in sh.RULES:
        ds, anc = _dataset(row, rule, pb, donors, lab, row.name)
        if ds is None:
            tabs.append(pd.DataFrame([dict(target=row.name, rule=rule, status="no anchor")]))
            continue
        if "p_pairing" not in gate_row:
            g = ds.gate(n_perm=run["n_perm"], seed=seed)
            gate_row.update(p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share,
                            passed=bool(g.passes))
        if not gate_row["passed"]:
            tabs.append(pd.DataFrame([dict(target=row.name, rule=rule, status="gate failed")]))
            continue
        c = ds.contrast("lab", levels=sh.GROUPS, companion=companion or None, blocks=blocks, n_perm=run["n_perm"],
                        n_boot=run["n_boot"], seed=seed)
        t = c.table.reset_index().assign(target=row.name, rule=rule, anchors=";".join(anc))
        t["primary"] = t.compartment.str.split(":").str[-1].isin(prim)
        tabs.append(t)
        if c.companion is not None and rule == "learned":
            comps.append(c.companion.assign(target_contrast=row.name))
    return gate_row, pd.concat(tabs, ignore_index=True), (pd.concat(comps) if comps else None)


def mde_one(row, k, n_mde, subsample, seed):
    """MDE of compartment k (learned anchors): `Dataset.contrast(mde=True)` restricted to one target."""
    limit_threads(1)
    pb, donors, lab, blocks = sh.load_contrast(row, subsample, seed)
    ds, _ = _dataset(row, "learned", pb, donors, lab, row.name)
    try:
        c = ds.contrast("lab", levels=sh.GROUPS, targets=[k], blocks=blocks, n_perm=0, n_boot=0, seed=seed,
                        mde=True, n_mde=n_mde)
        m = c.mde.iloc[0]
        return dict(target=row.name, compartment=k, mde=m.get("mde", np.nan), c_pooled=m.get("c_pooled", np.nan),
                    curve=m["curve"].to_json(orient="records") if "curve" in m and m["curve"] is not None else "")
    except Exception as e:  # reported in the table; the run goes on
        return dict(target=row.name, compartment=k, mde=np.nan, error=f"{type(e).__name__}: {e}")


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--mde", action="store_true", help="also the minimum detectable effect (slow)")
    ap.add_argument("--only", nargs="+", default=None, help="targets (default: all 36 contrasts)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    tg = sh.targets("contrast")
    if args.only:
        tg = tg[tg.name.isin(args.only)]
    rows = list(tg.itertuples())
    res = Parallel(n_jobs=args.jobs)(delayed(one_target)(r, run, args.subsample, args.seed) for r in rows)
    gates = pd.DataFrame([g for g, _, _ in res])
    tab = pd.concat([t for _, t, _ in res], ignore_index=True)
    ok = tab.p.notna() if "p" in tab else pd.Series(False, index=tab.index)
    tab["p_holm_primary"] = np.nan
    for rule in sh.RULES:
        m = (tab.rule == rule) & tab.primary.fillna(False).astype(bool) & ok
        tab.loc[m, "p_holm_primary"] = sh.holm_within(tab[m], "p")
    comp = pd.concat([c for _, _, c in res if c is not None], ignore_index=True) if any(c is not None for *_, c in res) \
        else pd.DataFrame()
    gates.to_csv(out / "a17_gates.csv", index=False)
    tab.to_csv(out / "a17_contrast.csv", index=False)
    comp.to_csv(out / "a17_contrast_companion.csv", index=False)
    print(gates.round(3).to_string(index=False))
    print(tab[tab.rule == "learned"].groupby("target").status.value_counts().unstack(fill_value=0).to_string())
    if args.mde:
        by = {r.name: r for r in rows}
        todo = tab[(tab.rule == "learned") & (tab.status == "ok") & tab.p.notna()][["target", "compartment"]]
        n_mde = 10 if run["quick"] else 100
        if run["quick"]:
            todo = todo.groupby("target").head(2)
        r = Parallel(n_jobs=args.jobs)(delayed(mde_one)(by[t], k, n_mde, args.subsample, args.seed)
                                       for t, k in todo.itertuples(index=False))
        pd.DataFrame(r).to_csv(out / "a17_mde.csv", index=False)


if __name__ == "__main__":
    main()
