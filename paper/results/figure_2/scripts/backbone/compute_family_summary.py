"""A17 healthy organisation by cell family: where does the host-locked clock sit? (stromal-vascular backbone)

Descriptive summary of the A17 re-analysis (docs/PLAN_internal_desync.md (archive), A17 step 8 and the note "A17 estimator vs
anchoring", 2026-10-13). For every gate-passing healthy fit (34 datasets): the point coupling of every compartment under
(i) the learned anchors, (ii) the fixed stromal-vascular / neuronal anchors (frozen tables of `common/anchors_a17.csv`),
(iii) the free host with the A17 corrections (`Dataset(anchors=None)`: each compartment on the host of all the others).
Per family (fibroblast_stromal, endothelial, mural, epithelial, neuron, glia, myeloid, lymphocyte, parenchyma): the median
coupling over compartments, the number of compartments and datasets, and a percentile interval of the median from a
bootstrap over datasets (2,000 draws; a dataset is the cluster). Spearman rank of the family medians across estimators.
Reading (A17 note): stroma, vascular and neurons above myeloid and lymphocytes under every estimator.

    python compute_family_summary.py [--subsample N] [--jobs 4] [--only NAME ...]
Outputs: paper/results/figure_2/output/backbone/family_couplings.csv (one row per healthy fit x compartment),
family_medians.csv (one row per variant x family x estimator; variant `all` = every healthy fit including the
Tabula Sapiens confirmation fit, as the legacy medians, `no_tabula` without it), family_dropped_targets.csv (target, reason: gate failed or error).
Legacy: legacy/desync/common/unified.py (a17_unified_healthy.csv; family medians of the PLAN note).
"""

from __future__ import annotations

import sys
import warnings
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("atlas_fits")))
sys.path.insert(0, str(group_dir("host_identity")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from _shared import FOLDER_4, Target, frozen, registry  # noqa: E402
from anchor_rule import family  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

N_PERM_GATE, N_BOOT_MED = 1000, 2000


def gate_seed(name: str) -> int:
    """Gate seed from a stable hash of the target name (not its registry row, which moves when rows are added)."""
    return zlib.crc32(name.encode()) % 2**31


def anchor_sets(name: str) -> dict:
    a = frozen()
    a = a[a.target == name]
    out = {}
    for rule, col in (("learned", "anchor"), ("fixed", "anchor_fixed")):
        comps = a.compartment[a[col].astype(bool)].tolist()
        ok = len(comps) >= 3 and not (rule == "learned" and str(a.status.iloc[0]).startswith("no anchor"))
        out[rule] = comps if ok else None
    return out


def one(name: str, args, n_perm: int) -> pd.DataFrame:
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        T = Target(name, args.subsample, args.seed)
        free = T.dataset(anchors=None)
        g = free.gate(n_perm=n_perm, seed=gate_seed(name))
        if not g.passes:
            return pd.DataFrame([dict(target=name, reason=f"gate failed (p_pairing {g.p_pairing:.3f}, "
                                                           f"p_acrophase {g.p_acrophase:.3f})")])
        cols = {"free_corrected": free.fit(n_boot=0, gate=g).table.coupling}
        for rule, anc in anchor_sets(name).items():
            anc = [a for a in anc or [] if a in set(T.pb.compartment)]
            if len(anc) >= 3:
                f = T.dataset(anchors=anc).fit(n_boot=0, gate=g).table
                cols[f"coupling_{rule}"], cols[f"anchor_{rule}"] = f.coupling, f.anchor
        t = pd.DataFrame(cols).rename(columns={"free_corrected": "coupling_free_corrected"})
        t = t.rename_axis("compartment").reset_index()
        return t.assign(target=name, family=t.compartment.map(family), tissue=T.tissue, study=T.study)
    except Exception as e:
        print("FAILED", name, type(e).__name__, e, flush=True)
        return pd.DataFrame([dict(target=name, reason=f"error: {type(e).__name__}: {e}")])


def medians(t: pd.DataFrame, rng, variant: str = "all") -> pd.DataFrame:
    rows = []
    for est in ("learned", "fixed", "free_corrected"):
        col = f"coupling_{est}"
        x = t.dropna(subset=[col])
        for fam, g in x.groupby("family"):
            ds = [d[col].to_numpy() for _, d in g.groupby("target")]
            boot = [np.median(np.concatenate([ds[j] for j in rng.integers(len(ds), size=len(ds))]))
                    for _ in range(N_BOOT_MED if len(ds) > 1 else 0)]
            rows.append(dict(variant=variant, estimator=est, family=fam, n_compartments=len(g), n_datasets=len(ds),
                             median=float(np.median(g[col])), lo=np.quantile(boot, 0.025) if boot else np.nan,
                             hi=np.quantile(boot, 0.975) if boot else np.nan))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER_4, run["quick"])
    reg = registry()
    names = [n for n in reg.name if reg.kind[n] == "healthy" and (args.only is None or n in args.only)]
    n_perm = 100 if run["quick"] else N_PERM_GATE
    t = pd.concat(Parallel(n_jobs=args.jobs)(delayed(one)(n, args, n_perm) for n in names), ignore_index=True)
    if "reason" not in t:
        t["reason"] = np.nan
    dropped = t.loc[t.reason.notna(), ["target", "reason"]]
    dropped.to_csv(out / "family_dropped_targets.csv", index=False)
    print(f"{len(dropped)} dropped targets (family_dropped_targets.csv)", flush=True)
    t = t[t.reason.isna()].drop(columns="reason")
    t.to_csv(out / "family_couplings.csv", index=False)
    m = pd.concat([medians(t, np.random.default_rng(args.seed)),
                   medians(t[t.target != "tabula_sapiens"], np.random.default_rng(args.seed), "no_tabula")])
    m.to_csv(out / "family_medians.csv", index=False)
    w = m[m.variant == "all"].pivot(index="family", columns="estimator", values="median")
    print(w.round(2).to_string(), flush=True)
    print("Spearman across families (learned vs free_corrected):",
          round(spearmanr(w.learned, w.free_corrected, nan_policy="omit").statistic, 2), flush=True)


if __name__ == "__main__":
    main()
