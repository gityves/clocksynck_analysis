"""A85: which cell types follow the timing cues rather than the local clock? (genome-wide Census atlases)

Pre-registered in docs/PLAN_internal_desync.md (archive) (A85, H85). Data: 9 atlases (`data/cue_units`; healthy donors; unit = donor x
tissue, majority assay; coarse classes; pseudobulks >= 20 cells; classes with >= 20 units). Per unit and class k: c_k = the core
clock state (GTEx waveform), H_-k = the mean core state of the other classes (reference-free host), Q_-k = the cue state (A81b
GTEx cue set, `figure_S07/a81b_cue_reference.csv`) of the summed counts of the other classes (>= 2 other classes); Poisson
variances as in `paired.states`; each state centred within class x stratum and scaled to unit variance. Complex least squares
c_k = a H_-k + b Q_-k; cue-following index |b| / (|a| + |b|). Donor-unit bootstrap (500, same draws for every class of an atlas).
H85: mean index of immune classes (Treg, CD4_T, CD8_T, T_other, NK, B, plasma, monocyte, macrophage_DC, mast) minus structural
ones (endothelial, mural_smooth_muscle, fibroblast_stromal) > 0, pooled over classes with inverse-variance weights (bootstrap SE),
one-sided bootstrap p < 0.05. Secondary (descriptive): corr(c_k, H_-k) vs corr(c_k, Q_-k) and corr(H_-k, Q_-k) per class.
Caveat (as registered): H and Q are correlated, so b is identified only through units whose clock and cues differ; a null is
uninformative unless the bootstrap SE is small.

    python compute_cue_following.py [--subsample N] [--seed S]     (N: units per atlas, 50 bootstrap draws)
Outputs (paper/results/figure_5/output/cancer/): a85_classes.csv, a85_verdict.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a85_cue_following.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "cancer"
W = 2 * np.pi / 24
IMMUNE = {"Treg", "CD4_T", "CD8_T", "T_other", "NK", "B", "plasma", "monocyte", "macrophage_DC", "mast"}
STRUCTURAL = {"endothelial", "mural_smooth_muscle", "fibroblast_stromal"}
N_BOOT, SEED = 500, 85
MIN_UNITS = 20


sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_out  # noqa: E402
def waveforms() -> dict:
    """Core (GTEx 14 genes) and cue sets (A81b cue reference of figure_S07, GTEx time-of-death frame)."""
    ref = pd.read_csv(group_out("gtex_cues") / "a81b_cue_reference.csv", index_col=0)
    x = ref[ref.set == "cue"]
    cue = cs.Waveform(tuple(canonical(g) or g for g in x.index), x.amp.to_numpy(float), x.acro_h.to_numpy(float) * W, None, None, "cue")
    return {"core": cs.gtex(), "cue": cue}


def state(x: pd.DataFrame, wf: cs.Waveform) -> np.ndarray:
    g = [y for y in wf.genes if y in x.columns]
    wf = wf.subset(g) if wf.name != "cue" else cs.Waveform(tuple(g), np.asarray(wf.amp)[[wf.genes.index(y) for y in g]],
                                                          np.asarray(wf.psi)[[wf.genes.index(y) for y in g]], None, None, "cue")
    L = projector(x, wf, counts=True)
    b, _ = states(x, L, wf, depth=x.depth.to_numpy(float))
    return b[:, 0] + 1j * b[:, 1]


def centre(z, strata):
    z = pd.Series(z)
    z = z - z.groupby(np.asarray(strata)).transform("mean")
    z = z.to_numpy()
    return z / np.sqrt(np.mean(np.abs(z) ** 2))


def design(t: pd.DataFrame, wfs, cue_g) -> dict:
    """Per class: aligned (c_k, H_-k, Q_-k) over units with >= 2 other classes."""
    t = t.assign(c=state(t, wfs["core"]))
    t["c"] = t.groupby("cls", group_keys=False).apply(lambda x: pd.Series(centre(x.c.to_numpy(), x.stratum), index=x.index))
    out = {}
    for k in t.cls.unique():
        own = t[t.cls == k].set_index("unit")
        rest = t[(t.cls != k) & t.unit.isin(own.index)]
        cnt = rest.groupby("unit").size()
        units = cnt[cnt >= 2].index
        if len(units) < MIN_UNITS:
            continue
        H = rest.groupby("unit").c.mean().loc[units]
        summed = rest.groupby("unit")[cue_g + ["depth"]].sum().loc[units]
        Q = centre(state(summed, wfs["cue"]), own.loc[units, "stratum"])
        out[k] = (own.loc[units, "c"].to_numpy(), centre(H.to_numpy(), own.loc[units, "stratum"]), Q)
    return out


def fit(c, H, Q):
    a, b = np.linalg.lstsq(np.column_stack([H, Q]), c, rcond=None)[0]
    return abs(b) / (abs(a) + abs(b))


def corr(x, y):
    return abs(np.vdot(y, x)) / np.sqrt(np.vdot(x, x).real * np.vdot(y, y).real)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = 50 if run["quick"] else N_BOOT
    cu = data._builder("cue_units")
    wfs = waveforms()
    core_g = list(wfs["core"].genes)
    cue_g = [g for g in wfs["cue"].genes if g not in core_g]
    genes = list(dict.fromkeys(core_g + list(wfs["cue"].genes)))
    rng = np.random.default_rng(SEED + args.seed)
    rows, boots = [], {}
    for name in (["lung_hlca", "heart_global"] if run["quick"] else cu.ATLASES):
        t = cu.table(name, genes)
        cue_present = [g for g in wfs["cue"].genes if g in t.columns]
        if run["quick"]:
            keep = np.random.default_rng(args.seed).permutation(t.unit.unique())[:args.subsample or 150]
            t = t[t.unit.isin(keep)]
        des = design(t, wfs, cue_present)
        draws = [rng.integers(0, 10 ** 9) for _ in range(n_boot)]
        for k, (c, H, Q) in des.items():
            n = len(c)
            bs = np.array([fit(c[i], H[i], Q[i]) for i in (np.random.default_rng(s).integers(0, n, n) for s in draws)])
            boots[(name, k)] = bs
            rows.append(dict(atlas=name, cls=k, units=n, index=fit(c, H, Q), se=float(bs.std()), corr_host=corr(c, H),
                             corr_cue=corr(c, Q), corr_host_cue=corr(H, Q),
                             group="immune" if k in IMMUNE else "structural" if k in STRUCTURAL else "other"))
        print(name, len(des), "classes", flush=True)
    d = pd.DataFrame(rows)
    w = 1 / np.maximum(d.se, 1e-3) ** 2
    pooled = lambda vals: {g: np.sum(w[d.group == g] * vals[d.group == g]) / np.sum(w[d.group == g]) for g in ("immune", "structural")}
    obs = pooled(d["index"].to_numpy())
    B = np.array([boots[(a, k)] for a, k in zip(d.atlas, d.cls)])
    diffs = np.array([(lambda p: p["immune"] - p["structural"])(pooled(B[:, r])) for r in range(n_boot)])
    v = dict(immune=obs["immune"], structural=obs["structural"], diff=obs["immune"] - obs["structural"],
             p_upper=float((1 + np.sum(diffs <= 0)) / (1 + n_boot)), lo=float(np.quantile(diffs, 0.025)),
             hi=float(np.quantile(diffs, 0.975)), classes=len(d), atlases=d.atlas.nunique())
    v["H85"] = "PASS" if (v["diff"] > 0 and v["p_upper"] < 0.05) else "FAIL"
    d.to_csv(out / "a85_classes.csv", index=False)
    pd.Series(v).to_csv(out / "a85_verdict.csv", header=False)
    print(d.groupby("group")[["index", "corr_host", "corr_cue"]].median().round(3).to_string(), "\n", pd.Series(v).to_string())


if __name__ == "__main__":
    main()
