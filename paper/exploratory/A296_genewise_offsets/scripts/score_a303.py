"""A303 (registered in `paper/docs/registered_experiments.md`, descriptive with a fixed criterion): is the adrenal offset
from the other organs in mouse and sheep?

Mouse (Wolff 2023, timed on CT): adrenal clock phase minus the circular mean of heart, kidney, lung, muscle, per age and
pooled (age as stratum); libraries resampled within CT (1000). Sheep atlas (untimed): adrenal offset minus the circular
mean of the non-brain tissues from `Dataset.fit`; donors resampled (500). Criterion: |offset| >= 2 h with the 90 %
interval excluding 0.

    python score_a303.py
Output (`../output/a303.csv`).
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import data  # noqa: E402

OUT = HERE.parents[1] / "output"
W = 2 * np.pi / 24
BRAIN_SHEEP = {"Brainstem", "Cerebellum", "Cerebral cortex", "Cerebrum", "Hippocampus", "Hypothalamus",
               "Medulla oblongata", "Optic chiasm", "Pineal"}


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def rel(off: pd.Series, k="adrenal"):
    others = off.drop(k)
    return float(wrap(off[k] - np.angle(np.mean(np.exp(1j * W * others.to_numpy()))) / W))


def mouse_offset(pb, donors):
    comps = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb, anchors=comps, donors=donors, min_donors=10, name="wolff")  # 12 libraries per age
    return rel(ds.timed_loadings("ct").clock_offset_h)


def resample_ct(pb, donors, rng):
    out_pb, out_d = [], []
    for ct, ids in donors.groupby("ct").groups.items():
        pick = rng.choice(np.asarray(ids), len(ids))
        for j, dn in enumerate(pick):
            new = f"{dn}_b{ct}_{j}"  # no "#": timed_loadings maps resections to donors by splitting at "#"
            out_pb.append(pb[pb.donor == dn].assign(donor=new))
            out_d.append(donors.loc[[dn]].set_axis([new]))
    return pd.concat(out_pb, ignore_index=True), pd.concat(out_d)


def main():
    warnings.filterwarnings("ignore")
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(303)
    rows = []
    pb, d = data.load("mouse_ageing_wolff")
    d = d.set_axis(d.index.astype(str))
    pb = pb[pb.compartment != "hypothalamus"].assign(donor=lambda x: x.donor.astype(str))
    for age in ["young", "aged", "old", "pooled"]:
        ids = d.index if age == "pooled" else d.index[d.age == age]
        p, dd = pb[pb.donor.isin(ids)].copy(), d.loc[ids].copy()
        if age == "pooled":
            p["stratum"] = p.donor.map(dd.age)
        est = mouse_offset(p, dd)
        bs = []
        for _ in range(1000):
            try:
                bp, bd = resample_ct(p, dd, rng)
                if age == "pooled":
                    bp["stratum"] = bp.donor.map(bd.age)
                bs.append(mouse_offset(bp, bd))
            except (ValueError, np.linalg.LinAlgError):
                continue
        bs = wrap(np.array(bs) - est) + est
        rows.append(dict(species="mouse", subset=age, offset_h=est, lo90=float(np.quantile(bs, 0.05)),
                         hi90=float(np.quantile(bs, 0.95)), n_boot=len(bs)))
        print(rows[-1], flush=True)
    pb, d = data.load("sheep")
    pb = pb[~pb.compartment.isin(BRAIN_SHEEP)]
    comps = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb, anchors=comps, donors=d, mixed_strata=True, name="sheep")  # as compute_sheep.py
    off = ds.fit(n_boot=0, gate=False).table.offset_h.dropna()
    est = rel(off, "Adrenal")
    bs = []
    for _ in range(500):
        try:
            o = ds.resample(np.random.default_rng(int(rng.integers(2 ** 31)))).fit(n_boot=0, gate=False).table.offset_h.dropna()
            bs.append(rel(o, "Adrenal"))
        except (ValueError, np.linalg.LinAlgError, KeyError):
            continue
    bs = wrap(np.array(bs) - est) + est
    rows.append(dict(species="sheep", subset="all (untimed)", offset_h=est, lo90=float(np.quantile(bs, 0.05)),
                     hi90=float(np.quantile(bs, 0.95)), n_boot=len(bs),
                     median_abs_offset_other_tissues=float(np.median(np.abs(wrap(off.drop("Adrenal") - np.angle(np.mean(np.exp(1j * W * off.to_numpy()))) / W))))))
    t = pd.DataFrame(rows)
    t["same_structure"] = (t.offset_h.abs() >= 2) & ((t.lo90 > 0) | (t.hi90 < 0))
    t.to_csv(OUT / "a303.csv", index=False)
    print(t.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
