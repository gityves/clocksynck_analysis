"""A296 (registered in `paper/docs/registered_experiments.md`): are the large organ offsets whole-clock properties or
carried by the signal-induced clock genes (PER1, PER2, NR1D1)?

GTEx fast deaths, A104 periphery rows, timed against the hour of death (`Dataset.timed_loadings("tod")`): each tissue's
offset relative to the mean over tissues, with the full waveform and with PER1, PER2, NR1D1 dropped; donor bootstrap of
the difference (200). H296.1: adrenal |reduced - full| <= 1 h and reduced offset >= +3 h; lung and skin secondary.
Descriptive: per-gene timing (`Dataset.gene_loadings`) for these tissues; baboon ADC, ADM, LUN, SKI (timed on ZT).
Aggregates only (GTEx time of death is controlled access).

    python score_a296.py [--n-boot 200]
Outputs (`../output/`): offsets.csv, tests.csv, genes.csv, baboon.csv.
"""
import argparse
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import data, gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
DROP = ["PER1", "PER2", "NR1D1"]
PRIMARY = "Adrenal Gland"
SECONDARY = ["Lung", "Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"]


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def offsets(ds, hours):
    return ds.timed_loadings(hours).offset_h


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-boot", type=int, default=200)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    pb, donors = gtex.load(min_donors=A104.MIN_DONORS)
    pb = pb[~pb.compartment.str.startswith("Brain")]
    comps = sorted(pb.compartment.unique())
    wf_full, wf_red = cs.gtex(), cs.gtex().drop(DROP)
    ds = {k: cs.Dataset(pb[["donor", "compartment", "depth"] + list(w.genes)], anchors=comps, donors=donors, waveform=w, name=k)
          for k, w in (("full", wf_full), ("reduced", wf_red))}
    off = pd.DataFrame({k: offsets(d, "tod") for k, d in ds.items()})
    off["diff"] = wrap(off.reduced - off.full)
    rng = np.random.default_rng(296)
    boots = []
    for b in range(a.n_boot):
        seed = int(rng.integers(2 ** 31))
        try:
            o = {k: offsets(d.resample(np.random.default_rng(seed)), "tod") for k, d in ds.items()}
            boots.append(wrap(o["reduced"] - o["full"]))
        except (ValueError, np.linalg.LinAlgError):
            continue
    B = pd.concat(boots, axis=1)
    off["diff_lo90"], off["diff_hi90"] = B.quantile(0.05, axis=1), B.quantile(0.95, axis=1)
    off.to_csv(OUT / "offsets.csv")
    rows = []
    for t in [PRIMARY] + SECONDARY:
        r = off.loc[t]
        ok = abs(r["diff"]) <= 1 and (t != PRIMARY or r.reduced >= 3)
        rows.append(dict(tissue=t, test="H296.1" if t == PRIMARY else "secondary", full_h=r.full, reduced_h=r.reduced,
                         diff_h=r["diff"], diff_lo90=r.diff_lo90, diff_hi90=r.diff_hi90, whole_clock=bool(ok)))
    tests = pd.DataFrame(rows)
    tests.to_csv(OUT / "tests.csv", index=False)
    print(off.round(2).sort_values("full").to_string())
    print(tests.round(2).to_string(index=False), flush=True)

    genes = []
    for t in [PRIMARY] + SECONDARY:
        try:
            g, summ = ds["full"].gene_loadings(t, n_perm=200, seed=296)
            genes.append(g.assign(tissue=t, conserved_waveform_p=summ.get("p", np.nan) if isinstance(summ, dict) else np.nan))
        except Exception as e:  # noqa: BLE001
            genes.append(pd.DataFrame([dict(tissue=t, error=f"{type(e).__name__}: {e}"[:200])]))
    pd.concat(genes).to_csv(OUT / "genes.csv")

    pbb, db = data.load("baboon")
    cb = sorted(pbb.compartment.unique())
    bb = {}
    for k, w in (("full", wf_full), ("reduced", wf_red)):
        g = [x for x in w.genes if x in pbb.columns]
        bb[k] = cs.Dataset(pbb[["donor", "compartment", "depth"] + g], anchors=cb, donors=db, waveform=w.subset(g),
                           name=f"baboon_{k}").timed_loadings("zt").offset_h
    bt = pd.DataFrame(bb)
    bt["diff"] = wrap(bt.reduced - bt.full)
    bt.loc[["ADC", "ADM", "LUN", "SKI"]].to_csv(OUT / "baboon.csv")
    print(bt.loc[["ADC", "ADM", "LUN", "SKI"]].round(2).to_string())
    print(f"baboon: median |offset| full {bt.full.abs().median():.2f} h over {len(bt)} tissues", flush=True)


if __name__ == "__main__":
    main()
