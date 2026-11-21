"""A328 part a scoring (registered in `paper/docs/registered_experiments.md`, A328, note 1 and power note): blind
replication of the fetal adrenal covariation in del Valle 2023 (32 whole fetal adrenals, 7-11.5 wpc).

Signed coherence C across fetuses (`common.development.l1_test`): log fractions of the exon counts residualised gene by
gene on the stage indicator and sex, 10,000 per-gene permutations within stage; 200 expression-matched and 200
responsive pseudo-clock sets (genes with mean count >= 1), same adjustment; leave-one-gene-out range; unsigned R.
Tests: H328.1 `dv:all` (14 genes), H328.2 `dv:core` (11-gene core waveform), each read on its own; descriptive
`dv:-PER1-NR1D1`, `dv:early` (7-8 wpc), `dv:late` (9-11.5 wpc). Seed 328.

    python score_a328a.py [--jobs 5] [--check]
Outputs (`../output/score_a/`): l1.csv, verdicts.csv (`_check`: 200 permutations, 10 sets).
"""
import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import development as dev  # noqa: E402

SRC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/sccc_mapping/delvalle2023_fetal_adrenal")
PREP = HERE.parents[1] / "output" / "prep"
OUT = HERE.parents[1] / "output" / "score_a"
SEED = 328
WF = cs.gtex()
WAVEFORMS = {"dv:all": WF, "dv:core": WF.drop(["PER1", "PER2", "NR1D1"]), "dv:-PER1-NR1D1": WF.drop(["PER1", "NR1D1"]),
             "dv:early": WF, "dv:late": WF}
PRIMARY = {"dv:all": "H328.1", "dv:core": "H328.2"}
S = {}


def cov(rows) -> np.ndarray:
    return dev.design(rows, age=None, factors=("stage", "sex"))


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fetus x gene exon counts (symbols, duplicates summed) and the sheet of `prep_a328.py`."""
    obs = pd.read_csv(PREP / "dv_sheet.csv")
    obs["stage"] = obs.stage.astype(str)
    exon = pd.read_csv(SRC / "counts_exon.tsv.gz", sep="\t", index_col=0)
    sym = exon.pop("symbol").astype(str).map(canonical).to_numpy()
    counts = exon[obs.run].groupby(sym).sum().T
    counts.index = obs.index
    return counts, obs


def rows(obs, name) -> np.ndarray:
    st = obs.stage.astype(float)
    return {"dv:early": st <= 8, "dv:late": st >= 9}.get(name, pd.Series(True, index=obs.index)).to_numpy()


def task(j):
    name = list(WAVEFORMS)[j]
    obs = S["obs"][rows(S["obs"], name)]
    r = dev.l1_test(S["counts"], obs, cov, WAVEFORMS[name], S["n_perm"], S["n_pseudo"], [SEED, 1, j], "stage")
    return dict(test=name, **r)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=5)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    counts, obs = load()
    S.update(counts=counts, obs=obs, n_perm=200 if a.check else 10000, n_pseudo=10 if a.check else 200)
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        l1 = pd.DataFrame(pool.map(task, range(len(WAVEFORMS))))
    sfx = "_check" if a.check else ""
    l1.to_csv(OUT / f"l1{sfx}.csv", index=False)
    v = l1[l1.test.isin(PRIMARY)].assign(hypothesis=lambda t: t.test.map(PRIMARY))
    v["verdict"] = np.where(v.p >= .05, "FAIL", np.where(v.clock_specific, "PASS", "not clock-specific"))
    v[["hypothesis", "test", "n", "C", "p", "matched_q95", "responsive_q95", "matched_pct", "responsive_pct",
       "verdict"]].to_csv(OUT / f"verdicts{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(l1.round(4).to_string(index=False), v[["hypothesis", "test", "verdict"]].to_string(index=False), sep="\n")


if __name__ == "__main__":
    main()
