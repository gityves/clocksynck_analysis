"""A360 (registered): is the person's clock volume set by the body's systemic state, or by a technical scale?

GTEx fast deaths, 25 peripheral tissues (A357 rows). Person volume log r_p: posterior-mean host radius from all tissues
(donors with >= 5 tissues). Regression on base (number of tissues, cos / sin time of death), technical (median RIN,
ischaemic time) and systemic (person GC and IEG states, A359 scores) blocks; partial R^2 of each block against 10,000
block permutations; standardised coefficients with donor-bootstrap 90 % intervals. Demographic block (age, sex, BMI)
descriptive. GTEx time of death and phenotypes are controlled access: only coefficients and R^2 are written.

    python a360.py
Output: ../output/a360_coefficients.csv, a360_tests.csv.
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP.parent))
sys.path.insert(0, str(EXP / "A143_brain_artefact" / "scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a359  # noqa: E402
import clocksync as cs  # noqa: E402
import compute_brain_artefact_checks as A143  # noqa: E402
from clocksync._core import host as HC  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
MIN_TISSUES = 5
BLOCKS = {"base": ["n_tissues", "tod_cos", "tod_sin"], "technical": ["rin", "ischaemic_h"], "systemic": ["gc", "ieg"],
          "demographic": ["age", "male", "bmi"]}


def rss(y, X):
    X = np.column_stack([np.ones(len(y)), X])
    return float(((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2).sum())


def partial_r2(df, block, given, rng, n=10_000):
    y = df.log_r.to_numpy()
    Xg = df[given].to_numpy()
    Xb = df[block].to_numpy()
    r_red = rss(y, Xg)
    pr2 = 1 - rss(y, np.column_stack([Xg, Xb])) / r_red
    null = np.empty(n)
    for k in range(n):
        null[k] = 1 - rss(y, np.column_stack([Xg, Xb[rng.permutation(len(y))]])) / r_red
    return pr2, float((1 + np.sum(null >= pr2)) / (n + 1))


def std_coefs(df, cols):
    Z = (df[cols] - df[cols].mean()) / df[cols].std(ddof=1)
    y = (df.log_r - df.log_r.mean()) / df.log_r.std(ddof=1)
    X = np.column_stack([np.ones(len(y)), Z.to_numpy()])
    return np.linalg.lstsq(X, y.to_numpy(), rcond=None)[0][1:]


def main():
    warnings.filterwarnings("ignore")
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    pb, wf, gc, ieg = a359.load()
    _, donors = gtex.load(min_donors=None, genes=list(wf.genes))
    comps = sorted(pb.compartment.unique())
    d = cs.Dataset(pb[a359.META + list(wf.genes)], anchors=comps, waveform=wf, name="a360").data
    h = HC.host(d, comps)
    per = pd.DataFrame({"log_r": np.log(np.hypot(h.f[:, 0], h.f[:, 1]))}, index=[str(r) for r in h.resections])
    per["n_tissues"] = pd.Series(d.resection.astype(str)).value_counts()
    for lab, genes in (("gc", gc), ("ieg", ieg)):
        _, z = a359.sample_scores(pb, d, genes)
        per[lab] = z.groupby(level=0).mean()
    att = A143.sample_attributes(pb)
    per["rin"] = pd.Series(att.SMRIN.to_numpy(float)).groupby(pb.donor.astype(str).to_numpy()).median()
    dn = donors.copy()
    dn.index = dn.index.astype(str)
    per["tod_cos"], per["tod_sin"] = np.cos(2 * np.pi * dn.tod / 24), np.sin(2 * np.pi * dn.tod / 24)
    per["ischaemic_h"], per["age"], per["male"] = dn.ischaemic_h, dn.age, dn.male
    per["bmi"] = gtex.phenotypes(["BMI"]).BMI.astype(float)
    per = per[per.n_tissues >= MIN_TISSUES]
    core = BLOCKS["base"] + BLOCKS["technical"] + BLOCKS["systemic"]
    df = per.dropna(subset=["log_r"] + core)
    rng = np.random.default_rng(360)
    r2_sys, p_sys = partial_r2(df, BLOCKS["systemic"], BLOCKS["base"] + BLOCKS["technical"], rng)
    r2_tech, p_tech = partial_r2(df, BLOCKS["technical"], BLOCKS["base"] + BLOCKS["systemic"], rng)
    coef = std_coefs(df, core)
    boot = np.array([std_coefs(df.iloc[rng.integers(0, len(df), len(df))], core) for _ in range(2000)])
    co = pd.DataFrame({"term": core, "std_coef": coef, "lo90": np.quantile(boot, 0.05, 0), "hi90": np.quantile(boot, 0.95, 0)})
    dfd = per.dropna(subset=["log_r"] + core + BLOCKS["demographic"])
    r2_dem, p_dem = partial_r2(dfd, BLOCKS["demographic"], core, rng)
    cod = pd.DataFrame({"term": core + BLOCKS["demographic"], "std_coef": std_coefs(dfd, core + BLOCKS["demographic"])})
    co = pd.concat([co.assign(model="registered"), cod.assign(model="with demographic (descriptive)")], ignore_index=True)
    co.to_csv(OUT / "a360_coefficients.csv", index=False)
    tests = pd.DataFrame([
        dict(test="H360.1 systemic", partial_r2=r2_sys, p=p_sys, passes=bool(p_sys < 0.05), donors=len(df)),
        dict(test="H360.2 technical", partial_r2=r2_tech, p=p_tech, passes=bool(p_tech < 0.05), donors=len(df)),
        dict(test="demographic (descriptive)", partial_r2=r2_dem, p=p_dem, donors=len(dfd)),
        dict(test="full model R2", partial_r2=1 - rss(df.log_r.to_numpy(), df[core].to_numpy())
             / float(((df.log_r - df.log_r.mean()) ** 2).sum()), donors=len(df))])
    tests.to_csv(OUT / "a360_tests.csv", index=False)
    print(co.round(3).to_string(index=False))
    print(tests.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
