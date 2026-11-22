"""A343 identity check and power stage (registered in `paper/docs/registered_experiments.md`, A343), before any GR or
clock-gene count of GSE320593 is read.

Identity: per run, log CPM of XIST and of the mean of RPS4Y1, DDX3Y, UTY, KDM5D (exon counts); sex = male if the
Y mean exceeds XIST. A patient code is one person if all its runs (both tissues, both biopsies) agree.

Power (registered rule): the primary units are patients whose non-surgical (N) and surgical (S) runs share
strandedness. Per run and gene, r = log((u + 0.5) / (s + 0.5)) (u = `counts_intron`, s = `counts_exon`). A random
four-gene set stands in for the GR set: genes outside the clock panel, `CLOCK_LIKE`, the responsive and broad
responsive families and the GR set itself, with mean N-run exon CPM >= 10 and intron-only length >= 1 kb. Per run,
R = mean r over the set minus the mean r of 50 control genes matched to each set gene (same exon-expression and
intron-length quantile bins, same exclusions); D = R(S) - R(N) per patient. Over 1,000 sets: SD of the patient-mean D
(the companion bar) and the median sign-flip standard error; MDE = 2.49 x the larger. Decidable if MDE <= 0.4.
H343.3: 200 pseudo core-loop waveforms (11 random non-target genes with mean N-run CPM >= 10 carrying the core
waveform; not matched to the clock genes, whose levels stay unread until scoring); per patient the S - N
displacement of the pseudo state; the 90 % bootstrap upper bound of |mean(DEX) - mean(NO_DEX)| (patients resampled
within arm); reachable if its median over sets is < 0.25 (state units: the GTEx adult daily amplitude is 1).

    python power_a343.py
Outputs (`../output/power/`): identity.csv, identity_runs.csv, power_gr.csv, power_core.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync import waveform as W  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE, RESPONSIVE_GENES  # noqa: E402
from clocksync._core.responsive_broad import RESPONSIVE_BROAD  # noqa: E402
from clocksync.genes import CLOCK_GENES  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/sccc_mapping/gse320593_liver_fat_dex")
OUT = HERE.parents[1] / "output" / "power"
SEED = 3430
GR = ("FKBP5", "TSC22D3", "ZBTB16", "KLF9")
TARGETS = set(GR) | {"PER1"} | set(CLOCK_GENES) | CLOCK_LIKE | set(RESPONSIVE_GENES) | set(RESPONSIVE_BROAD)
CORE_DROP = ["PER1", "PER2", "NR1D1"]
Y = ("RPS4Y1", "DDX3Y", "UTY", "KDM5D")
TISSUE = {"L": "Liver", "A": "Adipose - Subcutaneous"}
N_SETS, N_CTRL, N_CORE, N_BOOT, N_FLIP = 1000, 50, 200, 2000, 2000


def samples() -> pd.DataFrame:
    s = pd.read_csv(DATA / "samples.tsv", sep="\t")
    t = s.geo_title.str.extract(r"^([LA])(\d+)_.*_(surgical|non-surgical)_(NO_DEX|DEX)$")
    t.columns = ["tis", "pid", "biopsy", "dex"]
    t["biopsy"] = t.biopsy.map({"surgical": "S", "non-surgical": "N"})
    t["run"], t["strand"] = s.run_accession, s.strandedness_used.str[:3]
    assert t.notna().all().all()
    return t


def load(name: str) -> pd.DataFrame:
    x = pd.read_csv(DATA / f"counts_{name}.tsv.gz", sep="\t")
    x = x.drop_duplicates("symbol").set_index("symbol").drop(columns="gene_id")
    return x


def identity(ex: pd.DataFrame, t: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sub = ex.loc[["XIST", *Y], t.run]
    cpm = sub / ex[t.run].sum() * 1e6
    r = t.copy()
    r["log_xist"] = np.log1p(cpm.loc["XIST"].to_numpy())
    r["log_y"] = np.log1p(cpm.loc[list(Y)].mean().to_numpy())
    r["sex"] = np.where(r.log_y > r.log_xist, "M", "F")
    p = r.groupby("pid").agg(runs=("run", "size"), tissues=("tis", "nunique"), sexes=("sex", "nunique"),
                             sex=("sex", "first")).reset_index()
    p["consistent"] = p.sexes == 1
    return r, p


def primary(t: pd.DataFrame, tis: str) -> pd.DataFrame:
    """One row per patient of tissue `tis` with N and S runs, flag for shared strandedness."""
    u = t[t.tis == tis].pivot(index=["pid", "dex"], columns="biopsy", values=["run", "strand"])
    u.columns = [f"{a}_{b}" for a, b in u.columns]
    u = u.reset_index().dropna()
    u["same"] = u.strand_N == u.strand_S
    return u


def bins(v: pd.Series, k: int = 10) -> pd.Series:
    return pd.Series(pd.qcut(v.rank(method="first"), k, labels=False), index=v.index)


def gr_power(ex, it, ann, t, tis, rng) -> list[dict]:
    pts = primary(t, tis)
    runs_n = pts.run_N.tolist()
    cpm = ex[runs_n].div(ex[runs_n].sum(), axis=1) * 1e6
    mean_cpm = cpm.mean(1)
    pool = [g for g in ex.index if g not in TARGETS and g in ann.index and mean_cpm.get(g, 0) >= 10
            and ann.loc[g, "intron_only_len_bp_unstranded"] >= 1000]
    eb = bins(np.log(mean_cpm[pool]))
    lb = bins(np.log(ann.loc[pool, "intron_only_len_bp_unstranded"].astype(float)))
    key = eb.astype(str) + "_" + lb.astype(str)
    by_key = {k: v.index.to_numpy() for k, v in key.groupby(key)}
    allruns = pts.run_N.tolist() + pts.run_S.tolist()
    r = np.log((it.loc[pool, allruns] + 0.5) / (ex.loc[pool, allruns] + 0.5))
    rows = []
    for arm in ("DEX", "NO_DEX"):
        for units in ("primary", "all"):
            p = pts[(pts.dex == arm) & (pts.same | (units == "all"))]
            if len(p) < 3:
                continue
            means, ses = [], []
            for _ in range(N_SETS):
                genes = rng.choice(pool, 4, replace=False)
                ctrl = np.concatenate([rng.choice(by_key[key[g]], N_CTRL) for g in genes])

                def R(run):
                    return r.loc[genes, run].mean() - r.loc[ctrl, run].mean()

                d = (R(p.run_S.tolist()).to_numpy() - R(p.run_N.tolist()).to_numpy())
                means.append(d.mean())
                flips = rng.choice([-1, 1], (N_FLIP, len(d))) @ d / len(d)
                ses.append(flips.std())
            sd, se = float(np.std(means)), float(np.median(ses))
            rows.append(dict(tissue=tis, arm=arm, units=units, patients=len(p), sd_sets=sd, se_signflip=se,
                             mde=2.49 * max(sd, se), decidable=2.49 * max(sd, se) <= 0.4))
    return rows


def core_power(ex, t, tis, rng) -> list[dict]:
    pts = primary(t, tis)
    wf = W.gtex(TISSUE[tis]).drop(CORE_DROP)
    runs = pts.run_N.tolist() + pts.run_S.tolist()
    depth = ex[runs].sum()
    cpm = ex[pts.run_N.tolist()].div(depth[pts.run_N.tolist()], axis=1).mean(1) * 1e6
    pool = np.array([g for g in ex.index if g not in TARGETS and cpm.get(g, 0) >= 10])
    rows = []
    for units in ("primary", "all"):
        p = pts[pts.same | (units == "all")]
        dex = (p.dex == "DEX").to_numpy()
        if dex.sum() < 3 or (~dex).sum() < 3:
            rows.append(dict(tissue=tis, units=units, dex=int(dex.sum()), no_dex=int((~dex).sum()), ub90_median=np.nan))
            continue
        ubs = []
        for j in range(N_CORE):
            pw = Waveform(tuple(rng.choice(pool, len(wf.genes), replace=False)), wf.amp, wf.psi, None, None, f"p{j}")
            g = list(pw.genes)
            x = ex.loc[g, runs].T
            L = projector(x, pw)
            b, _ = states(x, L, pw, depth=depth[runs].to_numpy())
            n = len(pts)
            disp = pd.DataFrame(b[n:] - b[:n], index=pts.index).loc[p.index].to_numpy()
            z = disp[:, 0] + 1j * disp[:, 1]
            zd, zn = z[dex], z[~dex]
            bs = [abs(rng.choice(zd, len(zd)).mean() - rng.choice(zn, len(zn)).mean()) for _ in range(500)]
            ubs.append(np.quantile(bs, 0.95))
        rows.append(dict(tissue=tis, units=units, dex=int(dex.sum()), no_dex=int((~dex).sum()),
                         ub90_median=float(np.median(ubs)), reachable=bool(np.median(ubs) < 0.25)))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t = samples()
    ex, it = load("exon"), load("intron")
    ann = pd.read_csv(DATA / "gene_annotation.tsv.gz", sep="\t").drop_duplicates("symbol").set_index("symbol")
    r, p = identity(ex, t)
    r.to_csv(OUT / "identity_runs.csv", index=False)
    p.to_csv(OUT / "identity.csv", index=False)
    two = p[p.tissues == 2]
    print(f"identity: {int(p.consistent.sum())}/{len(p)} codes consistent; two-tissue codes {len(two)}, "
          f"consistent {int(two.consistent.sum())}", flush=True)
    gr = pd.DataFrame([row for tis in ("A", "L") for row in gr_power(ex, it, ann, t, tis, rng)])
    gr.to_csv(OUT / "power_gr.csv", index=False)
    core = pd.DataFrame([row for tis in ("A", "L") for row in core_power(ex, t, tis, rng)])
    core.to_csv(OUT / "power_core.csv", index=False)
    pd.set_option("display.width", 200)
    print(gr.round(3).to_string(index=False), core.round(3).to_string(index=False), sep="\n")


if __name__ == "__main__":
    main()
