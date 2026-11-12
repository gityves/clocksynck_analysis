"""A247: do glioma cells share a clock state with the patient's own non-malignant brain cells, and does MGMT follow it?

Registered in paper/docs/registered_experiments.md (A247). Extended GBmap droplet pseudobulks (`prep_a247.py`): patient x class
(malignant, myeloid, oligodendrocyte = mature only, T_NK), >= 100 cells, stratum = study | suspension (called "study" below; a
study enters a test with >= 4 patients holding the pseudobulks the test needs).

State z: complex clock state of a pseudobulk (`clocksync.paired.projector / states`, GTEx waveform, counts), projector fitted on all
rows (>= 100 cells) of that class in the stratum, centred within stratum over the patients of the test.
Tests (kind): P247.0 myeloid vs oligodendrocyte, H247.1 malignant vs myeloid, H247.2 malignant vs oligodendrocyte, descriptive T_NK vs
myeloid (pair: |complex correlation| of the two states, A227 `raw_stat`, pooled with patient-count weights; null: rows of the second class
permuted within stratum); H247.3 y = log2(CPM + 1) of MGMT in the malignant pseudobulk (>= 10 MGMT UMI) against the myeloid state,
H247.4a the same against the malignant state (assoc: A229 `stat`, null: y permuted within stratum); H247.4b methylation call against the
myeloid state, T = |mean z(methylated) - mean z(unmethylated)| per stratum holding both calls, pooled with patient-count weights, labels
permuted within stratum (label).
Clock specificity: 200 expression-matched and 200 responsive pseudo-clock sets (`clocksync.paired.pseudo_waveforms`), the same statistic
pooled per set; clock-specific = above the q95 of both. `:corrected`: the expected foreign clock-gene counts (marker share a of the other
class, A241 `contamination`) are subtracted from both classes of a pair and the pairing test is rerun.

  power        eligibility, the null of every test from SHUFFLED rows only (the identity pairing, the unshuffled readout or the unshuffled
               labels are never evaluated; no observed value is computed, stored or printed), MDE = 2.49 x SD, scored = >= 40 patients and
               MDE <= 0.25; marker-based foreign-RNA share a per pair and stratum.   -> eligibility.csv, power.csv, foreign_rna.csv
  calibration  healthy-brain pass rate of P247.0 (microglia vs oligodendrocyte of control donors of a cortex atlas, thinned to the GBmap
               clock-gene counts of the P247.0 patients), false-positive rate with mis-paired donors; observed statistics on the healthy
               atlas only.                                                            -> calibration.csv
  tests        the registered statistics, Holm, companions, `:corrected`, descriptives (IMPLEMENTED; run on real pairings only after the
               power note is committed). -> tests.csv, pooled.csv, per_study.csv, pseudo.csv, matched_genes.csv, descriptive.csv,
               foreign_clock_share.csv
  --placebo    code test: before anything else the rows of the second class (pairs), the y vector or the labels are shuffled once within
               stratum (seed 247) and everything is written to output/placebo/.

Choices not fixed by the registration:
- Projector per class and study (stratum), fitted on all rows of the class in the stratum, not only on the patients of the test.
  (A227 fitted one projector on both classes of a cohort; the |complex correlation| is invariant to the frame, the gene weights differ.)
- Pseudo-clock sets: one set of 200 per test and kind, drawn once (reference = the second class's pseudobulks of the test patients, genes
  with mean count >= 1 or waveform genes) and applied to every study; pooled per set with the same patient weights (A227 drew per cohort).
- Responsive sets repeat genes only if the responsive pool is smaller than the waveform (`pseudo_waveforms` warns).
- H247.4b pooling: patient-count weights over strata with >= 4 patients and >= 1 of each call; modulus per stratum (strata have their own
  frame). Its MDE is on the scale of z (not a correlation): the rule MDE <= 0.25 is applied literally and `mde_over_z_rms` gives it relative
  to the root-mean-square spread of the states.
- Foreign RNA: a = (marker counts of class y in the class-x pseudobulk / depth_x) / (the same in the class-y pseudobulk), clipped to [0, 1],
  per patient with both pseudobulks; a not defined (no marker count in y) is set to 0 in the correction. Both classes of a pair are
  corrected, each for the other's RNA, from the uncorrected counts. T_NK has no registered markers: T_NK vs myeloid is not corrected.
- Study unit of the pooled statistic = stratum; descriptives on suspension (cells / nuclei) and stage (primary / recurrent incl. neo-aPD1;
  unknown excluded) restrict the patients and pool the strata that keep >= 4.
- Calibration: control donors (`healthy` of the atlas builder) with both classes; one atlas at a time; thinning by binomial draws with
  p = min(1, target / count) per gene, depth scaled by the thinned share of the panel; no matched pseudo sets (the atlases on disk carry a
  45-gene panel, not the genome).

    python compute_a247.py --stage power|calibration|tests [--jobs J<=3] [--n-perm 10000] [--n-pseudo 200] [--placebo] [--smoke] [--out DIR]
"""
from __future__ import annotations

import importlib.util
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
for p in ("A227_tme_compartment_pairing", "A229_endothelial_gate", "A149_paired_tumour_normal", "A151_fibroblast_pairing",
          "A226_stroma_pairing_confirm"):
    sys.path.insert(0, str(PAPER / "exploratory" / p / "scripts"))
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_a227 as A  # noqa: E402
import compute_a229 as R  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import projector, pseudo_waveforms, states  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
SEED = 247
WF = cs.gtex()
CLOCK = list(WF.genes)
MIN_N, MIN_PATIENTS, MDE_K, MAX_MDE, MIN_STUDY = 4, 40, A.MDE_K, A.MAX_MDE, 12
MIN_MGMT = 10
MARKERS = {"myeloid": ["PTPRC", "CSF1R", "C1QA", "C1QB", "TYROBP"], "oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP"],
           "malignant": ["SOX2", "PTPRZ1", "EGFR"]}
# test -> (hypothesis, role, kind, first class / state class, second class, family)
TESTS = {
    "P247.0": ("P247.0", "precondition", "pair", "myeloid", "oligodendrocyte", None),
    "H247.1": ("H247.1", "primary", "pair", "malignant", "myeloid", "pairing"),
    "H247.2": ("H247.2", "primary", "pair", "malignant", "oligodendrocyte", "pairing"),
    "H247.3": ("H247.3", "primary", "assoc", "malignant", "myeloid", None),
    "H247.4a": ("H247.4a", "secondary", "assoc", "malignant", "malignant", "h4"),
    "H247.4b": ("H247.4b", "secondary", "label", "malignant", "myeloid", "h4"),
    "TNK_myeloid": (None, "descriptive", "pair", "T_NK", "myeloid", None),
}
NAMES = list(TESTS)
FOREIGN = {"P247.0": [("oligodendrocyte", "myeloid"), ("myeloid", "oligodendrocyte")],
           "H247.1": [("malignant", "myeloid"), ("myeloid", "malignant")],
           "H247.2": [("malignant", "oligodendrocyte"), ("oligodendrocyte", "malignant")],
           "TNK_myeloid": [("T_NK", "myeloid")]}  # (contaminated class, source class)
G: dict = {}  # shared with the workers (fork)


# ---- tables ----------------------------------------------------------------------------------------------------------

class Tab:
    """Pseudobulk table indexed for lookups: `rows(cls, stratum)` = DataFrame (index donor; depth, counts)."""

    def __init__(self, df: pd.DataFrame):
        self.df = df
        self.cache = {k: g.set_index("donor") for k, g in df.groupby(["compartment", "stratum"], sort=True)}

    def rows(self, cls, stratum):
        return self.cache.get((cls, stratum))

    def strata(self, cls):
        return sorted(s for c, s in self.cache if c == cls)


def load_table(path=OUT / "pseudobulks.parquet"):
    df = pd.read_parquet(path)
    pat = pd.read_csv(OUT / "patients.csv").set_index("donor")
    return df, pat


def shuffle_class(df: pd.DataFrame, cls: str, seed) -> pd.DataFrame:
    """Placebo: the rows of class `cls` are reassigned to the donors holding that class, once, within stratum."""
    rng, out = np.random.default_rng(seed), df.copy()
    for _, g in df[df.compartment == cls].groupby("stratum"):
        out.loc[g.index, "donor"] = g.donor.to_numpy()[rng.permutation(len(g))]
    return out


# ---- states and units ----------------------------------------------------------------------------------------------------

def state(tab: Tab, cls: str, stratum: str, pts, wf=WF):
    """Complex state of the rows `pts` of a class, projector fitted on all rows of the class in the stratum, centred over `pts`."""
    allr = tab.rows(cls, stratum)
    L = projector(allr, wf, counts=True)
    x = allr.loc[pts]
    z = S.complex_state(states(x, L, wf, depth=x.depth.to_numpy(float))[0])
    return z - z.mean()


def mgmt_y(tab: Tab, stratum: str):
    """log2(CPM + 1) of MGMT of the malignant pseudobulks with >= MIN_MGMT UMI (Series over donors)."""
    m = tab.rows("malignant", stratum)
    if m is None:
        return None
    m = m[m["MGMT"] >= MIN_MGMT]
    return np.log2(m["MGMT"] / m.depth * 1e6 + 1)


def make_designs(name, tab: Tab, pat: pd.DataFrame, subset=None, placebo=None):
    """Patients (and readout / labels) of a test per stratum; strata with < MIN_N patients (or one call missing) are not eligible.
    Rows: dict(stratum, pts, [y | lab], eligible, n). `subset`: boolean Series over patients (descriptive restrictions)."""
    _, _, kind, c1, c2, _ = TESTS[name]
    ok = set(pat.index if subset is None else pat.index[subset.reindex(pat.index).fillna(False).to_numpy()])
    rng = np.random.default_rng(placebo)
    out = []
    strata = sorted({s for c in (c1, c2) for s in tab.strata(c)})
    for s in strata:
        if kind == "pair":
            a, b = tab.rows(c1, s), tab.rows(c2, s)
            pts = sorted((set(a.index) & set(b.index) & ok) if a is not None and b is not None else [])
            d = dict(stratum=s, pts=pts, n=len(pts), eligible=len(pts) >= MIN_N)
        elif kind == "assoc":
            y = mgmt_y(tab, s)
            st = tab.rows(c2, s)
            pts = sorted((set(y.index) & set(st.index) & ok) if y is not None and st is not None else [])
            yy = y.loc[pts].to_numpy(float) if pts else np.array([])
            if placebo is not None and len(yy):
                yy = yy[rng.permutation(len(yy))]
            d = dict(stratum=s, pts=pts, n=len(pts), y=yy, eligible=len(pts) >= MIN_N)
        else:
            st = tab.rows(c2, s)
            call = pat.mgmt_call
            pts = sorted(p for p in (st.index if st is not None else []) if p in ok and call.get(p) in ("methylated", "unmethylated"))
            lab = np.array([call[p] == "methylated" for p in pts], bool)
            if placebo is not None and len(lab):
                lab = lab[rng.permutation(len(lab))]
            d = dict(stratum=s, pts=pts, n=len(pts), lab=lab,
                     eligible=len(pts) >= MIN_N and 0 < lab.sum() < len(lab))
        out.append(d)
    return out


def make_units(name, designs, tab: Tab, wf=WF):
    """States of the eligible designs on waveform `wf`: list of units (dict with kind, n, stratum, states)."""
    _, _, kind, c1, c2, _ = TESTS[name]
    units = []
    for d in designs:
        if not d["eligible"]:
            continue
        s, pts = d["stratum"], d["pts"]
        if kind == "pair":
            units.append(dict(kind="pair", stratum=s, n=d["n"], za=state(tab, c1, s, pts, wf), zb=state(tab, c2, s, pts, wf)))
        elif kind == "assoc":
            z = state(tab, c2, s, pts, wf)
            u, _ = R.make_unit(dict(stratum=np.array([s] * len(pts)), z=z), y=d["y"])
            if u:
                units.append(dict(u, kind="assoc", stratum=s))
        else:
            z = state(tab, c2, s, pts, wf)
            units.append(dict(kind="label", stratum=s, n=d["n"], z=z, lab=d["lab"], rms=float(np.sqrt(np.mean(abs(z) ** 2)))))
    return units


def unit_stat(u) -> float:
    if u["kind"] == "pair":
        return A.raw_stat(u["za"], u["zb"])
    if u["kind"] == "assoc":
        return R.stat(u)
    return float(abs(u["z"][u["lab"]].mean() - u["z"][~u["lab"]].mean()))


def unit_null(u, n_perm, rng) -> np.ndarray:
    if u["kind"] == "pair":
        return A.perm_null_raw(u, n_perm, rng)
    if u["kind"] == "assoc":
        return R.perm_null(u, n_perm, rng)
    perm = np.argsort(rng.random((n_perm, u["n"])), axis=1)
    L = u["lab"][perm]
    z = u["z"][None, :]
    return np.abs((L * z).sum(1) / L.sum(1) - (~L * z).sum(1) / (~L).sum(1))


def seed_of(name, stratum, extra=0):
    return [SEED, NAMES.index(name), sum(map(ord, stratum)), extra]


def pooled_null(name, units, n_perm):
    w = np.array([u["n"] for u in units], float)
    return w, A.pool_null(w, [unit_null(u, n_perm, np.random.default_rng(seed_of(name, u["stratum"]))) for u in units])


def pooled_obs(units) -> float:
    w = np.array([u["n"] for u in units], float)
    return float(np.sum(w * [unit_stat(u) for u in units]) / w.sum())


def pval(null, obs) -> float:
    return float((1 + np.sum(null >= obs)) / (1 + len(null)))


# ---- foreign RNA ---------------------------------------------------------------------------------------------------------

def marker_share(rows, cls):
    mk = [canonical(g) for g in MARKERS[cls] if canonical(g) in rows.columns]
    return rows[mk].sum(axis=1) / rows.depth


def foreign_a(tab: Tab, x: str, y: str, stratum: str) -> pd.Series:
    """Per patient with both pseudobulks: share of class y's marker counts found in the class-x pseudobulk relative to class y's own
    (clipped to [0, 1]; NaN where class y has no marker count)."""
    rx, ry = tab.rows(x, stratum), tab.rows(y, stratum)
    if rx is None or ry is None:
        return pd.Series(dtype=float)
    pts = sorted(set(rx.index) & set(ry.index))
    sy = marker_share(ry.loc[pts], y)
    return (marker_share(rx.loc[pts], y) / sy.where(sy > 0)).clip(0, 1)


def correct_pair(df: pd.DataFrame, tab: Tab, c1: str, c2: str) -> pd.DataFrame:
    """Table with the expected foreign clock-gene counts removed from both classes of a pair (floored at 0, rounded)."""
    out = df.copy()
    for s in sorted(set(tab.strata(c1)) & set(tab.strata(c2))):
        for x, y in ((c1, c2), (c2, c1)):
            if y not in MARKERS:
                continue
            rx, ry = tab.rows(x, s), tab.rows(y, s)
            a = foreign_a(tab, x, y, s).fillna(0)
            if a.empty:
                continue
            exp = a.to_numpy()[:, None] * rx.loc[a.index].depth.to_numpy()[:, None] * (
                ry.loc[a.index, CLOCK].to_numpy(float) / ry.loc[a.index].depth.to_numpy()[:, None])
            idx = df.index[(df.compartment == x) & (df.stratum == s) & df.donor.isin(a.index)]
            order = df.loc[idx, "donor"].map({d: i for i, d in enumerate(a.index)}).to_numpy()
            out.loc[idx, CLOCK] = np.round(np.clip(df.loc[idx, CLOCK].to_numpy(float) - exp[order], 0, None))
    return out


def foreign_table(name, tab: Tab, pat) -> pd.DataFrame:
    """Marker-based share a per direction and stratum (markers only; no clock gene is read)."""
    rows = []
    for x, y in FOREIGN[name]:
        allv = []
        for s in sorted(set(tab.strata(x)) & set(tab.strata(y))):
            a = foreign_a(tab, x, y, s).dropna()
            allv.append(a)
            rows.append(dict(test=name, contaminated=x, source=y, stratum=s, n=len(a), a_median=a.median(), a_q75=a.quantile(.75), a_max=a.max()))
        a = pd.concat(allv).dropna() if allv else pd.Series(dtype=float)
        rows.append(dict(test=name, contaminated=x, source=y, stratum="ALL", n=len(a), a_median=a.median(), a_q75=a.quantile(.75), a_max=a.max()))
    return pd.DataFrame(rows)


# ---- tab selection of a test ---------------------------------------------------------------------------------------------

def tab_of(name, args):
    """(Tab of the test, placebo seed or None, uncorrected DataFrame): under --placebo the second class is shuffled first."""
    df, kind = G["df"], TESTS[name][2]
    if args.placebo and kind == "pair":
        df = shuffle_class(df, TESTS[name][4], [SEED, 99, NAMES.index(name)])
    return Tab(df), df


def pl_seed(name, args):
    return [SEED, 98, NAMES.index(name)] if args.placebo and TESTS[name][2] != "pair" else None


# ---- power -------------------------------------------------------------------------------------------------------------------

def stage_power(args):
    df, pat = G["df"], G["pat"]
    elig, rows, fr = [], [], []
    base = Tab(df)
    for name in NAMES:
        hyp, role, kind, c1, c2, _ = TESTS[name]
        tab, _ = tab_of(name, args)
        ds = make_designs(name, tab, pat, placebo=pl_seed(name, args))
        for d in ds:
            elig.append(dict(test=name, stratum=d["stratum"], patients=d["n"], eligible=d["eligible"],
                             **({"methylated": int(d["lab"].sum()), "unmethylated": int((~d["lab"]).sum())} if kind == "label" else {})))
        us = make_units(name, ds, tab)
        if name in FOREIGN:
            fr.append(foreign_table(name, base, pat))
        if not us:
            rows.append(dict(test=name, hypothesis=hyp, role=role, patients=0, studies=0, scored=False))
            continue
        w, nul = pooled_null(name, us, args.n_perm)
        sd, n = float(nul.std()), int(w.sum())
        row = dict(test=name, hypothesis=hyp, role=role, patients=n, studies=len(us), study_list=";".join(u["stratum"] for u in us),
                   null_mean=float(nul.mean()), null_sd=sd, mde=MDE_K * sd, n_perm=args.n_perm,
                   scored=bool(n >= MIN_PATIENTS and MDE_K * sd <= MAX_MDE))
        if kind == "label":
            row["mde_over_z_rms"] = MDE_K * sd / float(np.sqrt(np.sum(w * [u["rms"] ** 2 for u in us]) / w.sum()))
        rows.append(row)
        print(f"{name}: {n} patients, {len(us)} studies, null {nul.mean():.4f} +- {sd:.4f}", flush=True)
    pd.DataFrame(elig).to_csv(args.out / "eligibility.csv", index=False)
    P = pd.DataFrame(rows)
    P.to_csv(args.out / "power.csv", index=False)
    pd.concat(fr).to_csv(args.out / "foreign_rna.csv", index=False)
    print(P.drop(columns="study_list", errors="ignore").round(4).to_string(), flush=True)


# ---- calibration ---------------------------------------------------------------------------------------------------------------

ATLASES = ["psychad_hbcc", "psychad_mssm", "seaad_dlpfc", "seaad_mtg"]


def atlas_pairs(ds):
    """Control donors of a cortex atlas with microglia and oligodendrocyte clock-gene counts: (micro, oligo) DataFrames (index donor)."""
    f = PAPER / "common" / "datasets" / ds / "build.py"
    spec = importlib.util.spec_from_file_location(f"build_{ds}", f)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    pb, donors = m.build()
    pb = pb.drop_duplicates(["donor", "compartment"])
    ctrl = set(donors.index[donors.healthy]) if "healthy" in donors else set(donors.index)
    mg, ol = (pb[pb.compartment == f"brain:{c}"].set_index("donor") for c in ("microglia", "oligodendrocyte"))
    d = sorted(set(mg.index) & set(ol.index) & ctrl)
    return mg.loc[d], ol.loc[d]


def thin(x: pd.DataFrame, target: np.ndarray, rng) -> pd.DataFrame:
    """Binomial thinning gene by gene of the clock-gene counts of x (n x genes) to `target` (n x genes); depth scaled by the thinned share."""
    c = x[CLOCK].to_numpy(float)
    p = np.where(c > 0, np.minimum(1.0, target / np.maximum(c, 1)), 0.0)
    new = rng.binomial(c.astype(np.int64), p).astype(float)
    out = x[CLOCK].copy()
    out[:] = new
    out["depth"] = x.depth.to_numpy(float) * new.sum(1) / np.maximum(c.sum(1), 1)
    return out


def calib_one(task):
    """One subsample: pass (p < 0.05 of the pooled P247.0 statistic) and the same with mis-paired donors."""
    ds, k, n_perm = task  # arms: thinned (pass), thinned mis-paired (false positive), native depth
    mg, ol = G["atlas"][ds]
    rng = np.random.default_rng([SEED, 7, k, sum(map(ord, ds))])
    slots = G["slots"]  # DataFrame: slot, stratum, targets
    pick = rng.permutation(len(mg))[:len(slots)]
    pick2 = pick[rng.permutation(len(pick))]  # the same donors, oligodendrocytes of other slots
    res = {}
    for label, po in (("pass", pick), ("mispaired", pick2), ("native", pick)):
        big = label == "native"  # atlas depth, no thinning: separates biology from depth
        m = thin(mg.iloc[pick], np.full_like(G["tm"], 1e12) if big else G["tm"], rng)
        o = thin(ol.iloc[po], np.full_like(G["to"], 1e12) if big else G["to"], rng)
        rows = []
        for cls, x in (("myeloid", m), ("oligodendrocyte", o)):
            t = x.reset_index(drop=True).assign(donor=slots.slot.to_numpy(), compartment=cls, stratum=slots.stratum.to_numpy())
            rows.append(t)
        df = pd.concat(rows, ignore_index=True)
        pat = pd.DataFrame(index=slots.slot)
        tab = Tab(df)
        us = make_units("P247.0", make_designs("P247.0", tab, pat), tab)
        w = np.array([u["n"] for u in us], float)
        nul = A.pool_null(w, [unit_null(u, n_perm, rng) for u in us])
        res[label] = pval(nul, pooled_obs(us))
    return ds, k, res["pass"], res["mispaired"], res["native"]


def stage_calibration(args):
    df, pat = G["df"], G["pat"]
    tab = Tab(df)
    ds0 = make_designs("P247.0", tab, pat)
    rows = [(d["stratum"], p) for d in ds0 if d["eligible"] for p in d["pts"]]
    slots = pd.DataFrame(rows, columns=["stratum", "patient"]).assign(slot=lambda t: [f"slot{i}" for i in range(len(t))])
    G["slots"] = slots
    for key, cls in (("tm", "myeloid"), ("to", "oligodendrocyte")):
        r = df[df.compartment == cls].set_index(["stratum", "donor"])
        G[key] = r.loc[list(zip(slots.stratum, slots.patient)), CLOCK].to_numpy(float)
    print(f"P247.0 design: {len(slots)} patients in {slots.stratum.nunique()} strata", flush=True)
    G["atlas"] = {}
    tasks, n_sub, out = [], (10 if args.smoke else 200), []
    for ds in ATLASES:
        try:
            mg, ol = atlas_pairs(ds)
        except Exception as e:  # noqa: BLE001
            out.append(dict(atlas=ds, note=f"not usable: {type(e).__name__}: {e}"))
            continue
        if len(mg) < len(slots):
            out.append(dict(atlas=ds, control_donors=len(mg), note=f"fewer control donors with both classes than P247.0 patients ({len(slots)})"))
            continue
        G["atlas"][ds] = (mg, ol)
        out.append(dict(atlas=ds, control_donors=len(mg), note=""))
        tasks += [(ds, k, args.n_perm) for k in range(n_sub)]
    t0, res = time.time(), []
    with Pool(args.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for i, r in enumerate(pool.imap_unordered(calib_one, tasks)):
            res.append(r)
            if (i + 1) % 100 == 0 or i + 1 == len(tasks):
                print(f"{i + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    from scipy.stats import binomtest
    R_ = pd.DataFrame(res, columns=["atlas", "k", "p_pass", "p_mispaired", "p_native"])
    for o in out:
        r = R_[R_.atlas == o["atlas"]]
        if not len(r):
            continue
        for key, col in (("pass_rate", "p_pass"), ("false_positive_rate", "p_mispaired"), ("pass_rate_native_depth", "p_native")):
            b = binomtest(int((r[col] < .05).sum()), len(r))
            ci = b.proportion_ci(.95, method="exact")
            o.update({key: b.proportion_estimate, f"{key}_lo": ci.low, f"{key}_hi": ci.high})
        o.update(subsamples=len(r), patients=len(slots), strata=slots.stratum.nunique(), n_perm=args.n_perm,
                 specificity="not computed: atlases on disk carry a 45-gene panel, no genome for matched pseudo sets")
    C = pd.DataFrame(out)
    C.to_csv(args.out / "calibration.csv", index=False)
    print(C.round(3).to_string(), flush=True)


# ---- tests ---------------------------------------------------------------------------------------------------------------------

def pseudo_task(t):
    name, kind, i = t
    tab, designs, w = G["tabs"][name], G["designs"][name], G["ws"][name][kind][i]
    try:
        us = make_units(name, designs, tab, w)
        return t, pooled_obs(us)
    except Exception:  # noqa: BLE001
        return t, np.nan


def pseudo_sets(name, tab: Tab, designs, n_pseudo):
    """Matched and responsive sets of one test (reference: the second class's pseudobulks of the test patients)."""
    c2 = TESTS[name][4]
    ref = pd.concat([tab.rows(c2, d["stratum"]).loc[d["pts"]] for d in designs if d["eligible"]])
    genes = [c for c in ref.columns if c not in ("compartment", "stratum", "study", "suspension", "n_cells", "depth")]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(CLOCK)]
    return {k: pseudo_waveforms(lm, n_pseudo, np.random.default_rng([SEED, NAMES.index(name), 44, j]), WF, k)
            for j, k in enumerate(("matched", "responsive"), 1)}


def run_test(name, designs, tab: Tab, n_perm):
    """Pooled observed statistic, null and p of one test on the given designs (None if no eligible study)."""
    us = make_units(name, designs, tab)
    if not us:
        return None
    w, nul = pooled_null(name, us, n_perm)
    obs = pooled_obs(us)
    out = dict(patients=int(w.sum()), studies=len(us), statistic=obs, null_mean=float(nul.mean()), null_sd=float(nul.std()), p=pval(nul, obs))
    if us[0]["kind"] == "assoc":
        c = np.sum(w * [R.cross(u["y"], u["z"]) / u["den"] for u in us]) / w.sum()
        out["implied_peak_h"] = R.peak_h(c)
    return out, us


def stage_tests(args):
    df, pat = G["df"], G["pat"]
    pw = pd.read_csv(OUT / "power.csv").set_index("test") if (OUT / "power.csv").exists() else None
    if pw is None and not args.placebo:
        raise SystemExit("run the power stage first")
    scored = {n: bool(pw.loc[n, "scored"]) if pw is not None and n in pw.index else args.placebo for n in NAMES}
    G.update(tabs={}, designs={}, ws={})
    tests, per, pseudo_rows = [], [], []
    for name in NAMES:
        tab, _ = tab_of(name, args)
        G["tabs"][name] = tab
        G["designs"][name] = make_designs(name, tab, pat, placebo=pl_seed(name, args))
        if not scored[name] and not (args.force or args.placebo):
            print(f"{name}: not scored (power note); --force to run anyway", flush=True)
            continue
        G["ws"][name] = pseudo_sets(name, tab, G["designs"][name], args.n_pseudo)
    run = [n for n in NAMES if n in G["ws"]]
    tasks = [(n, k, i) for n in run for k in ("matched", "responsive") for i in range(args.n_pseudo)]
    ps, t0 = {}, time.time()
    with Pool(args.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for j, (t, v) in enumerate(pool.imap_unordered(pseudo_task, tasks, chunksize=4)):
            ps[t] = v
            if (j + 1) % 200 == 0 or j + 1 == len(tasks):
                print(f"pseudo {j + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    for name in run:
        hyp, role, kind, c1, c2, fam = TESTS[name]
        tab, designs = G["tabs"][name], G["designs"][name]
        res = run_test(name, designs, tab, args.n_perm)
        if res is None:
            continue
        r, us = res
        row = dict(test=name, hypothesis=hyp, role=role, scored=scored[name], **r)
        for k in ("matched", "responsive"):
            v = np.array([ps[(name, k, i)] for i in range(args.n_pseudo)])
            v = v[np.isfinite(v)]
            row.update({f"{k}_median": float(np.median(v)), f"{k}_q95": float(np.quantile(v, .95)),
                        f"p_{k}": float((1 + np.sum(v >= r["statistic"])) / (1 + len(v))), f"n_sets_{k}": len(v)})
            pseudo_rows += [dict(test=name, kind=k, set=i, value=x) for i, x in enumerate(v)]
        row["clock_specific"] = bool(r["statistic"] > row["matched_q95"] and r["statistic"] > row["responsive_q95"])
        if kind == "pair" and name != "TNK_myeloid":
            tc = Tab(correct_pair(tab.df, tab, c1, c2))
            rc = run_test(name, make_designs(name, tc, pat), tc, args.n_perm)
            if rc:
                row.update(statistic_corrected=rc[0]["statistic"], p_corrected=rc[0]["p"])
        tests.append(row)
        for u in us:
            per.append(dict(test=name, study=u["stratum"], patients=u["n"], statistic=unit_stat(u), at_least_12=u["n"] >= MIN_STUDY))
    T = pd.DataFrame(tests)
    T["p_holm"] = np.nan
    for fam in ("pairing", "h4"):
        m = [t for t in T.test if TESTS[t][5] == fam]
        if m:
            T.loc[T.test.isin(m), "p_holm"] = A.holm(T[T.test.isin(m)].p)
    ok = T.clock_specific & (T.p_holm.fillna(T.p) < .05)
    if "p_corrected" in T:  # H247.1-2 must also hold after the foreign-RNA correction
        ok &= ~T.test.isin(["H247.1", "H247.2"]) | (T.p_corrected < .05)
    T["pass"] = np.where(T.role == "descriptive", np.nan, ok).astype(float)
    T.to_csv(args.out / "tests.csv", index=False)
    pd.DataFrame(per).to_csv(args.out / "per_study.csv", index=False)
    pd.DataFrame(pseudo_rows).to_csv(args.out / "pseudo.csv", index=False)
    T.to_csv(args.out / "pooled.csv", index=False)  # extended below by the descriptive restrictions
    descriptives(args, T)
    print(T.round(4).to_string(), flush=True)


def descriptives(args, T):
    """Registered descriptives: restrictions (cells / nuclei, primary / recurrent, H247.3 within methylation status), MGMT percentile among
    200 expression-matched malignant genes, foreign clock share, coupling c of the malignant class on an anchored host."""
    df, pat = G["df"], G["pat"]
    pooled = [T.assign(subset="all")]
    subsets = {"cells": pat.suspension == "cell", "nuclei": pat.suspension == "nucleus", "primary": pat.stage_group == "primary",
               "recurrent": pat.stage_group == "recurrent"}
    for name in ("P247.0", "H247.1", "H247.2", "H247.3"):
        for sub, m in subsets.items():
            tab = G["tabs"][name]
            r = run_test(name, make_designs(name, tab, pat, subset=m, placebo=pl_seed(name, args)), tab, args.n_perm)
            if r:
                pooled.append(pd.DataFrame([dict(test=name, subset=sub, **r[0])]))
    for sub in ("methylated", "unmethylated"):
        tab = G["tabs"]["H247.3"]
        r = run_test("H247.3", make_designs("H247.3", tab, pat, subset=pat.mgmt_call == sub, placebo=pl_seed("H247.3", args)), tab, args.n_perm)
        if r:
            pooled.append(pd.DataFrame([dict(test="H247.3", subset=f"within_{sub}", **r[0])]))
    pd.concat(pooled, ignore_index=True).to_csv(args.out / "pooled.csv", index=False)
    desc = []
    # MGMT percentile among expression-matched malignant genes (y), against the myeloid state of H247.3
    if "H247.3" in G["ws"] and (T.test == "H247.3").any():
        tab, designs = G["tabs"]["H247.3"], G["designs"]["H247.3"]
        elig = [d for d in designs if d["eligible"]]
        genes = [c for c in df.columns if c not in ("donor", "compartment", "stratum", "study", "suspension", "n_cells", "depth")]
        excl = set(CLOCK) | {canonical("MGMT")}
        lc = {d["stratum"]: np.log2(tab.rows("malignant", d["stratum"]).loc[d["pts"], genes].div(
            tab.rows("malignant", d["stratum"]).loc[d["pts"]].depth, axis=0).to_numpy(float) * 1e6 + 1) for d in elig}
        w = np.array([d["n"] for d in elig], float)
        mean = pd.Series(np.sum([wi * lc[d["stratum"]].mean(0) for wi, d in zip(w, elig)], 0) / w.sum(), index=genes)
        cand = (mean.drop([g for g in excl if g in mean.index]) - mean["MGMT"]).abs().sort_values().index[:400]
        pick = np.random.default_rng([SEED, 7]).choice(cand, size=200, replace=False)
        rs = []
        for g in pick:
            gi, st = genes.index(g), []
            for d in elig:
                z = state(tab, "myeloid", d["stratum"], d["pts"])
                u, _ = R.make_unit(dict(stratum=np.array([d["stratum"]] * d["n"]), z=z), y=lc[d["stratum"]][:, gi])
                if u:
                    st.append((u["n"], R.stat(u)))
            rs.append(np.average([s for _, s in st], weights=[n for n, _ in st]) if st else np.nan)
        rs = np.array(rs)
        obs = float(T.loc[T.test == "H247.3", "statistic"].iloc[0])
        pd.DataFrame(dict(gene=pick, mean_log2cpm=mean[pick].to_numpy(), statistic=rs, mgmt_statistic=obs,
                          mgmt_mean_log2cpm=mean["MGMT"])).to_csv(args.out / "matched_genes.csv", index=False)
        desc.append(dict(item="mgmt_percentile_among_matched", value=float(np.mean(rs[np.isfinite(rs)] <= obs) * 100)))
    pd.DataFrame(desc).to_csv(args.out / "descriptive.csv", index=False)
    clock_share(args)
    coupling(args, T)


def clock_share(args):
    """Expected foreign share of the clock-gene counts per pair, stratum and direction (the correction of the `:corrected` reruns)."""
    rows = []
    for name, dirs in FOREIGN.items():
        if name == "TNK_myeloid":
            continue
        if name not in G["tabs"]:
            continue
        tab = G["tabs"][name]
        for x, y in dirs:
            for s in sorted(set(tab.strata(x)) & set(tab.strata(y))):
                a = foreign_a(tab, x, y, s).fillna(0)
                rx, ry = tab.rows(x, s).loc[a.index], tab.rows(y, s).loc[a.index]
                exp = a.to_numpy()[:, None] * rx.depth.to_numpy()[:, None] * ry[CLOCK].to_numpy(float) / ry.depth.to_numpy()[:, None]
                share = exp.sum(0) / rx[CLOCK].sum(0).to_numpy()
                rows.append(dict(test=name, contaminated=x, source=y, stratum=s, n=len(a), clock_share_median=float(np.median(share)),
                                 clock_share_max=float(share.max()), confounded=bool(np.median(share) > .5)))
    pd.DataFrame(rows).to_csv(args.out / "foreign_clock_share.csv", index=False)


def coupling(args, T):
    """Descriptive: `clocksync` coupling of the malignant class on a host anchored on myeloid and oligodendrocyte cells; only if the
    gate over the three classes and both its pseudo-set companions pass (pairing_z above the q95 of 200 matched and 200 responsive sets)."""
    df, pat = G["df"], G["pat"]
    three = df[df.compartment.isin(["malignant", "myeloid", "oligodendrocyte"])]
    keep = three.groupby("donor").compartment.transform("nunique") == 3
    t3 = three[keep].reset_index(drop=True)
    if t3.donor.nunique() < 12:
        pd.DataFrame([dict(item="coupling", note=f"{t3.donor.nunique()} patients with the three classes")]).to_csv(args.out / "coupling.csv", index=False)
        return
    n_perm = min(args.n_perm, 1000)
    base = t3[["donor", "compartment", "stratum", "n_cells", "depth", *CLOCK]]
    ds = cs.Dataset(base, anchors=["myeloid", "oligodendrocyte"], waveform=WF, min_donors=4, name="a247")
    g = ds.gate(n_perm=n_perm, seed=SEED)
    genes = [c for c in t3.columns if c not in ("donor", "compartment", "stratum", "study", "suspension", "n_cells", "depth")]
    frac = t3[genes].div(t3.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(t3[genes].mean() >= 1) | frac.columns.isin(CLOCK)]
    row = dict(item="coupling", patients=t3.donor.nunique(), gate_p_pairing=g.p_pairing, gate_p_acrophase=g.p_acrophase, gate_pairing_z=g.pairing_z)
    for j, k in enumerate(("matched", "responsive"), 1):
        zs = []
        for w in pseudo_waveforms(lm, min(args.n_pseudo, 200), np.random.default_rng([SEED, 55, j]), WF, k):
            try:
                zs.append(cs.Dataset(t3[["donor", "compartment", "stratum", "n_cells", "depth", *w.genes]], anchors=None, waveform=w,
                                     min_donors=4).gate(n_perm=min(n_perm, 200), seed=SEED).pairing_z)
            except Exception:  # noqa: BLE001
                pass
        row[f"{k}_q95"] = float(np.quantile(zs, .95)) if zs else np.nan
    row["companions_pass"] = bool(g.passes and g.pairing_z > row["matched_q95"] and g.pairing_z > row["responsive_q95"])
    if row["companions_pass"]:
        fit = ds.fit(n_boot=200, seed=SEED, targets=["malignant"])
        row.update({f"fit_{k}": v for k, v in fit.table.loc["malignant"].items()})
    pd.DataFrame([row]).to_csv(args.out / "coupling.csv", index=False)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "calibration", "tests"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--placebo", action="store_true", help="code test: shuffle the second class / y / labels within stratum, write to output/placebo")
    ap.add_argument("--force", action="store_true", help="tests: also run tests flagged not scored")
    ap.add_argument("--smoke", action="store_true", help="200 permutations, 5 pseudo sets, 10 calibration subsamples")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    assert args.jobs <= 3, "at most 3 workers for A247"
    cli.setup(args)
    if args.placebo and args.stage != "tests" and args.stage != "power":
        raise SystemExit("--placebo applies to the power and tests stages")
    args.out = args.out or (OUT / "placebo" if args.placebo else OUT)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.smoke:
        args.n_perm, args.n_pseudo = 200, 5
        if args.stage == "calibration":
            args.n_perm = 50
    G["df"], G["pat"] = load_table()
    {"power": stage_power, "calibration": stage_calibration, "tests": stage_tests}[args.stage](args)


if __name__ == "__main__":
    main()
