"""A230: is the epithelial clock offset from the patient's stromal host different in inflamed vs non-inflamed intestine?

Registered in paper/docs/registered_experiments.md (A230). Cohorts: gse266616 (blind primary, Crohn's disease, inflamed vs
adjacent-to-inflamed / non-inflamed at the same site) and taurus (non-blind secondary, inflamed vs non-inflamed).
Tables: `paper/common/datasets/ibd_gse266616/build.py` (marker classes, >= 50 cells, median UMI >= 1,500, donor = sample).

  --stage power   Eligibility, cells / UMI per class and arm, host-only gate per arm (a precondition), offset MDE =
                  2.49 x SD of offset_diff_h under within-patient label swaps, coupling MDE (`clocksync.power.mde`).
                  NO OBSERVED CONTRAST IS EVALUATED OR STORED. The swap null is obtained by running `Dataset.contrast` on a
                  RANDOM within-patient swap of the arm labels (a draw of the null) and keeping only `Contrast.null`: permuting
                  a swapped labelling within patients has the same distribution as permuting the true one, so the null SD is
                  that of the registered test, and the contrast of the true labelling is never computed. The "observed" row
                  of that call belongs to the random swap and is discarded before anything is stored or printed.
  --stage tests   The registered contrast (H230.1-H230.3) with the true arms. Run only after the dated power note.

    python compute_a230.py --stage power --cohort gse266616|taurus [--jobs 3] [--n-perm 10000] [--n-rep 100]
    python compute_a230.py --stage tests --cohort gse266616|taurus [--out DIR]

Choices not fixed by the registration (all are data-independent)
- Pseudobulk tables (one per cohort) are read from `output/data/<cohort>/{pseudobulk,donors}.csv` (`build.py <cohort> <dir>`).
- Anchors: gut:fibroblast, gut:endothelial, + gut:mural if present in >= 60 % of all samples with an arm. A sample is
  eligible if it has epithelium and >= 2 of the anchors; a pair (patient-site; TAURUS: patient, each
  inflamed sample against each non-inflamed one is not formed: the patient's inflamed vs non-inflamed samples are the arms)
  is eligible if every sample of the patient-site pair qualifies; a patient is eligible with >= 1 eligible pair. TAURUS: the
  patient's samples are eligible one by one and the patient needs >= 1 eligible sample in each arm.
- Permutation count (tests): all within-patient label arrangements if <= 10,000 (drawn at random, `Dataset.contrast` has no
  enumeration, as many draws as arrangements), else 10,000; seed 0. Power stage: always 10,000 random swaps (SD precision). Gate: 1,000 permutations per arm.
- `targets` = [epithelium, fibroblast] when the companion is estimable (registered: epithelium + companion): `Contrast.module`
  needs the companion's draws; the epithelial rows are unchanged by it. The companion, an anchor, is read on its
  leave-one-out host, which needs >= 3 anchors: with fibroblast + endothelial only (no mural) it has no estimate, the
  companion clause of H230.1 / H230.2 is then not estimable (`companion_estimable` False; epithelium-only MDE reported).
- Offset difference of epithelium minus companion (module statistic) is wrapped to [-12, 12) h. The offset MDE flag uses the
  epithelial offset_diff_h SD as registered; the module SD is reported next to it.
- H230.1 passes if p_offset (epithelium, two-sided) < 0.05 and the module p < 0.05. H230.2: d* change of epithelium
  (`p_star`) and its companion difference (`Contrast.companion`, same permutations); Holm over H230.2 and H230.3 in
  `--stage tests` when both cohorts' score files exist (H230.2 gse266616, H230.3 = H230.1 rule in taurus, p = larger of the two p values).
"""

from __future__ import annotations

import argparse
import math
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import power, stats as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parent / "output"
EPI, COMP, FIB, ENDO, MURAL = "gut:epithelial", "gut:fibroblast", "gut:fibroblast", "gut:endothelial", "gut:mural"
LEVELS = ("reference", "inflamed")
MDE_K = 2.49
MIN_PATIENTS, MAX_OFFSET_MDE, MURAL_SHARE, N_PERM_MAX, N_GATE, SEED = 6, 6.0, 0.6, 10_000, 1000, 0
say = lambda *a: print(*a, flush=True)  # noqa: E731


# ---- data and eligibility -----------------------------------------------------------------------------------------
def load(cohort: str, data: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = (data or OUT / "data") / cohort
    pb = pd.read_csv(d / "pseudobulk.csv", keep_default_na=False)
    don = pd.read_csv(d / "donors.csv", index_col=0, keep_default_na=False)
    return pb, don


def eligible(pb: pd.DataFrame, don: pd.DataFrame, cohort: str) -> dict:
    """Eligible samples / patients and the anchor set. `don`: donor table (index sample; `arm`, `patient`, `kept`)."""
    don = don[(don.arm != "") & don.kept.astype(bool)].copy()
    have = pb[pb.donor.isin(don.index)].groupby("donor").compartment.agg(set).reindex(don.index).map(lambda s: s if isinstance(s, set) else set())
    share = float(have.map(lambda s: MURAL in s).mean())  # mural: share of all samples with an arm
    anchors = [FIB, ENDO] + ([MURAL] if share >= MURAL_SHARE else [])
    don["sample_ok"] = have.map(lambda s: EPI in s and len(set(anchors) & s) >= 2)
    unit = "patient_site" if cohort == "gse266616" else "patient"
    g = don.groupby(unit)
    if cohort == "gse266616":  # a pair needs both of its samples
        pair_ok = g.sample_ok.all() & (g.size() == 2)
        keep = don.index[don[unit].map(pair_ok)]
    else:  # every eligible sample; the patient needs both arms
        keep = don.index[don.sample_ok]
    don = don.loc[keep]
    pat_ok = don.groupby("patient").arm.nunique() == 2
    don = don[don.patient.map(pat_ok)]
    return dict(donors=don, anchors=anchors, mural_share=share, n_patients=don.patient.nunique(),
                pb=pb[pb.donor.isin(don.index)].reset_index(drop=True))


def n_arrangements(lab: pd.Series, blocks: pd.Series) -> int:
    """Number of distinct within-block label arrangements."""
    tot = 1
    for _, g in lab.groupby(blocks.reindex(lab.index)):
        c = g.value_counts()
        tot *= math.factorial(len(g)) // math.prod(math.factorial(v) for v in c)
    return tot


def dataset(pb, don, anchors, name):
    return cs.Dataset(pb, anchors=anchors, donors=don, waveform="human", name=name)


def companion_ok(ds) -> bool:
    """The companion (fibroblast, an anchor) is read on its leave-one-out host, which needs >= 3 anchors present
    (`clocksync._core.host._hosts`); with the two fixed anchors (no mural) it has no estimate."""
    return COMP in ds.anchors and len(ds.anchors) >= 3


def targets_of(ds) -> list:
    return [EPI, COMP] if companion_ok(ds) else [EPI]


def n_perm_for(lab, blocks) -> tuple[int, int]:
    n = n_arrangements(lab, blocks)
    return (n if n <= N_PERM_MAX else N_PERM_MAX), n


# ---- power -----------------------------------------------------------------------------------------------------------
def gates(ds, lab) -> dict:
    """Host-only gate (anchors among themselves) in each arm: the precondition `require_host_gate_each_arm` applies."""
    out = {}
    for g in LEVELS:
        try:
            gt = ds.subset(donors=lab.index[lab == g], compartments=ds.anchors, anchors=ds.anchors).gate(n_perm=N_GATE, seed=SEED)
            out[g] = gt.to_series().to_dict() | dict(passes=bool(gt.passes))
        except ValueError as e:
            out[g] = dict(passes=False, error=str(e))
    return out


def swap_null(ds, lab, blocks, n_perm: int) -> pd.DataFrame:
    """Null draws of the registered statistics under within-patient swaps; the identity labelling is never evaluated.
    Columns: offset_epi, offset_comp, offset_module (epithelium - companion, wrapped), d_star_epi, d_star_comp, d_epi."""
    rng = np.random.default_rng(SEED)
    pseudo = S.permute_within(lab, blocks, rng)  # one random swap: stands in for the observed labelling
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = ds.contrast(pseudo, levels=LEVELS, targets=targets_of(ds), companion=COMP if companion_ok(ds) else None, blocks=blocks, n_perm=n_perm, n_boot=0,
                        seed=SEED + 1, require_host_gate_each_arm=False)
    null = c.null  # draws only; c.table / c.observed belong to the random swap and are dropped here
    del c
    off = null["offset_diff_h"]
    nan = pd.Series(np.nan, index=off.index)
    ok = COMP in off
    return pd.DataFrame({"offset_epi": off[EPI], "offset_comp": off[COMP] if ok else nan,
                         "offset_module": S.wrap(off[EPI] - off[COMP]) if ok else nan, "d_star_epi": null["d_star"][EPI],
                         "d_star_comp": null["d_star"][COMP] if ok else nan, "d_epi": null["d"][EPI]})


def stage_power(cohort: str, a) -> None:
    OUT.mkdir(exist_ok=True)
    pb, don = load(cohort)
    e = eligible(pb, don, cohort)
    d, anchors, pbe = e["donors"], e["anchors"], e["pb"]
    say(f"{cohort}: {e['n_patients']} eligible patients, {len(d)} samples; anchors {anchors} (mural share {e['mural_share']:.2f})")
    d.to_csv(OUT / f"eligible_samples_{cohort}.csv")
    # cells and UMI per class and arm (descriptive, over eligible samples)
    pbe = pbe.assign(arm=pbe.donor.map(d.arm), umi_per_cell=pbe.depth / pbe.n_cells)
    cells = pbe.groupby(["compartment", "arm"]).agg(samples=("donor", "nunique"), median_cells=("n_cells", "median"),
                                                    min_cells=("n_cells", "min"), median_mean_umi_per_cell=("umi_per_cell", "median"))
    cells.to_csv(OUT / f"cells_{cohort}.csv")
    med_umi = d.groupby("arm").median_umi.median()
    row = dict(cohort=cohort, patients=e["n_patients"], samples=len(d), pairs=int(d.patient_site.nunique()) if cohort == "gse266616" else np.nan,
               anchors="+".join(anchors), mural_share=e["mural_share"], median_umi_reference=med_umi.get("reference"),
               median_umi_inflamed=med_umi.get("inflamed"), n_reference=int((d.arm == "reference").sum()),
               n_inflamed=int((d.arm == "inflamed").sum()))
    if e["n_patients"] < 4:
        pd.DataFrame([row | dict(note="too few eligible patients")]).to_csv(OUT / f"power_{cohort}.csv", index=False)
        return
    ds = dataset(pbe, d, anchors, cohort)
    lab, blocks = d.arm.astype(str), d.patient.astype(str)
    lab.index = blocks.index = lab.index.astype(str)
    # host-only gate per arm (precondition)
    gt = gates(ds, lab)
    pd.DataFrame(gt).T.to_csv(OUT / f"gate_{cohort}.csv")
    for g in LEVELS:
        say(f"gate {g}: {gt[g]}")
    # offset MDE from the swap null
    n_perm, n_all = n_perm_for(lab, blocks)
    n_perm = a.n_perm or N_PERM_MAX  # the null SD is estimated from 10,000 random swaps even when fewer arrangements exist
    say(f"within-patient arrangements: {n_all}; null draws: {n_perm}")
    null = swap_null(ds, lab, blocks, n_perm)
    null.to_csv(OUT / f"swap_null_{cohort}.csv.gz", index=False)
    sd = null.std(ddof=1)
    row |= dict(companion_estimable=companion_ok(ds), arrangements=n_all, n_null=int(null.offset_epi.notna().sum()), sd_offset_epi=sd.offset_epi,
                sd_offset_comp=sd.offset_comp, sd_offset_module=sd.offset_module, offset_mde_h=MDE_K * sd.offset_epi,
                offset_mde_module_h=MDE_K * sd.offset_module, sd_d_star_epi=sd.d_star_epi, d_star_mde=MDE_K * sd.d_star_epi,
                sd_d_star_comp=sd.d_star_comp)
    # coupling MDE (clocksync.power.mde; permutes labels and damps the epithelial host-locked part: no observed contrast)
    if a.n_rep:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = power.mde(ds.study, lab, anchors, EPI, n_rep=a.n_rep, n_perm=200, seed=SEED, case="inflamed",
                          control="reference", within_stratum=False, blocks=blocks, jobs=a.jobs)
        m["curve"].to_csv(OUT / f"coupling_mde_curve_{cohort}.csv", index=False)
        row |= dict(coupling_mde_d=m["mde"], c_pooled_epi=m["c_pooled"])
    row |= dict(gate_pass_reference=gt["reference"]["passes"], gate_pass_inflamed=gt["inflamed"]["passes"],
                gate_p_pairing_reference=gt["reference"].get("p_pairing"), gate_p_pairing_inflamed=gt["inflamed"].get("p_pairing"),
                gate_p_acrophase_reference=gt["reference"].get("p_acrophase"), gate_p_acrophase_inflamed=gt["inflamed"].get("p_acrophase"))
    row |= dict(flag_patients=row["patients"] >= MIN_PATIENTS, flag_offset_mde=row["offset_mde_h"] <= MAX_OFFSET_MDE,
                flag_gate=bool(row["gate_pass_reference"] and row["gate_pass_inflamed"]))
    row["scored"] = bool(row["flag_patients"] and row["flag_offset_mde"] and row["flag_gate"])
    pd.DataFrame([row]).to_csv(OUT / f"power_{cohort}.csv", index=False)
    say(pd.Series(row).to_string())


# ---- tests -----------------------------------------------------------------------------------------------------------
def holm(p: dict) -> dict:
    items = sorted(p.items(), key=lambda kv: kv[1])
    out, run = {}, 0.0
    for i, (k, v) in enumerate(items):
        run = max(run, min(1.0, (len(items) - i) * v))
        out[k] = run
    return out


def stage_tests(cohort: str, a) -> None:
    out = Path(a.out) if a.out else OUT
    out.mkdir(parents=True, exist_ok=True)
    pb, don = load(cohort)
    e = eligible(pb, don, cohort)
    d, anchors, pbe = e["donors"], e["anchors"], e["pb"]
    ds = dataset(pbe, d, anchors, cohort)
    lab, blocks = d.arm.astype(str), d.patient.astype(str)
    lab.index = blocks.index = lab.index.astype(str)
    n_perm, n_all = n_perm_for(lab, blocks)
    n_perm = a.n_perm or n_perm
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = ds.contrast(groups=lab, levels=LEVELS, targets=targets_of(ds), companion=COMP if companion_ok(ds) else None, blocks=blocks, n_perm=n_perm,
                        n_boot=a.n_boot, seed=SEED, require_host_gate_each_arm=True, n_gate=N_GATE)
    c.table.to_csv(out / f"contrast_{cohort}.csv")
    pd.DataFrame(c.gate).T.to_csv(out / f"tests_gate_{cohort}.csv")
    if c.companion is not None:
        c.companion.to_csv(out / f"companion_{cohort}.csv", index=False)
    if c.null is None:
        pd.DataFrame([dict(cohort=cohort, status="not estimable (host-only gate failed)")]).to_csv(out / f"score_{cohort}.csv", index=False)
        holm_files(out)
        return
    t = c.table.loc[EPI]
    if companion_ok(ds):
        mod = c.module({"epi_minus_fib": lambda x: float(S.wrap(x[EPI] - x[COMP]))}, scales=("offset_diff_h",))
        mod.to_csv(out / f"module_{cohort}.csv", index=False)
        comp = c.companion.iloc[0]
    else:  # companion not estimable: the specificity clause cannot be evaluated
        mod = pd.DataFrame(dict(value=[np.nan], p=[np.nan]))
        comp = pd.Series(dict(d_star_minus_companion=np.nan, p_star=np.nan))
    row = dict(cohort=cohort, patients=e["n_patients"], samples=len(d), n_perm=n_perm, arrangements=n_all,
               offset_diff_h=t.offset_diff_h, p_offset=t.p_offset, module_value=mod.value[0], module_p=mod.p[0],
               d=t.d, d_star=t.d_star, p_star=t.p_star, d_star_minus_companion=comp.d_star_minus_companion,
               p_star_companion=comp.p_star, companion_estimable=companion_ok(ds), status=t.status)
    row["H230.1_rule"] = bool(row["p_offset"] < 0.05 and row["module_p"] < 0.05)  # False when the companion is not estimable
    pd.DataFrame([row]).to_csv(out / f"score_{cohort}.csv", index=False)
    holm_files(out)


def holm_files(out: Path) -> None:
    """Holm over H230.2 (gse266616 d*) and H230.3 (taurus) once both score files exist; a test that is not estimable drops out."""
    other = {k: out / f"score_{k}.csv" for k in ("gse266616", "taurus")}
    if all(p.exists() for p in other.values()):
        g, tr = (pd.read_csv(other[k]).iloc[0] for k in ("gse266616", "taurus"))
        ps = {"H230.2": g.get("p_star", np.nan), "H230.3": max(tr.get("p_offset", np.nan), tr.get("module_p", np.nan))}
        ps = {k: float(v) for k, v in ps.items() if np.isfinite(v)}  # a test that is not estimable drops out of the family
        h = holm(ps)
        pd.DataFrame([dict(test=k, holm_p=v) for k, v in h.items()]).to_csv(out / "holm.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "tests"], required=True)
    ap.add_argument("--cohort", choices=["gse266616", "taurus"], required=True)
    ap.add_argument("--jobs", type=int, default=1, help="workers of the coupling MDE (at most 3 today)")
    ap.add_argument("--n-perm", type=int, default=0, help="override the permutation count (smoke runs)")
    ap.add_argument("--n-rep", type=int, default=100, help="replicates per damping level of the coupling MDE (0: skip)")
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--out", default=None, help="tests: output directory (default output/)")
    a = ap.parse_args()
    limit_threads(1 if a.jobs > 1 else 4)
    (stage_power if a.stage == "power" else stage_tests)(a.cohort, a)
