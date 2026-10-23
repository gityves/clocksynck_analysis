"""A27 part 2, the CHIRAL DIP arms: D5b, D6 (DIP arm), D6b (DIP parts), D6c part 1 (post hoc, exploratory).

Post hoc, exploratory (docs/PLAN_internal_desync.md (archive), A27 part 2; no hypothesis ids). CHIRAL's donor internal phase (DIP, hours
in clock time; `data/chiral_DIP.csv` of the repository, controlled per-donor data, in memory only) replaces the clock-gene host as
the internal time of ventilated donors. Estimands as the host arm (`compute_internal_external.d5_tissue`, `compute_glucocorticoid`):
  D5b  per tissue and gene in ventilated donors, partial F of the DIP (cos, sin) given the time of death and the reverse, counts
       of significant genes against permutation nulls, internal-only genes, phase against fast deaths (frame rotation); agreement
       of the host phase with the DIP and with the time of death (fast deaths, ventilated); overlap of the internal-only gene sets
       of the host and DIP arms (needs diag5_genes.csv.gz of compute_internal_external.py)
  D6   a priori gene classes (`compute_glucocorticoid.classes`) among fast-death rhythmic genes against matched random sets (same
       tissue and p_fast decile, 2,000 draws), on the DIP gene table; the amplitude ratio is comparable here (the DIP is in hours)
  D6b  GC score: R^2 on host f, DIP and time of death; host-DIP agreement, also after removing the GC direction from f and from
       the DIP, and by GC tertile (host-DIP, host-time of death)
  D6c  part 1: GC score profile along the DIP (4 h bins), the DIP phase by GC tertile; fast deaths and ventilated
Aggregates only.

    python compute_dip_arm.py [--subsample N] [--jobs 4]
Outputs (paper/results/figure_3/output/gtex_cues/): diag5_dip_{agreement,summary,overlap}.csv, diag5_dip_genes.csv.gz,
diag6_classes_dip.csv, diag6b_gc_dip.csv, diag6c_gc_dip_profile.csv
Legacy: legacy/desync/validation/glucocorticoid/scripts/{diag5_dip,diag6_classes (DIP arm),diag6b_gc (DIP),diag6c_gc_zeitgeber (part 1)}.py
Local: the DIP table, the DIP-as-host regression (`compute_internal_external.d5_tissue`).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _genes as G  # noqa: E402
import _shared as S  # noqa: E402
import compute_glucocorticoid as GC  # noqa: E402
import compute_internal_external as IE  # noqa: E402
from _shared import W  # noqa: E402

FOLDER = G.FOLDER
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import PAPER_DATA  # noqa: E402
DIP_FILE = PAPER_DATA / "chiral_DIP.csv"


def read_dip() -> pd.Series:
    """CHIRAL DIP (h) by donor (`GTEX-<id>`)."""
    d = pd.read_csv(DIP_FILE, sep=";", decimal=",")
    return pd.Series(d.hour.to_numpy(float), index=("GTEX-" + d["SUBJ.ID"].astype(str)).to_numpy())


def r2(D, y):
    return float(1 - ((y - D @ np.linalg.lstsq(D, y, rcond=None)[0]) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def agreement(pb, don, dip) -> pd.DataFrame:
    rows = []
    for lab, hd in (("fast deaths", (1, 2)), ("ventilator", (0,))):
        host = S.host_states(pb[pb.donor.isin(don.index[don.hardy.isin(hd)])], list(S.TISSUES), name=lab)
        s = host.index.intersection(dip.index)
        phi = np.angle(host.loc[s, "x"] + 1j * host.loc[s, "y"])
        d = W * dip.loc[s].to_numpy()
        R, Rm = np.exp(1j * (phi - d)).mean(), np.exp(1j * (phi + d)).mean()
        so = host.index.intersection(don.index[don.tod.notna()])
        ph_t = np.angle(host.loc[so, "x"] + 1j * host.loc[so, "y"])
        Rt = np.exp(1j * (ph_t - W * don.tod.loc[so].to_numpy(float))).mean()
        rows.append(dict(group=lab, n_hosts=len(host), n_with_time_of_death=len(so), n_with_dip=len(s),
                         R_host_vs_time_of_death=abs(Rt), R_host_vs_dip=abs(R), R_host_vs_dip_mirrored=abs(Rm)))
    return pd.DataFrame(rows)


def d6b_dip(pb, don, gc, dip) -> pd.DataFrame:
    rows = []
    for lab, hd in (("fast deaths", (1, 2)), ("ventilator", (0,))):
        host = S.host_states(pb[pb.donor.isin(don.index[don.hardy.isin(hd)])], list(S.TISSUES), name=lab)
        s = host.index.intersection(dip.index).intersection(gc.index)
        f, g = host.loc[s, ["x", "y"]].to_numpy(), gc.loc[s].to_numpy()
        d = W * dip.loc[s].to_numpy()
        t = W * don.tod.reindex(s).to_numpy(float)
        ok = np.isfinite(t)
        ones = np.ones(len(g))
        phi = np.angle(f[:, 0] + 1j * f[:, 1])
        gz = g - g.mean()
        fr = f - np.outer(gz, np.linalg.lstsq(gz[:, None], f - f.mean(0), rcond=None)[0][0])
        cd = np.c_[np.cos(d), np.sin(d)]
        cdr = cd - np.outer(gz, np.linalg.lstsq(gz[:, None], cd - cd.mean(0), rcond=None)[0][0])
        terc = pd.qcut(g, 3, labels=False)
        row = dict(group=lab, n=len(s), r2_gc_on_host=r2(np.c_[ones, f], g),
                   r2_gc_on_dip=r2(np.c_[ones, np.cos(d), np.sin(d)], g),
                   r2_gc_on_time_of_death=r2(np.c_[ones[ok], np.cos(t[ok]), np.sin(t[ok])], g[ok]),
                   R_host_dip=float(abs(np.exp(1j * (phi - d)).mean())),
                   R_host_dip_after_removing_gc=float(abs(np.exp(1j * (np.angle(fr[:, 0] + 1j * fr[:, 1])
                                                                       - np.angle(cdr[:, 0] + 1j * cdr[:, 1]))).mean())))
        for k in range(3):
            row[f"R_host_dip_gc_tertile{k + 1}"] = float(abs(np.exp(1j * (phi[terc == k] - d[terc == k])).mean()))
            m = (terc == k) & ok
            row[f"R_host_tod_gc_tertile{k + 1}"] = float(abs(np.exp(1j * (phi[m] - t[m])).mean()))
        rows.append(row)
    return pd.DataFrame(rows)


def d6c_profile(don, gc, dip) -> pd.DataFrame:
    rows = []
    for lab, hd in (("fast deaths", (1, 2)), ("ventilator", (0,))):
        s = pd.DataFrame(dict(gc=gc, dip=dip.reindex(gc.index))).dropna()
        s = s[s.index.isin(don.index[don.hardy.isin(hd)])]
        bins = (s.dip % 24 // 4).astype(int)
        for b, x in s.groupby(bins):
            rows.append(dict(group=lab, quantity="gc_by_dip_bin_4h", bin_start_h=4 * b, n=len(x), mean_gc=x.gc.mean()))
        terc = pd.qcut(s.gc, 3, labels=False)
        for k, x in s.groupby(terc):
            z = np.exp(1j * W * x.dip).mean()
            rows.append(dict(group=lab, quantity="dip_phase_by_gc_tertile", bin_start_h=k + 1, n=len(x),
                             mean_phase_h=float(np.angle(z) / W % 24), R=float(abs(z))))
    return pd.DataFrame(rows)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    dip = read_dip()
    n_perm, n_draw = (10, 100) if run["quick"] else (50, 2000)
    CK = GC.CK
    pb, don = gtex.load(tissues=S.TISSUES, hardy=(0, 1, 2), genes=CK + list(S.GC), depth="library_size", require_tod=False)
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[:args.subsample]
        pb, don = pb[pb.donor.isin(keep)].reset_index(drop=True), don.loc[don.index.isin(keep)]
    pd.set_option("display.width", 250)

    ag = agreement(pb[["donor", "compartment", "depth", *CK]], don, dip)
    ag.to_csv(out / "diag5_dip_agreement.csv", index=False)
    print(ag.round(3).to_string(index=False), flush=True)

    # D5b: DIP as the internal time of ventilated donors
    vd = set(don.index[don.hardy == 0]) & set(dip.index)
    fv = pd.DataFrame(dict(x=np.cos(W * dip.loc[sorted(vd)]), y=np.sin(W * dip.loc[sorted(vd)])), index=sorted(vd))
    targets = [t for t in S.TISSUES if t != "Pituitary"]
    res = Parallel(n_jobs=args.jobs)(delayed(IE.d5_tissue)(t, fv, 50 + i, n_perm) for i, t in enumerate(targets))
    summ = pd.DataFrame([r[0] for r in res])
    gen = pd.concat([r[1] for r in res], ignore_index=True)
    summ.to_csv(out / "diag5_dip_summary.csv", index=False)
    gen.to_csv(out / "diag5_dip_genes.csv.gz", index=False)
    print(summ.round(2).to_string(index=False), flush=True)
    hp = out / "diag5_genes.csv.gz"
    if hp.exists():
        h = pd.read_csv(hp)
        m = h.merge(gen, on=["tissue", "gene"], suffixes=("_host", "_dip"))
        m = m[~m.frame_host]
        rows = []
        for t, d in m.groupby("tissue"):
            a_, b_ = (d.p_int_host < 1e-3) & (d.p_ext_host > 0.1), (d.p_int_dip < 1e-3) & (d.p_ext_dip > 0.1)
            both = a_ & b_
            dphi = np.exp(1j * W * (d.int_h_host[both] - d.int_h_dip[both])).mean() if both.sum() else np.nan
            rows.append(dict(tissue=t, host_int_only=int(a_.sum()), dip_int_only=int(b_.sum()), overlap=int(both.sum()),
                             jaccard=both.sum() / max((a_ | b_).sum(), 1),
                             phase_agree_R=abs(dphi) if both.sum() else np.nan,
                             phase_offset_h=float(np.angle(dphi) / W % 24) if both.sum() else np.nan))
        pd.DataFrame(rows).to_csv(out / "diag5_dip_overlap.csv", index=False)
    else:
        print("diag5_genes.csv.gz absent (run compute_internal_external.py): overlap not written")

    # D6, DIP arm
    g = gen[gen.tissue.isin(GC.TARGETS)].drop(columns=["ext_h"]).rename(columns={"amp_int": "amp_int"})
    cls = GC.class_tests(g, n_draw, 0).assign(arm="dip")
    cls.to_csv(out / "diag6_classes_dip.csv", index=False)
    print(cls[["cls", "n", "tissues", "retained", "retained_null", "retained_p", "phase_kept", "phase_kept_p", "amp_rel",
               "amp_rel_p", "ext", "ext_p"]].round(3).to_string(index=False), flush=True)

    # D6b, D6c part 1
    gc = GC.gc_scores(pb)
    cols = ["donor", "compartment", "depth", *CK]
    b = d6b_dip(pb[cols], don, gc, dip)
    b.to_csv(out / "diag6b_gc_dip.csv", index=False)
    pr = d6c_profile(don, gc, dip)
    pr.to_csv(out / "diag6c_gc_dip_profile.csv", index=False)
    print(b.round(3).to_string(index=False), "\n", pr.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
