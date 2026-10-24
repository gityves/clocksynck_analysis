"""A111 / A112 / A129 (WP3, feeding and genetic controls): scoring.

Registered in paper/docs/registered_experiments.md (section "A111, A112, A129 (WP3 ...)" and its WP3 note). Order:
  1. A111 timed truth and informativeness (Manella 2021: Alb-Cre controls ad libitum (AL) vs daytime feeding (DF), liver /
     WAT / lung; Xin 2021 replication: night- (NRF) vs day-restricted (DRF) feeding, liver / heart / kidney / VAT).
     Truth change = timed clock-plane phase (DF - AL) per tissue; informative if the liver shift relative to the other
     tissues is >= 3 h (|.|); otherwise H111.1 / H111.2 are reported as uninformative.
  2. H111.1 offsets: Fisher-Lee correlation of the untimed change in tissue offsets with the timed change over tissues,
     exact permutation p over all tissue permutations (one-sided). H111.2 coupling: liver coupling to the non-liver host,
     DF minus AL (d < 0, one-sided, labels permuted within ZT), direction of the timed liver-vs-rest R^2 change; with
     specificity companions on the same permutations.
  3. H111.3 (= A129 primary): Manella BLKO vs Alb-Cre, ad libitum: liver coupling to the WAT + lung host at the common
     spread d*, one-sided, labels permuted within ZT (1,000); companions WAT on lung and lung on WAT (same permutations).
  4. A112 (Greenwell 2019): lock ratio arrhythmic / ad libitum of the mouse orthologues of the GTEx A81b cue and locked sets.
  5. Descriptive clock-loss readouts (BLKO liver offset; WAT-lung coherence; VIP-Cre Bmal1 KO; whole-body Bmal1 KO).

    python compute_feeding_controls.py [--subsample N] [--only a111 a112 clockloss]
`--subsample N` is a quick end-to-end check (fewer permutations / bootstrap draws, outputs in output/subsample/; the
mice are all kept: the designs are small and balanced). Outputs: paper/exploratory/A111_feeding_controls/output/.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.dataset import _labels  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

import _common as C  # noqa: E402

warnings.filterwarnings("ignore", message="no anchor compartments given")  # VIP-Cre: two tissues, free host
FOLDER = "A111_feeding_controls"
CUE_REF = HERE.parents[2] / "results" / "figure_3" / "output" / "gtex_cues" / "a81b_cue_reference.csv"
INFORMATIVE_H = 3.0


def sizes(quick: bool) -> dict:
    return dict(n_perm=50 if quick else 1000, n_boot=20 if quick else 200, n_gate=100 if quick else 1000,
                n_timed=200 if quick else 2000, n_boot_timed=20 if quick else 200, n_gene_perm=1000 if quick else 10000)


# ---------------------------------------------------------------------------------------------------------------------
# A111 shared pieces
# ---------------------------------------------------------------------------------------------------------------------

def arm_dataset(pb, donors, arm_col, arm, anchors, wf, name):
    ids = donors.index[donors[arm_col] == arm]
    return cs.Dataset(pb[pb.donor.isin(ids)], anchors=anchors, donors=donors, waveform=wf, name=name)


def timed_block(ds, hours, tag, n_timed, n_boot_timed, seed) -> tuple[pd.DataFrame, pd.Series, list]:
    tl = ds.timed_loadings(hours)
    ht = ds.timed(hours, n_perm=n_timed, seed=seed)
    return tl.assign(**tag), ht, C.boot_timed(ds, hours, n_boot_timed, seed)


def shift_stats(tl0, tl1, target):
    """Timed change (case minus reference arm) per tissue: phase shift (h), the target's shift relative to the other
    tissues, and the change of the target's timed R^2 (coupling to known time) relative to the others'."""
    comps = list(tl0.index)
    d_off = pd.Series(S.wrap(tl1.clock_offset_h.reindex(comps).to_numpy() - tl0.clock_offset_h.to_numpy()), index=comps)
    d_r2 = tl1.coupling.reindex(comps) - tl0.coupling
    return d_off, C.rel_to_rest(d_off, target), float(d_r2[target] - d_r2.drop(target).mean()), d_r2


def informativeness(dataset, variant, contrast_name, tl0, tl1, b0, b1, target, extra=None):
    d_off, shift, r2_rel, d_r2 = shift_stats(tl0, tl1, target)
    bs, br = [], []
    for x0, x1 in zip(b0, b1):
        if target in x0.index and target in x1.index:
            _, s_, r_, _ = shift_stats(x0, x1.reindex(x0.index), target)
            bs.append(shift + S.wrap(s_ - shift))
            br.append(r_)
    lo, hi = (np.quantile(bs, 0.025), np.quantile(bs, 0.975)) if len(bs) > 10 else (np.nan, np.nan)
    rows = [dict(dataset=dataset, variant=variant, contrast=contrast_name, tissue=k, shift_h=d_off[k],
                 shift_rel_mean_h=float(C.circ_centre(d_off.to_numpy())[i]), d_timed_r2=d_r2[k])
            for i, k in enumerate(d_off.index)]
    summ = dict(dataset=dataset, variant=variant, contrast=contrast_name, target=target, shift_vs_rest_h=shift,
                shift_lo=lo, shift_hi=hi, informative=bool(abs(shift) >= INFORMATIVE_H),
                d_timed_r2_vs_rest=r2_rel, r2_lo=np.quantile(br, 0.025) if len(br) > 10 else np.nan,
                r2_hi=np.quantile(br, 0.975) if len(br) > 10 else np.nan, n_boot=len(bs), **(extra or {}))
    return pd.DataFrame(rows), summ, d_off


def gate_row(ds, label, n_gate, **kw):
    g = ds.gate(n_perm=n_gate, seed=0, **kw)
    return dict(set=label, scope=("all" if not kw else "->".join(["+".join(kw["targets"]), "+".join(kw["hosts"])])),
                clock_share=g.clock_share, null_mean=g.share_null_mean, p_pairing=g.p_pairing,
                p_acrophase=g.p_acrophase, passes=bool(g.passes), n_donors=g.n_donors)


def untimed_offsets(ds) -> pd.Series:
    f = ds.fit(n_boot=0, gate=False)
    return f.table.offset_h


def h111_1(untimed0, untimed1, d_off_truth, label) -> dict:
    comps = list(untimed0.index)
    du = pd.Series(S.wrap(untimed1.reindex(comps).to_numpy() - untimed0.to_numpy()), index=comps)
    du_c, tr_c = C.circ_centre(du.to_numpy()), C.circ_centre(d_off_truth.reindex(comps).to_numpy())
    r = C.exact_circ_corr(du_c, tr_c)
    detail = {f"untimed_{k}": v for k, v in zip(comps, du_c)} | {f"timed_{k}": v for k, v in zip(comps, tr_c)}
    return dict(contrast=label, **r, **detail)


def public_contrast(ds, group, levels, targets, companion, blocks, sz, seed=0):
    return ds.contrast(group, levels=levels, targets=targets, companion=companion, blocks=blocks, n_perm=sz["n_perm"],
                       n_boot=sz["n_boot"], seed=seed)


def one_sided(c, target):
    """Contrast row of `target` with one-sided (less) p of d and d_star from the kept permutation draws."""
    r = c.table.loc[target].to_dict()
    r["p_less_d"] = S.perm_p(c.null["d"][target], c.observed["d"][target], "less")
    r["p_less_d_star"] = S.perm_p(c.null["d_star"][target], c.observed["d_star"][target], "less")
    return r


def engine_for(ds, group, levels, blocks, tissue_hosts: dict, sz, seed=0):
    """Contrast engine on the donors of `group`; tissue_hosts: key -> (target, host spec) with host spec = list of
    anchor tissues (>= 2) or a single tissue name."""
    donors = ds.donors
    lab, levels = _labels(ds._donor_values(group), levels, set(ds.table.resection))
    blk = ds._blocks(blocks, lab)
    d = ds.data.take(np.flatnonzero(np.isin(ds.data.resection, lab.index)))
    keys = {}
    for name, (target, spec) in tissue_hosts.items():
        h = C.multi_host(d, spec, target) if isinstance(spec, (list, tuple)) else C.pseudo_host(d, spec)
        keys[name] = (target, h)
    return C.contrast_engine(d, keys, lab, blk, levels, sz["n_perm"], sz["n_boot"], seed), d, lab


def check_against_public(res, contrast, key, target):
    a, b = res["null"][key]["d"].to_numpy(), contrast.null["d"][target].to_numpy()
    ok = np.allclose(a, b, equal_nan=True)
    e = res["obs"].loc[key, "d"] - contrast.observed["d"][target]
    print(f"   engine == public contrast for {target}: null identical {ok}, obs diff {e:.2e}", flush=True)
    return dict(null_identical=bool(ok), obs_diff=float(e))


def pair_coherence(ds, a, b, group, levels, blocks, sz, seed=0, tag=""):
    """|C| and offset (a minus b, h) of two tissues per group, difference case - reference, bootstrap interval and
    permutation p (labels within blocks), descriptive."""
    lab, levels = _labels(ds._donor_values(group), levels, set(ds.table.resection))
    blk = ds._blocks(blocks, lab)
    d = ds.data.take(np.flatnonzero(np.isin(ds.data.resection, lab.index)))
    both = pd.Index(np.intersect1d(d.donors(a), d.donors(b)))
    lab = lab.reindex(lab.index.intersection(both))
    blk = None if blk is None else blk.reindex(lab.index)

    def stat(labels):
        z = [C.pair_corr(d, a, b, labels.index[labels == g]) for g in levels]
        off = S.wrap(np.angle(z[1] / z[0]) * C.HOURS) if abs(z[0]) > 0 and abs(z[1]) > 0 else np.nan
        return abs(z[0]), abs(z[1]), off

    o = stat(lab)
    rng = np.random.default_rng(seed)
    nd, no = [], []
    for _ in range(sz["n_perm"]):
        r = stat(S.permute_within(lab, blk, rng))
        nd.append(r[1] - r[0])
        no.append(r[2])
    brng = np.random.default_rng(seed + 1)
    don = {g: np.array(lab.index[lab == g]) for g in levels}
    ia = pd.Series(np.flatnonzero(d.compartment == a), index=d.resection[d.compartment == a])
    ib = pd.Series(np.flatnonzero(d.compartment == b), index=d.resection[d.compartment == b])
    boot = []
    for _ in range(sz["n_boot"]):
        m = []
        for g in levels:
            dr = brng.choice(don[g], len(don[g]))
            x, y = ia[dr].to_numpy(), ib[dr].to_numpy()
            m.append(abs(S.ccorr(d.beta[x], d.beta[y], np.trace(d.V[x], axis1=1, axis2=2),
                                 np.trace(d.V[y], axis1=1, axis2=2))))
        boot.append(m[1] - m[0])
    boot = np.array(boot, float)
    return dict(pair=f"{a}-{b}", tag=tag, group=group, ref=levels[0], case=levels[1], n_ref=len(don[levels[0]]),
                n_case=len(don[levels[1]]), C_ref=o[0], C_case=o[1], dC=o[1] - o[0],
                dC_lo=np.nanquantile(boot, 0.025), dC_hi=np.nanquantile(boot, 0.975),
                p_two=S.perm_p(np.array(nd), o[1] - o[0]), p_less=S.perm_p(np.array(nd), o[1] - o[0], "less"),
                offset_ref_to_case_h=o[2], p_offset=S.perm_p(np.array(no), o[2], tie=0.0))


# ---------------------------------------------------------------------------------------------------------------------
# A111 Manella (+ H111.3 / A129, clock-loss readouts on the same data)
# ---------------------------------------------------------------------------------------------------------------------

def run_manella(out: Path, sz: dict, results: dict):
    pb, donors = data.load("manella2021")
    qc = C.gene_qc(pb, "manella2021")
    results["qc"].append(qc)
    fail = sorted(qc[~qc.kept].gene.unique())
    print("Manella QC: genes failing median count >= 5 in some tissue:", qc[~qc.kept][["compartment", "gene"]].values.tolist(),
          flush=True)
    # sensitivity: CIART kept (only the genes absent from a whole genotype x tissue file are dropped)
    absent = sorted(qc[qc.n_missing > 0].gene.unique())
    variants = {"primary": C.waveform_without(fail), "keep_failing": C.waveform_without(absent)}
    anchors = ["WAT", "Lung"]
    arms = ["AlbCre_AL", "AlbCre_DF", "BLKO_AL", "BLKO_DF"]
    for vname, wf in variants.items():
        print(f"== Manella, variant {vname} ({len(wf.genes)} genes)", flush=True)
        ds = {a: arm_dataset(pb, donors, "arm", a, anchors, wf, a) for a in arms}
        tl, boots = {}, {}
        for a in arms:
            t, ht, boots[a] = timed_block(ds[a], "zt", dict(dataset="manella2021", variant=vname, arm=a), sz["n_timed"],
                                          sz["n_boot_timed"], 0)
            tl[a] = t
            results["timed"].append(t.reset_index())
            results["host_timed"].append(pd.DataFrame([ht.to_dict() | dict(dataset="manella2021", variant=vname, arm=a)]))
            g = [gate_row(ds[a], a, sz["n_gate"]),
                 gate_row(ds[a], a, sz["n_gate"], targets=["Liver"], hosts=anchors),
                 gate_row(ds[a], a, sz["n_gate"], targets=["WAT"], hosts=["Lung"])]
            results["gates"].append(pd.DataFrame(g).assign(dataset="manella2021", variant=vname))
        # ---- informativeness (primary: Alb-Cre controls AL -> DF; BLKO descriptive) -------------------------------------
        trs, truth = {}, {}
        for gen in ("AlbCre", "BLKO"):
            a0, a1 = f"{gen}_AL", f"{gen}_DF"
            rows, summ, d_off = informativeness("manella2021", vname, f"{gen}: AL->DF", tl[a0], tl[a1], boots[a0], boots[a1],
                                                "Liver")
            results["truth_rows"].append(rows)
            results["informativeness"].append(pd.DataFrame([summ]))
            truth[gen] = (summ, d_off)
        s_albcre, d_off = truth["AlbCre"]
        informative = s_albcre["informative"]
        print(f"   timed liver shift vs WAT+lung (AlbCre AL->DF): {s_albcre['shift_vs_rest_h']:+.2f} h "
              f"[{s_albcre['shift_lo']:+.2f}, {s_albcre['shift_hi']:+.2f}] -> informative {informative}; "
              f"BLKO {truth['BLKO'][0]['shift_vs_rest_h']:+.2f} h", flush=True)
        # ---- H111.1 ------------------------------------------------------------------------------------------------------
        u0, u1 = untimed_offsets(ds["AlbCre_AL"]), untimed_offsets(ds["AlbCre_DF"])
        r1 = h111_1(u0, u1, d_off, "AlbCre: AL->DF")
        r1.update(dataset="manella2021", variant=vname, informative=informative,
                  verdict=("uninformative" if not informative else
                           "PASS" if (r1["r_fl"] > 0 and r1["p_exact"] < 0.05) else "FAIL"))
        results["h111_1"].append(pd.DataFrame([r1]))
        print(f"   H111.1 r_FL {r1['r_fl']:+.3f}, exact p {r1['p_exact']:.3f} (min attainable {r1['p_min']:.3f}) -> {r1['verdict']}",
              flush=True)
        # ---- H111.2: Alb-Cre, AL vs DF -----------------------------------------------------------------------------------
        pool = cs.Dataset(pb[pb.donor.isin(donors.index[donors.genotype == "AlbCre"])], anchors=anchors, donors=donors,
                          waveform=wf, name="AlbCre")
        c = public_contrast(pool, "feeding", ("AL", "DF"), ["Liver"], None, "zt", sz)
        res, d_, lab = engine_for(pool, "feeding", ("AL", "DF"), "zt",
                                  {"Liver": ("Liver", anchors), "WAT": ("WAT", "Lung"), "Lung": ("Lung", "WAT")}, sz)
        chk = check_against_public(res, c, "Liver", "Liver")
        summ = C.summarise_engine(res, ["WAT", "Lung"], "Liver").assign(dataset="manella2021", variant=vname,
                                                                          contrast="AlbCre: AL->DF (H111.2)")
        pub = one_sided(c, "Liver")
        gate_ok = bool(pd.concat(results["gates"][-4:]).query("set == 'AlbCre_AL' and scope == 'Liver->WAT+Lung'").passes.iloc[0])
        summ.loc[summ.key == "Liver", ["public_d", "public_p_less_d", "public_d_star", "public_p_less_d_star"]] = [
            pub["d"], pub["p_less_d"], pub["d_star"], pub["p_less_d_star"]]
        v = ("uninformative" if not informative else
             "PASS" if (pub["d"] < 0 and pub["p_less_d"] < 0.05 and truth["AlbCre"][0]["d_timed_r2_vs_rest"] < 0) else "FAIL")
        summ["verdict"] = np.where(summ.key == "Liver", v, "")
        summ["gate_liver_AlbCre_AL_passes"] = gate_ok
        summ["timed_dR2_liver_vs_rest"] = truth["AlbCre"][0]["d_timed_r2_vs_rest"]
        summ["check_null_identical"] = chk["null_identical"]
        results["h111_2"].append(summ)
        print(f"   H111.2 liver d {pub['d']:+.3f} (p_less {pub['p_less_d']:.3f}), d* {pub['d_star']:+.3f} "
              f"(p_less {pub['p_less_d_star']:.3f}); timed dR2 vs rest {truth['AlbCre'][0]['d_timed_r2_vs_rest']:+.3f} -> {v}",
              flush=True)
        # ---- H111.3: AL arms, BLKO vs Alb-Cre ----------------------------------------------------------------------------
        for feed, role in (("AL", "H111.3 (primary)"), ("DF", "descriptive")):
            pool = cs.Dataset(pb[pb.donor.isin(donors.index[donors.feeding == feed])], anchors=anchors, donors=donors,
                              waveform=wf, name=feed)
            c = public_contrast(pool, "genotype", ("AlbCre", "BLKO"), ["Liver"], None, "zt", sz)
            res, d_, lab = engine_for(pool, "genotype", ("AlbCre", "BLKO"), "zt",
                                      {"Liver": ("Liver", anchors), "WAT": ("WAT", "Lung"), "Lung": ("Lung", "WAT")}, sz)
            chk = check_against_public(res, c, "Liver", "Liver")
            summ = C.summarise_engine(res, ["WAT", "Lung"], "Liver").assign(dataset="manella2021", variant=vname,
                                                                              contrast=f"{feed}: BLKO vs AlbCre ({role})")
            pub = one_sided(c, "Liver")
            summ.loc[summ.key == "Liver", ["public_d", "public_p_less_d", "public_d_star", "public_p_less_d_star",
                                           "public_offset_diff_h", "public_p_offset"]] = [
                pub["d"], pub["p_less_d"], pub["d_star"], pub["p_less_d_star"], pub["offset_diff_h"], pub["p_offset"]]
            spec = {}
            for comp in ("WAT", "Lung"):
                row = summ[summ.key == "Liver"].iloc[0]
                spec[comp] = bool(row["d_star_minus_" + comp] < 0)
            summ["liver_drop_exceeds_both_companions"] = bool(all(spec.values()))
            if feed == "AL":
                gg = pd.concat(results["gates"][-4:]).query("scope == 'Liver->WAT+Lung' or scope == 'WAT->Lung'")
                gl = {r.set: r.passes for r in gg.itertuples() if r.scope == "Liver->WAT+Lung"}
                gw = {r.set: r.passes for r in gg.itertuples() if r.scope == "WAT->Lung"}
                alt = (not gl["BLKO_AL"]) and gl["AlbCre_AL"] and gw["BLKO_AL"]
                v = "PASS" if (pub["d_star"] < 0 and pub["p_less_d_star"] < 0.05) else "FAIL"
                summ["verdict"] = np.where(summ.key == "Liver", v, "")
                summ["alt_gate_criterion_met"] = bool(alt)
                summ["gate_liver_BLKO_AL_passes"] = bool(gl["BLKO_AL"])
                summ["gate_liver_AlbCre_AL_passes"] = bool(gl["AlbCre_AL"])
                summ["gate_WAT_Lung_BLKO_AL_passes"] = bool(gw["BLKO_AL"])
                print(f"   H111.3 liver d* {pub['d_star']:+.3f} (p_less {pub['p_less_d_star']:.3f}); d {pub['d']:+.3f} "
                      f"(p_less {pub['p_less_d']:.3f}); companions d* WAT {res['obs'].loc['WAT','d_star']:+.3f}, "
                      f"Lung {res['obs'].loc['Lung','d_star']:+.3f}; alt gate criterion {bool(alt)} -> {v}", flush=True)
            summ["check_null_identical"] = chk["null_identical"]
            results["h111_3"].append(summ)
            # clock-loss: liver offset to the host, BLKO vs AlbCre, untimed; timed from the arms
            gen0, gen1 = "AlbCre", "BLKO"
            d_off_t, shift_t, _, _ = shift_stats(tl[f"{gen0}_{feed}"], tl[f"{gen1}_{feed}"], "Liver")
            results["clockloss"].append(pd.DataFrame([dict(
                dataset="manella2021", variant=vname, feeding=feed, readout="liver offset to WAT+lung host, BLKO minus AlbCre",
                untimed_h=pub["offset_diff_h"], untimed_p=pub["p_offset"], timed_h=shift_t,
                timed_wat_h=d_off_t["WAT"], timed_lung_h=d_off_t["Lung"], timed_liver_h=d_off_t["Liver"])]))
            if vname == "primary":
                results["coherence"].append(pd.DataFrame([pair_coherence(pool, "WAT", "Lung", "genotype", ("AlbCre", "BLKO"),
                                                                       "zt", sz, tag=f"Manella {feed}")]))
        if vname == "primary":
            for gen in ("AlbCre", "BLKO"):
                pool = cs.Dataset(pb[pb.donor.isin(donors.index[donors.genotype == gen])], anchors=anchors, donors=donors,
                                  waveform=wf, name=gen)
                c = public_contrast(pool, "feeding", ("AL", "DF"), ["Liver"], None, "zt", sz)
                r = one_sided(c, "Liver")
                results["h111_2_descr"].append(pd.DataFrame([dict(dataset="manella2021", genotype=gen, **{
                    k: r[k] for k in ("c_AL", "c_DF", "d", "lo", "hi", "d_star", "p_less_d", "p_less_d_star", "offset_diff_h",
                                      "p_offset")})]))


# ---------------------------------------------------------------------------------------------------------------------
# A111 replication: Xin 2021
# ---------------------------------------------------------------------------------------------------------------------

def run_xin(out: Path, sz: dict, results: dict):
    pb, donors = data.load("xin2021")
    qc = C.gene_qc(pb, "xin2021")
    results["qc"].append(qc)
    fail = sorted(qc[~qc.kept].gene.unique())
    wf = C.waveform_without(fail)
    print(f"== Xin 2021 ({len(wf.genes)} genes; QC failures {fail})", flush=True)
    anchors = ["heart", "kidney", "VAT"]
    arms = ["NRF", "DRF"]
    ds = {a: arm_dataset(pb, donors, "arm", a, anchors, wf, a) for a in arms}
    tl, boots = {}, {}
    for a in arms:
        t, ht, boots[a] = timed_block(ds[a], "zt_fold", dict(dataset="xin2021", variant="primary", arm=a), sz["n_timed"],
                                      sz["n_boot_timed"], 0)
        tl[a] = t
        results["timed"].append(t.reset_index())
        results["host_timed"].append(pd.DataFrame([ht.to_dict() | dict(dataset="xin2021", variant="primary", arm=a)]))
        results["gates"].append(pd.DataFrame([gate_row(ds[a], a, sz["n_gate"]),
                                              gate_row(ds[a], a, sz["n_gate"], targets=["liver"], hosts=anchors)])
                                .assign(dataset="xin2021", variant="primary"))
    rows, summ, d_off = informativeness("xin2021", "primary", "NRF->DRF", tl["NRF"], tl["DRF"], boots["NRF"], boots["DRF"],
                                        "liver")
    results["truth_rows"].append(rows)
    results["informativeness"].append(pd.DataFrame([summ]))
    informative = summ["informative"]
    print(f"   timed liver shift vs rest (NRF->DRF): {summ['shift_vs_rest_h']:+.2f} h [{summ['shift_lo']:+.2f}, "
          f"{summ['shift_hi']:+.2f}] -> informative {informative}", flush=True)
    u0, u1 = untimed_offsets(ds["NRF"]), untimed_offsets(ds["DRF"])
    r1 = h111_1(u0, u1, d_off, "NRF->DRF")
    r1.update(dataset="xin2021", variant="primary", informative=informative,
              verdict=("uninformative" if not informative else "PASS" if (r1["r_fl"] > 0 and r1["p_exact"] < 0.05) else "FAIL"))
    results["h111_1"].append(pd.DataFrame([r1]))
    print(f"   H111.1 r_FL {r1['r_fl']:+.3f}, exact p {r1['p_exact']:.4f} (min {r1['p_min']:.4f}) -> {r1['verdict']}", flush=True)
    pool = cs.Dataset(pb, anchors=anchors, donors=donors, waveform=wf, name="xin")
    comps = ["liver", "heart", "kidney", "VAT"]
    c = pool.contrast("arm", levels=("NRF", "DRF"), targets=comps, companion=["heart", "kidney", "VAT"], blocks="zt_fold",
                      n_perm=sz["n_perm"], n_boot=sz["n_boot"], seed=0)
    r = one_sided(c, "liver")
    ct = c.companion
    v = ("uninformative" if not informative else
         "PASS" if (r["d"] < 0 and r["p_less_d"] < 0.05 and summ["d_timed_r2_vs_rest"] < 0) else "FAIL")
    row = dict(dataset="xin2021", variant="primary", contrast="NRF->DRF (H111.2)", **{k: r[k] for k in (
        "donors_NRF", "donors_DRF", "c_NRF", "c_DRF", "d", "lo", "hi", "d_star", "p_less_d", "p_less_d_star", "offset_diff_h",
        "p_offset")}, timed_dR2_liver_vs_rest=summ["d_timed_r2_vs_rest"], verdict=v)
    for k in ("d_minus_companion", "p", "d_star_minus_companion", "p_star"):
        row["companion_" + k] = ct.iloc[0][k]
    row["gate_liver_NRF_passes"] = bool(pd.concat(results["gates"][-2:]).query("set == 'NRF' and scope == 'liver->heart+kidney+VAT'").passes.iloc[0])
    results["h111_2"].append(pd.DataFrame([row]))
    other = c.table.loc[["heart", "kidney", "VAT"], ["c_NRF", "c_DRF", "d", "d_star"]].assign(
        p_less_d=[S.perm_p(c.null["d"][k], c.observed["d"][k], "less") for k in ["heart", "kidney", "VAT"]])
    other.rename_axis("tissue").reset_index().assign(dataset="xin2021").to_csv(out / "a111_xin_other_tissues.csv", index=False)
    print(f"   H111.2 liver d {r['d']:+.3f} (p_less {r['p_less_d']:.3f}), d* {r['d_star']:+.3f} (p_less {r['p_less_d_star']:.3f}); "
          f"timed dR2 vs rest {summ['d_timed_r2_vs_rest']:+.3f}; companion d-d_comp {ct.iloc[0]['d_minus_companion']:+.3f} -> {v}",
          flush=True)


# ---------------------------------------------------------------------------------------------------------------------
# A112: Greenwell 2019, lock ratio of the mouse orthologues of the A81b cue and locked sets
# ---------------------------------------------------------------------------------------------------------------------

def run_a112(out: Path, sz: dict, results: dict):
    mod = data._builder("greenwell2019")
    pb, donors = data.load("greenwell2019")
    cnt, meta = mod.counts()
    qc = C.gene_qc(pb, "greenwell2019")
    results["qc"].append(qc)
    fail = sorted(qc[~qc.kept].gene.unique())
    wf = C.waveform_without(fail)
    print(f"== Greenwell 2019 (liver; clock genes kept {len(wf.genes)}, dropped {fail})", flush=True)
    # host: the liver's own clock-plane state of every library (one projector over all 72 libraries)
    x = pb.set_index("donor")
    L = projector(x, wf, counts=True)
    beta, _ = states(x, L, wf, depth=x.depth.to_numpy(float))
    host = pd.DataFrame(beta, index=x.index, columns=["x", "y"])
    meta = meta.reset_index().set_index("donor")
    arm_of = {"ad_libitum": "LB", "arrhythmic": "AR", "arrhythmic_unadapted": "DR", "night_restricted": "NR"}
    meta["arm2"] = meta.arm.map(arm_of)
    meta = meta.loc[host.index]
    host["arm"], host["zt"] = meta.arm2.to_numpy(), meta.zt.to_numpy(float)
    # descriptive timed state per arm
    rows = []
    for a, h in host.groupby("arm"):
        t = h.zt.to_numpy() * C.W
        D = np.c_[np.ones(len(h)), np.cos(t), np.sin(t)]
        bx = np.linalg.lstsq(D, h.x.to_numpy(), rcond=None)[0]
        by = np.linalg.lstsq(D, h.y.to_numpy(), rcond=None)[0]
        amp = complex((bx[1] + by[2]) / 2, (by[1] - bx[2]) / 2)  # state = a e^{i W t}
        rows.append(dict(arm=a, n=len(h), state_amplitude=abs(amp), state_acrophase_h=float(np.mod(np.angle(amp), 2 * np.pi) / C.W)))
    pd.DataFrame(rows).to_csv(out / "a112_host_state_by_arm.csv", index=False)
    # orthologues of the A81b sets (upper-cased symbol, the A48 convention)
    ref = pd.read_csv(CUE_REF, index_col=0)
    sym = pd.Series(cnt.index.to_numpy(), index=[s.upper() for s in cnt.index])
    sym = sym[~sym.index.duplicated()]
    orth = pd.DataFrame(dict(set=ref.set, nt=ref.nt, human_lock=ref.lock), index=ref.index)
    orth["mouse_symbol"] = [sym.get(g) for g in orth.index]
    orth["mapped"] = orth.mouse_symbol.notna()
    # expression: median raw count >= 5 over the libraries of the two arms compared (per contrast), LB + AR here
    dep = pb.set_index("donor").depth.reindex(meta.index)
    lib_cols = meta.lib.to_numpy()
    cnt_d = cnt[lib_cols]
    cnt_d.columns = meta.index
    cpm = cnt_d.div(dep, axis=1) * 1e6
    y = np.log(cpm + 0.5)
    tw = host.zt.to_numpy() * C.W

    def loadings(Y: pd.DataFrame, arm: str):
        i = np.flatnonzero(host.arm.to_numpy() == arm)
        Yi = Y.iloc[:, i].to_numpy(float)
        Yi = Yi - Yi.mean(1, keepdims=True)
        f = host[["x", "y"]].to_numpy()[i]
        f = f - f.mean(0)
        B_u = np.linalg.lstsq(f, Yi.T, rcond=None)[0]          # untimed: coefficient on the host state
        t = np.c_[np.cos(tw[i]), np.sin(tw[i])]
        t = t - t.mean(0)
        B_t = np.linalg.lstsq(t, Yi.T, rcond=None)[0]          # timed: cosinor on ZT
        return np.hypot(*B_u), np.hypot(*B_t)

    gene_rows, res_rows, core_rows = [], [], []
    clock_syms = {g: sym.get(g) for g in wf.genes}
    for other in ("AR", "DR", "NR"):
        i_pair = np.flatnonzero(host.arm.isin(["LB", other]).to_numpy())
        med = cnt_d.iloc[:, i_pair].median(axis=1)
        expr = med >= C.EXPR_MIN
        au_o, at_o = loadings(y, other)
        au_l, at_l = loadings(y, "LB")
        with np.errstate(divide="ignore", invalid="ignore"):  # genes with no counts in an arm are not expressed
            lr = pd.DataFrame(dict(log_ratio_untimed=np.log(au_o / au_l), log_ratio_timed=np.log(at_o / at_l),
                                   amp_untimed_LB=au_l, amp_timed_LB=at_l, expressed=expr.to_numpy()), index=cnt.index)
        core = lr.loc[[s for s in clock_syms.values() if s in lr.index]]
        core_rows.append(dict(contrast=f"{other} vs LB", n_core=len(core), median_log_ratio_untimed=core.log_ratio_untimed.median(),
                              median_log_ratio_timed=core.log_ratio_timed.median()))
        o = orth[orth.mapped].copy()
        o = o.join(lr, on="mouse_symbol")
        o["expressed"] = o.expressed.fillna(False).astype(bool)
        use = o[o.expressed]
        if other == "AR":
            orth_out = orth.copy()
            orth_out["expressed_AR_LB"] = orth_out.mouse_symbol.map(lambda s: bool(lr.expressed.get(s, False)) if isinstance(s, str) else False)
            orth_out.to_csv(out / "a112_orthologues.csv")
            gene_tab = use.assign(contrast="AR vs LB")
            gene_tab.to_csv(out / "a112_genes_AR_vs_LB.csv.gz")
        rng = np.random.default_rng(112)
        for kind in ("untimed", "timed"):
            col = f"log_ratio_{kind}"
            cue, lock = use.loc[use.set == "cue", col].to_numpy(), use.loc[use.set == "locked", col].to_numpy()
            stat = float(np.median(cue) - np.median(lock))
            allv = np.r_[cue, lock]
            null = np.empty(sz["n_gene_perm"])
            for j in range(len(null)):
                p = rng.permutation(allv)
                null[j] = np.median(p[:len(cue)]) - np.median(p[len(cue):])
            brng = np.random.default_rng(1)
            boot = np.array([np.median(brng.choice(cue, len(cue))) - np.median(brng.choice(lock, len(lock)))
                             for _ in range(2000)])
            res_rows.append(dict(contrast=f"{other} vs LB", version=kind, primary=bool(other == "AR" and kind == "untimed"),
                                 n_cue_mapped=int(((o.set == "cue")).sum()), n_locked_mapped=int(((o.set == "locked")).sum()),
                                 n_cue=len(cue), n_locked=len(lock), median_cue=float(np.median(cue)),
                                 median_locked=float(np.median(lock)), diff_cue_minus_locked=stat,
                                 lo=float(np.quantile(boot, 0.025)), hi=float(np.quantile(boot, 0.975)),
                                 p_less=S.perm_p(null, stat, "less"), p_two=S.perm_p(null, stat),
                                 core_median=core[col].median(),
                                 verdict=("PASS" if (stat < 0 and S.perm_p(null, stat, "less") < 0.05) else "FAIL")
                                 if (other == "AR" and kind == "untimed") else ""))
    r = pd.DataFrame(res_rows)
    r.to_csv(out / "a112_h112_1.csv", index=False)
    pd.DataFrame(core_rows).to_csv(out / "a112_core_clock.csv", index=False)
    print(r.round(4).to_string(index=False), flush=True)
    print(f"   orthologues mapped: cue {int((orth.set=='cue').sum() and orth[orth.set=='cue'].mapped.sum())}/40, "
          f"locked {int(orth[orth.set=='locked'].mapped.sum())}/40", flush=True)


# ---------------------------------------------------------------------------------------------------------------------
# Descriptive clock-loss readouts: VIP-Cre Bmal1 KO (SCN knockout; liver-heart) and whole-body Bmal1 KO (Jouffe)
# ---------------------------------------------------------------------------------------------------------------------

def run_other_clockloss(out: Path, sz: dict, results: dict):
    # ---- VIP-Cre Bmal1 KO: liver and heart, CT0 / CT12, WT vs VIPKO within feeding ---------------------------------
    pb, donors = data.load("vipko_gse296764")
    qc = C.gene_qc(pb, "vipko_gse296764")
    results["qc"].append(qc)
    wf = C.waveform_without(sorted(qc[~qc.kept].gene.unique()))
    sub = pb[pb.compartment.isin(["liver", "heart"])]
    print(f"== VIP-Cre KO (liver + heart; mice with both: {sub.groupby('donor').compartment.nunique().eq(2).sum()})", flush=True)
    for feed in ("AL", "TRF"):
        ids = donors.index[donors.feeding == feed]
        ds = cs.Dataset(sub[sub.donor.isin(ids)], anchors=None, donors=donors, waveform=wf, name=f"vipko_{feed}")
        res, d_, lab = engine_for(ds, "genotype", ("WT", "VIPKO"), "ct", {"liver|heart": ("liver", "heart"),
                                                                           "heart|liver": ("heart", "liver")}, sz)
        t = C.summarise_engine(res).assign(dataset="vipko_gse296764", feeding=feed, role="descriptive")
        results["clockloss_coupling"].append(t)
        results["coherence"].append(pd.DataFrame([pair_coherence(ds, "liver", "heart", "genotype", ("WT", "VIPKO"), "ct", sz,
                                                                 tag=f"VIP-Cre {feed}")]))
        results["gates"].append(pd.DataFrame([gate_row(ds, f"vipko_{feed}", sz["n_gate"])]).assign(dataset="vipko_gse296764",
                                                                                                  variant="primary"))
        print(t[["key", "n_ref", "n_case", "c_ref", "c_case", "d", "d_star", "p_less_d_star", "offset_diff_h"]].round(3)
              .to_string(index=False), flush=True)
    # ---- whole-body Bmal1 KO (Jouffe): coupling only, single ZT -----------------------------------------------------
    pb, donors = data.load("jouffe2022_bmal1ko")
    qc = C.gene_qc(pb, "jouffe2022_bmal1ko")
    results["qc"].append(qc)
    wf = C.waveform_without(sorted(qc[~qc.kept].gene.unique()))
    print("== Jouffe 2022 whole-body Bmal1 KO: paired mice per arm", flush=True)
    have = pb.assign(v=1).pivot_table(index="donor", columns="compartment", values="v", aggfunc="size").notna()
    have["arm"] = (donors.bmal + "_" + donors.leptin).reindex(have.index)
    rows = []
    for a, b in (("liver", "muscle"), ("liver", "WAT"), ("muscle", "WAT")):
        n = have[have[a] & have[b]].groupby("arm").size()
        for arm in sorted(have.arm.unique()):
            rows.append(dict(dataset="jouffe2022_bmal1ko", pair=f"{a}-{b}", arm=arm, n_mice_with_both=int(n.get(arm, 0)),
                             estimable=bool(n.get(arm, 0) >= 5)))
    t = pd.DataFrame(rows)
    t.to_csv(out / "a111_jouffe_pair_counts.csv", index=False)
    ok = t.groupby("pair").estimable.apply(lambda x: x.reset_index(drop=True))
    print(t.pivot(index="pair", columns="arm", values="n_mice_with_both").to_string(), flush=True)
    print("   coupling contrast KO vs WT needs >= 5 mice with both tissues in each arm (Dataset min_donors): not met in any "
          "pair, so no coupling is estimated (descriptive readout not available)", flush=True)


# ---------------------------------------------------------------------------------------------------------------------

def save(results: dict, out: Path):
    names = dict(qc="a111_qc_genes", timed="a111_timed_truth", host_timed="a111_host_vs_time", gates="a111_gates",
                 truth_rows="a111_truth_shifts", informativeness="a111_informativeness", h111_1="a111_h111_1",
                 h111_2="a111_h111_2", h111_2_descr="a111_h111_2_descriptive", h111_3="a111_h111_3",
                 clockloss="a111_clockloss_liver_offset", clockloss_coupling="a111_clockloss_coupling",
                 coherence="a111_pair_coherence")
    for k, f in names.items():
        if results.get(k):
            pd.concat(results[k], ignore_index=True).to_csv(out / f"{f}.csv", index=False)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="*", default=["a111", "a112", "clockloss"])
    args = ap.parse_args(argv)
    run = cli.setup(args)
    sz = sizes(run["quick"])
    out = exploratory_dir(FOLDER, run["quick"])
    print("sizes", sz, "out", out, flush=True)
    results = {k: [] for k in ("qc", "timed", "host_timed", "gates", "truth_rows", "informativeness", "h111_1", "h111_2",
                               "h111_2_descr", "h111_3", "clockloss", "clockloss_coupling", "coherence")}
    if "a111" in args.only:
        run_manella(out, sz, results)
        save(results, out)
        run_xin(out, sz, results)
        save(results, out)
    if "clockloss" in args.only:
        run_other_clockloss(out, sz, results)
        save(results, out)
    if "a112" in args.only:
        run_a112(out, sz, results)
    save(results, out)
    print("done", flush=True)


if __name__ == "__main__":
    main()
