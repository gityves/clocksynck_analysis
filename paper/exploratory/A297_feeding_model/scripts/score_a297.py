"""A297 (registered in `paper/docs/registered_experiments.md`): does removing the feeding rhythm do what the additive cue
model, calibrated on inverted feeding, predicts?

Liver clock state per arm: one projector per dataset and tissue over all its arms (clock genes passing the A111 QC);
per arm, OLS of the state on cos, sin of ZT gives the complex amplitude A (state = mu + A e^{iWt}). Additive cues:
A_AL = d + f, A_inv = d - f, so rho = A_inv / A_AL and the removal prediction A_removed / A_AL = (1 + rho) / 2.
Calibration: Manella Alb-Cre controls (AL vs DF; Xin NRF vs DRF as a second calibration). Test: Greenwell arrhythmic
(adapted) vs ad libitum; night-restricted vs ad libitum as the internal control. Mouse bootstrap within arm (2000).

    python score_a297.py
Outputs (`../output/`): arms.csv, tests.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[2] / "A111_feeding_controls" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _common as C  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from common import data  # noqa: E402

OUT = HERE.parents[1] / "output"
W = 2 * np.pi / 24
N_BOOT = 2000


def liver_states(name, tissue, arm_col):
    pb, donors = data.load(name)
    qc = C.gene_qc(pb, name)
    wf = C.waveform_without(sorted(qc[~qc.kept].gene.unique()))
    x = pb[pb.compartment == tissue].set_index("donor")
    L = projector(x[list(wf.genes)], wf, counts=True)
    b, _ = states(x, L, wf, depth=x.depth.to_numpy(float))
    d = donors.reindex(x.index)
    return pd.DataFrame(dict(z=b[:, 0] + 1j * b[:, 1], zt=d.zt.to_numpy(float), arm=d[arm_col].astype(str).to_numpy()),
                        index=x.index), list(wf.genes)


def amplitude(s: pd.DataFrame) -> complex:
    t = s.zt.to_numpy() * W
    D = np.c_[np.ones(len(s)), np.exp(1j * t), np.exp(-1j * t)]  # z = mu + A e^{iWt} + B e^{-iWt}
    coef, *_ = np.linalg.lstsq(D, s.z.to_numpy(), rcond=None)
    return complex(coef[1])


def boot_amp(s: pd.DataFrame, rng) -> np.ndarray:
    """Resample libraries within each ZT (keeps the design)."""
    out = []
    groups = [g.index.to_numpy() for _, g in s.groupby("zt")]
    for _ in range(N_BOOT):
        idx = np.concatenate([rng.choice(g, len(g)) for g in groups])
        out.append(amplitude(s.loc[idx]))
    return np.array(out)


def ratio_row(label, a1, a0, b1, b0):
    r, rb = a1 / a0, b1 / b0
    return dict(contrast=label, amp_ratio=abs(r), amp_lo90=float(np.quantile(abs(rb), 0.05)),
                amp_hi90=float(np.quantile(abs(rb), 0.95)), shift_h=float(np.angle(r) / W),
                shift_lo90=float(np.quantile(np.angle(rb / r) / W, 0.05) + np.angle(r) / W),
                shift_hi90=float(np.quantile(np.angle(rb / r) / W, 0.95) + np.angle(r) / W)), r, rb


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(297)
    arms, rows = [], []
    amps, boots = {}, {}
    for name, tissue, col, keep in (("manella2021", "Liver", "arm", ("AlbCre_AL", "AlbCre_DF")),
                                    ("xin2021", "liver", "arm", ("NRF", "DRF")),
                                    ("greenwell2019", "liver", "arm", ("ad_libitum", "arrhythmic", "night_restricted"))):
        s, genes = liver_states(name, tissue, col)
        for a in keep:
            g = s[s.arm == a]
            amps[(name, a)] = amplitude(g)
            boots[(name, a)] = boot_amp(g, rng)
            arms.append(dict(dataset=name, arm=a, n=len(g), genes=len(genes), amplitude=abs(amps[(name, a)]),
                             acrophase_h=float(np.mod(np.angle(amps[(name, a)]), 2 * np.pi) / W)))
    pd.DataFrame(arms).to_csv(OUT / "arms.csv", index=False)

    k = lambda n, a: (amps[(n, a)], boots[(n, a)])  # noqa: E731
    r_man, rho_m, rho_mb = ratio_row("Manella DF / AL (rho, calibration)", *k("manella2021", "AlbCre_DF")[:1],
                                     k("manella2021", "AlbCre_AL")[0], k("manella2021", "AlbCre_DF")[1], k("manella2021", "AlbCre_AL")[1])
    r_xin, rho_x, rho_xb = ratio_row("Xin DRF / NRF (rho, second calibration)", k("xin2021", "DRF")[0], k("xin2021", "NRF")[0],
                                     k("xin2021", "DRF")[1], k("xin2021", "NRF")[1])
    rows += [r_man, r_xin]
    for lab, rho, rb in (("prediction from Manella: removal / AL = (1 + rho) / 2", rho_m, rho_mb),
                         ("prediction from Xin: removal / AL = (1 + rho) / 2", rho_x, rho_xb)):
        p, pb = (1 + rho) / 2, (1 + rb) / 2
        rows.append(dict(contrast=lab, amp_ratio=abs(p), amp_lo90=float(np.quantile(abs(pb), 0.05)),
                         amp_hi90=float(np.quantile(abs(pb), 0.95)), shift_h=float(np.angle(p) / W),
                         shift_lo90=float(np.quantile(np.angle(pb / p) / W, 0.05) + np.angle(p) / W),
                         shift_hi90=float(np.quantile(np.angle(pb / p) / W, 0.95) + np.angle(p) / W)))
    obs, _, _ = ratio_row("observed Greenwell arrhythmic / AL", k("greenwell2019", "arrhythmic")[0], k("greenwell2019", "ad_libitum")[0],
                          k("greenwell2019", "arrhythmic")[1], k("greenwell2019", "ad_libitum")[1])
    ctl, _, _ = ratio_row("control Greenwell night-restricted / AL", k("greenwell2019", "night_restricted")[0],
                          k("greenwell2019", "ad_libitum")[0], k("greenwell2019", "night_restricted")[1], k("greenwell2019", "ad_libitum")[1])
    rows += [obs, ctl]
    t = pd.DataFrame(rows)
    pred = t.iloc[2]
    inside = pred.amp_lo90 <= obs["amp_ratio"] <= pred.amp_hi90
    side = "inside" if inside else ("above" if obs["amp_ratio"] > pred.amp_hi90 else "below")
    t["H297_1"] = ""
    t.loc[t.contrast.str.startswith("observed"), "H297_1"] = f"{'PASS' if inside else 'FAIL'} ({side} the Manella prediction interval)"
    t.to_csv(OUT / "tests.csv", index=False)
    print(pd.DataFrame(arms).round(3).to_string(index=False))
    print(t.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
