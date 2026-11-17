"""A293 Part 2 (registered in `paper/docs/registered_experiments.md`): do stromal-vascular cell classes keep the same phase
order in every organ?

Human untimed pair offsets (A272, A274 pairs.csv), oriented as the first class named minus the second; readable pairs
|r| > matched q95; per atlas, circular mean weighted by |r|^2. Reference sign: the untimed offset of the same class
pair in mouse colon (A273, arm C-DC) and aorta (A279, aligned arm); the class pair is tested only if both agree.
H293.3: endothelial minus fibroblast has the mouse sign in more atlases than chance (one-sided sign test).

    python part2_phase_order.py
Outputs (`../output/part2/`): atlases.csv, tests.csv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
OUT = HERE.parents[1] / "output" / "part2"
W = 2 * np.pi / 24
CLASS_PAIRS = [("endothelial", "fibroblast_stromal", "H293.3"), ("endothelial", "mural", "secondary"),
               ("mural", "fibroblast_stromal", "secondary")]
MOUSE = {"colon": ("A273_pairwise_offset_truth", "C-DC", {"endothelial": "endothelial", "fibroblast_stromal": "fibroblast",
                                                          "mural": "smooth_muscle"}),
         "aorta": ("A279_aorta_pairwise_truth", "aligned", {"endothelial": "EC", "fibroblast_stromal": "Fibroblast",
                                                            "mural": "SMC"})}


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def mouse_offset(tissue, c1, c2):
    folder, arm, names = MOUSE[tissue]
    p = pd.read_csv(EXP / folder / "output" / "pairs.csv")
    p = p[p.arm == arm]
    a, b = names[c1], names[c2]
    r = p[(p.a == a) & (p.b == b)]
    if len(r):
        return float(r.untimed_h.iloc[0]), float(r.timed_h.iloc[0])
    r = p[(p.a == b) & (p.b == a)]
    return float(-r.untimed_h.iloc[0]), float(-r.timed_h.iloc[0])


def human_pairs():
    out = []
    for f in ("A272_backbone_pairwise", "A274_celltype_map"):
        p = pd.read_csv(EXP / f / "output" / "pairs.csv").assign(source=f)
        out.append(p[["source", "atlas", "a", "b", "family_a", "family_b", "r_abs", "offset_h", "matched_q95"]])
    return pd.concat(out, ignore_index=True)


def circ_mean(h, w):
    return float(np.angle(np.sum(w * np.exp(1j * W * np.asarray(h)))) / W)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    hp = human_pairs()
    rows, tests = [], []
    rng = np.random.default_rng(293)
    for c1, c2, label in CLASS_PAIRS:
        m = {t: mouse_offset(t, c1, c2) for t in MOUSE}
        signs = {t: np.sign(v[0]) for t, v in m.items()}
        sel = hp[((hp.family_a == c1) & (hp.family_b == c2)) | ((hp.family_a == c2) & (hp.family_b == c1))].copy()
        sel["oriented_h"] = np.where(sel.family_a == c1, sel.offset_h, -sel.offset_h)
        sel = sel[sel.r_abs > sel.matched_q95]
        per = []
        for atlas, g in sel.groupby("atlas"):
            mu = circ_mean(g.oriented_h, g.r_abs.to_numpy() ** 2)
            per.append(dict(class_pair=f"{c1}-{c2}", atlas=atlas, readable_pairs=len(g), offset_h=mu))
        per = pd.DataFrame(per)
        rows.append(per)
        t = dict(class_pair=f"{c1}-{c2}", test=label, mouse_colon_untimed_h=m["colon"][0], mouse_colon_timed_h=m["colon"][1],
                 mouse_aorta_untimed_h=m["aorta"][0], mouse_aorta_timed_h=m["aorta"][1], atlases=len(per))
        if signs["colon"] != signs["aorta"] or not len(per):
            t.update(tested=False, reason="mouse signs disagree" if signs["colon"] != signs["aorta"] else "no readable pair")
        else:
            s = signs["colon"]
            k = int((np.sign(per.offset_h) == s).sum())
            boot = [circ_mean(per.offset_h.to_numpy()[i], np.ones(len(per))) for i in
                    (rng.integers(0, len(per), len(per)) for _ in range(2000))]
            t.update(tested=True, mouse_sign=int(s), atlases_with_mouse_sign=k,
                     p_sign_one_sided=float(binomtest(k, len(per), 0.5, alternative="greater").pvalue),
                     human_circular_mean_h=circ_mean(per.offset_h, np.ones(len(per))),
                     boot_lo90=float(np.quantile(wrap(boot), 0.05)), boot_hi90=float(np.quantile(wrap(boot), 0.95)))
            t["passes"] = bool(label == "H293.3" and t["p_sign_one_sided"] < 0.05)
        tests.append(t)
    pd.concat(rows).to_csv(OUT / "atlases.csv", index=False)
    tests = pd.DataFrame(tests)
    tests.to_csv(OUT / "tests.csv", index=False)
    print(pd.concat(rows).round(3).to_string(index=False))
    print(tests.round(3).T.to_string())


if __name__ == "__main__":
    main()
