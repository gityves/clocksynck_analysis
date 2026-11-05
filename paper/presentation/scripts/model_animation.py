"""Slide video: the host-referenced factor model, built step by step in the clock plane.

    beta_ik = mu_k + z_k f_i + e_ik,   z_k = |z_k| exp(i phi_k),   e_ik ~ N(0, sigma_k^2 I + V_ik)

One colour per patient, one marker per compartment, on a clock face (midnight on top, the day turning clockwise).
The host state is centred over patients (sum_i f_i = 0), so mu_k is the mean of compartment k's points: with sampling
in hospital hours it is not the centre of the cycle. Fibroblasts (the reference, z_F = 1) appear first, sampled in
hospital hours; then the terms are introduced one at a time: the compartment mean mu_k, the patient's true clock
(phase theta_i and amplitude A_i around the centre of the cycle, which is unknown, so these are never estimated), the
host state f_i the model does estimate (the patient's deviation from the cohort mean, f_i = A_i exp(i theta_i) - mean:
its angle is not the sampling time, as patient 4 shows), the compartment loading z_k (endothelial, hepatocyte: rotated
by phi_k, scaled by |z_k|), the compartment noise sigma_k^2, and the full equation. All values are illustrative.

Output (paper/presentation/figures/): model_animation.mp4 (1920 x 864, the slide below its title) and one still per
step (model_step<k>.png, the last frame of the step). Runs in about a minute.

    python scripts/model_animation.py [--fps 30] [--stills-only]
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch

import paperfig as pf

OUT = Path(__file__).resolve().parents[1] / "figures"
WIDTH, HEIGHT = pf.SLIDE_WIDTH, 152.0  # mm: the slide below the title band
PX_W = 1920

# ---- the illustrative model --------------------------------------------------------------------------------
HOURS = np.array([8.5, 10.0, 11.5, 13.0, 15.0, 17.0])  # sampling hour of each patient (unknown to the method)
AMP = np.array([1.0, 1.0, 1.0, 0.5, 1.0, 1.0])  # A_i: patient 4 has a weaker clock
THETA = 2 * np.pi * HOURS / 24  # clock phase of each patient (radians, 0 = midnight)
CLOCK = 1.35 * AMP * np.exp(1j * THETA)  # clock state (fibroblast units); its angle is the clock phase
F = CLOCK - CLOCK.mean()  # host state f_i, centred over patients: what the model estimates
PAT = pf.OKABE_ITO[:6]


def D(w):
    """Clock face: angle a (clock phase) drawn at pi/2 - a, i.e. midnight on top and the day turning clockwise."""
    return 1j * np.conj(w)


COMPS = {  # centre of the cycle (on screen), z_k, marker, noise sd
    "F": dict(name="fibroblast", c=2.3 + 1.85j, z=1.0 + 0j, marker="o", sd=0.08),
    "E": dict(name="endothelial", c=6.7 + 1.3j, z=0.75 * np.exp(1j * np.pi / 4), marker="^", sd=0.09),  # 3 h ahead
    "H": dict(name="hepatocyte", c=3.4 - 1.7j, z=1.15 * np.exp(-1j * np.pi / 3), marker="s", sd=0.2),  # 4 h behind
}
for _c in COMPS.values():  # mu_k: the mean of the compartment's points
    _c["mu"] = _c["c"] + D(_c["z"] * CLOCK.mean())
NOISE = {k: np.random.default_rng(i).normal(0, c["sd"], (len(F), 2)) @ np.array([1, 1j]) for i, (k, c) in enumerate(COMPS.items())}
INK, MUTED, HIGHLIGHT = "#222222", pf.GREY, "#fff1a8"
VIEW = np.array([-0.8, 8.6, -3.6, 3.45])  # xlim, ylim of the plane
ZOOM = np.array([0.05, 4.55, 0.12, 3.495])  # the fibroblasts (same aspect), for the patient clock in steps 3-4

STEPS = [  # caption, sub-caption, seconds of transition, seconds of hold
    ("Fibroblasts of six patients", "each sampled once, between 08 and 18 h", 3.0, 2.5),
    (r"$\mu_k$: the compartment mean", "constant offsets (baseline, platform) only move it;\nnot the centre of the cycle: sampling covers 08-18 h only", 1.5, 3.5),
    (r"$\theta_i$, $A_i$: each patient's clock", "phase (clockwise from the top) and amplitude\naround the centre of the cycle;"
     "\nthe centre is unknown, so they are never estimated", 2.0, 3.5),
    (r"$f_i$: what the model estimates", r"the patient's deviation from the mean $\mu_\mathrm{F}$;" + "\n"
     "its angle is not the sampling time:\nsee patient 4 (weaker clock), sampled at 13 h", 1.5, 4.5),
    ("Endothelial cells: the same patients", r"moved to $\mu_\mathrm{E}$, rotated by $\varphi_\mathrm{E}$, scaled by $|z_\mathrm{E}|$", 4.0, 3.0),
    ("Hepatocytes", r"rotated by $\varphi_\mathrm{H}$, scaled by $|z_\mathrm{H}|$" + "\n" + r"fibroblasts are the reference: $z_\mathrm{F} = 1$", 3.0, 3.0),
    (r"$\sigma_k^2$: each compartment's own noise", "hepatocytes follow the shared time less closely", 1.5, 3.0),
    ("The model", "fitted without any sampling time: absolute time cancels", 1.5, 5.0),
]


def ease(u):
    u = np.clip(u, 0.0, 1.0)
    return u * u * (3 - 2 * u)


def on(step, s, u, start=0.0, end=1.0):
    """0 before ``step``, easing in over [start, end] of its transition, 1 after."""
    return 0.0 if s < step else 1.0 if s > step else ease((u - start) / (end - start))


def visible(step_in, step_out, s, u):
    """Fades in at ``step_in``, out at ``step_out`` (None: stays)."""
    a = on(step_in, s, u)
    return a if step_out is None else a * (1 - on(step_out, s, u, 0.0, 0.4))


def loading(k, s, u):
    """Centre and complex loading of compartment k while it is being built (copy of the fibroblast pattern moved,
    rotated, then scaled), or its final values."""
    c, step = COMPS[k], {"F": 0, "E": 4, "H": 5}[k]
    if s > step or k == "F":
        return c["mu"], c["z"]
    if s < step:
        return COMPS["F"]["mu"], 1.0 + 0j
    mv, rot, sc = ease(u / 0.35), ease((u - 0.35) / 0.35), ease((u - 0.7) / 0.3)
    mu = COMPS["F"]["mu"] + mv * (c["mu"] - COMPS["F"]["mu"])
    return mu, (1 + sc * (abs(c["z"]) - 1)) * np.exp(1j * rot * np.angle(c["z"]))


def xy(z):
    z = np.atleast_1d(z)
    return np.column_stack([z.real, z.imag])


# ---- figure -------------------------------------------------------------------------------------------------
def build():
    lay = pf.Layout("PD;PE", width=WIDTH, height=HEIGHT, letters=False, col_ratios=[1.75, 1], row_heights=[1, 2.2],
                    polar="D", blank="PE", aspect={"P": 0.75}, anchor="C", margins=(3, 3, 3, 3))
    ax, dial, eqa = lay.axes["P"], lay.axes["D"], lay.axes["E"]
    ax.set(xlim=VIEW[:2], ylim=VIEW[2:])
    for a, b, lab, off in (((-0.6, 0), (8.0, 0), r"$\beta_1$", (0, -16)), ((0, -3.4), (0, 3.0), r"$\beta_2$", (-22, -10))):
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.0, mutation_scale=12))
        ax.annotate(lab, xy=b, xytext=off, textcoords="offset points", color=MUTED, ha="center")
    A = {"ax": ax}
    for k, c in COMPS.items():
        A[k, "circle"] = ax.add_patch(Circle((c["c"].real, c["c"].imag), abs(c["z"] * CLOCK[0]), fill=False, ec="#c8c8c8", lw=1.0, ls=(0, (4, 3))))
        A[k, "mu"] = ax.add_patch(FancyArrowPatch((0, 0), (c["mu"].real, c["mu"].imag), arrowstyle="-|>", mutation_scale=14,
                                                  color=INK, lw=1.4, shrinkA=0, shrinkB=4))
        A[k, "pts"] = ax.scatter(*xy(c["mu"] + D(F)).T, s=110, marker=c["marker"], c=PAT, edgecolors="white", linewidths=1.2, zorder=5)
        A[k, "noise"] = [ax.add_patch(Circle((p.real, p.imag), c["sd"], fill=False, ec=PAT[i], lw=0.8, ls=":", zorder=4))
                         for i, p in enumerate(c["mu"] + D(c["z"] * F))]
        A[k, "links"] = [ax.plot([], [], color=PAT[i], lw=0.8, ls=":", zorder=2)[0] for i in range(len(F))]
    mF, cF = COMPS["F"]["mu"], COMPS["F"]["c"]
    A["rays"] = [ax.plot(*xy([mF, mF + D(f)]).T, color=PAT[i], lw=1.2, zorder=3)[0] for i, f in enumerate(F)]  # f_i
    A["true"] = [ax.plot(*xy([cF, cF + D(w)]).T, color=PAT[i], lw=1.2, zorder=3)[0] for i, w in enumerate(CLOCK)]
    A["centre"] = ax.plot([cF.real], [cF.imag], marker="+", ms=9, mew=1.2, color=MUTED)[0]
    A["ref"] = ax.plot([cF.real, cF.real], [cF.imag, cF.imag + 0.9], color=MUTED, lw=0.8, ls="--")[0]  # 0 h: up
    a1 = 90 - np.degrees(THETA[0])  # screen angle of patient 1, clockwise from the top
    A["theta_arc"] = ax.add_patch(Arc((cF.real, cF.imag), 0.6, 0.6, theta1=a1, theta2=90, color=PAT[0], lw=1.4))
    mid = np.radians((a1 + 90) / 2)
    theta_anchor = cF + 0.3 * np.exp(1j * mid)
    for k in ("E", "H"):
        c = COMPS[k]
        A[k, "phi_arc"] = ax.add_patch(Arc((c["mu"].real, c["mu"].imag), 0.8, 0.8, color=INK, lw=1.2))
        A[k, "ghost"] = ax.plot([], [], color=PAT[0], lw=1.2, ls="--", zorder=3)[0]  # patient 1 before the rotation
        A[k, "ray"] = ax.plot([], [], color=PAT[0], lw=1.4, zorder=3)[0]  # and after
    # text labels of the plane, placed clear of every final mark, arrow and circle
    final = np.concatenate([c["mu"] + D(c["z"] * F) + d * NOISE[k] for k, c in COMPS.items() for d in (0, 1)])
    segs = [[(0, 0), (c["mu"].real, c["mu"].imag)] for c in COMPS.values()]
    for c in COMPS.values():
        ring = c["c"] + abs(c["z"] * CLOCK[0]) * np.exp(1j * np.linspace(0, 2 * np.pi, 33))
        segs += [[(a.real, a.imag), (b.real, b.imag)] for a, b in zip(ring[:-1], ring[1:])]
    for k, c in COMPS.items():
        segs += [[(c["mu"].real, c["mu"].imag), (p.real, p.imag)] for p in c["mu"] + D(c["z"] * F)]  # rays and arcs region
    segs += [[(cF.real, cF.imag), (p.real, p.imag)] for p in cF + D(CLOCK)]  # true clock rays
    segs.append([(cF.real, cF.imag), (cF.real, cF.imag + 0.9)])  # 0 h reference of theta
    kw = dict(obstacles=xy(final), segments=np.array(segs), max_dist=14, leader_from=4, clearance=2.2)  # markers: ~1.9 mm radius
    mus = [0.5 * c["mu"] for c in COMPS.values()]
    A["mu_labs"] = pf.place_labels(ax, *xy(mus).T, [rf"$\mu_\mathrm{{{k}}}$" for k in COMPS], prefer=[2.4, 1.7, -2.4],
                                   fontsize=15, color=INK, **kw)
    tops = [c["c"] + 1j * abs(c["z"] * CLOCK[0]) * (1 if k != "H" else -1) for k, c in COMPS.items()]
    A["names"] = pf.place_labels(ax, *xy(tops).T, [c["name"] for c in COMPS.values()], prefer=[np.pi / 2, np.pi / 2, -np.pi / 2],
                                 fontsize=13, color=MUTED, **kw)
    A["theta_lab"] = pf.place_labels(ax, *xy(theta_anchor).T, [r"$\theta_1$"], prefer=[mid], fontsize=14, color=PAT[0], **kw)
    phi_anchor = []
    for k in ("E", "H"):  # middle of the final phase arc
        c = COMPS[k]
        a0, ang = np.angle(D(F[0])), -np.angle(c["z"])
        phi_anchor.append(c["mu"] + 0.4 * np.exp(1j * (a0 + ang / 2)))
    A["phi_labs"] = pf.place_labels(ax, *xy(phi_anchor).T, [r"$\varphi_\mathrm{E}$", r"$\varphi_\mathrm{H}$"],
                                    prefer=[np.angle(a - COMPS[k]["mu"]) for a, k in zip(phi_anchor, "EH")], fontsize=14,
                                    color=INK, **kw)
    A["amp_lab"] = pf.place_labels(ax, *xy(cF + D(CLOCK[3])).T, [r"$A_4$"], prefer=[np.pi], fontsize=14, color=PAT[3], **kw)
    A["f_lab"] = pf.place_labels(ax, *xy(mF + D(F[3])).T, [r"$f_4$"], prefer=[np.angle(D(F[3]))], fontsize=14,
                                 color=PAT[3], **kw)
    A["z_labs"] = {k: pf.place_labels(ax, *xy(COMPS[k]["c"] + abs(COMPS[k]["z"] * CLOCK[0])).T,
                                      [rf"$|z_\mathrm{{{k}}}| = {abs(COMPS[k]["z"]):.2f}$"], prefer=[0.0], fontsize=13, color=INK, **kw)
                   for k in ("E", "H")}

    # sampling dial
    dial.set_theta_zero_location("N"); dial.set_theta_direction(-1)
    dial.set_thetagrids([0, 90, 180, 270], ["0 h", "6", "12", "18"]); dial.set_yticks([]); dial.set_ylim(0, 1)
    dial.bar(np.radians(15 * 8), 1, width=np.radians(150), align="edge", color="#eeeeee", zorder=0)
    A["ticks"] = [dial.plot([t, t], [0, 0.95], color=PAT[i], lw=2.4)[0] for i, t in enumerate(THETA)]
    dial.set_title("sampling hour (unknown to the method)", fontsize=11, color=MUTED, pad=4)
    A["dial"] = dial

    # captions and equations, one set per step (faded in and out)
    A["cap"] = [(eqa.text(0.0, 0.97, t, fontsize=16, color=INK, va="top", transform=eqa.transAxes),
                 eqa.text(0.0, 0.84, st, fontsize=12, color=MUTED, va="top", transform=eqa.transAxes, wrap=True))
                for t, st, _, _ in STEPS]
    eqs = [
        ([r"\beta_{i\mathrm{F}}", "=", "?"], None),
        ([r"\beta_{i\mathrm{F}}", "=", r"\mu_\mathrm{F}", "+", "?"], 2),
        ([r"\beta_{i\mathrm{F}}", "=", r"c_\mathrm{F}", "+", r"A_i\, e^{i\theta_i}"], 4),
        ([r"\beta_{i\mathrm{F}}", "=", r"\mu_\mathrm{F}", "+", r"f_i"], 4),
        ([r"\beta_{i\mathrm{E}}", "=", r"\mu_\mathrm{E}", "+", r"|z_\mathrm{E}|\, e^{i\varphi_\mathrm{E}}", r"f_i"], 4),
        ([r"\beta_{ik}", "=", r"\mu_k", "+", r"z_k", r"f_i"], 4),
        ([r"\beta_{ik}", "=", r"\mu_k", "+", r"z_k f_i", "+", r"e_{ik}"], 6),
    ]
    A["eq"] = []
    for terms, hl in eqs:
        e = pf.equation(eqa, (0.0, 0.5), terms, ha="left", size=24)
        box = None
        if hl is not None:
            x0, y0, x1, y1 = e.boxes[hl]
            box = eqa.add_patch(FancyBboxPatch((x0 - 2, y0 - 2), x1 - x0 + 4, y1 - y0 + 4, transform=e.transform,
                                               boxstyle="round,pad=0,rounding_size=3", fc=HIGHLIGHT, ec="none", zorder=0))
        A["eq"].append((e, box))
    e3 = eqa.text(0.0, 0.26, r"$c_\mathrm{F}$: centre of the cycle (unknown)", fontsize=13, color=MUTED, transform=eqa.transAxes)
    A["eq"][2] = (*A["eq"][2], e3, None)
    e4 = pf.equation(eqa, (0.0, 0.26), [r"f_i = A_i\, e^{i\theta_i} - \overline{A\, e^{i\theta}}"], ha="left", size=18)
    t4 = eqa.text(0.0, 0.1, "the centre cancels; the time does not come back", fontsize=12, color=MUTED, transform=eqa.transAxes)
    A["eq"][3] = (*A["eq"][3], e4, None, t4, None)
    # step 7 adds the noise distribution, step 8 the labelled model
    e7 = pf.equation(eqa, (0.0, 0.26), [r"e_{ik} \sim \mathcal{N}(0,\ ", r"\sigma_k^2", r"I + V_{ik})"], ha="left", size=18)
    x0, y0, x1, y1 = e7.boxes[1]
    b7 = eqa.add_patch(FancyBboxPatch((x0 - 2, y0 - 2), x1 - x0 + 4, y1 - y0 + 4, transform=e7.transform,
                                      boxstyle="round,pad=0,rounding_size=3", fc=HIGHLIGHT, ec="none", zorder=0))
    A["eq"][-1] = (A["eq"][-1][0], A["eq"][-1][1], e7, b7)
    e8 = pf.equation(eqa, (0.5, 0.5), [r"\beta_{ik}", "=", (r"\mu_k", "centre", INK), "+",
                                       dict(tex=r"z_k", label="compartment\n" + r"$|z_k|,\ \varphi_k$", color=INK, side="above"),
                                       (r"f_i", "host state\n(centred)", INK), "+",
                                       (r"e_{ik}", "noise\n" + r"$\sigma_k^2$", INK)], size=24, gap=12, label_size=14, label_sep=6, marks="none")
    t8 = eqa.text(0.5, 0.08, r"reference: $z_\mathrm{F} = 1$" + "\nreadouts per compartment: coupling, $|z_k|$, $\\varphi_k$",
                  ha="center", va="bottom", fontsize=12, color=MUTED, transform=eqa.transAxes)
    A["eq8"] = (e8, t8)
    fig = lay.finalize(check=False)  # every step's text is laid over the others here; each step is checked when rendered
    return fig, A


def eq_artists(e):
    if not hasattr(e, "texts"):  # a plain Text
        return [e]
    out = list(e.texts)
    for _, lab, arts in e.labels:
        out += [lab, *arts]
    return out


def update(A, s, u):
    """Set every artist for step s (0-based) at progress u in [0, 1] of its transition."""
    pts = {}
    for k, c in COMPS.items():
        mu, z = loading(k, s, u)
        p = mu + D(z * F)
        noisy = on(6, s, u)
        p_obs = p + noisy * NOISE[k]
        pts[k] = p_obs
        A[k, "pts"].set_offsets(xy(p_obs))
        if k == "F":  # appear one patient at a time in step 1
            alpha = np.array([np.clip((u * 7 - i), 0, 1) if s == 0 else 1.0 for i in range(len(F))])
        else:
            alpha = np.full(len(F), on({"E": 4, "H": 5}[k], s, u, 0.0, 0.15))
        A[k, "pts"].set_alpha(None)
        A[k, "pts"].set_facecolors([(*_rgb(PAT[i]), a) for i, a in enumerate(alpha)])
        A[k, "pts"].set_edgecolors([(1, 1, 1, a) for a in alpha])
        step = {"F": 1, "E": 4, "H": 5}[k]
        A[k, "mu"].set_alpha(on(step, s, u) if k == "F" else on(step, s, u, 0.0, 0.35))
        j = list(COMPS).index(k)
        _fade(A["mu_labs"], j, on(step, s, u) if k == "F" else on(step, s, u, 0.3, 0.4))
        _fade(A["names"], j, on(0 if k == "F" else {"E": 4, "H": 5}[k], s, u, 0.0, 0.3))
        A[k, "circle"].set_alpha(on(step, s, u) if k == "F" else on(step, s, u, 0.7, 1.0))
        A[k, "pts"].set_visible(alpha.max() > 0.01)
        for n, q in zip(A[k, "noise"], p_obs):  # centred on the observed points
            n.set_center((q.real, q.imag))
            n.set_alpha(0.9 * on(6, s, u))
        if k != "F":  # dotted links patient by patient while the compartment is built
            la = 0.6 * visible(step, step + 1, s, u)
            for i, ln in enumerate(A[k, "links"]):
                ln.set_data(*xy([pts["F"][i], p_obs[i]]).T)
                ln.set_alpha(la)
            # phase arc from patient 1's fibroblast direction to its rotated direction, then its scale label
            a0 = np.degrees(np.angle(D(F[0])))  # patient 1's fibroblast direction on screen
            ang = -np.degrees(np.angle(z))  # a phase advance turns clockwise on the clock face
            arc = A[k, "phi_arc"]
            arc.set_center((mu.real, mu.imag))
            arc.theta1, arc.theta2 = (a0, a0 + ang) if ang >= 0 else (a0 + ang, a0)
            arc.set_alpha(visible(step, 6, s, u) * (on(step, s, u, 0.35, 0.5)))
            _fade(A["phi_labs"], "EH".index(k), visible(step, 6, s, u) * on(step, s, u, 0.7, 0.8))
            A[k, "ghost"].set_data(*xy([mu, mu + D(F[0])]).T)
            A[k, "ray"].set_data(*xy([mu, mu + D(z * F[0])]).T)
            A[k, "ghost"].set_alpha(arc.get_alpha())
            A[k, "ray"].set_alpha(arc.get_alpha())
            _fade(A["z_labs"][k], 0, visible(step, 6, s, u) * on(step, s, u, 0.85, 1.0))
    for i, (r, t) in enumerate(zip(A["rays"], A["true"])):  # true clock in step 3 (faint in 4), f_i in step 4
        t.set_alpha(visible(2, 4, s, u) * (1 - 0.75 * on(3, s, u)))
        r.set_alpha(visible(3, 4, s, u) * (1.0 if i == 3 else 0.6))
    A["true"][3].set_alpha(visible(2, 4, s, u) * (1 - 0.4 * on(3, s, u)))  # patient 4's true direction stays visible
    A["true"][3].set_linestyle("-" if s < 3 else "--")
    A["rays"][3].set_linewidth(2.2)
    for key in ("ref", "theta_arc", "centre"):
        A[key].set_alpha(visible(2, 3, s, u))
    _fade(A["theta_lab"], 0, visible(2, 3, s, u))
    _fade(A["amp_lab"], 0, visible(2, 3, s, u))
    _fade(A["f_lab"], 0, visible(3, 4, s, u))
    for i, t in enumerate(A["ticks"]):
        t.set_alpha(np.clip(u * 7 - i, 0, 1) if s == 0 else 1.0)
    for j, (t, st) in enumerate(A["cap"]):
        a = 1.0 if j == s else 0.0
        if j == s and s > 0:
            a = ease(u / 0.3)
        t.set_alpha(a)
        st.set_alpha(a)
    for j, item in enumerate(A["eq"]):
        a = (ease(u / 0.3) if s > 0 else 1.0) if j == s else 0.0
        for e in item[0::2]:
            for art in eq_artists(e):
                art.set_alpha(a)
        for b in item[1::2]:
            if b is not None:
                b.set_alpha(a)
    e8, t8 = A["eq8"]
    a8 = on(7, s, u, 0.0, 0.4)
    for art in eq_artists(e8) + [t8]:
        art.set_alpha(a8)
    for ax in A["dial"].figure.axes:  # fully faded artists are hidden (and skipped by the collision check)
        for art in ax.texts + ax.patches + ax.lines:
            a = art.get_alpha()
            art.set_visible(a is None or a > 0.01)
    lim = VIEW + ease(visible(2, 4, s, u)) * (ZOOM - VIEW)  # zoom on the fibroblasts for the patient clock
    if not np.array_equal(lim, A.get("lim")):
        A["lim"] = lim
        A["ax"].set(xlim=lim[:2], ylim=lim[2:])
        pf.place_all(A["ax"].figure)  # labels follow the zoom (placed among the texts visible now)


def _fade(labelset, j, a):
    """Alpha of label j of a paperfig LabelSet, with its leader line."""
    t = labelset.texts[j]
    t.set_alpha(a)
    t.arrow_patch.set_alpha(a)


def _rgb(hexcol):
    h = hexcol.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--stills-only", action="store_true")
    args = ap.parse_args()
    pf.use_style("slide")
    fig, A = build()
    fig.set_dpi(PX_W / (WIDTH / 25.4))
    pf.place_all(fig)  # at the video resolution
    w, h = (int(round(v)) for v in fig.canvas.get_width_height())
    w, h = w - w % 2, h - h % 2
    OUT.mkdir(exist_ok=True)
    enc = None
    if not args.stills_only:
        enc = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgba",
                                "-s", f"{w}x{h}", "-r", str(args.fps), "-i", "-", "-c:v", "libx264", "-pix_fmt",
                                "yuv420p", "-crf", "18", str(OUT / "model_animation.mp4")], stdin=subprocess.PIPE)
    for s, (_, _, t_in, t_hold) in enumerate(STEPS):
        n_in, n_hold = int(t_in * args.fps), int(t_hold * args.fps)
        frames = [(i / n_in) for i in range(n_in)] + [1.0] * n_hold if enc else [1.0]
        for u in frames:
            update(A, s, u)
            if enc:
                fig.canvas.draw()
                buf = np.asarray(fig.canvas.buffer_rgba())[:h, :w]
                enc.stdin.write(np.ascontiguousarray(buf).tobytes())
        update(A, s, 1.0)
        fig.savefig(OUT / f"model_step{s + 1}.png", dpi=fig.dpi)
        issues = pf.collisions(fig)
        print(f"step {s + 1}: {len(issues)} collision(s)" + "".join(f"\n  {m}" for _, m in issues), flush=True)
    if enc:
        enc.stdin.close()
        enc.wait()
        print("wrote", OUT / "model_animation.mp4")


if __name__ == "__main__":
    main()
