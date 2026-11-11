"""Slide video (method, slide M4a): step 1 of the fit, the host by EM on the anchor compartments.

    beta_ik = mu_k + z_k f_i + e_ik,   f_i ~ N(0, Sigma_f),   e_ik ~ N(0, sigma_k^2 I + V_ik),   z_F = 1

Left, the data (clock plane): three anchor compartments of six patients (one colour per patient, one marker per
compartment, as in the model video). Right, the host plane: one point f_i per patient, shared by its compartments.
Steps, text on screen only:
 1. the data and what is known (beta_ik, V_ik) or unknown (mu_k, z_k, sigma_k^2, f_i, Sigma_f);
 2. a starting guess: mu_k the compartment mean, z_k = 1, sigma_k^2 and Sigma_f wide;
 3. E-step (1/2): each compartment gives a guess of the patient's host state, z_k^-1 (beta_ik - mu_k);
 4. E-step (2/2): the guesses are combined, weighted by their precision and pulled toward 0 by Sigma_f: posterior
    mean f_i_hat and covariance S_i;
 5. M-step (1/2): mu_k and z_k by least squares of beta_ik on the host: the host pattern is shifted, turned and scaled
    onto each compartment;
 6. M-step (2/2): sigma_k^2 (what is left beyond the Poisson V_ik) and Sigma_f (spread of the host; drawn at 1.5 SD, S_i at 2 SD);
 7. repeat until nothing changes: the guesses of each patient's compartments come to agree;
 8. the result: f_i_hat and S_i, the regressor of step 2.

The E and M steps are the package's own (`clocksync._core.factor._EM.estep / mstep`; SQUAREM acceleration not
shown), from the starting guess of step 2 (not the package's initialisation). V_ik is isotropic here, so the E-step is
exactly the precision-weighted mean written on screen; the package uses the 2 x 2 matrix form. Patients, loadings and
noise are illustrative (patients and fibroblast/endothelial loadings as in model_animation.py; smooth muscle 0.9 at
+1 h).

Output (paper/presentation/figures/em_animation/, 1920 px wide): em_animation.mp4 and, for slides, em_01.mp4 (the
transition of step 1), em_02.png (its final frame), ... one clip and one still per step.

    python scripts/em_animation.py [--fps 30] [--stills-only]
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch
from PIL import Image

import paperfig as pf
from clocksync._core import factor

OUT = Path(__file__).resolve().parents[1] / "figures" / "em_animation"
WIDTH, HEIGHT = pf.SLIDE_WIDTH, 152.0  # mm: the slide below the title band
PX_W = 1920
H24 = 24 / (2 * np.pi)
N_ITER = 20

# ---- the illustrative data -----------------------------------------------------------------------------------
HOURS = np.array([8.5, 10.0, 11.5, 13.0, 15.0, 17.0])  # sampling hour of each patient (unknown to the method)
AMP = np.array([1.0, 1.0, 1.0, 0.85, 1.0, 1.0])
CLOCK = 1.35 * AMP * np.exp(1j * 2 * np.pi * HOURS / 24)
F_TRUE = CLOCK - CLOCK.mean()
PAT = pf.OKABE_ITO[:6]


def D(w):
    """Clock face: angle a drawn at pi/2 - a (midnight on top, the day turning clockwise)."""
    return 1j * np.conj(w)


COMPS = {  # key: name, marker, loading, noise sd, centre on screen
    "F": ("fibroblast", "o", 1.0 + 0j, 0.08, -3.3 + 0.35j),
    "E": ("endothelial", "^", 0.75 * np.exp(1j * np.pi / 4), 0.09, 0.0 + 0.35j),
    "S": ("smooth muscle", "s", 0.9 * np.exp(1j * np.pi / 12), 0.12, 3.3 + 0.35j),
}
KEYS = list(COMPS)
V_ISO = 0.004  # Poisson variance per coordinate (known)
INK, MUTED, FAINT, HIGHLIGHT = "#222222", pf.GREY, "#c8c8c8", "#fff1a8"


def simulate():
    rng = np.random.default_rng(3)
    beta, s, k = [], [], []
    for j, (key, (_, _, z, sd, centre)) in enumerate(COMPS.items()):
        mu = np.conj(centre / 1j)  # model-space centre whose display is `centre`
        beta += list(mu + z * F_TRUE + rng.normal(0, sd, 6) + 1j * rng.normal(0, sd, 6))
        s += list(range(6))
        k += [j] * 6
    return np.array(beta), np.array(s), np.array(k)


def trace(beta, s, k):
    """Parameters after each M-step and the E-step moments at each parameter set, from the step-2 guess."""
    b = np.c_[beta.real, beta.imag]
    V = np.repeat(V_ISO * np.eye(2)[None], len(b), 0)
    em = factor._EM(b, V, s, k, k, 6, 3, 3, 0, 1e-8)
    p = (np.ones(3, complex), np.array([b[k == j].mean(0) for j in range(3)]), np.full(3, 0.3), 0.5 * np.eye(2))
    params, moments = [], []
    for _ in range(N_ITER + 2):
        z, mu, psi, sig = p
        _, W, m, C = em.estep(z, mu, psi, sig)
        params.append(dict(z=z, mu=mu[:, 0] + 1j * mu[:, 1], psi=psi, sig=sig))
        moments.append(dict(m=m[:, 0] + 1j * m[:, 1], C=C))
        p = em.mstep(z, mu, psi, W, m, C)
    return params, moments


# ---- helpers -------------------------------------------------------------------------------------------------
def ease(u):
    u = np.clip(u, 0.0, 1.0)
    return u * u * (3 - 2 * u)


def on(step, s, u, start=0.0, end=1.0):
    """0 before `step`, easing in over [start, end] of its transition, 1 after."""
    return 0.0 if s < step else 1.0 if s > step else float(ease((u - start) / (end - start)))


def lerp(a, b, t):
    return a + (b - a) * t


def blend(states, x):
    """State at fractional index x: every array of the dict interpolated between its neighbours."""
    i = int(np.clip(np.floor(x), 0, len(states) - 2))
    t = float(np.clip(x - i, 0, 1))
    return {key: lerp(states[i][key], states[i + 1][key], t) for key in states[i]}


def xy(w):
    w = np.atleast_1d(w)
    return np.column_stack([w.real, w.imag])


def ellipse_params(cov, n_sd=2.0):
    """Width, height and angle (degrees) of the n_sd ellipse of a model-space 2 x 2 covariance, on the clock face."""
    c = np.array([[cov[1, 1], cov[0, 1]], [cov[0, 1], cov[0, 0]]])  # D swaps the two coordinates
    w, v = np.linalg.eigh(c)
    return 2 * n_sd * np.sqrt(max(w[1], 0)), 2 * n_sd * np.sqrt(max(w[0], 0)), np.degrees(np.arctan2(v[1, 1], v[0, 1]))


STEPS = [  # caption, sub-caption, highlighted terms of the model, formula of the step, seconds in, seconds hold
    ("Three anchor compartments of six patients", "one point per pseudobulk; its Poisson noise V$_{ik}$ is known",
     [0], r"known: $\beta_{ik},\ V_{ik}$        unknown: $\mu_k,\ z_k,\ \sigma_k^2$ (per compartment),   "
          r"$f_i$ (per patient),   $\Sigma_f$", 2.5, 3.0),
    ("Start from a guess", r"$\mu_k$: the compartment mean;  $z_k = 1$ (the arrow: 0 h, unturned);  "
                           r"$\sigma_k^2$ and $\Sigma_f$ wide",
     [2, 4], r"$\mu_k = \overline{\beta}_{\cdot k}, \quad z_k = 1$", 2.0, 3.0),
    ("E-step: each compartment guesses the patient's host state",
     "undo the compartment: subtract its centre, turn and scale back",
     [5], r"$\hat g_{ik} = z_k^{-1}\,(\beta_{ik} - \mu_k)$", 2.5, 3.0),
    ("E-step: combine the guesses", r"weighted by their precision, pulled toward 0 by the host spread $\Sigma_f$",
     [5], r"$\hat f_i = S_i \sum_k w_{ik}\,\hat g_{ik}, \quad S_i = (\Sigma_f^{-1} + \sum_k w_{ik})^{-1}, "
          r"\quad w_{ik} = |z_k|^2 / (\sigma_k^2 + v_{ik})$", 2.5, 3.5),
    ("M-step: place the host pattern on each compartment",
     r"$\mu_k$ shifts it, $z_k$ turns and scales it, to best fit the points",
     [2, 4], r"$(\mu_k, z_k) = \arg\min \sum_i \mathbb{E}\,\|\beta_{ik} - \mu_k - z_k f_i\|^2_{W_{ik}}$",
     3.0, 3.0),
    ("M-step: the unshared variance and the host spread", r"$\sigma_k^2$: what is left around the fit, beyond $V_{ik}$;  "
                                              r"$\Sigma_f$: the spread of the host states",
     [7], r"$\mathbb{E}\,\|e_{ik}\|^2 = 2\sigma_k^2 + \mathrm{tr}\,V_{ik}, \quad "
          r"\Sigma_f = \frac{1}{n}\sum_i (\hat f_i \hat f_i^{\top} + S_i)$", 2.0, 3.0),
    ("Repeat E and M until nothing changes", "the guesses of each patient's compartments come to agree",
     [], r"the likelihood increases at every iteration", 7.0, 2.0),
    ("Result: the host state of every patient", r"$\hat f_i$ with its uncertainty $S_i$: "
                                                "the regressor of step 2, where every compartment is read",
     [5], r"$\hat f_i,\ S_i$", 2.0, 4.0),
]


# ---- figure --------------------------------------------------------------------------------------------------
def build(beta, s, k, params, moments):
    lay = pf.Layout("CC;PH;EE", width=WIDTH, height=HEIGHT, letters=False, col_ratios=[2.35, 1],
                    row_heights=[0.3, 1, 0.3], blank="CE", anchor="C", margins=(3, 3, 3, 3), wgap=10)
    cap, ax, hx, eqa = lay.axes["C"], lay.axes["P"], lay.axes["H"], lay.axes["E"]
    A = {"ax": ax, "hx": hx}
    for a, lim, title in ((ax, ((-5.0, 5.0), (-1.3, 1.55)), r"data: the clock plane, $\beta_{ik}$"),
                          (hx, ((-1.55, 1.55), (-1.55, 1.55)), r"host plane, $f_i$")):
        a.set(xlim=lim[0], ylim=lim[1], aspect="equal", xticks=[], yticks=[])
        for sp in a.spines.values():
            sp.set_visible(False)
        a.set_title(title, fontsize=13, color=MUTED, pad=4)
    hx.axhline(0, color="#eeeeee", lw=0.8, zorder=0)
    hx.axvline(0, color="#eeeeee", lw=0.8, zorder=0)
    hx.add_patch(Circle((0, 0), 1.52, fill=False, ec="#e6e6e6", lw=0.8, zorder=0))

    # data plane
    disp = D(beta)
    for j, key in enumerate(KEYS):
        name, marker, *_ = COMPS[key]
        i = k == j
        A[key, "pts"] = ax.scatter(*xy(disp[i]).T, s=95, marker=marker, c=[PAT[q] for q in s[i]], edgecolors="white",
                                   linewidths=1.0, zorder=5)
        A[key, "pred"] = ax.scatter(*xy(disp[i]).T, s=95, marker=marker, facecolors="none",
                                    edgecolors=[PAT[q] for q in s[i]], linewidths=1.4, zorder=4)
        A[key, "res"] = [ax.plot([], [], color=PAT[q], lw=0.9, zorder=3)[0] for q in s[i]]
        A[key, "sig"] = [ax.add_patch(Circle((0, 0), 0.1, fill=False, ec=PAT[q], lw=0.8, ls=":", zorder=3))
                         for q in s[i]]
        A[key, "mu"] = ax.scatter([0], [0], s=160, marker="+", color=INK, linewidths=2.0, zorder=6)
        A[key, "z"] = ax.add_patch(FancyArrowPatch((0, 0), (0, 1), arrowstyle="-|>", mutation_scale=14, color=INK,
                                                   lw=1.6, shrinkA=0, shrinkB=0, zorder=6))
        A[key, "mu_lab"] = ax.text(0, 0, rf"$\mu_\mathrm{{{key}}}$", fontsize=14, color=INK, ha="right", va="bottom",
                                   zorder=7)
        A[key, "z_lab"] = ax.text(0, 0, rf"$z_\mathrm{{{key}}}$", fontsize=14, color=INK, ha="left", va="bottom",
                                  zorder=7)
        A[key, "name"] = ax.text(COMPS[key][4].real, -1.25, name, fontsize=13, color=MUTED, ha="center", va="bottom")

    # host plane
    A["guess"] = {key: hx.scatter(np.zeros(6), np.zeros(6), s=55, marker=COMPS[key][1], c=PAT, edgecolors="white",
                                  linewidths=0.6, zorder=4) for key in KEYS}
    A["fhat"] = hx.scatter(np.zeros(6), np.zeros(6), s=110, marker="o", c=PAT, edgecolors=INK, linewidths=1.2, zorder=5)
    A["S"] = [hx.add_patch(Ellipse((0, 0), 0.1, 0.1, fill=False, ec=PAT[q], lw=1.0, zorder=3)) for q in range(6)]
    A["Sigma"] = hx.add_patch(Ellipse((0, 0), 1, 1, fill=False, ec=MUTED, lw=1.2, ls=(0, (5, 3)), zorder=2))
    A["Sigma_lab"] = hx.text(0, 0, r"$\Sigma_f$", fontsize=14, color=MUTED, ha="left", va="bottom")
    A["f_lab"] = hx.text(0, 0, r"$\hat f_1,\ S_1$", fontsize=14, color=PAT[0], ha="left", va="bottom", zorder=6)
    A["iter"] = hx.text(0.0, 0.0, "", transform=hx.transAxes, fontsize=13, color=MUTED, ha="left", va="bottom")

    # captions
    A["cap"] = [(cap.text(0.0, 0.95, c, fontsize=19, color=INK, va="top", transform=cap.transAxes),
                 cap.text(0.0, 0.42, sc, fontsize=13, color=MUTED, va="top", transform=cap.transAxes))
                for c, sc, *_ in STEPS]
    # model with highlighted terms, and the formula of the step
    terms = [r"\beta_{ik}", "=", r"\mu_k", "+", r"z_k", r"f_i", "+", r"e_{ik}"]
    A["eq"] = pf.equation(eqa, (0.0, 0.72), terms, ha="left", size=22)
    A["boxes"] = []
    for x0, y0, x1, y1 in A["eq"].boxes:
        A["boxes"].append(eqa.add_patch(FancyBboxPatch((x0 - 2, y0 - 2), x1 - x0 + 4, y1 - y0 + 4,
                                                       transform=A["eq"].transform, boxstyle="round,pad=0,rounding_size=3",
                                                       fc=HIGHLIGHT, ec="none", zorder=0)))
    A["form"] = [eqa.text(0.0, 0.08, f, fontsize=16, color=INK, va="bottom", transform=eqa.transAxes)
                 for *_, f, _, _ in STEPS]
    fig = lay.finalize(check=False)  # every step's text is laid over the others here; each step is checked when rendered
    return fig, A


def state(s, u, params, moments):
    """Displayed parameters and host moments for step s at progress u."""
    P0, P1 = params[0], params[1]
    if s <= 3:
        return P0, moments[0]
    if s == 4:  # mu and z move to the M-step values; sigma and Sigma_f follow in step 5
        t = ease(u)
        return dict(z=lerp(P0["z"], P1["z"], t), mu=lerp(P0["mu"], P1["mu"], t), psi=P0["psi"], sig=P0["sig"]), moments[0]
    if s == 5:
        t = ease(u)
        return dict(z=P1["z"], mu=P1["mu"], psi=lerp(P0["psi"], P1["psi"], t), sig=lerp(P0["sig"], P1["sig"], t)), moments[0]
    x = N_ITER * (ease(u) if s == 6 else 1.0)  # iterations: parameters run one M-step ahead of the moments
    return blend(params[1:], x), blend(moments, x)


def update(A, s, u, beta, ss, kk, params, moments):
    P, Mo = state(s, u, params, moments)
    ax, hx = A["ax"], A["hx"]
    f, C = Mo["m"], Mo["C"]
    a_init = on(1, s, u, 0.0, 0.5)
    a_pred = on(4, s, u, 0.0, 0.3) * (1 - on(6, s, u, 0.0, 0.3))
    a_sig = on(5, s, u, 0.0, 0.5)
    a_guess = (on(2, s, u, 0.0, 0.5) * (1 - on(3, s, u, 0.5, 1.0)) + on(6, s, u, 0.0, 0.15)) * (1 - on(7, s, u, 0.0, 0.5))
    a_f = on(3, s, u, 0.4, 0.9)
    for j, key in enumerate(KEYS):
        i = kk == j
        A[key, "pts"].set_alpha(1.0 if s > 0 else float(ease(u / 0.6)))
        mu, z = P["mu"][j], P["z"][j]
        dm = D(mu)
        A[key, "mu"].set_offsets([[dm.real, dm.imag]])
        A[key, "mu"].set_alpha(a_init)
        tip = D(mu + 0.75 * z)
        A[key, "z"].set_positions((dm.real, dm.imag), (tip.real, tip.imag))
        A[key, "z"].set_alpha(a_init)
        A[key, "mu_lab"].set_position((dm.real - 0.1, dm.imag + 0.06))
        A[key, "mu_lab"].set_alpha(a_init * (key == "F"))  # one compartment labelled: the marks mean the same on all
        A[key, "z_lab"].set_position((tip.real + 0.06, tip.imag + 0.02))
        A[key, "z_lab"].set_alpha(a_init * (key == "F"))
        pred = D(mu + z * f[ss[i]])
        A[key, "pred"].set_offsets(xy(pred))
        A[key, "pred"].set_alpha(a_pred)
        obs = D(beta[i])
        for ln, a, b in zip(A[key, "res"], pred, obs):
            ln.set_data([a.real, b.real], [a.imag, b.imag])
            ln.set_alpha(0.7 * a_pred)
        r = 2 * np.sqrt(P["psi"][j] + V_ISO)
        for c, p in zip(A[key, "sig"], pred):
            c.set_center((p.real, p.imag))
            c.set_radius(r)
            c.set_alpha(a_sig)
        # guesses: z_k^-1 (beta_ik - mu_k); in step 3 they slide into f_hat
        g = D((beta[i] - mu) / z)
        if s == 3:
            g = lerp(g, D(f[ss[i]]), ease((u - 0.1) / 0.6))
        A["guess"][key].set_offsets(xy(g))
        A["guess"][key].set_alpha(a_guess)
    df = D(f)
    A["fhat"].set_offsets(xy(df))
    A["fhat"].set_alpha(a_f)
    for e, p, c in zip(A["S"], df, C):
        w, h, ang = ellipse_params(c)
        e.set_center((p.real, p.imag))
        e.width, e.height, e.angle = w, h, ang
        e.set_alpha(a_f)
    w, h, ang = ellipse_params(P["sig"], n_sd=1.5)
    A["Sigma"].width, A["Sigma"].height, A["Sigma"].angle = w, h, ang
    A["Sigma"].set_alpha(on(3, s, u, 0.4, 0.9))
    A["Sigma_lab"].set_position((0.5 * w * 0.72, -0.5 * h - 0.02))
    A["Sigma_lab"].set_alpha(A["Sigma"].get_alpha())
    A["f_lab"].set_position((df[0].real + 0.12, df[0].imag + 0.05))
    A["f_lab"].set_alpha(a_f * (1 if s in (3, 7) else 0.0))
    if s >= 6:
        n = int(round(N_ITER * (ease(u) if s == 6 else 1.0)))
        A["iter"].set_text(f"iteration {max(n, 1)}")
    A["iter"].set_alpha(on(6, s, u, 0.0, 0.1))
    # captions, highlighted terms, formula
    for j, texts in enumerate(A["cap"]):
        a = (float(ease(u / 0.25)) if s > 0 else 1.0) if j == s else 0.0
        for t in texts:
            t.set_alpha(a)
    hl = STEPS[s][2]
    for j, b in enumerate(A["boxes"]):
        b.set_alpha(1.0 if j in hl else 0.0)
    for j, t in enumerate(A["form"]):
        t.set_alpha((float(ease(u / 0.25)) if s > 0 else 1.0) if j == s else 0.0)
    for a in A["ax"].figure.axes:  # fully faded artists are hidden (and skipped by the collision check)
        for art in a.texts + a.patches + a.lines + a.collections:
            al = art.get_alpha()
            art.set_visible(al is None or al > 0.01)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--stills-only", action="store_true")
    args = ap.parse_args()
    pf.use_style("slide")
    beta, s, k = simulate()
    params, moments = trace(beta, s, k)
    fig, A = build(beta, s, k, params, moments)
    fig.set_dpi(PX_W / (WIDTH / 25.4))
    w, h = (int(round(v)) for v in fig.canvas.get_width_height())
    w, h = w - w % 2, h - h % 2
    OUT.mkdir(parents=True, exist_ok=True)

    def encoder(path):
        return subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgba",
                                 "-s", f"{w}x{h}", "-r", str(args.fps), "-i", "-", "-c:v", "libx264", "-pix_fmt",
                                 "yuv420p", "-crf", "18", str(path)], stdin=subprocess.PIPE)

    def frame():
        fig.canvas.draw()
        return np.ascontiguousarray(np.asarray(fig.canvas.buffer_rgba())[:h, :w]).tobytes()

    full = None if args.stills_only else encoder(OUT / "em_animation.mp4")
    for st, (*_, t_in, t_hold) in enumerate(STEPS):
        n_in, n_hold = int(t_in * args.fps), int(t_hold * args.fps)
        clip = None if args.stills_only else encoder(OUT / f"em_{2 * st + 1:02d}.mp4")
        for i in range(0 if args.stills_only else n_in + 1):
            update(A, st, i / n_in, beta, s, k, params, moments)
            buf = frame()
            clip.stdin.write(buf)
            if i < n_in:
                full.stdin.write(buf)
        if clip:
            clip.stdin.close()
            clip.wait()
        update(A, st, 1.0, beta, s, k, params, moments)
        buf = frame()
        for _ in range(0 if args.stills_only else n_hold):
            full.stdin.write(buf)
        Image.frombytes("RGBA", (w, h), buf).save(OUT / f"em_{2 * st + 2:02d}.png")
        issues = pf.collisions(fig)
        print(f"step {st + 1}: {len(issues)} collision(s)" + "".join(f"\n  {m}" for _, m in issues), flush=True)
    if full:
        full.stdin.close()
        full.wait()
    fin = params[N_ITER + 1]
    print("final: |z|", np.round(np.abs(fin["z"]), 3), "offset h", np.round(np.angle(fin["z"]) * H24, 2),
          "psi", np.round(fin["psi"], 4))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
