"""Slide video (method, slide M1): how one pseudobulk's 14 clock genes become one cosine, hence one point beta.

Each clock gene g follows y_g = m_g + A_g cos(phi - psi_g) (log expression; baseline m_g, amplitude A_g and acrophase
psi_g from the GTEx waveform, clocksync.gtex()). Steps:
 1. the 14 genes over the day; a sample is taken at one internal phase phi (cursor on the right, hand on the dial);
 2. phi turns: each gene is high when phi is near its acrophase (dial: point size = cos(phi - psi_g));
 3. remove the baseline m_g; 4. divide by the amplitude A_g: every gene reads cos(phi - psi_g);
 5. place each gene at its acrophase: the 14 values trace one cosine of psi, peaking at phi;
 6. phi turns: the cosine follows;
 7. a real sample is noisy (Poisson counts plus gene noise): the weighted least-squares cosine gives phi_hat and
    lam_hat, i.e. beta;
 8. weights: after division by A_g, low-amplitude and low-count genes are noisiest; the fit weighs gene g by
    w_g A_g^2 (w_g = 1 / (1/c_g + tau^2), the package's gene weight), against the unweighted fit;
 9. in data the baseline is the compartment's mean over donors; with donors sampled 10-18 h that mean still
    oscillates: cbar cos(14 h - psi), cbar = sin(D)/D = 0.83 for a uniform window of half-width D = 4 h;
10. subtracting it spreads the apparent peaks: donors at 10 / 14 / 18 h read 6.6 / 14.0 / 21.4 h, with relative
    heights 0.93 / 0.17 / 0.93 (dial: hand = the donor's vector minus the mean vector, as in the clock plane).
Fitting (y_g - m_g) / A_g with weights w_g A_g^2 is the package's projection (WLS of y_g on A_g (cos psi_g, sin psi_g)
with weights w_g): the fits use clocksync's projector (`_core.plane.projector`), applied to y - m. The sample is
simulated (depth DEPTH, gene noise TAU_SIM, seed SEED); illustrative, not fitted to data.

Output (paper/presentation/figures/projection_animation/, 1920 x 864, the slide below its title): the full video
projection_animation.mp4, and for slides: proj_01.mp4 (transition of step 1), proj_02.png (its final frame), ...

    python scripts/projection_animation.py [--fps 30] [--stills-only]
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image

import clocksync as cs
import paperfig as pf
from clocksync._core.plane import projector

OUT = Path(__file__).resolve().parents[1] / "figures" / "projection_animation"
WIDTH, HEIGHT = pf.SLIDE_WIDTH, 152.0  # mm: the slide below the title band
PX_W = 1920
H24 = 24 / (2 * np.pi)
PHI0 = 15.0  # internal phase of the sample (h)
DEPTH, TAU_SIM, TAU2, SEED = 4e5, 0.05, 0.05, 47  # seed: weighting visibly helps (phase)
WINDOW = (10.0, 18.0)  # hospital hours of steps 9-10
DONORS = {10.0: "#CC79A7", 14.0: "#D55E00", 18.0: "#009E73"}  # donor phase (h): colour
HALF = (WINDOW[1] - WINDOW[0]) / 2 / H24
CBAR = np.sin(HALF) / HALF  # mean of cos over a uniform window, relative to its centre
SHOW = {"ARNTL": pf.OKABE_ITO[0], "NR1D1": pf.OKABE_ITO[1], "PER2": pf.OKABE_ITO[2]}
INK, MUTED, FAINT, FIT = "#222222", pf.GREY, "#b4b4b4", pf.OKABE_ITO[4]
CLIP = 2.7  # display limit of the scaled values (noisy low-amplitude genes go beyond)

WF = cs.gtex()
GENES = list(WF.genes)
PSI = np.mod(WF.psi, 2 * np.pi)  # radians
A = np.asarray(WF.amp)
M = WF.mesor - WF.mesor.mean()
COL = [SHOW.get(g, FAINT) for g in GENES]

STEPS = [  # caption, sub-caption, equation, seconds of transition, seconds of hold
    ("Each clock gene oscillates around its own baseline", "log expression over the day: small rhythms on very different baselines",
     r"$y_g(t) = m_g + A_g \cos(t - \psi_g)$", 3.0, 2.5),
    ("A sample is taken at one internal phase φ", "each gene is high when φ is close to its peak time ψ$_g$",
     r"$y_g = m_g + A_g \cos(\varphi - \psi_g)$", 6.0, 2.0),
    ("Remove the baseline", "in data: each gene's mean over the compartment's donors",
     r"$y_g - m_g = A_g \cos(\varphi - \psi_g)$", 2.0, 2.5),
    ("Divide by the amplitude", "every gene now swings between −1 and +1",
     r"$(y_g - m_g)\,/\,A_g = \cos(\varphi - \psi_g)$", 2.0, 2.5),
    ("Place each gene at its peak time ψ$_g$", "the 14 values trace one cosine, whose peak is the internal phase φ",
     r"$\cos(\varphi - \psi)$ as a function of $\psi$", 3.0, 3.0),
    ("The cosine follows the internal phase", "two numbers describe the sample: peak φ and height λ",
     r"$\beta = \lambda\,(\cos\varphi,\ \sin\varphi)$", 6.0, 2.0),
    ("A real sample is noisy: fit the best cosine", r"weighted least squares over the 14 genes gives $\hat\varphi$ and $\hat\lambda$",
     r"$\hat\beta = \arg\min_\beta \sum_g w_g A_g^2 \left(\tilde y_g - \lambda\cos(\varphi - \psi_g)\right)^2$", 3.0, 3.0),
    ("Noisy genes count less", "after dividing by $A_g$, low-amplitude and low-count genes are noisiest (point area = weight)",
     r"weight $= w_g A_g^2$,   $w_g = 1 / (1/\bar c_g + \tau^2)$", 3.0, 4.0),
    ("In data, the baseline is the mean over donors", "donors sampled 10-18 h: their mean still oscillates (dashed)",
     r"mean $= m_g + 0.83\,A_g \cos(14\,\mathrm{h} - \psi_g)$", 3.0, 3.0),
    ("Subtracting it spreads the peaks", "8 h of sampling become apparent peaks from 6.6 to 21.4 h; the middle donor flattens",
     r"$\cos(\varphi_i - \psi) - 0.83\cos(14\,\mathrm{h} - \psi)$", 4.0, 4.0),
]
YLIM = {"raw": (-2.1, 2.4), "centred": (-0.7, 0.7), "unit": (-1.45, 1.45), "noisy": (-3.0, 3.0)}


def ease(u):
    u = np.clip(u, 0.0, 1.0)
    return u * u * (3 - 2 * u)


def on(step, s, u, start=0.0, end=1.0):
    """0 before ``step``, easing in over [start, end] of its transition, 1 after."""
    return 0.0 if s < step else 1.0 if s > step else ease((u - start) / (end - start))


def sample():
    """One simulated pseudobulk at PHI0: scaled values (y - m) / A, the gene weights (w A^2) and both fits."""
    rng = np.random.default_rng(SEED)
    log_rate = WF.mesor + A * np.cos(PHI0 / H24 - PSI) + rng.normal(0, TAU_SIM, len(A))
    counts = rng.poisson(DEPTH * np.exp(log_rate))
    y = np.log(np.maximum(counts, 0.5) / DEPTH) - WF.mesor
    cbar = DEPTH * np.exp(WF.mesor)  # expected count at the baseline (the package uses the compartment's median)
    w = 1.0 / (1.0 / cbar + TAU2)
    fits = {}
    for key, wt in (("weighted", w), ("unweighted", np.ones_like(w))):
        b = projector(WF, wt) @ y  # (cos, sin) coefficients
        fits[key] = (np.hypot(*b), np.mod(np.arctan2(b[1], b[0]), 2 * np.pi) * H24)
    return y / A, w * A ** 2, fits


def phase_at(s, u):
    if s == 0:
        return 6.0 + (PHI0 - 6.0) * ease(u)
    if s in (1, 5):
        return PHI0 + 24.0 * ease(u)
    return PHI0


def ylim_at(s, u):
    lerp = lambda a, b, t: tuple(np.add(np.multiply(1 - t, a), np.multiply(t, b)))
    if s <= 1:
        return YLIM["raw"]
    if s == 2:
        return lerp(YLIM["raw"], YLIM["centred"], ease(u))
    if s == 3:
        return lerp(YLIM["centred"], YLIM["unit"], ease(u))
    if s == 6:
        return lerp(YLIM["unit"], YLIM["noisy"], ease(u / 0.4))
    if s == 8:
        return lerp(YLIM["noisy"], YLIM["unit"], ease(u / 0.4))
    return YLIM["noisy"] if s == 7 else YLIM["unit"]


def build(data):
    scaled, weight, fits = data
    lay = pf.Layout("CC;DP", width=WIDTH, height=HEIGHT, letters=False, col_ratios=[1, 2.3], row_heights=[0.42, 1],
                    polar="D", blank="C", anchor="C", margins=(3, 3, 3, 3), wgap=14)
    cap, dial, ax = lay.axes["C"], lay.axes["D"], lay.axes["P"]
    S = {}
    # dial: genes at their acrophase, the internal phase as a hand
    dial.set_theta_zero_location("N"); dial.set_theta_direction(-1)
    dial.set_thetagrids([0, 90, 180, 270], ["0 h", "6", "12", "18"], color=MUTED)
    dial.set_yticks([]); dial.set_ylim(0, 1.25); dial.grid(False)
    dial.spines["polar"].set_color(MUTED)
    S["dial_pts"] = dial.scatter(PSI, np.ones_like(PSI), s=40, c=COL, edgecolor="white", linewidth=0.8, zorder=3)
    S["hand"] = dial.plot([0, 0], [0, 0.85], color=INK, lw=2.4, solid_capstyle="round", zorder=4)[0]
    S["hand_w"] = dial.plot([0, 0], [0, 0.85], color=FIT, lw=2.0, ls=(0, (4, 2)), zorder=4)[0]
    S["window"] = dial.bar(np.mean(WINDOW) / H24, 1.25, width=2 * HALF, color="#eeeeee", zorder=0)[0]
    S["donor_hands"] = {d: dial.plot([0, 0], [0, 0.85], color=c, lw=2.4, solid_capstyle="round", zorder=4)[0]
                        for d, c in DONORS.items()}
    S["donor_labs"] = {d: dial.text(0, 0, f"{d:.0f} h", ha="center", va="center", fontsize=12, color=c)
                       for d, c in DONORS.items()}
    S["dial_labs"] = [dial.text(PSI[GENES.index(g)], 1.17, g, ha="center", va="center", fontsize=11, color=c)
                      for g, c in SHOW.items()]
    S["phi_lab"] = dial.text(0, 0.55, "φ", fontsize=15, color=INK, ha="center", va="center")
    # main axes
    t = np.linspace(0, 24, 241)
    S["curves"] = [ax.plot(t, np.zeros_like(t), color=c, lw=2.0 if c != FAINT else 1.0, zorder=2)[0] for c in COL]
    S["cursor"] = ax.axvline(PHI0, color=INK, lw=1.2, ls=(0, (4, 3)), zorder=1)
    S["psi_curve"] = ax.plot(t, np.zeros_like(t), color=INK, lw=2.2, zorder=2)[0]
    S["fit_w"] = ax.plot(t, np.zeros_like(t), color=FIT, lw=2.6, zorder=3)[0]
    S["fit_u"] = ax.plot(t, np.zeros_like(t), color=MUTED, lw=1.8, ls=(0, (5, 3)), zorder=3)[0]
    S["donor_curves"] = {d: ax.plot(t, np.zeros_like(t), color=c, lw=2.2, zorder=3)[0] for d, c in DONORS.items()}
    S["mean_curve"] = ax.plot(t, np.zeros_like(t), color=INK, lw=2.0, ls=(0, (5, 3)), zorder=3)[0]
    S["dots"] = ax.scatter(np.zeros(len(A)), np.zeros(len(A)), s=60, c=COL, edgecolor="white", linewidth=0.9, zorder=4)
    S["labels"] = {g: ax.text(0, 0, g, fontsize=11, color=c, ha="left", va="center", zorder=5) for g, c in SHOW.items()}
    S["unw_lab"] = ax.text(0, 0, "unweighted", fontsize=11, color=MUTED, ha="left", va="top", zorder=5)
    ax.axhline(0, color="#e2e2e2", lw=0.8, zorder=0)
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])
    S["xlab_t"] = ax.text(0.5, -0.13, "time of day (h)", transform=ax.transAxes, ha="center", va="top", fontsize=13)
    S["xlab_p"] = ax.text(0.5, -0.13, "peak time of each gene, ψ (h)", transform=ax.transAxes, ha="center", va="top",
                          fontsize=13)
    S["ylab"] = [ax.text(-0.075, 0.5, lab, transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=13)
                 for lab in ("log expression (centred)", "log expression − baseline", "(log expression − baseline) / amplitude",
                             "(log expression − donor mean) / amplitude")]
    # captions
    S["cap"] = [(cap.text(0.0, 0.92, c, fontsize=19, color=INK, va="top", transform=cap.transAxes),
                 cap.text(0.0, 0.52, sc, fontsize=13, color=MUTED, va="top", transform=cap.transAxes),
                 cap.text(1.0, 0.52, eq, fontsize=17, color=INK, va="top", ha="right", transform=cap.transAxes))
                for c, sc, eq, _, _ in STEPS]
    fig = lay.finalize(check=False)  # every step's text is laid over the others here
    return fig, S, ax, dial


def update(S, ax, data, s, u):
    scaled, weight, fits = data
    phi = phase_at(s, u)
    p = phi / H24
    bm = 1 - on(2, s, u)  # baseline kept
    da = on(3, s, u)  # amplitude divided out
    sp = on(4, s, u)  # spread to the acrophases
    nz = on(6, s, u, 0.0, 0.4)  # noise
    ws = on(7, s, u, 0.0, 0.6)  # point area = weight
    amp = A ** (1 - da)
    clean = bm * M + amp * np.cos(p - PSI)
    y = (1 - nz) * clean + nz * np.clip(scaled, -CLIP, CLIP)
    x = (1 - sp) * phi % 24 + sp * PSI * H24 if sp < 1 else PSI * H24
    if sp == 0:
        x = np.full(len(A), phi % 24)
    ax.set_ylim(*ylim_at(s, u))
    gone = 1 - on(8, s, u, 0.0, 0.4)  # steps 9-10 keep only the donors' cosines
    # day curves: fade out while the genes spread
    t = np.linspace(0, 24, 241)
    ca = 1 - on(4, s, u, 0.0, 0.4)
    for i, ln in enumerate(S["curves"]):
        ln.set_ydata(bm * M[i] + amp[i] * np.cos(t / H24 - PSI[i]))
        ln.set_alpha(ca * (1.0 if COL[i] != FAINT else 0.8))
    S["cursor"].set_xdata([phi % 24] * 2)
    S["cursor"].set_alpha(ca)
    S["psi_curve"].set_ydata(np.cos(p - t / H24))
    S["psi_curve"].set_alpha(on(4, s, u, 0.6, 1.0) * (1 - 0.65 * nz) * gone)
    # dots: area = weight in step 8
    size = (1 - ws) * 70 + ws * np.maximum(260 * weight / weight.max(), 10)
    S["dots"].set_offsets(np.column_stack([x, y]))
    S["dots"].set_sizes(size)
    S["dots"].set_alpha(gone)
    for g, txt in S["labels"].items():
        i = GENES.index(g)
        txt.set_position((x[i] + 0.35, y[i] + (0.06 if sp < 1 else 0.0)))
        txt.set_alpha(on(0, s, u, 0.6, 1.0) * gone)
    # fits
    lw, pw = fits["weighted"]
    lu, pu = fits["unweighted"]
    S["fit_w"].set_ydata(lw * np.cos(pw / H24 - t / H24))
    S["fit_w"].set_alpha(on(6, s, u, 0.6, 1.0) * gone)
    S["fit_u"].set_ydata(lu * np.cos(pu / H24 - t / H24))
    S["fit_u"].set_alpha(on(7, s, u, 0.4, 1.0) * gone)
    S["unw_lab"].set_position(((pu + 1.0) % 24, lu * 1.0 + 0.45))
    S["unw_lab"].set_alpha(on(7, s, u, 0.6, 1.0) * gone)
    # steps 9-10: donors in hospital hours, minus their mean
    da9, k = on(8, s, u, 0.4, 1.0), on(9, s, u)
    centre = np.mean(WINDOW) / H24
    S["mean_curve"].set_ydata((1 - k) * CBAR * np.cos(centre - t / H24))
    S["mean_curve"].set_alpha(da9 * (1 - 0.7 * k))
    S["window"].set_alpha(da9)
    for d, ln in S["donor_curves"].items():
        ln.set_ydata(np.cos(d / H24 - t / H24) - k * CBAR * np.cos(centre - t / H24))
        ln.set_alpha(da9)
        v = np.exp(1j * d / H24) - k * CBAR * np.exp(1j * centre)
        S["donor_hands"][d].set_data([np.angle(v)] * 2, [0, 0.85 * abs(v)])
        S["donor_hands"][d].set_alpha(da9)
        S["donor_labs"][d].set_position((np.angle(v), 0.85 * abs(v) + 0.17))
        S["donor_labs"][d].set_alpha(da9)
    # dial
    S["hand"].set_xdata([p, p])
    S["hand"].set_alpha(gone)
    S["phi_lab"].set_alpha(gone)
    S["phi_lab"].set_position((p + 0.25, 0.55))
    resp = on(1, s, u, 0.0, 0.1)
    S["dial_pts"].set_sizes(40 + resp * 160 * (1 + np.cos(p - PSI)) / 2)
    S["hand_w"].set_xdata([pw / H24] * 2)
    S["hand_w"].set_alpha(on(6, s, u, 0.6, 1.0) * gone)
    S["dial_pts"].set_alpha(gone)
    for txt in S["dial_labs"]:
        txt.set_alpha(gone)
    # axis labels and captions
    S["xlab_t"].set_alpha(1 - sp)
    S["xlab_p"].set_alpha(sp)
    ya = [1 - on(2, s, u), on(2, s, u) * (1 - on(3, s, u)), on(3, s, u) * (1 - on(9, s, u)), on(9, s, u)]
    for txt, a in zip(S["ylab"], ya):
        txt.set_alpha(a)
    for j, texts in enumerate(S["cap"]):
        a = (ease(u / 0.25) if s > 0 else 1.0) if j == s else 0.0
        for txt in texts:
            txt.set_alpha(a)
    for art in ax.texts + ax.lines + S["cap"][0][0].axes.texts + S["dial_pts"].axes.lines + S["dial_pts"].axes.texts:
        a = art.get_alpha()
        art.set_visible(a is None or a > 0.01)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--stills-only", action="store_true")
    args = ap.parse_args()
    pf.use_style("slide")
    data = sample()
    fig, S, ax, dial = build(data)
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

    full = None if args.stills_only else encoder(OUT / "projection_animation.mp4")
    for s, (*_, t_in, t_hold) in enumerate(STEPS):
        n_in, n_hold = int(t_in * args.fps), int(t_hold * args.fps)
        clip = None if args.stills_only else encoder(OUT / f"proj_{2 * s + 1:02d}.mp4")
        for i in range(0 if args.stills_only else n_in + 1):
            update(S, ax, data, s, i / n_in)
            buf = frame()
            clip.stdin.write(buf)
            if i < n_in:
                full.stdin.write(buf)
        if clip:
            clip.stdin.close()
            clip.wait()
        update(S, ax, data, s, 1.0)
        buf = frame()
        for _ in range(0 if args.stills_only else n_hold):
            full.stdin.write(buf)
        Image.frombytes("RGBA", (w, h), buf).save(OUT / f"proj_{2 * s + 2:02d}.png")
    if full:
        full.stdin.close()
        full.wait()
    lw, pw = data[2]["weighted"]
    lu, pu = data[2]["unweighted"]
    print(f"weighted fit: phi {pw:.1f} h, lam {lw:.2f}; unweighted: phi {pu:.1f} h, lam {lu:.2f}; true {PHI0} h, 1")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
