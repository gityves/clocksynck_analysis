# The body's clock as driven Hopf oscillators: model statement (draft, 2026-11-18)

Candidate take-home model of the paper (user, 2026-11-18). It is an organiser: a quantitative multi-tissue fit was
tried and failed (A310). Model numbers come from `exploratory/A311_hopf/output/hopf_grid.csv`, data numbers from the
results entries cited.

## Model

Each compartment k (a cell class in an organ) of a person, in the frame of that person's internal time theta(t):

    dz_k/dt = (mu_k + i nu_k) z_k - |z_k|^2 z_k + F_k + noise

- **mu_k**, distance to the Hopf bifurcation. mu > 0: a self-sustained cycle of amplitude about sqrt(mu). mu < 0: a damped
  responder that oscillates only while driven.
- **nu_k**, detuning from the drive (intrinsic period mismatch); with F_k it sets the lock lag.
- **F_k**, the drive from the person's internal time, common to all compartments up to its strength.
- **Organ:** the mixture of its compartments, sum_c (share of clock-gene mRNA) x z_c.

This is the generic form of any oscillator near the onset of oscillation (normal-form theory), the form used for
amplitude-phase clock models (Abraham 2010; Bordyugov 2015). Its natural units are mu ~ F^{2/3} and gain ~ F^{1/3}.

## The three measured coordinates as model quantities (exact for the locked state)

With x = |w*|^2 solving x((mu - x)^2 + nu^2) = F^2 and w* = -F / (mu + i nu - x):

| coordinate (clocksync) | model |
|---|---|
| gain (relative amplitude on the host) | sqrt(x): rises with mu; above the bifurcation it barely depends on F (saturation), below it it is proportional to F |
| lag (offset) | arg w*: set by nu / F, growing with mu; strong cycles with detuning can lose lock (no locked state at mu 6, nu 0.5) |
| departure 2 psi | stationary response of the cycle to person-level perturbations (Lyapunov equation at w*) |
| **orientation of the departure** | **share tangential minus radial: 0.015 at mu -4, 0.500 at mu 0, 0.904 at mu +4 (nu 0)**: a cycle above the bifurcation holds its radius and lets its phase move; a damped responder has no preferred direction |
| coupling S / (S + 2 psi) | follows: a fixed-size perturbation turns a weak cycle more; over mu -4 to +4 the model's slope of log departure on log gain^2 is 0.332 (observed 0.223, interval 0.044 to 0.332, A294b) |

## The paper's messages in the model

| message | in the model | status |
|---|---|---|
| One internal time per person, advancing with elapsed hours | one theta(t) drives every F_k | tested (A251-A256) |
| Time carried by stromal-vascular, epithelial, parenchymal cells and macrophages, not tissue lymphocytes | mu > 0 vs mu <= 0 (lymphocytes respond through hormone-induced genes, A294; a weak core loop is open) | tested for gain; the mu sign is an interpretation; mouse lymphocyte clocks are reported (Druzd 2017) |
| Cell classes share the organ's phase; composition moves gain, not phase | one drive, small detuning: one lag; gain set by mu_c and composition | tested (A273, A279, A295); class phase scatter 0.9 h, conserved class gain within a factor of about 1.4, not within errors (A310) |
| Organs differ 85-fold in gain over a nearly constant departure | gain from mu and composition; departure weakly rising with gain | tested descriptively (A291, A294b); one noise intensity does not fit every tissue (A310) |
| The organ arrangement needs a living brain, not the day or meals | F from the central drive; brain death F -> 0: organs free-run at their own nu and drift apart | consistent (A293 Part 1, descriptive; A132 tertiles in the drift direction, not significant) |
| Peripheral clocks are self-sustained | above the bifurcation, amplitude saturates: removing one drive keeps it | consistent (A297: 0.72 kept; the saturating model fits the removal amplitude, A310) |
| Fetal organs: a driven state before a running clock | mu < 0, F > 0: the drive reaches the organ before a loop runs; development crosses the bifurcation | consistent (A243-A269, A307, A309); the emergence of the cellular clock during differentiation is shown in vitro (Yagita 2010; Umemura 2017) |
| Proliferating tumours keep the patient's time with less gain | mu decreases, the lock is kept | consistent (A304) |
| Psoriasis delays the lesion clock | nu changes | consistent (A183, A308) |
| Chronic conditions change nothing | no parameter moves | tested (A164-A306) |

## Predictions not yet measured

1. **Orientation (GTEx; A311): NOT SUPPORTED.** Predicted: departures mostly tangential, more so in high-gain tissues.
   Observed: median share -0.03 (90 % -0.064 to -0.007), no significant gradient (rho 0.30, p 0.064), powered for a
   0.16 share. The GTEx departure is not the phase-dominated response of a cycle above its bifurcation; it reads as
   variation added to the state (A311, A302).
2. **Gain against drive (A312): sign supported, size small.** Lower-gain GTEx tissues follow the persons' shared
   amplitude more steeply (rho -0.34, p 0.037), but every tissue follows it nearly in proportion (elasticities 0.83
   to 1.77; low/high-gain ratio 1.11, against 3.4 for mu -1 vs +1). Adult organ gains are at most weakly a matter of
   distance to the bifurcation.
3. **Tumours in single cells.** A fall of mu means lower per-cell amplitude (damping), not cells drifting out of phase.
4. **Brain death.** Organ disorder grows with time since the drive was lost (not blind: A132 seen).

## Outside the model

The brain-body split (A263, A289), the GTEx adrenal offset (read as an effect of death), waveform and harmonics.
Coupling between cells (Webb 2012; Bernard 2007) can raise a compartment's effective mu; the model does not separate
cell-autonomous from coupling-generated cycles.
