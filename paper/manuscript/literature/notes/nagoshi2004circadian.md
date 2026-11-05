# nagoshi2004circadian
**Ref:** Nagoshi, Saini, Bauer, Laroche, Naef, Schibler, Cell 2004. DOI 10.1016/j.cell.2004.11.015. Type: biology (with modelling). Lab: Schibler (Naef co-author).
**Source:** full text

## Question and design
Do fibroblasts carry cell-autonomous, self-sustained oscillators, and does damping of the population signal reflect
single-cell damping or dephasing? NIH3T3 reporter lines: Rev-erba-Venus-NLS-PEST (time-lapse, 30 min frames, 3 days;
42 and 26 nuclei in Figs 1D, 2B) and Bmal1-Luc (luminometer, 19 days). Co-culture with 20-fold excess of non-luminescent
feeder cells; dexamethasone pulses for phase response curves; proliferating versus confluent cells. **Sampling:** timed
by the experimenter (shock at t = 0), continuous in vitro recording, no human or tissue sampling.

## Method (if any)
Three models of the population envelope (pure damping, pure dephasing of sustained oscillators, both), fitted to the
detrended luminescence peaks/troughs for t < 10 days, finite populations simulated (5000 oscillators).

## Key results
- Single cells oscillate robustly: period Gaussian, mean 28.2 h, SD 2.9 h after serum shock (Fig. 1E); about 25 h at serum >= 5%.
- Without a shock, cells are out of phase: "individual cells displayed robust but out-of-phase circadian VNP accumulation." (Fig. 2)
  Serum therefore synchronises existing oscillators instead of starting them.
- Population luminescence amplitude falls with half-life 3.2 days (cell loss); envelope fits favour dephasing
  (sigma of frequency 0.93 h) over damping; small oscillations persist in finite populations (Fig. 3).
- A 15 min dexamethasone pulse resets every phase: phase shifts up to +/-12 h, slope of the phase transition curve near zero (type-zero reset) (Fig. 4).
- Co-culture with Rat-1 (6 h advanced), Per1-/- MEFs (about 4 h shorter period) or Per1/Per2 double knockout fibroblasts does not change
  phase or period of the reporter: "NIH3T3 cells do not communicate with each other with regard to circadian rhythm generation."
- Clock runs through division, daughters resume the mother's rhythm; doubling times 60, 34, 28 h at 5, 10, 20% serum;
  cytokinesis is gated into three windows (about 6, 12, 19 h after the VNP peak) and division shifts circadian phase (Figs 6, 7).

## Not addressed / limits
- One immortalised line, reporter-based, serum or dexamethasone synchronised; fluorescence imaging limited to ~3 days
  by photodamage; the "no coupling" conclusion applies to dispersed confluent fibroblasts, not to tissue.
- Does not test response to systemic cues from an intact animal (the authors say so).
- The dephasing model is a fit to population envelopes; single-cell tracks are short.

## Link to us
- **Intro claim(s):** single-cell clocks are cell-autonomous and self-sustained; population amplitude reflects
  synchrony (Internal desynchrony, within tissue). Our refs summary (including the clock continuing through the cell cycle)
  is supported; the "independently phased" statement is supported via Fig. 2.
- **Relation to clocksync:** background. Shows that a bulk (pseudobulk) clock readout confounds single-cell amplitude
  with cell-to-cell dispersion; our |z| is such a mixture.

## Ideas
- Within-compartment dispersion bounds the shared-clock fraction in a pseudobulk: simulate from their sigma_f (0.83-0.93 h
  on a 24.5 h period) how much dephasing alone lowers |z|.
- Cell-cycle gating (Fig. 6E) predicts a cell-cycle contribution to clock-gene pseudobulks in proliferating compartments (tumours): control with CoPhaser cell-cycle phase.
- Missing in our list if relevant: Pando 2002 (MEFs respond to systemic cues), Bieler 2014 (cell cycle - clock coupling in single cells).
