# reddy2002differential
**Ref:** Reddy, Field, Maywood, Hastings, J Neurosci 2002. DOI 10.1523/JNEUROSCI.22-17-07326.2002. Type: biology. Lab: Hastings.
**Source:** full text

## Question and design
Where does inertia in re-entrainment to a shifted light cycle arise? Male CD1 mice (6 wk), 12L:12DR (dim red), in situ hybridisation of
mPer1, mPer2 and mCry1 mRNA in the whole SCN, plus mPER/mCRY immunostaining. Exp 1: 6 h advance of lights-on; hourly kills for 18 h on
days 1, 3 and 8 after the shift (3 mice/hour) and on a baseline day. Exp 2: 6 h delay; hourly kills for 44 h. Wheel running in 5 extra mice.
**Sampling:** timed, hourly, SCN only (no peripheral tissue).

## Key results
- Baseline phases (ZT): "mPer1 mRNA peaked around ZT06, and mPer2 and mCry1 mRNA around lights-off (ZT12)" (Results).
- Behaviour after a 6 h advance: activity onset advanced by 4.17 +/- 0.25 h on day 3. "Stable re-entrainment was only achieved after 1 week." (Results)
- Advance: mPer1 and mPer2 are induced by the new light and return to the control cycle within 2-3 h; mCry1 is not acutely affected. "Dissociation of mPer and mCry1 expression was again apparent on the third day of the advance." (Results); mCry1 peak at day 3 is intermediate between baseline and day 8.
- "By day 8 the rhythm of mCry1 expression had re-entrained to the advanced lighting cycle" and re-synchronised with mPer2 (three-way ANOVA gene x time x day interaction, p stated as highly significant).
- Delay of 6 h: behaviour completed within 2 cycles, mPer and mCry1 shift together; "The phase delay was complete by the second night" (Results) for mPer1, mPer2 by the second day and night, mCry1 matched the projected 6 h delay on the second cycle.
- Cell counts: advance, mPER2-IR 618.7 +/- 74.5 vs 160.3 +/- 25.8 cells per section in controls (n = 3); delay, mCRY1-IR 389 +/- 62 vs 98 +/- 12 (p < 0.01).
- Model: inertia sits within the SCN in slow mCry1 response; yet the authors state that peripheral inertia is also real: "In transgenic rats rPer1::Luc expression in peripheral tissues readjusts more slowly, indicating that inertia lies in the periphery not in the" SCN (Discussion, citing Yamazaki 2000).

## Not addressed / limits
- Not about SCN subregions (shell vs core): signals are measured in the whole SCN, regional remarks are only about where induction is strongest (ventral SCN).
- Single time series of pooled mice (3 per hour), no peripheral tissues, no individual-animal covariance, no inference on tissue offsets.
- Gene-level desynchrony (Per vs Cry) inside the same oscillator, so a tissue's "phase" is gene dependent during transients.

## Link to us
- **Intro claim(s):** CLAIM WRONG as currently worded. refs_circadian_desynchrony.md says SCN subregions (shell vs core) re-entrain at different rates (also line 106, SCN subregions). The paper shows dissociation between mPer and mCry1 expression within the SCN, not between subregions. Use de la Iglesia 2004 for subregions and Reddy for gene-level dissociation.
- **Relation to clocksync:** background. Shows that during transients different clock genes of one tissue sit at different phases, relevant to how a gene-phase model (common phase per compartment) can fail.

## Ideas
- Stress-test the gene-phase coherence gate D with gene-specific partial shifts (Per genes shifted, Cry/Bmal1 not) within a compartment; scrambled-phase simulations exist (METHOD.md, D 4.70 vs 0.42) but a Per-vs-Cry specific dissociation is not documented.
- Reference we lack: Davidson et al. 2009 (Curr Biol, visualising jet lag in the mouse SCN with PER2::LUC; shell-core resetting differences).
