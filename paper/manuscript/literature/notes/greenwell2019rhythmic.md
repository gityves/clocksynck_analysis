# greenwell2019rhythmic
**Ref:** Greenwell, Trott, Beytebiere, Pao, Bosley, Beach, Finegan, Hernandez, Menet, Cell Rep 2019 (27:649-657). DOI 10.1016/j.celrep.2019.03.064. Type: biology. Lab: Menet (Texas A&M).
**Source:** full text

## Question and design
Mouse liver under three feeding paradigms in LD12:12 for 5 weeks: night-restricted (NR), ad libitum (LB) and arrhythmic (AR, 1/8th of daily intake every 3 h with a programmable 8-compartment dispenser). Does the rhythm of food intake (RFI) drive liver gene expression independently of the hepatic clock? Liver mRNA-seq, plus a public Bmal1-/- NR dataset (Atger 2015, GSE73554). Behaviour and temperature by telemetry; insulin tolerance test. Data: GSE118967.
**Sampling:** timed. Livers every 4 h for 24 h (6 time points), n = 3 per paradigm and time point. Bmal1-/- reanalysis: n = 4 per time point (WT), n = 2 (KO).

## Method (if any)
Rhythmicity by four statistical programs combined (MetaCycle among them); categories by rhythmic (R) or arrhythmic (A) status in NR, LB, AR: RRR, RRA, RAA, AAR. Differential rhythmicity with DODR (phase and amplitude, HANOVA); relative amplitude rAMP from MetaCycle.

## Key results
- RFI is flattened while behaviour stays rhythmic. Light-phase food share: 48.9% (AR), 23.8% (LB), 0% (NR); "mice fed AR remain behaviorally rhythmic and still exhibit a rhythmic drive to feed when active at night". Activity and drinking unaffected; body temperature normal in AR and 1 C lower by day in NR.
- Most cycling genes need RFI. "more than 70% of the cycling mouse liver transcriptome loses rhythmicity under arrhythmic feeding" (Summary); number of rhythmic genes tracks the amplitude of RFI. 324 RRR genes; RRA plus RAA 2,229 genes dampen in step with RFI.
- The core clock keeps oscillating without a feeding rhythm: "all core clock genes were in the RRR category", "the core molecular clock in the liver is not affected by RFI manipulation". Almost all clock genes are affected by Bmal1 knockout but not by RFI manipulation. Per2 is the exception: decreased amplitude and a phase advance of 2.8 h from LB under AR.
- Phase under feeding manipulation is small: most rhythmic genes are "phase-advanced by 1–2 h" in NR and AR relative to LB. In clock-deficient NR mice, RFI-driven genes keep a phase distribution close to WT, "globally phase-advanced by 1–2 h", i.e. feeding alone can set phase, the clock delays it by 1-2 h.
- Four BMAL1 targets (rate-limiting enzymes) lose rhythm under AR; Cpt1a and Alas1 cycle in NR Bmal1-/- mice, driven by RFI not the clock.

## Not addressed / limits
- No feeding inversion arm; no other tissue; no amplitude number for the core-clock genes beyond "small decrease" (RRR) and Per2 "decreased"; time resolution 4 h, n = 3.
- Core clock genes are called oscillating by rhythmicity category, not by a joint amplitude estimate; clock-gene amplitude ratios must be recomputed from GSE118967.
- AR mice still eat most of their food at night (light-phase share 48.9%), so RFI is dampened rather than abolished; activity and drinking stay rhythmic, so the clock may be held by a central or behavioural drive other than feeding (not measured: SCN, glucocorticoid).
- No test of the brake mechanism and no cue-interaction model.

## Link to us
- **Intro claim(s):** supports that the liver clock is self-sustained: clock genes keep rhythm, amplitude and phase (1-2 h) when feeding rhythm is removed, while the output transcriptome (70% of cyclers) follows feeding. Supports removal keeps amplitude, inversion flips phase (point 3 of the brief) for the liver.
- **Relation to clocksync:** data we could use (GSE118967 NR, LB, AR; GSE73554 NR WT vs Bmal1-/-). Feeds A297 (Greenwell AR vs LB: amplitude 0.724, shift 1.57 h are our numbers, not in the paper; the paper's 1-2 h and Per2 2.8 h are consistent with the shift).

## Ideas
- Recompute a joint clock-gene (Bmal1, Per, Cry, Rev-erb, Dbp) amplitude and phase for AR, NR, LB from GSE118967 and compare with the Manella DF inversion in the same units; the saturation reading predicts amplitude retained, phase near 0 for AR.
- Check that AR amplitude 0.724 is not driven by Per2 (the one clock gene with decreased amplitude).
- A tissue with higher b (WAT, kidney, heart) under arrhythmic feeding would discriminate self-sustained from feeding-driven: no such dataset here.
- Atger 2015 (GSE73554) gives liver of Bmal1-/- under NR: the clock-free reference for b near 1 in the liver.
