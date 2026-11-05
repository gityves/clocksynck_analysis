# damiola2000restricted
**Ref:** Damiola, Le Minh, Preitner, Kornmann, Fleury-Olela, Schibler, Genes Dev 2000. DOI 10.1101/gad.183500. Type: biology. Lab: Schibler.
**Source:** full text

## Question and design
Does time of feeding set the phase of peripheral clocks, and is the SCN affected? Mice on 12:12 LD (lights on 6 a.m.)
or in DD fed only in the light phase (6 a.m.-6 p.m.) or only in the dark phase for 6-9 days; sacrifice every 4 h over
one day (single profiles, animal numbers not given in the text). RNase protection for clock genes (Per1-3, Cry1, Rev-erba,
Dbp, Cyp2a5) in liver; in situ hybridisation for Per1/Per2 in SCN; Dbp in kidney, heart, pancreas; telemetric body temperature.
**Sampling:** timed; 4 h sampling over one full 24 h cycle per condition, plus 1/3/7 day entrainment series.

## Method (if any)
None; phase = time of peak mRNA of each profile.

## Key results
- Daytime feeding inverts liver phase: "the phases of all examined mRNA accumulation profiles differ by 8-12 h between mice fed during the day or during the night." (Results, Fig. 1)
- SCN unaffected, in LD and DD: "feeding time has little if any influence on circadian Per1 and Per2 expression in the SCN." (Fig. 1D, 3C)
- About 110 random genes (ADDER display): the large majority are constant; only circadian transcripts change phase (Fig. 2).
- Shifts are gradual (not immediate) over about a week, so they are clock-dependent, with transients on reversal and on food withdrawal (Figs 4, 5).
- Tissue kinetics differ: Dbp phase inverted in liver, kidney, heart and pancreas after 6 days of daytime feeding; on day 4 liver is not yet fully
  reversed and the other tissues show intermediate patterns (Fig. 6). "after 1 wk of daytime feeding, the phases of circadian gene expression are similar in all examined peripheral tissues."
- Body temperature profile changes strongly with daytime feeding (transient drop of up to 3 degrees C, Fig. 7).

## Not addressed / limits
- Mouse only; tissue phases from Dbp in three of four tissues (clock genes only in liver); in heart and pancreas amplitude is low on day 3 (cell-type mixture not resolved, as the authors note).
- Whole-tissue RNA: cannot separate cell types; no sex, age or genotype comparison.
- Mechanism of the feeding signal left open; SCN gene expression is Per1/Per2 only.
- No quantification of replicate variability.

## Link to us
- **Intro claim(s):** peripheral oscillators can be uncoupled from the SCN by feeding (Internal desynchrony).
  refs entry (inverts phase in liver, kidney, heart, pancreas; SCN unchanged) is supported but the multi-tissue part is Dbp only.
- **Relation to clocksync:** background. Provides tissue-specific kinetics: transient inter-organ offsets arise after any
  change of feeding time, which clocksync should read as phase offset arg z.

## Ideas
- Day-4 tissue-specific lags (liver ahead of kidney, heart, pancreas) as a benchmark for offsets between tissues in a
  simulation.
- Fasting and feeding state before surgery (hospital hours) could produce offsets; compare with brown2002rhythms,
  where temperature drives the same inversion.
