# kiessling2010adrenal
**Ref:** Kiessling, Eichele, Oster, J Clin Invest 2010. DOI 10.1172/JCI41192. Type: biology. Lab: Oster.
**Source:** full text

## Question and design
How fast do clock genes in different tissues re-entrain after a 6 h advance (also an 8 h delay for the drug experiment), and does the adrenal
clock (via glucocorticoids) regulate it? Mice, LD 12:12, wheel running (n = 9). Clock genes Per1, Per2, Dbp, Nr1d1, Arntl in SCN and
somatosensory cortex (in situ hybridisation) and in adrenal, kidney, liver, pancreas (qPCR), n = 3 animals per time point, sampled on day 0 and
days 2, 3, 4 and 8 or 12 after the shift; full daily profile at each day, then sine fit for peak time. Adrenal clock-deficient hosts
(Per2/Cry1 double mutant adrenal graft into adrenalectomised WT), fecal corticosterone at 4 h intervals, and timed metyrapone (MET) injections.
**Sampling:** timed, 24 h profile on each of 5-6 days; readout is PS50 (days to reach half of the shift), not hours of offset.

## Key results
- Behaviour: "This transition was completed after 8–9 days"; PS50 4.0 +/- 0.1 days (Fig. 1B).
- Within-tissue heterogeneity: "Unexpectedly, strong heterogeneity of entrainment kinetics was found not only between different organs, but also within the molecular clockwork of each tissue." (Abstract)
- SCN: Per genes fast, Dbp/Nr1d1/Arntl slower, closer to behaviour; full re-entrainment between day 3 (Per2) and day 8 (Nr1d1).
- Liver and kidney: "Per2 shifted fastest in kidney and liver, whereas Arntl shifted slowest" (Arntl PS50 4.1 +/- 0.2 days in liver, 3.5 +/- 0.2 in adrenal).
- Pancreas inverted: "In stark contrast, Per1 and Per2 were both slow to reset in the pancreas, whereas Nr1d1 adjusted fast"; PS50 4.8 +/- 0.7 days for Per1 and Per2, 3.5 +/- 0.4 for Nr1d1.
- Three groups: "Our data distinguish 3 different resetting groups" (CNS: SCN, cortex; peripheral kidney and liver; pancreas), adrenal intermediate.
- Adrenal clock removal speeds up, not slows: PS50 of activity "28.4% less in animals lacking an adrenal clock" (4.95 vs 3.55 days, p < 0.0001); corticosterone PS50 2.47 vs 1.57 days (-36%). Corticosterone shifts precede activity (r2 0.65-0.87).
- MET at day vs night phase shifted the corticosterone peak by about 1 h (ZT17 to ZT16 or ZT18) and either accelerated or slowed behavioural re-entrainment, also after an 8 h delay.
- Conclusion: "jet lag evokes a global desynchronization of clock gene expression rhythms that gradually returns to the robust alignment typical for the entrained state"

## Not addressed / limits
- Kinetics are per gene per tissue (PS50 in days), not an inter-tissue phase offset in hours at a given time; no tissue is referenced to another within animal (pooled, 3 animals per time point).
- Only advance for the multi-tissue analysis; delay only for behaviour and MET; males; mouse; no human or disease.
- Adrenal graft is a manipulation of the clock, not evidence of the actual adrenal role in intact animals beyond the MET timing.

## Link to us
- **Intro claim(s):** CLAIM PARTLY WRONG. refs_circadian_desynchrony.md says adrenal glucocorticoids needed for rapid resynchronisation of peripheral clocks. The paper finds the adrenal clock (via GC rhythms) is a brake: ablation accelerates re-entrainment by 28%, and MET timing can speed or slow it; peripheral clock resynchronisation itself was not tested with GC manipulation. "Tissues differ" is supported (gene by tissue kinetics).
- **Relation to clocksync:** background. Direct evidence that during recovery, genes of one tissue are at different phases, which challenges a single phase per compartment in transient states.

## Ideas
- Gene-wise partial shift test for the coherence gate D (see reddy2002differential): Arntl is the slowest and Per2 the fastest gene in liver (PS50 values in Fig. 3, not in the text). Not registered (no jet-lag hit in docs).
- Pancreas is the one tissue ordered differently from SCN/liver; a candidate target for tissue offset or coupling in GTEx pancreas (check A-numbers for pancreas in results_preliminary.md before registering).
