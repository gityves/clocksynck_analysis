# wefers2018circadian
**Ref:** Wefers J et al., PNAS 2018;115(30):7789-7794. DOI 10.1073/pnas.1722295115. Type: clinical (human misalignment protocol). Lab: Schrauwen (Maastricht), with Scheer, Kalsbeek.
**Source:** full text

## Question and design
Randomised crossover, n = 14 healthy lean young men (22.4 ± 2.8 y). Control (3 d aligned) vs 3.5 d misalignment induced by a 12 h rapid shift of the behavioural cycle in dim light (4 lx), in a respiration chamber. Read-outs: hyperinsulinemic euglycemic clamp, indirect calorimetry, core body temperature, plasma metabolites, vastus lateralis biopsies (qPCR of clock genes, Affymetrix microarray, mitochondrial respirometry).
**Sampling:** timed but sparse: "two skeletal muscle biopsies were obtained, at 8 AM and 8 PM" in each condition (7 AM / 7 PM in the text for qPCR), i.e. two time points per arm, 12 h apart. No 24 h time series.

## Method (if any)
None; paired t tests; multilevel PLS-DA and pathway enrichment on the microarray. No phase estimate is computed (only two time points).

## Key results
- Insulin sensitivity: muscle nonoxidative glucose disposal fell "23.7 ± 2.4 vs. 18.4 ± 1.4 mg/kg per minute" (control vs misalignment, P = 0.024; units in abstract vs results text differ, mg vs umol); hepatic insulin sensitivity unaffected.
- Different outputs shift differently (internal desynchrony): core body temperature "decreased in the late biological evening despite subjects remaining awake, indicating that the shifted behavioral cycle was misaligned relative to the central circadian clock".
- Muscle clock lags the behaviour: "mRNA expression levels of the core clock genes kept the same diurnal pattern as in the control condition, albeit blunted for PER and CRY, indicating that the core molecular clock in skeletal muscle indeed was misaligned relative to the inverted behavioral cycle." (BMAL1, CRY1, PER2 by qPCR).
- Authors' conclusion: "the molecular clock gene pathway had not adapted to the behavioral rhythm". So the clock genes did not shift by 12 h within 3.5 d; there is no measured phase shift of the muscle clock in hours.
- Transcriptome: most diurnally expressed genes tended to reverse with the behaviour while PPAR and fatty acid oxidation genes were up-regulated at 7 AM to 7 PM contrasts; mitochondrial state 3 respiration higher at 7 PM in misalignment.

## Not addressed / limits
- Only two biopsies per condition: figure caption states "the lack of difference in mRNA expression levels between 7 AM and 7 PM does not provide evidence that the measured parameter lacks 24 h rhythmicity." No amplitude or phase in hours.
- Only one peripheral tissue (muscle); other tissues (blood, adipose, SCN) are inferred from body temperature or earlier papers; no direct comparison of tissue clocks in the same person.
- n = 14 men; short protocol; no sleep quality data ("lack of information on sleep quality").
- No individual-level offsets reported.

## Link to us
- **Intro claim(s):** human misalignment alters the muscle clock and metabolism. Our refs line 76 (altered clock-gene expression phase shifts) is overstated/wrong: the clock genes did not shift; the paper shows the muscle clock failed to follow the shifted behaviour within 3.5 d, i.e. muscle and behaviour became misaligned. The grouped line 102 (two compartments in the same volunteers, muscle vs blood) is partly wrong: here the second reference is core body temperature, not a blood clock read-out.
- **Relation to clocksync:** background (supports that peripheral clock lags forced behavioural shift). The muscle microarray is not likely to give us a usable clock-gene phase.

## Ideas
- Test whether clocksync's cross-compartment offsets would show a muscle-behaviour lag if blood and muscle clock-gene data from the same donors under such protocols exist (e.g. Skene 2018 / Kervezee 2018 PBMC series).
- GTEx skeletal muscle vs blood offsets by donor time of death (day vs night deaths) should reproduce a muscle-lag signature: a pre-specifiable prediction.
- Missing reference: Morris 2015 PNAS (endogenous circadian system and misalignment affect glucose tolerance separately), cited as ref 3 here.
