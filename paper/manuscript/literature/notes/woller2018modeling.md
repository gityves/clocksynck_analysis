# woller2018modeling
**Ref:** Woller, Gonze, Sci Rep 2018 (8:13641). DOI 10.1038/s41598-018-31804-9. Type: model. Lab: Gonze (ULB, Brussels).
**Source:** full text

## Question and design
Model of the mouse pancreatic beta-cell clock coupled to a glucose-insulin module, forced by two periodic signals, an SCN-driven cue (strength parameter cneur) and a nutrient (glucose) cue. Calibrated on pancreas data from Mukherji 2015 (night-fed, ad libitum); then the food cue is shifted by 12 h while the SCN cue stays fixed. Question: why do animals fed in the rest phase develop hypoinsulinemia and hyperglycemia although "everything is simply inverted"? 13 ODEs plus 2 for food; fit by Hooke and Jeeves (COPASI), solved in MATLAB.
**Sampling:** none new; the model is calibrated on published 24 h profiles (pancreas, mouse); no GEO data in the paper.

## Method (if any)
Local clock entrained by two cues: glucose acts on Per transcription; SCN cues (neural, taken as rhythmic activity) act on Rev-Erb. Phase shifts of each clock gene after inversion of the food cue, as a function of cneur. Comparison with two alternative architectures: no central clock (SCN rhythm replaced by its mean) and no clock at all.

## Key results
- Conflicting cues produce a gene-specific partial shift: "Per expression shifts by 11.63 h while the phase shift of the other clock genes (Bmal1, Rev-Erb and Cry mRNAs) lies between 7 and 9 h". Experimental motivation: some genes shift by 12 h while others "only shift by 8 h", in liver and pancreas (Mukherji data).
- Dependence on the SCN weight. "In the absence of SCN-driven cues, the whole clock shifts by 12 h." As SCN-driven cues strengthen (cneur from 0 to 1.55; default 1.4336), Per keeps shifting by "almost 12 h regardless of the strength of the SCN cue", while Bmal1 shifts less. The phase shift difference between Bmal1 and Per grows with SCN strength (Fig. 5A).
- Amplitude only appears in outputs: the partial shift misaligns glucose and clock-controlled exocytosis cues, giving reduced insulin secretion; Fig. 7 draws amplitude as arrow length for clock genes and metabolic variables. In the "No central clock" architecture "all genes shift by 12 h" and insulin is the same in night and day feeding. The model gives no number for clock-gene amplitude after inversion.
- Complete uncoupling fails: "the hypothesis of a complete uncoupling from the SCN fails to reproduce the differential phase shift in clock gene expression as well as the metabolic syndrome."
- Small shifts: phase differences between clock components are only slightly changed for a short advance or delay in food intake.
- Generality argued: liver, intestinal epithelial cells and PVN also show differential shifts in Mukherji; skin clock shifts "only partially or not at all" (Takahashi, ref 32).

## Not addressed / limits
- Cues are additive-style inputs to different genes (glucose on Per, SCN on Rev-Erb); the cell is a single oscillator with no amplitude-dependent saturation analysis; no population of oscillators or phase coherence.
- No removal of one cue (food or SCN) with amplitude output for the clock genes; amplitude after inversion not tabulated; no per-tissue weight beyond cneur.
- Glucocorticoids are left out on purpose: "the bimodal profile does not seem to affect much Per time profile", so the brake is not modelled.
- One tissue (pancreas), fitted to a single study; parameters fitted jointly with many free values.

## Link to us
- **Intro claim(s):** supports that conflicting light and food cues distort rather than simply invert the clock, with a tissue-specific SCN weight setting how far each component shifts; it matches the idea that the shift depends on a weight b (here cneur, where b = 0 gives a full 12 h shift). It does not predict amplitude loss of the clock genes from cancelling cues, nor retention under removal.
- **Relation to clocksync:** background and model competitor: a mechanistic two-cue model. Listed in the reference list of manella2021liverclock.

## Ideas
- Run this model across cneur and report clock-gene amplitude (not given in the paper): does amplitude dip at intermediate SCN weight (cancellation) as in the human muscle case (shift 2.7 h, amplitude 0.23), or stay flat?
- Add food removal (flat nutrient cue) to see phase and amplitude after removal: a direct test of the saturation reading versus the additive model failing in A297.
- Data to test per-tissue weight: Mukherji 2015 (pancreas, liver, intestinal epithelium, PVN, day vs night feeding; accessions not in this text, to retrieve); Manella 2021 GSE159135 (liver, WAT, lung).
